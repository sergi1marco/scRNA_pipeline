# A Comprehensive and Reproducible R Pipeline for Single-Cell RNA Sequencing Data Analysis

# Set RETICULATE_CONDA_PATH and disable auto-installation
Sys.setenv(RETICULATE_CONDA_PATH = "C:/ProgramData/miniconda3/condabin/conda.bat")
Sys.setenv(RETICULATE_PYTHON_FALLBACK = "FALSE")
options(reticulate.conda_auto_install = FALSE)
options(reticulate.install_miniforge = FALSE)

# Point reticulate to the Python executable in your new, space-free environment
reticulate::use_python(python = "C:/miniconda3/envs/r_magic_env/python.exe", required = TRUE)

# Global options
options(future.globals.maxSize = 8000 * 1024^2)
options(scipen = 100)

# --- 0. Library Loading ---
library(Seurat)
library(SeuratObject)
library(SingleCellExperiment) # Explicitly added this to the installation script for robustness
library(scater)
library(scran)
library(DropletUtils)
library(BiocSingular)
library(PCAtools)
library(SingleR)
library(celldex)
library(dittoSeq)
library(TENxIO)
library(singleCellTK)
library(MAST)
library(limma)
library(edgeR)
library(ComplexHeatmap)
library(tidyverse)
library(cowplot)
library(patchwork)
library(ggrepel)
library(viridis)
library(data.table)
library(ggraph)
library(igraph)
library(RColorBrewer)
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(fgsea)
library(msigdbr)
library(circlize)
library(plotly)
library(reticulate)
library(remotes)
library(EnhancedVolcano)
library(nebula)
library(presto)
library(sceasy)
library(zellkonverter)

message("\n--- All required libraries loaded successfully. ---")

# --- 1. Interactive Directory Selection & Pipeline Start Decision ---

# Function to select a directory
select_data_directory <- function() {
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    message("Using RStudio API for directory selection...")
    data_path <- rstudioapi::selectDirectory(
      caption = "Select the working directory (or where your 10x Genomics data / saved Seurat objects are)",
      path = getwd()
    )
  } else {
    message("Using base R for directory selection...")
    data_path <- choose.dir(caption = "Select the working directory (or where your 10x Genomics data / saved Seurat objects are)")
  }
  
  if (is.null(data_path) || data_path == "") {
    stop("Directory selection cancelled. Please select a valid directory to proceed.")
  }
  return(data_path)
}

data_dir <- select_data_directory()
setwd(data_dir)

message(paste("Working directory set to:", getwd()))

seurat_object <- NULL
start_step <- 0

saved_seurat_files <- list.files(data_dir, pattern = "\\.rds$", full.names = TRUE, recursive = TRUE)
saved_seurat_files <- saved_seurat_files[grepl("seurat_object", basename(saved_seurat_files))]

if (length(saved_seurat_files) > 0) {
  message("\n--- Found existing Seurat object(s): ---")
  for (i in seq_along(saved_seurat_files)) {
    message(paste0(i, ": ", saved_seurat_files[i]))
  }
  message("------------------------------------------------------------------")
  
  choice_made <- FALSE
  while (!choice_made) {
    user_choice_file <- readline(prompt = "Enter the number of the Seurat object to load, or '0' to start from raw data: ")
    user_choice_file <- as.integer(user_choice_file)
    
    if (user_choice_file == 0) {
      message("Proceeding to import raw data.")
      break
    } else if (user_choice_file > 0 && user_choice_file <= length(saved_seurat_files)) {
      seurat_object <- readRDS(saved_seurat_files[user_choice_file])
      message(paste("Loaded Seurat object from:", saved_seurat_files[user_choice_file]))
      
      message("\n--- Pipeline Steps: ---")
      message("1: Quality Control and Filtering")
      message("2: Normalization and Scaling")
      message("3: Dimensionality Reduction (PCA and UMAP)")
      message("4: Visualization (Initial UMAPs, etc.)")
      message("5: Find Marker Genes for Clusters")
      message("6: Cell Type Annotation & CAF Identification (SingleR)")
      message("7: Differential Expression Analysis")
      message("8: Malignancy Detection (scMalignantFinder via reticulate)")
      message("9: Interactive Gene Expression Analysis & Visualization")
      message("10: Final Data Saving")
      message("-----------------------")
      
      user_choice_step <- readline(prompt = "From which step would you like to resume the pipeline (1-10)? ")
      start_step <- as.integer(user_choice_step)
      
      if (is.na(start_step) || start_step < 1 || start_step > 10) {
        warning("Invalid step chosen. Please enter a number between 1 and 10. Reloading Seurat object and re-asking.")
        seurat_object <- NULL
        next
      } else {
        step_names <- c("QC", "Normalization", "DimRed", "Visualization", "Markers", "Annotation_CAF", "DE", "scMalignantFinder_Detection", "Interactive_Gene_Analysis", "Final_Save")
        message(paste("Starting pipeline from step:", step_names[start_step]))
        choice_made <- TRUE
      }
    } else {
      message("Invalid file choice. Please enter a valid number or '0'." )
    }
  }
} else {
  message("\nNo existing Seurat object found. Proceeding to import raw data.")
}

current_seurat_object <- seurat_object

# --- 2. Data Acquisition and Initial Processing ---
if (is.null(seurat_object)) {
  # Function to import 10x Genomics data
  import_10x_data <- function(data_directory) {
    message(paste("Scanning directory for 10x Genomics data files:", data_directory))
    
    h5_files <- list.files(data_directory, pattern = "\\.h5$", full.names = TRUE, recursive = TRUE)
    mtx_dirs <- c()
    
    sub_dirs <- list.dirs(data_directory, recursive = TRUE)
    for (s_dir in sub_dirs) {
      if (file.exists(file.path(s_dir, "matrix.mtx.gz")) &&
          file.exists(file.path(s_dir, "features.tsv.gz")) &&
          file.exists(file.path(s_dir, "barcodes.tsv.gz"))) {
        mtx_dirs <- c(mtx_dirs, s_dir)
      }
    }
    
    selected_path <- NULL
    data_type <- NULL
    
    if (length(h5_files) > 0) {
      if (length(h5_files) > 1) {
        warning("Multiple .h5 files found. Using the first one detected (alphabetically): ", h5_files[1])
      }
      selected_path <- h5_files[1]
      data_type <- "h5"
    } else if (length(mtx_dirs) > 0) {
      if (length(mtx_dirs) > 1) {
        warning("Multiple MTX data directories found. Using the first one detected (alphabetically): ", mtx_dirs[1])
      }
      selected_path <- mtx_dirs[1]
      data_type <- "mtx"
    } else {
      stop("No 10x Genomics H5 file or MTX data directory found. Please ensure your data is correctly placed.")
    }
    
    if (data_type == "h5") {
      message(paste("Importing data from H5 file:", selected_path))
      seurat_obj_raw <- Read10X_h5(selected_path)
      if (is.list(seurat_obj_raw)) {
        if ("Gene Expression" %in% names(seurat_obj_raw)) {
          seurat_obj_raw <- seurat_obj_raw[["Gene Expression"]]
          message("Extracted 'Gene Expression' assay from H5 file.")
        } else if (length(seurat_obj_raw) == 1) {
          seurat_obj_raw <- seurat_obj_raw[[1]]
          message(paste("Extracted the only assay (", names(seurat_obj_raw)[1], ") from H5 file.", sep=""))
        } else {
          stop("H5 file contains multiple data types, but 'Gene Expression' was not found. Please select the correct assay.")
        }
      }
    } else if (data_type == "mtx") {
      message(paste("Importing data from MTX folder:", selected_path))
      seurat_obj_raw <- Read10X(data.dir = selected_path)
    }
    
    seurat_object <- CreateSeuratObject(counts = seurat_obj_raw, project = "SingleCellAnalysis")
    message("Seurat object created successfully.")
    
    return(seurat_object)
  }
  
  seurat_object <- import_10x_data(data_directory = data_dir)
}

current_seurat_object <- seurat_object

# --- 3. Quality Control and Filtering ---
if (start_step <= 1) {
  perform_qc_and_filter <- function(seurat_obj, min_umis = 500, max_umis = 100000,
                                    max_percent_mt = 10, min_features = 200, min_cells = 3) {
    message("Performing Quality Control and Filtering...")
    
    seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
    
    pdf("QC_Plots.pdf", width = 10, height = 8)
    print(VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1))
    print(FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.mt"))
    print(FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA"))
    dev.off()
    message("QC plots saved to 'QC_Plots.pdf'. Review to refine thresholds.")
    
    seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > min_features &
                           nCount_RNA > min_umis &
                           nCount_RNA < max_umis &
                           percent.mt < max_percent_mt)
    
    message(paste("QC and filtering complete. Remaining cells:", ncol(seurat_obj)))
    if (ncol(seurat_obj) == 0) {
      stop("No cells remaining after QC. Adjust filtering thresholds.")
    }
    
    return(seurat_obj)
  }
  
  current_seurat_object <- perform_qc_and_filter(current_seurat_object,
                                                 min_umis = 500,
                                                 max_umis = 100000,
                                                 max_percent_mt = 10,
                                                 min_features = 200,
                                                 min_cells = 3)
} else {
  message("Skipping Quality Control and Filtering.")
}

# --- 4. Normalization and Scaling ---
if (start_step <= 2) {
  perform_normalization_and_scaling <- function(seurat_obj) {
    message("Normalizing and Scaling data...")
    seurat_obj <- SCTransform(seurat_obj, vars.to.regress = "percent.mt", verbose = FALSE)
    message("Normalization and scaling complete.")
    return(seurat_obj)
  }
  current_seurat_object <- perform_normalization_and_scaling(current_seurat_object)
} else {
  message("Skipping Normalization and Scaling.")
}

# --- 5. Dimensionality Reduction (PCA and UMAP) ---
if (start_step <= 3) {
  perform_dimensionality_reduction <- function(seurat_obj, dims = 1:30, resolution = 0.5) {
    message("Performing Dimensionality Reduction (PCA and UMAP)...")
    seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj))
    
    pdf("PCA_ElbowPlot.pdf", width = 8, height = 6)
    print(ElbowPlot(seurat_obj))
    dev.off()
    message("PCA complete. See 'PCA_ElbowPlot.pdf' to select appropriate number of dimensions.")
    
    seurat_obj <- FindNeighbors(seurat_obj, dims = dims)
    seurat_obj <- FindClusters(seurat_obj, resolution = resolution)
    
    seurat_obj <- RunUMAP(seurat_obj, dims = dims)
    message("UMAP and Clustering complete.")
    return(seurat_obj)
  }
  current_seurat_object <- perform_dimensionality_reduction(current_seurat_object, dims = 1:30, resolution = 0.8)
} else {
  message("Skipping Dimensionality Reduction.")
}

# --- 6. Visualization ---
if (start_step <= 4) {
  perform_visualization <- function(seurat_obj) {
    message("Generating visualizations...")
    
    pdf("UMAP_Clusters.pdf", width = 8, height = 7)
    print(DimPlot(seurat_obj, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend())
    dev.off()
    message("UMAP plot by clusters saved to 'UMAP_Clusters.pdf'.")
    
    marker_genes <- c("CD3D", "CD8A", "MS4A1", "CD14", "FCGR3A", "NKG7", "PPBP", "S100A8")
    marker_genes_present <- marker_genes[marker_genes %in% rownames(GetAssayData(seurat_obj, assay = "SCT", slot = "data"))]
    
    if (length(marker_genes_present) > 0) {
      pdf("UMAP_FeaturePlots.pdf", width = 12, height = 10)
      print(FeaturePlot(seurat_obj, features = marker_genes_present, reduction = "umap", ncol = 3, slot = "data"))
      dev.off()
      message("UMAP feature plots for example marker genes saved to 'UMAP_FeaturePlots.pdf'.")
    } else {
      message("No example marker genes found. Skipping FeaturePlot generation.")
    }
    message("All standard visualizations complete.")
    return(seurat_obj)
  }
  current_seurat_object <- perform_visualization(current_seurat_object)
} else {
  message("Skipping Visualization.")
}

# --- 7. Find Marker Genes for Clusters ---
if (start_step <= 5) {
  find_cluster_markers <- function(seurat_obj) {
    message("Finding marker genes for each cluster...")
    
    cluster_markers <- FindAllMarkers(seurat_obj,
                                      assay = "SCT",
                                      min.pct = 0.25,
                                      logfc.threshold = 0.25,
                                      only.pos = TRUE,
                                      verbose = FALSE)
    
    if (nrow(cluster_markers) == 0) {
      message("WARNING: No differentially expressed genes found. Adjust 'min.pct' or 'logfc.threshold'.")
      top10_markers <- data.frame(cluster=character(), gene=character(), avg_log2FC=numeric())
      write.csv(cluster_markers, "All_Cluster_Markers.csv", row.names = FALSE)
      write.csv(top10_markers, "Top10_Cluster_Markers.csv", row.names = FALSE)
      return(list(all_markers = cluster_markers, top10_markers = top10_markers))
    }
    
    top10_markers <- cluster_markers %>%
      group_by(cluster) %>%
      slice_max(n = 10, order_by = avg_log2FC, with_ties = FALSE)
    
    write.csv(cluster_markers, "All_Cluster_Markers.csv", row.names = FALSE)
    write.csv(top10_markers, "Top10_Cluster_Markers.csv", row.names = FALSE)
    message("Marker gene finding complete. 'All_Cluster_Markers.csv' and 'Top10_Cluster_Markers.csv' saved.")
    return(list(all_markers = cluster_markers, top10_markers = top10_markers))
  }
  current_seurat_object@misc$all_markers <- find_cluster_markers(current_seurat_object)$all_markers
  current_seurat_object@misc$top10_markers <- find_cluster_markers(current_seurat_object)$top10_markers
} else {
  message("Skipping Find Marker Genes.")
  if (file.exists("All_Cluster_Markers.csv")) {
    current_seurat_object@misc$all_markers <- read.csv("All_Cluster_Markers.csv")
  }
  if (file.exists("Top10_Cluster_Markers.csv")) {
    current_seurat_object@misc$top10_markers <- read.csv("Top10_Cluster_Markers.csv")
  }
}
# --- 8. Cell Type Annotation & CAF Identification (SingleR and custom markers) ---
if (start_step <= 6) {
  perform_cell_type_annotation_and_caf <- function(seurat_obj) {
    message("Performing cell type annotation using SingleR...")
    
    message("\n--- SINGLE R: REFERENCE SELECTION ---")
    message("Choose a reference dataset for cell type annotation:")
    message("1: HumanPrimaryCellAtlasData")
    message("2: NovershternHematopoieticStemCellAtlas")
    message("3: BlueprintEncodeData")
    message("4: MonacoImmuneData")
    message("5: HPCA + BlueprintEncode (Combined)")
    
    ref_choice <- readline(prompt = "Enter your choice (1-5): ")
    ref_choice <- as.integer(ref_choice)
    
    ref <- NULL
    tryCatch({
      if (ref_choice == 1) {
        ref <- celldex::HumanPrimaryCellAtlasData()
        message("HumanPrimaryCellAtlasData reference loaded.")
        csv_file_path <- "HPCA_tissue_cell_types.csv"
        
        if (file.exists(csv_file_path)) {
          message(paste0("Reading cell type mapping from ", csv_file_path, "..."))
          tissue_cell_types_map <- data.table::fread(csv_file_path)
          
          if (!("Tissue" %in% colnames(tissue_cell_types_map)) || !("CellTypes" %in% colnames(tissue_cell_types_map))) {
            warning(paste0("CSV file '", csv_file_path, "' must contain 'Tissue' and 'CellTypes' columns. Using full HPCA data."))
            final_cell_types_to_filter <- unique(ref$label.main)
          } else {
            available_tissues <- unique(tissue_cell_types_map$Tissue)
            
            message("\nAvailable tissues from HPCA_tissue_cell_types.csv:")
            for (i in seq_along(available_tissues)) {
              message(paste0(i, ": ", available_tissues[i]))
            }
            
            user_tissue_choice_idx <- as.numeric(readline(prompt = "Enter the number corresponding to your tissue of interest: "))
            
            if (!is.na(user_tissue_choice_idx) && user_tissue_choice_idx %in% seq_along(available_tissues)) {
              chosen_tissue <- available_tissues[user_tissue_choice_idx]
              message(paste0("You selected: ", chosen_tissue))
              
              cell_types_string <- tissue_cell_types_map[Tissue == chosen_tissue]$CellTypes
              retrieved_cell_types <- unlist(strsplit(cell_types_string, split = ",", fixed = TRUE))
              retrieved_cell_types <- trimws(retrieved_cell_types)
              
              if (!("PBMCs" %in% retrieved_cell_types)) {
                retrieved_cell_types <- c(retrieved_cell_types, "PBMCs")
                message("Added 'PBMCs' to the cell type list for filtering.")
              }
              
              final_cell_types_to_filter <- retrieved_cell_types[retrieved_cell_types %in% unique(ref$label.main)]
              
              if (length(final_cell_types_to_filter) > 0) {
                ref <- ref[, ref$label.main %in% final_cell_types_to_filter]
                message(paste0("HumanPrimaryCellAtlasData filtered to include: ", paste(final_cell_types_to_filter, collapse = ", ")))
              } else {
                warning("No matching cell types found. Using full HPCA data for annotation.")
                final_cell_types_to_filter <- unique(ref$label.main)
              }
            } else {
              warning("Invalid tissue selection. Using full HumanPrimaryCellAtlasData for annotation.")
              final_cell_types_to_filter <- unique(ref$label.main)
            }
          }
        } else {
          warning(paste0("CSV file '", csv_file_path, "' not found. Using full HumanPrimaryCellAtlasData for annotation."))
          final_cell_types_to_filter <- unique(ref$label.main)
        }
      } else if (ref_choice == 2) {
        ref <- celldex::NovershternHematopoieticStemCellAtlas()
        message("NovershternHematopoieticStemCellAtlas reference loaded.")
      } else if (ref_choice == 3) {
        ref <- celldex::BlueprintEncodeData()
        message("BlueprintEncodeData reference loaded.")
      } else if (ref_choice == 4) {
        ref <- celldex::MonacoImmuneData()
        message("MonacoImmuneData reference loaded.")
      } else if (ref_choice == 5) {
        ref_hpca <- celldex::HumanPrimaryCellAtlasData()
        ref_bp <- celldex::BlueprintEncodeData()
        ref <- c(ref_hpca, ref_bp)
        message("HPCA and BlueprintEncodeData references combined.")
      } else {
        stop("Invalid reference choice. Using HumanPrimaryCellAtlasData as default.")
        ref <- celldex::HumanPrimaryCellAtlasData()
      }
    }, error = function(e) {
      stop("Failed to load selected celldex reference. Error: ", e$message)
    })
    
    sce <- as.SingleCellExperiment(seurat_obj)
    predictions <- SingleR(test = sce,
                           ref = ref,
                           labels = ref$label.main)
    
    seurat_obj$singler_labels <- predictions$pruned.labels
    message("Initial cell type annotation complete.")
    
    seurat_obj$combined_cell_types <- seurat_obj$singler_labels
    message("Initialized 'combined_cell_types' with SingleR labels.")
    
    message("Attempting to identify Cancer-Associated Fibroblasts (CAFs)...")
    
    label_counts_per_cluster <- as.data.frame(table(seurat_obj$seurat_clusters, seurat_obj$singler_labels))
    colnames(label_counts_per_cluster) <- c("cluster", "label", "count")
    
    most_frequent_labels <- label_counts_per_cluster %>%
      group_by(cluster) %>%
      slice_max(order_by = count, n = 1) %>%
      ungroup()
    
    fibroblast_clusters <- most_frequent_labels$cluster[grepl("fibroblast|stromal|mesenchymal", most_frequent_labels$label, ignore.case = TRUE)]
    fibroblast_clusters <- as.character(fibroblast_clusters)
    
    if (length(fibroblast_clusters) > 0) {
      message(paste("Potential fibroblast/stromal clusters:", paste(fibroblast_clusters, collapse = ", ")))
      
      # --- START OF MODIFIED CAF MARKER RETRIEVAL ---
      caf_marker_file <- "CAF_markers_all.csv" # Define the CAF marker file name
      
      if (!file.exists(caf_marker_file)) {
        stop(paste0("CAF marker file '", caf_marker_file, "' not found in the working directory. Please provide it to use custom CAF markers."))
      }
      
      caf_markers_df <- read.csv(caf_marker_file)
      
      if (!("cancer_type" %in% colnames(caf_markers_df)) || !("caf_markers" %in% colnames(caf_markers_df))) {
        stop(paste0("The file '", caf_marker_file, "' must contain 'cancer_type' and 'caf_markers' columns."))
      }
      
      available_caf_cancer_types <- unique(caf_markers_df$cancer_type)
      
      if (length(available_caf_cancer_types) == 0) {
        stop("No cancer types found in 'CAF_markers_all.csv'. Cannot proceed with custom CAF markers.")
      }
      
      message("\nAvailable Cancer Types for CAF Markers from 'CAF_markers_all.csv':")
      for (i in seq_along(available_caf_cancer_types)) {
        message(paste0(i, ": ", available_caf_cancer_types[i]))
      }
      
      caf_cancer_type_choice_idx <- as.integer(readline(prompt = "Enter the number corresponding to the cancer type for CAF markers: "))
      
      if (is.na(caf_cancer_type_choice_idx) || !(caf_cancer_type_choice_idx %in% seq_along(available_caf_cancer_types))) {
        stop("Invalid cancer type choice for CAF markers. Please run the script again and select a valid number.")
      }
      
      chosen_caf_cancer_type <- available_caf_cancer_types[caf_cancer_type_choice_idx]
      message(paste0("You selected for CAF markers: ", chosen_caf_cancer_type))
      
      selected_caf_markers_string <- caf_markers_df$caf_markers[caf_markers_df$cancer_type == chosen_caf_cancer_type]
      
      if (length(selected_caf_markers_string) == 0) {
        stop(paste0("No CAF markers found for cancer type '", chosen_caf_cancer_type, "' in the CSV file."))
      }
      
      # Split the comma-separated string into a vector of gene names
      custom_caf_markers <- unlist(strsplit(selected_caf_markers_string, split = ",\\s*"))
      custom_caf_markers <- toupper(trimws(custom_caf_markers)) # Ensure uppercase and no leading/trailing spaces
      
      message(paste0("Using custom CAF markers for '", chosen_caf_cancer_type, "': ", paste(custom_caf_markers, collapse = ", ")))
      
      # Use the dynamically loaded custom_caf_markers
      caf_markers_present <- custom_caf_markers[custom_caf_markers %in% rownames(seurat_obj)] #
      # --- END OF MODIFIED CAF MARKER RETRIEVAL ---
      
      # The rest of the CAF scoring and classification logic remains the same
      if (length(caf_markers_present) > 0) {
        seurat_obj <- AddModuleScore(seurat_obj, features = list(CAF_Signature = caf_markers_present), name = "CAF_Score")
        message("Added CAF signature score.")
        
        seurat_obj$is_CAF <- "Non_CAF"
        
        if (length(fibroblast_clusters) > 0 && length(caf_markers_present) > 0) {
          cells_in_fib_clusters <- colnames(seurat_obj)[seurat_obj$seurat_clusters %in% fibroblast_clusters]
          
          if (length(cells_in_fib_clusters) > 0) {
            caf_score_for_fibroblasts <- seurat_obj$CAF_Score1[cells_in_fib_clusters]
            caf_threshold_for_classification_fib <- quantile(caf_score_for_fibroblasts, 0.75, na.rm = TRUE)
            
            seurat_obj$is_CAF[seurat_obj$seurat_clusters %in% fibroblast_clusters & seurat_obj$CAF_Score1 > caf_threshold_for_classification_fib] <- "CAF"
            
            message(paste("Cells classified as CAF (score >", round(caf_threshold_for_classification_fib, 2), "):", sum(seurat_obj$is_CAF == "CAF")))
          } else {
            message("No cells in identified fibroblast clusters for detailed CAF classification.")
          }
        }
        
        seurat_obj$is_CAF <- factor(seurat_obj$is_CAF, levels = c("Non_CAF", "CAF"))
        
        seurat_obj$combined_cell_types <- seurat_obj$singler_labels
        if ("CAF" %in% levels(seurat_obj$is_CAF) && any(seurat_obj$is_CAF == "CAF")) {
          seurat_obj$combined_cell_types[seurat_obj$is_CAF == "CAF"] <- "Cancer-Associated Fibroblast (CAF)"
        }
        seurat_obj$combined_cell_types <- factor(seurat_obj$combined_cell_types)
        message("Updated 'combined_cell_types' with integrated CAF status.")
        
        if ("CAF_Score1" %in% colnames(seurat_obj@meta.data)) {
          pdf("UMAP_CAF_Score.pdf", width = 8, height = 7)
          print(FeaturePlot(seurat_obj, features = "CAF_Score1", reduction = "umap", pt.size = 0.5, slot = "data", order = TRUE))
          dev.off()
          message("UMAP plot of CAF signature score saved to 'UMAP_CAF_Score.pdf'.")
        }
        
        message("Review 'UMAP_CAF_Score.pdf' and individual marker expression to confirm CAF populations.")
        
      } else {
        message("No CAF markers present in data after filtering. Skipping detailed CAF scoring and classification.")
      }
    } else {
      message("No clusters identified as fibroblast/stromal. Skipping detailed CAF analysis.")
    }
    
    pdf("UMAP_CellTypeAnnotation.pdf", width = 10, height = 8)
    print(DimPlot(seurat_obj,
                  reduction = "umap",
                  group.by = "combined_cell_types",
                  label = TRUE,
                  repel = TRUE) +
            ggtitle("UMAP by Cell Type Annotation with Integrated CAF Status"))
    dev.off()
    message("Cell type annotation complete. UMAP plot by predicted cell types saved to 'UMAP_CellTypeAnnotation.pdf'.")
    
    return(seurat_obj)
  }
  current_seurat_object <- perform_cell_type_annotation_and_caf(current_seurat_object)
} else {
  message("Skipping Cell Type Annotation & CAF Identification.")
}
# --- INTERMEDIATE DATA SAVING (after Visualization and before Marker Gene Finding) ---
if (start_step < 7) { # Only execute if starting from a step earlier than (or equal to) step 4 (Visualization)
  message("\n--- INTERMEDIATE DATA SAVING (after Visualization and before Marker Gene Finding) ---")
  output_filename <- "normalised_seurat_object.rds"
  saveRDS(current_seurat_object, file = output_filename)
  message(paste0("Seurat object saved to '", output_filename, "' after visualization."))
} else {
  message("Skipping intermediate data saving as pipeline started from step 5 or later.")
}

# --- 9. Differential Expression Analysis ---
if (start_step <= 7) {
  perform_differential_expression <- function(seurat_obj, ident1, ident2, assay = "SCT") {
    message(paste("Performing differential expression analysis between cluster", ident1, "and cluster", ident2, "..."))
    
    if (!as.character(ident1) %in% levels(seurat_obj) || !as.character(ident2) %in% levels(seurat_obj)) {
      warning(paste("One or both specified identities (", ident1, ",", ident2, ") not found. Skipping DE analysis.", sep=""))
      return(NULL)
    }
    
    DefaultAssay(seurat_obj) <- assay
    seurat_obj <- SetIdent(seurat_obj, value = "seurat_clusters")
    
    de_results <- FindMarkers(seurat_obj,
                              ident.1 = ident1,
                              ident.2 = ident2,
                              min.pct = 0.25,
                              logfc.threshold = 0.25,
                              verbose = FALSE)
    
    output_file <- paste0("DE_Cluster_", ident1, "_vs_Cluster_", ident2, ".csv")
    write.csv(de_results, output_file, row.names = TRUE)
    message(paste("Differential expression analysis complete. Results saved to", output_file))
    
    if (requireNamespace("EnhancedVolcano", quietly = TRUE) && !is.null(de_results) && nrow(de_results) > 0) {
      pdf(paste0("VolcanoPlot_Cluster_", ident1, "_vs_Cluster_", ident2, ".pdf"), width = 10, height = 10)
      print(EnhancedVolcano(de_results,
                            lab = rownames(de_results),
                            x = 'avg_log2FC',
                            y = 'p_val_adj',
                            title = paste('Cluster', ident1, 'vs Cluster', ident2),
                            pCutoff = 0.05,
                            FCcutoff = 0.5,
                            pointSize = 2.0,
                            labSize = 4.0))
      dev.off()
      message(paste("Volcano plot saved to VolcanoPlot_Cluster_", ident1, "_vs_Cluster_", ident2, ".pdf", sep=""))
    } else {
      message("Skipping Volcano plot: EnhancedVolcano not installed, no DE results, or no significant genes.")
    }
    
    return(de_results)
  }
} else {
  message("Skipping Differential Expression Analysis.")
}
# --- 10. Malignancy Detection (scMalignantFinder or Custom Markers) ---
if (start_step <= 8) {
  message("\n--- MALIGNANCY DETECTION METHOD SELECTION ---")
  message("Choose a method for malignancy annotation:")
  message("1: Custom Marker List (from 'malignancy_markers_all.csv')")
  message("2: scMalignantFinder (Automated Classifier)")
  
  malignancy_method_choice <- readline(prompt = "Enter your choice (1 or 2): ")
  malignancy_method_choice <- as.integer(malignancy_method_choice)
  
  if (is.na(malignancy_method_choice) || !(malignancy_method_choice %in% c(1, 2))) {
    warning("Invalid choice. Defaulting to scMalignantFinder.")
    malignancy_method_choice <- 2
  }
  
  if (malignancy_method_choice == 1) {
    message("\n--- MALIGNANCY DETECTION: CUSTOM MARKER LIST ---")
    marker_file <- "malignancy_markers_all.csv"
    
    if (!file.exists(marker_file)) {
      stop(paste0("Malignancy marker file '", marker_file, "' not found in the working directory. Please provide it to use custom markers."))
    }
    
    malignancy_markers_df <- read.csv(marker_file)
    
    if (!("cancer_type" %in% colnames(malignancy_markers_df)) || !("malignancy_markers" %in% colnames(malignancy_markers_df))) {
      stop(paste0("The file '", marker_file, "' must contain 'cancer_type' and 'malignancy_markers' columns."))
    }
    
    available_cancer_types <- unique(malignancy_markers_df$cancer_type)
    
    if (length(available_cancer_types) == 0) {
      stop("No cancer types found in 'malignancy_markers_all.csv'. Cannot proceed with custom markers.")
    }
    
    message("\nAvailable Cancer Types from 'malignancy_markers_all.csv':")
    for (i in seq_along(available_cancer_types)) {
      message(paste0(i, ": ", available_cancer_types[i]))
    }
    
    cancer_type_choice_idx <- as.integer(readline(prompt = "Enter the number corresponding to the cancer type you want to analyze: "))
    
    if (is.na(cancer_type_choice_idx) || !(cancer_type_choice_idx %in% seq_along(available_cancer_types))) {
      stop("Invalid cancer type choice. Please run the script again and select a valid number.")
    }
    
    chosen_cancer_type <- available_cancer_types[cancer_type_choice_idx]
    message(paste0("You selected: ", chosen_cancer_type))
    
    selected_markers_string <- malignancy_markers_df$malignancy_markers[malignancy_markers_df$cancer_type == chosen_cancer_type]
    
    if (length(selected_markers_string) == 0) {
      stop(paste0("No markers found for cancer type '", chosen_cancer_type, "' in the CSV file."))
    }
    
    # Split the comma-separated string into a vector of gene names
    custom_malignancy_markers <- unlist(strsplit(selected_markers_string, split = ",\\s*"))
    custom_malignancy_markers <- toupper(trimws(custom_malignancy_markers)) # Ensure uppercase and no leading/trailing spaces
    
    message(paste0("Using custom malignancy markers: ", paste(custom_malignancy_markers, collapse = ", ")))
    
    # Check which markers are actually present in the Seurat object
    markers_present <- custom_malignancy_markers[custom_malignancy_markers %in% rownames(GetAssayData(current_seurat_object, assay = "SCT", slot = "data"))]
    
    if (length(markers_present) == 0) {
      warning("None of the specified custom malignancy markers were found in the Seurat object's gene list. Cannot calculate a score.")
      current_seurat_object$scm_malignancy <- "Undetermined"
      current_seurat_object$scm_scores <- 0 # Default score
    } else {
      message(paste0("Found ", length(markers_present), " of ", length(custom_malignancy_markers), " custom malignancy markers in data."))
      
      # Calculate a module score for the custom markers
      current_seurat_object <- AddModuleScore(current_seurat_object, features = list(Malignancy_Signature = markers_present), name = "Custom_Malignancy_Score")
      
      # You might need to adjust the score column name based on AddModuleScore output (e.g., Custom_Malignancy_Score1)
      score_col_name <- paste0("Custom_Malignancy_Score1") # AddModuleScore appends "1" by default
      
      if (!(score_col_name %in% colnames(current_seurat_object@meta.data))) {
        stop(paste0("Expected score column '", score_col_name, "' not found after AddModuleScore."))
      }
      
      # Determine a threshold for malignancy based on the score
      # This is a heuristic; you might want to adjust it based on your data.
      # For simplicity, let's say cells in the top quartile of the score are malignant.
      malignancy_threshold <- quantile(current_seurat_object@meta.data[[score_col_name]], 0.75, na.rm = TRUE)
      message(paste0("Using a malignancy score threshold of > ", round(malignancy_threshold, 3), " to classify cells as Malignant."))
      
      current_seurat_object$scm_malignancy <- ifelse(current_seurat_object@meta.data[[score_col_name]] > malignancy_threshold, "Malignant", "Normal")
      current_seurat_object$scm_scores <- current_seurat_object@meta.data[[score_col_name]]
      
      current_seurat_object$scm_malignancy <- factor(current_seurat_object$scm_malignancy, levels = c("Normal", "Malignant"))
      message("Malignancy detection by custom markers complete. Malignancy status and scores added.")
      
      # Plot 1: UMAP colored by custom marker malignancy status
      output_filename_malig <- "UMAP_Malignancy_Status_CustomMarkers.pdf"
      pdf(output_filename_malig, width = 8, height = 7)
      plot_malig_status <- DimPlot(current_seurat_object,
                                   reduction = "umap",
                                   group.by = "scm_malignancy",
                                   label = TRUE,
                                   repel = TRUE,
                                   pt.size = 0.5) +
        ggtitle(paste0("UMAP by Malignancy Status (Custom Markers: ", chosen_cancer_type, ")"))
      print(plot_malig_status)
      dev.off()
      message(paste0("UMAP plot by custom marker inferred malignancy saved to '", output_filename_malig, "'."))
      
      # --- NEW INTEGRATED PLOT FOR CUSTOM MARKERS (with Legend) ---
      message("Generating integrated UMAP plot for Cell Type and Malignancy Status (Custom Markers) with new rules...")
      
      # Create a new metadata column for display, implementing the new rules
      # 1. Start with original cell labels (singler_labels)
      current_seurat_object$display_cell_type_custom <- as.character(current_seurat_object$singler_labels)
      
      # 2. First, assign "Cancer-Associated Fibroblast (CAF)" to all cells identified as CAF.
      #    This ensures CAFs are correctly labeled.
      current_seurat_object$display_cell_type_custom[current_seurat_object$is_CAF == "CAF"] <- "Cancer-Associated Fibroblast (CAF)"
      
      # 3. Now, assign "Tumour Cell" to all malignant cells.
      #    This will *override* any previous "CAF" label if the cell is both malignant and a CAF,
      #    thereby prioritizing the "Tumour Cell" designation for malignant cells.
      current_seurat_object$display_cell_type_custom[current_seurat_object$scm_malignancy == "Malignant"] <- "Tumour Cell"
      
      # Reorder levels for plotting: Tumour Cell, CAF, then others alphabetically
      all_unique_labels_display <- unique(current_seurat_object$display_cell_type_custom)
      special_labels_order <- c("Tumour Cell", "Cancer-Associated Fibroblast (CAF)") # Prioritize these
      remaining_labels_sorted <- sort(setdiff(all_unique_labels_display, special_labels_order))
      
      # Combine, ensuring only labels that actually exist are included and in desired order
      ordered_levels_display <- c(intersect(special_labels_order, all_unique_labels_display), remaining_labels_sorted)
      
      current_seurat_object$display_cell_type_custom <- factor(current_seurat_object$display_cell_type_custom, levels = ordered_levels_display)
      
      output_filename_integrated_custom <- paste0("UMAP_CellType_Integrated_CustomMarkers_", chosen_cancer_type, ".pdf")
      pdf(output_filename_integrated_custom, width = 10, height = 8)
      plot_integrated_custom <- DimPlot(current_seurat_object,
                                        reduction = "umap",
                                        group.by = "display_cell_type_custom",
                                        label = TRUE,
                                        repel = TRUE,
                                        pt.size = 0.5) +
        ggtitle(paste0("UMAP by Cell Type and Malignancy Status (Custom Markers: ", chosen_cancer_type, ")")) +
        guides(color = guide_legend(override.aes = list(size = 3))) # Increase legend point size
      print(plot_integrated_custom)
      dev.off()
      message(paste0("Integrated UMAP plot for custom marker malignancy saved to '", output_filename_integrated_custom, "'."))
      
      # Plot 3: UMAP colored by the continuous custom malignancy score
      output_filename_score <- "UMAP_Custom_Malignancy_Score.pdf"
      pdf(output_filename_score, width = 8, height = 7)
      plot_malig_score <- FeaturePlot(current_seurat_object,
                                      features = score_col_name,
                                      reduction = "umap",
                                      pt.size = 0.5,
                                      order = TRUE,
                                      cols = c("lightgrey", "darkred")) +
        ggtitle(paste0("UMAP by Custom Malignancy Score (", chosen_cancer_type, ")"))
      print(plot_malig_score)
      dev.off()
      message(paste0("UMAP plot showing continuous custom malignancy score saved to '", output_filename_score, "'."))
    }
    
  } else if (malignancy_method_choice == 2) {
    message("\n--- MALIGNANCY DETECTION (scMalignantFinder) ---")
    
    # Ensure reticulate is configured and the Python module is available
    if (!reticulate::py_module_available("scMalignantFinder")) {
      message("scMalignantFinder Python package not found. Please ensure the 'scmalignant' conda environment is activated and accessible.")
      stop("scMalignantFinder Python package missing. Cannot proceed with malignancy detection.")
    }
    
    scm_classifier_module <- reticulate::import("scMalignantFinder.classifier")
    
    message("Converting Seurat object to AnnData object for scMalignantFinder (this may take a moment for large datasets)...")
    sce_obj <- as.SingleCellExperiment(current_seurat_object)
    adata_obj <- scater::to_AnnData(sce_obj)
    rm(sce_obj) # Remove intermediate object to free memory
    gc() # Explicit garbage collection
    
    message("Initializing scMalignantFinder model and performing prediction (this can be computationally intensive)...")
    malignancy_model <- scm_classifier_module$scMalignantFinder(
      test_input = adata_obj,
      pretrain_dir = NULL,
      train_h5ad_path = NULL,
      feature_path = NULL,
      model_method = "LogisticRegression",
      norm_type = TRUE,
      n_thread = 7L # Ensure it's an integer
    )
    
    result_adata_py <- malignancy_model$predict()
    message("scMalignantFinder prediction complete. Retrieving results from Python...")
    
    current_seurat_object$scm_malignancy <- result_adata_py$obs$scMalignantFinder_prediction
    current_seurat_object$scm_scores <- result_adata_py$obs$malignancy_probability
    message("Malignancy status and scores successfully added to Seurat object metadata.")
    
    current_seurat_object$scm_malignancy <- factor(current_seurat_object$scm_malignancy, levels = c("Normal", "Malignant"))
    
    message("Generating UMAP with scMalignantFinder malignancy status...")
    
    output_filename_malig <- "UMAP_Malignancy_Status_scMalignantFinder.pdf"
    pdf(output_filename_malig, width = 8, height = 7)
    plot_malig_status <- DimPlot(current_seurat_object,
                                 reduction = "umap",
                                 group.by = "scm_malignancy",
                                 label = TRUE,
                                 repel = TRUE,
                                 pt.size = 0.5) +
      ggtitle("UMAP by scMalignantFinder Malignancy Status")
    print(plot_malig_status)
    dev.off()
    message(paste0("UMAP plot by scMalignantFinder inferred malignancy saved to '", output_filename_malig, "'."))
    
    message("Generating integrated UMAP plot for Cell Type and Malignancy Status (scMalignantFinder) with new rules...")
    
    # Create a new metadata column for display, implementing the new rules
    # 1. Start with original cell labels (singler_labels)
    current_seurat_object$display_cell_type_scm <- as.character(current_seurat_object$singler_labels)
    
    # 2. First, assign "Cancer-Associated Fibroblast (CAF)" to all cells identified as CAF.
    #    This ensures CAFs are correctly labeled.
    current_seurat_object$display_cell_type_scm[current_seurat_object$is_CAF == "CAF"] <- "Cancer-Associated Fibroblast (CAF)"
    
    # 3. Now, assign "Tumour Cell" to all malignant cells.
    #    This will *override* any previous "CAF" label if the cell is both malignant and a CAF,
    #    thereby prioritizing the "Tumour Cell" designation for malignant cells.
    current_seurat_object$display_cell_type_scm[current_seurat_object$scm_malignancy == "Malignant"] <- "Tumour Cell"
    
    # Reorder levels for plotting: Tumour Cell, CAF, then others alphabetically
    all_unique_labels_display_scm <- unique(current_seurat_object$display_cell_type_scm)
    special_labels_order_scm <- c("Tumour Cell", "Cancer-Associated Fibroblast (CAF)") # Prioritize these
    remaining_labels_sorted_scm <- sort(setdiff(all_unique_labels_display_scm, special_labels_order_scm))
    
    # Combine, ensuring only labels that actually exist are included and in desired order
    ordered_levels_display_scm <- c(intersect(special_labels_order_scm, all_unique_labels_display_scm), remaining_labels_sorted_scm)
    
    current_seurat_object$display_cell_type_scm <- factor(current_seurat_object$display_cell_type_scm, levels = ordered_levels_display_scm)
    
    output_filename_integrated_scm <- "UMAP_CellType_Integrated_scMalignantFinder.pdf"
    pdf(output_filename_integrated_scm, width = 10, height = 8)
    plot_integrated_scm <- DimPlot(current_seurat_object,
                                   reduction = "umap",
                                   group.by = "display_cell_type_scm",
                                   label = TRUE,
                                   repel = TRUE,
                                   pt.size = 0.5) +
      ggtitle("UMAP by Cell Type and Malignancy Status (scMalignantFinder)") +
      guides(color = guide_legend(override.aes = list(size = 3))) # Increase legend point size
    print(plot_integrated_scm)
    dev.off()
    message(paste0("Integrated UMAP plot for scMalignantFinder malignancy saved to '", output_filename_integrated_scm, "'."))
    
    output_filename_score <- "UMAP_scMalignantFinder_Score.pdf"
    pdf(output_filename_score, width = 8, height = 7)
    plot_malig_score <- FeaturePlot(current_seurat_object,
                                    features = "scm_scores",
                                    reduction = "umap",
                                    pt.size = 0.5,
                                    order = TRUE,
                                    cols = c("lightgrey", "darkred")) +
      ggtitle("UMAP by scMalignantFinder Malignancy Score")
    print(plot_malig_score)
    dev.off()
    message(paste0("UMAP plot showing continuous scMalignantFinder score saved to '", output_filename_score, "'."))
    message("All scMalignantFinder plots generated and saved.")
    message("scMalignantFinder detection step completed.")
  }
}
# --- 11. Interactive Gene Expression Analysis & Visualization ---
if (start_step <= 9) {
  message("\n--- INTERACTIVE GENE EXPRESSION ANALYSIS ---")
  
  target_gene <- toupper(readline(prompt = "Enter the gene name to analyze (e.g., GAPDH, ACTA2, FAP): "))
  
  if (!target_gene %in% rownames(current_seurat_object[["SCT"]])) {
    warning(paste0("Gene '", target_gene, "' not found. Skipping gene expression analysis."))
  } else {
    message(paste0("Analyzing expression of gene: ", target_gene))
    
    gene_expression_values <- GetAssayData(current_seurat_object, assay = "SCT", layer = "data")[target_gene, ]
    
    if (target_gene == "FAP") {
      message("\n--- FAP EXPRESSION STATISTICS ---")
      message(paste0("Minimum FAP expression: ", round(min(gene_expression_values), 4)))
      message(paste0("Maximum FAP expression: ", round(max(gene_expression_values), 4)))
      message(paste0("Mean FAP expression: ", round(mean(gene_expression_values), 4)))
      message(paste0("Median FAP expression: ", round(median(gene_expression_values), 4)))
      
      message("FAP expression quantiles:")
      print(quantile(gene_expression_values, c(0, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 1), na.rm = TRUE))
      
      pdf(paste0("Histogram_", target_gene, "_Expression.pdf"), width = 8, height = 6)
      hist(gene_expression_values,
           breaks = 50,
           xlab = paste0(target_gene, " Expression (SCT Normalized)"),
           main = paste0("Distribution of ", target_gene, " Expression"),
           col = "lightblue", border = "black")
      abline(v = mean(gene_expression_values), col = "red", lty = 2, lwd = 2)
      text(x = mean(gene_expression_values), y = max(hist(gene_expression_values, plot=FALSE)$counts) * 0.9,
           labels = paste0("Mean: ", round(mean(gene_expression_values), 2)), col = "red", pos = 4)
      dev.off()
      message(paste0("Histogram of FAP expression saved to 'Histogram_FAP_Expression.pdf'. Review this plot to help set your threshold."))
    }
    
    threshold_input <- readline(prompt = paste0("Enter a threshold value for '", target_gene, "' to consider cells positive (e.g., 0, 1, 2): "))
    threshold_value <- as.numeric(threshold_input)
    
    if (is.na(threshold_value)) {
      warning("Invalid threshold value. Skipping gene expression analysis.")
    } else {
      message(paste0("Using threshold: ", threshold_value, " for positive cells."))
      
      current_seurat_object$gene_expression_status <- ifelse(
        gene_expression_values > threshold_value,
        paste0(target_gene, "_Positive"),
        paste0(target_gene, "_Negative")
      )
      current_seurat_object$gene_expression_status <- factor(current_seurat_object$gene_expression_status,
                                                             levels = c(paste0(target_gene, "_Negative"), paste0(target_gene, "_Positive")))
      
      num_positive_cells <- sum(current_seurat_object$gene_expression_status == paste0(target_gene, "_Positive"))
      message(paste0("Number of ", target_gene, "_Positive cells: ", num_positive_cells))
      message(paste0("Number of ", target_gene, "_Negative cells: ", sum(current_seurat_object$gene_expression_status == paste0(target_gene, "_Negative"))))
      
      plot_title <- paste0("UMAP of ", target_gene, " Expression (Threshold > ", threshold_value, ")")
      output_filename <- paste0("UMAP_", target_gene, "_Positive_Negative_Plot.pdf")
      
      pdf(output_filename, width = 8, height = 7)
      print(DimPlot(current_seurat_object,
                    reduction = "umap",
                    group.by = "gene_expression_status",
                    cols = c("lightgrey", "red"),
                    label = FALSE,
                    repel = TRUE) +
              ggtitle(plot_title))
      dev.off()
      message(paste0("UMAP plot for '", target_gene, "' expression saved to '", output_filename, "'."))
      
      output_filename_feature <- paste0("UMAP_", target_gene, "_FeaturePlot.pdf")
      pdf(output_filename_feature, width = 8, height = 7)
      print(FeaturePlot(current_seurat_object,
                        features = target_gene,
                        reduction = "umap",
                        pt.size = 0.5,
                        order = TRUE) +
              ggtitle(paste0("UMAP of Continuous ", target_gene, " Expression")))
      dev.off()
      message(paste0("UMAP feature plot for '", target_gene, "' continuous expression saved to '", output_filename_feature, "'."))
    }
  }
}

# --- 12. Final Data Saving ---
if (start_step <= 10) {
  saveRDS(seurat_object, file = "seurat_object_final.rds")
} else {
}
# --- END OF PIPELINE ---