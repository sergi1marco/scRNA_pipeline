# A Comprehensive and Reproducible R Pipeline for Single-Cell RNA Sequencing Data Analysis

# Correcting the future.globals.maxSize option to use 1024^2 for MB to Bytes conversion.
# 8000 * 1024^2 means 8 GB. Adjust as needed for your system's RAM.
options(future.globals.maxSize = 8000 * 1024^2) # Set to 8 GB (adjust as needed)
options(scipen = 100) # Added as per inferCNV warning for hclust generation

# --- 0. Library Loading (Assumes all packages are pre-installed) ---
# All necessary packages must be installed and available in your R environment.
# If any package is missing, this script will stop with an error.

library(Seurat)
library(SeuratObject)
library(SingleCellExperiment)
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
library(infercnv) # Added for malignancy detection
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
library(org.Hs.eg.db) # For gene annotation, useful for enrichment
library(AnnotationDbi) # Explicitly load for select function
library(fgsea) # For GSEA
library(msigdbr) # For gene set collections (e.g., hallmark, oncogenic signatures)
library(circlize)
library(plotly)
library(reticulate) # Required for Python integration (e.g., if inferCNV uses it)
library(remotes) # Still needed if you explicitly call remotes::install_github somewhere else
library(EnhancedVolcano)
library(nebula)
library(presto) # For faster marker finding

message("\n--- All required libraries loaded successfully. ---")

# --- End of Library Loading ---


# --- 1. Cross-Platform Interactive Directory Selection & Pipeline Start Decision ---

# Function to select a directory interactively based on OS and R environment
select_data_directory <- function() {
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    message("Using RStudio API for directory selection...")
    data_path <- rstudioapi::selectDirectory(
      caption = "Select the working directory (or where your 10x Genomics data / saved Seurat objects are)",
      path = getwd() # Start in current working directory
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

# Set working directory to the selected data directory
data_dir <- select_data_directory()
setwd(data_dir)

message(paste("Working directory set to:", getwd()))

# Initialize seurat_object to NULL
seurat_object <- NULL
start_step <- 0 # Default to start from the very beginning (raw data)

# Check for existing Seurat objects
saved_seurat_files <- list.files(data_dir, pattern = "\\.rds$", full.names = TRUE, recursive = TRUE)
saved_seurat_files <- saved_seurat_files[grepl("seurat_object", basename(saved_seurat_files))] # Filter for files with "seurat_object" in name

if (length(saved_seurat_files) > 0) {
  message("\n--- Found existing Seurat object(s) in the working directory: ---")
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
      break # Exit loop, start_step remains 0
    } else if (user_choice_file > 0 && user_choice_file <= length(saved_seurat_files)) {
      seurat_object <- readRDS(saved_seurat_files[user_choice_file])
      message(paste("Loaded Seurat object from:", saved_seurat_files[user_choice_file]))
      
      # Ask user for desired starting step
      message("\n--- Pipeline Steps: ---")
      message("1: Quality Control and Filtering")
      message("2: Normalization and Scaling")
      message("3: Dimensionality Reduction (PCA and UMAP)")
      message("4: Visualization (Initial UMAPs, etc.)")
      message("5: Find Marker Genes for Clusters")
      message("6: Malignancy Detection (inferCNV Setup and Run)")
      message("7: Cell Type Annotation & CAF Identification (SingleR)")
      message("8: Differential Expression Analysis") # As per your list - functionality for this step would need to be added
      message("9: InferCNV Results Integration & Visualization")
      message("10: Interactive Gene Expression Analysis & Visualization") # New module
      message("11: Final Data Saving") # Renumbered
      message("-----------------------")
      
      user_choice_step <- readline(prompt = "From which step would you like to resume the pipeline (1-10)? ") # Updated max step
      start_step <- as.integer(user_choice_step)
      
      if (is.na(start_step) || start_step < 1 || start_step > 10) { # Updated max step
        warning("Invalid step chosen. Please enter a number between 1 and 10. Reloading Seurat object and re-asking.") # Updated max step
        seurat_object <- NULL # Reset to force re-load or new data
        next # Go back to start of while loop
      } else {
        step_names <- c("QC", "Normalization", "DimRed", "Visualization", "Markers", "InferCNV_Detection", "Annotation_CAF", "DE", "InferCNV_Integration_Vis", "Final_Save") # Updated step names
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

# Assign to a common variable for pipeline flow regardless of loading source
current_seurat_object <- seurat_object

# --- End of Directory Selection & Pipeline Start Decision ---


# --- 2. Data Acquisition and Initial Processing (Conditional) ---

# This section only runs if seurat_object is NULL (i.e., no saved object was loaded or user chose to start from raw)
if (is.null(seurat_object)) { # Corrected from is.is.null
  # Function to import 10x Genomics data (H5 or MTX)
  import_10x_data <- function(data_directory) {
    message(paste("Scanning directory for 10x Genomics data files:", data_directory))
    
    h5_files <- list.files(data_directory, pattern = "\\.h5$", full.names = TRUE, recursive = TRUE)
    mtx_dirs <- c()
    
    # Find directories that contain the three MTX files
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
    
    # Prioritize H5 files if available
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
      stop("No 10x Genomics H5 file or MTX data directory found in the selected folder or its subfolders. Please ensure your data is correctly placed.")
    }
    
    # Import data based on determined type
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
          stop("H5 file contains multiple data types, but 'Gene Expression' was not found and there are multiple options. Please adjust the script to select the correct assay from the H5 file (e.g., seurat_obj_raw[['VDJ_T']]) based on your data.")
        }
      }
    } else if (data_type == "mtx") {
      message(paste("Importing data from MTX folder:", selected_path))
      seurat_obj_raw <- Read10X(data.dir = selected_path)
    }
    
    # Create Seurat object
    seurat_object <- CreateSeuratObject(counts = seurat_obj_raw, project = "SingleCellAnalysis")
    message("Seurat object created successfully from 10x Genomics data.")
    
    return(seurat_object)
  }
  
  # Call the import_10x_data function with the selected data_dir.
  seurat_object <- import_10x_data(data_directory = data_dir)
}

# Assign to a common variable for pipeline flow regardless of loading source
current_seurat_object <- seurat_object

# --- End of Data Acquisition and Initial Processing ---


# --- 3. Quality Control and Filtering ---
if (start_step <= 1) {
  # Function to perform QC metrics calculation and initial filtering
  perform_qc_and_filter <- function(seurat_obj, min_umis = 500, max_umis = 100000,
                                    max_percent_mt = 10, min_features = 200, min_cells = 3) {
    message("Performing Quality Control and Filtering...")
    
    # Calculate mitochondrial percentage
    seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
    
    # Generate QC diagnostic plots
    pdf("QC_Plots.pdf", width = 10, height = 8)
    print(VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1))
    print(FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.mt"))
    print(FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA"))
    dev.off()
    message("Generating QC diagnostic plots (saved to 'QC_Plots.pdf').")
    message("Review 'QC_Plots.pdf' to refine filtering thresholds. Applying default/user-specified thresholds now.")
    
    # Apply filtering thresholds
    seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > min_features &
                           nCount_RNA > min_umis &
                           nCount_RNA < max_umis &
                           percent.mt < max_percent_mt)
    
    message(paste("QC and filtering complete. Remaining cells:", ncol(seurat_obj)))
    if (ncol(seurat_obj) == 0) {
      stop("No cells remaining after quality control and filtering. Please review 'QC_Plots.pdf' and adjust filtering thresholds (min_umis, max_umis, max_percent_mt, min_features).")
    }
    
    return(seurat_obj)
  }
  
  # Run QC and filtering
  current_seurat_object <- perform_qc_and_filter(current_seurat_object,
                                                 min_umis = 500,
                                                 max_umis = 100000,
                                                 max_percent_mt = 10,
                                                 min_features = 200,
                                                 min_cells = 3)
} else {
  message("Skipping Quality Control and Filtering (resuming from later step).")
}

# --- 4. Normalization and Scaling ---
if (start_step <= 2) {
  # Function for normalization and scaling
  perform_normalization_and_scaling <- function(seurat_obj) {
    message("Normalizing and Scaling data...")
    seurat_obj <- SCTransform(seurat_obj, vars.to.regress = "percent.mt", verbose = FALSE)
    message("Normalization and scaling complete.")
    return(seurat_obj)
  }
  # Run Normalization and Scaling
  current_seurat_object <- perform_normalization_and_scaling(current_seurat_object)
} else {
  message("Skipping Normalization and Scaling (resuming from later step).")
}


# --- 5. Dimensionality Reduction (PCA and UMAP) ---
if (start_step <= 3) {
  # Function for dimensionality reduction
  perform_dimensionality_reduction <- function(seurat_obj, dims = 1:30, resolution = 0.5) {
    message("Performing Dimensionality Reduction (PCA and UMAP)...")
    seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj))
    
    pdf("PCA_ElbowPlot.pdf", width = 8, height = 6)
    print(ElbowPlot(seurat_obj))
    dev.off()
    message("PCA complete. See 'PCA_ElbowPlot.pdf' to select appropriate number of dimensions.")
    
    seurat_obj <- FindNeighbors(seurat_obj, dims = dims)
    seurat_obj <- FindClusters(seurat_obj, resolution = resolution)
    
    # Using R-native UWOT for UMAP (default in Seurat v4+)
    seurat_obj <- RunUMAP(seurat_obj, dims = dims)
    message("UMAP and Clustering complete.")
    return(seurat_obj)
  }
  # Run Dimensionality Reduction
  current_seurat_object <- perform_dimensionality_reduction(current_seurat_object, dims = 1:30, resolution = 0.8)
} else {
  message("Skipping Dimensionality Reduction (resuming from later step).")
}

# --- 6. Visualization ---
if (start_step <= 4) {
  # Function for visualization
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
      # IMPORTANT: Using 'slot = "data"' as confirmed by args(Seurat::FeaturePlot) output
      print(FeaturePlot(seurat_obj, features = marker_genes_present, reduction = "umap", ncol = 3, slot = "data"))
      dev.off()
      message("UMAP feature plots for example marker genes saved to 'UMAP_FeaturePlots.pdf'.")
    } else {
      message("No example marker genes found in the dataset. Skipping FeaturePlot generation.")
    }
    message("All standard visualizations complete.")
    return(seurat_obj)
  }
  # Run Visualization
  current_seurat_object <- perform_visualization(current_seurat_object)
} else {
  message("Skipping Visualization (resuming from later step).")
}


# --- 7. Find Marker Genes for Clusters ---
if (start_step <= 5) {
  # Function to find cluster marker genes
  find_cluster_markers <- function(seurat_obj) {
    message("Finding marker genes for each cluster...")
    
    # The global option for future.globals.maxSize is set at the top of the script
    # options(future.globals.maxSize = 4000 * 1024^2) # No need to set it again here
    
    cluster_markers <- FindAllMarkers(seurat_obj,
                                      assay = "SCT",
                                      min.pct = 0.25,
                                      logfc.threshold = 0.25,
                                      only.pos = TRUE,
                                      verbose = FALSE)
    
    if (nrow(cluster_markers) == 0) {
      message("WARNING: No differentially expressed genes found for any cluster with the current thresholds.")
      message("Consider adjusting 'min.pct' or 'logfc.threshold' in FindAllMarkers if you expect marker genes.")
      # Create empty data frames to avoid errors later
      top10_markers <- data.frame(cluster=character(), gene=character(), avg_log2FC=numeric())
      write.csv(cluster_markers, "All_Cluster_Markers.csv", row.names = FALSE) # Will be an empty file
      write.csv(top10_markers, "Top10_Cluster_Markers.csv", row.names = FALSE) # Will be an empty file
      return(list(all_markers = cluster_markers, top10_markers = top10_markers))
    }
    
    # Proceed only if there are actual marker genes
    top10_markers <- cluster_markers %>%
      group_by(cluster) %>%
      slice_max(n = 10, order_by = avg_log2FC, with_ties = FALSE) # Added with_ties = FALSE for deterministic results
    
    write.csv(cluster_markers, "All_Cluster_Markers.csv", row.names = FALSE)
    write.csv(top10_markers, "Top10_Cluster_Markers.csv", row.names = FALSE)
    message("Marker gene finding complete. 'All_Cluster_Markers.csv' and 'Top10_Cluster_Markers.csv' saved.")
    return(list(all_markers = cluster_markers, top10_markers = top10_markers))
  }
  # Run Marker Gene Finding
  # Store markers_list globally or pass it to inferCNV
  markers_list <- find_cluster_markers(current_seurat_object)
  current_seurat_object@misc$all_markers <- markers_list$all_markers
  current_seurat_object@misc$top10_markers <- markers_list$top10_markers
  
} else {
  message("Skipping Find Marker Genes (resuming from later step).")
  # Attempt to load markers if resuming and they might be needed later (e.g., for auto-inferCNV)
  if (file.exists("All_Cluster_Markers.csv")) {
    current_seurat_object@misc$all_markers <- read.csv("All_Cluster_Markers.csv")
  }
  if (file.exists("Top10_Cluster_Markers.csv")) {
    current_seurat_object@misc$top10_markers <- read.csv("Top10_Cluster_Markers.csv")
  }
}

# --- 8. Malignancy Detection (inferCNV) ---
# This step is placed after normalization/scaling and before cell type annotation.
# It requires identifying a 'normal' cell population for reference.
if (start_step <= 6) {
  message("Starting Malignancy Detection using inferCNV...")
  
  # Add a flag to indicate if 'g' prefix is needed for cluster names in inferCNV
  # Based on the recurring warning "First group.by variable `seurat_clusters` starts with a number, appending `g`..."
  # We assume 'g' prefix is necessary for inferCNV's parsing and internal consistency.
  needs_g_prefix <- TRUE
  
  # Initialize normal_cell_groups
  normal_cell_groups <- NULL
  
  # >>> INSERT THIS LINE EXACTLY HERE <<<
  normal_cell_groups <- c("g2", "g4", "g5", "g8", "g17", "g21") # YOUR MANUALLY IDENTIFIED NORMAL GROUPS
  message(paste("Manually set normal cell groups for inferCNV:", paste(normal_cell_groups, collapse = ", ")))
  # >>> END INSERTION <<<
  
  # Strategy 1: Automatic selection based on SingleR labels (if available and run)
  if ("singler_labels" %in% names(current_seurat_object@meta.data)) {
    message("\n--- INFERCNV: Attempting automatic selection using SingleR annotations. ---")
    # Identify common normal cell types from SingleR labels
    normal_cell_type_keywords <- c("t cell", "b cell", "macrophage", "dendritic", "endothelial", "fibroblast", "immune", "myeloid", "lymphoid")
    
    # Get clusters associated with these normal keywords
    potential_normal_clusters <- character(0)
    for (keyword in normal_cell_type_keywords) {
      matching_labels <- unique(current_seurat_object$singler_labels[grepl(keyword, current_seurat_object$singler_labels, ignore.case = TRUE)])
      if (length(matching_labels) > 0) {
        # Map these labels back to Seurat cluster IDs
        matching_cluster_ids <- unique(current_seurat_object$seurat_clusters[current_seurat_object$singler_labels %in% matching_labels])
        potential_normal_clusters <- c(potential_normal_clusters, as.character(matching_cluster_ids))
      }
    }
    
    # Remove duplicates and ensure they are valid clusters
    potential_normal_clusters <- unique(potential_normal_clusters)
    
    # If g-prefix is needed, apply it to potential_normal_clusters nowC:\Users\Sergi Marco\OneDrive - University of Glasgow\SMARCO\AVACTA project\Databases\scRNAseq\GSE125449\CancerSCEM\inferCNV_output\BayesNetOutput.HMMi6.leiden.hmm_mode-subclusters
    if (needs_g_prefix && !any(grepl("^g", potential_normal_clusters))) {
      potential_normal_clusters <- paste0("g", potential_normal_clusters)
    }
    
    # Filter against actual Seurat levels (which should now also be g-prefixed if needs_g_prefix is TRUE)
    potential_normal_clusters <- potential_normal_clusters[potential_normal_clusters %in% levels(current_seurat_object$seurat_clusters)]
    
    if (length(potential_normal_clusters) > 0) {
      # Filter out clusters that were potentially classified as CAFs if CAF_Score1 exists
      if ("CAF_Score1" %in% colnames(current_seurat_object@meta.data)) {
        fibroblast_clusters_with_high_caf <- character(0)
        # To match "g" prefixed cluster IDs for filtering:
        seurat_clusters_char <- as.character(current_seurat_object$seurat_clusters)
        if (needs_g_prefix && !any(grepl("^g", seurat_clusters_char))) {
          seurat_clusters_char <- paste0("g", seurat_clusters_char)
        }
        
        for (clus_id in potential_normal_clusters) {
          if (any(grepl("fibroblast", current_seurat_object$singler_labels[seurat_clusters_char == clus_id], ignore.case = TRUE))) {
            median_caf_score_in_cluster <- median(current_seurat_object$CAF_Score1[seurat_clusters_char == clus_id], na.rm = TRUE)
            if (!is.na(median_caf_score_in_cluster) && median_caf_score_in_cluster > quantile(current_seurat_object$CAF_Score1, 0.75, na.rm = TRUE)) {
              fibroblast_clusters_with_high_caf <- c(fibroblast_clusters_with_high_caf, clus_id)
            }
          }
        }
        potential_normal_clusters <- setdiff(potential_normal_clusters, fibroblast_clusters_with_high_caf)
        if (length(fibroblast_clusters_with_high_caf) > 0) {
          message(paste("Automatically excluded fibroblast clusters due to high CAF score (based on SingleR labels):", paste(fibroblast_clusters_with_high_caf, collapse = ", ")))
        }
      }
      
      if (length(potential_normal_clusters) > 0) {
        normal_cell_groups <- potential_normal_clusters
        message(paste("Automatically selected normal clusters based on SingleR annotations:", paste(normal_cell_groups, collapse = ", ")))
      } else {
        message("Could not find suitable normal clusters automatically from SingleR annotations.")
      }
    } else {
      message("SingleR annotations did not yield clear normal cell populations.")
    }
  }
  
  
  # Strategy 2: Automatic selection based on immune markers (default fallback if SingleR not sufficient)
  # This section now runs automatically if normal_cell_groups is still NULL
  if (is.null(normal_cell_groups)) {
    message("\n--- INFERCNV: Attempting automatic selection using immune markers (as SingleR labels were not sufficient). ---")
    
    # Common immune cell markers (expand as needed for specificity)
    immune_markers <- c("CD3D", "CD19", "CD14", "PTPRC", "LYZ", "NKG7", "MS4A1")
    immune_markers_present <- immune_markers[immune_markers %in% rownames(current_seurat_object)]
    
    if (length(immune_markers_present) == 0) {
      # If even immune markers are missing, then we must fall back to manual input.
      warning("No common immune markers found in your dataset to attempt automatic normal cluster selection.")
      normal_clusters_input <- readline(prompt = paste0("AUTOMATIC SELECTION FAILED. Please manually enter normal cluster IDs separated by commas (e.g., '0,2,5'). Current clusters: ", paste(levels(current_seurat_object$seurat_clusters), collapse = ", "), ": "))
      normal_cell_groups <- unlist(strsplit(normal_clusters_input, ","))
      normal_cell_groups <- trimws(normal_cell_groups)
      
      # If g-prefix is needed, apply it to manually entered clusters
      if (needs_g_prefix && !any(grepl("^g", normal_cell_groups))) {
        normal_cell_groups <- paste0("g", normal_cell_groups)
      }
      if (length(normal_cell_groups) == 0 || any(!normal_cell_groups %in% levels(current_seurat_object$seurat_clusters))) {
        stop("Invalid or no normal cluster IDs provided. Cannot proceed with inferCNV.")
      }
      message(paste("Using manually provided normal cell clusters for inferCNV:", paste(normal_cell_groups, collapse = ", ")))
      
    } else {
      # Find clusters with highest average expression of immune markers
      # AverageExpression will return column names as the actual cluster names (e.g., "g0", "g1" if Seurat changed them)
      avg_immune_expr <- AverageExpression(current_seurat_object, features = immune_markers_present, assays = "SCT", group.by = "seurat_clusters")$SCT
      
      if (is.null(avg_immune_expr) || nrow(avg_immune_expr) == 0) {
        warning("Could not calculate average immune marker expression for automatic selection.")
        normal_clusters_input <- readline(prompt = paste0("AUTOMATIC SELECTION FAILED. Please manually enter normal cluster IDs separated by commas (e.g., '0,2,5'). Current clusters: ", paste(levels(current_seurat_object$seurat_clusters), collapse = ", "), ": "))
        normal_cell_groups <- unlist(strsplit(normal_clusters_input, ","))
        normal_cell_groups <- trimws(normal_cell_groups)
        # If g-prefix is needed, apply it to manually entered clusters
        if (needs_g_prefix && !any(grepl("^g", normal_cell_groups))) {
          normal_cell_groups <- paste0("g", normal_cell_groups)
        }
        if (length(normal_cell_groups) == 0 || any(!normal_cell_groups %in% levels(current_seurat_object$seurat_clusters))) {
          stop("Invalid or no normal cluster IDs provided. Cannot proceed with inferCNV.")
        }
        message(paste("Using manually provided normal cell clusters for inferCNV:", paste(normal_cell_groups, collapse = ", ")))
      } else {
        # Sum average expression for all immune markers per cluster
        cluster_immune_scores <- colSums(avg_immune_expr)
        
        # Identify top cluster(s) with highest immune score. Consider top 1-3.
        # Take clusters whose immune score is above the 75th percentile of all cluster immune scores.
        threshold_score <- quantile(cluster_immune_scores, 0.75)
        candidate_normal_clusters <- names(cluster_immune_scores[cluster_immune_scores >= threshold_score])
        
        # Ensure at least one cluster is selected if possible, or just the very highest one
        if (length(candidate_normal_clusters) == 0 && length(cluster_immune_scores) > 0) {
          candidate_normal_clusters <- names(which.max(cluster_immune_scores))
        }
        
        if (length(candidate_normal_clusters) > 0) {
          normal_cell_groups <- candidate_normal_clusters
          message(paste("Automatically selected normal clusters based on immune markers:", paste(normal_cell_groups, collapse = ", ")))
          message("WARNING: This is an automated best-guess. Please review the 'UMAP_Clusters.pdf' and marker gene lists to confirm these clusters are truly normal before proceeding with inferCNV interpretation.")
        } else {
          # If still no clusters identified, fall back to explicit manual prompt
          normal_clusters_input <- readline(prompt = paste0("AUTOMATIC SELECTION FAILED. Please manually enter normal cluster IDs separated by commas (e.g., '0,2,5'). Current clusters: ", paste(levels(current_seurat_object$seurat_clusters), collapse = ", "), ": "))
          normal_cell_groups <- unlist(strsplit(normal_clusters_input, ","))
          normal_cell_groups <- trimws(normal_cell_groups)
          # If g-prefix is needed, apply it to manually entered clusters
          if (needs_g_prefix && !any(grepl("^g", normal_cell_groups))) {
            normal_cell_groups <- paste0("g", normal_cell_groups)
          }
          if (length(normal_cell_groups) == 0 || any(!normal_cell_groups %in% levels(current_seurat_object$seurat_clusters))) {
            stop("Invalid or no normal cluster IDs provided. Cannot proceed with inferCNV.")
          }
          message(paste("Using manually provided normal cell clusters for inferCNV:", paste(normal_cell_groups, collapse = ", ")))
        }
      }
    }
  }
  
  
  # Create annotations file for inferCNV
  # This maps cell barcodes to their cluster IDs for inferCNV's grouping.
  annotations_df <- data.frame(
    cell_id = colnames(current_seurat_object),
    cell_group = as.character(current_seurat_object$seurat_clusters)
  )
  
  # Apply 'g' prefix to annotations_df$cell_group if needed and not already present
  if (needs_g_prefix && !any(grepl("^g", annotations_df$cell_group))) {
    annotations_df$cell_group <- paste0("g", annotations_df$cell_group)
    message("Adjusted cell_group names in inferCNV annotations with 'g' prefix.")
  }
  
  # Save annotations to a temporary file
  infercnv_annotations_file <- "infercnv_cell_annotations.txt"
  write.table(annotations_df, infercnv_annotations_file, sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
  
  # Get gene order file (required by inferCNV)
  gene_symbols <- rownames(current_seurat_object)
  gene_order_df <- tryCatch({
    # Query for CHR and CHRLOC (chromosome location, often the start or a single point)
    tmp_df <- AnnotationDbi::select(org.Hs.eg.db, keys = gene_symbols, columns = c("CHR", "CHRLOC"), keytype = "SYMBOL")
    
    # Filter out genes without CHR or CHRLOC
    tmp_df <- tmp_df[!is.na(tmp_df$CHR) & !is.na(tmp_df$CHRLOC), ]
    
    # Create required columns for inferCNV: SYMBOL, CHR, TXSTART, TXEND
    # Use CHRLOC as approximate start/end for inferCNV sorting.
    gene_order_res <- data.frame(
      SYMBOL = tmp_df$SYMBOL,
      CHR = tmp_df$CHR,
      TXSTART = as.integer(tmp_df$CHRLOC), # Use CHRLOC as start
      TXEND = as.integer(tmp_df$CHRLOC) + 1 # Use CHRLOC + 1 as a dummy end for minimal range
    ) %>% distinct(SYMBOL, .keep_all = TRUE) # Ensure unique gene symbols
    
    # Match the order of genes in the Seurat object
    gene_order_res <- gene_order_res[match(rownames(current_seurat_object), gene_order_res$SYMBOL), ]
    na.omit(gene_order_res) # Remove genes not found in the annotation or with NA values
    
  }, error = function(e) {
    warning("Failed to retrieve gene order from org.Hs.eg.db using CHR and CHRLOC. Error: ", e$message, " Attempting fallback gene order creation.")
    NULL # Return NULL to trigger fallback
  })
  
  # Fallback if gene_order_df is NULL or empty after initial attempt
  if (is.null(gene_order_df) || nrow(gene_order_df) == 0 || !all(c("SYMBOL", "CHR", "TXSTART", "TXEND") %in% colnames(gene_order_df))) {
    warning("Final fallback: Creating a minimal gene order file. InferCNV might use its internal ordering if this is insufficient.")
    gene_order_df <- data.frame(
      SYMBOL = rownames(current_seurat_object),
      CHR = "Unknown", # Placeholder for chromosome
      TXSTART = 1,     # Placeholder for start
      TXEND = 2        # Placeholder for end
    )
    gene_order_df <- gene_order_df[!duplicated(gene_order_df$SYMBOL), ] # Ensure unique symbols
  }
  
  infercnv_gene_order_file <- "infercnv_gene_order.txt"
  write.table(gene_order_df, infercnv_gene_order_file, sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
  
  
  # Create the inferCNV object
  infercnv_obj <- CreateInfercnvObject(
    raw_counts_matrix = GetAssayData(current_seurat_object, assay = "RNA", slot = "counts"), # Use raw counts
    annotations_file = infercnv_annotations_file,
    gene_order_file = infercnv_gene_order_file,
    delim = "\t",
    ref_group_names = normal_cell_groups # CORRECT: ref_group_names is for CreateInfercnvObject
  )
  
  # Run inferCNV
  message("Running inferCNV... This may take a long time depending on data size.")
  infercnv_obj <- infercnv::run(infercnv_obj,
                                cutoff = 0.1, # Min avg read count to include gene
                                out_dir = "inferCNV_output", # Output directory
                                cluster_by_groups = TRUE, # Clusters cells within defined groups
                                denoise = TRUE, # Denoises the expression data
                                HMM = TRUE, # Apply Hidden Markov Model for CNV prediction
                                num_threads = parallel::detectCores() - 1, # Use most available cores
                                plot_steps = FALSE, # Set to TRUE for debugging/visualization of intermediate steps
                                HMM_type = "i6" # Corrected from "hspike" to "i6"
  )
  message("inferCNV analysis complete. Results saved in 'inferCNV_output' directory.")
  
  # Optional: Integrate inferCNV results back into Seurat object
  message("Please review the 'inferCNV_output' directory for CNV heatmaps and classifications.")
  message("You may use inferCNV results to annotate malignant cells in a later step or refine clustering.")
  
} else {
  message("Skipping Malignancy Detection (inferCNV) (resuming from later step).")
}
# --- 9. Cell Type Annotation & CAF Identification (using SingleR and custom markers) ---
if (start_step <= 7) { # This is now step 7
  
  # Function for cell type annotation and CAF identification
  perform_cell_type_annotation_and_caf <- function(seurat_obj) {
    message("Performing general cell type annotation using SingleR...")
    
    # --- Reference Selection for SingleR ---
    message("\n--- SINGLE R: REFERENCE SELECTION ---")
    message("Choose a reference dataset for cell type annotation:")
    message("1: HumanPrimaryCellAtlasData (General human cell types)")
    message("2: NovershternHematopoieticStemCellAtlas (Hematopoietic cells)")
    message("3: BlueprintEncodeData (Bulk RNA-seq for immune/stromal cells)")
    message("4: MonacoImmuneData (Comprehensive immune cell types)")
    message("5: HPCA + BlueprintEncode (Combined, good for broader coverage)")
    
    ref_choice <- readline(prompt = "Enter your choice (1-5): ")
    ref_choice <- as.integer(ref_choice)
    
    ref <- NULL
    tryCatch({
      if (ref_choice == 1) {
        ref <- celldex::HumanPrimaryCellAtlasData()
        message("HumanPrimaryCellAtlasData reference loaded.")
        # --- New CSV-based Tissue/Cell Type Filtering Logic ---
        csv_file_path <- "HPCA_tissue_cell_types.csv"
        
        if (file.exists(csv_file_path)) {
          message(paste0("Reading cell type mapping from ", csv_file_path, "..."))
          tissue_cell_types_map <- data.table::fread(csv_file_path) # Uses data.table::fread for efficiency
          
          # Validate CSV columns
          if (!("Tissue" %in% colnames(tissue_cell_types_map)) || !("CellTypes" %in% colnames(tissue_cell_types_map))) {
            warning(paste0("CSV file '", csv_file_path, "' must contain 'Tissue' and 'CellTypes' columns. Skipping tissue-specific filtering and using full HPCA data."))
            final_cell_types_to_filter <- unique(ref$label.main) # Fallback: No filtering
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
              
              # Retrieve cell types for the chosen tissue
              cell_types_string <- tissue_cell_types_map[Tissue == chosen_tissue]$CellTypes
              retrieved_cell_types <- unlist(strsplit(cell_types_string, split = ",", fixed = TRUE))
              retrieved_cell_types <- trimws(retrieved_cell_types) # Remove any leading/trailing spaces
              
              # Add PBMCs to the list if not already present
              if (!("PBMCs" %in% retrieved_cell_types)) {
                retrieved_cell_types <- c(retrieved_cell_types, "PBMCs")
                message("Added 'PBMCs' to the cell type list for filtering.")
              }
              
              # Filter the HPCA reference based on the combined list of cell types
              # Only use cell types that actually exist in HPCA's label.main
              final_cell_types_to_filter <- retrieved_cell_types[retrieved_cell_types %in% unique(ref$label.main)]
              
              if (length(final_cell_types_to_filter) > 0) {
                ref <- ref[, ref$label.main %in% final_cell_types_to_filter]
                message(paste0("HumanPrimaryCellAtlasData filtered to include: ", paste(final_cell_types_to_filter, collapse = ", ")))
              } else {
                warning("None of the cell types retrieved from CSV (including PBMCs) were found in HPCA's label.main. Using full HPCA data for annotation.")
                final_cell_types_to_filter <- unique(ref$label.main) # Fallback: No filtering
              }
              
            } else {
              warning("Invalid tissue selection. Using full HumanPrimaryCellAtlasData for annotation.")
              final_cell_types_to_filter <- unique(ref$label.main) # Fallback: No filtering
            }
          }
        } else {
          warning(paste0("CSV file '", csv_file_path, "' not found in the working directory. Using full HumanPrimaryCellAtlasData for annotation."))
          final_cell_types_to_filter <- unique(ref$label.main) # Fallback: No filtering
        }
        # --- End New CSV-based Tissue/Cell Type Filtering Logic ---
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
        ref <- c(ref_hpca, ref_bp) # Combine references if needed
        message("HumanPrimaryCellAtlasData and BlueprintEncodeData references combined.")
      } else {
        stop("Invalid reference choice. Using HumanPrimaryCellAtlasData as default.")
        ref <- celldex::HumanPrimaryCellAtlasData() # Fallback
      }
    }, error = function(e) {
      stop("Failed to load selected celldex reference. Please check your internet connection or celldex installation. Error: ", e$message)
    })
    
    sce <- as.SingleCellExperiment(seurat_obj)
    predictions <- SingleR(test = sce,
                           ref = ref,
                           labels = ref$label.main)
    
    seurat_obj$singler_labels <- predictions$pruned.labels
    message("Initial cell type annotation complete.")
    
    # Initialize combined_cell_types before CAF identification.
    # This ensures it always exists for plotting, even if CAF analysis is skipped.
    seurat_obj$combined_cell_types <- seurat_obj$singler_labels
    message("Initialized 'combined_cell_types' with SingleR labels.")
    
    # --- CAF Identification ---
    message("Attempting to identify Cancer-Associated Fibroblasts (CAFs)...")
    
    # Identify clusters predominantly composed of fibroblast/stromal cells based on SingleR annotation
    # Count occurrences of SingleR labels per cluster
    label_counts_per_cluster <- as.data.frame(table(seurat_obj$seurat_clusters, seurat_obj$singler_labels))
    colnames(label_counts_per_cluster) <- c("cluster", "label", "count")
    
    # Find the most frequent label for each cluster
    most_frequent_labels <- label_counts_per_cluster %>%
      group_by(cluster) %>%
      slice_max(order_by = count, n = 1) %>%
      ungroup()
    
    # Identify clusters where the most frequent label is fibroblast/stromal/mesenchymal
    fibroblast_clusters <- most_frequent_labels$cluster[grepl("fibroblast|stromal|mesenchymal", most_frequent_labels$label, ignore.case = TRUE)]
    fibroblast_clusters <- as.character(fibroblast_clusters) # Ensure character vector
    
    if (length(fibroblast_clusters) > 0) {
      message(paste("Potential fibroblast/stromal clusters based on initial annotation:", paste(fibroblast_clusters, collapse = ", ")))
      
      # Define common CAF markers and normal fibroblast markers
      # These lists can be expanded based on specific cancer types.
      caf_markers <- c("ACTA2", "FAP", "PDGFRA", "PDGFRB", "COL1A1", "COL1A2", "FN1", "VIM", "TWIST1", "SNAI1") # General CAF markers
      normal_fib_markers <- c("DCN", "LUM", "CD34", "THY1") # General normal fibroblast markers
      
      # Add module score for CAF signature to all cells
      caf_markers_present <- caf_markers[caf_markers %in% rownames(seurat_obj)]
      if (length(caf_markers_present) > 0) {
        seurat_obj <- AddModuleScore(seurat_obj, features = list(CAF_Signature = caf_markers_present), name = "CAF_Score")
        message("Added CAF signature score to all cells.")
        
        # Initialize is_CAF for all cells as "Non_CAF"
        seurat_obj$is_CAF <- "Non_CAF"
        
        # Proceed with CAF classification only if there are fibroblast clusters and CAF markers
        if (length(fibroblast_clusters) > 0 && length(caf_markers_present) > 0) {
          # Identify cells belonging to the fibroblast clusters
          cells_in_fib_clusters <- colnames(seurat_obj)[seurat_obj$seurat_clusters %in% fibroblast_clusters]
          
          if (length(cells_in_fib_clusters) > 0) {
            # Determine a dynamic threshold based on the CAF_Score1 distribution *within* identified fibroblast cells
            # For instance, cells in the top 75th percentile of CAF_Score1 among fibroblasts
            caf_score_for_fibroblasts <- seurat_obj$CAF_Score1[cells_in_fib_clusters]
            caf_threshold_for_classification_fib <- quantile(caf_score_for_fibroblasts, 0.75, na.rm = TRUE)
            
            # Classify cells as CAF if they are in a fibroblast cluster AND their CAF_Score1 is above the threshold
            seurat_obj$is_CAF[seurat_obj$seurat_clusters %in% fibroblast_clusters & seurat_obj$CAF_Score1 > caf_threshold_for_classification_fib] <- "CAF"
            
            message(paste("Cells classified as CAF (within fibroblast/stromal clusters) based on CAF_Score1 >", round(caf_threshold_for_classification_fib, 2), ":", sum(seurat_obj$is_CAF == "CAF")))
          } else {
            message("No cells found in identified fibroblast clusters to perform detailed CAF classification.")
          }
        }
        
        # Ensure 'is_CAF' is a factor with preferred order for plotting
        seurat_obj$is_CAF <- factor(seurat_obj$is_CAF, levels = c("Non_CAF", "CAF"))
        
        # Create a combined cell type label for plotting, integrating CAF status
        seurat_obj$combined_cell_types <- seurat_obj$singler_labels
        # If any cells were classified as CAF, update their labels
        if ("CAF" %in% levels(seurat_obj$is_CAF) && any(seurat_obj$is_CAF == "CAF")) {
          seurat_obj$combined_cell_types[seurat_obj$is_CAF == "CAF"] <- "Cancer-Associated Fibroblast (CAF)" # Assign a unique CAF label
        }
        seurat_obj$combined_cell_types <- factor(seurat_obj$combined_cell_types)
        message("Updated 'combined_cell_types' with integrated CAF status.")
        
        # Visualize CAF score on UMAP
        if ("CAF_Score1" %in% colnames(seurat_obj@meta.data)) {
          pdf("UMAP_CAF_Score.pdf", width = 8, height = 7)
          print(FeaturePlot(seurat_obj, features = "CAF_Score1", reduction = "umap", pt.size = 0.5, slot = "data", order = TRUE))
          dev.off()
          message("UMAP plot of CAF signature score saved to 'UMAP_CAF_Score.pdf'.")
        }
        
        message("Review 'UMAP_CAF_Score.pdf' and individual marker expression to confirm CAF populations.")
        
      } else {
        message("No CAF markers present in the Seurat object. Skipping detailed CAF scoring and classification.")
      }
      
    } else {
      message("No clusters identified as fibroblast/stromal based on initial SingleR annotation. Skipping detailed CAF analysis.")
    }
    
    # Plot 1: UMAP by Automated Cell Type Annotation (UMAP_CellTypeAnnotation.pdf)
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
  
  # Run Cell Type Annotation & CAF Identification (this call is inside the 'if (start_step <= 7)' block)
  current_seurat_object <- perform_cell_type_annotation_and_caf(current_seurat_object)
  
} else { # This brace correctly closes the 'if (start_step <= 7)' block, and the 'else' immediately follows it.
  message("Skipping Cell Type Annotation & CAF Identification (resuming from later step).")
}
# --- 10. Interactive Gene Expression Analysis & Visualization ---
# This module allows you to query a specific gene, define a positivity threshold,
# and visualize cells expressing that gene above the threshold on a UMAP plot.
if (start_step <= 10) { # Adjust step number as needed
  message("\n--- INTERACTIVE GENE EXPRESSION ANALYSIS ---")
  
  # 1. Interactively ask for a gene name
  target_gene <- toupper(readline(prompt = "Enter the gene name to analyze (e.g., GAPDH, ACTA2, FAP): "))
  
  # Check if the gene exists in the SCT assay's features
  if (!target_gene %in% rownames(current_seurat_object[["SCT"]])) {
    warning(paste0("Gene '", target_gene, "' not found in the SCT assay features. Skipping gene expression analysis."))
  } else {
    message(paste0("Analyzing expression of gene: ", target_gene))
    
    # Get gene expression values from the SCT assay's 'data' layer
    gene_expression_values <- GetAssayData(current_seurat_object, assay = "SCT", layer = "data")[target_gene, ]
    
    # *** ADDED SECTION: FAP EXPRESSION STATISTICS ***
    if (target_gene == "FAP") { # Check specifically for FAP
      message("\n--- FAP EXPRESSION STATISTICS ---")
      message(paste0("Minimum FAP expression: ", round(min(gene_expression_values), 4)))
      message(paste0("Maximum FAP expression: ", round(max(gene_expression_values), 4)))
      message(paste0("Mean FAP expression: ", round(mean(gene_expression_values), 4)))
      message(paste0("Median FAP expression: ", round(median(gene_expression_values), 4)))
      
      # Quantiles can be very helpful
      message("FAP expression quantiles:")
      print(quantile(gene_expression_values, c(0, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 1), na.rm = TRUE))
      
      # For a visual representation, a histogram is also very useful.
      pdf(paste0("Histogram_", target_gene, "_Expression.pdf"), width = 8, height = 6)
      hist(gene_expression_values,
           breaks = 50, # Adjust number of breaks for better visualization
           xlab = paste0(target_gene, " Expression (SCT Normalized)"),
           main = paste0("Distribution of ", target_gene, " Expression"),
           col = "lightblue", border = "black")
      abline(v = mean(gene_expression_values), col = "red", lty = 2, lwd = 2)
      text(x = mean(gene_expression_values), y = max(hist(gene_expression_values, plot=FALSE)$counts) * 0.9,
           labels = paste0("Mean: ", round(mean(gene_expression_values), 2)), col = "red", pos = 4)
      dev.off()
      message(paste0("Histogram of FAP expression saved to 'Histogram_FAP_Expression.pdf'. Review this plot to help set your threshold."))
    }
    # *** END ADDED SECTION ***
    
    # 2. Ask for a threshold value to consider cells positive
    threshold_input <- readline(prompt = paste0("Enter a threshold value for '", target_gene, "' to consider cells positive (e.g., 0, 1, 2): "))
    threshold_value <- as.numeric(threshold_input)
    
    if (is.na(threshold_value)) {
      warning("Invalid threshold value entered. Please enter a number. Skipping gene expression analysis.")
    } else {
      message(paste0("Using threshold: ", threshold_value, " for positive cells."))
      
      # 3. Generate positive vs. negative population
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
      
      # 4. Plot the positive vs negative population on a UMAP plot
      plot_title <- paste0("UMAP of ", target_gene, " Expression (Threshold > ", threshold_value, ")")
      output_filename <- paste0("UMAP_", target_gene, "_Positive_Negative_Plot.pdf")
      
      pdf(output_filename, width = 8, height = 7)
      print(DimPlot(current_seurat_object,
                    reduction = "umap",
                    group.by = "gene_expression_status",
                    cols = c("lightgrey", "red"), # Negative cells in grey, positive in red
                    label = FALSE,
                    repel = TRUE) +
              ggtitle(plot_title))
      dev.off()
      message(paste0("UMAP plot for '", target_gene, "' expression saved to '", output_filename, "'."))
      
      # Optional: Also plot the continuous expression for context
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
# --- 11. Differential Expression Analysis (Optional, Example per cluster) ---
if (start_step <= 8) { # Now step 8
  # Function to perform differential expression (e.g., between two clusters)
  perform_differential_expression <- function(seurat_obj, ident1, ident2, assay = "SCT") {
    message(paste("Performing differential expression analysis between cluster", ident1, "and cluster", ident2, "..."))
    
    if (!as.character(ident1) %in% levels(seurat_obj) || !as.character(ident2) %in% levels(seurat_obj)) {
      warning(paste("One or both specified identities (", ident1, ",", ident2, ") not found in Seurat object's current active identity (seurat_clusters). Skipping DE analysis.", sep=""))
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
  
  # Example usage for DE analysis (uncomment and modify as needed)
  # Replace '0' and '1' with actual cluster numbers found in your data (e.g., from UMAP_Clusters.pdf)
  # de_results_cluster_0_vs_1 <- perform_differential_expression(current_seurat_object, ident1 = 0, ident2 = 1, assay = "SCT")
} else {
  message("Skipping Differential Expression Analysis (resuming from later step).")
}

# Automated inferCNV Integration and Visualization
# Place this code after the "8. Differential Expression Analysis" section
# and before the "10. Final Data Saving" section in your script.

# --- Prompt to choose inferCNV object file ---
message("\nPlease select your saved inferCNV object file (e.g., 'final.infercnv_obj.rds').")
message("This will open a file browser window.")

# Use file.choose() to open a file selection dialog
selected_infercnv_obj_path <- tryCatch(
  {
    file.choose()
  },
  error = function(e) {
    NULL # Return NULL if an error occurs (e.g., user cancels)
  }
)

loaded_infercnv_obj <- NULL
if (is.null(selected_infercnv_obj_path) || !file.exists(selected_infercnv_obj_path)) {
  stop("ERROR: No inferCNV object file selected or file does not exist. Please run the script again and select the 'final.infercnv_obj.rds' file.")
} else {
  loaded_infercnv_obj <- readRDS(selected_infercnv_obj_path)
  message(paste("Loaded inferCNV object from:", selected_infercnv_obj_path))
}

# --- Automated Malignancy Classification using CNV Deviation Score ---
# Extract the inferCNV expression matrix (transformed and denoised)
cnv_matrix <- loaded_infercnv_obj@expr.infercnv

# Calculate a 'CNV Deviation Score' for each cell (sum of absolute deviations from baseline 1.0)
cell_cnv_scores <- colSums(abs(cnv_matrix - 1), na.rm = TRUE)

# Add this score to the Seurat object's metadata
current_seurat_object$infercnv_score <- cell_cnv_scores[colnames(current_seurat_object)] # Ensure order matches

# Define 'Normal' cell groups based on what was used in inferCNV run.
# If `normal_cell_groups` is not defined (e.g., when resuming from a later step),
# you might need to manually set it here based on your inferCNV reference clusters.
# Example: normal_cell_groups_for_exclusion <- c("g0", "g1", "g2")
normal_cell_groups_for_exclusion <- if (exists("normal_cell_groups") && !is.null(normal_cell_groups)) {
  temp_groups <- normal_cell_groups
  if (any(grepl("^g", levels(current_seurat_object$seurat_clusters))) && !any(grepl("^g", temp_groups))) {
    temp_groups <- paste0("g", temp_groups)
  }
  temp_groups
} else {
  warning("`normal_cell_groups` variable not found. Assuming no specific normal groups are explicitly excluded from automated malignancy classification. Review 'UMAP_InferCNV_Score.pdf' carefully.")
  character(0) # Exclude no groups if variable not found
}

# Classify cells: "Malignant" if CNV score is in the top X percentile AND not in normal groups
# Adjust the quantile (e.g., 0.9 for top 10%, 0.75 for top 25%) based on your data and expected malignancy rate.
cnv_score_threshold <- quantile(current_seurat_object$infercnv_score, 0.9, na.rm = TRUE) # Top 10% by default

current_seurat_object$infercnv_malignancy <- ifelse(
  current_seurat_object$infercnv_score > cnv_score_threshold &
    !(as.character(current_seurat_object$seurat_clusters) %in% normal_cell_groups_for_exclusion), # Exclude known normal cells
  "Malignant",
  "Non-Malignant"
)
message(paste("Automated classification: Cells with CNV score >", round(cnv_score_threshold, 2), "(top 10%) and not in specified normal clusters are classified as 'Malignant'."))


# Ensure the factor levels are ordered for consistent plotting
current_seurat_object$infercnv_malignancy <- factor(
  current_seurat_object$infercnv_malignancy,
  levels = c("Non-Malignant", "Malignant")
)

# --- Visualize Malignant Cells on UMAP ---
message("Generating UMAP visualizations with inferred malignancy status...")

# Plot 1: UMAP colored by Malignancy status
pdf("UMAP_Malignancy_Status.pdf", width = 8, height = 7)
print(DimPlot(current_seurat_object,
              reduction = "umap",
              group.by = "infercnv_malignancy",
              label = FALSE,
              pt.size = 0.5,
              cols = c("Non-Malignant" = "grey", "Malignant" = "red") # Customize colors
) + ggtitle("UMAP by Inferred Malignancy Status (Automated)"))
dev.off()
message("UMAP plot showing inferred malignancy saved to 'UMAP_Malignancy_Status.pdf'.")


# Plot 2: UMAP combining Cell Type Annotation and Malignancy Status (split view)
pdf("UMAP_CellType_SplitByMalignancy.pdf", width = 16, height = 8)
print(DimPlot(current_seurat_object,
              reduction = "umap",
              group.by = "singler_labels",
              label = TRUE,
              repel = TRUE,
              pt.size = 0.5,
              split.by = "infercnv_malignancy"
) + ggtitle("UMAP by Cell Type and Malignancy Status (Automated)") + NoLegend())
dev.off()
message("UMAP plot by cell type, split by malignancy, saved to 'UMAP_CellType_SplitByMalignancy.pdf'.")

# Plot 3: UMAP colored by the continuous inferCNV score
pdf("UMAP_InferCNV_Score.pdf", width = 8, height = 7)
print(FeaturePlot(current_seurat_object,
                  features = "infercnv_score",
                  reduction = "umap",
                  pt.size = 0.5,
                  order = TRUE,
                  cols = c("lightgrey", "darkred")) + # Color scale for continuous score
        ggtitle("UMAP by InferCNV Deviation Score"))
dev.off()
message("UMAP plot showing continuous inferCNV deviation score saved to 'UMAP_InferCNV_Score.pdf'.")

# --- 12. Final Data Saving ---

message("Saving final processed Seurat object...")
saveRDS(current_seurat_object, "final_seurat_object.rds")
message("Final Seurat object saved as 'final_seurat_object.rds'.")

# --- End of Pipeline ---

run_sc_pipeline()