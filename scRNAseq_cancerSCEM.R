# =============================================================================
# COMPREHENSIVE SINGLE-CELL RNA-SEQ ANALYSIS PIPELINE
# =============================================================================
# Following Luecken & Theis (2019) and scRNA-seq best practices
# Adapted to introduce changes according to CancerSCEM pipeline
# =============================================================================

# GLOBAL SETTINGS AND CONFIGURATION
options(
  future.globals.maxSize = 8000 * 1024^2, # 8GB
  Seurat.object.assay.version = "v5"
)
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

# =============================================================================
# PACKAGE INSTALLATION AND MANAGEMENT
# =============================================================================
install_required_packages <- function() {
  cat("Installing required packages...\n")
  core_packages <- c(
    "Seurat", "SeuratObject", "dplyr", "ggplot2", "patchwork",
    "Matrix", "R.utils", "reticulate", "hdf5r", "devtools", "openxlsx", "readxl",
    "tibble"
  )
  bioc_packages <- c(
    "scater", "scran", "BiocParallel", "AnnotationDbi", "org.Hs.eg.db",
    "SingleR", "celldex", "scDblFinder", "SingleCellExperiment"
  )
  specialized_packages <- c(
    "harmony", "fastMNN", "batchelor", "slingshot", "tradeSeq",
    "monocle3", "CellChat", "nichenetr", "velocyto.R", "SeuratWrappers"
  )
  install_if_missing <- function(packages, use_bioc = FALSE) {
    for (pkg in packages) {
      if (!requireNamespace(pkg, quietly = TRUE)) {
        tryCatch({
          if (use_bioc) {
            if (!requireNamespace("BiocManager", quietly = TRUE)) {
              install.packages("BiocManager")
            }
            BiocManager::install(pkg, dependencies = TRUE, update = FALSE)
          } else {
            install.packages(pkg, dependencies = TRUE)
          }
          cat(paste("✓ Installed:", pkg, "\n"))
        }, error = function(e) {
          cat(paste("✗ Failed to install:", pkg, "\n"))
        })
      } else {
        cat(paste("✓ Already installed:", pkg, "\n"))
      }
    }
  }
  install_if_missing(core_packages)
  install_if_missing(bioc_packages, use_bioc = TRUE)
  install_if_missing(specialized_packages)
  # Install CellChat from GitHub if needed
  if (!requireNamespace("CellChat", quietly = TRUE)) {
    tryCatch({
      devtools::install_github("sqjin/CellChat")
      cat("✓ Installed CellChat from GitHub\n")
    }, error = function(e) {
      cat("✗ Failed to install CellChat from GitHub\n")
    })
  }
  # Install sc-type from GitHub if needed
  if (!requireNamespace("scType", quietly = TRUE)) {
    cat("Note: sc-type scripts are sourced directly later, not installed as a package.\n")
  }
}

load_required_packages <- function() {
  cat("Loading required packages...\n")
  suppressPackageStartupMessages({
    library(Seurat)
    library(dplyr)
    library(ggplot2)
    library(patchwork)
    library(Matrix)
    library(hdf5r)
    library(openxlsx)
    library(readxl)
    library(tibble)
    library(SingleR)
    library(celldex)
    library(AnnotationDbi)
    library(org.Hs.eg.db)
    library(scDblFinder)
    library(SingleCellExperiment)
  })
  cat("✓ Packages loaded successfully\n")
}
# =============================================================================
# INTERACTIVE WORKSPACE (PROJECT) DIRECTORY SELECTION (tcltk-based)
# =============================================================================

#' Select a workspace/project directory interactively
#' @return Character string with the selected directory path
select_workspace_directory <- function() {
  if (capabilities("tcltk") && interactive()) {
    workspace_dir <- tcltk::tk_choose.dir(default = getwd(), caption = "Select your workspace/project directory")
    if (is.na(workspace_dir) || workspace_dir == "") {
      stop("No directory selected. Please rerun and select a valid workspace/project directory.")
    }
    return(normalizePath(workspace_dir))
  } else {
    stop("Interactive directory selection is only supported in interactive R sessions with tcltk.")
  }
}

# =============================================================================
# INTERACTIVE DATA DIRECTORY SELECTION (tcltk-based)
# =============================================================================

#' Select a data directory interactively (if different from workspace)
#' @return Character string with the selected directory path
select_data_directory <- function() {
  if (capabilities("tcltk") && interactive()) {
    data_dir <- tcltk::tk_choose.dir(default = getwd(), caption = "Select your working/data directory")
    if (is.na(data_dir) || data_dir == "") {
      stop("No directory selected. Please rerun and select a valid working/data directory.")
    }
    return(normalizePath(data_dir))
  } else {
    stop("Interactive directory selection is only supported in interactive R sessions with tcltk.")
  }
}

# =============================================================================
# DATA LOADING FUNCTION
# =============================================================================
load_data <- function(data_path) {
  if (dir.exists(data_path)) {
    # Look for .rds or .RData files first
    seurat_files <- list.files(data_path, pattern = "\\.(rds|RData)$", full.names = TRUE, ignore.case = TRUE)
    h5_files <- list.files(data_path, pattern = "\\.h5$", full.names = TRUE, ignore.case = TRUE)
    mtx_files <- list.files(data_path, pattern = "^.*matrix.*\\.mtx\\.gz$", full.names = TRUE, ignore.case = TRUE)
    all_files <- c(seurat_files, h5_files, mtx_files)
    all_files <- all_files[file.info(all_files)$isdir == FALSE]
    if (length(all_files) == 0) stop("No .rds, .RData, .h5, or matrix .mtx.gz data files found in the directory.")
    data_path <- all_files[1]
    cat("Auto-selected data file:", data_path, "\n")
  }
  if (!file.exists(data_path)) stop("Selected file does not exist: ", data_path)
  
  if (grepl("\\.rds$", data_path, ignore.case = TRUE)) {
    seurat_obj <- readRDS(data_path)
    cat(paste("✓ Loaded Seurat object with", ncol(seurat_obj), "cells and", nrow(seurat_obj), "genes\n"))
    return(seurat_obj)
  } else if (grepl("\\.RData$", data_path, ignore.case = TRUE)) {
    load(data_path)
    # Try to find the Seurat object in the environment
    seurat_objs <- Filter(function(x) inherits(get(x), "Seurat"), ls())
    if (length(seurat_objs) == 0) stop("No Seurat object found in the loaded .RData file.")
    seurat_obj <- get(seurat_objs[1])
    cat(paste("✓ Loaded Seurat object from .RData with", ncol(seurat_obj), "cells and", nrow(seurat_obj), "genes\n"))
    return(seurat_obj)
  } else if (grepl("\\.h5$", data_path, ignore.case = TRUE)) {
    data_matrix <- Read10X_h5(data_path)
    if (is.list(data_matrix)) {
      if ("Gene Expression" %in% names(data_matrix)) {
        data_matrix <- data_matrix$`Gene Expression`
      } else {
        data_matrix <- data_matrix[[1]]
      }
    }
    seurat_obj <- CreateSeuratObject(
      counts = data_matrix,
      project = "scRNA_analysis",
      min.cells = 3,
      min.features = 200
    )
    cat(paste("✓ Loaded", ncol(seurat_obj), "cells and", nrow(seurat_obj), "genes\n"))
    return(seurat_obj)
  } else if (grepl("matrix.*\\.mtx\\.gz$", data_path, ignore.case = TRUE)) {
    file_dir <- dirname(data_path)
    all_files_in_dir <- list.files(file_dir, full.names = TRUE, ignore.case = TRUE)
    barcodes_file <- grep("barcodes.*\\.tsv(\\.gz)?$", all_files_in_dir, value = TRUE, ignore.case = TRUE)
    features_file <- grep("features.*\\.tsv(\\.gz)?$", all_files_in_dir, value = TRUE, ignore.case = TRUE)
    genes_file <- grep("genes.*\\.tsv(\\.gz)?$", all_files_in_dir, value = TRUE, ignore.case = TRUE)
    if (length(features_file) == 0 && length(genes_file) > 0) {
      features_file <- genes_file
    }
    if (length(barcodes_file) == 0 || length(features_file) == 0) {
      stop("Could not find matching barcodes or features/genes file in the selected directory.")
    }
    cat("Found barcodes file:", barcodes_file[1], "\n")
    cat("Found features/genes file:", features_file[1], "\n")
    
    mtx <- Matrix::readMM(gzfile(data_path))
    barcodes <- readLines(gzfile(barcodes_file[1]))
    features <- read.delim(gzfile(features_file[1]), header = FALSE)
    gene_symbols <- features[,2]
    if (any(duplicated(gene_symbols))) {
      cat("Warning: Duplicate gene symbols found. Making unique...\n")
      gene_symbols <- make.unique(gene_symbols)
    }
    rownames(mtx) <- gene_symbols
    colnames(mtx) <- barcodes
    
    seurat_obj <- CreateSeuratObject(
      counts = mtx,
      project = "scRNA_analysis",
      min.cells = 3,
      min.features = 200
    )
    cat(paste("✓ Loaded", ncol(seurat_obj), "cells and", nrow(seurat_obj), "genes\n"))
    return(seurat_obj)
  } else {
    stop("Selected file is not a supported data file: ", data_path)
  }
}
# =============================================================================
# QUALITY CONTROL FUNCTIONS
# =============================================================================
calculate_qc_metrics <- function(seurat_obj) {
  cat("Calculating QC metrics...\n")
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  seurat_obj[["percent.ribo"]] <- PercentageFeatureSet(seurat_obj, pattern = "^RP[SL]")
  seurat_obj[["percent.hb"]] <- PercentageFeatureSet(seurat_obj, pattern = "^HB[^(P)]")
  seurat_obj$log10GenesPerUMI <- log10(seurat_obj$nFeature_RNA) / log10(seurat_obj$nCount_RNA)
  return(seurat_obj)
}

visualize_qc_metrics <- function(seurat_obj, output_dir = "plots") {
  cat("Creating QC visualizations...\n")
  p1 <- VlnPlot(seurat_obj,
                features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                ncol = 3, pt.size = 0.1) +
    patchwork::plot_annotation(title = "QC Metrics Distribution")
  p2 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.mt") +
    geom_smooth(method = "lm", se = FALSE, color = "red")
  p3 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
    geom_smooth(method = "lm", se = FALSE, color = "red")
  p4 <- ggplot(seurat_obj@meta.data, aes(x = log10GenesPerUMI)) +
    geom_histogram(bins = 50, fill = "skyblue", alpha = 0.7) +
    geom_vline(xintercept = 0.8, color = "red", linetype = "dashed") +
    labs(title = "Genes per UMI", x = "log10(Genes per UMI)", y = "Frequency") +
    theme_minimal()
  qc_plot <- (p1) / (p2 | p3 | p4)
  
  print(qc_plot) # This will display the plot in the RStudio Plots pane
  
  if (!is.null(output_dir)) {
    ggsave(file.path(output_dir, "qc_metrics.pdf"), qc_plot,
           width = 16, height = 12, dpi = 300)
  }
  return(qc_plot)
}

filter_cells <- function(seurat_obj,
                         min_features = NULL, max_features = NULL,
                         min_counts = NULL, max_counts = NULL,
                         max_mt = 20, min_complexity = 0.8) {
  cat("Filtering cells based on QC metrics...\n")
  if (is.null(min_features)) {
    median_features <- median(seurat_obj$nFeature_RNA)
    mad_features <- mad(seurat_obj$nFeature_RNA)
    min_features <- max(200, median_features - 3 * mad_features)
  }
  if (is.null(max_features)) {
    median_features <- median(seurat_obj$nFeature_RNA)
    mad_features <- mad(seurat_obj$nFeature_RNA)
    max_features <- median_features + 3 * mad_features
  }
  if (is.null(min_counts)) {
    median_counts <- median(seurat_obj$nCount_RNA)
    mad_counts <- mad(seurat_obj$nCount_RNA)
    min_counts <- max(500, median_counts - 3 * mad_counts)
  }
  if (is.null(max_counts)) {
    median_counts <- median(seurat_obj$nCount_RNA)
    mad_counts <- mad(seurat_obj$nCount_RNA)
    max_counts <- median_counts + 3 * mad_counts
  }
  cat("Filtering criteria:\n")
  cat(paste(" Features per cell:", min_features, "-", max_features, "\n"))
  cat(paste(" UMI per cell:", min_counts, "-", max_counts, "\n"))
  cat(paste(" Max mitochondrial %:", max_mt, "\n"))
  cat(paste(" Min complexity:", min_complexity, "\n"))
  cells_before <- ncol(seurat_obj)
  seurat_obj <- subset(seurat_obj,
                       subset = nFeature_RNA >= min_features &
                         nFeature_RNA <= max_features &
                         nCount_RNA >= min_counts &
                         nCount_RNA <= max_counts &
                         percent.mt <= max_mt &
                         log10GenesPerUMI >= min_complexity)
  cells_after <- ncol(seurat_obj)
  cat(paste("✓ Filtered from", cells_before, "to", cells_after, "cells\n"))
  cat(paste(" Removed", cells_before - cells_after, "cells (",
            round((cells_before - cells_after) / cells_before * 100, 1), "%)\n"))
  return(seurat_obj)
}

# =============================================================================
# DOUBLET REMOVAL (Using scDblFinder)
# =============================================================================
perform_doublet_removal <- function(seurat_obj, output_dir = "plots") {
  cat("\nPerforming doublet removal with scDblFinder...\n")
  
  # Ensure scDblFinder and SingleCellExperiment are loaded
  if (!requireNamespace("scDblFinder", quietly = TRUE) || !requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    stop("scDblFinder or SingleCellExperiment package is not loaded. Please ensure you run install_required_packages() and load_required_packages().")
  }
  suppressPackageStartupMessages(library(scDblFinder))
  suppressPackageStartupMessages(library(SingleCellExperiment))
  
  # Check if 'RNA' assay exists
  if (!"RNA" %in% names(seurat_obj@assays)) {
    stop("RNA assay not found in Seurat object. Cannot proceed with scDblFinder.")
  }
  
  # Explicitly extract the 'counts' layer from the Seurat v5 object.
  # This avoids the deprecated slot access and uses the correct 'layer' argument.
  counts_matrix <- tryCatch({
    GetAssayData(seurat_obj, assay = "RNA", layer = "counts")
  }, error = function(e) {
    stop(paste("Failed to retrieve 'counts' layer from RNA assay:", e$message,
               "\nPlease ensure your Seurat object has raw counts in the 'RNA' assay's 'counts' layer (often the default after Read10X or CreateSeuratObject)."))
  })
  
  # Create SingleCellExperiment object using the extracted counts matrix
  # and transfer the metadata (colData) directly from the Seurat object.
  sce_obj <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = counts_matrix),
    colData = seurat_obj@meta.data
  )
  
  # Run scDblFinder
  tryCatch({
    # scDblFinder will perform its own internal normalization/processing.
    # If your samples are merged, you might want to specify 'samples' argument, e.g., 'samples = seurat_obj$orig.ident'
    sce_obj <- scDblFinder(sce_obj, verbose = TRUE)
    
    # Extract doublet information and add to Seurat object
    seurat_obj$doublet_class <- sce_obj$scDblFinder.class
    seurat_obj$doublet_score <- sce_obj$scDblFinder.score
    
    cat(paste("scDblFinder classifications: ", table(seurat_obj$doublet_class)[["doublet"]], "doublets identified.\n"))
    
    # Visualize doublets
    # Ensure UMAP is computed if not already (scDblFinder itself doesn't need UMAP, but for visualization)
    if (!"umap" %in% names(seurat_obj@reductions)) {
      cat("UMAP reduction not found for visualization. Please ensure dimensionality reduction (Step 11) is run if you want to see UMAP plots.\n")
      # Skipping UMAP visualization if reduction not present, to avoid error.
      # You can add a RunUMAP here if you want it to always plot, but it's better to ensure pipeline order.
    } else {
      p_doublets <- DimPlot(seurat_obj, reduction = "umap", group.by = "doublet_class",
                            cols = c("singlet" = "grey", "doublet" = "red"),
                            pt.size = 0.5) + # Smaller points for clarity
        ggtitle("scDblFinder Classification (UMAP)") +
        theme_minimal() +
        theme(plot.title = element_text(hjust = 0.5))
      
      if (!is.null(output_dir)) {
        ggsave(file.path(output_dir, "scdblfinder_umap.pdf"), p_doublets, width = 8, height = 6)
      }
      print(p_doublets)
    }
    
    # Remove doublets
    cells_before_doublet_removal <- ncol(seurat_obj)
    seurat_obj <- subset(seurat_obj, subset = doublet_class == "singlet")
    cells_after_doublet_removal <- ncol(seurat_obj)
    
    cat(paste("✓ Removed", cells_before_doublet_removal - cells_after_doublet_removal, "doublet cells.\n"))
    
  }, error = function(e) {
    warning(paste("scDblFinder failed:", e$message, "\nSkipping doublet removal and retaining all cells."))
    seurat_obj$doublet_class <- "not_run"
    seurat_obj$doublet_score <- NA
  })
  
  return(seurat_obj)
}

# =============================================================================
# NORMALIZATION AND FEATURE SELECTION
# =============================================================================
normalize_data <- function(seurat_obj, method = "SCTransform") {
  cat(paste("Normalizing data using", method, "...\n"))
  if (method == "SCTransform") {
    seurat_obj <- SCTransform(seurat_obj,
                              vars.to.regress = "percent.mt",
                              verbose = FALSE,
                              variable.features.n = 3000)
  } else if (method == "LogNormalize") {
    seurat_obj <- NormalizeData(seurat_obj,
                                normalization.method = "LogNormalize",
                                scale.factor = 10000,
                                verbose = FALSE)
    seurat_obj <- FindVariableFeatures(seurat_obj,
                                       selection.method = "vst",
                                       nfeatures = 2000,
                                       verbose = FALSE)
    seurat_obj <- ScaleData(seurat_obj,
                            vars.to.regress = "percent.mt",
                            verbose = FALSE)
  } else {
    stop("Unsupported normalization method. Choose 'SCTransform' or 'LogNormalize'.")
  }
  cat("✓ Normalization complete\n")
  return(seurat_obj)
}

# =============================================================================
# DIMENSIONALITY REDUCTION AND CLUSTERING
# =============================================================================
perform_dimensionality_reduction <- function(seurat_obj, n_pcs = 50) {
  cat("Performing dimensionality reduction...\n")
  # Use the default assay (SCT if SCTransform was run, otherwise RNA)
  seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj),
                       npcs = n_pcs, verbose = FALSE)
  elbow_plot <- ElbowPlot(seurat_obj, ndims = n_pcs) +
    ggtitle("PCA Elbow Plot") +
    theme_minimal()
  
  # Determine optimal number of dimensions for UMAP/tSNE
  # This could be user-defined or determined from the elbow plot
  n_dims <- 30 # Default, or dynamically choose based on elbow_plot analysis
  
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:n_dims,
                        verbose = FALSE, seed.use = 42)
  seurat_obj <- RunTSNE(seurat_obj, dims = 1:n_dims,
                        verbose = FALSE, seed.use = 42)
  cat("✓ Dimensionality reduction complete\n")
  return(list(seurat_obj = seurat_obj, elbow_plot = elbow_plot))
}

find_clusters <- function(seurat_obj, resolutions = c(0.1, 0.3, 0.5, 0.8, 1.0)) {
  cat("Finding clusters...\n")
  # Use the default assay's PCA results for FindNeighbors
  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30, verbose = FALSE) # dims should match those used for UMAP/tSNE
  seurat_obj <- FindClusters(seurat_obj,
                             resolution = resolutions,
                             verbose = FALSE)
  
  # Set default resolution for Idents
  default_res_col <- paste0(DefaultAssay(seurat_obj), "_snn_res.", resolutions[3])
  if (default_res_col %in% colnames(seurat_obj@meta.data)) {
    Idents(seurat_obj) <- default_res_col
  } else {
    warning("Default resolution not found. Using first available resolution.")
    # Find the first column that looks like a clustering resolution
    Idents(seurat_obj) <- grep(paste0("^", DefaultAssay(seurat_obj), "_snn_res\\."), colnames(seurat_obj@meta.data), value = TRUE)[1]
  }
  cat(paste("✓ Found clusters at resolutions:", paste(resolutions, collapse = ", "), "\n"))
  cat(paste("✓ Default resolution set to:", resolutions[3], "\n"))
  return(seurat_obj)
}

visualize_clusters <- function(seurat_obj, output_dir = "plots") {
  cat("Creating cluster visualizations...\n")
  p1 <- DimPlot(seurat_obj, reduction = "umap", label = TRUE,
                pt.size = 0.5, label.size = 6) +
    ggtitle("Clusters (UMAP)") +
    theme_minimal() +
    NoLegend()
  p2 <- DimPlot(seurat_obj, reduction = "tsne", label = TRUE,
                pt.size = 0.5, label.size = 6) +
    ggtitle("Clusters (t-SNE)") +
    theme_minimal() +
    NoLegend()
  p3 <- FeaturePlot(seurat_obj, features = "nFeature_RNA",
                    reduction = "umap", pt.size = 0.5) +
    scale_color_viridis_c() +
    ggtitle("nFeature_RNA")
  p4 <- FeaturePlot(seurat_obj, features = "percent.mt",
                    reduction = "umap", pt.size = 0.5) +
    scale_color_viridis_c() +
    ggtitle("Mitochondrial %")
  cluster_plot <- (p1) / (p2 | p3 | p4)
  if (!is.null(output_dir)) {
    ggsave(file.path(output_dir, "clusters.pdf"), cluster_plot,
           width = 16, height = 12, dpi = 300)
  }
  print(cluster_plot)
  return(cluster_plot)
}
# =============================================================================
# ENHANCED CELL TYPE ANNOTATION WITH MALIGNANCY & CAF DETECTION
# (Adapted for CancerSCEM alignment: SingleR added, marker lists adjusted for example)
# =============================================================================
get_cancerscem_markers <- function() {
  # Load HGNC helper if needed
  if (!requireNamespace("HGNChelper", quietly = TRUE)) {
    install.packages("HGNChelper")
  }
  suppressPackageStartupMessages(library(HGNChelper))
  
  # Define markers and validate symbols
  markers <- list(
    "Astrocyte" = c("AGXT2L1", "GFAP", "ALDOC", "SLC1A3", "AGT", "ALDH1L1"),
    "B cell" = c("CD19", "MS4A1", "BANK1", "BLK", "IRF8", "ABCB4", "ABCB9", "AFF4", "AIDA", "AIM2"),
    "Endothelial cell" = c("VWF", "PECAM1", "CDH5", "VEGFA", "FLT1", "ECSCR", "ACYP1", "ADGRL2", "SELE", "ICAM1"),
    "Epithelial cell" = c("CDH1", "MYLK", "ANKRD30A", "ABCA13", "ABCB10", "ADGB", "SFTPB", "SFTPC"),
    "Erythrocyte" = c("ALAS2", "CA1", "HBB", "HBE1", "HBA1", "HBG1", "GYPA"),
    "Fibroblast" = c("COL1A1", "COL3A1", "THY1", "NECTIN1", "FAP", "PTPN13", "C5AR2", "LRP1"),
    "GMP" = c("CD38", "KIT", "ADK", "CD123", "ALDH4A1", "ANXA1", "AP3S1", "APLP2", "APPL1", "AREG", "ASPM", "CDKN3", "CLSPN", "DEPDC7", "MCM10", "MUCB2", "SDC4", "RMI2"),
    "HSC" = c("CD34", "ITGA5", "PROM1", "CD105", "VCAM1", "CD164", "THY1", "KIT", "ACE", "CMAH", "ABCG2", "CD41", "ALDH1A1", "BMI1"),
    "Macro/Mono/DC" = c("CD68", "CD14", "MRC1", "BHLHE40", "CD93", "CREM", "CSF1R", "CCL18", "ICAM4", "ACPP", "ACSL3", "ADGRE2", "ADGRE3", "CD209", "CD83", "CD1A"),
    "Malignant cell" = c("EPCAM", "FOLH1", "KLK3", "KRT8", "KRT18", "KRT19"),
    "Mast cell" = c("SLC18A2", "ADIRF", "ASIC4", "BACE2", "ENPP3", "CADPS", "CAPN3", "CDK15", "CMA1", "GCSAML", "MAML1", "MAOB", "CAVIN2"),
    "Melanocyte" = c("MLANA", "PMEL", "DCT", "KIT", "MITF", "TYR", "MC1R", "OCA2", "BCL2"),
    "Muscle cell" = c("MYH2", "MYL2", "ACTA1", "CKM", "MYOM3", "CRYAB", "CMD", "APLN", "HOMG4"),
    "Myeloid cell" = c("PTPRC", "CD14", "AIF1", "TYROBP", "CD163"),
    "Neural progenitor cell" = c("SOX2", "NESTIN", "DCX", "NES1", "PAX6", "CENPF", "UBE2C", "TMN2", "MKI67"),
    "Neuron" = c("STMN2", "RBFOX3", "MAP2", "TUBB3", "CSF3", "DLG4", "ENO2"),
    "Neutrophil" = c("ADGRG3", "CXCL8", "FCGR3B", "MNDA", "USP10", "CSF3R", "ANXA3", "AQP9", "BTNL8", "LGALS13", "G0S2", "NFE4", "IL5RA"),
    "NK cell" = c("FCGR3A", "KLRB1", "KLRD1", "NKG7", "XCL1", "XCL2", "NCR3", "NCR1", "CD247", "GZMB", "KLRC1", "KLRK1"),
    "Oligodendrocyte" = c("MOG", "OLIG1", "OLIG2", "PDGFRA", "PLP1", "MBP", "MAG", "SOX10"),
    "Oligodendrocyte precursor cell" = c("PDGFRALPHA", "OLIG2", "CSPG4", "OLIG1", "VCAN", "LHFPL3", "GPR17", "APOD"),
    "Pericyte" = c("CSPG4", "PDGFRB", "MCAM", "RGS5", "ALPHA-SMA", "KCNJ8", "PDGFRALPHA", "THY1"),
    "Plasma cell" = c("MZB1", "BRSK1", "AC026202.3", "JSRP1", "LINC00582", "PARM1", "TAS1R3"),
    "Plasmacytoid dendritic cell" = c("CLE4C", "IL3RA", "LILRA4", "GZMB", "JCHAIN", "IRF7", "TCF4", "NRP1", "IRF8"),
    "Progenitor" = c("CD38", "CASR", "ALDH", "CAR", "KDR", "MME", "FLT3", "CD90", "CD123"),
    "Retinal ganglion cell" = c("RBPMS", "NESTIN", "PACA", "GAP43", "MAP2", "ATOH1", "NEFM", "NEFH"),
    "Smooth muscle cell" = c("TAGLN", "MYH11", "DES", "CNN1", "RGS5", "MYL9", "NOTCH3", "SYNPO2", "PLN"),
    "T cell" = c("CD3D", "CD3G", "CD3E"),
    "Treg" = c("IL2RA", "FOXP3", "CTLA4", "LRRC32", "TNFRSF18", "TNFRSF4")
  )
  # Validate and correct gene symbols
  for (i in seq_along(markers)) {
    markers[[i]] <- unique(toupper(markers[[i]])) # Convert to uppercase
    validated <- HGNChelper::checkGeneSymbols(markers[[i]])
    markers[[i]] <- validated$Suggested.Symbol
  }
  
  return(markers)
}

# =============================================================================
# MODIFIED annotate_cell_types function with proper scType setup
# =============================================================================
annotate_cell_types <- function(seurat_obj, tissue = "Liver", output_dir = getwd()) {
  cat("Performing automated cell type annotation...\n")
  
  # Initialize celltype_auto column to 'Unassigned'
  seurat_obj$celltype_auto <- "Unassigned"
  
  # Get CancerSCEM markers
  cancerscem_markers <- get_cancerscem_markers()
  
  # Add module scores for each cell type
  for (celltype in names(cancerscem_markers)) {
    markers <- cancerscem_markers[[celltype]]
    valid_markers <- markers[markers %in% rownames(seurat_obj)]
    
    if (length(valid_markers) > 0) {
      seurat_obj <- AddModuleScore(
        object = seurat_obj,
        features = list(valid_markers),
        name = paste0("score_", celltype),
        assay = DefaultAssay(seurat_obj)
      )
    }
  }
  
  # Get all score columns
  score_cols <- grep("^score_", colnames(seurat_obj@meta.data), value = TRUE)
  
  # For each cell, assign the cell type with the highest score
  if (length(score_cols) > 0) {
    score_matrix <- seurat_obj@meta.data[, score_cols, drop = FALSE]
    max_scores <- apply(score_matrix, 1, max)
    best_types <- apply(score_matrix, 1, function(x) {
      if (max(x) > 0) {
        gsub("^score_", "", score_cols[which.max(x)])
      } else {
        "Unassigned"
      }
    })
    
    seurat_obj$celltype_auto <- best_types
  }
  
  return(seurat_obj)
}
# =============================================================================
# DETECTION OF CAFS
# =============================================================================
detect_CAFs <- function(seurat_obj, output_dir = getwd()) {
  cat("\nPerforming Cancer-Associated Fibroblast (CAF) detection...\n")
  
  # Define CAF markers
  caf_markers <- c("COL1A1", "COL3A1", "ACTA2", "PDGFRA", "PDGFRB", "THY1", "POSTN")
  
  # Filter to markers present in the dataset
  valid_caf_markers <- caf_markers[caf_markers %in% rownames(seurat_obj)]
  cat("Using CAF markers found in data:", paste(valid_caf_markers, collapse = ", "), "\n")
  
  # Calculate CAF score
  if (length(valid_caf_markers) == 0) {
    warning("No valid CAF markers found in dataset! Setting all CAF scores to 0.")
    seurat_obj$caf_score <- 0
  } else {
    seurat_obj <- AddModuleScore(
      object = seurat_obj,
      features = list(valid_caf_markers),
      name = "caf_score",
      assay = DefaultAssay(seurat_obj)
    )
    # Clean up score column naming
    if ("caf_score1" %in% colnames(seurat_obj@meta.data)) {
      seurat_obj$caf_score <- seurat_obj@meta.data$caf_score1
      seurat_obj@meta.data$caf_score1 <- NULL
    }
  }
  
  # Assign CAF calls with automatic threshold
  caf_threshold <- quantile(seurat_obj$caf_score, probs = 0.9, na.rm = TRUE)
  seurat_obj$CAF_call <- ifelse(seurat_obj$caf_score > caf_threshold, "CAF", "Non-CAF")
  
  # IMPLEMENT PRIORITY ANNOTATION SYSTEM:
  # 1. First use CancerSCEM annotation (celltype_auto)
  # 2. Then override with CAF identity if applicable
  seurat_obj$celltype_priority <- ifelse(
    seurat_obj$CAF_call == "CAF",
    "CAF",
    as.character(seurat_obj$celltype_auto)
  )
  
  # Visualization
  if ("umap" %in% names(seurat_obj@reductions)) {
    # Plot showing the priority annotations
    p_priority <- DimPlot(seurat_obj, 
                          group.by = "celltype_priority",
                          reduction = "umap",
                          label = TRUE,
                          repel = TRUE) +
      ggtitle("Cell Types (Priority: CAF > CancerSCEM)") +
      theme_minimal()
    
    ggsave(file.path(output_dir, "celltype_priority_umap.pdf"), p_priority, width = 10, height = 8)
    print(p_priority)
  }
  
  cat("✓ CAF detection complete\n")
  return(seurat_obj)
}
# =============================================================================
# DIFFERENTIAL EXPRESSION ANALYSIS AND MARKER FINDING
# =============================================================================
find_cluster_markers <- function(seurat_obj, output_dir = "plots", min_pct = 0.25, logfc_threshold = 0.25) {
  cat("\nFinding cluster markers...\n")
  # Find all markers for each cluster
  seurat_obj.markers <- FindAllMarkers(seurat_obj,
                                       only.pos = TRUE,
                                       min.pct = min_pct,
                                       logfc.threshold = logfc_threshold,
                                       verbose = FALSE)
  
  # Save top markers
  if (!is.null(output_dir) && !is.null(seurat_obj.markers) && nrow(seurat_obj.markers) > 0) {
    write.xlsx(seurat_obj.markers, file.path(output_dir, "all_cluster_markers.xlsx"))
    
    top_n_markers <- seurat_obj.markers %>%
      group_by(cluster) %>%
      top_n(n = 10, wt = avg_log2FC)
    
    write.xlsx(top_n_markers, file.path(output_dir, "top10_cluster_markers.xlsx"))
    
    # Visualize top 10 markers per cluster
    if (nrow(top_n_markers) > 0) {
      cat("Creating heatmap of top markers...\n")
      # Handle potential empty clusters or clusters with too few cells for plotting
      valid_clusters_for_heatmap <- unique(top_n_markers$cluster)
      valid_cells_for_heatmap <- WhichCells(seurat_obj, idents = valid_clusters_for_heatmap)
      
      if (length(valid_cells_for_heatmap) > 0) {
        heatmap_plot <- DoHeatmap(subset(seurat_obj, cells = valid_cells_for_heatmap),
                                  features = top_n_markers$gene,
                                  group.by = "seurat_clusters",
                                  angle = 90,
                                  size = 3) +
          scale_fill_gradientn(colors = c("blue", "white", "red")) +
          theme(axis.text.y = element_text(size = 8))
        
        ggsave(file.path(output_dir, "top_cluster_markers_heatmap.pdf"), heatmap_plot,
               width = 12, height = 18, limitsize = FALSE)
        print(heatmap_plot)
      } else {
        warning("No valid cells found in selected clusters to plot heatmap for top markers. Skipping heatmap generation.")
      }
    } else {
      warning("No top markers found to visualize. Skipping heatmap generation.")
    }
  } else {
    warning("No markers found or output directory not specified. Skipping marker saving and heatmap generation.")
  }
  cat("✓ Cluster marker finding complete\n")
  return(seurat_obj.markers)
}

# =============================================================================
# GENE EXPRESSION ANALYSIS
# =============================================================================
analyze_gene_expression <- function(seurat_obj, gene, output_dir = "plots") {
  cat(paste0("\nAnalyzing expression of gene: ", gene, "...\n"))
  
  if (!gene %in% rownames(seurat_obj)) {
    warning(paste0("Gene '", gene, "' not found in the dataset. Skipping analysis for this gene."))
    return(NULL)
  }
  
  # UMAP FeaturePlot
  p_umap_gene <- FeaturePlot(seurat_obj, features = gene, reduction = "umap",
                             pt.size = 0.5, order = TRUE) +
    scale_color_viridis_c() +
    ggtitle(paste0(gene, " Expression (UMAP)")) +
    theme_minimal()
  
  # Violin Plot by cluster
  p_vln_gene <- VlnPlot(seurat_obj, features = gene, pt.size = 0.1) +
    ggtitle(paste0(gene, " Expression by Cluster")) +
    theme_minimal() +
    NoLegend()
  
  # Combine plots
  gene_plot <- p_umap_gene + p_vln_gene
  
  if (!is.null(output_dir)) {
    ggsave(file.path(output_dir, paste0(gene, "_expression.pdf")), gene_plot,
           width = 12, height = 6, dpi = 300)
  }
  print(gene_plot)
  cat(paste0("✓ Analysis for gene '", gene, "' complete.\n"))
  return(gene_plot)
}


# =============================================================================
# MAIN PIPELINE EXECUTION FUNCTION (MODIFIED for interactive menu)
# =============================================================================
run_sc_pipeline <- function(start_num = NULL) { # Make start_num optional
  # Create plots directory
  plots_dir <- "plots"
  if (!dir.exists(plots_dir)) {
    dir.create(plots_dir)
  }
  
  # Initialize seurat_obj and workspace_dir outside if blocks
  seurat_obj <- NULL
  workspace_dir <- NULL
  
  # Attempt to load existing Seurat object at the very beginning
  cat("\nAttempting to load a previously saved Seurat object...\n")
  potential_clustered_rds_path <- file.path(getwd(), "seurat_obj_after_clustering.rds")
  potential_final_rds_path <- file.path(getwd(), "seurat_obj_processed.rds")
  
  if (file.exists(potential_clustered_rds_path)) {
    cat("Found 'seurat_obj_after_clustering.rds'. Loading...\n")
    seurat_obj <- readRDS(potential_clustered_rds_path)
    cat("✓ Loaded Seurat object after clustering.\n")
    workspace_dir <- getwd()
  } else if (file.exists(potential_final_rds_path)) {
    cat("Found 'seurat_obj_processed.rds'. Loading...\n")
    seurat_obj <- readRDS(potential_final_rds_path)
    cat("✓ Loaded previous final Seurat object.\n")
    workspace_dir <- getwd()
  } else {
    cat("No previously saved Seurat object found. You may need to start from Step 4: Load data.\n")
  }
  
  # Define pipeline steps for interactive menu (UPDATED NUMBERS)
  pipeline_steps <- list(
    "Install required packages",
    "Load required packages",
    "Select workspace directory",
    "Select data directory and load data",
    "Calculate QC metrics",
    "Visualize QC metrics (pre-filtering)",
    "Filter cells",
    "Visualize QC metrics (post-filtering)",
    "Normalize data",
    "Perform doublet removal",
    "Dimensionality reduction (PCA, UMAP, tSNE)",
    "Find and visualize clusters",
    "Save Seurat object after clustering", # NEW STEP 13
    "Automated Cell Type Annotation (SingleR/scType)",
    "Malignancy and CAF detection & Enhanced UMAP",
    "Gene expression analysis",
    "Save processed Seurat object"
  )
  
  # Interactive menu for choosing starting step
  if (is.null(start_num)) {
    repeat {
      cat("\n--- Choose a starting step for the pipeline ---\n")
      for (i in seq_along(pipeline_steps)) {
        # Indicate if an object is already loaded and which step it covers implicitly
        status_msg <- ""
        if (!is.null(seurat_obj)) {
          if (i <= 4) status_msg <- " (Choosing this step will overwrite the currently loaded object.)"
          else if (i == 5 && !"nCount_RNA" %in% colnames(seurat_obj@meta.data)) status_msg <- " (Loaded object might not have QC metrics calculated. Consider starting from Step 5 or earlier.)"
          else if (i == 7 && (!"nFeature_RNA" %in% colnames(seurat_obj@meta.data) || !"nCount_RNA" %in% colnames(seurat_obj@meta.data))) status_msg <- " (Loaded object might not have QC metrics. Consider starting from Step 5 or earlier.)"
          else if (i == 9 && DefaultAssay(seurat_obj) != "SCT" && !"RNA" %in% names(seurat_obj@assays) ) status_msg <- " (Loaded object might not be normalized. Consider starting from Step 9 or earlier.)"
          else if (i == 10 && !"doublet_class" %in% colnames(seurat_obj@meta.data)) status_msg <- " (Loaded object might not have doublet info. Consider starting from Step 10 or earlier.)"
          else if (i == 11 && is.null(seurat_obj@reductions$pca)) status_msg <- " (Loaded object might not have PCA. Consider starting from Step 11 or earlier.)"
          else if (i == 12 && is.null(seurat_obj@meta.data$seurat_clusters)) status_msg <- " (Loaded object might not have clusters. Consider starting from Step 12 or earlier.)"
          else if (i > 13 && file.exists(potential_final_rds_path) && !file.exists(potential_clustered_rds_path)) status_msg <- " (Loaded object is the final processed one, might be ahead of this step.)"
          else if (i > 13 && file.exists(potential_clustered_rds_path)) status_msg <- " (Loaded object is processed up to clustering.)"
        }
        cat(paste0(i, ": ", pipeline_steps[[i]], status_msg, "\n"))
      }
      cat("----------------------------------------------\n")
      
      user_input <- readline(prompt = "Enter the number of the step to start from (e.g., 1 for beginning): ")
      selected_step <- as.integer(user_input)
      
      if (!is.na(selected_step) && selected_step >= 1 && selected_step <= length(pipeline_steps)) {
        start_num <- selected_step
        cat(paste0("Starting pipeline from Step ", start_num, ": ", pipeline_steps[[start_num]], "\n"))
        break
      } else {
        cat("Invalid input. Please enter a valid step number.\n")
      }
    }
  }
  
  # Ensure seurat_obj and workspace_dir are available for subsequent steps
  # If starting at a step > 4 and no object was loaded initially, this will stop execution.
  if (start_num > 4 && is.null(seurat_obj)) {
    stop("Seurat object is NULL. It must be loaded in Step 4 or from a saved file to proceed.")
  }
  
  # Step 1: Install packages
  if (start_num <= 1) {
    cat("\n=== Step 1: Installing required packages ===\n")
    install_required_packages()
  }
  
  # Step 2: Load packages
  if (start_num <= 2) {
    cat("\n=== Step 2: Loading required packages ===\n")
    load_required_packages()
  }
  
  # Step 3: Select workspace directory
  if (start_num <= 3) {
    cat("\n=== Step 3: Selecting workspace directory ===\n")
    # Only ask for workspace if it wasn't already set by loading a saved object
    if (is.null(workspace_dir) || workspace_dir == "") {
      workspace_dir <- select_workspace_directory()
    }
    cat("Selected workspace directory:", workspace_dir, "\n")
    setwd(workspace_dir) # Set working directory
  }
  
  # Step 4: Select data directory and load data
  if (start_num <= 4) {
    cat("\n=== Step 4: Selecting data directory and loading data ===\n")
    data_dir <- select_data_directory()
    cat("Selected data directory:", data_dir, "\n")
    seurat_obj <- load_data(data_dir)
    cat(paste("Initial cell count:", ncol(seurat_obj), "\n"))
  }
  
  # Check if seurat_obj is available before proceeding with steps that require it
  if (is.null(seurat_obj) && start_num > 4) {
    stop("Seurat object is NULL. It must be loaded in Step 4 or from a saved file to proceed.")
  }
  
  # Step 5: Calculate QC metrics
  if (start_num <= 5) {
    cat("\n=== Step 5: Calculating QC metrics ===\n")
    seurat_obj <- calculate_qc_metrics(seurat_obj)
  }
  
  # Step 6: Visualize QC metrics (pre-filtering)
  if (start_num <= 6) {
    cat("\n=== Step 6: Visualizing QC metrics (pre-filtering) ===\n")
    visualize_qc_metrics(seurat_obj, output_dir = plots_dir)
  }
  
  # Step 7: Filter cells
  if (start_num <= 7) {
    cat("\n=== Step 7: Filtering cells ===\n")
    # Using default parameters here, adjust as needed
    seurat_obj <- filter_cells(seurat_obj, max_mt = 20, min_complexity = 0.8)
    cat(paste("Cells after filtering:", ncol(seurat_obj), "\n"))
  }
  
  # Step 8: Visualize QC metrics (post-filtering)
  if (start_num <= 8) {
    cat("\n=== Step 8: Visualizing QC metrics (post-filtering) ===\n")
    visualize_qc_metrics(seurat_obj, output_dir = plots_dir)
  }
  
  # Step 9: Normalize data
  if (start_num <= 9) {
    cat("\n=== Step 9: Normalizing data ===\n")
    seurat_obj <- normalize_data(seurat_obj, method = "SCTransform")
  }
  
  # Step 10: Perform doublet removal
  if (start_num <= 10) {
    cat("\n=== Step 10: Performing doublet removal ===\n")
    if (ncol(seurat_obj) < 500) { # Add cell count check for reliable scDblFinder
      warning("Too few cells (less than 500) after filtering to reliably perform scDblFinder. Skipping doublet removal.")
      seurat_obj$doublet_class <- "not_run"
      seurat_obj$doublet_score <- NA
    } else {
      seurat_obj <- perform_doublet_removal(seurat_obj, output_dir = plots_dir)
    }
    cat(paste("Cells after doublet removal:", ncol(seurat_obj), "\n"))
  }
  
  # Step 11: Dimensionality reduction (PCA, UMAP, tSNE)
  if (start_num <= 11) {
    cat("\n=== Step 11: Performing dimensionality reduction ===\n")
    dim_red_results <- perform_dimensionality_reduction(seurat_obj)
    seurat_obj <- dim_red_results$seurat_obj
    
    if (!is.null(plots_dir)) {
      ggsave(file.path(plots_dir, "pca_elbow_plot.pdf"), dim_red_results$elbow_plot,
             width = 8, height = 6, dpi = 300)
      print(dim_red_results$elbow_plot)
    }
  }
  
  # Step 12: Find and visualize clusters
  if (start_num <= 12) {
    cat("\n=== Step 12: Finding and visualizing clusters ===\n")
    seurat_obj <- find_clusters(seurat_obj, resolutions = c(0.1, 0.3, 0.5, 0.8, 1.0))
    visualize_clusters(seurat_obj, output_dir = plots_dir)
  }
  
  # NEW Step 13: Save Seurat object after clustering
  if (start_num <= 13) {
    cat("\n=== Step 13: Saving Seurat object after clustering ===\n")
    # Ensure workspace_dir is set
    if (is.null(workspace_dir) || workspace_dir == "") {
      warning("Workspace directory not explicitly set. Using current working directory for saving.")
      workspace_dir <- getwd()
    }
    saveRDS(seurat_obj, file = file.path(workspace_dir, "seurat_obj_after_clustering.rds"))
    cat("✓ Seurat object after clustering saved to seurat_obj_after_clustering.rds\n")
  }
  
  # Step 14: Automated Cell Type Annotation (now only scType + CancerSCEM)
  if (start_num <= 14) {
    cat("\n=== Step 14: Performing automated cell type annotation ===\n")
    seurat_obj <- annotate_cell_types(seurat_obj, tissue = "Liver", output_dir = plots_dir)
  }
  # Step 15: CAF detection and enhanced UMAP
  if (start_num <= 15) {
    cat("\n=== Step 15: Performing CAF detection ===\n")
    seurat_obj <- detect_CAFs(seurat_obj, output_dir = plots_dir)
    
    # Create enhanced UMAP using priority annotations
    p_umap_enhanced <- DimPlot(seurat_obj, 
                               group.by = "celltype_priority", 
                               reduction = "umap", 
                               label = TRUE, 
                               repel = TRUE, 
                               label.size = 3) +
      ggtitle("Cell Types (Priority: CAF > CancerSCEM)") +
      theme_minimal() +
      theme(legend.position = "right")
    
    ggsave(file.path(plots_dir, "umap_celltypes_enhanced.pdf"), p_umap_enhanced, width = 10, height = 8)
    print(p_umap_enhanced)
  }
  
  # Step 16: Gene expression analysis (was 15)
  if (start_num <= 16) {
    repeat {
      cat("\n=== Step 16: Performing gene expression analysis ===\n")
      gene_of_interest <- readline(prompt = "Enter the gene symbol you want to analyze (or 'done' to finish): ")
      
      if (tolower(gene_of_interest) == "done") {
        cat("Finished gene expression analysis.\n")
        break
      }
      
      if (gene_of_interest == "") {
        gene_of_interest <- "FAP"
        cat("No gene entered. Defaulting to FAP.\n")
      }
      
      # Run gene analysis
      gene_analysis <- analyze_gene_expression(seurat_obj, gene = gene_of_interest, output_dir = plots_dir)
      
      # Ask if user wants to analyze another gene
      another_gene <- readline(prompt = "Would you like to analyze another gene? (y/n): ")
      if (tolower(another_gene) %in% c("n", "no")) {
        break
      }
    }
  }
  
  # Step 17: Save processed Seurat object (was 16)
  if (start_num <= 17) {
    cat("\n=== Step 17: Saving processed object ===\n")
    # Ensure workspace_dir is set
    if (is.null(workspace_dir) || workspace_dir == "") {
      warning("Workspace directory not explicitly set. Using current working directory for saving.")
      workspace_dir <- getwd()
    }
    saveRDS(seurat_obj, file = file.path(workspace_dir, "seurat_obj_processed.rds"))
    cat("✓ Processed Seurat object saved to seurat_obj_processed.rds\n")
  }
  
  cat("\nPipeline execution complete!\n")
}
run_sc_pipeline()