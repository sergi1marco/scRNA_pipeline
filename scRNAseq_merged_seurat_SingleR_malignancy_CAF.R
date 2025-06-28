# =============================================================================
# COMPREHENSIVE SINGLE-CELL RNA-SEQ ANALYSIS PIPELINE
# =============================================================================
# Following Luecken & Theis (2019) and scRNA-seq best practices
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
    "Matrix", "R.utils", "reticulate", "hdf5r", "devtools", "openxlsx"
  )
  bioc_packages <- c(
    "scater", "scran", "BiocParallel"
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
    library(celldex)
    library(scuttle)
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
  seurat_obj <- RunPCA(seurat_obj, features = NULL,
                       npcs = n_pcs, verbose = FALSE)
  elbow_plot <- ElbowPlot(seurat_obj, ndims = n_pcs) +
    ggtitle("PCA Elbow Plot") +
    theme_minimal()
  n_dims <- 30
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:n_dims,
                        verbose = FALSE, seed.use = 42)
  seurat_obj <- RunTSNE(seurat_obj, dims = 1:n_dims,
                        verbose = FALSE, seed.use = 42)
  cat("✓ Dimensionality reduction complete\n")
  return(list(seurat_obj = seurat_obj, elbow_plot = elbow_plot))
}

find_clusters <- function(seurat_obj, resolutions = c(0.1, 0.3, 0.5, 0.8, 1.0)) {
  cat("Finding clusters...\n")
  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30, verbose = FALSE)
  seurat_obj <- FindClusters(seurat_obj,
                             resolution = resolutions,
                             verbose = FALSE)
  default_res <- paste0("SCT_snn_res.", resolutions[3]) # Usually 0.5
  if (default_res %in% colnames(seurat_obj@meta.data)) {
    Idents(seurat_obj) <- default_res
  } else {
    warning("Default resolution not found. Using available resolution.")
    Idents(seurat_obj) <- grep("^SCT_snn_res", colnames(seurat_obj@meta.data), value = TRUE)[1]
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
  cluster_plot <- (p1 | p2) / (p3 | p4)
  if (!is.null(output_dir)) {
    ggsave(file.path(output_dir, "clusters.pdf"), cluster_plot,
           width = 16, height = 12, dpi = 300)
  }
  return(cluster_plot)
}
# =============================================================================
# COMPLETE CELL TYPE ANNOTATION MODULE (SINGLER + MALIGNANCY + CAF)
# =============================================================================

#' Execute Complete Cell Typing Pipeline
#' @param seurat_obj Processed Seurat object with clustering
#' @param tissue_type Tissue type ("Liver" by default)
#' @param output_dir Directory to save results
#' @param min_confidence Minimum SingleR confidence score (0-1)
complete_cell_typing <- function(seurat_obj, tissue_type = "Liver", output_dir = "celltype_results", min_confidence = 0.4, min_cells = 24) {
  
  # Create output directory
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  cat("\n===== STEP 1: SINGLER ANNOTATION =====\n")
  seurat_obj <- run_singler_annotation(seurat_obj, tissue_type, output_dir, min_confidence)
  
  cat("\n===== STEP 2: MALIGNANCY DETECTION =====\n")
  seurat_obj <- detect_malignancy(seurat_obj, tissue_type, output_dir)
  
  cat("\n===== STEP 3: CAF DETECTION=====\n")
  seurat_obj <- detect_CAFs(seurat_obj, output_dir)
  
  cat("\n===== STEP 4: FINAL ANNOTATIONS =====\n")
  seurat_obj <- create_final_annotations(seurat_obj, output_dir)
  
  cat("\n===== STEP 5: IMMUNE CELL SUBTYPING =====\n")
  if (sum(seurat_obj$celltype_final == "Immune system cells") >= min_cells) {
    seurat_obj <- refine_immune_subtypes(
      seurat_obj,
      output_dir = file.path(output_dir, "immune_subtyping"),
      min_cells = 24  # Minimum cells required
    )
  } else {
    warning("Skipping immune subtyping: <", min_cells, " immune cells detected")
    seurat_obj$immune_subtype <- NA
  }
  
  cat("\n===== STEP 6: CAF SUBTYPING =====\n")
  if (sum(seurat_obj$celltype_final == "CAF") >= min_cells) {
    seurat_obj <- subtype_CAFs(seurat_obj, output_dir, min_cells)
  } else {
    warning("Skipping CAF subtyping: <", min_cells, " CAF cells detected")
    seurat_obj$caf_subtype <- NA
  }
  
  # Create enhanced annotations combining main types and subtypes
  seurat_obj$celltype_enhanced <- seurat_obj$celltype_final
  
  # Update with immune subtypes if available
  if ("immune_subtype" %in% colnames(seurat_obj@meta.data)) {
    immune_cells <- !is.na(seurat_obj$immune_subtype)
    seurat_obj$celltype_enhanced[immune_cells] <- paste(
      "Immune", 
      seurat_obj$immune_subtype[immune_cells]
    )
  }
  
  # Update with CAF subtypes if available
  if ("caf_subtype" %in% colnames(seurat_obj@meta.data)) {
    caf_cells <- !is.na(seurat_obj$caf_subtype)
    seurat_obj$celltype_enhanced[caf_cells] <- paste(
      "CAF", 
      seurat_obj$caf_subtype[caf_cells]
    )
  }
  
  cat("\n===== STEP 7: VISUALIZE FINAL RESULTS =====\n")
  p_enhanced <- DimPlot(seurat_obj, group.by = "celltype_enhanced", 
                        reduction = "umap", label = TRUE, repel = TRUE, 
                        pt.size = 0.7, label.size = 4) +
    ggtitle("Enhanced Cell Type Annotation") +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold", size = 16))
  
  ggsave(file.path(output_dir, "enhanced_annotation.pdf"), p_enhanced, width = 12, height = 8)
  print(p_enhanced)
  
  # Save final metadata
  final_meta <- data.frame(
    Cell = colnames(seurat_obj),
    Cluster = seurat_obj$seurat_clusters,
    SingleR_Type = seurat_obj$singler_type,
    SingleR_Score = seurat_obj$singler_score,
    Malignancy_Score = seurat_obj$malignancy_score,
    CAF_Score = seurat_obj$CAF_score,
    Final_Type = seurat_obj$celltype_final,
    Immune_Subtype = seurat_obj$immune_subtype,
    CAF_Subtype = seurat_obj$caf_subtype,
    Enhanced_Type = seurat_obj$celltype_enhanced,
    stringsAsFactors = FALSE
  )
  write.csv(final_meta, file.path(output_dir, "final_annotations.csv"), row.names = FALSE)
  
  return(seurat_obj)
}
# =============================================================================
# SINGLER ANNOTATION CORE FUNCTION
# =============================================================================

run_singler_annotation <- function(seurat_obj, tissue_type, output_dir, min_confidence) {
  
  cat("\n=== RUNNING SINGLER ANNOTATION ===\n")
  
  # Load required packages silently
  suppressPackageStartupMessages({
    library(SingleR)
    library(celldex)
    library(scuttle)
  })
  
  # Set default assay and convert to SCE
  DefaultAssay(seurat_obj) <- "SCT"
  sce <- as.SingleCellExperiment(seurat_obj)
  
  # Load references
  hpca <- HumanPrimaryCellAtlasData()
  monaco <- MonacoImmuneData()
  
  # Run SingleR (cluster-level)
  pred_hpca <- SingleR(
    test = sce,
    ref = hpca,
    labels = hpca$label.main,
    clusters = seurat_obj$seurat_clusters,
    assay.type.test = "logcounts"
  )
  
  pred_monaco <- SingleR(
    test = sce,
    ref = monaco,
    labels = monaco$label.fine,
    clusters = seurat_obj$seurat_clusters,
    assay.type.test = "logcounts"
  )
  
  # Create consensus annotations
  consensus_df <- data.frame(
    cluster = rownames(pred_hpca),
    hpca_type = pred_hpca$labels,
    hpca_score = apply(pred_hpca$scores, 1, max),
    monaco_type = pred_monaco$labels,
    monaco_score = apply(pred_monaco$scores, 1, max),
    stringsAsFactors = FALSE
  )
  
  # Determine final type (Monaco for immune, HPCA for others)
  consensus_df$final_type <- ifelse(
    grepl("T cell|B cell|Macrophage|Monocyte|NK cell", consensus_df$hpca_type),
    consensus_df$monaco_type,
    consensus_df$hpca_type
  )
  
  consensus_df$final_score <- ifelse(
    grepl("T cell|B cell|Macrophage|Monocyte|NK cell", consensus_df$hpca_type),
    consensus_df$monaco_score,
    consensus_df$hpca_score
  )
  
  # Clean cell type names
  consensus_df$final_type <- gsub("_", " ", consensus_df$final_type)
  consensus_df$final_type <- gsub(" cells$", "", consensus_df$final_type)
  
  # Map to cells
  seurat_obj$singler_type <- consensus_df$final_type[match(as.character(seurat_obj$seurat_clusters), consensus_df$cluster)]
  seurat_obj$singler_score <- consensus_df$final_score[match(as.character(seurat_obj$seurat_clusters), consensus_df$cluster)]
  
  # Filter low confidence
  seurat_obj$singler_type[seurat_obj$singler_score < min_confidence] <- "Low_confidence"
  
  # Save results
  write.csv(consensus_df, file.path(output_dir, "singler_annotations.csv"), row.names = FALSE)
  
  # Plot
  p <- DimPlot(seurat_obj, group.by = "singler_type", label = TRUE, repel = TRUE) +
    ggtitle("SingleR Cell Type Annotation")
  ggsave(file.path(output_dir, "singler_annotation.pdf"), p, width = 10, height = 8)
  
  cat("SingleR identified", length(unique(seurat_obj$singler_type)), "cell types\n")
  print(table(seurat_obj$singler_type))
  
  return(seurat_obj)
}
# =============================================================================
# MALIGNANCY DETECTION
# =============================================================================

detect_malignancy <- function(seurat_obj, tissue_type, output_dir) {
  
  cat("\n=== DETECTING MALIGNANT CELLS ===\n")
  
  # Liver-specific malignancy markers
  markers <- switch(tolower(tissue_type),
                    "liver" = c("KRT19", "EPCAM", "MUC1", "SPP1", "CEACAM5"),
                    c("EPCAM", "KRT19", "MUC1") # Default markers
  )
  
  # Use only markers present in data
  markers <- markers[markers %in% rownames(seurat_obj)]
  
  if (length(markers) == 0) {
    warning("No malignancy markers found in dataset!")
    seurat_obj$malignancy_score <- 0
  } else {
    seurat_obj <- AddModuleScore(seurat_obj,
                                 features = list(markers),
                                 name = "malignancy",
                                 ctrl = min(100, ncol(seurat_obj)/3))
    seurat_obj$malignancy_score <- seurat_obj$malignancy1
    seurat_obj$malignancy1 <- NULL
  }
  
  # Plot
  p <- FeaturePlot(seurat_obj, "malignancy_score", pt.size = 0.5) +
    scale_color_viridis_c() +
    ggtitle("Malignancy Score")
  ggsave(file.path(output_dir, "malignancy_scores.pdf"), p, width = 8, height = 6)
  
  return(seurat_obj)
}
# ===========================================
# IMMUNE SUBTYPING MODULE
# ===========================================
refine_immune_subtypes <- function(seurat_obj, output_dir = getwd(), min_cells = 24) {
  
  # 1. Setup and validation ---------------------------------------------------
  if (!"singler_type" %in% colnames(seurat_obj@meta.data)) {
    stop("Run SingleR annotation first using run_singler_annotation()")
  }
  
  # Create output directory
  output_dir <- normalizePath(output_dir, winslash = "/", mustWork = FALSE)
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  # 2. Detect immune cells more aggressively ----------------------------------
  immune_terms <- c("T cell", "B cell", "NK cell", "Monocyte", "Macrophage",
                    "Dendritic", "Plasma cell", "Granulocyte", "HSC", "Progenitor",
                    "CD4", "CD8", "Treg", "Memory", "Naive", "ILC", "MAIT")
  
  immune_cells <- grepl(paste(immune_terms, collapse = "|"), 
                        seurat_obj$singler_type, ignore.case = TRUE)
  
  if (sum(immune_cells) < min_cells) {
    warning("Only ", sum(immune_cells), " immune cells found (minimum ", min_cells, " required)")
    seurat_obj$immune_subtype <- NA
    return(seurat_obj)
  }
  
  cat("Found", sum(immune_cells), "immune cells for subtyping\n")
  print(table(seurat_obj$singler_type[immune_cells]))
  
  # 3. Enhanced immune cell processing ----------------------------------------
  immune_obj <- subset(seurat_obj, cells = colnames(seurat_obj)[immune_cells])
  
  # Re-process with higher sensitivity
  immune_obj <- NormalizeData(immune_obj, verbose = FALSE)
  immune_obj <- FindVariableFeatures(immune_obj, nfeatures = 3000, verbose = FALSE)
  immune_obj <- ScaleData(immune_obj, verbose = FALSE)
  immune_obj <- RunPCA(immune_obj, npcs = 30, verbose = FALSE)
  
  # Determine optimal clusters
  immune_obj <- FindNeighbors(immune_obj, dims = 1:20, verbose = FALSE)
  immune_obj <- FindClusters(immune_obj, resolution = c(0.2, 0.5, 0.8), verbose = FALSE)
  
  # Use highest resolution that produces reasonable clusters
  if ("RNA_snn_res.0.8" %in% colnames(immune_obj@meta.data)) {
    Idents(immune_obj) <- "RNA_snn_res.0.8"
  } else if ("RNA_snn_res.0.5" %in% colnames(immune_obj@meta.data)) {
    Idents(immune_obj) <- "RNA_snn_res.0.5"
  } else {
    Idents(immune_obj) <- "RNA_snn_res.0.2"
  }
  
  immune_obj <- RunUMAP(immune_obj, dims = 1:20, verbose = FALSE)
  
  # 4. Robust annotation with multiple references -----------------------------
  immune_sce <- as.SingleCellExperiment(immune_obj)
  
  # Try both Monaco and HPCA references
  tryCatch({
    monaco_ref <- MonacoImmuneData()
    pred_monaco <- SingleR(test = immune_sce, ref = monaco_ref, labels = monaco_ref$label.fine)
    immune_obj$monaco_labels <- pred_monaco$labels
  }, error = function(e) warning("Monaco reference failed: ", e$message))
  
  tryCatch({
    hpca_ref <- HumanPrimaryCellAtlasData()
    pred_hpca <- SingleR(test = immune_sce, ref = hpca_ref, labels = hpca_ref$label.main)
    immune_obj$hpca_labels <- pred_hpca$labels
  }, error = function(e) warning("HPCA reference failed: ", e$message))
  
  # Create consensus annotation
  if (all(c("monaco_labels", "hpca_labels") %in% colnames(immune_obj@meta.data))) {
    immune_obj$immune_subtype <- ifelse(
      grepl("^T cell|^B cell|^NK cell", immune_obj$hpca_labels),
      immune_obj$monaco_labels,
      immune_obj$hpca_labels
    )
  } else if ("monaco_labels" %in% colnames(immune_obj@meta.data)) {
    immune_obj$immune_subtype <- immune_obj$monaco_labels
  } else if ("hpca_labels" %in% colnames(immune_obj@meta.data)) {
    immune_obj$immune_subtype <- immune_obj$hpca_labels
  } else {
    stop("Failed to annotate using both references")
  }
  
  # 5. Visualization and output -----------------------------------------------
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  # Plot UMAP with subtypes
  p <- DimPlot(immune_obj, group.by = "immune_subtype", label = TRUE) +
    ggtitle("Immune Cell Subtypes")
  ggsave(file.path(output_dir, "immune_subtypes_umap.pdf"), p, width = 10, height = 8)
  
  # Plot composition
  comp_plot <- immune_obj@meta.data %>%
    count(immune_subtype) %>%
    ggplot(aes(x = reorder(immune_subtype, n), y = n, fill = immune_subtype)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    labs(title = "Immune Cell Composition", x = "", y = "Cell Count") +
    theme_minimal()
  ggsave(file.path(output_dir, "immune_composition.pdf"), comp_plot, width = 8, height = 6)
  
  # 6. Integrate back to main object -----------------------------------------
  seurat_obj$immune_subtype <- NA
  seurat_obj$immune_subtype[colnames(immune_obj)] <- immune_obj$immune_subtype
  
  # Update enhanced annotations
  if (!"celltype_enhanced" %in% colnames(seurat_obj@meta.data)) {
    seurat_obj$celltype_enhanced <- seurat_obj$singler_type
  }
  
  immune_idx <- which(!is.na(seurat_obj$immune_subtype))
  seurat_obj$celltype_enhanced[immune_idx] <- paste(
    seurat_obj$celltype_enhanced[immune_idx],
    paste0("(", seurat_obj$immune_subtype[immune_idx], ")")
  )
  
  cat("\n===== Final Immune Subtype Distribution =====\n")
  print(table(seurat_obj$immune_subtype))
  
  return(seurat_obj)
}
# =============================================================================
# CAF DETECTION 
# =============================================================================
detect_CAFs <- function(seurat_obj, output_dir = getwd()) {
  
  cat("\n=== IDENTIFYING CAFs ===\n")
  
  # Universal CAF markers
  caf_markers <- c("ACTA2", "PDGFRA", "PDGFRB", "COL1A1", "COL3A1", "COL1A2", "C1R", "SERPINF1", "C1S")
  caf_markers <- caf_markers[caf_markers %in% rownames(seurat_obj)]
  
  if (length(caf_markers) == 0) {
    warning("No CAF markers found in dataset!")
    seurat_obj$CAF_score <- 0
    seurat_obj$is_CAF <- "Non_CAF"
  } else {
    seurat_obj <- AddModuleScore(
      seurat_obj,
      features = list(caf_markers),
      name = "CAF",
      ctrl = min(100, ncol(seurat_obj)/3)
    )
    seurat_obj$CAF_score <- seurat_obj$CAF1
    seurat_obj$CAF1 <- NULL
    
    # Identify top 10% as CAFs
    caf_threshold <- quantile(seurat_obj$CAF_score, 0.9)
    seurat_obj$is_CAF <- ifelse(seurat_obj$CAF_score > caf_threshold, "CAF", "Non_CAF")
    
    # Create plots directory if it doesn't exist
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    
    # Plot 1: CAF Scores (FeaturePlot)
    p1 <- FeaturePlot(seurat_obj, "CAF_score", pt.size = 0.5) +
      scale_color_viridis_c(option = "viridis") +
      ggtitle("CAF Signature Score") +
      theme(plot.title = element_text(face = "bold", size = 14))
    
    ggsave(file.path(output_dir, "caf_scores.pdf"), p1, width = 8, height = 6)
    ggsave(file.path(output_dir, "caf_scores.png"), p1, width = 8, height = 6, dpi = 300)
    
    # Plot 2: CAF Identification (DimPlot)
    p2 <- DimPlot(seurat_obj, 
                  group.by = "is_CAF", 
                  cols = c("CAF" = "red", "Non_CAF" = "gray90"),
                  pt.size = 0.5) +
      ggtitle("CAF Identification") +
      theme(plot.title = element_text(face = "bold", size = 14))
    
    ggsave(file.path(output_dir, "caf_identification.pdf"), p2, width = 8, height = 6)
    ggsave(file.path(output_dir, "caf_identification.png"), p2, width = 8, height = 6, dpi = 300)
    
    # Print summary
    cat("CAF detection results:\n")
    cat("-", sum(seurat_obj$is_CAF == "CAF"), "cells identified as CAFs\n")
    cat("-", sum(seurat_obj$is_CAF == "Non_CAF"), "cells as Non-CAF\n")
    cat("- Threshold score:", round(caf_threshold, 3), "\n")
  }
  
  return(seurat_obj)
}
# =============================================================================
# FINAL ANNOTATION INTEGRATION (UPDATED WITH IMMUNE SUBTYPES)
# =============================================================================

create_final_annotations <- function(seurat_obj, output_dir) {
  
  cat("\n=== CREATING FINAL ANNOTATIONS ===\n")
  
  # ---------------------------------------------------------------------------
  # 1. CREATE BASE ANNOTATIONS
  # ---------------------------------------------------------------------------
  
  # Start with SingleR annotations
  seurat_obj$celltype_final <- seurat_obj$singler_type
  
  # Override with malignancy calls (top 5% scores)
  malignant_cells <- seurat_obj$malignancy_score > quantile(seurat_obj$malignancy_score, 0.95, na.rm = TRUE)
  seurat_obj$celltype_final[malignant_cells] <- "Malignant"
  
  # Override with CAF calls
  seurat_obj$celltype_final[seurat_obj$is_CAF == "CAF"] <- "CAF"
  
  # ---------------------------------------------------------------------------
  # 2. INTEGRATE IMMUNE SUBTYPES (IF AVAILABLE)
  # ---------------------------------------------------------------------------
  
  if ("immune_subtype" %in% colnames(seurat_obj@meta.data)) {
    # Only update if immune_subtype exists and isn't NA
    immune_cells <- !is.na(seurat_obj$immune_subtype)
    
    if (any(immune_cells)) {
      # Preserve the broad category (e.g., "T cells") while adding subtype info
      seurat_obj$celltype_final[immune_cells] <- paste(
        seurat_obj$celltype_final[immune_cells],
        paste0("(", seurat_obj$immune_subtype[immune_cells], ")")
      )
      
      # Clean up any cases where we might have double parentheses
      seurat_obj$celltype_final <- gsub("\\(\\(", "(", seurat_obj$celltype_final)
    }
  }
  
  # ---------------------------------------------------------------------------
  # 3. SAVE COMPREHENSIVE RESULTS
  # ---------------------------------------------------------------------------
  
  # Prepare metadata with all annotation layers
  final_meta <- data.frame(
    Cell = colnames(seurat_obj),
    Cluster = seurat_obj$seurat_clusters,
    SingleR_Type = seurat_obj$singler_type,
    SingleR_Score = seurat_obj$singler_score,
    Malignancy_Score = seurat_obj$malignancy_score,
    CAF_Score = seurat_obj$CAF_score,
    Immune_Subtype = if ("immune_subtype" %in% colnames(seurat_obj@meta.data)) seurat_obj$immune_subtype else NA,
    Final_Type = seurat_obj$celltype_final,
    stringsAsFactors = FALSE
  )
  
  # Save with error handling for paths
  tryCatch({
    write.csv(final_meta, file.path(output_dir, "final_annotations.csv"), row.names = FALSE)
    cat("Saved final annotations to:", file.path(output_dir, "final_annotations.csv"), "\n")
  }, error = function(e) {
    warning("Failed to save annotations: ", e$message)
  })
  
  # ---------------------------------------------------------------------------
  # 4. VISUALIZATION
  # ---------------------------------------------------------------------------
  
  # Create color palette (consistent colors for same cell types across plots)
  celltypes <- unique(seurat_obj$celltype_final)
  color_palette <- setNames(
    viridis::viridis(length(celltypes)),
    celltypes
  )
  
  p <- DimPlot(seurat_obj, 
               group.by = "celltype_final",
               cols = color_palette,
               label = TRUE, 
               repel = TRUE,
               label.size = 3,
               pt.size = 0.7) +
    ggtitle("Final Cell Type Annotations") +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold", size = 14))
  
  # Save plot with path handling
  tryCatch({
    ggsave(file.path(output_dir, "final_annotations.pdf"), p, width = 12, height = 8)
    cat("Saved final annotation plot to:", file.path(output_dir, "final_annotations.pdf"), "\n")
  }, error = function(e) {
    warning("Failed to save final plot: ", e$message)
  })
  
  # ---------------------------------------------------------------------------
  # 5. SUMMARY STATISTICS
  # ---------------------------------------------------------------------------
  
  cat("\nFinal cell type distribution:\n")
  print(sort(table(seurat_obj$celltype_final), decreasing = TRUE))
  
  if ("immune_subtype" %in% colnames(seurat_obj@meta.data)) {
    cat("\nImmune subtype distribution:\n")
    print(sort(table(seurat_obj$immune_subtype), decreasing = TRUE))
  }
  
  return(seurat_obj)
}
# =============================================================================
# ANALYSIS OF GENE EXPRESSION WITH VIRIDIS GRADIENT (SAVES ALL PLOTS & TABLES)
# =============================================================================

get_user_input <- function(prompt) {
  if (interactive()) {
    return(readline(prompt))
  } else {
    cat(prompt)
    return(readLines("stdin", n = 1))
  }
}

analyze_gene_expression <- function(seurat_obj, gene = NULL, gene_threshold = NULL, output_dir = getwd()) {
  # Initialize empty results list
  results <- list()
  
  # 1. Prompt user for gene if not provided
  if (is.null(gene) || gene == "") {
    gene <- get_user_input("Enter the gene symbol you want to analyze (e.g., FAP): ")
    if (gene == "") {
      gene <- "FAP"
      cat("No gene entered. Defaulting to FAP.\n")
    }
  }
  
  cat("Analyzing", gene, "expression...\n")
  
  # Check if gene exists (case-insensitive)
  gene_upper <- toupper(gene)
  matches <- which(toupper(rownames(seurat_obj)) == gene_upper)
  if (length(matches) == 0) {
    warning(paste("Gene ('", gene, "') not found in the dataset. Showing empty expression plot.", sep = ""))
    
    # Create empty UMAP plot to show absence of gene
    p_empty <- FeaturePlot(
      seurat_obj,
      features = sample(rownames(seurat_obj), 1),  # Use random gene just to show the UMAP
      reduction = "umap",
      pt.size = 0.5,
      cols = c("grey90", "grey90")  # All cells grey
    ) +
      ggtitle(paste(gene, "not found in dataset")) +
      theme_minimal()
    
    print(p_empty)
    pdf(file.path(output_dir, paste0("umap_", gene, "_not_found.pdf")), width = 8, height = 6)
    print(p_empty)
    dev.off()
    
    return(list(
      gene_found = FALSE,
      umap_plot = p_empty
    ))
  }
  
  # Get correct case of gene
  gene <- rownames(seurat_obj)[matches[1]]
  if (gene != gene_upper) {
    cat(paste("Note: Using case-matched gene", gene, "\n"))
  }
  
  # 2. Get gene expression vector (normalized data slot)
  gene_expr <- as.numeric(GetAssayData(seurat_obj, slot = "data")[gene, ])
  
  # 3. Show summary statistics to user
  cat("\nSummary statistics for", gene, "expression:\n")
  cat("Mean:    ", mean(gene_expr), "\n")
  cat("Median:  ", median(gene_expr), "\n")
  cat("Min:     ", min(gene_expr), "\n")
  cat("Max:     ", max(gene_expr), "\n")
  cat("Quantiles:\n")
  print(quantile(gene_expr))
  cat("\n")
  
  # 4. Prompt user for threshold, using summary as guidance
  if (is.null(gene_threshold) || is.na(gene_threshold)) {
    threshold_input <- get_user_input(
      paste0("Enter the threshold for ", gene, " positivity (e.g., 0.0005): ")
    )
    if (threshold_input == "") {
      gene_threshold <- 0.0005
      cat("No threshold entered. Defaulting to 0.0005\n")
    } else {
      gene_threshold <- as.numeric(threshold_input)
      if (is.na(gene_threshold)) {
        gene_threshold <- 0.0005
        cat("Invalid input. Defaulting to 0.0005\n")
      }
    }
  }
  cat("Using threshold:", gene_threshold, "\n")
  
  # 5. Store expression and positivity in metadata
  seurat_obj[[paste0(gene, "_expression")]] <- gene_expr
  seurat_obj[[paste0(gene, "_positive")]] <- gene_expr > gene_threshold
  
  # 6. Average gene expression per cluster
  avg_gene_expression <- AggregateExpression(seurat_obj,
                                             features = gene,
                                             return.seurat = FALSE)
  avg_gene_expression <- avg_gene_expression[[gene]]
  cat("Average", gene, "expression per cluster:\n")
  print(avg_gene_expression)
  
  # --- Calculate and save gene+ percentages by cell type/subtype ---
  library(dplyr)
  
  # Helper function to calculate percentages by group
  percent_table <- function(meta, group_col, gene_col) {
    meta %>%
      group_by(.data[[group_col]]) %>%
      summarise(
        n_total = n(),
        n_positive = sum(.data[[gene_col]], na.rm = TRUE),
        percent_positive = round(100 * n_positive / n_total, 2)
      ) %>%
      arrange(desc(percent_positive))
  }
  
  meta <- seurat_obj@meta.data
  gene_pos_col <- paste0(gene, "_positive")
  
  # a) By cell subtype
  subtype_table <- NULL
  if ("celltype_enhanced" %in% colnames(meta)) {
    subtype_table <- percent_table(meta, "celltype_enhanced", gene_pos_col)
    print(subtype_table)
    write.csv(subtype_table, file = file.path(output_dir, paste0(gene, "_positive_by_celltype.csv")), row.names = FALSE)
  }
  
  # b) By CAF subtype
  caf_table <- NULL
  if ("caf_subtype" %in% colnames(meta)) {
    caf_table <- percent_table(meta, "caf_subtype", gene_pos_col)
    print(caf_table)
    write.csv(caf_table, file = file.path(output_dir, paste0(gene, "_positive_by_CAF_subtype.csv")), row.names = FALSE)
  }
  
  # c) By immune cell subtype
  immune_table <- NULL
  if ("immune_subtype" %in% colnames(meta)) {
    immune_table <- percent_table(meta, "immune_subtype", gene_pos_col)
    print(immune_table)
    write.csv(immune_table, file = file.path(output_dir, paste0(gene, "_positive_by_immune_subtype.csv")), row.names = FALSE)
  }
  
  # 7. Violin plot (no dots, stylised)
  vln_plot <- VlnPlot(
    seurat_obj,
    features = gene,
    pt.size = 0,
    group.by = "celltype_enhanced"
  ) +
    ggtitle(paste(gene, "Expression per Annotated Cell Type")) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_line(color = "grey90"),
      plot.title = element_text(face = "bold", hjust = 0.5)
    ) +
    scale_fill_viridis_d(option = "C") +
    ylab("Normalized Expression") +
    xlab("Cell Type")
  
  print(vln_plot)
  pdf(file.path(output_dir, paste0("vlnplot_", gene, ".pdf")), width = 8, height = 6)
  print(vln_plot)
  dev.off()
  
  # 8. UMAP plot: continuous gene expression (viridis gradient)
  viridis_pal <- viridis::viridis(100)
  p_umap_gene_gradient <- FeaturePlot(
    seurat_obj,
    features = gene,
    reduction = "umap",
    pt.size = 0.5,
    cols = viridis_pal
  ) +
    ggtitle(paste(gene, "Expression (UMAP, viridis gradient)")) +
    theme_minimal()
  print(p_umap_gene_gradient)
  pdf(file.path(output_dir, paste0("umap_", gene, "_gradient.pdf")), width = 8, height = 6)
  print(p_umap_gene_gradient)
  dev.off()
  
  # 9. UMAP plot: gene positive vs negative (binary)
  p_umap_gene_binary <- DimPlot(seurat_obj, group.by = paste0(gene, "_positive"), reduction = "umap", pt.size = 0.5) +
    scale_color_manual(values = c("FALSE" = "grey80", "TRUE" = "red")) +
    ggtitle(paste0(gene, "+ Cells (UMAP, threshold > ", gene_threshold, ")")) +
    theme_minimal()
  print(p_umap_gene_binary)
  pdf(file.path(output_dir, paste0("umap_", gene, "_positive.pdf")), width = 8, height = 6)
  print(p_umap_gene_binary)
  dev.off()
  
  cat("✓", gene, "expression analysis complete\n")
  
  return(list(
    gene_found = TRUE,
    avg_expression = avg_gene_expression,
    vln_plot = vln_plot,
    umap_gene_gradient = p_umap_gene_gradient,
    umap_gene_binary = p_umap_gene_binary,
    subtype_table = subtype_table,
    caf_table = caf_table,
    immune_table = immune_table
  ))
}

# =============================================================================
# MAIN PIPELINE EXECUTION (MODULAR, INTERACTIVE, AND COMPLETE)
# =============================================================================

# 1. Install and load required packages
install_required_packages()
load_required_packages()

# 2. Choose and set up the workspace/project directory
project_dir <- select_workspace_directory()

# Optionally, create output subdirectories (plots, results, etc.)
plots_dir <- file.path(project_dir, "plots")
if (!dir.exists(plots_dir)) {
  dir.create(plots_dir, recursive = TRUE)
}

# 3. Detect available data files in the workspace directory
data_files <- list.files(project_dir, full.names = TRUE)
if (length(data_files) == 0) {
  stop("No data files found. Please add data to the selected folder.")
}

# 4. Automated data file selection: prefer .rds, then .RData, then .h5, then matrix.mtx.gz
rds_file <- data_files[grepl("\\.rds$", data_files, ignore.case = TRUE)]
rdata_file <- data_files[grepl("\\.RData$", data_files, ignore.case = TRUE)]
h5_file <- data_files[grepl("\\.h5$", data_files, ignore.case = TRUE)]
mtx_file <- data_files[grepl("matrix.*\\.mtx\\.gz$", data_files, ignore.case = TRUE)]

if (length(rds_file) > 0) {
  data_path <- rds_file[1]
  cat("Auto-selected Seurat .rds file:", data_path, "\n")
} else if (length(rdata_file) > 0) {
  data_path <- rdata_file[1]
  cat("Auto-selected Seurat .RData file:", data_path, "\n")
} else if (length(h5_file) > 0) {
  data_path <- h5_file[1]
  cat("Auto-selected HDF5 file:", data_path, "\n")
} else if (length(mtx_file) > 0) {
  data_path <- mtx_file[1]
  cat("Auto-selected matrix.mtx.gz file:", data_path, "\n")
} else {
  stop("No suitable data file (.rds, .RData, .h5 or matrix.mtx.gz) found in the selected folder.")
}

# 5. Load data and create Seurat object
seurat_obj <- load_data(data_path)

# NEW: Interactive step selection
if (inherits(seurat_obj, "Seurat")) {
  cat("\nLoaded Seurat object with", ncol(seurat_obj), "cells and", nrow(seurat_obj), "genes\n")
  
  # Define pipeline steps with descriptions
  pipeline_steps <- list(
    "6" = "Calculate QC metrics",
    "7" = "Visualize QC metrics",
    "8" = "Filter cells",
    "9" = "Normalize data",
    "10" = "Dimensionality reduction",
    "11" = "Visualize clusters",
    "12" = "Find clusters",
    "13" = "Enhanced cell type annotation",
    "14" = "Visualize enhanced annotation",
    "15" = "Gene expression analysis",
    "16" = "Save processed object"
  )
  
  # Display available steps
  cat("\nPipeline steps available:\n")
  for (step in names(pipeline_steps)) {
    cat(step, ":", pipeline_steps[[step]], "\n")
  }
  
  # Get user input for starting step
  start_step <- readline(prompt = "Enter the step number to start from (or press Enter to run full pipeline): ")
  
  # If user provides input, validate and set start point
  if (start_step != "") {
    if (!start_step %in% names(pipeline_steps)) {
      stop("Invalid step number. Please enter a number between 6 and 16.")
    }
    cat("Starting pipeline from step", start_step, ":", pipeline_steps[[step]], "\n")
  } else {
    cat("Running full pipeline from the beginning.\n")
    start_step <- "6"  # Default to first analysis step
  }
} else {
  stop("Loaded object is not a Seurat object. Please check your input data.")
}

# Function to execute pipeline from selected step
execute_pipeline <- function(seurat_obj, start_step) {
  # Convert start_step to numeric for comparison
  start_num <- as.numeric(start_step)
  
  # Step 6: Calculate QC metrics
  if (start_num <= 6) {
    cat("\n=== Step 6: Calculating QC metrics ===\n")
    seurat_obj <- calculate_qc_metrics(seurat_obj)
  }
  
  # Step 7: Visualize QC metrics
  if (start_num <= 7) {
    cat("\n=== Step 7: Visualizing QC metrics ===\n")
    visualize_qc_metrics(seurat_obj, output_dir = plots_dir)
  }
  
  # Step 8: Filter cells
  if (start_num <= 8) {
    cat("\n=== Step 8: Filtering cells ===\n")
    seurat_obj <- filter_cells(seurat_obj)
  }
  
  # Step 9: Normalize data
  if (start_num <= 9) {
    cat("\n=== Step 9: Normalizing data ===\n")
    seurat_obj <- normalize_data(seurat_obj)
  }
  
  # Step 10: Dimensionality reduction
  if (start_num <= 10) {
    cat("\n=== Step 10: Performing dimensionality reduction ===\n")
    dr_results <- perform_dimensionality_reduction(seurat_obj)
    seurat_obj <- dr_results$seurat_obj
  }
  
  # Step 11: Visualize clusters
  if (start_num <= 11) {
    cat("\n=== Step 11: Visualizing clusters ===\n")
    visualize_clusters(seurat_obj, output_dir = plots_dir)
  }
  
  # Step 12: Find clusters
  if (start_num <= 12) {
    cat("\n=== Step 12: Finding clusters ===\n")
    seurat_obj <- find_clusters(seurat_obj)
    # Preliminary save after finding clusters
    saveRDS(seurat_obj, file = file.path(project_dir, "seurat_obj_processed.rds"))
  }
  
  # Step 13: Enhanced cell type annotation
  if (start_num <= 13) {
    cat("\n=== Step 13: Performing cell type annotation ===\n")
    seurat_obj <- complete_cell_typing(
      seurat_obj = seurat_obj,
      tissue_type = "Liver", 
      output_dir = file.path(plots_dir, "cell_typing"),
      min_confidence = 0.4
    )
    # Save intermediate results
    saveRDS(seurat_obj, file = file.path(project_dir, "seurat_obj_annotated.rds"))
  }
  
  if (start_num <= 14) {
    cat("\n=== Step 14: Visualizing cell type annotation ===\n")
    
    # Create a subdirectory for annotation plots
    annotation_plot_dir <- file.path(plots_dir, "annotations")
    if (!dir.exists(annotation_plot_dir)) dir.create(annotation_plot_dir)
    
    # Final UMAP - use enhanced annotations if available
    group_by_col <- ifelse("celltype_enhanced" %in% colnames(seurat_obj@meta.data),
                           "celltype_enhanced", "celltype_final")
    
    p_final <- DimPlot(seurat_obj, 
                       group.by = group_by_col,
                       label = TRUE,
                       repel = TRUE,
                       label.size = 5) +
      ggtitle("Final Cell Type Annotation") +
      theme_minimal()
    
    ggsave(file.path(annotation_plot_dir, "final_annotations.pdf"), p_final, width = 12, height = 8)
    print(p_final)
  }
  
  # Step 15: Gene expression analysis
  if (start_num <= 15) {
    repeat {
      cat("\n=== Step 15: Performing gene expression analysis ===\n")
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
  
  # Step 16: Save processed Seurat object
  if (start_num <= 16) {
    cat("\n=== Step 16: Saving processed object ===\n")
    saveRDS(seurat_obj, file = file.path(project_dir, "seurat_obj_processed.rds"))
    cat("Saved processed Seurat object to:", file.path(project_dir, "seurat_obj_processed.rds"), "\n")
  }
  
  return(seurat_obj)
}

# Execute the pipeline from the selected step
seurat_obj <- execute_pipeline(seurat_obj, start_step)

cat("\n=== Pipeline execution complete ===\n")