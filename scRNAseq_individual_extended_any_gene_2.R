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
    h5_files <- list.files(data_path, pattern = "\\.h5$", full.names = TRUE, ignore.case = TRUE)
    mtx_files <- list.files(data_path, pattern = "^.*matrix.*\\.mtx\\.gz$", full.names = TRUE, ignore.case = TRUE)
    all_files <- c(h5_files, mtx_files)
    all_files <- all_files[file.info(all_files)$isdir == FALSE]
    if (length(all_files) == 0) stop("No .h5 or matrix .mtx.gz data files found in the directory.")
    data_path <- all_files[1]
    cat("Auto-selected data file:", data_path, "\n")
  }
  
  if (!file.exists(data_path)) stop("Selected file does not exist: ", data_path)
  
  if (grepl("\\.h5$", data_path, ignore.case = TRUE)) {
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
# ENHANCED CELL TYPE ANNOTATION WITH MALIGNANCY & CAF DETECTION
# =============================================================================

# CORE SC-TYPE FUNCTIONS
setup_sctype <- function() {
  cat("Setting up scType functions...\n")
  
  # Install and load HGNChelper for checkGeneSymbols
  if (!requireNamespace("HGNChelper", quietly = TRUE)) install.packages("HGNChelper")
  suppressPackageStartupMessages(library(HGNChelper))
  
  # Load AnnotationDbi and org.Hs.eg.db
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  if (!requireNamespace("AnnotationDbi", quietly = TRUE)) BiocManager::install("AnnotationDbi")
  if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) BiocManager::install("org.Hs.eg.db")
  
  suppressPackageStartupMessages({
    library(AnnotationDbi)
    library(org.Hs.eg.db)
  })
  
  sctype_dir <- file.path(tempdir(), "sctype_scripts")
  dir.create(sctype_dir, showWarnings = FALSE, recursive = TRUE)
  
  sctype_scripts <- c(
    "gene_sets_prepare.R",
    "sctype_score_.R"
  )
  
  base_url <- "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/"
  
  for (script in sctype_scripts) {
    script_path <- file.path(sctype_dir, script)
    if (!file.exists(script_path)) {
      download.file(
        url = paste0(base_url, script),
        destfile = script_path,
        mode = "wb"
      )
    }
    source(script_path)
    cat(paste("✓ Sourced:", script, "\n"))
  }
  
  marker_file <- file.path(sctype_dir, "ScTypeDB_full.xlsx")
  if (!file.exists(marker_file)) {
    download.file(
      "https://github.com/IanevskiAleksandr/sc-type/raw/master/ScTypeDB_full.xlsx",
      destfile = marker_file,
      mode = "wb"
    )
    cat("✓ Downloaded marker database\n")
  }
  
  return(list(
    scripts_dir = sctype_dir,
    marker_file = marker_file
  ))
}

annotate_cell_types <- function(seurat_obj, tissue = "Liver", output_dir = getwd()) {
  cat("Performing automated cell type annotation with scType...\n")
  
  sctype_resources <- setup_sctype()
  marker_file <- sctype_resources$marker_file
  
  gs_list <- gene_sets_prepare(marker_file, tissue)
  
  if (length(gs_list) == 0) {
    stop(paste("No markers found for tissue:", tissue))
  }
  
  es.max <- sctype_score(
    scRNAseqData = as.matrix(GetAssayData(seurat_obj, slot = "data")), 
    scaled = TRUE,
    gs = gs_list$gs_positive,
    gs2 = gs_list$gs_negative
  )
  
  clusters <- seurat_obj$seurat_clusters
  
  cL_resutls <- do.call(
    "rbind", 
    lapply(unique(clusters), function(cl) {
      es.max.cl <- sort(rowSums(es.max[, names(clusters)[clusters == cl], drop = FALSE]), decreasing = TRUE)
      head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl), 10)
    })
  )
  
  cL_resutls_top <- cL_resutls %>%
    dplyr::group_by(cluster) %>%
    dplyr::top_n(1, scores)
  
  cluster_ids <- as.character(seurat_obj$seurat_clusters)
  celltype_vec <- setNames(as.character(cL_resutls_top$type), as.character(cL_resutls_top$cluster))
  celltype_auto <- celltype_vec[cluster_ids]
  names(celltype_auto) <- NULL
  seurat_obj$celltype_auto <- celltype_auto
  
  cat("✓ Cell type annotation with scType complete\n")
  
  # UMAP plot
  p_umap_celltypes <- DimPlot(seurat_obj, group.by = "celltype_auto", reduction = "umap", 
                              label = TRUE, repel = TRUE, label.size = 3) +
    ggtitle("Cell Types (scType)") + 
    theme_minimal() +
    theme(legend.position = "right")
  
  print(p_umap_celltypes)
  ggsave(file.path(output_dir, "scType_annotation.pdf"), p_umap_celltypes, width = 10, height = 8)
  
  return(seurat_obj)
}

# MALIGNANCY & CAF DETECTION FUNCTIONS

detect_malignancy <- function(seurat_obj, tissue_type = "Liver", output_dir = getwd()) {
  cat("\nPerforming malignancy detection...\n")
  
  # Define tissue-specific cancer markers
  marker_db <- list(
    "Liver" = c("AFP", "GPC3", "SPP1", "MAGEB2", "AKR1B10"),
    "Pancreas" = c("MSLN", "CEACAM5", "TFF1", "CLDN18", "MUC5AC"),
    "Breast" = c("ESR1", "ERBB2", "KRT5", "KRT17", "MUC1"),
    "Lung" = c("NKX2-1", "NAPSA", "SFTPB", "CEACAM5"),
    "Immune system" = c("CD274", "PDCD1LG2", "CTLA4", "IDO1")  # Immune checkpoint markers
  )
  
  # Select markers or use default
  if (tissue_type %in% names(marker_db)) {
    malignancy_markers <- marker_db[[tissue_type]]
    cat("Using tissue-specific markers for:", tissue_type, "\n")
  } else {
    warning("Using generic epithelial cancer markers. Consider adding tissue-specific markers.")
    malignancy_markers <- c("EPCAM", "KRT19", "MUC1", "CEACAM5")
  }
  
  # Filter to markers present in dataset
  valid_markers <- malignancy_markers[malignancy_markers %in% rownames(seurat_obj)]
  cat("Using malignancy markers:", paste(valid_markers, collapse = ", "), "\n")
  
  if (length(valid_markers) == 0) {
    warning("No valid malignancy markers found in dataset!")
    seurat_obj$malignancy_score <- 0
  } else {
    # Calculate malignancy score
    seurat_obj <- AddModuleScore(
      object = seurat_obj,
      features = list(valid_markers),
      name = "malignancy",
      ctrl = min(100, floor(ncol(seurat_obj)/3))
    )
    seurat_obj$malignancy_score <- seurat_obj$malignancy1
    seurat_obj$malignancy1 <- NULL  # Clean up temporary column
  }
  
  # Visualize
  p_malignancy <- FeaturePlot(seurat_obj, "malignancy_score", reduction = "umap", pt.size = 0.7) +
    scale_color_viridis_c(option = "magma", direction = -1) + 
    ggtitle("Malignancy Score") +
    theme(plot.title = element_text(face = "bold", size = 14))
  
  ggsave(file.path(output_dir, "malignancy_score.pdf"), p_malignancy, width = 8, height = 6)
  print(p_malignancy)
  
  return(seurat_obj)
}

detect_CAFs <- function(seurat_obj, output_dir = getwd()) {
  cat("\nDetecting Cancer-Associated Fibroblasts (CAFs)...\n")
  
  # Comprehensive CAF markers
  caf_markers <- c("ACTA2", "FAP", "PDGFRA", "PDGFRB", "THY1", "TAGLN", 
                   "COL1A1", "COL1A2", "COL3A1", "MMP2", "S100A4")
  
  valid_markers <- caf_markers[caf_markers %in% rownames(seurat_obj)]
  cat("Using CAF markers:", paste(valid_markers, collapse = ", "), "\n")
  
  if (length(valid_markers) == 0) {
    warning("No valid CAF markers found in dataset!")
    seurat_obj$CAF_score <- 0
    seurat_obj$is_CAF <- "Non-CAF"
  } else {
    # Calculate CAF score
    seurat_obj <- AddModuleScore(
      object = seurat_obj,
      features = list(valid_markers),
      name = "CAF",
      ctrl = min(100, floor(ncol(seurat_obj)/3))
    )
    
    seurat_obj$CAF_score <- seurat_obj$CAF1
    seurat_obj$CAF1 <- NULL  # Clean up temporary column
    
    # Identify CAF-enriched cells (top 10% scoring cells)
    caf_threshold <- quantile(seurat_obj$CAF_score, 0.90)
    seurat_obj$is_CAF <- ifelse(seurat_obj$CAF_score > caf_threshold, "CAF", "Non-CAF")
  }
  
  # Visualize
  p_caf <- FeaturePlot(seurat_obj, "CAF_score", reduction = "umap", pt.size = 0.7) +
    scale_color_viridis_c(option = "viridis") + 
    ggtitle("CAF Score") +
    theme(plot.title = element_text(face = "bold", size = 14))
  
  p_caf_type <- DimPlot(seurat_obj, group.by = "is_CAF", reduction = "umap", pt.size = 0.7) +
    scale_color_manual(values = c("CAF" = "darkorange", "Non-CAF" = "lightgray")) +
    ggtitle("CAF Identification") +
    theme(plot.title = element_text(face = "bold", size = 14))
  
  ggsave(file.path(output_dir, "caf_score.pdf"), p_caf, width = 8, height = 6)
  ggsave(file.path(output_dir, "caf_identification.pdf"), p_caf_type, width = 8, height = 6)
  
  print(p_caf)
  print(p_caf_type)
  
  return(seurat_obj)
}

# MASTER ANNOTATION FUNCTION

enhanced_annotation <- function(seurat_obj, tissue = "Liver", output_dir = getwd()) {
  # Create output directory if needed
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  # Step 1: scType annotation
  cat("\n===== STEP 1: BASIC CELL TYPE ANNOTATION =====\n")
  seurat_obj <- annotate_cell_types(seurat_obj, tissue, output_dir)
  
  # Step 2: Malignancy detection
  cat("\n===== STEP 2: MALIGNANCY DETECTION =====\n")
  seurat_obj <- detect_malignancy(seurat_obj, tissue, output_dir)
  
  # Step 3: CAF detection
  cat("\n===== STEP 3: CAF DETECTION =====\n")
  seurat_obj <- detect_CAFs(seurat_obj, output_dir)
  
  # Step 4: Integrate annotations
  cat("\n===== STEP 4: INTEGRATING ANNOTATIONS =====\n")
  seurat_obj$celltype_enhanced <- seurat_obj$celltype_auto
  
  # Override with malignancy calls (top 5% malignant cells)
  if ("malignancy_score" %in% colnames(seurat_obj@meta.data)) {
    malignancy_threshold <- quantile(seurat_obj$malignancy_score, 0.95)
    malignant_cells <- seurat_obj$malignancy_score > malignancy_threshold
    seurat_obj$celltype_enhanced[malignant_cells] <- "Malignant"
    cat(paste("✓ Identified", sum(malignant_cells), "malignant cells\n"))
  }
  
  # Override with CAF calls (preserve malignancy status)
  if ("is_CAF" %in% colnames(seurat_obj@meta.data)) {
    caf_cells <- seurat_obj$is_CAF == "CAF" & seurat_obj$celltype_enhanced != "Malignant"
    seurat_obj$celltype_enhanced[caf_cells] <- "CAF"
    cat(paste("✓ Identified", sum(caf_cells), "CAF cells\n"))
  }
  
  # Step 5: Visualize final results
  p_enhanced <- DimPlot(seurat_obj, group.by = "celltype_enhanced", 
                        reduction = "umap", label = TRUE, repel = TRUE, 
                        pt.size = 0.7, label.size = 4) +
    ggtitle("Enhanced Cell Type Annotation") +
    theme_minimal() +
    theme(plot.title = element_text(face = "bold", size = 16))
  
  ggsave(file.path(output_dir, "enhanced_annotation.pdf"), p_enhanced, width = 12, height = 8)
  print(p_enhanced)
  
  # Save metadata
  write.csv(seurat_obj@meta.data, file.path(output_dir, "cell_metadata.csv"))
  
  cat("\n===== ANNOTATION COMPLETE =====\n")
  cat(paste("Results saved to:", output_dir, "\n"))
  
  return(seurat_obj)
}
# =============================================================================
# ANALYSIS OF GENE EXPRESSION WITH VIRIDIS GRADIENT (SAVES ALL PLOTS)
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
  # 1. Prompt user for gene if not provided
  if (is.null(gene) || gene == "") {
    gene <- get_user_input("Enter the gene symbol you want to analyze (e.g., FAP): ")
    if (gene == "") {
      gene <- "FAP"
      cat("No gene entered. Defaulting to FAP.\n")
    }
  }
  cat("Analyzing", gene, "expression...\n")
  if (!(gene %in% rownames(seurat_obj))) {
    stop(paste("Gene ('", gene, "') not found in the dataset.", sep = ""))
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
  
  # 7. Violin plot
  vln_plot <- VlnPlot(
    seurat_obj,
    features = gene,
    pt.size = 0.5,
    group.by = "celltype_enhanced"
  ) +
    ggtitle(paste(gene, "Expression per Annotated Cell Type")) +
    theme_minimal()
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
    avg_expression = avg_gene_expression,
    vln_plot = vln_plot,
    umap_gene_gradient = p_umap_gene_gradient,
    umap_gene_binary = p_umap_gene_binary
  ))
}
# =============================================================================
# MAIN PIPELINE EXECUTION (MODULAR, INTERACTIVE, AND COMPLETE)
# =============================================================================

# 1. Install and load required packages (ensure these functions are defined elsewhere)
install_required_packages()
load_required_packages()

# 2. Choose and set up the workspace/project directory
project_dir <- select_workspace_directory()

# Optionally, create output subdirectories (plots, results, etc.)
plots_dir <- file.path(project_dir, "plots")
if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)

# 3. Detect available data files in the workspace directory
data_files <- list.files(project_dir, full.names = TRUE)
if (length(data_files) == 0) stop("No data files found. Please add data to the selected folder.")

# 4. Automated data file selection: prefer .h5, then matrix.mtx.gz
h5_file <- data_files[grepl("\\.h5$", data_files, ignore.case = TRUE)]
mtx_file <- data_files[grepl("matrix.*\\.mtx\\.gz$", data_files, ignore.case = TRUE)]

if (length(h5_file) > 0) {
  data_path <- h5_file[1]
  cat("Auto-selected HDF5 file:", data_path, "\n")
} else if (length(mtx_file) > 0) {
  data_path <- mtx_file[1]
  cat("Auto-selected matrix.mtx.gz file:", data_path, "\n")
} else {
  stop("No suitable data file (.h5 or matrix.mtx.gz) found in the selected folder.")
}

# 5. Load data and create Seurat object
seurat_obj <- load_data(data_path)

# 6. Calculate QC metrics (custom function)
seurat_obj <- calculate_qc_metrics(seurat_obj)

# 7. Visualize QC metrics (custom function, saves plots)
visualize_qc_metrics(seurat_obj, output_dir = plots_dir)

# 8. Filter cells (custom function)
seurat_obj <- filter_cells(seurat_obj)

# 9. Normalize data (custom function)
seurat_obj <- normalize_data(seurat_obj)

# 10. Dimensionality reduction (custom function, returns list)
dr_results <- perform_dimensionality_reduction(seurat_obj)
seurat_obj <- dr_results$seurat_obj

# 11. Visualize clusters (custom function, saves plots)
visualize_clusters(seurat_obj, output_dir = plots_dir)

# 12. Find clusters (custom function)
seurat_obj <- find_clusters(seurat_obj)

# 13. Enhanced cell type annotation (modular, with scType + CAF/malignant logic)
# This step should create a metadata column, e.g., 'celltype_enhanced'
seurat_obj <- enhanced_annotation(seurat_obj, tissue = "Liver", output_dir = plots_dir) # Change tissue as needed

# 14. Visualize enhanced annotation (CAFs, Malignant, etc.) and save as PDF
p_umap_enhanced <- DimPlot(
  seurat_obj,
  group.by = "celltype_enhanced",   # Make sure this matches your enhanced annotation column name
  reduction = "umap",
  label = TRUE,
  repel = TRUE,
  label.size = 5
) + ggtitle("Enhanced Cell Type Annotation (UMAP)") + theme_minimal()

ggsave(file.path(plots_dir, "umap_celltype_enhanced.pdf"), p_umap_enhanced, width = 8, height = 6)
print(p_umap_enhanced)

# 15. Prompt for gene and run gene expression analysis (SAVES ALL PLOTS AS PDF)
gene_of_interest <- readline(prompt = "Enter the gene symbol you want to analyze (e.g., FAP): ")
if (gene_of_interest == "") {
  gene_of_interest <- "FAP"
  cat("No gene entered. Defaulting to FAP.\n")
}
gene_analysis <- analyze_gene_expression(seurat_obj, gene = gene_of_interest, output_dir = plots_dir)

# 16. Save the processed Seurat object
saveRDS(seurat_obj, file = file.path(project_dir, "seurat_obj_processed.rds"))
