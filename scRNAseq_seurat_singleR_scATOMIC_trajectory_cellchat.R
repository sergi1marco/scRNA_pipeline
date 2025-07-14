# A Comprehensive and Reproducible R Pipeline for Single-Cell RNA Sequencing Data Analysis

# Set RETICULATE_CONDA_PATH and disable auto-installation
options(reticulate.conda_auto_install = FALSE)
options(reticulate.install_miniforge = FALSE)

# Point reticulate to the Python executable for Rmagic/scATOMIC environment
reticulate::use_python(python = "C:/miniconda3/envs/r_magic_env/python.exe", required = TRUE)

# Global options
options(future.globals.maxSize = 8000 * 1024^2)
options(scipen = 100)

# --- 0. Library Loading ---
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
library(MAST)
library(limma)
library(edgeR)
library(ComplexHeatmap)
library(tidyverse) # Includes ggplot2, dplyr, etc.
library(cowplot)
library(patchwork)
library(ggrepel)
library(viridis)
library(data.table)
library(ggraph)
library(igraph)
library(RColorBrewer)
library(clusterProfiler)
library(org.Hs.eg.db) # For gene ID conversion
library(AnnotationDbi) # For mapIds
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
library(Rmagic)
library(scATOMIC)
library(cutoff.scATOMIC) # scATOMIC uses this
library(grid) # For unit() functions, essential for custom plots
library(ggforce) # For geom_arc_bar in pie charts
library(CellChat) # Make sure this is loaded and updated!
library(scales) # ADDED: For generating color palettes (used in CellChat color assignment)
library(grDevices) # ADDED: For rainbow() function

message("\n--- All required libraries loaded successfully. ---")
message("Verifying Python environment for scATOMIC:")
py_config() # Print Python configuration to confirm correct environment


# --- Custom Helper Function for UMAP Plot Styling with Direct Arrows ---
add_umap_arrows_and_theme_void <- function(base_plot, seurat_obj_for_coords,
                                           arrow_length_prop = 0.05, arrow_offset_prop = 0.02,
                                           umap_dim_names = c("UMAP1", "UMAP2"),
                                           show_legend = TRUE) {
  
  if (!inherits(base_plot, "ggplot")) {
    stop("Input 'base_plot' must be a ggplot object.")
  }
  
  umap_coords <- Embeddings(seurat_obj_for_coords, reduction = "umap")
  min_x <- min(umap_coords[,1], na.rm = TRUE)
  max_x <- max(umap_coords[,1], na.rm = TRUE)
  min_y <- min(umap_coords[,2], na.rm = TRUE)
  max_y <- max(umap_coords[,2], na.rm = TRUE)
  
  range_x <- max_x - min_x
  range_y <- max_y - min_y
  
  effective_arrow_length <- min(range_x, range_y) * arrow_length_prop
  effective_arrow_offset <- min(range_x, range_y) * arrow_offset_prop
  
  start_x_base <- min_x + effective_arrow_offset
  start_y_base <- min_y + effective_arrow_offset
  
  end_x_umap1 <- start_x_base + effective_arrow_length
  end_y_umap1 <- start_y_base
  
  end_x_umap2 <- start_x_base
  end_y_umap2 <- start_y_base + effective_arrow_length
  
  arrow_style <- ggplot2::arrow(length = grid::unit(0.08, "inches"), type = "closed", ends = "last")
  
  custom_plot <- base_plot +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 14),
      legend.position = if (show_legend) "right" else "none",
      legend.box.just = "left",
      legend.direction = "vertical",
      legend.key.size = grid::unit(10, "pt"),
      legend.background = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(t = 5, r = 5, b = 15, l = 15, unit = "pt")
    ) +
    geom_segment(aes(x = start_x_base, y = start_y_base, xend = end_x_umap1, yend = end_y_umap1),
                 arrow = arrow_style, color = "black", linewidth = 0.6) +
    annotate("text", x = end_x_umap1 + effective_arrow_offset * 0.75,
             y = end_y_umap1, label = umap_dim_names[1],
             hjust = 0, vjust = 0.5, size = 3, fontface = "bold") +
    geom_segment(aes(x = start_x_base, y = start_y_base, xend = end_x_umap2, yend = end_y_umap2),
                 arrow = arrow_style, color = "black", linewidth = 0.6) +
    annotate("text", x = end_x_umap2,
             y = end_y_umap2 + effective_arrow_offset * 1.5,
             label = umap_dim_names[2],
             hjust = 0.5, vjust = 0, size = 3, fontface = "bold", angle = 90) +
    coord_cartesian(xlim = c(min_x - range_x * 0.05, max_x + range_x * 0.05),
                    ylim = c(min_y - range_y * 0.05, max_y + range_y * 0.05),
                    expand = FALSE,
                    clip = "off")
  
  return(custom_plot)
}

# New helper function for gene symbol conversion
convert_gene_symbols <- function(seurat_obj) {
  message("Starting gene symbol conversion...")
  current_genes <- rownames(seurat_obj)
  original_gene_count <- length(current_genes)
  
  # Try to remove version numbers from Ensembl-like IDs if present
  # Pattern for ENSG/ENSMBL IDs (e.g., ENSG00000123456.1 -> ENSG00000123456)
  stripped_genes <- gsub("\\.[0-9]+$", "", current_genes)
  
  # Map all stripped genes to SYMBOL using org.Hs.eg.db
  # use_names_as_keys = TRUE is important to map original names
  mapped_symbols <- AnnotationDbi::mapIds(
    org.Hs.eg.db,
    keys = stripped_genes,
    column = "SYMBOL",
    keytype = "ENSEMBL", # Assuming ENSEMBL-like IDs are the primary unmapped type
    multiVals = "first" # Take the first symbol if multiple map
  )
  
  # Handle cases where original genes were already symbols, or couldn't be mapped via ENSEMBL
  # Fill NA values in mapped_symbols with original gene names
  mapped_symbols[is.na(mapped_symbols)] <- current_genes[is.na(mapped_symbols)]
  
  # Create a data frame for mapping
  gene_map_df <- data.frame(
    original_gene = current_genes,
    mapped_symbol = mapped_symbols,
    stringsAsFactors = FALSE
  )
  
  # Filter out genes that could not be mapped to any reasonable symbol (e.g., if they were just random strings)
  # A simple check: if the mapped_symbol is still exactly the same as the original_gene AND it's not a known symbol
  # This is tricky; for now, we'll keep everything that got *some* mapping or was already present.
  # The strict filtering happens when CellChat queries its database.
  
  # Identify unique new symbols
  unique_new_symbols <- unique(gene_map_df$mapped_symbol)
  message(paste0("Number of unique new gene symbols after mapping: ", length(unique_new_symbols)))
  
  # Create a matrix for the new counts/data, summing up duplicates
  new_counts <- matrix(0, nrow = length(unique_new_symbols), ncol = ncol(seurat_obj))
  rownames(new_counts) <- unique_new_symbols
  colnames(new_counts) <- colnames(seurat_obj)
  
  new_data <- matrix(0, nrow = length(unique_new_symbols), ncol = ncol(seurat_obj))
  rownames(new_data) <- unique_new_symbols
  colnames(new_data) <- colnames(seurat_obj)
  
  # Aggregate counts and data by new symbol
  for (symbol in unique_new_symbols) {
    original_genes_for_symbol <- gene_map_df$original_gene[gene_map_df$mapped_symbol == symbol]
    
    # Ensure original_genes_for_symbol are present in the Seurat object's assay
    valid_original_genes <- intersect(original_genes_for_symbol, rownames(GetAssayData(seurat_obj, assay = "RNA", slot = "counts")))
    
    if (length(valid_original_genes) > 0) {
      if (length(valid_original_genes) == 1) {
        new_counts[symbol, ] <- GetAssayData(seurat_obj, assay = "RNA", slot = "counts")[valid_original_genes, ]
        new_data[symbol, ] <- GetAssayData(seurat_obj, assay = "RNA", slot = "data")[valid_original_genes, ]
      } else {
        new_counts[symbol, ] <- colSums(GetAssayData(seurat_obj, assay = "RNA", slot = "counts")[valid_original_genes, ])
        new_data[symbol, ] <- colSums(GetAssayData(seurat_obj, assay = "RNA", slot = "data")[valid_original_genes, ])
      }
    } else {
      # If no valid original genes map to this symbol (shouldn't happen if mapped_symbols are from original_genes)
      # But as a safeguard, set to zero.
      new_counts[symbol, ] <- 0
      new_data[symbol, ] <- 0
    }
  }
  
  # Replace the RNA assay data with the new, symbol-mapped data
  seurat_obj[["RNA"]] <- CreateAssayObject(counts = new_counts)
  seurat_obj[["RNA"]]@data <- new_data # Manually set the data slot
  
  message(paste0("Gene symbol conversion complete. Original genes: ", original_gene_count,
                 ", Genes after conversion: ", nrow(seurat_obj[["RNA"]]), " (unique symbols)."))
  
  return(seurat_obj)
}


# Function for scATOMIC-based malignancy detection and annotation
perform_scatomic_malignancy_detection <- function(seurat_obj) {
  message("Converting Seurat object to matrix for scATOMIC (using RNA assay raw counts)...")
  if (!"RNA" %in% names(seurat_obj@assays)) {
    stop("RNA assay not found. scATOMIC requires raw counts from 'RNA' assay.")
  }
  input_matrix <- GetAssayData(seurat_obj, assay = "RNA", layer = "counts")
  
  message("Prompting for known cancer type for scATOMIC (optional but recommended)...")
  cancer_type_input <- readline(prompt = "Enter known cancer type (e.g., 'Lung_Adenocarcinoma', 'Breast_Ductal_Carcinoma') or leave blank for general: ")
  known_cancer_type_param <- if (nchar(cancer_type_input) > 0) cancer_type_input else NULL
  
  message("Running scATOMIC for malignancy prediction (this can be computationally intensive)...")
  scatomic_predictions <- scATOMIC::run_scATOMIC(
    rna_counts = input_matrix,
    imputation = TRUE,
    mc.cores = 7,
    confidence_cutoff = TRUE
  )
  message("scATOMIC prediction complete. Generating summary matrix...")
  
  final_scatomic_annotations <- scATOMIC::create_summary_matrix(
    prediction_list = scatomic_predictions,
    raw_counts = input_matrix,
    use_CNVs = FALSE,
    modify_results = TRUE,
    known_cancer_type = known_cancer_type_param,
    normal_tissue = FALSE,
    low_res_mode = FALSE
  )
  message("scATOMIC summary matrix created.")
  
  if (!"scATOMIC_pred" %in% colnames(final_scatomic_annotations)) {
    stop("scATOMIC::create_summary_matrix output 'scATOMIC_pred' column not found.")
  }
  
  common_cells_scatomic <- intersect(rownames(final_scatomic_annotations), colnames(seurat_obj))
  
  if (length(common_cells_scatomic) == 0) {
    stop("No common cells found between scATOMIC results and Seurat object. Cannot merge.")
  }
  
  seurat_obj_subsetted <- subset(seurat_obj, cells = common_cells_scatomic)
  
  predictions_for_seurat_raw <- final_scatomic_annotations[colnames(seurat_obj_subsetted), "scATOMIC_pred"]
  
  problematic_string <- "source(\"~/.active−rstudio−document\", echo = TRUE)"
  predictions_for_seurat <- ifelse(predictions_for_seurat_raw == problematic_string, "Unknown_scATOMIC_Cell", predictions_for_seurat_raw)
  
  seurat_obj_subsetted <- AddMetaData(seurat_obj_subsetted, metadata = predictions_for_seurat, col.name = "scATOMIC_prediction")
  
  seurat_obj_subsetted$scATOMIC_malignancy <- seurat_obj_subsetted$scATOMIC_prediction
  if(problematic_string %in% levels(seurat_obj_subsetted$scATOMIC_malignancy)) {
    levels(seurat_obj_subsetted$scATOMIC_malignancy) <- levels(seurat_obj_subsetted$scATOMIC_malignancy)[levels(seurat_obj_subsetted$scATOMIC_malignancy) != problematic_string]
    message("Removed problematic 'source(...)' level from scATOMIC_malignancy.")
  }
  seurat_obj_subsetted$scATOMIC_malignancy <- factor(seurat_obj_subsetted$scATOMIC_malignancy,
                                                     levels = sort(unique(seurat_obj_subsetted$scATOMIC_malignancy)))
  
  
  message("scATOMIC cell type predictions added to Seurat object metadata.")
  
  output_filename_scatomic_umap <- "UMAP_scATOMIC_Predicted_CellTypes.pdf"
  pdf(output_filename_scatomic_umap, width = 10, height = 8)
  plot_scatomic_umap <- DimPlot(seurat_obj_subsetted,
                                reduction = "umap",
                                group.by = "scATOMIC_prediction",
                                label = TRUE,
                                repel = TRUE,
                                pt.size = 0.5) +
    ggplot2::ggtitle("UMAP by scATOMIC Predicted Cell Types") +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 3), ncol = 1))
  print(add_umap_arrows_and_theme_void(plot_scatomic_umap, seurat_obj_subsetted, show_legend = TRUE))
  dev.off()
  message(paste0("UMAP plot by scATOMIC predicted cell types saved to '", output_filename_scatomic_umap, "'."))
  
  message("scATOMIC detection step completed.")
  
  return(seurat_obj_subsetted)
}

# --- SCP-like Plotting Functions (from SCP_plots.R) ---

`%||%` <- function(a, b) {
  if (!is.null(a)) a else b
}

col2hex <- function(color_vector) {
  sapply(color_vector, function(color) {
    if (is.character(color) && length(color) == 1 && nchar(color) %in% c(7, 9) && substr(color, 1, 1) == "#") {
      return(color)
    }
    tryCatch({
      rgb_val <- grDevices::col2rgb(color)
      grDevices::rgb(rgb_val[1], rgb_val[2], rgb_val[3], maxColorValue = 255)
    }, error = function(e) {
      warning(paste("Invalid color specified:", color, "- returning black."), call. = FALSE)
      return("#000000")
    })
  }, USE.NAMES = FALSE)
}

DefaultReduction <- function(object, assay = NULL, min_dim = 2, pattern = NULL) {
  if (!is.null(pattern)) {
    if (pattern %in% Seurat::Reductions(object)) {
      return(pattern)
    } else {
      available_reductions <- Seurat::Reductions(object)
      matched_reduction <- grep(pattern, available_reductions, ignore.case = TRUE, value = TRUE)
      if (length(matched_reduction) > 0) {
        return(matched_reduction[1])
      }
    }
  }
  
  reductions <- Seurat::Reductions(object)
  preferred_reductions <- c("umap", "tsne", "pca")
  
  for (red in preferred_reductions) {
    if (red %in% reductions) {
      if (ncol(Seurat::Embeddings(object, reduction = red)) >= min_dim) {
        return(red)
      }
    }
  }
  
  for (red in reductions) {
    if (ncol(Seurat::Embeddings(object, reduction = red)) >= min_dim) {
      return(red)
    }
  }
  
  stop("No suitable dimensionality reduction found with at least ", min_dim, " dimensions.")
}

palette_list_internal <- list(
  "Paired" = c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928"),
  "Spectral" = c("#9E0142", "#D53E4F", "#F46D43", "#FDAE61", "#FEE08B", "#FFFFBF", "#E6F598", "#ABDDA4", "#66C2A5", "#3288BD", "#5E4FA2"),
  "jet" = c("#00007F", "#0000FF", "#007FFF", "#00FFFF", "#7FFF7F", "#FFFF00", "#FF7F00", "#FF0000", "#7F0000"),
  "Dark2" = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666"),
  "Set1" = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999"),
  "Greys" = c("#FFFFFF", "#F0F0F0", "#D9D9D9", "#BDBDBD", "#969696", "#737373", "#525252", "#252525", "#000000"),
  "RdYlBu" = c("#A50026", "#D73027", "#F46D43", "#FDAE61", "#FEE090", "#FFFFBF", "#E0F3F8", "#ABDDA4", "#66C2A5", "#3288BD", "#5E4FA2")
)
attr(palette_list_internal$Spectral, "type") <- "continuous"
attr(palette_list_internal$jet, "type") <- "continuous"
attr(palette_list_internal$Paired, "type") <- "discrete"
attr(palette_list_internal$Dark2, "type") <- "discrete"
attr(palette_list_internal$Set1, "type") <- "discrete"
attr(palette_list_internal$Greys, "type") <- "continuous"
attr(palette_list_internal$RdYlBu, "type") <- "continuous"

palette_scp <- function(x, n = 100, palette = "Paired", palcolor = NULL, type = "auto",
                        matched = FALSE, reverse = FALSE, NA_keep = FALSE, NA_color = "grey80") {
  palette_list <- palette_list_internal
  if (missing(x)) {
    x <- 1:n
    type <- "continuous"
  }
  if (!palette %in% names(palette_list)) {
    stop("The palette name (", palette, ") is invalid! You can check the available palette names with 'show_palettes()'. Or pass palette colors via the 'palcolor' parameter.")
  }
  if (is.list(palcolor)) {
    palcolor <- unlist(palcolor)
  }
  if (all(palcolor == "")) {
    palcolor <- palette_list[[palette]]
  }
  if (is.null(palcolor) || length(palcolor) == 0) {
    palcolor <- palette_list[[palette]]
  }
  if (!is.null(names(palcolor))) {
    if (all(x %in% names(palcolor))) {
      palcolor <- palcolor[intersect(names(palcolor), x)]
    }
  }
  pal_n <- length(palcolor)
  
  if (!type %in% c("auto", "discrete", "continuous")) {
    stop("'type' must be one of 'auto','discrete' and 'continuous'.")
  }
  if (type == "auto") {
    if (is.numeric(x)) {
      type <- "continuous"
    } else {
      type <- "discrete"
    }
  }
  
  if (type == "discrete") {
    if (!is.factor(x)) {
      x <- factor(x, levels = unique(x))
    }
    n_x <- nlevels(x)
    if (isTRUE(attr(palcolor, "type") == "continuous")) {
      color <- grDevices::colorRampPalette(palcolor)(n_x)
    } else {
      color <- ifelse(rep(n_x, n_x) <= pal_n,
                      palcolor[1:n_x],
                      grDevices::colorRampPalette(palcolor)(n_x)
      )
    }
    names(color) <- levels(x)
    if (any(is.na(x))) {
      color <- c(color, stats::setNames(NA_color, "NA"))
    }
    if (isTRUE(matched)) {
      color <- color[x]
      color[is.na(color)] <- NA_color
    }
  } else if (type == "continuous") {
    if (!is.numeric(x)) {
      stop("'x' must be of numeric type when using continuous color palettes.")
    }
    
    color_gradient_steps <- grDevices::colorRampPalette(palcolor)(n)
    
    if (isTRUE(reverse)) {
      color_gradient_steps <- rev(color_gradient_steps)
    }
    
    return(color_gradient_steps)
  }
  
  if (isTRUE(reverse)) {
    color <- rev(color)
  }
  if (!isTRUE(NA_keep)) {
    color <- color[names(color) != "NA"]
  }
  return(color)
}


CellDimPlot <- function(srt, group.by, reduction = NULL, dims = c(1, 2), split.by = NULL, cells = NULL,
                        show_na = FALSE, show_stat = FALSE,
                        pt.size = NULL, pt.alpha = 1, palette = "Paired", palcolor = NULL, bg_color = "grey80",
                        label = FALSE, label.size = 4, label.fg = "white", label.bg = "black", label.bg.r = 0.1,
                        label_repel = FALSE, label_repulsion = 20, label_point_size = 1, label_point_color = "black", label_segment_color = "black",
                        stat.by = NULL, stat_plot_type = c("pie"), stat_plot_size = 0.05, stat_plot_label = FALSE, stat_plot_label_size = 3,
                        aspect.ratio = 1, title = NULL, subtitle = NULL, xlab = NULL, ylab = NULL,
                        legend.position = "right", legend.direction = "vertical",
                        combine = TRUE, nrow = NULL, ncol = NULL, byrow = TRUE, force = FALSE, seed = 11) {
  set.seed(seed)
  stat_plot_type <- match.arg(stat_plot_type)
  
  if (is.null(split.by)) {
    split.by <- "All.groups"
    srt@meta.data[[split.by]] <- factor("")
  }
  
  cols_to_check <- unique(c(group.by, split.by))
  if (!is.null(stat.by)) {
    cols_to_check <- unique(c(cols_to_check, stat.by))
  }
  
  for (col_name in cols_to_check) {
    if (!col_name %in% colnames(srt@meta.data)) {
      stop(paste0("Metadata column '", col_name, "' is not in the meta.data of srt object."))
    }
    if (!is.factor(srt@meta.data[[col_name]])) {
      srt@meta.data[[col_name]] <- factor(srt@meta.data[[col_name]], levels = unique(srt@meta.data[[col_name]]))
    }
    if (isTRUE(show_na) && any(is.na(srt@meta.data[[col_name]]))) {
      raw_levels <- unique(c(levels(srt@meta.data[[col_name]]), "NA"))
      srt@meta.data[[col_name]] <- as.character(srt@meta.data[[col_name]])
      srt@meta.data[[col_name]][is.na(srt@meta.data[[col_name]])] <- "NA"
      srt@meta.data[[col_name]] <- factor(srt@meta.data[[col_name]], levels = raw_levels)
    }
  }
  
  if (is.null(reduction)) {
    reduction <- DefaultReduction(srt)
  } else {
    reduction <- DefaultReduction(srt, pattern = reduction)
  }
  if (!reduction %in% names(srt@reductions)) {
    stop(paste0(reduction, " is not in the srt reduction names."))
  }
  
  reduction_key <- Seurat::Key(srt@reductions[[reduction]])
  dat_dim <- Seurat::Embeddings(srt, reduction = reduction)
  colnames(dat_dim) <- paste0(reduction_key, dims)
  
  dat_meta_orig <- srt@meta.data[, unique(c(group.by, split.by, stat.by)), drop = FALSE]
  
  common_cells <- intersect(rownames(dat_dim), rownames(dat_meta_orig))
  if (length(common_cells) == 0) {
    stop("No common cells found across reduction data and metadata. Check Seurat object integrity.")
  }
  
  dat_dim_ordered <- dat_dim[common_cells, , drop = FALSE]
  dat_meta_ordered <- dat_meta_orig[common_cells, , drop = FALSE]
  dat_use <- cbind(dat_dim_ordered, dat_meta_ordered)
  
  if (!is.null(cells)) {
    dat_use <- dat_use[intersect(rownames(dat_use), cells), , drop = FALSE]
  }
  if (is.null(pt.size)) {
    pt.size <- min(3000 / nrow(dat_use), 0.5)
  }
  
  plist <- list()
  xlab <- xlab %||% paste0(reduction_key, dims[1])
  ylab <- ylab %||% paste0(reduction_key, dims[2])
  
  comb <- expand.grid(split = levels(dat_use[[split.by]]), group = group.by, stringsAsFactors = FALSE)
  rownames(comb) <- paste0(comb[["split"]], ":", comb[["group"]])
  
  stat_by_levels <- if (!is.null(stat.by)) levels(dat_use[[stat.by]]) else NULL
  stat_by_colors <- if (!is.null(stat.by)) palette_scp(stat_by_levels, palette = "Set1", type = "discrete") else NULL
  
  plist <- lapply(setNames(rownames(comb), rownames(comb)), function(i) {
    legend_list <- list()
    
    g <- comb[i, "group"]
    s <- comb[i, "split"]
    
    dat_filtered <- dat_use
    if (s != "") {
      cells_mask <- dat_filtered[[split.by]] != s
      dat_filtered[[g]][cells_mask] <- NA
      if (!is.null(stat.by)) dat_filtered[[stat.by]][cells_mask] <- NA
    }
    
    dat_filtered <- dat_filtered %>% filter(!is.na(.data[[g]]))
    
    labels_tb <- table(dat_filtered[[g]])
    labels_tb <- labels_tb[labels_tb != 0]
    colors <- palette_scp(levels(dat_use[[g]]), palette = palette, palcolor = palcolor, NA_keep = TRUE)
    
    label_use <- names(labels_tb)
    
    p <- ggplot2::ggplot(dat_filtered, ggplot2::aes(x = .data[[paste0(reduction_key, dims[1])]], y = .data[[paste0(reduction_key, dims[2])]], color = .data[[g]])) +
      ggplot2::geom_point(size = pt.size, alpha = pt.alpha) +
      ggplot2::scale_color_manual(
        name = paste0(g, ":"),
        values = colors[names(labels_tb)],
        labels = label_use,
        na.value = bg_color,
        guide = ggplot2::guide_legend(
          title.hjust = 0,
          order = 1,
          override.aes = list(size = 4, alpha = 1)
        )
      ) +
      ggplot2::theme_void() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 14),
        legend.position = legend.position,
        legend.box.just = "left",
        legend.direction = "vertical",
        legend.key.size = grid::unit(10, "pt"),
        legend.background = ggplot2::element_blank(),
        plot.margin = ggplot2::margin(t = 5, r = 5, b = 15, l = 15, unit = "pt")
      ) +
      ggplot2::labs(x = xlab, y = ylab, title = title, subtitle = subtitle)
    
    label_df <- dat_filtered %>%
      dplyr::group_by(.data[[g]]) %>%
      dplyr::summarise(
        x = stats::median(.data[[paste0(reduction_key, dims[1])]], na.rm = TRUE),
        y = stats::median(.data[[paste0(reduction_key, dims[2])]], na.rm = TRUE),
        .groups = 'drop'
      ) %>%
      as.data.frame()
    colnames(label_df)[colnames(label_df) == g] <- "label"
    label_df[["label_text"]] <- as.character(label_df[["label"]])
    
    if (isTRUE(label)) {
      if (isTRUE(show_stat)) {
        label_df[["count_info"]] <- vapply(as.character(label_df[["label"]]), function(grp) {
          count_val <- labels_tb[grp]
          if (is.na(count_val)) { return("") } else { return(paste0(" (n=", count_val, ")")) }
        }, FUN.VALUE = character(1), USE.NAMES = FALSE)
        label_df[["label_text"]] <- paste0(label_df[["label_text"]], label_df[["count_info"]])
      }
      
      if (isTRUE(label_repel)) {
        p <- p + ggplot2::annotate(
          geom = "point", x = label_df[["x"]], y = label_df[["y"]],
          color = label_point_color, size = label_point_size
        ) + ggrepel::geom_text_repel(
          data = label_df, ggplot2::aes(x = .data[["x"]], y = .data[["y"]], label = .data[["label_text"]]),
          fontface = "bold", min.segment.length = 0, segment.color = label_segment_color,
          point.size = label_point_size, max.overlaps = 100, force = label_repulsion,
          color = label.fg, bg.color = label.bg, bg.r = label.bg.r, size = label.size, inherit.aes = FALSE
        )
      } else {
        p <- p + ggrepel::geom_text_repel(
          data = label_df, ggplot2::aes(x = .data[["x"]], y = .data[["y"]], label = .data[["label_text"]]),
          fontface = "bold",
          point.size = NA, max.overlaps = 100, force = 0,
          color = label.fg, bg.color = label.bg, bg.r = label.bg.r, size = label.size, inherit.aes = FALSE
        )
      }
    }
    
    if (!is.null(stat.by) && isTRUE(label)) {
      pie_data_per_cluster <- dat_filtered %>%
        dplyr::group_by(.data[[g]], .data[[stat.by]]) %>%
        dplyr::summarise(count = n(), .groups = 'drop') %>%
        dplyr::group_by(.data[[g]]) %>%
        dplyr::mutate(percentage = count / sum(count)) %>%
        dplyr::ungroup()
      
      umap_range_x <- max(dat_dim[, 1], na.rm = TRUE) - min(dat_dim[, 1], na.rm = TRUE)
      umap_range_y <- max(dat_dim[, 2], na.rm = TRUE) - min(dat_dim[, 2], na.rm = TRUE)
      
      size_factor_x <- umap_range_x * stat_plot_size
      size_factor_y <- umap_range_y * stat_plot_size
      
      for (cluster_label in unique(pie_data_per_cluster[[g]])) {
        current_cluster_data <- pie_data_per_cluster %>%
          dplyr::filter(.data[[g]] == cluster_label) %>%
          dplyr::arrange(.data[[stat.by]])
        
        centroid <- label_df %>%
          dplyr::filter(.data[["label"]] == cluster_label) %>%
          dplyr::select(x, y)
        
        if (nrow(centroid) == 0) next
        
        pie_chart <- ggplot(current_cluster_data, aes(x = 1, y = percentage, fill = .data[[stat.by]])) +
          geom_bar(width = 1, stat = "identity", color = "white", linewidth = 0.5) +
          coord_polar("y", start = 0) +
          scale_fill_manual(values = stat_by_colors, drop = FALSE) + # Use all defined levels
          theme_void() +
          theme(
            legend.position = "none",
            plot.margin = ggplot2::margin(0,0,0,0)
          )
        
        pie_grob <- ggplot2::ggplotGrob(pie_chart)
        
        p <- p + annotation_custom(
          grob = pie_grob,
          xmin = centroid$x - size_factor_x / 2,
          xmax = centroid$x + size_factor_x / 2,
          ymin = centroid$y - size_factor_y / 2,
          ymax = centroid$y + size_factor_y / 2
        )
      }
      
      stat_legend_plot <- ggplot(data.frame(var = stat_by_levels, dummy = 1), aes(x = dummy, y = var, fill = var)) +
        geom_tile() +
        scale_fill_manual(values = stat_by_colors, drop = FALSE) +
        guides(fill = guide_legend(title = stat.by, override.aes = list(size = 3))) +
        theme_void() +
        theme(legend.position = "right")
      
      stat_legend_grob <- cowplot::get_legend(stat_legend_plot)
      legend_list <- c(legend_list, list(stat_legend_grob))
    }
    
    p_final <- p
    
    if (length(legend_list) > 0) {
      legend_combined <- cowplot::plot_grid(plotlist = legend_list, ncol = 1, align = "v", axis = "l")
      p_final <- patchwork::wrap_plots(p + ggplot2::theme(legend.position = "none"), legend_combined, ncol = 2, widths = c(3,1))
    }
    
    return(p_final)
  })
  
  if (isTRUE(combine)) {
    if (length(plist) > 1) {
      plot <- patchwork::wrap_plots(plotlist = plist, nrow = nrow, ncol = ncol, byrow = byrow)
    } else {
      plot <- plist[[1]]
    }
    return(plot)
  } else {
    return(plist)
  }
}

FeatureDimPlot <- function(srt, variable, reduction = NULL, dims = c(1, 2), split.by = NULL, cells = NULL,
                           show_na = FALSE, show_stat = FALSE,
                           pt.size = NULL, pt.alpha = 1, palette = "Greys", palcolor = NULL, bg_color = "grey80",
                           label = FALSE, label.size = 4, label.fg = "white", label.bg = "black", label.bg.r = 0.1,
                           label_insitu = FALSE, label_repel = FALSE, label_repulsion = 20,
                           label_point_size = 1, label_point_color = "black", label_segment_color = "black",
                           cells.highlight = NULL, cols.highlight = "black", sizes.highlight = 1, alpha.highlight = 1, stroke.highlight = 0.5,
                           aspect.ratio = 1, title = NULL, subtitle = NULL, xlab = NULL, ylab = NULL,
                           legend.position = "right", legend.direction = "vertical",
                           theme_use = "theme_scp", theme_args = list(),
                           combine = TRUE, nrow = NULL, ncol = NULL, byrow = TRUE, force = FALSE, seed = 11) {
  set.seed(seed)
  
  if (is.null(split.by)) {
    split.by <- "All.groups"
    srt@meta.data[[split.by]] <- factor("")
  }
  for (i in unique(c(split.by))) {
    if (!i %in% colnames(srt@meta.data)) {
      stop(paste0(i, " is not in the meta.data of srt object."))
    }
    if (!is.factor(srt@meta.data[[i]])) {
      srt@meta.data[[i]] <- factor(srt@meta.data[[i]], levels = unique(srt@meta.data[[i]]))
    }
    if (isTRUE(show_na) && any(is.na(srt@meta.data[[i]]))) {
      raw_levels <- unique(c(levels(srt@meta.data[[i]]), "NA"))
      srt@meta.data[[i]] <- as.character(srt@meta.data[[i]])
      srt@meta.data[[i]][is.na(srt@meta.data[[i]])] <- "NA"
      srt@meta.data[[i]] <- factor(srt@meta.data[[i]], levels = raw_levels)
    }
  }
  
  if (!variable %in% colnames(srt@meta.data)) {
    stop(paste0("Variable '", variable, "' not found in srt@meta.data."))
  }
  if (!is.numeric(srt@meta.data[[variable]])) {
    stop(paste0("Variable '", variable, "' must be numeric for FeatureDimPlot (continuous scale)."))
  }
  
  if (is.null(reduction)) {
    reduction <- DefaultReduction(srt)
  } else {
    reduction <- DefaultReduction(srt, pattern = reduction)
  }
  if (!reduction %in% names(srt@reductions)) {
    stop(paste0(reduction, " is not in the srt reduction names."))
  }
  
  dat_meta_orig <- srt@meta.data[, unique(c(split.by, variable)), drop = FALSE]
  
  reduction_key <- Seurat::Key(srt@reductions[[reduction]])
  dat_dim <- Seurat::Embeddings(srt, reduction = reduction)
  colnames(dat_dim) <- paste0(reduction_key, dims)
  
  common_cells <- intersect(rownames(dat_dim), rownames(dat_meta_orig))
  
  if (length(common_cells) == 0) {
    stop("No common cells found across reduction data and metadata. Check Seurat object integrity.")
  }
  
  dat_dim_ordered <- dat_dim[common_cells, , drop = FALSE]
  dat_meta_ordered <- dat_meta_orig[common_cells, , drop = FALSE]
  
  dat_use <- cbind(dat_dim_ordered, dat_meta_ordered)
  
  if (!is.null(cells)) {
    dat_use <- dat_use[intersect(rownames(dat_use), cells), , drop = FALSE]
  }
  if (is.null(pt.size)) {
    pt.size <- min(3000 / nrow(dat_use), 0.5)
  }
  
  plist <- list()
  xlab <- xlab %||% paste0(reduction_key, dims[1])
  ylab <- ylab %||% paste0(reduction_key, dims[2])
  
  comb <- expand.grid(split = levels(dat_use[[split.by]]), variable = variable, stringsAsFactors = FALSE)
  rownames(comb) <- paste0(comb[["split"]], ":", comb[["variable"]])
  plist <- lapply(setNames(rownames(comb), rownames(comb)), function(i) {
    legend_list <- list()
    
    f <- comb[i, "variable"]
    s <- comb[i, "split"]
    
    dat <- dat_use
    cells_mask <- dat[[split.by]] != s
    dat[[f]][cells_mask] <- NA
    
    feature_data_for_plotting <- as.numeric(dat[[f]])
    colors <- palette_scp(feature_data_for_plotting, palette = palette, palcolor = palcolor, type = "continuous", NA_keep = TRUE)
    
    p <- ggplot2::ggplot(dat, ggplot2::aes(x = .data[[paste0(reduction_key, dims[1])]], y = .data[[paste0(reduction_key, dims[2])]], color = as.numeric(.data[[f]]))) +
      ggplot2::geom_point(size = pt.size, alpha = pt.alpha)
    
    p <- p + ggplot2::scale_color_gradientn(
      colors = colors,
      na.value = bg_color,
      guide = ggplot2::guide_colorbar(
        title = f,
        order = 1
      )
    ) + ggplot2::labs(color = f)
    
    p_base <- p
    
    if (length(legend_list) > 0) {
      legend_list <- legend_list[!sapply(legend_list, is.null)]
      legend_base <- cowplot::get_legend(p_base)
      if (legend.direction == "vertical") {
        legend <- do.call(cowplot::plot_grid, c(list(plotlist = c(list(legend_base), legend_list), ncol = 1), list(align = "v", axis = "l", rel_heights = c(1, rep(0.5, length(legend_list))))))
      } else {
        legend <- do.call(cowplot::plot_grid, c(list(plotlist = c(list(legend_base), legend_list), nrow = 1), list(align = "h", axis = "t", rel_widths = c(1, rep(0.5, length(legend_list))))))
      }
      p <- patchwork::wrap_plots(p + ggplot2::theme(legend.position = "none"), legend, ncol = 2, widths = c(3,1))
    }
    return(p)
  })
  
  if (isTRUE(combine)) {
    if (length(plist) > 1) {
      plot <- patchwork::wrap_plots(plotlist = plist, nrow = nrow, ncol = ncol, byrow = byrow)
    } else {
      plot <- plist[[1]]
    }
    return(plot)
  } else {
    return(plist)
  }
}

add_umap_arrows_for_scp <- function(plot, srt_object, reduction_name, show_legend = TRUE) {
  reduction_key <- Seurat::Key(srt_object@reductions[[reduction_name]])
  embeddings <- Seurat::Embeddings(srt_object, reduction = reduction_name)
  
  x_min <- min(embeddings[, 1], na.rm = TRUE)
  x_max <- max(embeddings[, 1], na.rm = TRUE)
  y_min <- min(embeddings[, 2], na.rm = TRUE)
  y_max <- max(embeddings[, 2], na.rm = TRUE)
  
  range_x <- x_max - x_min
  range_y <- y_max - y_min
  
  arrow_length_prop <- 0.05
  arrow_offset_prop <- 0.02
  
  effective_arrow_length <- min(range_x, range_y) * arrow_length_prop
  effective_arrow_offset <- min(range_x, range_y) * arrow_offset_prop
  
  start_x_base <- x_min + effective_arrow_offset
  start_y_base <- y_min + effective_arrow_offset
  
  end_x_umap1 <- start_x_base + effective_arrow_length
  end_y_umap1 <- start_y_base
  
  end_x_umap2 <- start_x_base
  end_y_umap2 <- start_y_base + effective_arrow_length
  
  arrow_style <- ggplot2::arrow(length = grid::unit(0.08, "inches"), type = "closed", ends = "last")
  
  plot_with_arrows <- plot +
    ggplot2::geom_segment(
      aes(x = start_x_base, y = start_y_base, xend = end_x_umap1, yend = end_y_umap1),
      arrow = arrow_style, color = "black", linewidth = 0.6, inherit.aes = FALSE
    ) +
    ggplot2::annotate("text", x = end_x_umap1 + effective_arrow_offset * 0.75,
                      y = end_y_umap1, label = paste0(reduction_key, "1"),
                      hjust = 0, vjust = 0.5, size = 3, fontface = "bold", color = "black") +
    ggplot2::geom_segment(
      aes(x = start_x_base, y = start_y_base, xend = end_x_umap2, yend = end_y_umap2),
      arrow = arrow_style, color = "black", linewidth = 0.6, inherit.aes = FALSE
    ) +
    ggplot2::annotate("text", x = end_x_umap2,
                      y = end_y_umap2 + effective_arrow_offset * 1.5,
                      label = paste0(reduction_key, "2"),
                      hjust = 0.5, vjust = 0, size = 3, fontface = "bold", color = "black", angle = 90) +
    ggplot2::coord_cartesian(xlim = c(x_min - range_x * 0.05, x_max + range_x * 0.05),
                             ylim = c(y_min - range_y * 0.05, y_max + range_y * 0.05),
                             expand = FALSE,
                             clip = "off") +
    ggplot2::theme(legend.position = ifelse(show_legend, "right", "none"))
  
  return(plot_with_arrows)
}


# --- 1. Interactive Directory Selection & Pipeline Start Decision ---

# Function to select a directory
select_data_directory <- function() {
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    message("Using RStudio API for directory selection...")
    data_path <- rstudioapi::selectDirectory(
      caption = "Select the working directory",
      path = getwd()
    )
  } else {
    message("Using base R for directory selection...")
    data_path <- choose.dir(caption = "Select the working directory")
  }
  
  if (is.null(data_path) || data_path == "") {
    stop("Directory selection cancelled. Select a valid directory.")
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
      message("3: Gene Symbol Conversion") # NEW STEP
      message("4: Dimensionality Reduction (PCA and UMAP)")
      message("5: Visualization")
      message("6: Find Marker Genes for Clusters")
      message("7: Cell Type Annotation & CAF Identification (SingleR)")
      message("8: Differential Expression Analysis")
      message("9: Malignancy Detection (scATOMIC Classifier)")
      message("10: Cell Cycle Scoring")
      message("11: Generate SCP-like Cell Cycle Plots")
      message("12: Interactive Gene Expression Analysis & Visualization")
      message("13: Trajectory inference and Pseudotime Analysis (Monocle3)")
      message("14: Cell-Cell Communication Analysis (CellChat)") # NEW STEP
      message("15: Final Data Saving") # NEW STEP
      message("-----------------------")
      
      user_choice_step <- readline(prompt = "From which step would you like to resume (1-15)? ") # Corrected prompt max step
      start_step <- as.integer(user_choice_step)
      
      if (is.na(start_step) || start_step < 1 || start_step > 15) { # Corrected max step
        warning("Invalid step. Enter a number between 1 and 15. Reloading object and re-asking.")
        seurat_object <- NULL
        next
      } else {
        step_names <- c(
          "QC", # 1
          "Normalization", # 2
          "Gene_Symbol_Conversion", # 3 (NEW)
          "DimRed", # 4 (was 3)
          "Visualization", # 5 (was 4)
          "Markers", # 6 (was 5)
          "Annotation_CAF", # 7 (was 6)
          "DE", # 8 (was 7)
          "scATOMIC_Detection", # 9 (was 8)
          "Cell_Cycle_Scoring", # 10 (was 9)
          "Generate_SCP_Plots", # 11 (was 10)
          "Interactive_Gene_Analysis", # 12 (was 11)
          "Trajectory_Inference", # 13 (was 12)
          "CellChat_Analysis", # 14 (was 13, NEW)
          "Final_Save" # 15 (was 14, NEW)
        )
        message(paste("Starting pipeline from step:", step_names[start_step]))
        choice_made <- TRUE
      }
    } else {
      message("Invalid file choice. Enter a valid number or '0'." )
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
        warning("Multiple .h5 files found. Using the first one detected: ", h5_files[1])
      }
      selected_path <- h5_files[1]
      data_type <- "h5"
    } else if (length(mtx_dirs) > 0) {
      if (length(mtx_dirs) > 1) {
        warning("Multiple MTX data directories found. Using the first one detected: ", mtx_dirs[1])
      }
      selected_path <- mtx_dirs[1]
      data_type <- "mtx"
    } else {
      stop("No 10x Genomics H5 file or MTX data directory found. Ensure your data is correctly placed.")
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
          stop("H5 file contains multiple data types, but 'Gene Expression' was not found.")
        }
      }
    } else if (data_type == "mtx") {
      message(paste("Importing data from MTX folder:", selected_path))
      seurat_obj_raw <- Read10X(data.dir = selected_path)
    }
    
    seurat_object <- CreateSeuratObject(counts = seurat_obj_raw, project = "SingleCellAnalysis")
    message("Seurat object created.")
    
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
    
    seurat_obj <- subset(seurat_obj,
                         subset = nFeature_RNA > min_features &
                           nFeature_RNA < max_umis &
                           percent.mt < max_percent_mt)
    seurat_obj <- SetAssayData(seurat_obj, layer = "counts", new.data = GetAssayData(seurat_obj, assay = "RNA", layer = "counts")[rowSums(GetAssayData(seurat_obj, assay = "RNA", layer = "counts") > 0) >= min_cells, ])
    message(paste0("QC and filtering complete. Remaining cells: ", ncol(seurat_obj)))
    return(seurat_obj)
  }
  current_seurat_object <- perform_qc_and_filter(current_seurat_object)
}

# --- 4. Normalization and Scaling ---
if (start_step <= 2) {
  perform_normalization_and_scaling <- function(seurat_obj) {
    message("Normalizing and Scaling data...")
    seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
    seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
    seurat_obj <- ScaleData(seurat_obj, features = rownames(seurat_obj))
    message("Normalization and scaling complete.")
    return(seurat_obj)
  }
  current_seurat_object <- perform_normalization_and_scaling(current_seurat_object)
}

# --- 5. Gene Symbol Conversion (NEW STEP) ---
if (start_step <= 3) {
  message("\n--- Step 5: Performing Gene Symbol Conversion ---")
  # Ensure org.Hs.eg.db is installed
  if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager")
    }
    BiocManager::install("org.Hs.eg.db", update = FALSE, ask = FALSE)
  }
  library(org.Hs.eg.db)
  
  current_seurat_object <- convert_gene_symbols(current_seurat_object)
  message("Gene symbol conversion step completed.")
}


# --- 6. Dimensionality Reduction (PCA and UMAP) ---
if (start_step <= 4) { # Step changed from 3 to 4
  perform_dimensionality_reduction <- function(seurat_obj, dims = 1:30) {
    message("Performing Dimensionality Reduction (PCA and UMAP)...")
    seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))
    print(VizDimLoadings(seurat_obj, dims = 1:2, reduction = "pca"))
    pdf("PCA_ElbowPlot.pdf", width = 8, height = 6)
    print(ElbowPlot(seurat_obj, ndims = 30))
    dev.off()
    message("PCA complete. See 'PCA_ElbowPlot.pdf' to select appropriate number of dimensions.")
    
    seurat_obj <- RunUMAP(seurat_obj, dims = dims)
    seurat_obj <- FindNeighbors(seurat_obj, dims = dims)
    seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
    message("UMAP and Clustering complete.")
    return(seurat_obj)
  }
  current_seurat_object <- perform_dimensionality_reduction(current_seurat_object)
}

# --- 7. Visualization ---
if (start_step <= 5) { # Step changed from 4 to 5
  perform_visualization <- function(seurat_obj) {
    message("Generating visualizations...")
    
    output_filename_cluster_umap <- "UMAP_Clusters.pdf"
    pdf(output_filename_cluster_umap, width = 10, height = 8)
    p_cluster <- DimPlot(seurat_obj, reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE) +
      ggtitle("UMAP by Clusters") +
      guides(color = guide_legend(override.aes = list(size = 3)))
    print(add_umap_arrows_and_theme_void(p_cluster, seurat_obj, show_legend = FALSE))
    dev.off()
    message(paste0("UMAP plot by clusters saved to '", output_filename_cluster_umap, "'."))
    
    if ("orig.ident" %in% colnames(seurat_obj@meta.data)) {
      output_filename_sample_umap <- "UMAP_Samples.pdf"
      pdf(output_filename_sample_umap, width = 10, height = 8)
      p_sample <- DimPlot(seurat_obj, reduction = "umap", group.by = "orig.ident", label = FALSE, repel = TRUE) +
        ggtitle("UMAP by Sample (orig.ident)") +
        guides(color = guide_legend(override.aes = list(size = 3)))
      print(add_umap_arrows_and_theme_void(p_sample, seurat_obj, show_legend = TRUE))
      dev.off()
      message(paste0("UMAP plot by samples saved to '", output_filename_sample_umap, "'."))
    } else {
      message("Warning: 'orig.ident' not found in metadata. Skipping UMAP by samples.")
    }
    
    message("Visualizations complete.")
  }
  perform_visualization(current_seurat_object)
}

# --- 8. Find Marker Genes for Clusters ---
if (start_step <= 6) { # Step changed from 5 to 6
  find_cluster_markers <- function(seurat_obj) {
    message("Finding marker genes for each cluster (this may take some time)...")
    seurat_obj.markers <- FindAllMarkers(seurat_obj,
                                         only.pos = TRUE,
                                         min.pct = 0.25,
                                         logfc.threshold = 0.25,
                                         test.use = "wilcox")
    
    top_n_markers <- seurat_obj.markers %>%
      group_by(cluster) %>%
      top_n(n = 10, wt = avg_log2FC)
    
    write.csv(seurat_obj.markers, "All_Cluster_Markers.csv", row.names = FALSE)
    write.csv(top_n_markers, "Top10_Cluster_Markers.csv", row.names = FALSE)
    
    message("Cluster markers identified and saved to 'All_Cluster_Markers.csv' and 'Top10_Cluster_Markers.csv'.")
    
    return(seurat_obj.markers)
  }
  cluster_markers <- find_cluster_markers(current_seurat_object)
}

# --- 9. Cell Type Annotation (SingleR) & CAF Identification ---
if (start_step <= 7) { # Step changed from 6 to 7
  perform_cell_type_annotation <- function(seurat_obj) {
    message("Performing cell type annotation using SingleR...")
    
    if (!requireNamespace("celldex", quietly = TRUE)) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
      }
      BiocManager::install("celldex")
    }
    library(celldex)
    
    hpca.se <- HumanPrimaryCellAtlasData()
    
    sce_obj <- as.SingleCellExperiment(seurat_obj)
    
    predictions <- SingleR(test = sce_obj,
                           ref = hpca.se,
                           labels = hpca.se$label.main)
    
    seurat_obj$singler_labels <- predictions$pruned.labels
    
    message("Cell type annotation complete. Adding to metadata and generating UMAP.")
    
    output_filename_singler_umap <- "UMAP_SingleR_CellTypes.pdf"
    pdf(output_filename_singler_umap, width = 10, height = 8)
    p_singler <- DimPlot(seurat_obj, reduction = "umap", group.by = "singler_labels", label = TRUE, repel = TRUE) +
      ggtitle("UMAP by SingleR Cell Types") +
      guides(color = guide_legend(override.aes = list(size = 3)))
    print(add_umap_arrows_and_theme_void(p_singler, seurat_obj, show_legend = TRUE))
    dev.off()
    message(paste0("UMAP plot with SingleR cell types saved to '", output_filename_singler_umap, "'."))
    
    message("Identifying Cancer-Associated Fibroblasts (CAFs) based on marker expression...")
    caf_markers <- c("ACTA2", "COL1A1", "DCN", "FAP", "PDGFRA", "THY1", "VIM")
    
    # Filter for markers present in the Seurat object after gene conversion
    valid_caf_markers <- intersect(caf_markers, rownames(seurat_obj))
    if (length(valid_caf_markers) == 0) {
      warning("No common CAF markers found after gene symbol conversion. Skipping CAF scoring.")
      seurat_obj$CAF_score <- NA
      seurat_obj$is_CAF <- "Unknown"
    } else {
      seurat_obj$CAF_score <- colMeans(x = GetAssayData(seurat_obj, assay = "RNA", layer = "data")[valid_caf_markers, , drop = FALSE])
      seurat_obj$is_CAF <- ifelse(seurat_obj$CAF_score > quantile(seurat_obj$CAF_score, 0.75), "CAF", "Non-CAF")
      seurat_obj$is_CAF <- factor(seurat_obj$is_CAF, levels = c("Non-CAF", "CAF"))
      
      output_filename_caf_umap <- "UMAP_CAF_Status.pdf"
      pdf(output_filename_caf_umap, width = 10, height = 8)
      p_caf <- DimPlot(seurat_obj, reduction = "umap", group.by = "is_CAF", label = FALSE, cols = c("lightgrey", "darkgreen")) +
        ggtitle("UMAP by CAF Status") +
        guides(color = guide_legend(override.aes = list(size = 3)))
      print(add_umap_arrows_and_theme_void(p_caf, seurat_obj, show_legend = TRUE))
      dev.off()
      message(paste0("UMAP plot for CAF status saved to '", output_filename_caf_umap, "'."))
    }
    
    return(seurat_obj)
  }
  current_seurat_object <- perform_cell_type_annotation(current_seurat_object)
}

# --- 10. Differential Expression Analysis ---
if (start_step <= 8) { # Step changed from 7 to 8
  perform_differential_expression <- function(seurat_obj) {
    message("Performing Differential Expression (DE) analysis between clusters...")
    
    de_results_cluster_0 <- FindMarkers(seurat_obj, ident.1 = 0, min.pct = 0.25, logfc.threshold = 0.25)
    write.csv(de_results_cluster_0, "DE_Cluster_0_vs_All.csv", row.names = TRUE)
    message("DE analysis for Cluster 0 vs. all other cells saved to 'DE_Cluster_0_vs_All.csv'.")
    
    pdf("VolcanoPlot_Cluster_0.pdf", width = 10, height = 10)
    print(EnhancedVolcano(de_results_cluster_0,
                          lab = rownames(de_results_cluster_0),
                          x = 'avg_log2FC',
                          y = 'p_val_adj',
                          title = 'Cluster 0 vs. All Other Cells',
                          pCutoff = 0.05,
                          FCcutoff = 0.5,
                          pointSize = 3.0,
                          labSize = 4.0))
    dev.off()
    message("Volcano plot for Cluster 0 DE genes saved to 'VolcanoPlot_Cluster_0.pdf'.")
    
    message("Differential Expression analysis complete.")
    return(de_results_cluster_0)
  }
  de_analysis_results <- perform_differential_expression(current_seurat_object)
}

# --- INTERMEDIATE DATA SAVING (After Differential Expression) ---
if (start_step <= 8) { # Step changed from 7 to 8
  message("\n--- INTERMEDIATE DATA SAVING (after Differential Expression) ---")
  output_filename_intermediate <- "seurat_object_after_DE.rds"
  saveRDS(current_seurat_object, file = output_filename_intermediate)
  message(paste0("Seurat object saved to '", output_filename_intermediate, "'."))
} else {
  message("Skipping intermediate data saving after DE (already past this step).")
}

# --- 11. Malignancy Detection (scATOMIC Classifier) ---
if (start_step <= 9) { # Step changed from 8 to 9
  message("Running scATOMIC malignancy detection...")
  current_seurat_object <- perform_scatomic_malignancy_detection(current_seurat_object)
  message("scATOMIC malignancy detection step completed.")
}

# --- 12. Cell Cycle Scoring ---
if (start_step <= 10) { # Step changed from 9 to 10
  message("Performing cell cycle scoring...")
  
  s.genes <- Seurat::cc.genes$s.genes
  g2m.genes <- Seurat::cc.genes$g2m.genes
  
  current_seurat_object <- CellCycleScoring(current_seurat_object,
                                            s.features = s.genes,
                                            g2m.features = g2m.genes,
                                            set.ident = FALSE)
  
  message("Cell cycle scoring complete. 'S.Score', 'G2M.Score', and 'Phase' added to metadata.")
  
  message("Ensuring 'Phase' metadata is a factor for consistent plotting...")
  phase_order <- c("G1", "S", "G2M")
  current_seurat_object$Phase <- factor(current_seurat_object$Phase, levels = phase_order)
  
  user_regress_choice <- tolower(readline(prompt = "Do you want to regress out cell cycle effects? (yes/no): "))
  
  if (user_regress_choice == "yes") {
    message("Regressing out cell cycle effects...")
    current_seurat_object <- ScaleData(current_seurat_object,
                                       vars.to.regress = c("S.Score", "G2M.Score"),
                                       features = rownames(current_seurat_object))
    message("Cell cycle effects regressed out. Consider re-running PCA and UMAP for optimal results.")
    
    user_re_run_dimred <- tolower(readline(prompt = "Cell cycle effects regressed. Re-run PCA, UMAP, and clustering now? (yes/no): "))
    if (user_re_run_dimred == "yes") {
      message("Re-running PCA, UMAP, and clustering after cell cycle regression...")
      current_seurat_object <- RunPCA(current_seurat_object, features = VariableFeatures(object = current_seurat_object))
      current_seurat_object <- RunUMAP(current_seurat_object, dims = 1:30)
      current_seurat_object <- FindNeighbors(current_seurat_object, dims = 1:30)
      current_seurat_object <- FindClusters(current_seurat_object, resolution = 0.5)
      message("Dimensionality reduction and clustering re-run.")
    } else {
      message("Skipping re-running dimensionality reduction and clustering.")
    }
    
  } else {
    message("Skipping cell cycle regression.")
  }
  
  message("Generating UMAP plots with cell cycle phases (standard Seurat plots)...")
  
  output_filename_cc_umap <- "UMAP_CellCycle_Phases.pdf"
  pdf(output_filename_cc_umap, width = 10, height = 8)
  p_cc_phase <- DimPlot(current_seurat_object, reduction = "umap", group.by = "Phase", label = FALSE, repel = TRUE) +
    ggplot2::ggtitle("UMAP by Cell Cycle Phase") +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 3)))
  print(add_umap_arrows_and_theme_void(p_cc_phase, current_seurat_object, show_legend = TRUE))
  dev.off()
  message(paste0("UMAP plot by cell cycle phases saved to '", output_filename_cc_umap, "'."))
  
  output_filename_s_score_umap <- "UMAP_S_Score.pdf"
  pdf(output_filename_s_score_umap, width = 10, height = 8)
  p_s_score <- FeaturePlot(current_seurat_object, features = "S.Score", reduction = "umap") +
    ggplot2::ggtitle("UMAP by S Phase Score")
  print(add_umap_arrows_and_theme_void(p_s_score, current_seurat_object, show_legend = TRUE))
  dev.off()
  message(paste0("UMAP plot by S.Score saved to '", output_filename_s_score_umap, "'."))
  
  output_filename_g2m_score_umap <- "UMAP_G2M_Score.pdf"
  pdf(output_filename_g2m_score_umap, width = 10, height = 8)
  p_g2m_score <- FeaturePlot(current_seurat_object, features = "G2M.Score", reduction = "umap") +
    ggplot2::ggtitle("UMAP by G2M Phase Score")
  print(add_umap_arrows_and_theme_void(p_g2m_score, current_seurat_object, show_legend = TRUE))
  dev.off()
  message(paste0("UMAP plot by G2M.Score saved to '", output_filename_g2m_score_umap, "'."))
  message("Cell cycle scoring step completed.")
}

# --- 13. Generate SCP-like Cell Cycle Plots ---
if (start_step <= 11) { # Step changed from 10 to 11
  message("Generating UMAP plots with cell cycle phases using SCP-like functions...")
  
  if (!("Phase" %in% colnames(current_seurat_object@meta.data) &&
        "S.Score" %in% colnames(current_seurat_object@meta.data) &&
        "G2M.Score" %in% colnames(current_seurat_object@meta.data))) {
    stop("Cell cycle scores (Phase, S.Score, G2M.Score) not found in Seurat object metadata. Please run step 10 first.")
  }
  
  output_filename_cc_phase_scp <- "UMAP_CellCycle_Phases_SCP.pdf"
  pdf(output_filename_cc_phase_scp, width = 10, height = 8)
  p_cc_phase_scp <- CellDimPlot(current_seurat_object,
                                group.by = "seurat_clusters",
                                stat.by = "Phase",
                                reduction = "umap",
                                label = TRUE,
                                pt.size = 0.5,
                                stat_plot_size = 0.08,
                                stat_plot_label = FALSE) +
    ggplot2::ggtitle("UMAP of Cell Cycle Phase Distribution within Clusters (SCP)")
  print(add_umap_arrows_for_scp(p_cc_phase_scp, current_seurat_object, reduction_name = "umap", show_legend = TRUE))
  dev.off()
  message(paste0("UMAP plot by cell cycle phases (SCP style with phase distribution) saved to '", output_filename_cc_phase_scp, "'."))
  
  
  output_filename_s_score_scp <- "UMAP_S_Score_SCP.pdf"
  pdf(output_filename_s_score_scp, width = 10, height = 8)
  p_s_score_scp <- FeatureDimPlot(current_seurat_object,
                                  variable = "S.Score",
                                  reduction = "umap",
                                  pt.size = 0.5,
                                  palette = "Spectral") +
    ggplot2::ggtitle("UMAP by S Phase Score (SCP)")
  print(add_umap_arrows_for_scp(p_s_score_scp, current_seurat_object, reduction_name = "umap", show_legend = TRUE))
  dev.off()
  message(paste0("UMAP plot by S.Score (SCP style) saved to '", output_filename_s_score_scp, "'."))
  
  
  output_filename_g2m_score_scp <- "UMAP_G2M_Score_SCP.pdf"
  pdf(output_filename_g2m_score_scp, width = 10, height = 8)
  p_g2m_score_scp <- FeatureDimPlot(current_seurat_object,
                                    variable = "G2M.Score",
                                    reduction = "umap",
                                    pt.size = 0.5,
                                    palette = "Spectral") +
    ggplot2::ggtitle("UMAP by G2M Phase Score (SCP)")
  print(add_umap_arrows_for_scp(p_g2m_score_scp, current_seurat_object, reduction_name = "umap", show_legend = TRUE))
  dev.off()
  message(paste0("UMAP plot by G2M.Score (SCP style) saved to '", output_filename_g2m_score_scp, "'."))
  
  message("All cell cycle plots generated using SCP-like functions and saved.")
}

# --- 14. Interactive Gene Expression Analysis & Visualization ---
if (start_step <= 12) { # Step changed from 11 to 12
  message("Starting interactive gene expression analysis...")
  
  interactive_gene_analysis <- function(current_seurat_object) {
    while (TRUE) {
      gene_input <- readline(prompt = "Enter gene symbol(s) (comma-separated, e.g., EPCAM, PTPRC) or 'done' to finish: ")
      if (tolower(gene_input) == "done") {
        message("Exiting interactive gene expression analysis.")
        break
      }
      
      target_genes <- str_trim(unlist(str_split(gene_input, ",")))
      target_genes <- intersect(target_genes, rownames(current_seurat_object))
      
      if (length(target_genes) == 0) {
        message("No valid genes found. Please try again.")
        next
      }
      
      message(paste0("Analyzing gene(s): ", paste(target_genes, collapse = ", ")))
      
      save_plots_choice <- tolower(readline(prompt = "Save plots to PDF (yes/no)? ")) == "yes"
      
      for (target_gene in target_genes) {
        plot_title_feature <- paste0("UMAP of Continuous ", target_gene, " Expression")
        if (save_plots_choice) {
          output_filename_feature <- paste0("UMAP_", target_gene, "_FeaturePlot.pdf")
          pdf(output_filename_feature, width = 8, height = 7)
        }
        
        p_feature <- FeaturePlot(current_seurat_object,
                                 features = target_gene,
                                 reduction = "umap",
                                 pt.size = 0.5,
                                 order = TRUE) +
          ggtitle(plot_title_feature)
        print(add_umap_arrows_and_theme_void(p_feature, current_seurat_object, show_legend = TRUE))
        
        if (save_plots_choice) {
          dev.off()
          message(paste0("UMAP feature plot for '", target_gene, "' continuous expression saved to '", output_filename_feature, "'."))
        }
        
        threshold_value_input <- readline(prompt = paste0("Enter expression threshold for ", target_gene, " (e.g., 0.5) or leave blank to skip binary plot: "))
        if (nchar(threshold_value_input) > 0) {
          threshold_value <- as.numeric(threshold_value_input)
          if (is.na(threshold_value)) {
            message("Invalid threshold value. Skipping binary plot for this gene.")
            next
          }
          
          gene_status_vector <- as.factor(
            ifelse(GetAssayData(current_seurat_object, assay = "RNA", layer = "data")[target_gene, ] > threshold_value,
                   paste0(target_gene, "_Positive"), paste0(target_gene, "_Negative"))
          )
          names(gene_status_vector) <- colnames(current_seurat_object)
          
          current_seurat_object <- AddMetaData(
            object = current_seurat_object,
            metadata = gene_status_vector,
            col.name = "gene_expression_status"
          )
          
          plot_title_binary <- paste0("UMAP of ", target_gene, " Expression (Threshold > ", threshold_value, ")")
          if (save_plots_choice) {
            output_filename_binary <- paste0("UMAP_", target_gene, "_Positive_Negative_Plot.pdf")
            pdf(output_filename_binary, width = 8, height = 7)
          }
          
          p_binary <- DimPlot(current_seurat_object,
                              reduction = "umap",
                              group.by = "gene_expression_status",
                              cols = c("lightgrey", "red"),
                              label = FALSE,
                              repel = TRUE) +
            ggtitle(plot_title_binary)
          print(add_umap_arrows_and_theme_void(p_binary, current_seurat_object, show_legend = TRUE))
          
          if (save_plots_choice) {
            dev.off()
            message(paste0("UMAP plot for '", target_gene, "' expression saved to '", output_filename_binary, "'."))
          }
        }
      }
    }
  }
  interactive_gene_analysis(current_seurat_object)
}
# --- 15. Trajectory Inference (Monocle3) ---
if (start_step <= 13) { # Step changed from 12 to 13
  
  message("Starting Trajectory Inference using Monocle3...")
  
  # Define helper function to load CDS if exists
  load_cds_if_exists <- function(filename) {
    if (file.exists(filename)) {
      message(paste0("Found existing Monocle3 CDS file: ", filename, ". Loading..."))
      return(readRDS(filename))
    } else {
      message(paste0("No existing CDS file found at: ", filename))
      return(NULL)
    }
  }
  
  # Initial CDS filename (you can adjust or add dynamic naming after root selection)
  cds_filename <- "monocle3_cds_object.rds"
  
  # Try to load existing CDS
  cds <- load_cds_if_exists(cds_filename)
  
  if (is.null(cds)) {
    # CDS not found — run full Monocle3 analysis
    
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager")
    }
    
    required_monocle_packages <- c(
      "BiocGenerics", "DelayedArray", "DelayedMatrixStats", "limma", "lme4",
      "S4Vectors", "SingleCellExperiment", "SummarizedExperiment", "sparseMatrixStats",
      "IRanges", "GenomicRanges", "ggplot2", "dplyr", "igraph", "patchwork", "Matrix"
    )
    
    for (pkg in required_monocle_packages) {
      if (!requireNamespace(pkg, quietly = TRUE)) {
        message(paste0("Installing missing package: ", pkg))
        if (pkg %in% c("monocle3", "DelayedArray", "DelayedMatrixStats", "SingleCellExperiment", "SummarizedExperiment", "sparseMatrixStats", "IRanges", "GenomicRanges")) {
          BiocManager::install(pkg, update = FALSE, ask = FALSE)
        } else {
          install.packages(pkg, dependencies = TRUE)
        }
      }
    }
    
    if (!requireNamespace("monocle3", quietly = TRUE)) {
      message("Installing monocle3...")
      BiocManager::install("monocle3", update = FALSE, ask = FALSE)
    }
    library(monocle3)
    message("Monocle3 and its dependencies loaded successfully.")
    
    message("Preprocessing 'scATOMIC_prediction' column for Monocle3...")
    
    current_seurat_object@meta.data$scATOMIC_prediction <- as.character(current_seurat_object@meta.data$scATOMIC_prediction)
    
    problematic_string_scatomic <- "source(\"~/.active−rstudio−document\", echo = TRUE)"
    unknown_category <- "Unknown_scATOMIC_Cell"
    unassigned_category <- "Unassigned_scATOMIC"
    
    current_seurat_object@meta.data$scATOMIC_prediction[current_seurat_object@meta.data$scATOMIC_prediction == problematic_string_scatomic] <- unknown_category
    current_seurat_object@meta.data$scATOMIC_prediction[is.na(current_seurat_object@meta.data$scATOMIC_prediction)] <- unassigned_category
    
    all_levels <- sort(unique(c(current_seurat_object@meta.data$scATOMIC_prediction, unknown_category, unassigned_category)))
    current_seurat_object$scATOMIC_prediction <- factor(current_seurat_object@meta.data$scATOMIC_prediction,
                                                        levels = all_levels)
    message("'scATOMIC_prediction' column processed and converted to factor.")
    
    message("Converting Seurat object to Monocle3 cell_data_set (CDS)...")
    
    # Ensure you are using the correct assay and layer for expression data
    expression_matrix <- GetAssayData(current_seurat_object, assay = "RNA", layer = "counts")
    cell_metadata <- current_seurat_object@meta.data
    gene_metadata <- data.frame(gene_short_name = rownames(expression_matrix), row.names = rownames(expression_matrix))
    
    if (!all(rownames(cell_metadata) == colnames(expression_matrix))) {
      warning("Cell metadata rownames do not match expression matrix column names. Attempting reorder.")
      cell_metadata <- cell_metadata[colnames(expression_matrix), ]
    }
    if (!all(rownames(gene_metadata) == rownames(expression_matrix))) {
      warning("Gene metadata rownames do not match expression matrix row names. Attempting reorder.")
      gene_metadata <- gene_metadata[rownames(expression_matrix), ]
    }
    
    cds <- new_cell_data_set(expression_matrix,
                             cell_metadata = cell_metadata,
                             gene_metadata = gene_metadata)
    
    message("Adding Seurat UMAP embeddings to CDS...")
    umap_coords <- Embeddings(current_seurat_object, reduction = "umap")
    colnames(umap_coords) <- c("UMAP_1", "UMAP_2")
    
    reducedDims(cds)$UMAP <- umap_coords
    
    message("CDS object created and Seurat UMAP embeddings transferred.")
    
    message("Running Monocle3 clustering...")
    cds <- cluster_cells(cds, reduction_method = "UMAP", k = 20, resolution = 1e-6)
    message("Monocle3 clustering complete.")
    
    message("Learning cell trajectory graph and ordering cells with Monocle3...")
    
    cds <- learn_graph(cds, use_partition = FALSE)
    
  } else {
    library(monocle3) # Ensure monocle3 is loaded if CDS was loaded
  }
  
  # --- Interactive Root Selection for Monocle3 ---
  unique_scatomic_labels <- levels(factor(colData(cds)$scATOMIC_prediction))
  message("\n--- Available scATOMIC Cell Types for Monocle3 Root Selection ---")
  print(unique_scatomic_labels)
  message("-----------------------------------------------------")
  
  selected_root_label <- readline(prompt = "Enter the full name of the scATOMIC cell type to use as the trajectory root (e.g., 'Normal Fibroblast' or 'Cholangiocarcinoma'). Type 'auto' for automatic detection: ")
  
  selected_root_label <- trimws(selected_root_label)
  
  if (tolower(selected_root_label) == "auto") {
    message("Automatically detecting root node(s) for trajectory inference.")
    cds <- order_cells(cds)
    root_label_for_filename <- "AutoRoot"
  } else {
    message(paste0("Attempting to set '", selected_root_label, "' as the root node for trajectory inference."))
    
    root_cells <- colData(cds) %>%
      as.data.frame() %>%
      dplyr::filter(scATOMIC_prediction == selected_root_label) %>%
      rownames()
    
    if (length(root_cells) == 0) {
      warning(paste0("No cells found with label '", selected_root_label, "'. Reverting to automatic root detection."))
      cds <- order_cells(cds, root_cells = NULL)
      root_label_for_filename <- "AutoRoot"
    } else {
      cds <- order_cells(cds, root_cells = root_cells)
      root_label_for_filename <- gsub("[^A-Za-z0-9_]", "", selected_root_label)
    }
  }
  
  message("Cells ordered along trajectory.")
  
  # Save cds with root label in filename
  cds_filename_root <- paste0("monocle3_cds_object_", root_label_for_filename, "Root.rds")
  saveRDS(cds, cds_filename_root)
  message(paste0("Monocle3 CDS object saved as '", cds_filename_root, "'."))
  
  # Plot and save trajectory plots to PDF (with print() to fix blank pdf issue)
  output_filename_monocle_traj <- paste0("UMAP_Monocle3_Trajectory_", root_label_for_filename, "Root.pdf")
  pdf(output_filename_monocle_traj, width = 10, height = 8)
  
  print(
    plot_cells(cds,
               color_cells_by = "pseudotime",
               label_groups_by_cluster = FALSE,
               label_leaves = TRUE,
               label_branch_points = TRUE,
               graph_label_size = 3) +
      ggplot2::ggtitle(paste0("Monocle3 Trajectory (Pseudotime - ", gsub("_", " ", root_label_for_filename), " Root)")) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))
  )
  
  print(
    plot_cells(cds,
               color_cells_by = "scATOMIC_prediction",
               label_groups_by_cluster = TRUE,
               label_leaves = TRUE,
               label_branch_points = TRUE,
               graph_label_size = 3) +
      ggplot2::ggtitle(paste0("Monocle3 Trajectory (scATOMIC Prediction - ", gsub("_", " ", root_label_for_filename), " Root)")) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))
  )
  
  dev.off()
  message(paste0("Monocle3 trajectory plots saved to '", output_filename_monocle_traj, "'."))
  
  message("Trajectory Inference step completed.")
  
}

# --- 16. Cell-Cell Communication Analysis (CellChat) ---
if (start_step <= 14) { # Step changed from 16 to 14 in previous output snippet, confirming this is the correct block
  
  message("Running CellChat analysis...")
  
  # Check if CellChat object already exists
  cellchat_file <- file.path(data_dir, "cellchat_object.Rdata")
  if (file.exists(cellchat_file)) {
    reuse <- tolower(readline(prompt = "CellChat object found. Do you want to reuse it? (yes/no): "))
    if (reuse == "yes") {
      load(cellchat_file)
      message("Loaded existing CellChat object.")
    } else {
      rm(cellchat)  # Remove any previously loaded object
    }
  }
  
  # If cellchat doesn't exist or was cleared, run new analysis
  if (!exists("cellchat")) {
    
    # Step 1: Prepare protein-coding gene list
    protein_coding_genes <- keys(org.Hs.eg.db, keytype = "SYMBOL")
    
    # Step 2: Extract normalized data
    data.input <- GetAssayData(current_seurat_object, assay = "RNA", slot = "data")
    data.input <- data.input[rownames(data.input) %in% protein_coding_genes, ]
    
    # Step 3: Extract metadata and sanitize cluster labels
    meta <- current_seurat_object@meta.data
    meta$seurat_clusters <- as.character(meta$seurat_clusters)
    meta$seurat_clusters[meta$seurat_clusters == "0"] <- "0_"
    meta$seurat_clusters <- factor(meta$seurat_clusters)
    
    # Convert scATOMIC_prediction to factor BEFORE dropping levels
    current_seurat_object$scATOMIC_prediction <- as.factor(current_seurat_object$scATOMIC_prediction)
    
    # Drop unused levels from the scATOMIC_prediction factor in your Seurat object
    # This ensures only levels with actual cells are considered when creating CellChat object
    current_seurat_object$scATOMIC_prediction <- droplevels(current_seurat_object$scATOMIC_prediction)
    
    # Optional: Verify the levels and counts after dropping unused levels
    message("Levels in current_seurat_object$scATOMIC_prediction after dropping unused levels:")
    print(levels(current_seurat_object$scATOMIC_prediction))
    message("Counts in current_seurat_object$scATOMIC_prediction after dropping unused levels:")
    print(table(current_seurat_object$scATOMIC_prediction))
    
    # Step 4: Create CellChat object
    # The 'group.by' argument will now use the cleaned scATOMIC_prediction levels
    cellchat <- createCellChat(object = data.input, meta = meta, group.by = "scATOMIC_prediction")
    cellchat@DB <- CellChatDB.human
    
    # Step 5: Process cell-cell communication (these steps remain as they are)
    cellchat <- subsetData(cellchat)
    cellchat <- identifyOverExpressedGenes(cellchat)
    cellchat <- identifyOverExpressedInteractions(cellchat)
    cellchat <- computeCommunProb(cellchat)
    # TEMPORARY FIX: Set min.cells to 0 to bypass potential filtering issues that lead to zero cell types
    cellchat <- filterCommunication(cellchat, min.cells = 0) # CHANGED from 10 to 0
    cellchat <- computeCommunProbPathway(cellchat)
    cellchat <- aggregateNet(cellchat)
    cellchat <- netAnalysis_computeCentrality(cellchat)
    
    # Save CellChat object for reuse
    save(cellchat, file = cellchat_file)
    message("New CellChat object created and saved.")
  }
  
  # --- CRITICAL DIAGNOSTIC: Check cellchat@idents levels before color generation ---
  message("Levels of cellchat@idents BEFORE color generation:")
  print(levels(cellchat@idents))
  message("Number of cell types BEFORE color generation:")
  print(length(levels(cellchat@idents)))
  # --- END CRITICAL DIAGNOSTIC ---
  
  # --- CRITICAL FIX: Generate cell_colors ALWAYS after cellchat object is loaded/created ---
  # This ensures 'cell_colors' is available for visualization regardless of reuse option.
  num_cell_types <- length(levels(cellchat@idents))
  # NEW FIX: Use grDevices::rainbow as an alternative to scales::hue_pal() to generate colors
  cell_colors <- grDevices::rainbow(num_cell_types) # CHANGED

  # Keep setting as an attribute on the idents factor, as this is a common practice and doesn't error.
  attr(cellchat@idents, "colors") <- cell_colors
  cellchat@colors <- cell_colors # ADD THIS LINE
  # --- END CRITICAL FIX ---
  
  
  # --- Visualization ---
  # NEW STRATEGY for netVisual_circle: Remove 'color.use' to bypass igraph error for now.
  # If this runs successfully, it confirms the issue is specifically with color application.
  pdf(file.path(data_dir, "cellchat_netVisual_circle.pdf"), width = 8, height = 8)
  print(netVisual_circle(cellchat@net$weight,
                         vertex.weight = as.numeric(table(cellchat@idents)),
                         weight.scale = TRUE,
                         label.edge = FALSE,
                         vertex.label.color = "black",
                         color.use = cell_colors)) # ADDED color.use to explicitly pass your generated colors
  dev.off()
  
  # ---- 1. Signaling Role Heatmap ----
  ht <- netAnalysis_signalingRole_heatmap(cellchat)
  ht@row_names_param$gp <- grid::gpar(fontsize = 3)
  ht@row_names_param$max_width <- grid::unit(6, "cm")
  
  pdf(file.path(data_dir, "cellchat_netAnalysis_signalingRole_heatmap.pdf"), width = 8, height = 12)
  ComplexHeatmap::draw(ht, heatmap_legend_side = "right")
  dev.off()
  
  # ---- 2. Bubble Plot for Selected Pathway ----
  message("Available pathways:")
  print(cellchat@netP$pathways)
  selected_pathway <- readline(prompt = "Enter the name of a pathway to visualize (or press Enter to skip): ")
  
  if (nzchar(selected_pathway) && selected_pathway %in% cellchat@netP$pathways) {
    pdf(file.path(data_dir, paste0("cellchat_netVisual_bubble_", selected_pathway, ".pdf")), width = 15, height = 12)
    # FIX: Removed 'color.use' as it's an unused argument for netVisual_bubble based on documentation
    p <- netVisual_bubble(cellchat, signaling = selected_pathway) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 3, vjust = 1, margin = ggplot2::margin(t = 2)),
        axis.text.y = ggplot2::element_text(size = 6)
      )
    print(p)
    dev.off()
    
    # ---- 3. Chord Plot for Selected Pathway ----
    strwidth <- function(x) {0.02}
    pdf(file.path(data_dir, paste0("cellchat_netVisual_chord_", selected_pathway, ".pdf")), width = 14, height = 14)
    circlize::circos.clear()
    circlize::circos.par(gap.degree = 5, track.margin = c(0.01, 0.01), start.degree = 90, track.height = 0.05)
    graphics::par(cex = 0.4, mar = c(1, 1, 1, 1))
    # FIX: Removed 'group.colors' as it's an unused argument for netVisual_chord_cell
    netVisual_chord_cell(cellchat, signaling = selected_pathway, lab.cex = 3) # Removed group.colors
    dev.off()
    
  } else if (nzchar(selected_pathway)) {
    warning("Selected pathway not found in available pathways.")
  } else {
    message("No pathway selected, skipping pathway visualizations.")
  }
  
  message("CellChat analysis and visualizations complete.")
}
# --- 17. Final Data Saving ---
if (start_step <= 15) { # Step changed from 14 to 15
  user_save_choice <- tolower(readline(prompt = "Do you want to save the final Seurat object? (yes/no): "))
  
  if (user_save_choice == "yes") {
    message("Saving final Seurat object...")
    saveRDS(current_seurat_object, "seurat_object_final.rds")
    message("Final Seurat object saved as 'seurat_object_final.rds'.")
    message("Pipeline complete. You can now close RStudio.")
  } else {
    message("Skipping final Seurat object save. Pipeline complete.")
  }
}