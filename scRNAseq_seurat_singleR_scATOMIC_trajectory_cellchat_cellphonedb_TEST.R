# A Comprehensive and Reproducible R Pipeline for Single-Cell RNA Sequencing Data Analysis

# Set RETICULATE_CONDA_PATH and disable auto-installation
options(reticulate.conda_auto_install = FALSE)
options(reticulate.install_miniforge = FALSE)

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
library(reticulate) # Already loaded, but good to ensure
library(remotes) # For installing from GitHub
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
library(tools) # For file_path_sans_ext, basename etc.
library(monocle3) # For trajectory analysis
library(scales) # ADDED: For generating color palettes (used in CellChat color assignment)
library(grDevices) # ADDED: For rainbow() function
library(ktplots)
library(scatterpie)

# --- NEW: Install ktplots if not already installed (from zktuong/ktplots GitHub) ---
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}
library(devtools)

if (!requireNamespace("ktplots", quietly = TRUE)) {
  message("Installing ktplots from GitHub (zktuong/ktplots)... This may take a moment.")
  devtools::install_github("zktuong/ktplots")
}
library(ktplots)


message("\n--- All required libraries loaded successfully. ---")

# --- Python Environment Selection at Start ---
cat("Please select the primary Python environment for this R session:\n")
cat("1. scATOMIC / Rmagic (C:/miniconda3/envs/r_magic_env/python.exe)\n")
cat("2. CellPhoneDB (C:/miniconda3/envs/cellphonedb_env/python.exe)\n")

env_choice <- as.integer(readline(prompt = "Enter your choice (1 or 2): "))

if (env_choice == 1) {
  python_env_path <- "C:/miniconda3/envs/r_magic_env/python.exe"
  assign("CELLPHONEDB_PYTHON_PATH", NULL, envir = .GlobalEnv)
  cat("\n--- Initializing R session with scATOMIC/Rmagic Python environment ---\n")
} else if (env_choice == 2) {
  python_env_path <- "C:/miniconda3/envs/cellphonedb_env/python.exe"
  assign("CELLPHONEDB_PYTHON_PATH", python_env_path, envir = .GlobalEnv)
  cat("\n--- Initializing R session with CellPhoneDB Python environment ---\n")
} else {
  stop("Invalid choice. Please restart your R session and enter 1 or 2.")
}

reticulate::use_python(python = python_env_path, required = TRUE)

# Verify the selected Python environment
cat("Verifying Python environment configuration:\n")
py_config()

# --- Custom Helper Function for UMAP Plot Styling with Direct Arrows ---
# This function adds UMAP dimension arrows and applies a theme_void for clean plots.
add_umap_arrows_and_theme_void <- function(base_plot, seurat_obj_for_coords,
                                           arrow_length_prop = 0.05, arrow_offset_prop = 0.02,
                                           umap_dim_names = c("UMAP1", "UMAP2"),
                                           show_legend = TRUE) {
  
  if (!inherits(base_plot, "ggplot")) {
    stop("Input 'base_plot' must be a ggplot object.")
  }
  
  # Extract UMAP coordinates
  umap_coords <- Embeddings(seurat_obj_for_coords, reduction = "umap")
  min_x <- min(umap_coords[,1], na.rm = TRUE)
  max_x <- max(umap_coords[,1], na.rm = TRUE)
  min_y <- min(umap_coords[,2], na.rm = TRUE)
  max_y <- max(umap_coords[,2], na.rm = TRUE)
  
  range_x <- max_x - min_x
  range_y <- max_y - min_y
  
  # Calculate effective arrow length and offset based on plot dimensions
  effective_arrow_length <- min(range_x, range_y) * arrow_length_prop
  effective_arrow_offset <- min(range_x, range_y) * arrow_offset_prop
  
  # Define arrow starting and ending points
  start_x_base <- min_x + effective_arrow_offset
  start_y_base <- min_y + effective_arrow_offset
  
  end_x_umap1 <- start_x_base + effective_arrow_length
  end_y_umap1 <- start_y_base
  
  end_x_umap2 <- start_x_base
  end_y_umap2 <- start_y_base + effective_arrow_length
  
  arrow_style <- ggplot2::arrow(length = grid::unit(0.08, "inches"), type = "closed", ends = "last")
  
  # Apply theme, add arrows and labels
  custom_plot <- base_plot +
    ggplot2::theme_void() + # Use theme_void for a clean background
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 14),
      legend.position = if (show_legend) "right" else "none",
      legend.box.just = "left",
      legend.direction = "vertical",
      legend.key.size = grid::unit(10, "pt"),
      legend.background = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(t = 5, r = 5, b = 15, l = 15, unit = "pt") # Adjust margins to make space for arrows
    ) +
    # UMAP1 arrow
    geom_segment(aes(x = start_x_base, y = start_y_base, xend = end_x_umap1, yend = end_y_umap1),
                 arrow = arrow_style, color = "black", linewidth = 0.6) +
    annotate("text", x = end_x_umap1 + effective_arrow_offset * 0.75,
             y = end_y_umap1, label = umap_dim_names[1],
             hjust = 0, vjust = 0.5, size = 3, fontface = "bold") +
    # UMAP2 arrow
    geom_segment(aes(x = start_x_base, y = start_y_base, xend = end_x_umap2, yend = end_y_umap2),
                 arrow = arrow_style, color = "black", linewidth = 0.6) +
    annotate("text", x = end_x_umap2,
             y = end_y_umap2 + effective_arrow_offset * 1.5,
             label = umap_dim_names[2],
             hjust = 0.5, vjust = 0, size = 3, fontface = "bold", angle = 90) +
    # Adjust plot limits to ensure arrows are visible and not clipped
    coord_cartesian(xlim = c(min_x - range_x * 0.05, max_x + range_x * 0.05),
                    ylim = c(min_y - range_y * 0.05, max_y + range_y * 0.05),
                    expand = FALSE,
                    clip = "off") # Important: allow drawing outside plot area
  
  return(custom_plot)
}

# Helper function for gene symbol conversion (for CellChat/CellphoneDB compatibility)
convert_gene_symbols <- function(seurat_obj) {
  message("Starting gene symbol conversion for CellChat/CellphoneDB compatibility...")
  current_genes <- rownames(seurat_obj)
  original_gene_count <- length(current_genes)
  
  # Try to remove version numbers from Ensembl-like IDs if present (e.g., ENSG00000123456.1 -> ENSG00000123456)
  stripped_genes <- gsub("\\.[0-9]+$", "", current_genes)
  
  # Map all stripped genes to SYMBOL using org.Hs.eg.db
  # CORRECTED: Changed keytype from "ENSEMBL" to "SYMBOL"
  mapped_symbols <- AnnotationDbi::mapIds(
    org.Hs.eg.db,
    keys = stripped_genes,
    column = "SYMBOL",
    keytype = "SYMBOL", # Changed from "ENSEMBL"
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
  
  # Identify unique new symbols
  unique_new_symbols <- unique(gene_map_df$mapped_symbol)
  message(paste0("Number of unique new gene symbols after initial mapping: ", length(unique_new_symbols)))
  
  # Create a new matrix for the RNA assay, aggregating duplicate gene symbols
  # Initialize with 0s
  new_counts <- matrix(0, nrow = length(unique_new_symbols), ncol = ncol(seurat_obj))
  rownames(new_counts) <- unique_new_symbols
  colnames(new_counts) <- colnames(seurat_obj)
  
  new_data <- matrix(0, nrow = length(unique_new_symbols), ncol = ncol(seurat_obj))
  rownames(new_data) <- unique_new_symbols
  colnames(new_data) <- colnames(seurat_obj)
  
  # Aggregate counts and data by new symbol
  for (symbol in unique_new_symbols) {
    # Find all original gene IDs that map to the current symbol
    original_genes_for_symbol <- gene_map_df$original_gene[gene_map_df$mapped_symbol == symbol]
    
    # Filter to only those present in the current Seurat assay
    valid_original_genes <- intersect(original_genes_for_symbol, rownames(GetAssayData(seurat_obj, assay = "RNA", slot = "counts")))
    
    if (length(valid_original_genes) > 0) {
      if (length(valid_original_genes) == 1) {
        new_counts[symbol, ] <- GetAssayData(seurat_obj, assay = "RNA", slot = "counts")[valid_original_genes, ]
        new_data[symbol, ] <- GetAssayData(seurat_obj, assay = "RNA", slot = "data")[valid_original_genes, ]
      } else {
        # Sum counts/data for genes that map to the same symbol
        new_counts[symbol, ] <- colSums(GetAssayData(seurat_obj, assay = "RNA", slot = "counts")[valid_original_genes, , drop=FALSE])
        new_data[symbol, ] <- colSums(GetAssayData(seurat_obj, assay = "RNA", slot = "data")[valid_original_genes, , drop=FALSE])
      }
    }
  }
  
  # Replace the RNA assay in the Seurat object
  # Create a new AssayObject from the aggregated counts
  new_rna_assay <- CreateAssayObject(counts = new_counts)
  # Manually set the 'data' slot as CreateAssayObject only initializes 'counts' and 'scale.data'
  new_rna_assay@data <- new_data 
  
  seurat_obj[["RNA"]] <- new_rna_assay
  
  message(paste0("Gene symbol conversion complete. Original genes: ", original_gene_count,
                 ", Genes after conversion: ", nrow(seurat_obj[["RNA"]]), " (unique symbols)."))
  
  return(seurat_obj)
}


# Function for scATOMIC-based malignancy detection and annotation
# Modified to handle problematic strings and use AddMetaData
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
    mc.cores = 7, # Adjust based on your system's cores
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
  
  # Ensure the 'scATOMIC_prediction_name' column (or similar) exists in the output
  if (!"scATOMIC_prediction_name" %in% colnames(final_scatomic_annotations)) {
    stop("scATOMIC::create_summary_matrix output did not contain 'scATOMIC_prediction_name' column. Check scATOMIC output structure.")
  }
  
  # Align cell barcodes between scATOMIC results and Seurat object
  common_cells_scatomic <- intersect(rownames(final_scatomic_annotations), colnames(seurat_obj))
  
  if (length(common_cells_scatomic) == 0) {
    stop("No common cells found between scATOMIC results and Seurat object. Cannot merge scATOMIC results.")
  }
  
  # Extract relevant predictions, ensuring order matches seurat_obj cells
  predictions_for_seurat_raw <- final_scatomic_annotations[common_cells_scatomic, "scATOMIC_prediction_name", drop = FALSE]
  malignancy_for_seurat_raw <- final_scatomic_annotations[common_cells_scatomic, "scATOMIC_malignancy_prediction", drop = FALSE]
  
  # Initialize new metadata columns in Seurat object to handle all cells
  seurat_obj$scATOMIC_prediction <- NA_character_
  seurat_obj$scATOMIC_malignancy <- NA_character_
  
  # Assign predictions to the correct cells in the Seurat object
  seurat_obj@meta.data[common_cells_scatomic, "scATOMIC_prediction"] <- predictions_for_seurat_raw
  seurat_obj@meta.data[common_cells_scatomic, "scATOMIC_malignancy"] <- malignancy_for_seurat_raw
  
  message("scATOMIC cell type predictions and malignancy scores added to Seurat object metadata.")
  
  # Convert to factor and handle potential problematic strings if they exist (from previous issues)
  # The problematic string might come from unparsed R code in some environments, ensure it's handled.
  problematic_string <- "source(\"~/.active−rstudio−document\", echo = TRUE)" # Keep for robustness
  
  if ("scATOMIC_prediction" %in% colnames(seurat_obj@meta.data)) {
    # Replace problematic string with a more descriptive label
    seurat_obj$scATOMIC_prediction[seurat_obj$scATOMIC_prediction == problematic_string] <- "Unknown_scATOMIC_Cell"
    # Convert to factor and re-level to remove NAs from factor levels if present
    seurat_obj$scATOMIC_prediction <- factor(seurat_obj$scATOMIC_prediction)
    seurat_obj$scATOMIC_prediction <- droplevels(seurat_obj$scATOMIC_prediction) # Remove unused levels
  }
  
  if ("scATOMIC_malignancy" %in% colnames(seurat_obj@meta.data)) {
    seurat_obj$scATOMIC_malignancy[seurat_obj$scATOMIC_malignancy == problematic_string] <- "Unknown_scATOMIC_Malignancy"
    seurat_obj$scATOMIC_malignancy <- factor(seurat_obj$scATOMIC_malignancy)
    seurat_obj$scATOMIC_malignancy <- droplevels(seurat_obj$scATOMIC_malignancy)
  }
  
  # Save UMAP plot by scATOMIC_prediction
  
  # Get counts for each scATOMIC_prediction group
  if (!is.null(seurat_object$scATOMIC_prediction)) { # Changed from seurat_obj to seurat_object
    group_counts <- table(seurat_object$scATOMIC_prediction) # Changed from seurat_obj to seurat_object
    
    cat("\nAvailable scATOMIC_prediction groups:\n")
    for (i in seq_along(group_counts)) {
      cat(sprintf("%d. %s (%d cells)\n", i, names(group_counts)[i], group_counts[i]))
    }
  } else {
    stop("scATOMIC_prediction labels not found in Seurat object. Please ensure scATOMIC analysis was run correctly.")
  }
  
  # Interactively get user input for cell types to plot
  selected_indices_input <- readline(prompt = "Enter the numbers of cell types to plot (e.g., '1,3,5' or '2-4'): ")
  
  # Parse the input
  selected_indices <- unique(unlist(sapply(strsplit(selected_indices_input, ",")[[1]], function(x) {
    if (grepl("-", x)) {
      range_vals <- as.numeric(strsplit(x, "-")[[1]])
      return(seq(min(range_vals), max(range_vals)))
    } else {
      return(as.numeric(x))
    }
  })))
  
  # Validate input indices
  if (any(selected_indices < 1) || any(selected_indices > length(group_counts))) {
    stop("Invalid selection. Please enter numbers corresponding to the listed groups.")
  }
  
  # Get the names of the selected cell types
  selected_cell_types <- names(group_counts)[selected_indices]
  cat("\nSelected cell types for UMAP plot:", paste(selected_cell_types, collapse = ", "), "\n")
  
  # Create a subset of the Seurat object for plotting
  seurat_object_filtered <- subset(seurat_object, subset = scATOMIC_prediction %in% selected_cell_types) # Changed from seurat_obj to seurat_object
  
  # Ensure that the factor levels are dropped for the filtered object to avoid displaying unselected levels
  seurat_object_filtered$scATOMIC_prediction <- factor(seurat_object_filtered$scATOMIC_prediction)
  seurat_object_filtered$scATOMIC_prediction <- droplevels(seurat_object_filtered$scATOMIC_prediction)
  
  # Generate the UMAP plot using the filtered object and selected groups
  plot_scatomic_umap <- DimPlot(seurat_object_filtered, # Changed from seurat_obj_filtered to seurat_object_filtered
                                reduction = "umap",
                                group.by = "scATOMIC_prediction",
                                label = TRUE,
                                repel = TRUE,
                                pt.size = 0.5) +
    ggtitle("UMAP of Selected scATOMIC Predicted Cell Types") +
    NoLegend() # Optionally remove legend if labels are sufficient, or keep if preferred
  
  # Call your custom function for arrows and theme (assuming add_umap_arrows_and_theme_void is defined elsewhere)
  plot_scatomic_umap <- add_umap_arrows_and_theme_void(plot_scatomic_umap, seurat_object_filtered) # Changed from seurat_obj_filtered to seurat_object_filtered
  
  # Save the plot
  plot_filename_umap <- file.path(output_dir, "UMAP_Selected_scATOMIC_Predicted_CellTypes.pdf")
  ggsave(plot_filename_umap, plot_scatomic_umap, width = 10, height = 8)
  cat(paste0("UMAP plot for selected cell types saved to: '", plot_filename_umap, "'\n"))
  
  message("scATOMIC detection step completed.")
  
  return(seurat_obj)
}

# --- SCP-like Plotting Functions (from SCP_plots.R) ---

# Operator to provide a default value if 'a' is NULL
`%||%` <- function(a, b) {
  if (!is.null(a)) a else b
}

# Function to convert R color names to hex codes
col2hex <- function(color_vector) {
  sapply(color_vector, function(color) {
    if (is.character(color) && length(color) == 1 && nchar(color) %in% c(7, 9) && substr(color, 1, 1) == "#") {
      return(color) # Already hex
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

# Function to determine a default dimensionality reduction if not specified
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
  preferred_reductions <- c("umap", "tsne", "pca") # Order of preference
  
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

# Internal list of predefined color palettes
palette_list_internal <- list(
  "Paired" = c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928"),
  "Spectral" = c("#9E0142", "#D53E4F", "#F46D43", "#FDAE61", "#FEE08B", "#FFFFBF", "#E6F598", "#ABDDA4", "#66C2A5", "#3288BD", "#5E4FA2"),
  "jet" = c("#00007F", "#0000FF", "#007FFF", "#00FFFF", "#7FFF7F", "#FFFF00", "#FF7F00", "#FF0000", "#7F0000"),
  "Dark2" = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666"),
  "Set1" = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999"),
  "Greys" = c("#FFFFFF", "#F0F0F0", "#D9D9D9", "#BDBDBD", "#969696", "#737373", "#525252", "#252525", "#000000"),
  "RdYlBu" = c("#A50026", "#D73027", "#F46D43", "#FDAE61", "#FEE090", "#FFFFBF", "#E0F3F8", "#ABDDA4", "#66C2A5", "#3288BD", "#5E4FA2")
)
# Assign types to palettes for `palette_scp` function
attr(palette_list_internal$Spectral, "type") <- "continuous"
attr(palette_list_internal$jet, "type") <- "continuous"
attr(palette_list_internal$Paired, "type") <- "discrete"
attr(palette_list_internal$Dark2, "type") <- "discrete"
attr(palette_list_internal$Set1, "type") <- "discrete"
attr(palette_list_internal$Greys, "type") <- "continuous"
attr(palette_list_internal$RdYlBu, "type") <- "continuous"

# Flexible color palette function
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
  if (is.list(palcolor)) { # Flatten list if palcolor is a list
    palcolor <- unlist(palcolor)
  }
  if (all(palcolor == "")) { # If palcolor is empty, use default from list
    palcolor <- palette_list[[palette]]
  }
  if (is.null(palcolor) || length(palcolor) == 0) { # If still empty, use default
    palcolor <- palette_list[[palette]]
  }
  if (!is.null(names(palcolor))) { # If named colors are provided, match to x
    if (all(x %in% names(palcolor))) {
      palcolor <- palcolor[intersect(names(palcolor), x)]
    }
  }
  pal_n <- length(palcolor) # Number of colors in the chosen palette
  
  if (!type %in% c("auto", "discrete", "continuous")) {
    stop("'type' must be one of 'auto','discrete' and 'continuous'.")
  }
  if (type == "auto") { # Automatically determine type based on x
    if (is.numeric(x)) {
      type <- "continuous"
    } else {
      type <- "discrete"
    }
  }
  
  if (type == "discrete") {
    if (!is.factor(x)) {
      x <- factor(x, levels = unique(x)) # Convert to factor if not already
    }
    n_x <- nlevels(x) # Number of unique levels in x
    if (isTRUE(attr(palcolor, "type") == "continuous")) {
      color <- grDevices::colorRampPalette(palcolor)(n_x) # Interpolate continuous palette for discrete levels
    } else {
      # Use colors directly, cycle if not enough, or interpolate if necessary
      color <- ifelse(rep(n_x, n_x) <= pal_n,
                      palcolor[1:n_x],
                      grDevices::colorRampPalette(palcolor)(n_x)
      )
    }
    names(color) <- levels(x) # Name colors by factor levels
    if (any(is.na(x))) {
      color <- c(color, stats::setNames(NA_color, "NA")) # Add NA color if NAs present
    }
    if (isTRUE(matched)) { # Return colors matched to original x values
      color <- color[x]
      color[is.na(color)] <- NA_color
    }
  } else if (type == "continuous") {
    if (!is.numeric(x)) {
      stop("'x' must be of numeric type when using continuous color palettes.")
    }
    
    color_gradient_steps <- grDevices::colorRampPalette(palcolor)(n) # Create gradient
    
    if (isTRUE(reverse)) {
      color_gradient_steps <- rev(color_gradient_steps)
    }
    
    return(color_gradient_steps)
  }
  
  if (isTRUE(reverse)) {
    color <- rev(color)
  }
  if (!isTRUE(NA_keep)) { # Remove NA color if not keeping NAs
    color <- color[names(color) != "NA"]
  }
  return(color)
}


# Custom DimPlot-like function for cell grouping (similar to SCP's CellDimPlot)
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
    srt@meta.data[[split.by]] <- factor("") # Create a dummy column if no splitting
  }
  
  cols_to_check <- unique(c(group.by, split.by))
  if (!is.null(stat.by)) {
    cols_to_check <- unique(c(cols_to_check, stat.by))
  }
  
  # Ensure all relevant metadata columns are factors and handle NAs if requested
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
  
  # Determine the reduction to use
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
  colnames(dat_dim) <- paste0(reduction_key, dims) # Rename columns for clarity
  
  # Combine dimensionality reduction data and metadata
  dat_meta_orig <- srt@meta.data[, unique(c(group.by, split.by, stat.by)), drop = FALSE]
  
  common_cells <- intersect(rownames(dat_dim), rownames(dat_meta_orig))
  if (length(common_cells) == 0) {
    stop("No common cells found across reduction data and metadata. Check Seurat object integrity.")
  }
  
  dat_dim_ordered <- dat_dim[common_cells, , drop = FALSE]
  dat_meta_ordered <- dat_meta_orig[common_cells, , drop = FALSE]
  dat_use <- cbind(dat_dim_ordered, dat_meta_ordered)
  
  # Subset cells if specified
  if (!is.null(cells)) {
    dat_use <- dat_use[intersect(rownames(dat_use), cells), , drop = FALSE]
  }
  # Set point size if not specified
  if (is.null(pt.size)) {
    pt.size <- min(3000 / nrow(dat_use), 0.5)
  }
  
  plist <- list() # List to store plots
  xlab <- xlab %||% paste0(reduction_key, dims[1])
  ylab <- ylab %||% paste0(reduction_key, dims[2])
  
  # Create all combinations of split.by and group.by for plotting
  comb <- expand.grid(split = levels(dat_use[[split.by]]), group = group.by, stringsAsFactors = FALSE)
  rownames(comb) <- paste0(comb[["split"]], ":", comb[["group"]])
  
  # Prepare colors for stat.by if a pie chart is requested
  stat_by_levels <- if (!is.null(stat.by)) levels(dat_use[[stat.by]]) else NULL
  stat_by_colors <- if (!is.null(stat.by)) palette_scp(stat_by_levels, palette = "Set1", type = "discrete") else NULL
  
  # Generate plots for each combination
  plist <- lapply(setNames(rownames(comb), rownames(comb)), function(i) {
    legend_list <- list() # For additional legends (e.g., stat.by)
    
    g <- comb[i, "group"] # Current group.by variable
    s <- comb[i, "split"] # Current split.by variable
    
    dat_filtered <- dat_use
    if (s != "") { # Filter if split.by is active
      cells_mask <- dat_filtered[[split.by]] != s
      dat_filtered[[g]][cells_mask] <- NA # Set group.by to NA for cells not in current split
      if (!is.null(stat.by)) dat_filtered[[stat.by]][cells_mask] <- NA
    }
    
    dat_filtered <- dat_filtered %>% filter(!is.na(.data[[g]])) # Remove NA cells
    
    labels_tb <- table(dat_filtered[[g]]) # Count cells per group
    labels_tb <- labels_tb[labels_tb != 0] # Remove empty groups
    colors <- palette_scp(levels(dat_use[[g]]), palette = palette, palcolor = palcolor, NA_keep = TRUE)
    
    label_use <- names(labels_tb) # Labels for legend and plot
    
    p <- ggplot2::ggplot(dat_filtered, ggplot2::aes(x = .data[[paste0(reduction_key, dims[1])]], y = .data[[paste0(reduction_key, dims[2])]], color = .data[[g]])) +
      ggplot2::geom_point(size = pt.size, alpha = pt.alpha) +
      ggplot2::scale_color_manual(
        name = paste0(g, ":"),
        values = colors[names(labels_tb)],
        labels = label_use,
        na.value = bg_color, # Color for NA cells
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
    
    # Prepare data for labels (centroids of clusters)
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
    
    # Add labels to plot if requested
    if (isTRUE(label)) {
      if (isTRUE(show_stat)) { # Add cell counts to labels
        label_df[["count_info"]] <- vapply(as.character(label_df[["label"]]), function(grp) {
          count_val <- labels_tb[grp]
          if (is.na(count_val)) { return("") } else { return(paste0(" (n=", count_val, ")")) }
        }, FUN.VALUE = character(1), USE.NAMES = FALSE)
        label_df[["label_text"]] <- paste0(label_df[["label_text"]], label_df[["count_info"]])
      }
      
      # Use ggrepel for non-overlapping labels
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
          point.size = NA, max.overlaps = 100, force = 0, # No point, no repulsion for static labels
          color = label.fg, bg.color = label.bg, bg.r = label.bg.r, size = label.size, inherit.aes = FALSE
        )
      }
    }
    
    # Add pie charts if stat.by is specified and label is TRUE
    if (!is.null(stat.by) && isTRUE(label) && stat_plot_type == "pie") {
      pie_data_per_cluster <- dat_filtered %>%
        dplyr::group_by(.data[[g]], .data[[stat.by]]) %>%
        dplyr::summarise(count = n(), .groups = 'drop') %>%
        dplyr::group_by(.data[[g]]) %>%
        dplyr::mutate(percentage = count / sum(count)) %>%
        dplyr::ungroup()
      
      # Calculate size of pie charts based on UMAP range
      umap_range_x <- max(dat_dim[, 1], na.rm = TRUE) - min(dat_dim[, 1], na.rm = TRUE)
      umap_range_y <- max(dat_dim[, 2], na.rm = TRUE) - min(dat_dim[, 2], na.rm = TRUE)
      
      size_factor_x <- umap_range_x * stat_plot_size
      size_factor_y <- umap_range_y * stat_plot_size
      
      # Loop through each cluster to create and add pie charts
      for (cluster_label in unique(pie_data_per_cluster[[g]])) {
        current_cluster_data <- pie_data_per_cluster %>%
          dplyr::filter(.data[[g]] == cluster_label) %>%
          dplyr::arrange(.data[[stat.by]])
        
        centroid <- label_df %>%
          dplyr::filter(.data[["label"]] == cluster_label) %>%
          dplyr::select(x, y)
        
        if (nrow(centroid) == 0 || nrow(current_cluster_data) == 0) next # Skip if no centroid or data
        
        # Create a tiny ggplot for the pie chart
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
        
        # Add the pie chart as a custom annotation at the cluster centroid
        p <- p + annotation_custom(
          grob = pie_grob,
          xmin = centroid$x - size_factor_x / 2,
          xmax = centroid$x + size_factor_x / 2,
          ymin = centroid$y - size_factor_y / 2,
          ymax = centroid$y + size_factor_y / 2
        )
      }
      
      # Create a separate legend for the pie chart variable
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
    
    # Combine main plot with any extra legends
    if (length(legend_list) > 0) {
      legend_list <- legend_list[!sapply(legend_list, is.null)] # Remove NULL legends
      if (length(legend_list) > 0) {
        legend_combined <- cowplot::plot_grid(plotlist = legend_list, ncol = 1, align = "v", axis = "l")
        p_final <- patchwork::wrap_plots(p + ggplot2::theme(legend.position = "none"), legend_combined, ncol = 2, widths = c(3,1))
      }
    }
    
    return(p_final)
  })
  
  # Combine individual plots if requested
  if (isTRUE(combine)) {
    if (length(plist) > 1) {
      plot <- patchwork::wrap_plots(plotlist = plist, nrow = nrow, ncol = ncol, byrow = byrow)
    } else {
      plot <- plist[[1]]
    }
    return(plot)
  } else {
    return(plist) # Return list of plots if not combined
  }
}

# Custom FeaturePlot-like function for continuous variables (similar to SCP's FeatureDimPlot)
FeatureDimPlot <- function(srt, variable, reduction = NULL, dims = c(1, 2), split.by = NULL, cells = NULL,
                           show_na = FALSE, show_stat = FALSE, # show_stat not relevant for FeatureDimPlot
                           pt.size = NULL, pt.alpha = 1, palette = "Greys", palcolor = NULL, bg_color = "grey80",
                           label = FALSE, label.size = 4, label.fg = "white", label.bg = "black", label.bg.r = 0.1,
                           label_insitu = FALSE, label_repel = FALSE, label_repulsion = 20, # label_insitu not implemented in this version
                           label_point_size = 1, label_point_color = "black", label_segment_color = "black",
                           cells.highlight = NULL, cols.highlight = "black", sizes.highlight = 1, alpha.highlight = 1, stroke.highlight = 0.5,
                           aspect.ratio = 1, title = NULL, subtitle = NULL, xlab = NULL, ylab = NULL,
                           legend.position = "right", legend.direction = "vertical",
                           theme_use = "theme_scp", theme_args = list(), # theme_use not implemented
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
    dat[[f]][cells_mask] <- NA # Set variable to NA for cells not in current split
    
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
    return(plist) # Return list of plots if not combined
  }
}

# Helper for SCP-like plots to add arrows
add_umap_arrows_for_scp <- function(plot, srt_object, reduction_name, show_legend = TRUE) {
  # This function is specifically tailored for SCP-like plotting functions which might return
  # ggplot objects that need an arrow overlay. It uses the Seurat object to get reduction embeddings.
  
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
  
  # Add arrows and labels using annotation_custom for precise placement
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
                             clip = "off") + # Allow drawing outside plot area
    ggplot2::theme(legend.position = ifelse(show_legend, "right", "none")) # Control legend display
  
  return(plot_with_arrows)
}


# --- PIPELINE START ---
# User chooses starting point or loads data.

# --- 1. Interactive Directory Selection & Pipeline Start Decision ---

# Function to select a directory (cross-platform RStudio/base R)
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

seurat_object <- NULL # Initialize seurat_object
start_step_numeric <- 0 # Initialize for actual step number

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
      break # Exit loop to proceed to data import
    } else if (user_choice_file > 0 && user_choice_file <= length(saved_seurat_files)) {
      seurat_object <- readRDS(saved_seurat_files[user_choice_file])
      message(paste("Loaded Seurat object from:", saved_seurat_files[user_choice_file]))
      
      message("\n--- Pipeline Steps: ---")
      message("1: Data Acquisition and Initial Processing (if starting fresh)")
      message("2: Quality Control and Filtering")
      message("3: Normalization and Scaling")
      message("4: Gene Symbol Conversion")
      message("5: Dimensionality Reduction (PCA and UMAP)")
      message("6: Visualization")
      message("7: Find Marker Genes for Clusters")
      message("8: Cell Type Annotation & CAF Identification (SingleR)")
      message("9: Differential Expression Analysis")
      message("10: Malignancy Detection (scATOMIC Classifier)")
      message("11: Cell Cycle Scoring")
      message("12: Generate SCP-like Cell Cycle Plots")
      message("13: Interactive Gene Expression Analysis & Visualization")
      message("14: Trajectory inference and Pseudotime Analysis (Monocle3)")
      message("15: Cell-Cell Communication Analysis (CellChat)")
      message("16: Cell-Cell Communication Analysis (CellphoneDB)") # NEWLY INTEGRATED STEP
      message("17: Final Data Saving")
      message("-----------------------")
      
      user_choice_step <- readline(prompt = "From which step would you like to resume (1-17)? ") # Updated max step
      start_step_numeric <- as.integer(user_choice_step)
      
      if (is.na(start_step_numeric) || start_step_numeric < 1 || start_step_numeric > 17) { # Updated max step
        warning("Invalid step. Enter a number between 1 and 17. Reloading object and re-asking.")
        seurat_object <- NULL # Reset seurat_object to force re-load or new start
        next # Continue to the next iteration of the while loop
      } else {
        step_names_map <- c(
          "Data_Acquisition", # 1
          "QC", # 2
          "Normalization", # 3
          "Gene_Symbol_Conversion", # 4
          "DimRed", # 5
          "Visualization", # 6
          "Markers", # 7
          "Annotation_CAF", # 8
          "DE", # 9
          "scATOMIC_Detection", # 10
          "Cell_Cycle_Scoring", # 11
          "Generate_SCP_Plots", # 12
          "Interactive_Gene_Analysis", # 13
          "Trajectory_Inference", # 14
          "CellChat_Analysis", # 15
          "CellphoneDB_Analysis", # 16
          "Final_Save" # 17
        )
        message(paste("Starting pipeline from step:", step_names_map[start_step_numeric]))
        choice_made <- TRUE
      }
    } else {
      message("Invalid file choice. Enter a valid number or '0'." )
    }
  }
} else {
  message("\nNo existing Seurat object found. Proceeding to import raw data.")
  # start_step_numeric will remain 0, leading to initial data acquisition
}

current_seurat_object <- seurat_object # Assign loaded object to working variable

# Initialize the pipeline step marker for conditional execution
current_pipeline_step_marker <- 1


# --- 2. Data Acquisition and Initial Processing ---
# This step only runs if no existing Seurat object was loaded (start_step_numeric is 0 or 1)
if (current_pipeline_step_marker >= start_step_numeric) {
  message("--- Step ", current_pipeline_step_marker, ": Data Acquisition and Initial Processing ---")
  # Function to import 10x Genomics data
  import_10x_data <- function(data_directory) {
    message(paste("Scanning directory for 10x Genomics data files:", data_directory))
    
    h5_files <- list.files(data_directory, pattern = "\\.h5$", full.names = TRUE, recursive = TRUE)
    mtx_dirs <- c()
    
    # Check for MTX format directories
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
      # Handle potential multi-assay H5 files
      if (is.list(seurat_obj_raw)) {
        if ("Gene Expression" %in% names(seurat_obj_raw)) {
          seurat_obj_raw <- seurat_obj_raw[["Gene Expression"]]
          message("Extracted 'Gene Expression' assay from H5 file.")
        } else if (length(seurat_obj_raw) == 1) {
          seurat_obj_raw <- seurat_obj_raw[[1]]
          message(paste("Extracted the only assay (", names(seurat_obj_raw)[1], ") from H5 file.", sep=""))
        } else {
          stop("H5 file contains multiple data types, but 'Gene Expression' was not found. Please specify which assay to use.")
        }
      }
    } else if (data_type == "mtx") {
      message(paste("Importing data from MTX folder:", selected_path))
      seurat_obj_raw <- Read10X(data.dir = selected_path)
    }
    
    seurat_object_out <- CreateSeuratObject(counts = seurat_obj_raw, project = "SingleCellAnalysis")
    message("Seurat object created from raw data.")
    
    return(seurat_object_out)
  }
  
  # Only run if current_seurat_object is NULL (meaning no existing object was loaded)
  if (is.null(current_seurat_object)) {
    current_seurat_object <- import_10x_data(data_directory = data_dir)
  } else {
    message("Seurat object already loaded, skipping raw data acquisition.")
  }
}
current_pipeline_step_marker <- current_pipeline_step_marker + 1


# --- 3. Quality Control and Filtering ---
if (current_pipeline_step_marker >= start_step_numeric) {
  message("--- Step ", current_pipeline_step_marker, ": Quality Control and Filtering ---")
  perform_qc_and_filter <- function(seurat_obj, min_umis = 500, max_umis = 100000,
                                    max_percent_mt = 10, min_features = 200, min_cells = 3) {
    message("Performing Quality Control and Filtering...")
    seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
    
    # Save QC plots to PDF
    pdf(file.path(data_dir, "QC_Plots.pdf"), width = 10, height = 8)
    print(VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1))
    print(FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.mt"))
    print(FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA"))
    dev.off()
    message("QC plots saved to 'QC_Plots.pdf'. Review to refine thresholds.")
    
    # Subset cells based on QC metrics
    seurat_obj <- subset(seurat_obj,
                         subset = nFeature_RNA > min_features &
                           nFeature_RNA < max_umis &
                           percent.mt < max_percent_mt)
    
    # Filter genes that are detected in very few cells
    # Note: For versions of Seurat where filtering happens directly on assay, this might need adjustment
    # For now, manually filtering counts and then setting.
    # This line was originally present and is kept for compatibility, assuming it correctly updates the RNA assay.
    seurat_obj <- SetAssayData(seurat_obj, layer = "counts", new.data = GetAssayData(seurat_obj, assay = "RNA", layer = "counts")[rowSums(GetAssayData(seurat_obj, assay = "RNA", layer = "counts") > 0) >= min_cells, , drop = FALSE])
    
    message(paste0("QC and filtering complete. Remaining cells: ", ncol(seurat_obj)))
    return(seurat_obj)
  }
  current_seurat_object <- perform_qc_and_filter(current_seurat_object)
}
current_pipeline_step_marker <- current_pipeline_step_marker + 1


# --- 4. Normalization and Scaling ---
if (current_pipeline_step_marker >= start_step_numeric) {
  message("--- Step ", current_pipeline_step_marker, ": Normalization and Scaling ---")
  perform_normalization_and_scaling <- function(seurat_obj) {
    message("Normalizing and Scaling data...")
    seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
    seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
    seurat_obj <- ScaleData(seurat_obj, features = rownames(seurat_obj)) # Scale all features
    message("Normalization and scaling complete.")
    return(seurat_obj)
  }
  current_seurat_object <- perform_normalization_and_scaling(current_seurat_object)
}
current_pipeline_step_marker <- current_pipeline_step_marker + 1


# --- 5. Gene Symbol Conversion (NEW STEP) ---
# This step ensures gene names are in SYMBOL format, crucial for CellChat and CellphoneDB.
if (current_pipeline_step_marker >= start_step_numeric) {
  message("--- Step ", current_pipeline_step_marker, ": Performing Gene Symbol Conversion ---")
  # Ensure org.Hs.eg.db is installed (required for gene ID mapping)
  if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager")
    }
    BiocManager::install("org.Hs.eg.db", update = FALSE, ask = FALSE)
  }
  library(org.Hs.eg.db) # Load the library
  
  current_seurat_object <- convert_gene_symbols(current_seurat_object)
  message("Gene symbol conversion step completed.")
}
current_pipeline_step_marker <- current_pipeline_step_marker + 1


# --- 6. Dimensionality Reduction (PCA and UMAP) ---
if (current_pipeline_step_marker >= start_step_numeric) {
  message("--- Step ", current_pipeline_step_marker, ": Dimensionality Reduction (PCA and UMAP) ---")
  perform_dimensionality_reduction <- function(seurat_obj, dims_to_use = 1:30) { # Renamed argument to avoid conflict
    message("Performing Dimensionality Reduction (PCA and UMAP)...")
    seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))
    # Visualize PCA loadings
    pdf(file.path(data_dir, "PCA_DimLoadings.pdf"), width = 10, height = 8)
    print(VizDimLoadings(seurat_obj, dims = 1:2, reduction = "pca"))
    dev.off()
    # Plot elbow curve for PCA dimension selection
    pdf(file.path(data_dir, "PCA_ElbowPlot.pdf"), width = 8, height = 6)
    print(ElbowPlot(seurat_obj, ndims = 30))
    dev.off()
    message("PCA complete. See 'PCA_DimLoadings.pdf' and 'PCA_ElbowPlot.pdf' to select appropriate number of dimensions.")
    
    # User input for number of dimensions for UMAP/clustering
    n_dims_umap_input <- readline(prompt = paste0("Enter the number of PCA dimensions to use for UMAP and clustering (e.g., 20, based on ElbowPlot, max ", length(Reductions(seurat_obj, slot = "pca")), "): "))
    n_dims_umap <- as.integer(n_dims_umap_input)
    
    if (is.na(n_dims_umap) || n_dims_umap < 1 || n_dims_umap > length(Reductions(seurat_obj, slot = "pca"))) {
      warning("Invalid number of dimensions. Using default 1:20.")
      n_dims_umap <- 20
    }
    
    message(paste0("Running UMAP with ", n_dims_umap, " PCA dimensions..."))
    seurat_obj <- RunUMAP(seurat_obj, dims = 1:n_dims_umap)
    message(paste0("Finding neighbors with ", n_dims_umap, " PCA dimensions..."))
    seurat_obj <- FindNeighbors(seurat_obj, dims = 1:n_dims_umap)
    
    # User input for clustering resolution
    res_input <- readline(prompt = "Enter clustering resolution (e.g., 0.5): ")
    resolution <- as.numeric(res_input)
    if (is.na(resolution) || resolution <= 0) {
      warning("Invalid resolution. Using default 0.5.")
      resolution <- 0.5
    }
    message(paste0("Finding clusters with resolution ", resolution, "..."))
    seurat_obj <- FindClusters(seurat_obj, resolution = resolution)
    message("UMAP and Clustering complete.")
    return(seurat_obj)
  }
  current_seurat_object <- perform_dimensionality_reduction(current_seurat_object)
}
current_pipeline_step_marker <- current_pipeline_step_marker + 1


# --- 7. Visualization ---
if (current_pipeline_step_marker >= start_step_numeric) {
  message("--- Step ", current_pipeline_step_marker, ": Generating Visualizations ---")
  perform_visualization <- function(seurat_obj) {
    message("Generating UMAP plots...")
    
    output_filename_cluster_umap <- file.path(data_dir, "UMAP_Clusters.pdf")
    pdf(output_filename_cluster_umap, width = 10, height = 8)
    p_cluster <- DimPlot(seurat_obj, reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE) +
      ggtitle("UMAP by Seurat Clusters") +
      guides(color = guide_legend(override.aes = list(size = 3)))
    print(add_umap_arrows_and_theme_void(p_cluster, seurat_obj, show_legend = FALSE))
    dev.off()
    message(paste0("UMAP plot by clusters saved to '", output_filename_cluster_umap, "'."))
    
    if ("orig.ident" %in% colnames(seurat_obj@meta.data)) {
      output_filename_sample_umap <- file.path(data_dir, "UMAP_Samples.pdf")
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
current_pipeline_step_marker <- current_pipeline_step_marker + 1


# --- 8. Find Marker Genes for Clusters ---
if (current_pipeline_step_marker >= start_step_numeric) {
  message("--- Step ", current_pipeline_step_marker, ": Finding Marker Genes for Clusters ---")
  find_cluster_markers <- function(seurat_obj) {
    message("Finding marker genes for each cluster (this may take some time)...")
    # Using presto::wilcoxauc for faster marker detection
    message("Using presto::wilcoxauc for differential gene expression analysis...")
    all_markers <- presto::wilcoxauc(seurat_obj, group_by = 'seurat_clusters', assay = 'data')
    
    top_n_markers <- all_markers %>%
      group_by(group) %>% # 'group' column from presto output corresponds to cluster
      top_n(n = 10, wt = auc) # Using AUC to rank top markers for presto output
    
    write.csv(all_markers, file.path(data_dir, "All_Cluster_Markers.csv"), row.names = FALSE)
    write.csv(top_n_markers, file.path(data_dir, "Top10_Cluster_Markers.csv"), row.names = FALSE)
    
    message("Cluster markers identified and saved to 'All_Cluster_Markers.csv' and 'Top10_Cluster_Markers.csv'.")
    
    return(all_markers) # Return the full markers table
  }
  cluster_markers <- find_cluster_markers(current_seurat_object)
}
current_pipeline_step_marker <- current_pipeline_step_marker + 1


# --- 9. Cell Type Annotation (SingleR) & CAF Identification ---
if (current_pipeline_step_marker >= start_step_numeric) {
  message("--- Step ", current_pipeline_step_marker, ": Cell Type Annotation (SingleR) & CAF Identification ---")
  perform_cell_type_annotation <- function(seurat_obj) {
    message("Performing cell type annotation using SingleR...")
    
    # Ensure celldex is installed for reference datasets
    if (!requireNamespace("celldex", quietly = TRUE)) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
      }
      BiocManager::install("celldex")
    }
    library(celldex)
    
    # Load HumanPrimaryCellAtlasData reference
    # It's recommended to download and save this large object once to avoid repeated downloads.
    message("Loading HumanPrimaryCellAtlasData for SingleR (may download if not cached)...")
    hpca.se <- HumanPrimaryCellAtlasData()
    
    # Convert Seurat object to SingleCellExperiment format for SingleR
    sce_obj <- as.SingleCellExperiment(seurat_obj)
    
    # Run SingleR annotation
    predictions <- SingleR(test = sce_obj,
                           ref = hpca.se,
                           labels = hpca.se$label.main, # Main cell type labels
                           assay.type.test = "logcounts") # Ensure correct assay is used
    
    seurat_obj$singler_labels <- predictions$pruned.labels
    
    message("Cell type annotation complete. Adding to metadata and generating UMAP.")
    
    output_filename_singler_umap <- file.path(data_dir, "UMAP_SingleR_CellTypes.pdf")
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
      # Calculate average expression of CAF markers
      seurat_obj$CAF_score <- colMeans(x = GetAssayData(seurat_obj, assay = "RNA", layer = "data")[valid_caf_markers, , drop = FALSE])
      # Define CAFs based on top 25% expression
      seurat_obj$is_CAF <- ifelse(seurat_obj$CAF_score > quantile(seurat_obj$CAF_score, 0.75, na.rm = TRUE), "CAF", "Non-CAF")
      seurat_obj$is_CAF <- factor(seurat_obj$is_CAF, levels = c("Non-CAF", "CAF"))
      
      output_filename_caf_umap <- file.path(data_dir, "UMAP_CAF_Status.pdf")
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
current_pipeline_step_marker <- current_pipeline_step_marker + 1


# --- 10. Differential Expression Analysis ---
if (current_pipeline_step_marker >= start_step_numeric) {
  message("--- Step ", current_pipeline_step_marker, ": Differential Expression Analysis ---")
  perform_differential_expression <- function(seurat_obj) {
    message("Performing Differential Expression (DE) analysis between clusters (e.g., Cluster 0 vs. All)...")
    
    # Ensure current identity is set to seurat_clusters
    Idents(seurat_obj) <- "seurat_clusters"
    
    de_results_cluster_0 <- FindMarkers(seurat_obj, ident.1 = "0", min.pct = 0.25, logfc.threshold = 0.25)
    write.csv(de_results_cluster_0, file.path(data_dir, "DE_Cluster_0_vs_All.csv"), row.names = TRUE)
    message("DE analysis for Cluster 0 vs. all other cells saved to 'DE_Cluster_0_vs_All.csv'.")
    
    # Generate Volcano Plot
    pdf(file.path(data_dir, "VolcanoPlot_Cluster_0.pdf"), width = 10, height = 10)
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
current_pipeline_step_marker <- current_pipeline_step_marker + 1


# --- INTERMEDIATE DATA SAVING (After Differential Expression) ---
# This intermediate save is good practice, especially before computationally intensive steps.
# This is NOT a user-selectable step number, so its condition is based on whether it should have been passed.
if (current_pipeline_step_marker > start_step_numeric && start_step_numeric <= (current_pipeline_step_marker - 1)) {
  message("\n--- INTERMEDIATE DATA SAVING (after Differential Expression) ---")
  output_filename_intermediate <- file.path(data_dir, "seurat_object_after_DE.rds")
  saveRDS(current_seurat_object, file = output_filename_intermediate)
  message(paste0("Seurat object saved to '", output_filename_intermediate, "'."))
} else if (start_step_numeric > (current_pipeline_step_marker - 1)) { # If starting from later step, explicitly skip message
  message("Skipping intermediate data saving after DE (already past this step).")
}
# Note: current_pipeline_step_marker is NOT incremented here as this is an auxiliary save, not a main numbered step.


# --- 11. Malignancy Detection (scATOMIC Classifier) ---
if (current_pipeline_step_marker >= start_step_numeric) {
  message("--- Step ", current_pipeline_step_marker, ": Malignancy Detection (scATOMIC Classifier) ---")
  # This step uses the custom function defined above.
  # It will prompt for cancer type input.
  current_seurat_object <- perform_scatomic_malignancy_detection(current_seurat_object)
  
  message("scATOMIC malignancy detection step completed.")
}
current_pipeline_step_marker <- current_pipeline_step_marker + 1


# --- 12. Cell Cycle Scoring ---
if (current_pipeline_step_marker >= start_step_numeric) {
  message("--- Step ", current_pipeline_step_marker, ": Cell Cycle Scoring ---")

  s.genes <- Seurat::cc.genes$s.genes
  g2m.genes <- Seurat::cc.genes$g2m.genes

  # --- DEBUGGING LINES START ---
  message(paste0("DEBUG: Length of s.genes: ", length(s.genes)))
  if (length(s.genes) > 0) {
    message(paste0("DEBUG: First 5 s.genes: ", paste(head(s.genes, 5), collapse = ", ")))
  } else {
    message("DEBUG: s.genes is empty!")
  }

  message(paste0("DEBUG: Length of g2m.genes: ", length(g2m.genes)))
  if (length(g2m.genes) > 0) {
    message(paste0("DEBUG: First 5 g2m.genes: ", paste(head(g2m.genes, 5), collapse = ", ")))
  } else {
    message("DEBUG: g2m.genes is empty!")
  }
  # --- DEBUGGING LINES END ---

  current_seurat_object <- CellCycleScoring(current_seurat_object,
                                            s.features = s.genes,
                                            g2m.features = g2m.genes,
                                            set.ident = FALSE) # Do not set identities

  message("Cell cycle scoring complete. 'S.Score', 'G2M.Score', and 'Phase' added to metadata.")

  message("Ensuring 'Phase' metadata is a factor for consistent plotting...")
  phase_order <- c("G1", "S", "G2M")
  current_seurat_object$Phase <- factor(current_seurat_object$Phase, levels = phase_order)

  user_regress_choice <- tolower(readline(prompt = "Do you want to regress out cell cycle effects (recommended if cell cycle is a major source of variation)? (yes/no): "))

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
      current_seurat_object <- RunUMAP(current_seurat_object, dims = 1:min(30, length(Reductions(current_seurat_object, slot = "pca"))))
      current_seurat_object <- FindNeighbors(current_seurat_object, dims = 1:min(30, length(Reductions(current_seurat_object, slot = "pca"))))
      current_seurat_object <- FindClusters(current_seurat_object, resolution = 0.5)
      message("Dimensionality reduction and clustering re-run.")
    } else {
      message("Skipping re-running dimensionality reduction and clustering.")
    }

  } else {
    message("Skipping cell cycle regression.")
  }

  message("Generating UMAP plots with cell cycle phases (standard Seurat plots)...")

  output_filename_cc_umap <- file.path(data_dir, "UMAP_CellCycle_Phases.pdf")
  pdf(output_filename_cc_umap, width = 10, height = 8)
  p_cc_phase <- DimPlot(current_seurat_object, reduction = "umap", group.by = "Phase", label = FALSE, repel = TRUE) +
    ggplot2::ggtitle("UMAP by Cell Cycle Phase") +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 3)))
  print(add_umap_arrows_and_theme_void(p_cc_phase, current_seurat_object, show_legend = TRUE))
  dev.off()
  message(paste0("UMAP plot by cell cycle phases saved to '", output_filename_cc_umap, "'."))

  output_filename_s_score_umap <- file.path(data_dir, "UMAP_S_Score.pdf")
  pdf(output_filename_s_score_umap, width = 10, height = 8)
  p_s_score <- FeaturePlot(current_seurat_object, features = "S.Score", reduction = "umap") +
    ggplot2::ggtitle("UMAP by S Phase Score")
  print(add_umap_arrows_and_theme_void(p_s_score, current_seurat_object, show_legend = TRUE))
  dev.off()
  message(paste0("UMAP plot by S.Score saved to '", output_filename_s_score_umap, "'."))

  output_filename_g2m_score_umap <- file.path(data_dir, "UMAP_G2M_Score.pdf")
  pdf(output_filename_g2m_score_umap, width = 10, height = 8)
  p_g2m_score <- FeaturePlot(current_seurat_object, features = "G2M.Score", reduction = "umap") +
    ggplot2::ggtitle("UMAP by G2M Phase Score")
  print(add_umap_arrows_and_theme_void(p_g2m_score, current_seurat_object, show_legend = TRUE))
  dev.off()
  message(paste0("UMAP plot by G2M.Score saved to '", output_filename_g2m_score_umap, "'."))
  message("Cell cycle scoring step completed.")
}
current_pipeline_step_marker <- current_pipeline_step_marker + 1


# --- 13. Generate SCP-like Cell Cycle Plots (with Custom Pie Legend) ---
if (current_pipeline_step_marker >= start_step_numeric) {
  message("--- Step ", current_pipeline_step_marker, ": Generating SCP-like Cell Cycle Plots (with Custom Pie Legend) ---")
  
  # --- Pre-check for required metadata columns ---
  if (!("Phase" %in% colnames(current_seurat_object@meta.data) &&
        "S.Score" %in% colnames(current_seurat_object@meta.data) &&
        "G2M.Score" %in% colnames(current_seurat_object@meta.data))) {
    stop("Cell cycle scores (Phase, S.Score, G2M.Score) not found in Seurat object metadata. Please run step ", current_pipeline_step_marker - 1, " (Cell Cycle Scoring) first.")
  }
  
  # --- Start of interactive Cell Type (scATOMIC_prediction) selection for plotting ---
  if (!("scATOMIC_prediction" %in% colnames(current_seurat_object@meta.data)) ||
      length(unique(current_seurat_object$scATOMIC_prediction)) == 0) {
    stop("Error: 'scATOMIC_prediction' column not found or is empty in Seurat object metadata. Cannot proceed with interactive cell type selection for plotting.")
  }
  
  current_cell_types_atomic <- levels(factor(current_seurat_object$scATOMIC_prediction))
  cell_type_counts_atomic <- table(current_seurat_object$scATOMIC_prediction)
  
  cat("\nAvailable cell types (from scATOMIC_prediction) in Seurat object for plotting:\n")
  if (length(current_cell_types_atomic) == 0) {
    stop("No cell types found in 'scATOMIC_prediction' for plotting. Please check your data.")
  }
  
  for (i in seq_along(current_cell_types_atomic)) {
    cat(sprintf("%d. %s (%d cells)\n", i, current_cell_types_atomic[i], cell_type_counts_atomic[current_cell_types_atomic[i]]))
  }
  
  selected_indices_input_cc_plot <- readline(prompt = "Enter the numbers of cell types to plot (e.g., '1,3,5' or '2-4' for all cell cycle plots): ")
  
  selected_indices_cc_plot <- unique(unlist(sapply(strsplit(selected_indices_input_cc_plot, ",")[[1]], function(x) {
    if (grepl("-", x)) {
      range_vals <- as.numeric(strsplit(x, "-")[[1]])
      return(seq(min(range_vals), max(range_vals)))
    } else {
      return(as.numeric(x))
    }
  })))
  
  selected_cell_types_for_plotting <- current_cell_types_atomic[selected_indices_cc_plot]
  
  if (length(selected_cell_types_for_plotting) == 0 || any(is.na(selected_cell_types_for_plotting))) {
    stop("No valid cell types selected or invalid input provided for plotting. Exiting plotting step.")
  }
  cat("\nSelected cell types for cell cycle plots:", paste(selected_cell_types_for_plotting, collapse = ", "), "\n")
  
  seurat_object_filtered_for_cc_plots <- subset(current_seurat_object,
                                                scATOMIC_prediction %in% selected_cell_types_for_plotting)
  seurat_object_filtered_for_cc_plots$scATOMIC_prediction <- factor(seurat_object_filtered_for_cc_plots$scATOMIC_prediction,
                                                                    levels = selected_cell_types_for_plotting)
  # --- End of interactive Cell Type (scATOMIC_prediction) selection for plotting ---
  
  # --- Start of Custom UMAP Plot with Pies as Custom Legend Panel ---
  
  # I. Preparation: Generate colors for main UMAP and custom legend icons
  # Using scales::hue_pal() for perceptually distinct categorical colors
  cell_type_colors <- scales::hue_pal()(length(selected_cell_types_for_plotting))
  names(cell_type_colors) <- selected_cell_types_for_plotting # Name the vector for easy lookup
  
  # 1. Main UMAP Plot (points colored by scATOMIC_prediction, NO PIES HERE)
  p_main_umap <- DimPlot(seurat_object_filtered_for_cc_plots,
                         reduction = "umap",
                         group.by = "scATOMIC_prediction",
                         pt.size = 0.5,
                         cols = cell_type_colors) + # Use named color vector
    ggplot2::ggtitle("UMAP of Selected Cell Types") +
    ggplot2::theme_void() + # Minimal theme for UMAP background
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 16)) +
    NoLegend() # We will create a custom legend panel
  
  # Add arrows to the UMAP plot
  p_main_umap <- add_umap_arrows_for_scp(p_main_umap, seurat_object_filtered_for_cc_plots, reduction_name = "umap", show_legend = FALSE) # Arrows legend handled by main legend if needed
  
  # Add cell type labels to the UMAP plot without overlap and no connecting lines
  # Calculate centroids for label positioning
  centroids <- as.data.frame(Seurat::Embeddings(object = seurat_object_filtered_for_cc_plots, reduction = "umap")) %>%
    # Rename columns by position to ensure they are UMAP_1 and UMAP_2
    dplyr::rename(UMAP_1 = 1, UMAP_2 = 2) %>% # Renames the 1st column to UMAP_1 and 2nd to UMAP_2
    cbind(scATOMIC_prediction = seurat_object_filtered_for_cc_plots$scATOMIC_prediction) %>%
    dplyr::group_by(scATOMIC_prediction) %>%
    dplyr::summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2), .groups = 'drop') # Added .groups = 'drop' for good practice
  
  p_main_umap <- p_main_umap +
    ggrepel::geom_text_repel(data = centroids,
                             aes(x = UMAP_1, y = UMAP_2, label = scATOMIC_prediction),
                             box.padding = 0.6, # Increased padding
                             point.padding = 0.6, # Increased padding
                             segment.color = NA, # No connecting lines
                             # max.overlaps = Inf, # Removed to allow ggrepel to prevent overlaps
                             min.segment.length = 0)
  
  # 2. Calculate cell cycle phase proportions for each scATOMIC_prediction group
  cell_cycle_proportions <- seurat_object_filtered_for_cc_plots@meta.data %>%
    group_by(scATOMIC_prediction, Phase) %>%
    summarise(count = n(), .groups = 'drop') %>%
    group_by(scATOMIC_prediction) %>%
    mutate(proportion = count / sum(count)) %>%
    ungroup() %>%
    dplyr::select(scATOMIC_prediction, Phase, proportion) %>%
    pivot_wider(names_from = Phase, values_from = proportion, values_fill = 0)
  
  cell_cycle_phase_cols <- levels(seurat_object_filtered_for_cc_plots@meta.data$Phase)
  for (col in cell_cycle_phase_cols) {
    if (!(col %in% colnames(cell_cycle_proportions))) {
      cell_cycle_proportions[[col]] <- 0
    }
  }
  
  cell_cycle_proportions_long <- cell_cycle_proportions %>%
    pivot_longer(cols = all_of(cell_cycle_phase_cols), names_to = "Phase", values_to = "Proportion") %>%
    mutate(Phase = factor(Phase, levels = c("G1", "S", "G2M")))
  
  # 3. Generate individual pie charts for each selected cell type, with color icon and title
  pie_plots_list <- list()
  for (cell_type in levels(cell_cycle_proportions_long$scATOMIC_prediction)) {
    current_pie_data <- cell_cycle_proportions_long %>%
      filter(scATOMIC_prediction == cell_type)
    
    # Create the pie chart itself (without title)
    p_pie_only <- ggplot(current_pie_data, aes(x = "", y = Proportion, fill = Phase)) +
      geom_bar(stat = "identity", width = 0.8, color = "black", linewidth = 0.2) + # Reduced pie radius
      coord_polar("y", start = 0) +
      scale_fill_manual(values = c("G1" = "steelblue", "S" = "firebrick", "G2M" = "darkgreen")) +
      ggplot2::theme_void() +
      NoLegend()
    
    # Create a panel with the color swatch (now a circle) and the cell type name
    p_color_and_title_panel <- ggplot() +
      geom_point(aes(x = 0.15, y = 0.5), shape = 21, size = 6, fill = cell_type_colors[cell_type], color = "black", stroke = 0.2) + # Changed to circle icon
      geom_text(aes(x = 0.3, y = 0.5, label = cell_type), hjust = 0, size = 3.5) + # Cell type name, increased font size
      ggplot2::theme_void() +
      ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1) # Set limits for panel
    
    # Combine the color/title panel with the pie chart
    combined_pie_item <- cowplot::plot_grid(p_color_and_title_panel, p_pie_only,
                                            rel_widths = c(0.6, 0.4), # Adjusted widths for better label space
                                            ncol = 2, align = "h")
    
    pie_plots_list[[cell_type]] <- combined_pie_item
  }
  
  # 4. Create a common legend for cell cycle phases
  dummy_legend_data <- data.frame(Phase = factor(c("G1", "S", "G2M"), levels = c("G1", "S", "G2M")),
                                  value = 1)
  p_legend_only <- ggplot(dummy_legend_data, aes(x = value, y = value, fill = Phase)) +
    geom_point(size = 5, shape = 21, color = "black") + # Increased icon size, set shape to circle with black border
    scale_fill_manual(values = c("G1" = "steelblue", "S" = "firebrick", "G2M" = "darkgreen"), name = "Cell Cycle Phase") +
    ggplot2::theme(legend.position = "right",
                   panel.background = ggplot2::element_blank(), # Remove background color
                   panel.grid.major = ggplot2::element_blank(), # Remove major grid lines
                   panel.grid.minor = ggplot2::element_blank()  # Remove minor grid lines
    )
  
  cell_cycle_legend <- cowplot::get_legend(p_legend_only)
  
  # 5. Arrange all plots using cowplot::plot_grid
  ncol_pies <- 2 # Force 2 columns
  pies_grid_panel <- cowplot::plot_grid(plotlist = pie_plots_list, ncol = ncol_pies, row_spacing = unit(-0.5, "lines")) # Reduced row spacing with negative value
  
  right_panel_layout <- cowplot::plot_grid(pies_grid_panel, cell_cycle_legend, ncol = 1,
                                           rel_heights = c(1, 0.2),
                                           align = "v", axis = "l")
  
  final_composite_plot <- cowplot::plot_grid(p_main_umap, right_panel_layout, nrow = 1,
                                             rel_widths = c(0.65, 0.35),
                                             align = "h", axis = "t")
  
  output_filename_custom_umap_pies_legend <- file.path(data_dir, "UMAP_CellTypes_Pies_CustomLegend.pdf")
  pdf(output_filename_custom_umap_pies_legend, width = 18, height = 10)
  print(final_composite_plot)
  dev.off()
  message(paste0("UMAP plot with Cell Cycle Pies in custom legend saved to '", output_filename_custom_umap_pies_legend, "'."))
  
  # --- End of Custom UMAP Plot with Pies as Custom Legend Panel ---
  
  # NEW PLOT: Stacked Bar Chart of Cell Cycle Phase Distribution per Cell Type (remains unchanged)
  output_filename_cc_phase_barchart <- file.path(data_dir, "CellCycle_Phase_Distribution_Per_CellType_SCP.pdf")
  pdf(output_filename_cc_phase_barchart, width = 10, height = 8)
  
  plot_data_bar <- seurat_object_filtered_for_cc_plots@meta.data %>%
    dplyr::group_by(scATOMIC_prediction, Phase) %>%
    dplyr::summarise(Count = n(), .groups = 'drop') %>%
    dplyr::group_by(scATOMIC_prediction) %>%
    dplyr::mutate(Percentage = Count / sum(Count) * 100)
  
  phase_order <- c("G1", "S", "G2M")
  plot_data_bar$Phase <- factor(plot_data_bar$Phase, levels = phase_order)
  
  p_cc_phase_barchart <- ggplot(plot_data_bar, aes(x = scATOMIC_prediction, y = Percentage, fill = Phase)) +
    geom_bar(stat = "identity", position = "fill", color = "black", linewidth = 0.2) +
    scale_fill_manual(values = c("G1" = "steelblue", "S" = "firebrick", "G2M" = "darkgreen")) +
    labs(x = "Cell Type (scATOMIC Prediction)",
         y = "Percentage of Cells",
         title = "Cell Cycle Phase Distribution within Selected Cell Types") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5, face = "bold"),
          legend.position = "right")
  
  print(p_cc_phase_barchart)
  dev.off()
  message(paste0("Cell cycle phase distribution bar chart for selected cell types saved to '", output_filename_cc_phase_barchart, "'."))
  
  
  # UMAP by S Phase Score for Selected Cell Types (theme adjustments for no background/grid)
  output_filename_s_score_scp <- file.path(data_dir, "UMAP_S_Score_SCP.pdf")
  pdf(output_filename_s_score_scp, width = 10, height = 8)
  p_s_score_scp <- FeatureDimPlot(seurat_object_filtered_for_cc_plots,
                                  variable = "S.Score",
                                  reduction = "umap",
                                  pt.size = 0.5,
                                  palette = "Spectral") +
    ggplot2::ggtitle("UMAP by S Phase Score for Selected Cell Types (SCP)") +
    ggplot2::theme(
      panel.background = ggplot2::element_blank(), # Remove background color
      panel.grid.major = ggplot2::element_blank(), # Remove major grid lines
      panel.grid.minor = ggplot2::element_blank(),  # Remove minor grid lines
      axis.title.x = ggplot2::element_blank(), # Remove x-axis title
      axis.title.y = ggplot2::element_blank(), # Remove y-axis title
      axis.text.x = ggplot2::element_blank(),  # Remove x-axis text (values/ticks)
      axis.text.y = ggplot2::element_blank(),   # Remove y-axis text (values/ticks)
      axis.ticks = ggplot2::element_blank()     # Remove axis ticks
    )
  print(add_umap_arrows_for_scp(p_s_score_scp, seurat_object_filtered_for_cc_plots, reduction_name = "umap", show_legend = TRUE))
  dev.off()
  message(paste0("UMAP plot by S.Score (SCP style for selected cell types) saved to '", output_filename_s_score_scp, "'."))
  
  
  # UMAP by G2M Phase Score for Selected Cell Types (theme adjustments for no background/grid)
  output_filename_g2m_score_scp <- file.path(data_dir, "UMAP_G2M_Score_SCP.pdf")
  pdf(output_filename_g2m_score_scp, width = 10, height = 8)
  p_g2m_score_scp <- FeatureDimPlot(seurat_object_filtered_for_cc_plots,
                                    variable = "G2M.Score",
                                    reduction = "umap",
                                    pt.size = 0.5,
                                    palette = "Spectral") +
    ggplot2::ggtitle("UMAP by G2M Phase Score for Selected Cell Types (SCP)") +
    ggplot2::theme(
      panel.background = ggplot2::element_blank(), # Remove background color
      panel.grid.major = ggplot2::element_blank(), # Remove major grid lines
      panel.grid.minor = ggplot2::element_blank(),  # Remove minor grid lines
      axis.title.x = ggplot2::element_blank(), # Remove x-axis title
      axis.title.y = ggplot2::element_blank(), # Remove y-axis title
      axis.text.x = ggplot2::element_blank(),  # Remove x-axis text (values/ticks)
      axis.text.y = ggplot2::element_blank(),   # Remove y-axis text (values/ticks)
      axis.ticks = ggplot2::element_blank()     # Remove axis ticks
    )
  print(add_umap_arrows_for_scp(p_g2m_score_scp, seurat_object_filtered_for_cc_plots, reduction_name = "umap", show_legend = TRUE))
  dev.off()
  message(paste0("UMAP plot by G2M.Score (SCP style for selected cell types) saved to '", output_filename_g2m_score_scp, "'."))
  
  # NEW PLOT: UMAP by Normalized Proliferation Score with Cell Type Labels
  message("--- Generating UMAP plot for Normalized Proliferation Score (Revised for Legend and Control) ---")
  
  # 1. Normalize S.Score and G2M.Score using min-max scaling (0-1 range)
  # Function to perform min-max normalization
  normalize_score <- function(x) {
    min_val <- min(x, na.rm = TRUE)
    max_val <- max(x, na.rm = TRUE)
    if (max_val == min_val) { # Handle cases where all values are the same
      return(rep(0.5, length(x))) # Return 0.5 or any constant if no variance
    }
    (x - min_val) / (max_val - min_val)
  }
  
  # Apply normalization
  seurat_object_filtered_for_cc_plots$Normalized_S_Score <-
    normalize_score(seurat_object_filtered_for_cc_plots$S.Score)
  
  seurat_object_filtered_for_cc_plots$Normalized_G2M_Score <-
    normalize_score(seurat_object_filtered_for_cc_plots$G2M.Score)
  
  # 2. Combine the normalized scores
  seurat_object_filtered_for_cc_plots$Normalized_Proliferation_Score <-
    seurat_object_filtered_for_cc_plots$Normalized_S_Score +
    seurat_object_filtered_for_cc_plots$Normalized_G2M_Score
  
  # Re-normalize the combined score to 0-1 range
  seurat_object_filtered_for_cc_plots$Normalized_Proliferation_Score <-
    normalize_score(seurat_object_filtered_for_cc_plots$Normalized_Proliferation_Score)
  
  # Extract UMAP coordinates and the new score into a data frame for direct ggplot2 use
  # CORRECTED: Ensure UMAP_1 and UMAP_2 column names are set correctly
  plot_data_proliferation_umap <- as.data.frame(Seurat::Embeddings(object = seurat_object_filtered_for_cc_plots, reduction = "umap")) %>%
    dplyr::rename(UMAP_1 = 1, UMAP_2 = 2) # Explicitly rename the first two columns to UMAP_1 and UMAP_2
  
  plot_data_proliferation_umap$Normalized_Proliferation_Score <- seurat_object_filtered_for_cc_plots$Normalized_Proliferation_Score
  
  output_filename_normalized_proliferation_score_umap <- file.path(data_dir, "UMAP_Normalized_Proliferation_Score_CellTypes.pdf")
  pdf(output_filename_normalized_proliferation_score_umap, width = 12, height = 10)
  
  p_normalized_proliferation_score_umap <- ggplot2::ggplot(plot_data_proliferation_umap,
                                                           ggplot2::aes(x = UMAP_1, y = UMAP_2, color = Normalized_Proliferation_Score)) +
    ggplot2::geom_point(size = 0.5, shape = 16) + # Use geom_point directly
    
    # Use ggplot2::scale_color_viridis_c for continuous data for better perceptual uniformity and automatic legend
    ggplot2::scale_color_viridis_c(option = "plasma", name = "Proliferation Score") + # Added name for clarity, legend should now appear
    
    ggplot2::ggtitle("UMAP: Normalized Proliferation Score with Cell Type Labels") +
    ggplot2::theme_void() + # Minimal theme
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 16),
                   legend.position = "right", # Explicitly set legend position
                   plot.margin = ggplot2::margin(t = 5, r = 20, b = 15, l = 15, unit = "pt") # Adjust right margin for legend
    ) +
    ggplot2::coord_cartesian(clip = "off") # Important: allow drawing outside plot area for arrows and labels
  
  # Add arrows to the UMAP plot (ensure add_umap_arrows_for_scp is available)
  p_normalized_proliferation_score_umap <- add_umap_arrows_for_scp(p_normalized_proliferation_score_umap,
                                                                   seurat_object_filtered_for_cc_plots, # Still need Seurat object for coordinates
                                                                   reduction_name = "umap",
                                                                   show_legend = TRUE)
  
  # Add cell type labels to the UMAP plot without overlap and no connecting lines
  # Reuse the 'centroids' data frame calculated earlier (assuming 'centroids' is available in scope)
  p_normalized_proliferation_score_umap <- p_normalized_proliferation_score_umap +
    ggrepel::geom_text_repel(data = centroids,
                             aes(x = UMAP_1, y = UMAP_2, label = scATOMIC_prediction),
                             box.padding = 0.6,
                             point.padding = 0.6,
                             segment.color = NA,
                             min.segment.length = 0,
                             size = 3.5, # Adjust label font size
                             color = "black", # Label color
                             bg.color = "white", # Background color for labels
                             bg.r = 0.15) # Radius for the background box
  
  print(p_normalized_proliferation_score_umap)
  dev.off()
  message(paste0("UMAP plot for Normalized Proliferation Score with cell type labels saved to '", output_filename_normalized_proliferation_score_umap, "'."))
  
  message("All cell cycle plots generated using SCP-like functions and saved for selected scATOMIC_prediction cell types.")
}
current_pipeline_step_marker <- current_pipeline_step_marker + 1

# --- 14. Interactive Gene Expression Analysis & Visualization ---
# This step allows interactive plotting of gene expression on UMAP,
# with gene analysis performed on all cells, and plotting filtered by user-selected cell types.
# It assumes 'current_seurat_object' is available from previous steps.
# It also assumes 'add_umap_arrows_and_theme_void' function is defined.

if (current_pipeline_step_marker >= start_step_numeric) {
  message("--- Step ", current_pipeline_step_marker, ": Interactive Gene Expression Analysis & Visualization ---")
  
  interactive_gene_analysis <- function(seurat_obj) {
    message("Starting interactive gene expression analysis. Plots will appear one by one. Close plot window to continue.")
    
    # Get available cell types once at the start, to list for user
    if (!is.null(seurat_obj$scATOMIC_prediction)) {
      all_group_counts <- table(seurat_obj$scATOMIC_prediction)
    } else {
      warning("scATOMIC_prediction labels not found in Seurat object. Cell type selection and labeling for plotting will be skipped.")
      all_group_counts <- NULL # Indicate that cell type selection is not possible
    }
    
    while (TRUE) {
      gene_input <- readline(prompt = "Enter gene symbol(s) (comma-separated, e.g., EPCAM, PTPRC) or 'done' to finish: ")
      if (tolower(gene_input) == "done") {
        message("Exiting interactive gene expression analysis.")
        break
      }
      
      target_genes <- stringr::str_trim(unlist(stringr::str_split(gene_input, ",")))
      # Check if gene exists in the ENTIRE Seurat object first
      target_genes <- intersect(target_genes, rownames(seurat_obj))
      
      if (length(target_genes) == 0) {
        message("No valid genes found in the Seurat object for the given input. Please try again.")
        next
      }
      
      message(paste0("Analyzing gene(s): ", paste(target_genes, collapse = ", ")))
      
      # --- Cell Type Selection for PLOTTING (inside gene loop) ---
      seurat_object_for_plotting <- seurat_obj # Default to full object for plotting
      
      if (!is.null(all_group_counts)) {
        cat("\nAvailable scATOMIC_prediction groups for plotting this gene:\n")
        for (i in seq_along(all_group_counts)) {
          cat(sprintf("%d. %s (%d cells)\n", i, names(all_group_counts)[i], all_group_counts[i]))
        }
        
        selected_indices_input_celltypes <- readline(prompt = "Enter the numbers of cell types to plot (e.g., '1,3,5' or '2-4'), or leave blank/type 'all' to plot all cells: ")
        
        if (tolower(stringr::str_trim(selected_indices_input_celltypes)) %in% c("", "all")) {
          message("Plotting all cell types for this gene.")
          seurat_object_for_plotting <- seurat_obj
        } else {
          selected_indices_celltypes <- unique(unlist(sapply(stringr::str_split(selected_indices_input_celltypes, ",")[[1]], function(x) {
            if (grepl("-", x)) {
              range_vals <- as.numeric(stringr::str_split(x, "-")[[1]])
              return(seq(min(range_vals), max(range_vals)))
            } else {
              return(as.numeric(x))
            }
          })))
          
          # Validate input indices
          if (any(selected_indices_celltypes < 1) || any(selected_indices_celltypes > length(all_group_counts)) || length(selected_indices_celltypes) == 0) {
            message("Invalid cell type selection. Plotting all cell types for this gene.")
            seurat_object_for_plotting <- seurat_obj # Fallback to all cells if selection is invalid
          } else {
            selected_cell_types <- names(all_group_counts)[selected_indices_celltypes]
            message(paste0("Plotting gene expression for cell types: ", paste(selected_cell_types, collapse = ", "), "\n"))
            seurat_object_for_plotting <- subset(seurat_obj, subset = scATOMIC_prediction %in% selected_cell_types)
            
            # Ensure that the factor levels are dropped for the filtered object to avoid displaying unselected levels
            seurat_object_for_plotting$scATOMIC_prediction <- factor(seurat_object_for_plotting$scATOMIC_prediction)
            seurat_object_for_plotting$scATOMIC_prediction <- droplevels(seurat_object_for_plotting$scATOMIC_prediction)
            
            if (ncol(seurat_object_for_plotting) == 0) {
              warning("No cells found for the selected cell types. Plotting all cells from original object for this gene.")
              seurat_object_for_plotting <- seurat_obj
            }
          }
        }
      }
      
      save_plots_choice <- tolower(readline(prompt = "Save plots to PDF (yes/no)? ")) == "yes"
      
      for (target_gene in target_genes) {
        # --- Prepare data for cell type labels ---
        # Calculate median UMAP coordinates for each scATOMIC_prediction group
        cluster_labels_df <- NULL
        if (!is.null(all_group_counts) && ncol(seurat_object_for_plotting) > 0) {
          cluster_label_data <- as.data.frame(Seurat::Embeddings(seurat_object_for_plotting, reduction = "umap"))
          # Ensure the column name matches UMAP_1, UMAP_2 if they are different
          colnames(cluster_label_data) <- c("UMAP_1", "UMAP_2")
          cluster_label_data$cluster <- seurat_object_for_plotting$scATOMIC_prediction
          
          cluster_labels_df <- cluster_label_data %>%
            dplyr::group_by(cluster) %>%
            dplyr::summarize(
              UMAP_1 = median(UMAP_1, na.rm = TRUE),
              UMAP_2 = median(UMAP_2, na.rm = TRUE)
            )
        }
        
        # Feature plot (continuous expression)
        plot_title_feature <- paste0("UMAP of Continuous ", target_gene, " Expression")
        if (save_plots_choice) {
          output_filename_feature <- file.path(data_dir, paste0("UMAP_", target_gene, "_FeaturePlot.pdf"))
          pdf(output_filename_feature, width = 8, height = 7)
        }
        
        p_feature <- Seurat::FeaturePlot(seurat_object_for_plotting, # Use the filtered object for plotting
                                         features = target_gene,
                                         reduction = "umap",
                                         pt.size = 0.5,
                                         order = TRUE) +
          ggplot2::ggtitle(plot_title_feature)
        
        # Add scATOMIC_prediction labels to FeaturePlot if available
        if (!is.null(cluster_labels_df) && nrow(cluster_labels_df) > 0) {
          p_feature <- p_feature +
            ggrepel::geom_text_repel(data = cluster_labels_df,
                                     ggplot2::aes(x = UMAP_1, y = UMAP_2, label = cluster),
                                     size = 3, # Adjust as needed
                                     color = "black",
                                     fontface = "plain",
                                     max.overlaps = Inf,
                                     segment.colour = "grey50", # Adjust as needed
                                     segment.size = 0.5, # Adjust as needed
                                     force_pull = 1) # <--- Add/Adjust this line (default is 0.5)
        }
        
        print(add_umap_arrows_and_theme_void(p_feature, seurat_object_for_plotting, show_legend = TRUE))
        
        if (save_plots_choice) {
          dev.off()
          message(paste0("UMAP feature plot for '", target_gene, "' continuous expression saved to '", output_filename_feature, "'."))
        }
        
        # Binary expression plot (positive/negative based on threshold)
        threshold_value_input <- readline(prompt = paste0("Enter expression threshold for ", target_gene, " (e.g., 0.5) or leave blank to skip binary plot: "))
        if (nchar(threshold_value_input) > 0) {
          threshold_value <- as.numeric(threshold_value_input)
          if (is.na(threshold_value)) {
            message("Invalid threshold value. Skipping binary plot for this gene.")
            next
          }
          
          # Create a binary metadata column for the currently selected plotting subset
          gene_status_vector <- as.factor(
            ifelse(Seurat::GetAssayData(seurat_object_for_plotting, assay = "RNA", layer = "data")[target_gene, ] > threshold_value,
                   paste0(target_gene, "_Positive"), paste0(target_gene, "_Negative"))
          )
          names(gene_status_vector) <- colnames(seurat_object_for_plotting)
          
          # Temporarily add metadata to the object used for plotting
          temp_seurat_obj_with_metadata <- Seurat::AddMetaData(
            object = seurat_object_for_plotting,
            metadata = gene_status_vector,
            col.name = "gene_expression_status" # Temporary column
          )
          
          plot_title_binary <- paste0("UMAP of ", target_gene, " Expression (Threshold > ", threshold_value, ")")
          if (save_plots_choice) {
            output_filename_binary <- file.path(data_dir, paste0("UMAP_", target_gene, "_Positive_Negative_Plot.pdf"))
            pdf(output_filename_binary, width = 8, height = 7)
          }
          
          p_binary <- Seurat::DimPlot(temp_seurat_obj_with_metadata, # Use object with temp metadata
                                      reduction = "umap",
                                      group.by = "gene_expression_status",
                                      cols = c("lightgrey", "red"), # Customize colors for binary groups
                                      label = FALSE, # Set to FALSE here because we are adding custom labels
                                      repel = TRUE) +
            ggplot2::ggtitle(plot_title_binary)
          
          # Add scATOMIC_prediction labels to Binary Plot if available
          if (!is.null(cluster_labels_df) && nrow(cluster_labels_df) > 0) {
            p_binary <- p_binary +
              ggrepel::geom_text_repel(data = cluster_labels_df,
                                       ggplot2::aes(x = UMAP_1, y = UMAP_2, label = cluster),
                                       size = 3, # Adjust as needed
                                       color = "black",
                                       fontface = "plain",
                                       max.overlaps = Inf,
                                       segment.colour = "grey50", # Adjust as needed
                                       segment.size = 0.5, # Adjust as needed
                                       force_pull = 1) # <--- Add/Adjust this line (default is 0.5)
          }
          
          print(add_umap_arrows_and_theme_void(p_binary, temp_seurat_obj_with_metadata, show_legend = TRUE))
          
          if (save_plots_choice) {
            dev.off()
            message(paste0("UMAP plot for '", target_gene, "' expression (binary) saved to '", output_filename_binary, "'."))
          }
          # No need to remove metadata from seurat_object_for_plotting as it was added to a temporary object
        }
      }
    }
    return(seurat_obj) # Return original object (no permanent changes to metadata)
  }
  current_seurat_object <- interactive_gene_analysis(current_seurat_object)
}
current_pipeline_step_marker <- current_pipeline_step_marker + 1

# --- 15. Trajectory Inference (Monocle3) ---
if (current_pipeline_step_marker >= start_step_numeric) {
  
  message("--- Step ", current_pipeline_step_marker, ": Starting Trajectory Inference using Monocle3 ---")
  
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
  cds_filename <- file.path(data_dir, "monocle3_cds_object.rds")
  
  # Try to load existing CDS
  # For this specific module, if we are filtering, we will *always* regenerate CDS
  # unless we implement more complex logic for saving/loading filtered CDS objects.
  # For simplicity, we'll assume a new CDS is generated if cell type filtering is applied.
  cds <- NULL # Reset CDS to force recreation if filtering is to be done
  
  if (is.null(cds)) {
    # CDS not found or forcing recreation — run full Monocle3 analysis
    message("Monocle3 CDS object not found or regeneration forced. Running full Monocle3 analysis.")
    
    # Ensure monocle3 and its dependencies are installed
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
    
    # Convert to character to handle replacement, then back to factor
    current_seurat_object@meta.data$scATOMIC_prediction <- as.character(current_seurat_object@meta.data$scATOMIC_prediction)
    
    # Define problematic strings and replacement categories
    problematic_string_scatomic <- "source(\"~/.active−rstudio−document\", echo = TRUE)" # From prior issues
    unknown_category <- "Unknown_scATOMIC_Cell"
    unassigned_category <- "Unassigned_scATOMIC"
    
    # Replace problematic string and NA values
    current_seurat_object@meta.data$scATOMIC_prediction[current_seurat_object@meta.data$scATOMIC_prediction == problematic_string_scatomic] <- unknown_category
    current_seurat_object@meta.data$scATOMIC_prediction[is.na(current_seurat_object@meta.data$scATOMIC_prediction)] <- unassigned_category
    
    # Convert back to factor with all expected levels
    all_levels <- sort(unique(c(current_seurat_object@meta.data$scATOMIC_prediction, unknown_category, unassigned_category)))
    current_seurat_object$scATOMIC_prediction <- factor(current_seurat_object@meta.data$scATOMIC_prediction,
                                                        levels = all_levels)
    message("'scATOMIC_prediction' column processed and converted to factor for Monocle3 compatibility.")
    
    # --- NEW STEP: Cell Type Selection for Monocle3 Analysis ---
    message("\n--- Cell Type Selection for Monocle3 Trajectory Analysis ---")
    unique_scatomic_labels_all <- levels(factor(current_seurat_object$scATOMIC_prediction))
    message("Available scATOMIC Cell Types in Seurat object:")
    
    # Calculate cell counts for each cell type
    cell_type_counts <- as.data.frame(table(current_seurat_object$scATOMIC_prediction))
    colnames(cell_type_counts) <- c("CellType", "Count")
    
    # Create a numbered list for display and mapping, including counts
    # Ensure the order of counts matches unique_scatomic_labels_all
    ordered_counts <- cell_type_counts$Count[match(unique_scatomic_labels_all, cell_type_counts$CellType)]
    
    numbered_cell_types <- data.frame(
      Number = 1:length(unique_scatomic_labels_all),
      CellType = unique_scatomic_labels_all,
      Count = ordered_counts # Add the count here
    )
    print(numbered_cell_types) # Print the numbered list with counts
    
    message("-----------------------------------------------------")
    message("\n\n") # Added empty lines for better visibility of the prompt
    
    # Update prompt to ask for numbers, including ranges
    selected_cell_types_input <- readline(prompt = "Enter comma-separated NUMBERS or RANGES (e.g., '1,3,4,6-10') corresponding to cell types to include, or type 'all' to include all cell types: ")
    selected_cell_types_input <- trimws(selected_cell_types_input)
    
    seurat_object_for_monocle <- NULL # Initialize
    
    if (tolower(selected_cell_types_input) == "all" || selected_cell_types_input == "") {
      message("Including all cell types in Monocle3 analysis.")
      seurat_object_for_monocle <- current_seurat_object
    } else {
      # Parse input numbers and ranges
      selected_tokens <- unlist(strsplit(selected_cell_types_input, ",", fixed = TRUE))
      selected_tokens <- trimws(selected_tokens)
      
      parsed_numbers <- c()
      for (token in selected_tokens) {
        if (grepl("-", token)) { # Check for range
          range_parts <- unlist(strsplit(token, "-", fixed = TRUE))
          if (length(range_parts) == 2) {
            start_num <- suppressWarnings(as.numeric(range_parts[1]))
            end_num <- suppressWarnings(as.numeric(range_parts[2]))
            if (!is.na(start_num) && !is.na(end_num) && start_num <= end_num) {
              parsed_numbers <- c(parsed_numbers, start_num:end_num)
            } else {
              warning(paste0("Invalid range format or non-numeric values in '", token, "'. Ignoring this part of the input."))
            }
          } else {
            warning(paste0("Invalid range format '", token, "'. Ignoring this part of the input."))
          }
        } else { # Single number
          single_num <- suppressWarnings(as.numeric(token))
          if (!is.na(single_num)) {
            parsed_numbers <- c(parsed_numbers, single_num)
          } else {
            warning(paste0("Non-numeric input '", token, "'. Ignoring this part of the input."))
          }
        }
      }
      
      # Remove duplicates and sort
      selected_numbers <- unique(parsed_numbers)
      selected_numbers <- sort(selected_numbers) # Optional: keep order for consistent messaging
      
      # Filter out out-of-range numbers
      valid_numbers <- selected_numbers[selected_numbers >= 1 & selected_numbers <= length(unique_scatomic_labels_all)]
      
      if (length(valid_numbers) == 0) {
        warning("No valid numbers or ranges entered or numbers/ranges are out of bounds. Including all cell types by default.")
        seurat_object_for_monocle <- current_seurat_object
      } else {
        # Map valid numbers back to cell type names
        selected_cell_types_vector <- unique_scatomic_labels_all[valid_numbers]
        
        message(paste0("Including the following cell types in Monocle3 analysis: ", paste(selected_cell_types_vector, collapse = ", ")))
        # Subset the Seurat object based on selected cell types
        seurat_object_for_monocle <- subset(current_seurat_object, subset = scATOMIC_prediction %in% selected_cell_types_vector)
        
        if (ncol(seurat_object_for_monocle) == 0) {
          stop("After filtering, no cells remain in the Seurat object for Monocle3 analysis. Please check your cell type selection.")
        }
      }
    }
    message("Cell type selection complete. Proceeding with Monocle3 on selected cells.")
    
    message("Converting Seurat object to Monocle3 cell_data_set (CDS)...")
    
    # Extract data for CDS object creation from the *filtered* Seurat object
    expression_matrix <- GetAssayData(seurat_object_for_monocle, assay = "RNA", layer = "counts")
    cell_metadata <- seurat_object_for_monocle@meta.data
    gene_metadata <- data.frame(gene_short_name = rownames(expression_matrix), row.names = rownames(expression_matrix))
    
    # Ensure dimensions match before creating CDS
    if (!all(rownames(cell_metadata) == colnames(expression_matrix))) {
      warning("Cell metadata rownames do not match expression matrix column names. Attempting reorder.")
      cell_metadata <- cell_metadata[colnames(expression_matrix), , drop = FALSE]
    }
    if (!all(rownames(gene_metadata) == rownames(expression_matrix))) {
      warning("Gene metadata rownames do not match expression matrix row names. Attempting reorder.")
      gene_metadata <- gene_metadata[rownames(expression_matrix), , drop = FALSE]
    }
    
    cds <- new_cell_data_set(expression_matrix,
                             cell_metadata = cell_metadata,
                             gene_metadata = gene_metadata)
    
    message("Adding Seurat UMAP embeddings to CDS...")
    # Transfer UMAP coordinates from *filtered* Seurat object to Monocle3 CDS
    umap_coords <- Embeddings(seurat_object_for_monocle, reduction = "umap")
    colnames(umap_coords) <- c("UMAP_1", "UMAP_2") # Monocle3 expects these column names
    
    reducedDims(cds)$UMAP <- umap_coords
    
    message("CDS object created and Seurat UMAP embeddings transferred.")
    
    message("Running Monocle3 clustering (using transferred UMAP)...")
    # Clustering using the transferred UMAP coordinates
    cds <- cluster_cells(cds, reduction_method = "UMAP", k = 20, resolution = 1e-6) # k and resolution can be adjusted
    message("Monocle3 clustering complete.")
    
    message("Learning cell trajectory graph with Monocle3...")
    cds <- learn_graph(cds, use_partition = FALSE) # Set use_partition = FALSE for a single trajectory
    message("Monocle3 graph learned.")
    
  } else {
    message("Reusing existing Monocle3 CDS object.")
    library(monocle3) # Ensure monocle3 is loaded if CDS was loaded
    # If a CDS was loaded, it would typically be the unfiltered one.
    # We would need additional logic here to re-filter if desired,
    # or ensure filtering happens before initial CDS saving.
    # For now, if loaded, we proceed with the loaded (potentially unfiltered) CDS.
  }
  
  # --- Interactive Root Selection for Monocle3 ---
  # Note: This root selection will now be based on the *filtered* CDS if cells were selected.
  unique_scatomic_labels_cds <- levels(factor(colData(cds)$scATOMIC_prediction))
  message("\n--- Available scATOMIC Cell Types for Monocle3 Root Selection (from current CDS) ---")
  
  # Calculate cell counts for each cell type in CDS
  cell_type_counts_cds <- as.data.frame(table(colData(cds)$scATOMIC_prediction))
  colnames(cell_type_counts_cds) <- c("CellType", "Count")
  
  # Create a numbered list for display and mapping, including counts
  ordered_counts_cds <- cell_type_counts_cds$Count[match(unique_scatomic_labels_cds, cell_type_counts_cds$CellType)]
  
  numbered_cell_types_cds <- data.frame( # Corrected from data.data.frame
    Number = 1:length(unique_scatomic_labels_cds),
    CellType = unique_scatomic_labels_cds,
    Count = ordered_counts_cds
  )
  print(numbered_cell_types_cds)
  message("-----------------------------------------------------")
  message("\n\n") # Added empty lines for better visibility of the prompt
  
  # Update prompt to ask for a number for root selection
  selected_root_input <- readline(prompt = "Enter the NUMBER corresponding to the scATOMIC cell type to use as the trajectory root, or type 'auto' for automatic detection: ")
  
  selected_root_input <- trimws(selected_root_input) # Remove leading/trailing whitespace
  
  selected_root_label <- NULL # Initialize the actual label to be used
  
  if (tolower(selected_root_input) == "auto") {
    message("Automatically detecting root node(s) for trajectory inference.")
    cds <- order_cells(cds) # Monocle3 will attempt to find a root automatically
    root_label_for_filename <- "AutoRoot"
  } else {
    # Parse the number input for root selection
    selected_root_number <- suppressWarnings(as.numeric(selected_root_input))
    
    if (!is.na(selected_root_number) && selected_root_number >= 1 && selected_root_number <= length(unique_scatomic_labels_cds)) {
      selected_root_label <- unique_scatomic_labels_cds[selected_root_number]
      message(paste0("Attempting to set '", selected_root_label, "' (Number ", selected_root_number, ") as the root node for trajectory inference."))
      
      # Identify cells belonging to the selected root type
      root_cells <- colData(cds) %>%
        as.data.frame() %>%
        dplyr::filter(scATOMIC_prediction == selected_root_label) %>%
        rownames()
      
      if (length(root_cells) == 0) {
        warning(paste0("No cells found with label '", selected_root_label, "' in the current CDS object. Reverting to automatic root detection."))
        cds <- order_cells(cds, root_cells = NULL) # Revert to auto if label not found
        root_label_for_filename <- "AutoRoot"
      } else {
        cds <- order_cells(cds, root_cells = root_cells)
        # Sanitize label for filename (remove special characters)
        root_label_for_filename <- gsub("[^A-Za-z0-9_]", "", selected_root_label)
      }
    } else {
      warning("Invalid number entered for root selection. Reverting to automatic root detection.")
      cds <- order_cells(cds) # Default to auto if invalid number
      root_label_for_filename <- "AutoRoot"
    }
  }
  
  message("Cells ordered along trajectory.")
  
  # Save cds with root label in filename
  # Append selected cell types to filename for clarity if filtering occurred
  if (!tolower(selected_cell_types_input) == "all" && selected_cell_types_input != "") {
    selected_types_short <- paste(head(selected_cell_types_vector, 2), collapse="_")
    if(length(selected_cell_types_vector) > 2) {
      selected_types_short <- paste0(selected_types_short, "_etc")
    }
    cds_filename_root <- file.path(data_dir, paste0("monocle3_cds_object_", root_label_for_filename, "Root_Filtered_", selected_types_short, ".rds"))
  } else {
    cds_filename_root <- file.path(data_dir, paste0("monocle3_cds_object_", root_label_for_filename, "Root.rds"))
  }
  
  saveRDS(cds, cds_filename_root)
  message(paste0("Monocle3 CDS object saved as '", cds_filename_root, "'."))
  
  # Plot and save trajectory plots to PDF (with print() to ensure plot rendering)
  if (!tolower(selected_cell_types_input) == "all" && selected_cell_types_input != "") {
    plot_filename_suffix <- paste0(gsub("_", " ", root_label_for_filename), "Root_Filtered_", selected_types_short)
  } else {
    plot_filename_suffix <- paste0(gsub("_", " ", root_label_for_filename), "Root")
  }
  output_filename_monocle_traj <- file.path(data_dir, paste0("UMAP_Monocle3_Trajectory_", plot_filename_suffix, ".pdf"))
  pdf(output_filename_monocle_traj, width = 10, height = 8)
  
  # Plot cells colored by pseudotime
  print(
    plot_cells(cds,
               color_cells_by = "pseudotime",
               label_groups_by_cluster = FALSE,
               label_leaves = TRUE,
               label_branch_points = TRUE,
               graph_label_size = 3) +
      ggplot2::ggtitle(paste0("Monocle3 Trajectory (Pseudotime - ", gsub("_", " ", plot_filename_suffix), ")")) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))
  )
  
  # Plot cells colored by scATOMIC prediction
  print(
    plot_cells(cds,
               color_cells_by = "scATOMIC_prediction",
               label_groups_by_cluster = TRUE,
               label_leaves = TRUE,
               label_branch_points = TRUE,
               graph_label_size = 3) +
      ggplot2::ggtitle(paste0("Monocle3 Trajectory (scATOMIC Prediction - ", gsub("_", " ", plot_filename_suffix), ")")) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))
  )
  
  dev.off()
  message(paste0("Monocle3 trajectory plots saved to '", output_filename_monocle_traj, "'."))
  
  message("Trajectory Inference step completed.")
  
}
current_pipeline_step_marker <- current_pipeline_step_marker + 1


# --- 16. Cell-Cell Communication Analysis (CellChat) ---
if (current_pipeline_step_marker >= start_step_numeric) {
  
  message("--- Step ", current_pipeline_step_marker, ": Running CellChat analysis ---")
  
  # Check if CellChat object already exists locally
  cellchat_file <- file.path(data_dir, "cellchat_object.Rdata")
  if (file.exists(cellchat_file)) {
    reuse <- tolower(readline(prompt = "CellChat object found. Do you want to reuse it? (yes/no): "))
    if (reuse == "yes") {
      load(cellchat_file)
      message("Loaded existing CellChat object.")
    } else {
      rm(list = c("cellchat"))  # Remove any previously loaded object before re-running
      if (exists("cellchat")) message("Failed to remove old cellchat object, proceeding anyway.")
    }
  }
  
  # If cellchat doesn't exist or was cleared, run new analysis
  if (!exists("cellchat")) {
    message("Running new CellChat analysis.")
    
    # Step 1: Prepare protein-coding gene list (ensure org.Hs.eg.db is loaded)
    # This assumes org.Hs.eg.db is available and gene symbols are in rownames.
    protein_coding_genes <- keys(org.Hs.eg.db, keytype = "SYMBOL")
    
    # Step 2: Extract normalized data for CellChat
    # CellChat expects gene symbols in rownames and cells in colnames.
    # The 'data' slot of the RNA assay contains log-normalized expression.
    data.input <- GetAssayData(current_seurat_object, assay = "RNA", slot = "data")
    
    # Filter for protein-coding genes (if desired, based on org.Hs.eg.db)
    data.input <- data.input[rownames(data.input) %in% protein_coding_genes, , drop = FALSE]
    
    # Step 3: Extract metadata and set 'scATOMIC_prediction' as the grouping variable
    meta <- current_seurat_object@meta.data
    
    # Ensure 'scATOMIC_prediction' is a factor and used as identity
    if (!("scATOMIC_prediction" %in% colnames(meta))) {
      stop("Metadata column 'scATOMIC_prediction' not found for CellChat grouping. Please check step 11.")
    }
    
    # Convert 'scATOMIC_prediction' to character, handle potential empty/NA values for CellChat robustness
    meta$cell_type_for_cellchat <- as.character(meta$scATOMIC_prediction)
    meta$cell_type_for_cellchat[is.na(meta$cell_type_for_cellchat) | meta$cell_type_for_cellchat == ""] <- "Unassigned"
    meta$cell_type_for_cellchat <- factor(meta$cell_type_for_cellchat)
    
    # Create CellChat object
    cellchat <- createCellChat(object = data.input, meta = meta, group.by = "cell_type_for_cellchat")
    
    # Set CellChat database (human)
    CellChatDB <- CellChatDB.human # or CellChatDB.mouse if applicable
    cellchat@DB <- CellChatDB
    
    # Step 4: Preprocess and compute communication
    cellchat <- subsetData(cellchat) # Filters out low expressed genes/cells
    
    # Parallel computing (adjust workers based on your CPU cores)
    future::plan("multisession", workers = 4) # Adjust 'workers' as needed
    
    cellchat <- identifyOverExpressedGenes(cellchat)
    cellchat <- identifyOverExpressedInteractions(cellchat)
    cellchat <- projectData(cellchat, CellChatDB)
    
    cellchat <- computeCommunProb(cellchat)
    cellchat <- filterCommunication(cellchat, min.cells = 10) # Filter out infrequent interactions
    cellchat <- computeCommunProbPathway(cellchat)
    cellchat <- aggregateNet(cellchat)
    cellchat <- netAnalysis_computeCentrality(cellchat) # Compute centrality measures
    
    # Save CellChat object for reuse
    save(cellchat, file = cellchat_file)
    message("New CellChat object created and saved as 'cellchat_object.Rdata'.")
  }
  
  # --- CellChat Visualization ---
  message("Generating CellChat visualizations...")
  
  # --- Start of interactive CellChat plotting selection block (with cell counts) ---
  
  # Get cell types (idents) from the existing CellChat object
  cellchat_cell_types <- levels(cellchat@idents)
  # Get cell counts for each cell type in the CellChat object
  cell_type_counts <- table(cellchat@idents)
  
  cat("\nAvailable cell types in CellChat object for plotting:\n")
  for (i in seq_along(cellchat_cell_types)) {
    # Display cell type name and its count
    cat(sprintf("%d. %s (%d cells)\n", i, cellchat_cell_types[i], cell_type_counts[cellchat_cell_types[i]]))
  }
  
  # Interactively get user input for cell types to include in plots
  selected_indices_input_cellchat_plot <- readline(prompt = "Enter the numbers of cell types to plot (e.g., '1,3,5' or '2-4'): ")
  
  # Parse the input
  selected_indices_cellchat_plot <- unique(unlist(sapply(strsplit(selected_indices_input_cellchat_plot, ",")[[1]], function(x) {
    if (grepl("-", x)) {
      range_vals <- as.numeric(strsplit(x, "-")[[1]])
      return(seq(min(range_vals), max(range_vals)))
    } else {
      return(as.numeric(x))
    }
  })))
  
  # Validate input indices
  if (any(selected_indices_cellchat_plot < 1) || any(selected_indices_cellchat_plot > length(cellchat_cell_types))) {
    stop("Invalid selection. Please enter numbers corresponding to the listed cell types.")
  }
  
  # Get the names of the selected cell types for plotting
  selected_cell_types_for_plotting <- cellchat_cell_types[selected_indices_cellchat_plot]
  cat("\nSelected cell types for CellChat plots:", paste(selected_cell_types_for_plotting, collapse = ", "), "\n")
  
  # --- NEW: Subset the CellChat object for plotting ---
  # This creates a new CellChat object containing only the selected cell types
  cellchat_filtered_for_plotting <- subsetCellChat(cellchat, idents.use = selected_cell_types_for_plotting)
  
  # Recompute aggregated networks for the subsetted object if needed by plots
  # (Often subsetCellChat handles this, but sometimes explicit re-aggregation is safer)
  # cellchat_filtered_for_plotting <- aggregateNet(cellchat_filtered_for_plotting)
  # cellchat_filtered_for_plotting <- netAnalysis_computeCentrality(cellchat_filtered_for_plotting)
  
  # --- End of interactive CellChat plotting selection block (with cell counts) ---
  
  
  # NetVisual Circle Plot (overall communication)
  pdf(file.path(data_dir, "cellchat_netVisual_circle_selected_types.pdf"), width = 8, height = 8)
  # Use the filtered CellChat object's net$weight
  print(netVisual_circle(cellchat_filtered_for_plotting@net$weight, # <-- MODIFICATION HERE (using filtered object's weight)
                         vertex.weight = as.numeric(table(cellchat_filtered_for_plotting@idents)), # <-- MODIFICATION HERE (using filtered object's idents)
                         weight.scale = TRUE, label.edge = FALSE,
                         title.name = "Overall Communication Strength (Selected Cell Types)"))
  dev.off()
  message("CellChat overall communication circle plot saved for selected cell types.")
  
  # Signaling Role Heatmap
  # Use the filtered CellChat object for netAnalysis_signalingRole_heatmap
  ht <- netAnalysis_signalingRole_heatmap(cellchat_filtered_for_plotting)
  ht@row_names_param$gp <- grid::gpar(fontsize = 4) # Keep or adjust fontsize (e.g., 2.5, 3.5)
  ht@row_names_param$max_width <- grid::unit(8, "cm") # Keep or adjust for longer names
  
  pdf(file.path(data_dir, "cellchat_netAnalysis_signalingRole_heatmap_selected_types.pdf"), width = 8, height = 80) # Keep or adjust PDF height as needed
  ComplexHeatmap::draw(ht,
                       heatmap_legend_side = "right",
                       heatmap_height = grid::unit(200, "cm")) # <-- MODIFIED: Using absolute unit (cm)
  dev.off()
  message("CellChat signaling role heatmap saved for selected cell types.")
  
  # Bubble Plot for Selected Pathway (Interactive)
  message("Available pathways for CellChat visualization:")
  # Get pathways from the filtered object, as some might be lost if no cells in selected group support them
  print(cellchat_filtered_for_plotting@netP$pathways) # <-- MODIFICATION HERE (using filtered object)
  selected_pathway <- readline(prompt = "Enter the name of a specific CellChat pathway to visualize (e.g., 'VEGF', 'TGFb') or press Enter to skip pathway-specific plots: ")
  
  # Check if selected pathway is in the pathways of the filtered object
  if (nzchar(selected_pathway) && selected_pathway %in% cellchat_filtered_for_plotting@netP$pathways) { # <-- MODIFICATION HERE
    message(paste0("Generating CellChat bubble plot for pathway: ", selected_pathway))
    pdf(file.path(data_dir, paste0("cellchat_netVisual_bubble_", selected_pathway, "_selected_types.pdf")), width = 15, height = 12) # Added suffix
    # Use the filtered CellChat object for netVisual_bubble
    p_bubble <- netVisual_bubble(cellchat_filtered_for_plotting, # <-- MODIFICATION HERE (using filtered object)
                                 signaling = selected_pathway,
                                 remove.isolate = FALSE) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 8, vjust = 1, margin = ggplot2::margin(t = 2)),
        axis.text.y = ggplot2::element_text(size = 8),
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
      ) +
      ggplot2::ggtitle(paste0("CellChat Communication Probabilities: ", selected_pathway, " Pathway (Selected Cell Types)")) # Added suffix
    print(p_bubble)
    dev.off()
    message(paste0("CellChat bubble plot for '", selected_pathway, "' saved for selected cell types."))
    
    # Chord Plot for Selected Pathway
    message(paste0("Generating CellChat chord plot for pathway: ", selected_pathway))
    pdf(file.path(data_dir, paste0("cellchat_netVisual_chord_", selected_pathway, "_selected_types.pdf")), width = 14, height = 14) # Added suffix
    circlize::circos.clear()
    circlize::circos.par(gap.degree = 5, track.margin = c(0.01, 0.01), start.degree = 90, track.height = 0.05)
    graphics::par(cex = 0.4, mar = c(1, 1, 1, 1)) # Adjust margins for the plot
    # Use the filtered CellChat object for netVisual_chord_cell
    netVisual_chord_cell(cellchat_filtered_for_plotting, # <-- MODIFICATION HERE (using filtered object)
                         signaling = selected_pathway,
                         lab.cex = 3) # lab.cex adjusts label size
    dev.off()
    message(paste0("CellChat chord plot for '", selected_pathway, "' saved for selected cell types."))
    
  } else if (nzchar(selected_pathway)) {
    warning(paste0("Selected CellChat pathway '", selected_pathway, "' not found in available pathways. Skipping pathway visualizations."))
  } else {
    message("No CellChat pathway selected, skipping pathway visualizations.")
  }
  
  message("CellChat analysis and visualizations complete.")
}
current_pipeline_step_marker <- current_pipeline_step_marker + 1

# --- CellPhoneDB Analysis and Plotting Module ---

cat("--- CellphoneDB Analysis ---\n")

cellphonedb_io_for_current_run <- file.path(getwd(), "current_dataset_cellphonedb_io")
cellphonedb_output_dir_default <- file.path(cellphonedb_io_for_current_run, "cellphonedb_output")

if (!dir.exists(cellphonedb_io_for_current_run)) {
  cat("Creating CellphoneDB I/O directory: '", cellphonedb_io_for_current_run, "'\n")
  dir.create(cellphonedb_io_for_current_run, recursive = TRUE)
}

cat("Please select the directory to save/load CellPhoneDB output files.\n")
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  selected_output_folder <- rstudioapi::selectDirectory(caption = "Select Directory for CellPhoneDB Output Files")
} else {
  message("rstudioapi not available. Please manually set cellphonedb_output_dir.")
  selected_output_folder <- readline(prompt = "Enter path to CellPhoneDB output folder: ")
}

if (!is.null(selected_output_folder) && selected_output_folder != "") {
  cellphonedb_output_dir <- selected_output_folder
} else {
  cat("No output directory selected or selection cancelled.\n")
  cat("Using default output directory: '", cellphonedb_output_dir_default, "'\n")
  cellphonedb_output_dir <- cellphonedb_output_dir_default
}
cellphonedb_output_dir <- trimws(cellphonedb_output_dir)

if (!dir.exists(cellphonedb_output_dir)) {
  cat("Creating CellphoneDB output directory: '", cellphonedb_output_dir, "'\n")
  dir.create(cellphonedb_output_dir, recursive = TRUE)
}

run_cellphonedb_analysis <- TRUE
means_file_out <- NULL
pvalues_file_out <- NULL
deconvoluted_file_out <- NULL

existing_means_files <- list.files(cellphonedb_output_dir, pattern = "statistical_analysis_means_.*\\.txt", full.names = TRUE)
existing_pvalues_files <- list.files(cellphonedb_output_dir, pattern = "statistical_analysis_pvalues_.*\\.txt", full.names = TRUE)
existing_deconvoluted_files <- list.files(cellphonedb_output_dir, pattern = "statistical_analysis_deconvoluted_.*\\.txt", full.names = TRUE)

if (length(existing_means_files) > 0 && length(existing_pvalues_files) > 0 && length(existing_deconvoluted_files) > 0) {
  latest_means_file <- sort(existing_means_files, decreasing = TRUE)[1]
  latest_pvalues_file <- sort(existing_pvalues_files, decreasing = TRUE)[1]
  latest_deconvoluted_file <- sort(existing_deconvoluted_files, decreasing = TRUE)[1]
  
  cat("\nCellPhoneDB output files found in '", cellphonedb_output_dir, "'.\n")
  cat(paste0("Most recent means file: '", basename(latest_means_file), "'\n"))
  cat(paste0("Most recent pvalues file: '", basename(latest_pvalues_file), "'\n"))
  cat(paste0("Most recent deconvoluted file: '", basename(latest_deconvoluted_file), "'\n"))
  user_reuse_output <- tolower(readline(prompt = "Do you want to use these existing outputs for plotting (yes/no)? "))
  
  if (user_reuse_output == "yes") {
    run_cellphonedb_analysis <- FALSE
    means_file_out <- latest_means_file
    pvalues_file_out <- latest_pvalues_file
    deconvoluted_file_out <- latest_deconvoluted_file
    cat("Using existing CellPhoneDB outputs. Skipping analysis generation.\n")
  } else {
    cat("You chose to run analysis again. Existing files will be overwritten by new output.\n")
    run_cellphonedb_analysis <- TRUE
  }
}

counts_file <- NULL
meta_file <- NULL
create_cellphonedb_inputs <- FALSE

user_choice_reuse_inputs <- tolower(readline(prompt = "Do you have existing CellPhoneDB input files (counts.txt and meta.txt) to reuse? (yes/no): "))

if (user_choice_reuse_inputs == "yes") {
  if (!requireNamespace("rstudioapi", quietly = TRUE) || !rstudioapi::isAvailable()) {
    cat("Warning: The 'rstudioapi' package is required for folder selection. Falling back to creating new input files.\n")
    create_cellphonedb_inputs <- TRUE
  } else {
    cat("Please select the directory containing your CellPhoneDB input files (counts.txt and meta.txt).\n")
    selected_input_folder <- rstudioapi::selectDirectory(caption = "Select Folder with CellphoneDB Input Files")
    
    if (!is.null(selected_input_folder) && selected_input_folder != "") {
      counts_file_candidate <- file.path(selected_input_folder, "counts.txt")
      meta_file_candidate <- file.path(selected_input_folder, "meta.txt")
      
      if (file.exists(counts_file_candidate) && file.exists(meta_file_candidate)) {
        counts_file <- counts_file_candidate
        meta_file <- meta_file_candidate
        cat("Reusing existing input files from: '", selected_input_folder, "'\n")
      } else {
        cat("Warning: 'counts.txt' or 'meta.txt' not found in selected directory. Falling back to creating new input files.\n")
        create_cellphonedb_inputs <- TRUE
      }
    } else {
      cat("No folder selected or selection cancelled. Creating new input files.\n")
      create_cellphonedb_inputs <- TRUE
    }
  }
} else {
  create_cellphonedb_inputs <- TRUE
}

if (create_cellphonedb_inputs) {
  counts_file <- file.path(cellphonedb_io_for_current_run, "counts.txt")
  meta_file <- file.path(cellphonedb_io_for_current_run, "meta.txt")
  
  cat("Generating counts.txt from Seurat object (RNA assay, 'data' slot)....\n")
  cellphonedb_counts <- as.matrix(Seurat::GetAssayData(current_seurat_object, assay = "RNA", slot = "data"))
  cellphonedb_counts_df <- cbind(gene = rownames(cellphonedb_counts), as.data.frame(cellphonedb_counts))
  write.table(cellphonedb_counts_df, file = counts_file, sep = "\t", quote = FALSE, row.names = FALSE)
  
  cat("Generating meta.txt from Seurat object metadata...\n")
  cellphonedb_meta <- data.frame(
    Cell = colnames(current_seurat_object),
    cell_type = current_seurat_object@meta.data$scATOMIC_prediction
  )
  write.table(cellphonedb_meta, file = meta_file, sep = "\t", quote = FALSE, row.names = FALSE)
  cat("New input files 'counts.txt' and 'meta.txt' generated in: '", cellphonedb_io_for_current_run, "'\n")
} else {
  cat("Using existing input files directly.\n")
}

if (run_cellphonedb_analysis) {
  cat("Preparing to run CellPhoneDB v5.0.1 analysis via Python API...\n")
  
  python_path_env <- Sys.getenv("CELLPHONEDB_PYTHON_PATH")
  python_path_opt <- getOption("CELLPHONEDB_PYTHON_PATH")
  
  if (is.null(python_path_opt) && python_path_env == "") {
    stop("Error: CELLPHONEDB_PYTHON_PATH is not set. Please set it via Sys.setenv('CELLPHONEDB_PYTHON_PATH' = '/path/to/your/python') or options(CELLPHONEDB_PYTHON_PATH = '/path/to/your/python') before running this block. This should point to the Python executable within your CellPhoneDB conda environment.")
  }
  
  if (python_path_env != "") {
    reticulate::use_python(python_path_env, required = TRUE)
    cat(paste0("Reticulate is using Python from system env var: ", reticulate::py_discover_config()$python, "\n"))
  } else if (!is.null(python_path_opt)) {
    reticulate::use_python(python_path_opt, required = TRUE)
    cat(paste0("Reticulate is using Python from R option: ", reticulate::py_discover_config()$python, "\n"))
  }
  
  cellphonedb_method_func <- NULL
  tryCatch({
    cellphonedb_method_func <- reticulate::import("cellphonedb.src.core.methods.cpdb_statistical_analysis_method")
    cat("CellPhoneDB statistical analysis method imported successfully.\n")
  }, error = function(e) {
    stop(paste0("Error: Could not import CellPhoneDB Python API. Ensure 'cellphonedb' is installed and verify the correct import path for v5.0.1. Details: ", e$message))
  })
  
  cat("Please select the CellPhoneDB database file (e.g., 'cellphonedb.zip' or 'cellphonedb.json').\n")
  cat("This is typically found in your CellPhoneDB conda environment, for example:\n")
  cat("  'C:/miniconda3/envs/cellphonedb_env/Lib/site-packages/cellphonedb/data/'\n")
  
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    cpdb_file_path <- rstudioapi::selectFile(
      caption = "Select CellPhoneDB Database File",
      path = "C:/miniconda3/envs/cellphonedb_env/Lib/site-packages/cellphonedb/data/",
      filter = "CellPhoneDB Database Files (*.zip, *.json)"
    )
  } else {
    message("rstudioapi not available. Please manually set cpdb_file_path.")
    cpdb_file_path <- readline(prompt = "Enter path to CellPhoneDB database file: ")
  }
  
  if (is.null(cpdb_file_path) || cpdb_file_path == "") {
    stop("No CellPhoneDB database file selected. CellPhoneDB analysis aborted.")
  }
  cpdb_file_path <- trimws(cpdb_file_path)
  cat("Selected CellPhoneDB database file: '", cpdb_file_path, "'\n")
  
  cat("Initiating CellPhoneDB statistical analysis...\n")
  tryCatch({
    meta_file_clean <- trimws(meta_file)
    counts_file_clean <- trimws(counts_file)
    
    cellphonedb_method_func$call(
      cpdb_file_path = normalizePath(cpdb_file_path, winslash = "/"),
      meta_file_path = normalizePath(meta_file_clean, winslash = "/"),
      counts_file_path = normalizePath(counts_file_clean, winslash = "/"),
      counts_data = 'hgnc_symbol',
      output_path = normalizePath(cellphonedb_output_dir, winslash = "/"),
      iterations = as.integer(1000),
      threads = as.integer(7),
      debug_seed = as.integer(-1),
      result_precision = as.integer(3)
    )
    cat("CellPhoneDB statistical analysis completed.\n")
  }, error = function(e) {
    cat("Error during CellPhoneDB statistical analysis: ", e$message, "\n")
    stop("CellPhoneDB analysis failed. Please check the Python environment, input files, and CellPhoneDB installation.")
  })
  
  means_files_after_run <- list.files(cellphonedb_output_dir, pattern = "statistical_analysis_means_.*\\.txt", full.names = TRUE)
  pvalues_files_after_files <- list.files(cellphonedb_output_dir, pattern = "statistical_analysis_pvalues_.*\\.txt", full.names = TRUE)
  deconvoluted_files_after_run <- list.files(cellphonedb_output_dir, pattern = "statistical_analysis_deconvoluted_.*\\.txt", full.names = TRUE)
  
  if (length(means_files_after_run) == 0 || length(pvalues_files_after_files) == 0 || length(deconvoluted_files_after_run) == 0) {
    stop("CellPhoneDB output files not found after analysis. Please check logs.")
  }
  means_file_out <- sort(means_files_after_run, decreasing = TRUE)[1]
  pvalues_file_out <- sort(pvalues_files_after_files, decreasing = TRUE)[1]
  deconvoluted_file_out <- sort(deconvoluted_files_after_run, decreasing = TRUE)[1]
  cat("CellPhoneDB analysis completed. New output files will be used for plotting.\n")
  
} else {
  cat("Skipping CellPhoneDB analysis. Proceeding with plotting using existing outputs.\n")
}

cat("\n--- Generating CellPhoneDB Plots ---\n")

if (!requireNamespace("ktplots", quietly = TRUE)) {
  install.packages("ktplots")
}
library(ktplots)

if (!requireNamespace("viridis", quietly = TRUE)) {
  install.packages("viridis")
}
library(viridis)

if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
  BiocManager::install("SingleCellExperiment")
}
library(SingleCellExperiment)

if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
}
if (!requireNamespace("tidyr", quietly = TRUE)) {
  install.packages("tidyr")
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
  BiocManager::install("ComplexHeatmap")
}
if (!requireNamespace("ggraph", quietly = TRUE)) {
  install.packages("ggraph")
}
if (!requireNamespace("igraph", quietly = TRUE)) {
  install.packages("igraph")
}
library(dplyr)
library(tidyr)
library(ggplot2)
library(ComplexHeatmap)
library(ggraph)
library(igraph)
library(circlize)

cat("Loading CellPhoneDB means, pvalues, and deconvoluted files into R...\n")
means_df <- read.delim(means_file_out, sep = "\t", header = TRUE, check.names = FALSE)
pvalues_df <- read.delim(pvalues_file_out, sep = "\t", header = TRUE, check.names = FALSE)
deconvoluted_df <- read.delim(deconvoluted_file_out, sep = "\t", header = TRUE, check.names = FALSE)
cat("Files loaded successfully.\n")

metadata_cols_to_exclude <- c(
  "id_cp_interaction", "interacting_pair", "partner_a", "partner_b",
  "gene_a", "gene_b", "secreted", "receptor_a", "receptor_b",
  "annotation_strategy", "is_integrin", "directionality", "classification"
)

cellphonedb_meta_for_mapping <- read.delim(meta_file, sep = "\t", header = TRUE, check.names = FALSE)

all_cell_types <- sort(unique(as.character(current_seurat_object@meta.data$scATOMIC_prediction)))
numeric_id_to_cell_name_map <- setNames(all_cell_types, as.character(0:(length(all_cell_types)-1)))

interaction_cols_raw <- setdiff(colnames(means_df), metadata_cols_to_exclude)

new_interaction_col_names <- sapply(interaction_cols_raw, function(col_name) {
  if (grepl("\\|", col_name)) {
    parts <- strsplit(col_name, "\\|")[[1]]
    id1 <- parts[1]
    id2 <- parts[2]
    name1 <- numeric_id_to_cell_name_map[id1]
    name2 <- numeric_id_to_cell_name_map[id2]
    if (!is.na(name1) && !is.na(name2)) {
      return(paste0(name1, "|", name2))
    } else {
      return(col_name)
    }
  } else {
    return(col_name)
  }
})

colnames(means_df)[match(interaction_cols_raw, colnames(means_df))] <- new_interaction_col_names
colnames(pvalues_df)[match(interaction_cols_raw, colnames(pvalues_df))] <- new_interaction_col_names

cat("\n--- Cell Type Selection for Plotting ---\n")
cat("Available Cell Types:\n")

num_cols <- 1
max_name_length <- max(nchar(all_cell_types))
col_width <- max_name_length + nchar(as.character(length(all_cell_types))) + 5

for (k in seq_along(all_cell_types)) {
  item_string <- sprintf("[%-3d] %-*s", k, col_width, all_cell_types[k])
  cat(item_string)
  if (k %% num_cols == 0 || k == length(all_cell_types)) {
    cat("\n")
  }
}

user_selection_raw <- readline(prompt = "Enter the numbers of cell types to include in plots (e.g., 1,5,10 or 4-8), or type 'all' to select all: ")

selected_indices <- numeric(0)
user_input_cleaned <- tolower(trimws(user_selection_raw))

if (user_input_cleaned == "all") {
  selected_indices <- seq_along(all_cell_types)
} else if (user_input_cleaned != "") {
  # Split by comma first to handle mixed input like "1,4-8,10"
  input_parts <- strsplit(user_input_cleaned, ",")[[1]]
  
  for (part in input_parts) {
    if (grepl("-", part)) {
      # It's a range like "4-8"
      range_parts <- as.numeric(strsplit(part, "-")[[1]])
      if (length(range_parts) == 2 && !is.na(range_parts[1]) && !is.na(range_parts[2])) {
        selected_indices <- c(selected_indices, seq(range_parts[1], range_parts[2]))
      } else {
        warning(paste0("Invalid range format: ", part, ". Skipping."))
      }
    } else {
      # It's a single number
      num <- as.numeric(part)
      if (!is.na(num)) {
        selected_indices <- c(selected_indices, num)
      } else {
        warning(paste0("Invalid number format: ", part, ". Skipping."))
      }
    }
  }
  selected_indices <- unique(selected_indices) # Remove duplicates
}

if (any(is.na(selected_indices)) || length(selected_indices) == 0 || any(selected_indices < 1) || any(selected_indices > length(all_cell_types))) {
  stop("Invalid selection. Please enter valid numbers from the list, separated by commas, or type 'all'.")
}

chosen_cell_types <- all_cell_types[selected_indices]
cat("Selected Cell Types for Plotting:\n")
print(chosen_cell_types)

if (length(chosen_cell_types) == 0) {
  stop("No cell types selected. Aborting plotting steps.")
}

filtered_interaction_cols <- c()
for (col_name in new_interaction_col_names) {
  if (grepl("\\|", col_name)) {
    parts <- strsplit(col_name, "\\|")[[1]]
    if (parts[1] %in% chosen_cell_types && parts[2] %in% chosen_cell_types) {
      filtered_interaction_cols <- c(filtered_interaction_cols, col_name)
    }
  }
}

means_df_filtered_for_plotting <- means_df %>% dplyr::select(all_of(c(metadata_cols_to_exclude, filtered_interaction_cols)))
pvalues_df_filtered_for_plotting <- pvalues_df %>% dplyr::select(all_of(c(metadata_cols_to_exclude, filtered_interaction_cols)))

interaction_cols_means_filtered <- setdiff(colnames(means_df_filtered_for_plotting), metadata_cols_to_exclude)
interaction_cols_pvalues_filtered <- setdiff(colnames(pvalues_df_filtered_for_plotting), metadata_cols_to_exclude)

means_df_processed <- means_df_filtered_for_plotting %>%
  dplyr::mutate(across(all_of(interaction_cols_means_filtered), as.numeric))

pvalues_df_processed <- pvalues_df_filtered_for_plotting %>%
  dplyr::mutate(across(all_of(interaction_cols_pvalues_filtered), as.numeric))

plot_df <- means_df_processed %>%
  tidyr::pivot_longer(
    cols = all_of(interaction_cols_means_filtered),
    names_to = "cell_pair",
    values_to = "mean"
  )

pvalues_long <- pvalues_df_processed %>%
  tidyr::pivot_longer(
    cols = all_of(interaction_cols_pvalues_filtered),
    names_to = "cell_pair",
    values_to = "pvalue"
  )

plot_df <- plot_df %>%
  dplyr::left_join(pvalues_long, by = c("id_cp_interaction", "interacting_pair", "cell_pair")) %>%
  dplyr::filter(!is.na(mean) & !is.na(pvalue)) %>%
  dplyr::mutate(neg_log10_pvalue = -log10(pvalue + .Machine$double.eps))

plot_df_significant <- plot_df %>% dplyr::filter(pvalue < 0.05)

cat("\n--- Interaction Selection for Plotting ---\n")
all_interactions <- sort(unique(plot_df_significant$interacting_pair))

selected_interactions <- NULL

if (length(all_interactions) == 0) {
  cat("No significant interactions found for the selected cell types. Skipping interaction selection.\n")
} else {
  cat("Available Interactions:\n")
  num_cols_interactions <- 1
  max_name_length_interactions <- max(nchar(all_interactions))
  col_width_interactions <- max_name_length_interactions + nchar(as.character(length(all_interactions))) + 5
  
  for (k_int in seq_along(all_interactions)) {
    item_string_int <- sprintf("[%-3d] %-*s", k_int, col_width_interactions, all_interactions[k_int])
    cat(item_string_int)
    if (k_int %% num_cols_interactions == 0 || k_int == length(all_interactions)) {
      cat("\n")
    }
  }
  
  user_interaction_selection_raw <- readline(prompt = "Enter the numbers of interactions to include (e.g., 1,5,10 or 4-8), or type 'all' to select all: ")
  
  selected_interaction_indices <- numeric(0)
  user_interaction_input_cleaned <- tolower(trimws(user_interaction_selection_raw))
  
  if (user_interaction_input_cleaned == "all") {
    selected_interaction_indices <- seq_along(all_interactions)
  } else if (user_interaction_input_cleaned != "") {
    input_parts_int <- strsplit(user_interaction_input_cleaned, ",")[[1]]
    
    for (part_int in input_parts_int) {
      if (grepl("-", part_int)) {
        range_parts_int <- as.numeric(strsplit(part_int, "-")[[1]])
        if (length(range_parts_int) == 2 && !is.na(range_parts_int[1]) && !is.na(range_parts_int[2])) {
          selected_interaction_indices <- c(selected_interaction_indices, seq(range_parts_int[1], range_parts_int[2]))
        } else {
          warning(paste0("Invalid range format: ", part_int, ". Skipping."))
        }
      } else {
        num_int <- as.numeric(part_int)
        if (!is.na(num_int)) {
          selected_interaction_indices <- c(selected_interaction_indices, num_int)
        } else {
          warning(paste0("Invalid number format: ", part_int, ". Skipping."))
        }
      }
    }
    selected_interaction_indices <- unique(selected_interaction_indices) # Remove duplicates
  }
  
  if (any(is.na(selected_interaction_indices)) || length(selected_interaction_indices) == 0 || any(selected_interaction_indices < 1) || any(selected_interaction_indices > length(all_interactions))) {
    cat("Invalid interaction selection. Proceeding with all significant interactions.\n")
    selected_interactions <- all_interactions
  } else {
    selected_interactions <- all_interactions[selected_interaction_indices]
    cat("Selected Interactions for Plotting:\n")
    print(selected_interactions)
  }
  
  if (length(selected_interactions) == 0) {
    cat("No interactions selected. Plots might be empty.\n")
  }
}

if (!is.null(selected_interactions) && length(selected_interactions) > 0) {
  plot_df_significant <- plot_df_significant %>%
    dplyr::filter(interacting_pair %in% selected_interactions)
  cat(paste0("Filtered plot_df_significant to ", nrow(plot_df_significant), " rows based on interaction selection.\n"))
} else {
  cat("No specific interactions selected, proceeding with all significant interactions for the chosen cell types (or none if none were found).\n")
}

cat("Generating Consolidated Dot Plot for all selected cell type interactions...\n")
plot_filename_dot_consolidated <- file.path(cellphonedb_output_dir, "cellphonedb_consolidated_dot_plot.pdf")

if (nrow(plot_df_significant) > 0) {
  tryCatch({
    all_cell_pair_combinations <- c()
    for (ct1 in chosen_cell_types) {
      for (ct2 in chosen_cell_types) {
        all_cell_pair_combinations <- c(all_cell_pair_combinations, paste0(ct1, "|", ct2))
      }
    }
    ordered_cell_pairs_in_data <- intersect(all_cell_pair_combinations, unique(plot_df_significant$cell_pair))
    plot_df_significant$cell_pair <- factor(plot_df_significant$cell_pair, levels = ordered_cell_pairs_in_data)
    plot_df_significant$interacting_pair <- factor(plot_df_significant$interacting_pair,
                                                   levels = unique(plot_df_significant$interacting_pair[order(plot_df_significant$interacting_pair)]))
    cpdb_consolidated_dot_plot <- ggplot(plot_df_significant, aes(x = cell_pair, y = interacting_pair)) +
      geom_point(aes(size = neg_log10_pvalue, color = mean)) +
      scale_size_area(max_size = 8, name = "-log10(p-value)", breaks = c(1.3, 2, 3, 5, 10)) +
      scale_color_viridis_c(option = "B", name = "Mean Expression") +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
        axis.text.y = element_text(size = 8),
        axis.title = element_blank(),
        panel.grid.major = element_line(color = "lightgray", linewidth = 0.2),
        panel.grid.minor = element_blank(),
        legend.position = "right",
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5)
      ) +
      labs(
        title = "Consolidated CellPhoneDB Significant Interactions",
        subtitle = "Dot size represents significance (-log10(p-value)), color represents mean expression"
      )
    ggsave(plot_filename_dot_consolidated, cpdb_consolidated_dot_plot, width = 14, height = 18)
    cat(paste0(" Saved: '", plot_filename_dot_consolidated, "'\n"))
  }, error = function(e) {
    cat(paste0(" Error generating Consolidated Dot Plot: ", e$message, "\n"))
  })
} else {
  cat(" No significant interactions found for the selected cell types to generate a consolidated dot plot.\n")
}
cat("Consolidated dot plot generation complete.\n")

cat("Generating Heatmap manually using ComplexHeatmap...\n")
plot_filename_heatmap <- file.path(cellphonedb_output_dir, "cellphonedb_manual_heatmap.pdf")

p_value_threshold <- 0.05
cat(paste0("Attempting to save manual heatmap to: '", plot_filename_heatmap, "'\n"))
tryCatch({
  means_df_for_heatmap <- means_df_filtered_for_plotting
  pvalues_df_for_heatmap <- pvalues_df_filtered_for_plotting
  if (!is.null(selected_interactions) && length(selected_interactions) > 0) {
    means_df_for_heatmap <- means_df_for_heatmap %>% dplyr::filter(interacting_pair %in% selected_interactions)
    pvalues_df_for_heatmap <- pvalues_df_for_heatmap %>% dplyr::filter(interacting_pair %in% selected_interactions)
  }
  metadata_cols_for_cpdb <- c("id_cp_interaction", "interacting_pair", "partner_a", "partner_b", "gene_a", "gene_b", "secreted", "receptor_a", "receptor_b", "annotation_strategy", "is_integrin", "directionality", "classification")
  metadata_cols_present <- intersect(metadata_cols_for_cpdb, colnames(means_df_for_heatmap))
  interaction_cols <- setdiff(colnames(means_df_for_heatmap), metadata_cols_present)
  
  if (length(interaction_cols) == 0) {
    stop("No interaction columns found in means_df_for_heatmap after excluding metadata columns. Cannot create heatmap.")
  }
  
  means_matrix <- as.matrix(means_df_for_heatmap[, interaction_cols])
  rownames(means_matrix) <- means_df_for_heatmap$interacting_pair
  pvalues_matrix <- as.matrix(pvalues_df_for_heatmap[, interaction_cols])
  rownames(pvalues_matrix) <- pvalues_df_for_heatmap$interacting_pair
  
  ordered_interaction_cols <- c()
  for (ct1 in chosen_cell_types) {
    for (ct2 in chosen_cell_types) {
      pair_col <- paste0(ct1, "|", ct2)
      if (pair_col %in% interaction_cols) {
        ordered_interaction_cols <- c(ordered_interaction_cols, pair_col)
      }
    }
  }
  ordered_interaction_cols <- intersect(ordered_interaction_cols, colnames(means_matrix))
  means_matrix <- means_matrix[, ordered_interaction_cols, drop = FALSE]
  pvalues_matrix <- pvalues_matrix[, ordered_interaction_cols, drop = FALSE]
  
  if (nrow(means_matrix) == 0 || ncol(means_matrix) == 0) {
    stop("Means matrix is empty after filtering. Cannot generate heatmap.")
  }
  if (nrow(pvalues_matrix) == 0 || ncol(pvalues_matrix) == 0) {
    stop("P-values matrix is empty after filtering. Cannot generate heatmap.")
  }
  
  col_fun <- circlize::colorRamp2(
    c(min(means_matrix, na.rm = TRUE), max(means_matrix, na.rm = TRUE)),
    c(viridis::viridis_pal(option = "A")(2)[1], viridis::viridis_pal(option = "A")(2)[2])
  )
  
  cell_fun = function(j, i, x, y, width, height, fill) {
    p_value <- pvalues_matrix[i, j]
    if (!is.na(p_value) && p_value < p_value_threshold) {
      grid::grid.text("*", x, y, gp = grid::gpar(fontsize = 16, col = "black"))
    }
  }
  
  cpdb_heatmap <- ComplexHeatmap::Heatmap(
    means_matrix,
    name = "Mean Expression",
    col = col_fun,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    column_names_rot = 90,
    heatmap_legend_param = list(title = "Mean Expression"),
    cell_fun = cell_fun,
    row_names_gp = grid::gpar(fontsize = 8),
    column_names_gp = grid::gpar(fontsize = 8),
    width = grid::unit(ncol(means_matrix) * 0.5, "cm"),
    height = grid::unit(nrow(means_matrix) * 0.5, "cm"),
    row_title = "Interacting Pair",
    column_title = "Cell Pair",
    column_title_gp = grid::gpar(fontsize = 12, fontface = "bold"),
    row_title_gp = grid::gpar(fontsize = 12, fontface = "bold")
  )
  
  pdf(plot_filename_heatmap, width = 8 + ncol(means_matrix) * 0.1, height = 8 + nrow(means_matrix) * 0.1)
  ComplexHeatmap::draw(cpdb_heatmap, heatmap_legend_side = "right", annotation_legend_side = "right")
  dev.off()
  cat(paste0(" Saved: '", plot_filename_heatmap, "'\n"))
  
}, error = function(e) {
  cat(paste0(" Error generating Heatmap: ", e$message, "\n"))
})
cat("Heatmap generation complete.\n")


cat("Generating Network Plot...\n")
plot_filename_network <- file.path(cellphonedb_output_dir, "cellphonedb_network_plot.pdf")

tryCatch({
  network_df <- plot_df_significant %>%
    tidyr::separate(cell_pair, into = c("source_cell", "target_cell"), sep = "\\|") %>%
    dplyr::group_by(source_cell, target_cell) %>%
    dplyr::summarise(num_interactions = n(), .groups = "drop") %>%
    dplyr::ungroup()
  
  if (nrow(network_df) == 0) {
    stop("No significant interactions to plot in the network. Network plot skipped.")
  }
  
  graph_obj <- igraph::graph_from_data_frame(network_df, directed = TRUE)
  
  node_interactions <- network_df %>%
    tidyr::pivot_longer(c(source_cell, target_cell), names_to = "role", values_to = "cell_type") %>%
    dplyr::group_by(cell_type) %>%
    dplyr::summarise(total_interactions = sum(num_interactions), .groups = "drop")
  
  igraph::V(graph_obj)$total_interactions <- node_interactions$total_interactions[match(igraph::V(graph_obj)$name, node_interactions$cell_type)]
  igraph::V(graph_obj)$total_interactions[is.na(igraph::V(graph_obj)$total_interactions)] <- 0
  
  cpdb_network_plot <- ggraph::ggraph(graph_obj, layout = "fr") +
    ggraph::geom_edge_fan(aes(edge_width = num_interactions, color = num_interactions),
                          arrow = grid::arrow(length = grid::unit(2, 'mm'), type = "closed"),
                          end_cap = ggraph::circle(5, 'mm')) +
    ggraph::scale_edge_width_continuous(range = c(0.5, 3), name = "No. of Significant Interactions") +
    ggraph::scale_edge_color_continuous(low = "yellow", high = "purple", name = "No. of Significant Interactions") +
    ggraph::geom_node_point(color = "lightblue", size = 10) +
    ggraph::geom_node_text(aes(label = name), repel = TRUE, size = 4) +
    theme_void() +
    labs(
      title = "Cell-Cell Communication Network (Significant Interactions)",
      subtitle = "Nodes are cell types; Edge width and color represent number of significant interactions"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position = "right"
    )
  
  ggplot2::ggsave(plot_filename_network, cpdb_network_plot, width = 12, height = 10)
  cat(paste0("Network plot saved to: '", plot_filename_network, "'\n"))
  
}, error = function(e) {
  cat(paste0("Error generating Network Plot: ", e$message, "\n"))
})


cat("Generating Chord Plot using ktplots::plot_cpdb4...\n")

plot_filename_chord <- file.path(cellphonedb_output_dir, "cellphonedb_chord_plot.pdf")

# Convert Seurat to SingleCellExperiment
sce_object <- as.SingleCellExperiment(current_seurat_object)

# --- CRITICAL FIX: Populate logcounts assay explicitly ---
logcounts(sce_object) <- as.matrix(Seurat::GetAssayData(
  current_seurat_object,
  assay = "RNA",
  slot = "data"
))

# Retrieve colData, modify the column, and reassign
current_coldata <- SingleCellExperiment::colData(sce_object)
current_coldata$scATOMIC_prediction <- as.factor(current_coldata$scATOMIC_prediction)
colData(sce_object) <- current_coldata

# Optional: Check if logcounts is correctly set
if (is.null(assay(sce_object, "logcounts"))) {
  stop("logcounts assay is NULL. Check RNA data transfer.")
}
# Generate all possible cell pairs as before
cell_pairs <- as.vector(outer(chosen_cell_types, chosen_cell_types, paste, sep = "|"))

# Get available pairs present in all three data frames
available_pairs <- intersect(
  intersect(colnames(means_df), colnames(pvalues_df)),
  colnames(deconvoluted_df)
)
available_pairs <- available_pairs[grepl("\\|", available_pairs)]

# Warn about missing pairs in deconvoluted_df
missing_in_deconv <- setdiff(cell_pairs, colnames(deconvoluted_df))
if(length(missing_in_deconv) > 0) {
  cat("Warning: The following cell pairs are missing in deconvoluted_df and will be skipped:\n")
  print(missing_in_deconv)
}

# Filter interaction_vector to only those present in means_df and available_pairs
interaction_vector <- as.character(plot_df_significant$interacting_pair)
interaction_vector <- interaction_vector[interaction_vector %in% means_df$interacting_pair]

# Only plot if valid pairs exist
valid_pairs <- cell_pairs[cell_pairs %in% available_pairs]
if(length(valid_pairs) == 0) stop("No valid cell pairs present in all data frames.")

tryCatch({
  cat(paste0("Attempting to save Chord Plot to: '", plot_filename_chord, "'\n"))
  cpdb_chord_plot <- ktplots::plot_cpdb4(
    scdata = sce_object,
    celltype_key = "scATOMIC_prediction",
    cell_type1 = "scATOMIC_prediction",
    cell_type2 = "scATOMIC_prediction",
    means = means_df,
    pvals = pvalues_df,
    deconvoluted = deconvoluted_df,
    interaction = interaction_vector
  )
  ggsave(plot_filename_chord, cpdb_chord_plot, width = 12, height = 10)
  cat(paste0("Chord plot saved to: '", plot_filename_chord, "'\n"))
}, error = function(e) {
  cat("Error generating Chord Plot: ", e$message, "\n")
})


cat("All CellPhoneDB analysis and plotting steps complete.\n")
# --- 17. Final Data Saving ---
if (current_pipeline_step_marker >= start_step_numeric) {
  message("--- Step ", current_pipeline_step_marker, ": Final Data Saving ---")
  user_save_choice <- tolower(readline(prompt = "Do you want to save the final Seurat object? (yes/no): "))
  
  if (user_save_choice == "yes") {
    message("Saving final Seurat object...")
    # Define filename for the final object
    final_seurat_object_filename <- file.path(data_dir, "final_seurat_object_pipeline.rds")
    saveRDS(current_seurat_object, file = final_seurat_object_filename)
    message(paste0("Final Seurat object saved to: '", final_seurat_object_filename, "'."))
    message("Pipeline complete. You can now close RStudio.")
  } else {
    message("Skipping final Seurat object save. Pipeline complete.")
  }
}
# No increment needed after the last step

message("\n--- Pipeline execution finished. ---")