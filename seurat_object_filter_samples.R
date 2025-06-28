library(Seurat)

# 1. Select and load the Seurat object (.rds)
cat("Please select your Seurat object (.rds) file to load.\n")
seurat_file <- file.choose()
seurat_obj <- readRDS(seurat_file)
cat("Seurat object loaded successfully!\n")

# 2. Show available sample IDs
cat("Unique sample IDs in $orig.ident:\n")
print(unique(seurat_obj$orig.ident))

# 3. Ask user for samples to exclude (comma-separated)
samples_input <- readline(prompt = "Enter sample IDs to exclude, separated by commas (e.g., 18P,23P,25P): ")
samples_to_exclude <- trimws(unlist(strsplit(samples_input, ",")))

# 4. Filter Seurat object
seurat_obj_filtered <- subset(seurat_obj, subset = !(orig.ident %in% samples_to_exclude))
cat("Filtered Seurat object now contains these samples:\n")
print(unique(seurat_obj_filtered$orig.ident))

# 5. Ask user for the output file name (no path, just the file name)
save_name <- readline(prompt = "Enter the file name to save the filtered Seurat object (e.g., seurat_obj_filtered.rds): ")

# 6. Save in the same directory as the original Seurat object
output_dir <- dirname(seurat_file)
save_path <- file.path(output_dir, save_name)
saveRDS(seurat_obj_filtered, file = save_path)
cat(paste("Filtered Seurat object saved as:", save_path, "\n"))
