#!/usr/bin/env Rscript
# Run MMUPHin batch correction across multiple gut microbiome cohorts.

suppressPackageStartupMessages({
  library(methods)
  library(mmuphin)
  library(vegan)
  library(ape)
  library(ggplot2)
})

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  opts <- list(
    feature_tables = NULL,
    metadata = NULL,
    output_dir = "data/processed/mmuphin"
  )

  for (arg in args) {
    if (grepl("^--feature-tables=", arg)) {
      opts$feature_tables <- sub("^--feature-tables=", "", arg)
    } else if (grepl("^--metadata=", arg)) {
      opts$metadata <- sub("^--metadata=", "", arg)
    } else if (grepl("^--output-dir=", arg)) {
      opts$output_dir <- sub("^--output-dir=", "", arg)
    }
  }

  if (is.null(opts$feature_tables) || is.null(opts$metadata)) {
    stop("Usage: Rscript run_mmuphin.R --feature-tables=table1.tsv,table2.tsv --metadata=meta.tsv --output-dir=...\n",
         call. = FALSE)
  }

  opts$feature_tables <- strsplit(opts$feature_tables, ",")[[1]]
  opts
}

load_feature_table <- function(path) {
  # Read a feature table with features as rows and sample IDs as columns.
  table <- read.table(path, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE, quote = "\"")
  as.matrix(table)
}

run_pcoa <- function(feature_table, metadata, title, output_file) {
  # Compute Bray-Curtis distance and create a PCoA plot colored by study.
  dist_mat <- vegdist(t(feature_table), method = "bray")
  pcoa_obj <- pcoa(dist_mat)
  coords <- as.data.frame(pcoa_obj$vectors[, 1:2])
  coords$sample_id <- rownames(coords)
  coords <- merge(coords, metadata, by.x = "sample_id", by.y = "sample_id", all.x = TRUE)

  p <- ggplot(coords, aes(x = Axis.1, y = Axis.2, color = study)) +
    geom_point(size = 3, alpha = 0.8) +
    ggtitle(title) +
    xlab(sprintf("PCoA 1 (%.1f%%)", 100 * pcoa_obj$values$Relative_eig[1])) +
    ylab(sprintf("PCoA 2 (%.1f%%)", 100 * pcoa_obj$values$Relative_eig[2])) +
    theme_minimal() +
    theme(legend.position = "right")

  ggsave(output_file, p, width = 8, height = 6)
}

main <- function() {
  opts <- parse_args()
  dir.create(opts$output_dir, recursive = TRUE, showWarnings = FALSE)

  # Load the provided metadata file.
  metadata <- read.table(opts$metadata, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  if (!"sample_id" %in% colnames(metadata) || !"study" %in% colnames(metadata)) {
    stop("Metadata must contain columns 'sample_id' and 'study'.")
  }

  # Load each feature table and combine them by sample ID.
  feature_tables <- lapply(opts$feature_tables, load_feature_table)
  merged_table <- do.call(cbind, feature_tables)

  # Ensure the combined feature table matches metadata samples.
  sample_ids <- colnames(merged_table)
  if (!all(sample_ids %in% metadata$sample_id)) {
    stop("Some sample IDs in the feature tables are missing from metadata.")
  }
  metadata <- metadata[match(sample_ids, metadata$sample_id), ]

  # Save the raw merged feature table for reference.
  raw_path <- file.path(opts$output_dir, "merged_feature_table_raw.tsv")
  write.table(merged_table, raw_path, sep = "\t", quote = FALSE, col.names = NA)

  # Visualize batch structure before correction.
  run_pcoa(merged_table, metadata, "PCoA before MMUPHin batch correction", file.path(opts$output_dir, "pcoa_before_correction.png"))

  # Run MMUPHin batch correction across studies using the study column.
  corrected_table <- adjust_batch(merged_table, metadata$study)

  # Visualize the corrected data to assess batch correction.
  run_pcoa(corrected_table, metadata, "PCoA after MMUPHin batch correction", file.path(opts$output_dir, "pcoa_after_correction.png"))

  # Save the batch-corrected feature table.
  corrected_path <- file.path(opts$output_dir, "batch_corrected_feature_table.tsv")
  write.table(corrected_table, corrected_path, sep = "\t", quote = FALSE, col.names = NA)

  message("Batch correction complete. Corrected feature table saved to:", corrected_path)
}

main()
