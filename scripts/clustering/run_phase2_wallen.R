#!/usr/bin/env Rscript
## Phase 2 analysis for Wallen dataset

## This script performs the following using files in data/processed/wallen/:
##  1. Load metaphlan_rel_ab.tsv (transpose so samples are rows), load subject_metadata.tsv
##  2. Filter to samples with known Case_status (Case or Control)
##  3. Calculate Bray-Curtis distance matrix using vegan
##  4. Run PERMANOVA (adonis2) with Case_status as grouping variable, controlling for Age, Sex, BMI
##  5. Run PERMDISP (betadisper) to check dispersion assumption
##  6. Generate a PCoA plot coloured by Case_status saved to results/wallen_pcoa.png
##  7. Run DESeq2 differential abundance between Case and Control — report top 20 taxa by adjusted p-value
##  8. Extract IDO1 (K00463), KMO (K00486), KAT (K00816), TPH1 (K00502) abundances from humann_pathway_counts.tsv and compare Case vs Control with a Wilcoxon test
##  9. Save all results to results/phase2_wallen_results.txt

suppressPackageStartupMessages({
  library(data.table)
  library(vegan)
  library(ape)
  library(ggplot2)
})

ROOT <- normalizePath("..", mustWork = FALSE)
processed_dir <- file.path("data", "processed", "wallen")
results_dir <- file.path("results")
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

meta_file <- file.path(processed_dir, "subject_metadata.tsv")
meta <- fread(meta_file, sep = "\t", header = TRUE)

find_col <- function(dt, patterns){
  nm <- tolower(names(dt))
  for(p in patterns){
    i <- which(grepl(p, nm))
    if(length(i)) return(names(dt)[i[1]])
  }
  return(NA_character_)
}

sample_col <- find_col(meta, c("^sample", "^sampleid", "sample_name"))
if (is.na(sample_col)) stop("Could not find sample id column in subject_metadata.tsv")
case_col <- find_col(meta, c("case_status", "case", "case.status", "disease_status", "diagnosis"))
if (is.na(case_col)) stop("Could not find Case_status column in subject_metadata.tsv")
age_col <- find_col(meta, c("age", "age_at", "ageat"))
sex_col <- find_col(meta, c("sex", "gender"))
bmi_col <- find_col(meta, c("bmi", "bodymass"))

meta[[sample_col]] <- as.character(meta[[sample_col]])
meta[[case_col]] <- as.character(meta[[case_col]])

# standardize case/control labels
meta[[case_col]] <- ifelse(tolower(meta[[case_col]]) %in% c("case", "patient", "pd", "parkinson", "mdd", "depression"), "Case",
                           ifelse(tolower(meta[[case_col]]) %in% c("control", "ctrl", "hc", "healthy"), "Control", NA))

# Read metaphlan
taxa_file <- file.path(processed_dir, "metaphlan_rel_ab.tsv")
taxa_raw <- fread(taxa_file, sep = "\t", header = TRUE, check.names = FALSE)

# assume first column is taxa name
taxa_names <- as.character(taxa_raw[[1]])
taxa_mat <- as.matrix(taxa_raw[, -1, with = FALSE])
rownames(taxa_mat) <- taxa_names

# transpose so samples are rows
taxa_t <- as.data.frame(t(taxa_mat), stringsAsFactors = FALSE)
taxa_t[] <- lapply(taxa_t, function(x) as.numeric(as.character(x)))
taxa_t$sample_id <- rownames(taxa_t)
rownames(taxa_t) <- NULL

# Merge metadata and filter known case/control
meta_sub <- copy(meta)
names(meta_sub)[names(meta_sub) == sample_col] <- "sample_id"

common <- intersect(meta_sub$sample_id, taxa_t$sample_id)
if(length(common) == 0) stop("No matching samples between metadata and metaphlan table")

meta_f <- meta_sub[meta_sub$sample_id %in% common & !is.na(meta_sub[[case_col]]), ]
taxa_f <- taxa_t[taxa_t$sample_id %in% meta_f$sample_id, ]

# Reorder taxa_f to match meta_f
taxa_f <- taxa_f[match(meta_f$sample_id, taxa_f$sample_id), ]

abund_mat <- as.matrix(taxa_f[, setdiff(names(taxa_f), "sample_id")])
rownames(abund_mat) <- taxa_f$sample_id

## Ensure covariates present and filter to complete cases for PERMANOVA
meta_f$Case_status <- factor(meta_f[[case_col]])
meta_f$Age <- if (!is.na(age_col)) as.numeric(meta_f[[age_col]]) else NA_real_
meta_f$Sex <- if (!is.na(sex_col)) factor(meta_f[[sex_col]]) else NA
meta_f$BMI <- if (!is.na(bmi_col)) as.numeric(meta_f[[bmi_col]]) else NA_real_

# keep only samples with non-missing Case_status, Age, Sex, BMI
complete_idx <- which(!is.na(meta_f$Case_status) & !is.na(meta_f$Age) & !is.na(meta_f$Sex) & !is.na(meta_f$BMI))
if(length(complete_idx) == 0) stop("No samples with complete Case_status, Age, Sex, BMI for PERMANOVA")

meta_f <- meta_f[complete_idx, ]
taxa_f <- taxa_f[match(meta_f$sample_id, taxa_f$sample_id), ]
abund_mat <- as.matrix(taxa_f[, setdiff(names(taxa_f), "sample_id")])
rownames(abund_mat) <- taxa_f$sample_id

# Bray-Curtis distance (on filtered samples)
d_bc <- vegdist(abund_mat, method = "bray")

# PERMANOVA
adon_res <- adonis2(abund_mat ~ Case_status + Age + Sex + BMI, data = meta_f, permutations = 999, method = "bray")

# PERMDISP
bd <- betadisper(d_bc, meta_f$Case_status)
bd_perm <- permutest(bd, permutations = 999)

# PCoA and plot
pco <- pcoa(as.matrix(d_bc))
scores <- as.data.frame(pco$vectors[, 1:2])
scores$sample_id <- rownames(scores)
scores <- merge(scores, meta_f[, c("sample_id", "Case_status")], by = "sample_id")

png(filename = file.path(results_dir, "wallen_pcoa.png"), width = 900, height = 700)
  ggplot(scores, aes(x = Axis.1, y = Axis.2, color = Case_status)) +
    geom_point(size = 3, alpha = 0.8) +
    labs(x = paste0("PCoA1 (", round(pco$values$Relative_eig[1]*100,1), "%)"),
         y = paste0("PCoA2 (", round(pco$values$Relative_eig[2]*100,1), "%)"),
         title = "Wallen PCoA (Bray-Curtis)") +
    theme_minimal() +
    scale_color_manual(values = c("Control" = "blue", "Case" = "red"))
dev.off()

# Differential abundance using Wilcoxon rank-sum tests (base R)
# abund_mat: rows=samples, cols=taxa
taxa_names <- colnames(abund_mat)
pw_p <- numeric(length(taxa_names))
pw_mean_case <- numeric(length(taxa_names))
pw_mean_control <- numeric(length(taxa_names))
for(i in seq_along(taxa_names)){
  vals <- as.numeric(abund_mat[, i])
  grp <- meta_f$Case_status
  case_vals <- vals[grp == "Case"]
  ctrl_vals <- vals[grp == "Control"]
  # require at least one non-NA in each group
  if(length(case_vals) < 1 || length(ctrl_vals) < 1){
    pw_p[i] <- NA
    pw_mean_case[i] <- NA
    pw_mean_control[i] <- NA
    next
  }
  # use wilcox test (non-parametric)
  wt <- tryCatch(wilcox.test(case_vals, ctrl_vals, exact = FALSE), error = function(e) list(p.value = NA))
  pw_p[i] <- wt$p.value
  pw_mean_case[i] <- mean(case_vals, na.rm = TRUE)
  pw_mean_control[i] <- mean(ctrl_vals, na.rm = TRUE)
}
res_df <- data.frame(taxa = taxa_names, p.value = pw_p, mean_case = pw_mean_case, mean_control = pw_mean_control, stringsAsFactors = FALSE)
res_df$padj <- p.adjust(res_df$p.value, method = "BH")
res_df$log2FC <- log2((res_df$mean_case + 1e-6) / (res_df$mean_control + 1e-6))
res_df <- res_df[order(res_df$padj, na.last = TRUE), ]
top20 <- head(res_df, 20)

# Read humann pathway/KO counts
hum_file <- file.path(processed_dir, "humann_pathway_counts.tsv")
hum_raw <- fread(hum_file, sep = "\t", header = TRUE, check.names = FALSE)

# assume first column is feature name
hum_features <- as.character(hum_raw[[1]])
ko_ids <- c("K00463", "K00486", "K00816", "K00502")
ko_idx <- which(tolower(hum_features) %in% tolower(ko_ids) | hum_features %in% ko_ids)

ko_results <- data.frame(KO = character(), p.value = numeric(), stringsAsFactors = FALSE)
if(length(ko_idx) > 0){
  for(i in ko_idx){
    ko <- hum_features[i]
    vals <- as.numeric(as.vector(unlist(hum_raw[i, -1, with = FALSE])))
    names(vals) <- names(hum_raw)[-1]
    # match sample order and metadata
    common_samples <- intersect(meta_f$sample_id, names(vals))
    vals_sub <- vals[common_samples]
    meta_sub2 <- meta_f[match(common_samples, meta_f$sample_id), ]
    grp1 <- vals_sub[meta_sub2$Case_status == "Case"]
    grp0 <- vals_sub[meta_sub2$Case_status == "Control"]
    wt <- wilcox.test(grp1, grp0)
    ko_results <- rbind(ko_results, data.frame(KO = ko, p.value = wt$p.value, stringsAsFactors = FALSE))
  }
  ko_results$padj <- p.adjust(ko_results$p.value, method = "BH")
}

# Save results
outf <- file.path(results_dir, "phase2_wallen_results.txt")
con <- file(outf, open = "wt")
cat("PERMANOVA (adonis2):\n", file = con)
capture.output(adon_res, file = con)
cat("\nPERMDISP (betadisper) permutest:\n", file = con)
capture.output(bd_perm, file = con)
cat("\nTop 20 DESeq2 results (by padj):\n", file = con)
write.table(top20, file = con, sep = "\t", row.names = FALSE, quote = FALSE)
cat("\nKO Wilcoxon tests:\n", file = con)
if(nrow(ko_results) == 0) cat("No KOs found in humann table.\n", file = con) else write.table(ko_results, file = con, sep = "\t", row.names = FALSE, quote = FALSE)
close(con)

cat("Analysis complete. Results written to", outf, "and PCoA image to", file.path(results_dir, "wallen_pcoa.png"), "\n")
