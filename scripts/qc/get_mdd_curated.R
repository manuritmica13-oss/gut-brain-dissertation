#!/usr/bin/env Rscript
## Download MDD datasets from curatedMetagenomicData and save tables

install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages('BiocManager', repos = 'https://cloud.r-project.org')
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}

if (!requireNamespace('BiocManager', quietly = TRUE)) {
  install.packages('BiocManager', repos = 'https://cloud.r-project.org')
}

install_if_missing('curatedMetagenomicData')
library(curatedMetagenomicData)

message('Fetching sample metadata from curatedMetagenomicData...')
sm <- tryCatch(curatedMetagenomicData::sampleMetadata(), error = function(e) {
  stop('Failed to retrieve sampleMetadata(): ', e$message)
})

dfsm <- as.data.frame(sm)

# search for depression / MDD keywords across all metadata fields
kw <- '(depress|major depressive|mdd)'
hit_rows <- apply(dfsm, 1, function(r) any(grepl(kw, as.character(r), ignore.case = TRUE)))
matched <- dfsm[hit_rows, , drop = FALSE]

if (nrow(matched) == 0) {
  message('No samples matching MDD/depression found in curatedMetagenomicData sampleMetadata().')
  quit(status = 0)
}

message(sprintf('Found %d samples matching MDD/depression in sampleMetadata()', nrow(matched)))

# try to determine dataset/study identifier column
possible_cols <- c('study_name', 'study', 'dataset', 'studyid', 'project')
dataset_col <- intersect(possible_cols, colnames(matched))
if (length(dataset_col) >= 1) {
  ds_col <- dataset_col[1]
  ds_names <- unique(as.character(matched[[ds_col]]))
} else {
  # fallback: parse rownames which often have the dataset prefix before ':'
  rn <- rownames(matched)
  ds_names <- unique(sub(':.*$', '', rn))
}

message(sprintf('Datasets identified (%d): %s', length(ds_names), paste(head(ds_names, 5), collapse = ', ')))

# prepare output directories
dir.create(file.path('data', 'raw'), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path('data', 'processed', 'mdd_curated'), recursive = TRUE, showWarnings = FALSE)

meta_out <- file.path('data', 'raw', 'curatedMetagenomicData_metadata.tsv')
write.table(matched, file = meta_out, sep = '\t', quote = FALSE, row.names = TRUE)
message('Wrote sample metadata to ', meta_out)

relab_list <- list()
path_list <- list()
meta_list <- list()

for (ds in ds_names) {
  message('Requesting dataset: ', ds)
  # attempt to download objects matching this dataset
  res <- tryCatch(
    curatedMetagenomicData(ds, dryrun = FALSE),
    error = function(e) {
      message('  failed to download ', ds, ': ', e$message);
      return(NULL)
    }
  )
  if (is.null(res)) next

  # res may be a list or a single object; coerce to list
  objs <- if (is.list(res) && length(res) > 1) res else list(res)

  for (obj in objs) {
    # extract sample metadata
    md <- NULL
    if (inherits(obj, 'ExpressionSet')) {
      md <- as.data.frame(Biobase::pData(obj))
      expr_mat <- Biobase::exprs(obj)
      assays_list <- list(exprs = expr_mat)
    } else if (inherits(obj, 'SummarizedExperiment')) {
      md <- as.data.frame(SummarizedExperiment::colData(obj))
      # get all assays
      assays_list <- SummarizedExperiment::assays(obj)
    } else {
      # unknown object type: try to coerce
      next
    }

    if (!is.null(md)) meta_list[[length(meta_list) + 1]] <- md

    # inspect assays and classify as taxa (metaphlan) or pathways (humann)
    for (an in names(assays_list)) {
      mat <- assays_list[[an]]
      if (is.null(dim(mat))) next
      rn <- rownames(mat)
      if (is.null(rn)) next
      # heuristics: taxa often contain '__' (metaphlan) or underscore names; pathways often contain '|' or ':' or 'PATH' words
      taxa_like <- mean(grepl('__|_', rn)) > 0.2
      pathway_like <- mean(grepl(':|\\||path|PATH|UNMAPPED', rn, ignore.case = TRUE)) > 0.01

      if (taxa_like && !pathway_like) {
        message('  found taxa assay: ', an, ' (', nrow(mat), ' features, ', ncol(mat), ' samples)')
        relab_list[[length(relab_list) + 1]] <- mat
      } else if (pathway_like && !taxa_like) {
        message('  found pathway assay: ', an, ' (', nrow(mat), ' features, ', ncol(mat), ' samples)')
        path_list[[length(path_list) + 1]] <- mat
      } else {
        # if uncertain, try both tests: if many rows with typical taxon patterns treat as taxa
        if (taxa_like) relab_list[[length(relab_list) + 1]] <- mat
        else if (pathway_like) path_list[[length(path_list) + 1]] <- mat
        else message('  assay ', an, ' could not be classified automatically (features sample: ', nrow(mat), ', ', ncol(mat), ')')
      }
    }
  }
}

merge_matrices <- function(mat_list) {
  if (length(mat_list) == 0) return(NULL)
  # convert each to data.frame with feature column
  df_list <- lapply(mat_list, function(m) {
    df <- as.data.frame(as.matrix(m))
    df$feature <- rownames(m)
    return(df)
  })
  # full join sequentially by feature
  merged <- df_list[[1]]
  for (i in seq(2, length(df_list))) {
    merged <- merge(merged, df_list[[i]], by = 'feature', all = TRUE)
  }
  rownames(merged) <- merged$feature
  merged$feature <- NULL
  # replace NA with 0
  merged[is.na(merged)] <- 0
  return(as.matrix(merged))
}

relab_merged <- merge_matrices(relab_list)
path_merged <- merge_matrices(path_list)

if (!is.null(relab_merged)) {
  dir.create(dirname(file.path('data', 'processed', 'mdd_curated', 'metaphlan_rel_ab.tsv')), recursive = TRUE, showWarnings = FALSE)
  write.table(relab_merged, file = file.path('data', 'processed', 'mdd_curated', 'metaphlan_rel_ab.tsv'), sep = '\t', quote = FALSE, col.names = NA)
  message('Wrote MetaPhlAn relative abundance matrix with ', nrow(relab_merged), ' features and ', ncol(relab_merged), ' samples')
} else {
  message('No MetaPhlAn-like assays found.')
}

if (!is.null(path_merged)) {
  write.table(path_merged, file = file.path('data', 'processed', 'mdd_curated', 'humann_pathway_ab.tsv'), sep = '\t', quote = FALSE, col.names = NA)
  message('Wrote HUMAnN pathway abundance matrix with ', nrow(path_merged), ' features and ', ncol(path_merged), ' samples')
} else {
  message('No HUMAnN-like assays found.')
}

# combine metadata and report counts of MDD vs controls
if (length(meta_list) > 0) {
  all_meta <- do.call(rbind, lapply(meta_list, function(x) { x }))
  # write combined metadata
  write.table(all_meta, file = file.path('data', 'raw', 'curatedMetagenomicData_metadata.tsv'), sep = '\t', quote = FALSE, row.names = TRUE)
  message('Wrote combined sample metadata to data/raw/curatedMetagenomicData_metadata.tsv')

  # try to find a column describing disease/condition
  cols <- colnames(all_meta)
  disease_col <- NULL
  for (c in cols) {
    if (any(grepl('disease|condition|diagnos', c, ignore.case = TRUE))) { disease_col <- c; break }
  }
  if (!is.null(disease_col)) {
    vals <- as.character(all_meta[[disease_col]])
    mdd_flag <- grepl(kw, vals, ignore.case = TRUE)
    ctrl_flag <- grepl('control|healthy|normal|none', vals, ignore.case = TRUE)
    n_mdd <- sum(mdd_flag, na.rm = TRUE)
    n_ctrl <- sum(ctrl_flag, na.rm = TRUE)
    message(sprintf('Samples with MDD/depression in metadata column %s: %d', disease_col, n_mdd))
    message(sprintf('Samples with control/healthy labels in metadata column %s: %d', disease_col, n_ctrl))
  } else {
    # fallback: use earlier matched sample list
    # count matched vs others in sampleMetadata
    matched_ids <- rownames(matched)
    n_mdd <- length(matched_ids)
    message(sprintf('Matched samples found (from sampleMetadata search): %d (controls count unknown)', n_mdd))
  }
} else {
  message('No sample-level metadata objects were extracted from the downloaded datasets.')
}

message('Done.')
