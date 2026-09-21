#!/usr/bin/env Rscript

# Audit-only, deterministic inspection of Hazzard2022. This script reads existing
# inputs/artifacts and writes evidence beside AUDIT.md; it never saves a production object.
suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
  library(digest)
})

audit_dir <- normalizePath("../audit/mancio_silva2022", mustWork = TRUE)
data_dir <- normalizePath("../counts/36525464", mustWork = TRUE)
runs <- c("SRR20710498", "SRR20710499", "SRR20710500", "SRR20710501")
suffix <- sub("^SRR20710", "", runs)
write_tsv <- function(x, name) {
  write.table(x, file.path(audit_dir, name), sep = "\t", quote = FALSE,
              row.names = FALSE, na = "NA")
}
`%notin%` <- function(x, y) !(x %in% y)

study <- readRDS("../hazzard2022.rds")
counts <- GetAssayData(study, assay = "RNA", layer = "counts")
meta <- study@meta.data
stopifnot(identical(colnames(counts), rownames(meta)), anyDuplicated(colnames(counts)) == 0L)

app <- readRDS("../../PlaViSca/data/cleaned_dataset.rds")$mr_data
final <- app[app$study_label == "Hazzard2022", , drop = FALSE]
stopifnot(anyDuplicated(rownames(final)) == 0L)

# MatrixMarket headers and barcode inventories, without loading the 6.8M-column raw
# matrices. All STARsolo raw barcode files are deliberately counted linewise.
inventory <- lapply(seq_along(runs), function(i) {
  root <- file.path(data_dir, paste0(runs[i], "_solo_out/Solo.out/GeneFull/raw"))
  con <- file(file.path(root, "matrix.mtx"), "rt")
  on.exit(close(con), add = TRUE)
  header <- character()
  repeat {
    z <- readLines(con, n = 1L)
    if (!length(z)) stop("No MatrixMarket dimension line: ", root)
    if (!startsWith(z, "%")) { header <- scan(text = z, quiet = TRUE); break }
  }
  barcodes <- readLines(file.path(root, "barcodes.tsv"))
  features <- read.delim(file.path(root, "features.tsv"), header = FALSE,
                         stringsAsFactors = FALSE)
  ids <- colnames(study)[meta$run_id == runs[i]]
  tokens <- sub("^[0-9]+_", "", ids)
  data.frame(run_id = runs[i], raw_features = header[1],
             raw_barcodes_matrix = header[2], raw_nonzero_entries = header[3],
             raw_barcode_lines = length(barcodes), raw_unique_barcodes = length(unique(barcodes)),
             raw_feature_lines = nrow(features), raw_unique_feature_ids = length(unique(features[[1]])),
             retained_cells = length(ids), retained_unique_barcodes = length(unique(tokens)),
             retained_barcodes_present_raw = sum(tokens %in% barcodes),
             retained_barcodes_absent_raw = sum(!tokens %in% barcodes),
             stringsAsFactors = FALSE)
})
inventory <- do.call(rbind, inventory)
write_tsv(inventory, "hazzard2022_raw_inventory.tsv")

# Exact record lineage through every artifact that exists locally. The central
# pv_all_studies.rds is absent and is represented explicitly rather than inferred.
final_ids <- rownames(final)
lineage <- do.call(rbind, lapply(seq_along(runs), function(i) {
  ids <- colnames(study)[meta$run_id == runs[i]]
  tok <- sub("^[0-9]+_", "", ids)
  data.frame(run_id = runs[i], source_barcode = tok,
             plavisca_barcode = ids,
             study_rds_cell = ifelse(ids %in% colnames(study), ids, NA_character_),
             pv_all_studies_cell = NA_character_,
             final_cell_id = ifelse(ids %in% final_ids, ids, NA_character_),
             mapping_status = ifelse(ids %in% final_ids, "one_to_one_study_to_deployed", "missing_deployed"),
             stringsAsFactors = FALSE)
}))
write_tsv(lineage, "hazzard2022_cell_lineage.tsv")

# Identifier integrity and barcode-token reuse across runs.
token_tab <- table(lineage$source_barcode)
collisions <- lineage[lineage$source_barcode %in% names(token_tab)[token_tab > 1L], ]
if (!nrow(collisions)) collisions <- lineage[FALSE, ]
write_tsv(collisions, "hazzard2022_barcode_token_reuse.tsv")

lineage_summary <- do.call(rbind, lapply(runs, function(run) {
  ids <- lineage$plavisca_barcode[lineage$run_id == run]
  data.frame(run_id = run, after_emptydrops = length(ids), study_rds = sum(ids %in% colnames(study)),
             pv_all_studies = NA_integer_, deployed_cleaned = sum(ids %in% final_ids),
             deployed_raw = NA_integer_, deployed_normalized = NA_integer_, deployed_scaled = NA_integer_,
             missing_study_to_cleaned = sum(!ids %in% final_ids),
             extra_cleaned_vs_study = sum(final_ids[final$run_id == run] %notin% ids),
             stringsAsFactors = FALSE)
}))
# Exact sparse-profile hashes followed by elementwise verification.
hashes <- vapply(seq_len(ncol(counts)), function(j) {
  p1 <- counts@p[j] + 1L; p2 <- counts@p[j + 1L]
  idx <- if (p2 >= p1) counts@i[p1:p2] else integer()
  val <- if (p2 >= p1) counts@x[p1:p2] else numeric()
  digest(list(idx, val), algo = "sha256", serialize = TRUE)
}, character(1))
groups <- split(seq_along(hashes), hashes)
groups <- groups[lengths(groups) > 1L]
dup_rows <- list(); k <- 0L
for (g in groups) {
  ref <- counts[, g[1], drop = FALSE]
  if (!all(vapply(g[-1], function(j) identical(ref, counts[, j, drop = FALSE]), logical(1)))) next
  k <- k + 1L
  dup_rows[[k]] <- data.frame(
    duplicate_group = k, cell_id = colnames(counts)[g], run_id = meta$run_id[g],
    source_barcode = sub("^[0-9]+_", "", colnames(counts)[g]),
    upstream_columns_separately_named = TRUE,
    upstream_status = "present_as_separate_STARsolo_barcode_columns_pending_raw_value_check",
    classification = "unresolved", stringsAsFactors = FALSE)
}
duplicates <- if (length(dup_rows)) do.call(rbind, dup_rows) else data.frame(
  duplicate_group = integer(), cell_id = character(), run_id = character(),
  source_barcode = character(), upstream_columns_separately_named = logical(),
  upstream_status = character(), classification = character())
write_tsv(duplicates, "hazzard2022_exact_duplicate_profiles.tsv")

# Metadata/stage/HPI distributions and field-level study-to-deployed concordance.
stage_cols <- intersect(c("run_id", "sample_type", "host_species", "parasite_stage",
                          "life_cycle_stage", "development_phase", "blood_stage",
                          "parasite_stages", "hour_post_invasion", "pred_gametocyte"), names(final))
stage_dist <- as.data.frame(table(final[stage_cols], useNA = "ifany"), stringsAsFactors = FALSE)
stage_dist <- stage_dist[stage_dist$Freq > 0L, , drop = FALSE]
write_tsv(stage_dist, "hazzard2022_deployed_annotation_distribution.tsv")

shared <- intersect(names(meta), names(final))
mm <- match(colnames(study), rownames(final))
metadata_concordance <- do.call(rbind, lapply(shared, function(field) {
  a <- as.character(meta[[field]]); b <- as.character(final[[field]][mm])
  eq <- (is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & a == b)
  data.frame(field = field, compared = length(a), mismatches = sum(!eq), stringsAsFactors = FALSE)
}))
write_tsv(metadata_concordance, "hazzard2022_study_deployed_metadata.tsv")

# Read deployed expression only after cheap checks. Flattened files use row names as cell IDs.
raw_df <- readRDS("../../PlaViSca/data/raw_df.rds")
norm_df <- readRDS("../../PlaViSca/data/normalize_df.rds")
scale_df <- readRDS("../../PlaViSca/data/scale_df.rds")
expr_cols <- intersect(gsub("-", "_", rownames(counts)), names(raw_df))
source_rows <- match(gsub("-", "_", expr_cols), gsub("-", "_", rownames(counts)))
raw_rows <- match(colnames(study), rownames(raw_df))
norm_rows <- match(colnames(study), rownames(norm_df))
scale_rows <- match(colnames(study), rownames(scale_df))
stopifnot(!anyNA(raw_rows), !anyNA(norm_rows), !anyNA(scale_rows))
lineage_summary$deployed_raw <- vapply(runs, function(run) sum(colnames(study)[meta$run_id == run] %in% rownames(raw_df)), integer(1))
lineage_summary$deployed_normalized <- vapply(runs, function(run) sum(colnames(study)[meta$run_id == run] %in% rownames(norm_df)), integer(1))
lineage_summary$deployed_scaled <- vapply(runs, function(run) sum(colnames(study)[meta$run_id == run] %in% rownames(scale_df)), integer(1))
write_tsv(lineage_summary, "hazzard2022_lineage_summary.tsv")

raw_mismatch <- 0; raw_max <- 0; norm_mismatch <- 0; norm_max <- 0
values <- 0
for (j in seq_along(expr_cols)) {
  src <- as.numeric(counts[source_rows[j], ])
  rv <- as.numeric(raw_df[raw_rows, expr_cols[j]])
  expected <- log1p(src / Matrix::colSums(counts) * 10000)
  nv <- as.numeric(norm_df[norm_rows, expr_cols[j]])
  d1 <- abs(src - rv); d2 <- abs(expected - nv)
  raw_mismatch <- raw_mismatch + sum(d1 != 0)
  norm_mismatch <- norm_mismatch + sum(d2 > 1e-12)
  raw_max <- max(raw_max, d1); norm_max <- max(norm_max, d2)
  values <- values + length(src)
}
expr_summary <- data.frame(
  comparison = c("study_counts_vs_deployed_raw", "formula_vs_deployed_normalized"),
  cells_compared = ncol(counts), genes_compared = length(expr_cols), values_compared = values,
  mismatch_count = c(raw_mismatch, norm_mismatch), maximum_absolute_difference = c(raw_max, norm_max),
  stringsAsFactors = FALSE)
write_tsv(expr_summary, "hazzard2022_expression_concordance.tsv")

# Cell totals are separately reported, and non-source app features are tested for nonzero values.
extra_raw <- setdiff(names(raw_df), c(expr_cols, names(final)))
extra_norm <- setdiff(names(norm_df), c(expr_cols, names(final)))
cell_totals <- data.frame(cell_id = colnames(study), run_id = meta$run_id,
                          source_total = as.numeric(Matrix::colSums(counts)),
                          deployed_raw_source_gene_total = rowSums(raw_df[raw_rows, expr_cols, drop = FALSE]),
                          source_gene_total_match = as.numeric(Matrix::colSums(counts)) == rowSums(raw_df[raw_rows, expr_cols, drop = FALSE]))
write_tsv(cell_totals, "hazzard2022_cell_totals.tsv")
extra_summary <- data.frame(artifact = c("raw_df", "normalize_df"),
                            extra_gene_columns = c(length(extra_raw), length(extra_norm)),
                            nonzero_extra_values = c(sum(as.matrix(raw_df[raw_rows, extra_raw, drop = FALSE]) != 0),
                                                     sum(as.matrix(norm_df[norm_rows, extra_norm, drop = FALSE]) != 0)))
write_tsv(extra_summary, "hazzard2022_extra_expression_features.tsv")

cat("Hazzard2022 audit inspection complete\n")
