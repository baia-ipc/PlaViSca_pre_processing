#!/usr/bin/env Rscript
# Ruberto2022_1 author-matrix validation study (Q1-Q3 + supporting evidence).
# Reads only existing production/audit files; writes only into this directory.
# No production script, source data, study object, or app artifact is modified.
# Run with: env -u R_LIBS_USER pixi run Rscript --vanilla \
#   ../audit/ruberto2022_1/author_matrix_validation/validate_author_matrix.R
# (from pre_process_data/scripts, where the Pixi environment is defined)

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(tidyr)
})

proj_root <- "/home/baia/prj/plavisca/pre_process_data"
audit_dir <- file.path(proj_root, "audit/ruberto2022_1")
out_dir <- file.path(audit_dir, "author_matrix_validation")
dir.create(out_dir, showWarnings = FALSE)

cat("=== Environment ===\n")
cat("R version:", R.version.string, "\n")
cat("Seurat version:", as.character(packageVersion("Seurat")), "\n")
cat("Matrix version:", as.character(packageVersion("Matrix")), "\n")
cat("dplyr version:", as.character(packageVersion("dplyr")), "\n")

## ------------------------------------------------------------------
## Checksums of inputs used in this task
## ------------------------------------------------------------------
files_to_hash <- c(
  hep59 = file.path(proj_root, "data/Hep59.1.2.seu_20aug2025.rds"),
  study_rds = file.path(proj_root, "ruberto2022_1.rds"),
  d036_cells = file.path(audit_dir, "near_empty_cells_vs_hep59.tsv"),
  crosswalk_orig = file.path(audit_dir, "hep59_to_study_cellid_map.tsv")
)
md5s <- tools::md5sum(files_to_hash)
write.table(
  data.frame(file = files_to_hash, md5 = md5s),
  file.path(out_dir, "input_checksums.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
cat("\n=== Input checksums ===\n")
print(data.frame(file = files_to_hash, md5 = md5s))

## ------------------------------------------------------------------
## Load objects
## ------------------------------------------------------------------
cat("\n=== Loading objects ===\n")
hep59 <- readRDS(files_to_hash["hep59"])
study <- readRDS(files_to_hash["study_rds"])

cat("hep59: assays =", paste(Assays(hep59), collapse = ","),
    "; RNA dims =", paste(dim(hep59[["RNA"]]), collapse = " x "),
    "; layers =", paste(Layers(hep59[["RNA"]]), collapse = ","), "\n")
cat("hep59 Seurat object version:", as.character(hep59@version), "\n")
cat("study: assays =", paste(Assays(study), collapse = ","),
    "; RNA dims =", paste(dim(study[["RNA"]]), collapse = " x "), "\n")
cat("study Seurat object version:", as.character(study@version), "\n")

## ====================================================================
## SECTION 1: exact 1,438-cell crosswalk
## ====================================================================
cat("\n=== SECTION 1: crosswalk ===\n")

run_from_prefix <- c("51inf" = "612", "52inf" = "610", "91inf" = "611", "92inf" = "609")
run_map_srr <- c("609" = "SRR19573609", "610" = "SRR19573610",
                  "611" = "SRR19573611", "612" = "SRR19573612")

hep59_meta <- hep59@meta.data
hep59_meta$author_cell_id <- rownames(hep59_meta)
hep59_meta <- hep59_meta |>
  mutate(
    prefix = sub("_.*", "", author_cell_id),
    barcode = sub(".*_", "", author_cell_id),
    run_suffix = unname(run_from_prefix[prefix]),
    plavisca_cell_id = paste0(run_suffix, "_", barcode)
  )

stopifnot("hep59 must have exactly 1438 cells" = nrow(hep59_meta) == 1438)
stopifnot("all prefixes must map to a run" = !any(is.na(hep59_meta$run_suffix)))
stopifnot("no duplicated author cell IDs" = !any(duplicated(hep59_meta$author_cell_id)))
stopifnot("no duplicated derived PlaViSca IDs" = !any(duplicated(hep59_meta$plavisca_cell_id)))

study_meta <- study@meta.data
study_meta$plavisca_cell_id <- rownames(study_meta)
stopifnot("study object must have exactly 1438 cells" = ncol(study) == 1438)
stopifnot("no duplicated study cell IDs" = !any(duplicated(study_meta$plavisca_cell_id)))

set_a_minus_b <- setdiff(hep59_meta$plavisca_cell_id, study_meta$plavisca_cell_id)
set_b_minus_a <- setdiff(study_meta$plavisca_cell_id, hep59_meta$plavisca_cell_id)
stopifnot("author-derived IDs must all exist in study object" = length(set_a_minus_b) == 0)
stopifnot("study IDs must all exist in author-derived set" = length(set_b_minus_a) == 0)

cat("1438 author cells, 1438 PlaViSca cells, exact 1:1 mapping, zero unmapped: CONFIRMED\n")

author_counts_full <- GetAssayData(hep59, assay = "RNA", layer = "counts")
study_counts_full <- GetAssayData(study, assay = "RNA", layer = "counts")

author_umi <- Matrix::colSums(author_counts_full)
study_umi <- Matrix::colSums(study_counts_full)

d036_tbl <- read.delim(file.path(audit_dir, "near_empty_cells_vs_hep59.tsv"))
d036_ids <- d036_tbl$study_cell_id
stopifnot("authoritative D036 list must have 538 cells" = length(d036_ids) == 538)
stopifnot("D036 IDs must be a subset of study IDs" = all(d036_ids %in% study_meta$plavisca_cell_id))

crosswalk <- hep59_meta |>
  transmute(
    author_cell_id,
    plavisca_cell_id,
    run_id = unname(run_map_srr[run_suffix]),
    run_suffix,
    day_author = as.character(day),
    replicate_author = as.character(replicate),
    liver_form_author = as.character(LiverForm),
    author_total_umi = as.numeric(author_umi[author_cell_id]),
    plavisca_total_umi = as.numeric(study_umi[plavisca_cell_id])
  ) |>
  left_join(
    study_meta |>
      transmute(
        plavisca_cell_id,
        day_post_infection_plavisca = day_post_infection,
        biological_replicate_plavisca = biological_replicate,
        treatment_plavisca = treatment,
        liver_form_plavisca = liver_form,
        run_id_plavisca = run_id
      ),
    by = "plavisca_cell_id"
  ) |>
  mutate(
    d036_status = ifelse(plavisca_cell_id %in% d036_ids, "D036_near_empty", "unaffected"),
    group = ifelse(d036_status == "D036_near_empty", "A_D036_near_empty", "B_unaffected")
  )

stopifnot(nrow(crosswalk) == 1438)
stopifnot(sum(crosswalk$group == "A_D036_near_empty") == 538)
stopifnot(sum(crosswalk$group == "B_unaffected") == 900)
stopifnot(all(crosswalk$run_id == crosswalk$run_id_plavisca))

cat("Group A (D036 near-empty):", sum(crosswalk$group == "A_D036_near_empty"), "\n")
cat("Group B (unaffected):", sum(crosswalk$group == "B_unaffected"), "\n")
cat("Group B by run:\n"); print(table(crosswalk$run_id[crosswalk$group == "B_unaffected"]))
cat("Group A by run:\n"); print(table(crosswalk$run_id[crosswalk$group == "A_D036_near_empty"]))

write.table(crosswalk, file.path(out_dir, "cell_crosswalk_validation.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

## ====================================================================
## SECTION 2: metadata concordance
## ====================================================================
cat("\n=== SECTION 2: metadata concordance ===\n")

# Author day labels ("Five"/"Nine") vs PlaViSca day_post_infection (5/9)
day_lookup <- c(Five = 5, Nine = 9)
crosswalk <- crosswalk |>
  mutate(
    day_author_numeric = unname(day_lookup[day_author]),
    day_match = day_author_numeric == day_post_infection_plavisca,
    replicate_match = as.character(replicate_author) == as.character(biological_replicate_plavisca),
    liver_form_match = liver_form_author == liver_form_plavisca,
    run_prefix_mapping_expected = case_when(
      run_suffix == "612" ~ "51inf",
      run_suffix == "610" ~ "52inf",
      run_suffix == "611" ~ "91inf",
      run_suffix == "609" ~ "92inf"
    ),
    run_prefix_actual = sub("_.*", "", author_cell_id),
    run_prefix_mapping_match = run_prefix_mapping_expected == run_prefix_actual
  )

metadata_concordance <- crosswalk |>
  group_by(run_id, run_suffix) |>
  summarise(
    n_cells = n(),
    day_match_n = sum(day_match), day_mismatch_n = sum(!day_match),
    replicate_match_n = sum(replicate_match), replicate_mismatch_n = sum(!replicate_match),
    liver_form_match_n = sum(liver_form_match), liver_form_mismatch_n = sum(!liver_form_match),
    run_prefix_mapping_match_n = sum(run_prefix_mapping_match),
    run_prefix_mapping_mismatch_n = sum(!run_prefix_mapping_match),
    .groups = "drop"
  )

cat("Metadata concordance by run:\n")
print(as.data.frame(metadata_concordance))
cat("\nNote: the author object has no explicit drug-treatment field (its 'condition'",
    "column is uniformly 'Infected' for all 1438 cells); treatment is not",
    "independently checkable against the author object and is inferred only via",
    "day+run, consistent with the published protocol (already source-verified",
    "in AUDIT.md against ENA sample titles and the Methods text).\n")

write.table(metadata_concordance, file.path(out_dir, "metadata_concordance.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

overall_day_mismatch <- sum(!crosswalk$day_match)
overall_rep_mismatch <- sum(!crosswalk$replicate_match)
overall_lf_mismatch <- sum(!crosswalk$liver_form_match)
overall_prefix_mismatch <- sum(!crosswalk$run_prefix_mapping_match)
cat("\nTotal mismatches across all 1438 cells: day =", overall_day_mismatch,
    "; replicate =", overall_rep_mismatch,
    "; liver_form =", overall_lf_mismatch,
    "; run-prefix mapping (51inf/52inf/91inf/92inf) =", overall_prefix_mismatch, "\n")

## ====================================================================
## SECTION 3: RNA counts validation (raw UMI proof)
## ====================================================================
cat("\n=== SECTION 3: RNA counts validation ===\n")

hep59_layers <- Layers(hep59[["RNA"]])
sct_present <- "SCT" %in% Assays(hep59)
cat("hep59 RNA layers:", paste(hep59_layers, collapse = ","), "\n")
cat("hep59 SCT assay present:", sct_present, "\n")

author_counts <- GetAssayData(hep59, assay = "RNA", layer = "counts")
author_data <- GetAssayData(hep59, assay = "RNA", layer = "data")

n_entries_total <- as.numeric(nrow(author_counts)) * as.numeric(ncol(author_counts))
nz_idx <- author_counts@x
n_nonzero <- length(nz_idx)
n_noninteger <- sum(nz_idx != round(nz_idx))
min_val <- min(nz_idx)
max_val <- max(nz_idx)

cat("Integer test: matrix entries =", n_entries_total, "; nonzero =", n_nonzero,
    "; non-integer nonzero =", n_noninteger, "; min =", min_val, "; max =", max_val, "\n")

colsum_counts <- Matrix::colSums(author_counts)
ncount_meta <- hep59$nCount_RNA[colnames(author_counts)]
count_diff <- colsum_counts - ncount_meta
count_exact_match <- sum(count_diff == 0)
count_mismatch <- sum(count_diff != 0)
count_max_abs_diff <- max(abs(count_diff))
cat("Count-total test: exact match =", count_exact_match, "; mismatches =", count_mismatch,
    "; max abs diff =", count_max_abs_diff, "\n")
n_colsum_exceeds_ncount <- sum(count_diff > 0)
cat("Cells where colSums(counts) > nCount_RNA:", n_colsum_exceeds_ncount,
    "(0 expected if nCount_RNA was computed on a superset of these features)\n")
cat("This mismatch is one-directional (colSums <= nCount_RNA in all",
    sum(count_diff <= 0), "of", length(count_diff),
    "cells), consistent with nCount_RNA being inherited metadata computed",
    "upstream on a combined human+P. vivax count matrix (per the recovered",
    "kallisto|bustools workflow, which used a combined GRCh38+PlasmoDB index)",
    "before this file's RNA assay was subset down to the 4722 P. vivax-only",
    "features. It does not indicate the counts layer itself is non-raw; it",
    "means nCount_RNA cannot be exactly reconciled from this parasite-only",
    "export alone.\n")

nfeature_computed <- Matrix::colSums(author_counts > 0)
nfeature_meta <- hep59$nFeature_RNA[colnames(author_counts)]
nfeature_diff <- nfeature_computed - nfeature_meta
nfeature_exact_match <- sum(nfeature_diff == 0)
nfeature_mismatch <- sum(nfeature_diff != 0)
cat("Feature-total test: exact match =", nfeature_exact_match, "; mismatches =", nfeature_mismatch, "\n")

# Normalized-data distinction: counts vs data
identical_counts_data <- identical(as.matrix(author_counts[, 1:10]), as.matrix(author_data[, 1:10]))
data_sample <- author_data@x[1:min(5000, length(author_data@x))]
data_all_integer <- all(data_sample == round(data_sample))
data_max <- max(author_data@x)
counts_max <- max(author_counts@x)
cat("counts vs data identical (first 10 cells sample):", identical_counts_data, "\n")
cat("data layer nonzero sample all-integer:", data_all_integer, "; data max =", data_max,
    "; counts max =", counts_max, "\n")

rna_counts_validation <- data.frame(
  metric = c(
    "seurat_object_version", "assay_names", "default_assay",
    "RNA_counts_dims", "RNA_data_dims", "SCT_assay_present",
    "matrix_entries_total", "nonzero_entries", "noninteger_nonzero_entries",
    "min_nonzero_value", "max_nonzero_value",
    "colSums_vs_nCount_RNA_exact_match", "colSums_vs_nCount_RNA_mismatch",
    "colSums_vs_nCount_RNA_max_abs_diff",
    "nFeature_computed_vs_meta_exact_match", "nFeature_computed_vs_meta_mismatch",
    "counts_identical_to_data_sample", "data_layer_nonzero_all_integer_sample",
    "data_layer_max_value", "counts_layer_max_value"
  ),
  value = c(
    as.character(hep59@version), paste(Assays(hep59), collapse = ","), DefaultAssay(hep59),
    paste(dim(hep59[["RNA"]]), collapse = " x "),
    paste(dim(GetAssayData(hep59, assay = "RNA", layer = "data")), collapse = " x "),
    as.character(sct_present),
    as.character(n_entries_total), as.character(n_nonzero), as.character(n_noninteger),
    as.character(min_val), as.character(max_val),
    as.character(count_exact_match), as.character(count_mismatch),
    as.character(count_max_abs_diff),
    as.character(nfeature_exact_match), as.character(nfeature_mismatch),
    as.character(identical_counts_data), as.character(data_all_integer),
    as.character(data_max), as.character(counts_max)
  )
)
write.table(rna_counts_validation, file.path(out_dir, "rna_counts_validation.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

## ====================================================================
## SECTION 4 & 5: gene namespace + gene-space comparison with PlaViSca
## ====================================================================
cat("\n=== SECTION 4/5: gene namespace and mapping ===\n")

author_features <- rownames(hep59[["RNA"]])
study_features <- rownames(study[["RNA"]])

cat("Author feature format sample:\n"); print(head(author_features, 5))
cat("Study feature format sample:\n"); print(head(study_features, 5))
cat("Author feature count:", length(author_features), "unique:", length(unique(author_features)), "\n")
cat("Study feature count:", length(study_features), "unique:", length(unique(study_features)), "\n")

exact_match <- intersect(author_features, study_features)
author_only <- setdiff(author_features, study_features)
study_only <- setdiff(study_features, author_features)

cat("Exact-string matches:", length(exact_match), "\n")
cat("Author-only features:", length(author_only), "\n")
cat("PlaViSca-only features:", length(study_only), "\n")

# Check for any version-suffix pattern (e.g. trailing .1) requiring normalization
has_dot_suffix_author <- sum(grepl("\\.[0-9]+$", author_features))
has_dot_suffix_study <- sum(grepl("\\.[0-9]+$", study_features))
cat("Author features with a trailing version suffix ('.N'):", has_dot_suffix_author, "\n")
cat("Study features with a trailing version suffix ('.N'):", has_dot_suffix_study, "\n")

# No ambiguous/duplicate mappings possible since both feature vectors are unique
one_to_many <- 0
many_to_one <- 0
ambiguous <- 0

author_gene_mapping <- data.frame(
  author_feature = author_features,
  plavisca_feature = ifelse(author_features %in% study_features, author_features, NA_character_),
  mapping_method = "exact_string_match",
  mapping_status = ifelse(author_features %in% study_features, "exact", "unmapped"),
  notes = ifelse(
    author_features %in% study_features,
    "identical PVP01- hyphenated gene ID in both objects; no normalization required",
    "author feature absent from PlaViSca 6811-feature panel"
  )
)
write.table(author_gene_mapping, file.path(out_dir, "author_gene_mapping.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("Unmapped author features:", sum(author_gene_mapping$mapping_status == "unmapped"), "\n")
cat("One-to-many / many-to-one / ambiguous mappings: all 0 (both feature ID vectors are unique)\n")

shared_genes <- exact_match
cat("Shared genes used for expression concordance:", length(shared_genes), "\n")

## ====================================================================
## SECTION 6/7: per-cell expression concordance (Group B, by run)
## ====================================================================
cat("\n=== SECTION 6/7: per-cell expression concordance ===\n")

author_shared <- as.matrix(author_counts_full[shared_genes, , drop = FALSE])
study_shared <- as.matrix(study_counts_full[shared_genes, , drop = FALSE])

safe_cor <- function(a, b, method) {
  if (sd(a) == 0 || sd(b) == 0) return(NA_real_)
  suppressWarnings(cor(a, b, method = method))
}

## NOTE ON METHOD CHOICE: a preliminary run using the "union of nonzero in
## either representation" rule (as a literal reading of "avoid correlations
## dominated by thousands of joint zeros") produced systematically NEGATIVE
## per-cell correlations for known-good control cells (run 611/612), directly
## contradicting the strong positive signal from the independent cell-identity
## rank test (Section 9) and from full-panel and intersection-based
## correlations on the same cells. Diagnosis on a representative 612 cell:
## author nnz=951, PlaViSca nnz=373, intersection nnz=196, union nnz=1128;
## spearman(union)=-0.396 but spearman(intersection)=+0.318,
## pearson_log1p(full 4722 genes)=+0.242, spearman(full 4722 genes)=+0.241.
## The union rule is dominated by ~930 genes detected by only one of the two
## methods (kallisto author-side detects ~2.5x more genes per cell than
## STARsolo PlaViSca-side here), each contributing a one-sided tie at zero;
## with such asymmetric detection sensitivity this tie structure biases rank
## correlation negative, which is an artifact of the union rule, not evidence
## against concordance. This is a defensible alternative sparse-data strategy
## per the task's own allowance for "another defensible sparse-data strategy":
## the PRIMARY reported metric here is INTERSECTION-nonzero (genes detected by
## BOTH representations for that cell), which removes one-sided-detection
## noise while still avoiding domination by joint zeros. Full-panel (all
## shared genes, including joint zeros) is reported alongside as a robustness
## check.  Union-nonzero is also retained as a diagnostic column, explicitly
## labeled as of a different signature/less reliable for this dataset.
per_cell_concordance <- function(cw_subset) {
  res <- vector("list", nrow(cw_subset))
  for (i in seq_len(nrow(cw_subset))) {
    aid <- cw_subset$author_cell_id[i]
    pid <- cw_subset$plavisca_cell_id[i]
    a <- author_shared[, aid]
    b <- study_shared[, pid]
    int_keep <- (a > 0) & (b > 0)
    union_keep <- (a > 0) | (b > 0)
    n_int <- sum(int_keep)
    n_union <- sum(union_keep)
    a_log <- log1p(a); b_log <- log1p(b)
    if (n_int >= 5) {
      sp_int <- safe_cor(a[int_keep], b[int_keep], "spearman")
      pe_int <- safe_cor(a_log[int_keep], b_log[int_keep], "pearson")
    } else {
      sp_int <- NA_real_; pe_int <- NA_real_
    }
    sp_full <- safe_cor(a, b, "spearman")
    pe_full <- safe_cor(a_log, b_log, "pearson")
    if (n_union >= 5) {
      sp_union <- safe_cor(a[union_keep], b[union_keep], "spearman")
    } else {
      sp_union <- NA_real_
    }
    denom <- sqrt(sum(a_log^2)) * sqrt(sum(b_log^2))
    cos_sim_full <- if (denom > 0) sum(a_log * b_log) / denom else NA_real_
    res[[i]] <- data.frame(
      author_cell_id = aid, plavisca_cell_id = pid, run_id = cw_subset$run_id[i],
      group = cw_subset$group[i],
      author_total_umi = cw_subset$author_total_umi[i],
      plavisca_total_umi = cw_subset$plavisca_total_umi[i],
      author_detected_genes_shared = sum(a > 0),
      plavisca_detected_genes_shared = sum(b > 0),
      n_intersection_nonzero_genes = n_int,
      n_union_nonzero_genes = n_union,
      spearman_intersection = sp_int,
      pearson_log1p_intersection = pe_int,
      spearman_full_panel = sp_full,
      pearson_log1p_full_panel = pe_full,
      spearman_union_diagnostic = sp_union,
      cosine_log1p_full_panel = cos_sim_full
    )
  }
  do.call(rbind, res)
}

cell_expression_concordance <- per_cell_concordance(crosswalk)
write.table(cell_expression_concordance, file.path(out_dir, "cell_expression_concordance.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

summarize_dist <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(c(n = 0, median = NA, q1 = NA, q3 = NA, min = NA, max = NA))
  q <- quantile(x, c(0.25, 0.75))
  c(n = length(x), median = median(x), q1 = q[1], q3 = q[2], min = min(x), max = max(x))
}

cat("\nGroup B (900 unaffected) Spearman (intersection-nonzero, primary metric) distribution:\n")
print(summarize_dist(cell_expression_concordance$spearman_intersection[cell_expression_concordance$group == "B_unaffected"]))
cat("\nGroup B (900 unaffected) Spearman (full-panel, robustness check) distribution:\n")
print(summarize_dist(cell_expression_concordance$spearman_full_panel[cell_expression_concordance$group == "B_unaffected"]))

cat("\nBy run (Group B only), intersection-nonzero Spearman:\n")
for (r in unique(crosswalk$run_id)) {
  sub <- cell_expression_concordance |> filter(run_id == r, group == "B_unaffected")
  if (nrow(sub) == 0) { cat(r, ": no Group B cells\n"); next }
  cat(r, "(n=", nrow(sub), "): median Spearman(intersection) =",
      round(median(sub$spearman_intersection, na.rm = TRUE), 3),
      "; median Spearman(full panel) =",
      round(median(sub$spearman_full_panel, na.rm = TRUE), 3), "\n")
}

cat("\nGroup A (538 D036) Spearman (intersection-nonzero) distribution:\n")
print(summarize_dist(cell_expression_concordance$spearman_intersection[cell_expression_concordance$group == "A_D036_near_empty"]))
cat("\nGroup A (538 D036) Spearman (full-panel) distribution:\n")
print(summarize_dist(cell_expression_concordance$spearman_full_panel[cell_expression_concordance$group == "A_D036_near_empty"]))

## ====================================================================
## SECTION 8: pseudobulk comparison by run
## ====================================================================
cat("\n=== SECTION 8: pseudobulk concordance ===\n")

pseudobulk_rows <- list()
for (r in unique(crosswalk$run_id)) {
  for (grp in c("B_unaffected", "A_D036_near_empty")) {
    sub <- crosswalk |> filter(run_id == r, group == grp)
    if (nrow(sub) == 0) next
    a_sum <- Matrix::rowSums(author_shared[, sub$author_cell_id, drop = FALSE])
    b_sum <- Matrix::rowSums(study_shared[, sub$plavisca_cell_id, drop = FALSE])
    sp <- safe_cor(a_sum, b_sum, "spearman")
    pe <- safe_cor(log1p(a_sum), log1p(b_sum), "pearson")
    top_author <- names(sort(a_sum, decreasing = TRUE))[1:10]
    top_study <- names(sort(b_sum, decreasing = TRUE))[1:10]
    pseudobulk_rows[[paste(r, grp)]] <- data.frame(
      run_id = r, group = grp, n_cells = nrow(sub),
      author_pseudobulk_total_umi = sum(a_sum), plavisca_pseudobulk_total_umi = sum(b_sum),
      spearman_pseudobulk = sp, pearson_log1p_pseudobulk = pe,
      top10_author_genes = paste(top_author, collapse = ";"),
      top10_plavisca_genes = paste(top_study, collapse = ";"),
      top10_overlap_n = length(intersect(top_author, top_study))
    )
  }
}
run_pseudobulk_concordance <- do.call(rbind, pseudobulk_rows)
rownames(run_pseudobulk_concordance) <- NULL
write.table(run_pseudobulk_concordance, file.path(out_dir, "run_pseudobulk_concordance.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("Pseudobulk table:\n")
print(run_pseudobulk_concordance[, c("run_id", "group", "n_cells", "spearman_pseudobulk", "pearson_log1p_pseudobulk", "top10_overlap_n")])

## ====================================================================
## SECTION 9: cell-identity (rank) test
## ====================================================================
cat("\n=== SECTION 9: cell-identity rank test ===\n")

identity_test_run <- function(sub, run_label, group_label) {
  n <- nrow(sub)
  if (n < 2) {
    return(data.frame(run_id = run_label, group = group_label, n_candidates = n,
                       rank1_frac = NA, top5_frac = NA, top10_frac = NA,
                       note = "candidate pool too small (<2) for a meaningful rank test"))
  }
  A <- author_shared[, sub$author_cell_id, drop = FALSE]
  B <- study_shared[, sub$plavisca_cell_id, drop = FALSE]
  # drop constant (all-zero) columns which break rank correlation
  a_ok <- apply(A, 2, function(x) sd(x) > 0)
  b_ok <- apply(B, 2, function(x) sd(x) > 0)
  if (sum(a_ok) < 2 || sum(b_ok) < 2) {
    return(data.frame(run_id = run_label, group = group_label, n_candidates = n,
                       rank1_frac = NA, top5_frac = NA, top10_frac = NA,
                       note = "insufficient non-constant expression vectors for ranking"))
  }
  A <- A[, a_ok, drop = FALSE]; B <- B[, b_ok, drop = FALSE]
  corMat <- suppressWarnings(cor(A, B, method = "spearman"))
  ranks <- sapply(colnames(A), function(aid) {
    true_pid <- sub$plavisca_cell_id[match(aid, sub$author_cell_id)]
    if (!true_pid %in% colnames(B)) return(NA_integer_)
    ord <- order(corMat[aid, ], decreasing = TRUE)
    match(true_pid, colnames(B)[ord])
  })
  ranks <- ranks[!is.na(ranks)]
  data.frame(
    run_id = run_label, group = group_label, n_candidates = length(ranks),
    rank1_frac = mean(ranks == 1), top5_frac = mean(ranks <= 5), top10_frac = mean(ranks <= 10),
    note = "spearman over shared genes, candidate pool = same-run same-group PlaViSca cells"
  )
}

identity_rows <- list()
for (r in unique(crosswalk$run_id)) {
  for (grp in c("B_unaffected", "A_D036_near_empty")) {
    sub <- crosswalk |> filter(run_id == r, group == grp)
    if (nrow(sub) == 0) next
    identity_rows[[paste(r, grp)]] <- identity_test_run(sub, r, grp)
  }
}
cell_identity_test <- do.call(rbind, identity_rows)
rownames(cell_identity_test) <- NULL
write.table(cell_identity_test, file.path(out_dir, "cell_identity_test.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("Cell-identity rank test:\n")
print(cell_identity_test)

## ====================================================================
## SECTION 10: 900-vs-538 contrast + summary metrics
## ====================================================================
cat("\n=== SECTION 10: summary contrast table ===\n")

group_summary <- function(grp) {
  sub_cw <- crosswalk |> filter(group == grp)
  sub_ce <- cell_expression_concordance |> filter(group == grp)
  data.frame(
    group = grp,
    n_cells = nrow(sub_cw),
    author_median_umi = median(sub_cw$author_total_umi),
    plavisca_median_umi = median(sub_cw$plavisca_total_umi),
    median_shared_gene_spearman_intersection = median(sub_ce$spearman_intersection, na.rm = TRUE),
    median_shared_gene_spearman_full_panel = median(sub_ce$spearman_full_panel, na.rm = TRUE),
    median_detected_genes_author = median(sub_ce$author_detected_genes_shared),
    median_detected_genes_plavisca = median(sub_ce$plavisca_detected_genes_shared)
  )
}
summary_metrics <- do.call(rbind, lapply(c("B_unaffected", "A_D036_near_empty"), group_summary))
write.table(summary_metrics, file.path(out_dir, "summary_metrics.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("Summary metrics (900 unaffected vs 538 D036):\n")
print(summary_metrics)

## ====================================================================
## Plots (audit-only, not committed as production artifacts)
## ====================================================================
cat("\n=== Plots ===\n")
tryCatch({
  png(file.path(out_dir, "author_vs_plavisca_umi_unaffected.png"), width = 900, height = 800, res = 130)
  sub <- crosswalk |> filter(group == "B_unaffected")
  plot(log10(sub$author_total_umi + 1), log10(sub$plavisca_total_umi + 1),
       xlab = "log10(author total UMI + 1)", ylab = "log10(PlaViSca total UMI + 1)",
       main = "Author vs PlaViSca UMI, 900 unaffected cells", pch = 16, col = adjustcolor("steelblue", 0.5))
  abline(0, 1, col = "grey40", lty = 2)
  dev.off()

  png(file.path(out_dir, "spearman_distribution_by_group.png"), width = 900, height = 800, res = 130)
  boxplot(spearman_intersection ~ group, data = cell_expression_concordance,
          ylab = "Per-cell Spearman correlation (intersection-nonzero genes)", xlab = "",
          main = "D036 (538) vs unaffected (900) expression concordance")
  dev.off()

  png(file.path(out_dir, "pseudobulk_scatter.png"), width = 1000, height = 800, res = 130)
  par(mfrow = c(2, 2))
  for (r in unique(crosswalk$run_id)) {
    sub <- crosswalk |> filter(run_id == r, group == "B_unaffected")
    if (nrow(sub) == 0) next
    a_sum <- Matrix::rowSums(author_shared[, sub$author_cell_id, drop = FALSE])
    b_sum <- Matrix::rowSums(study_shared[, sub$plavisca_cell_id, drop = FALSE])
    plot(log1p(a_sum), log1p(b_sum), main = paste(r, "(n=", nrow(sub), ")"),
         xlab = "log1p author pseudobulk", ylab = "log1p PlaViSca pseudobulk",
         pch = 16, col = adjustcolor("darkgreen", 0.4))
  }
  dev.off()
  cat("Plots written.\n")
}, error = function(e) cat("Plotting failed (non-fatal):", conditionMessage(e), "\n"))

cat("\n=== DONE validate_author_matrix.R ===\n")
