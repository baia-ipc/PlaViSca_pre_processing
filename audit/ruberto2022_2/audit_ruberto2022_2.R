#!/usr/bin/env Rscript
# Ruberto2022_2 full forensic audit — read-only inspection script.
# Writes TSV evidence into pre_process_data/audit/ruberto2022_2/. Does not
# write to any production input, script, or RDS artifact.
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(tidyverse)
  library(digest)
})

audit_dir <- "../audit/ruberto2022_2"
dir.create(audit_dir, showWarnings = FALSE, recursive = TRUE)

log <- function(...) cat(sprintf(...), "\n")

## ---------------------------------------------------------------------
## A. Source object provenance
## ---------------------------------------------------------------------
log("=== A. Source object provenance ===")
src <- readRDS("../data/6_PvSPZ.BS.combined.rds")

src_meta <- src@meta.data
write.table(
  data.frame(field = colnames(src_meta)),
  file.path(audit_dir, "source_object_metadata_fields.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)

stage_tab <- as.data.frame(table(src$Stage, useNA = "always"))
colnames(stage_tab) <- c("Stage", "n_cells")
write.table(stage_tab, file.path(audit_dir, "source_stage_counts.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE)

sample_tab <- as.data.frame(table(src$Sample, src$Stage, useNA = "always"))
colnames(sample_tab) <- c("Sample", "Stage", "n_cells")
write.table(sample_tab, file.path(audit_dir, "source_sample_stage_counts.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE)

assay_info <- do.call(rbind, lapply(Assays(src), function(a) {
  data.frame(
    assay = a,
    n_features = nrow(src[[a]]),
    counts_range_min = suppressWarnings(min(GetAssayData(src, assay = a, layer = "counts")@x, na.rm = TRUE)),
    counts_range_max = suppressWarnings(max(GetAssayData(src, assay = a, layer = "counts")@x, na.rm = TRUE)),
    counts_integer_like = all(GetAssayData(src, assay = a, layer = "counts")@x ==
      round(GetAssayData(src, assay = a, layer = "counts")@x))
  )
}))
write.table(assay_info, file.path(audit_dir, "source_object_assay_summary.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE)

writeLines(
  c(
    paste("source_object_class:", class(src)),
    paste("source_object_seurat_version:", as.character(src@version)),
    paste("n_cells_total:", ncol(src)),
    paste("default_assay:", DefaultAssay(src)),
    paste("commands:", paste(names(src@commands), collapse = ";")),
    paste("n_sporozoite:", sum(src$Stage == "Sporozoite", na.rm = TRUE)),
    paste("n_bloodstage:", sum(src$Stage == "BloodStage", na.rm = TRUE))
  ),
  file.path(audit_dir, "source_object_provenance_summary.txt")
)

## ---------------------------------------------------------------------
## B. Lineage reconstruction
## ---------------------------------------------------------------------
log("=== B. Lineage reconstruction ===")

spz_cell_ids <- WhichCells(src, expression = Stage == "Sporozoite")

lineage <- data.frame(source_cell_id = spz_cell_ids) %>%
  mutate(
    source_barcode = sub("_.*", "", source_cell_id),
    source_sample = sub("^[^_]+_", "", source_cell_id),
    run_id = case_when(
      source_sample == "Case_909" ~ "ERR5087438",
      source_sample == "Case_922" ~ "ERR5087439",
      source_sample == "Case_923" ~ "ERR5087440",
      TRUE ~ NA_character_
    ),
    run_suffix = case_when(
      source_sample == "Case_909" ~ "438",
      source_sample == "Case_922" ~ "439",
      source_sample == "Case_923" ~ "440",
      TRUE ~ NA_character_
    ),
    plavisca_barcode = paste0(run_suffix, "_", source_barcode)
  )

r2 <- readRDS("../ruberto2022_2.rds")
r2_ids <- colnames(r2)

cd <- readRDS("../../PlaViSca/data/cleaned_dataset.rds")
mr <- cd$mr_data
deployed_ids <- rownames(mr)[mr$study_label == "Ruberto2022_2"]

lineage <- lineage %>%
  mutate(
    final_cell_id = plavisca_barcode,
    in_ruberto2022_2_rds = final_cell_id %in% r2_ids,
    in_deployed_cleaned_dataset = final_cell_id %in% deployed_ids,
    mapping_status = case_when(
      is.na(run_id) ~ "unmapped_source_sample",
      !in_ruberto2022_2_rds ~ "missing_from_study_rds",
      !in_deployed_cleaned_dataset ~ "missing_from_deployed",
      TRUE ~ "matched"
    )
  )

dup_final <- lineage %>%
  count(final_cell_id) %>%
  filter(n > 1)

write.table(lineage, file.path(audit_dir, "lineage_table.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE)

extra_in_r2 <- setdiff(r2_ids, lineage$final_cell_id)
extra_in_deployed <- setdiff(deployed_ids, lineage$final_cell_id)

writeLines(
  c(
    paste("n_source_sporozoite_cells:", nrow(lineage)),
    paste("n_unique_source_cell_ids:", length(unique(lineage$source_cell_id))),
    paste("n_unique_final_cell_ids_in_lineage:", length(unique(lineage$final_cell_id))),
    paste("n_matched:", sum(lineage$mapping_status == "matched")),
    paste("n_missing_from_study_rds:", sum(lineage$mapping_status == "missing_from_study_rds")),
    paste("n_missing_from_deployed:", sum(lineage$mapping_status == "missing_from_deployed")),
    paste("n_unmapped_source_sample:", sum(lineage$mapping_status == "unmapped_source_sample")),
    paste("n_duplicate_final_ids_in_lineage:", nrow(dup_final)),
    paste("n_ruberto2022_2_rds_cells:", length(r2_ids)),
    paste("n_extra_cells_in_r2_rds_not_in_lineage:", length(extra_in_r2)),
    paste("n_deployed_ruberto2022_2_cells:", length(deployed_ids)),
    paste("n_extra_cells_in_deployed_not_in_lineage:", length(extra_in_deployed)),
    paste("r2_rds_ids_equal_deployed_ids_asset:", setequal(r2_ids, deployed_ids))
  ),
  file.path(audit_dir, "lineage_summary.txt")
)

per_run <- lineage %>%
  group_by(run_id) %>%
  summarise(
    n_source = n(),
    n_matched = sum(mapping_status == "matched"),
    .groups = "drop"
  )
r2_per_run <- as.data.frame(table(r2$run_id))
colnames(r2_per_run) <- c("run_id", "n_ruberto2022_2_rds")
deployed_per_run <- as.data.frame(table(mr$run_id[mr$study_label == "Ruberto2022_2"]))
colnames(deployed_per_run) <- c("run_id", "n_deployed")

per_run_full <- per_run %>%
  left_join(r2_per_run, by = "run_id") %>%
  left_join(deployed_per_run, by = "run_id")
write.table(per_run_full, file.path(audit_dir, "per_run_cell_counts.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE)

log("Lineage reconstruction complete.")

## ---------------------------------------------------------------------
## C. STARsolo raw-matrix membership and raw expression reconstruction
## ---------------------------------------------------------------------
log("=== C. STARsolo raw matrix reconstruction ===")

# base-R GFF3 parsing (avoids the broken rtracklayer/GenomeInfoDbData chain
# in this Pixi environment); extracts only the rRNA ID attribute, matching
# the production script's rRNA_pv$ID construction exactly.
gff_lines <- readLines("ref/PlasmoDB-68_PvivaxP01.gff")
gff_lines <- gff_lines[!startsWith(gff_lines, "#")]
gff_split <- strsplit(gff_lines, "\t", fixed = TRUE)
gff_split <- gff_split[lengths(gff_split) >= 9]
gff_type <- vapply(gff_split, `[`, character(1), 3)
gff_attr <- vapply(gff_split, `[`, character(1), 9)
rRNA_attr <- gff_attr[gff_type == "rRNA"]
rRNA_id <- sub(";.*", "", sub("^ID=", "", rRNA_attr))
rRNA_pv <- data.frame(ID = gsub("\\.1", "", rRNA_id))

runs <- c("ERR5087438", "ERR5087439", "ERR5087440")
recon_list <- list()
membership_rows <- list()

for (run in runs) {
  path <- file.path("../counts/35926062", paste0(run, "_solo_out/Solo.out/GeneFull/raw"))
  bc <- readLines(file.path(path, "barcodes.tsv"))
  feat <- read.delim(file.path(path, "features.tsv"), header = FALSE, stringsAsFactors = FALSE)
  mm <- readMM(file.path(path, "matrix.mtx"))
  rownames(mm) <- feat$V1
  colnames(mm) <- bc

  # membership check: which selected barcodes for this run exist in the raw whitelist
  sel <- lineage %>% filter(run_id == run)
  in_whitelist <- sel$source_barcode %in% bc
  membership_rows[[run]] <- data.frame(
    run_id = run,
    n_selected = nrow(sel),
    n_in_raw_whitelist = sum(in_whitelist),
    n_missing_from_raw_whitelist = sum(!in_whitelist)
  )

  # subset to selected barcodes present, remove rRNA
  keep_bc <- sel$source_barcode[in_whitelist]
  sub_mm <- mm[, keep_bc, drop = FALSE]
  sub_mm <- sub_mm[!(rownames(sub_mm) %in% rRNA_pv$ID), , drop = FALSE]
  colnames(sub_mm) <- paste0(sel$run_suffix[in_whitelist], "_", keep_bc)
  recon_list[[run]] <- sub_mm
  rm(mm)
  gc()
}

membership_tab <- do.call(rbind, membership_rows)
write.table(membership_tab, file.path(audit_dir, "starsolo_raw_membership.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE)

# combine reconstructed raw counts across runs (union of gene rows)
all_genes <- Reduce(union, lapply(recon_list, rownames))
recon_full <- do.call(cbind, lapply(recon_list, function(m) {
  missing_genes <- setdiff(all_genes, rownames(m))
  if (length(missing_genes) > 0) {
    extra <- Matrix(0, nrow = length(missing_genes), ncol = ncol(m), sparse = TRUE)
    rownames(extra) <- missing_genes
    m <- rbind(m, extra)
  }
  m[all_genes, , drop = FALSE]
}))

log("Reconstructed STARsolo raw (rRNA-removed) matrix: %d genes x %d cells", nrow(recon_full), ncol(recon_full))

## ---------------------------------------------------------------------
## D. Expression concordance: STARsolo reconstruction vs ruberto2022_2.rds
## ---------------------------------------------------------------------
log("=== D. Expression concordance: STARsolo -> ruberto2022_2.rds ===")

r2_counts <- GetAssayData(r2, assay = "RNA", layer = "counts")
rownames(r2_counts) <- gsub("-", "_", rownames(r2_counts))

common_cells_1 <- intersect(colnames(recon_full), colnames(r2_counts))
common_genes_1 <- intersect(rownames(recon_full), rownames(r2_counts))

log("Common cells (STARsolo recon vs study RDS): %d / %d", length(common_cells_1), ncol(r2_counts))
log("Common genes (STARsolo recon vs study RDS): %d / %d", length(common_genes_1), nrow(r2_counts))

a1 <- recon_full[common_genes_1, common_cells_1]
b1 <- r2_counts[common_genes_1, common_cells_1]
diff1 <- a1 - b1
n_mismatch_1 <- sum(diff1@x != 0)
max_abs_diff_1 <- if (length(diff1@x) > 0) max(abs(diff1@x)) else 0

genes_only_in_recon <- setdiff(rownames(recon_full), rownames(r2_counts))
genes_only_in_r2 <- setdiff(rownames(r2_counts), rownames(recon_full))

writeLines(
  c(
    "STARsolo reconstruction vs ruberto2022_2.rds",
    paste("cells_compared:", length(common_cells_1)),
    paste("genes_compared:", length(common_genes_1)),
    paste("values_compared:", length(common_genes_1) * length(common_cells_1)),
    paste("n_nonzero_mismatches:", n_mismatch_1),
    paste("max_abs_difference:", max_abs_diff_1),
    paste("genes_only_in_starsolo_recon_n:", length(genes_only_in_recon)),
    paste("genes_only_in_ruberto2022_2_rds_n:", length(genes_only_in_r2)),
    paste("cells_missing_from_recon:", length(setdiff(colnames(r2_counts), colnames(recon_full)))),
    paste("cells_missing_from_r2_rds:", length(setdiff(colnames(recon_full), colnames(r2_counts))))
  ),
  file.path(audit_dir, "expression_concordance_starsolo_vs_studyrds.txt")
)

## ---------------------------------------------------------------------
## E. Expression concordance: ruberto2022_2.rds vs deployed raw_df / normalize_df
## ---------------------------------------------------------------------
log("=== E. Expression concordance: study RDS -> deployed raw_df/normalize_df ===")

raw_df <- readRDS("../../PlaViSca/data/raw_df.rds")
r2_rows_in_raw <- rownames(raw_df)[rownames(raw_df) %in% colnames(r2_counts)]
log("Ruberto2022_2 rows found in raw_df: %d / %d", length(r2_rows_in_raw), ncol(r2_counts))

gene_cols <- setdiff(colnames(raw_df), colnames(mr))  # gene columns = not a metadata column of mr_data
# more robust: gene columns are those matching feature ids (start with PVP01 or similar) and not in known metadata
meta_cols_raw <- c(
  "study_pmid", "run_id", "study_label", "num_srr", "pub_year", "geographic_location",
  "sc_technology", "sequencer", "host_species", "host_id", "sample_type", "parasite_stage",
  "strain", "day_post_infection", "treatment", "biological_replicate", "barcode",
  "liver_form", "goegraphic_location", "sequncer", "refine_state", "hour_post_invasion",
  "life_cycle_stage", "seurat_clusters", "pred_gametocyte",
  "umap_u_1", "umap_u_2", "umap_u_3", "umap_i_1", "umap_i_2", "umap_i_3",
  "pca_u_1", "pca_u_2", "pca_u_3", "pca_i_1", "pca_i_2", "pca_i_3"
)
gene_cols <- setdiff(colnames(raw_df), meta_cols_raw)
log("Detected %d gene columns in raw_df", length(gene_cols))

raw_gene_mat <- as.matrix(raw_df[r2_rows_in_raw, gene_cols])
# transpose to genes x cells for comparison, rename genes underscore already
storage.mode(raw_gene_mat) <- "numeric"

common_genes_2 <- intersect(gene_cols, rownames(r2_counts))
common_cells_2 <- r2_rows_in_raw
log("Common genes (study RDS vs deployed raw_df): %d", length(common_genes_2))

a2 <- as.matrix(r2_counts[common_genes_2, common_cells_2])
b2 <- t(raw_gene_mat[common_cells_2, common_genes_2])
diff2 <- a2 - b2
n_mismatch_2 <- sum(diff2 != 0)
max_abs_diff_2 <- if (length(diff2) > 0) max(abs(diff2)) else 0

genes_only_in_r2_not_deployed <- setdiff(rownames(r2_counts), gene_cols)
genes_only_in_deployed_not_r2 <- setdiff(gene_cols, rownames(r2_counts))

writeLines(
  c(
    "ruberto2022_2.rds vs deployed raw_df.rds",
    paste("cells_compared:", length(common_cells_2)),
    paste("genes_compared:", length(common_genes_2)),
    paste("values_compared:", length(common_genes_2) * length(common_cells_2)),
    paste("n_mismatches:", n_mismatch_2),
    paste("max_abs_difference:", max_abs_diff_2),
    paste("genes_only_in_study_rds_not_deployed_n:", length(genes_only_in_r2_not_deployed)),
    paste("genes_only_in_deployed_not_study_rds_n:", length(genes_only_in_deployed_not_r2))
  ),
  file.path(audit_dir, "expression_concordance_studyrds_vs_deployed_raw.txt")
)

## ---------------------------------------------------------------------
## F. Normalized-expression formula check
## ---------------------------------------------------------------------
log("=== F. Normalized expression formula check ===")

normalize_df <- readRDS("../../PlaViSca/data/normalize_df.rds")
norm_gene_mat <- as.matrix(normalize_df[r2_rows_in_raw, gene_cols])
storage.mode(norm_gene_mat) <- "numeric"

cell_totals <- rowSums(raw_gene_mat[r2_rows_in_raw, gene_cols])
expected_norm <- log1p(sweep(raw_gene_mat[r2_rows_in_raw, gene_cols], 1, cell_totals, "/") * 10000)
diff3 <- norm_gene_mat - expected_norm
max_abs_diff_3 <- max(abs(diff3), na.rm = TRUE)
n_mismatch_3 <- sum(abs(diff3) > 1e-8, na.rm = TRUE)

writeLines(
  c(
    "Deployed normalize_df.rds vs log1p(count / cell_total * 10000)",
    paste("cells_compared:", length(r2_rows_in_raw)),
    paste("genes_compared:", length(gene_cols)),
    paste("values_compared:", length(r2_rows_in_raw) * length(gene_cols)),
    paste("n_mismatches_gt_1e-8:", n_mismatch_3),
    paste("max_abs_difference:", max_abs_diff_3)
  ),
  file.path(audit_dir, "normalization_formula_check.txt")
)

## ---------------------------------------------------------------------
## G. Duplicate-profile forensics
## ---------------------------------------------------------------------
log("=== G. Duplicate-profile forensics ===")

hash_col <- function(m, j) {
  col <- m[, j]
  nz <- which(col != 0)
  digest(list(rownames(m)[nz], col[nz]), algo = "sha256")
}

hashes <- vapply(seq_len(ncol(r2_counts)), function(j) hash_col(r2_counts, j), character(1))
names(hashes) <- colnames(r2_counts)

dup_groups <- split(names(hashes), hashes)
dup_groups <- dup_groups[lengths(dup_groups) > 1]

log("Number of exact-duplicate raw-profile groups: %d", length(dup_groups))

if (length(dup_groups) > 0) {
  dup_df <- do.call(rbind, lapply(seq_along(dup_groups), function(i) {
    ids <- dup_groups[[i]]
    data.frame(
      group_id = i,
      cell_id = ids,
      run_id = sub("_.*", "", ids)
    )
  }))
  dup_df$cross_run <- ave(dup_df$run_id, dup_df$group_id, FUN = function(x) length(unique(x)) > 1)
  write.table(dup_df, file.path(audit_dir, "duplicate_profile_groups.tsv"),
    sep = "\t", row.names = FALSE, quote = FALSE)
} else {
  writeLines("no_exact_duplicate_raw_profiles_found", file.path(audit_dir, "duplicate_profile_groups.tsv"))
}

n_cross_run <- if (length(dup_groups) > 0) {
  sum(sapply(dup_groups, function(g) length(unique(sub("_.*", "", g))) > 1))
} else {
  0
}
n_within_run <- length(dup_groups) - n_cross_run

writeLines(
  c(
    paste("n_cells_hashed:", length(hashes)),
    paste("n_duplicate_groups:", length(dup_groups)),
    paste("n_within_run_duplicate_groups:", n_within_run),
    paste("n_cross_run_duplicate_groups:", n_cross_run),
    paste("n_cells_in_duplicate_groups:", sum(lengths(dup_groups)))
  ),
  file.path(audit_dir, "duplicate_profile_summary.txt")
)

# reused barcode tokens across runs (barcode part only, regardless of expression)
bc_tokens <- sub("^[0-9]+_", "", colnames(r2_counts))
run_tokens <- sub("_.*", "", colnames(r2_counts))
tok_tab <- data.frame(barcode = bc_tokens, run = run_tokens)
reused <- tok_tab %>% count(barcode) %>% filter(n > 1)
write.table(reused, file.path(audit_dir, "reused_barcode_tokens_across_runs.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE)
log("Barcode tokens reused across >=2 runs: %d", nrow(reused))

## ---------------------------------------------------------------------
## H. HPI / stage / gametocyte audit on deployed data
## ---------------------------------------------------------------------
log("=== H. HPI / stage / gametocyte audit (deployed cleaned_dataset.rds) ===")

r2_mr <- mr[mr$study_label == "Ruberto2022_2", ]

hpi_tab <- as.data.frame(table(r2_mr$hour_post_invasion, useNA = "always"))
colnames(hpi_tab) <- c("hour_post_invasion", "n_cells")
write.table(hpi_tab, file.path(audit_dir, "deployed_hpi_distribution.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE)

stage_fields <- intersect(
  c("parasite_stage", "life_cycle_stage", "development_phase", "blood_stage", "liver_form", "refine_state", "pred_gametocyte"),
  colnames(r2_mr)
)
stage_summary <- lapply(stage_fields, function(f) as.data.frame(table(r2_mr[[f]], useNA = "always")))
names(stage_summary) <- stage_fields
for (f in stage_fields) {
  colnames(stage_summary[[f]]) <- c(f, "n_cells")
  write.table(stage_summary[[f]], file.path(audit_dir, paste0("deployed_", f, "_distribution.tsv")),
    sep = "\t", row.names = FALSE, quote = FALSE)
}

writeLines(
  c(
    paste("n_deployed_ruberto2022_2_cells:", nrow(r2_mr)),
    paste("n_nonNA_hpi:", sum(!is.na(r2_mr$hour_post_invasion))),
    paste("n_study_rds_cells:", ncol(r2)),
    paste("n_cells_lost_between_study_rds_and_deployed:", ncol(r2) - nrow(r2_mr)),
    paste("stage_fields_present:", paste(stage_fields, collapse = ";"))
  ),
  file.path(audit_dir, "hpi_stage_gametocyte_summary.txt")
)

log("=== Audit script complete ===")
