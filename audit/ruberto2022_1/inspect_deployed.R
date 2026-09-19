#!/usr/bin/env Rscript
# Ruberto2022_1: lineage/concordance check against the deployed PlaViSca app data.
# Reads only existing production/app files; writes only into this audit directory.

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
})

audit_dir <- "../audit/ruberto2022_1"
proj_root <- ".."

cat("=== Load deployed cleaned_dataset.rds ===\n")
cleaned <- readRDS(file.path(proj_root, "..", "PlaViSca/data/cleaned_dataset.rds"))
cat("names(cleaned):", paste(names(cleaned), collapse = ", "), "\n")
mr <- cleaned$mr_data
cat("mr_data dim:", paste(dim(mr), collapse = " x "), "\n")
cat("study_label values:\n")
print(table(mr$study_label, useNA = "ifany"))

rub <- mr[mr$study_label == "Ruberto2022_1", ]
cat("\nRuberto2022_1 rows in deployed mr_data:", nrow(rub), "\n")

cat("\n=== Load ruberto2022_1.rds study object for comparison ===\n")
study <- readRDS(file.path(proj_root, "ruberto2022_1.rds"))
study_cells <- colnames(study)

cat("Deployed Ruberto rownames vs study object barcodes:\n")
cat("deployed rows:", nrow(rub), "; study cells:", length(study_cells), "\n")
cat("deployed rownames NOT in study barcodes:", sum(!rownames(rub) %in% study_cells), "\n")
cat("study barcodes NOT in deployed rownames:", sum(!study_cells %in% rownames(rub)), "\n")

cat("\n=== Deployed metadata field concordance ===\n")
common_cells <- intersect(rownames(rub), study_cells)
cat("common cells:", length(common_cells), "\n")
rub_c <- rub[common_cells, ]
study_meta <- study@meta.data[common_cells, ]

cat("\ndeployed study_pmid table:\n"); print(table(rub_c$study_pmid, useNA = "ifany"))
cat("\ndeployed pub_year table:\n"); print(table(rub_c$pub_year, useNA = "ifany"))
cat("\ndeployed run_id table:\n"); print(table(rub_c$run_id, useNA = "ifany"))
cat("\ndeployed day_post_infection table:\n"); print(table(rub_c$day_post_infection, useNA = "ifany"))
cat("\ndeployed treatment table:\n"); print(table(rub_c$treatment, useNA = "ifany"))
cat("\ndeployed geographic_location table:\n"); print(table(rub_c$geographic_location, useNA = "ifany"))
cat("\ndeployed sc_technology table:\n"); print(table(rub_c$sc_technology, useNA = "ifany"))
cat("\ndeployed sequencer table:\n"); print(table(rub_c$sequencer, useNA = "ifany"))
cat("\ndeployed host_id table:\n"); print(table(rub_c$host_id, useNA = "ifany"))
cat("\ndeployed strain table:\n"); print(table(rub_c$strain, useNA = "ifany"))
cat("\ndeployed sample_type table:\n"); print(table(rub_c$sample_type, useNA = "ifany"))

cat("\n=== Final stage/life-cycle field(s) present ===\n")
stage_cols <- grep("life_cycle|stage|liver|schizont|hypnozoite|refine|pred_gameto", colnames(mr), ignore.case = TRUE, value = TRUE)
print(stage_cols)
for (col in stage_cols) {
  cat("\n", col, "table for Ruberto2022_1 rows:\n")
  print(table(rub_c[[col]], useNA = "ifany"))
}

cat("\n=== HPI field ===\n")
hpi_cols <- grep("hour_post|hpi", colnames(mr), ignore.case = TRUE, value = TRUE)
print(hpi_cols)
for (col in hpi_cols) {
  cat("\n", col, "for Ruberto2022_1 rows, NA count:", sum(is.na(rub_c[[col]])), "of", nrow(rub_c), "\n")
  print(table(rub_c[[col]], useNA = "ifany"))
}

cat("\n=== barcode field ===\n")
if ("barcode" %in% colnames(mr)) {
  cat("NA count:", sum(is.na(rub_c$barcode)), "of", nrow(rub_c), "\n")
} else {
  cat("No 'barcode' column present in deployed mr_data.\n")
}

cat("\n=== liver_form field presence in deployed data ===\n")
cat("'liver_form' %in% colnames(mr):", "liver_form" %in% colnames(mr), "\n")

cat("\n=== Compare study liver_form to deployed final stage, cell by cell ===\n")
liver_form_study <- study_meta$liver_form
final_stage_col <- intersect(c("life_cycle_stage", "parasite_stages"), colnames(mr))
cat("Available final stage columns:", paste(final_stage_col, collapse = ","), "\n")
for (col in final_stage_col) {
  cat("\nCross-tab study liver_form x deployed", col, ":\n")
  print(table(liver_form_study, rub_c[[col]], useNA = "ifany"))
}

cat("\n=== Expression concordance for the near-empty cells in the deployed raw_df ===\n")
raw_df <- readRDS(file.path(proj_root, "..", "PlaViSca/data/raw_df.rds"))
cat("raw_df dim:", paste(dim(raw_df), collapse = " x "), "\n")
study_umi <- Matrix::colSums(GetAssayData(study, layer = "counts"))
near_empty <- names(study_umi)[study_umi <= 1]
near_empty_in_raw_df <- intersect(near_empty, rownames(raw_df))
cat("near-empty study cells also present in deployed raw_df:", length(near_empty_in_raw_df), "of", length(near_empty), "\n")
if (length(near_empty_in_raw_df) > 0) {
  gene_cols <- setdiff(colnames(raw_df), colnames(mr))
  sub <- raw_df[near_empty_in_raw_df, gene_cols, drop = FALSE]
  sums <- rowSums(as.matrix(sub))
  cat("deployed raw_df row-sums (gene columns only) for these cells:\n")
  print(summary(sums))
}

write.table(
  rub_c,
  file.path(audit_dir, "deployed_ruberto1_metadata.tsv"),
  sep = "\t", quote = FALSE, row.names = TRUE
)

cat("\n=== DONE inspect_deployed.R ===\n")
