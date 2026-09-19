#!/usr/bin/env Rscript
# Ruberto2022_1 forensic audit inspection.
# Reads only existing production files; writes only into this audit directory.
# No production script, source data, study object or app object is modified.

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(digest)
})

audit_dir <- "../audit/ruberto2022_1"
proj_root <- ".."

cat("=== 1. Checksums of production inputs ===\n")
files <- c(
  hep59 = file.path(proj_root, "data/Hep59.1.2.seu_20aug2025.rds"),
  study_rds = file.path(proj_root, "ruberto2022_1.rds"),
  gff = "ref/PlasmoDB-68_PvivaxP01.gff"
)
sha <- sapply(files, function(f) tools::md5sum(f))
write.table(
  data.frame(file = files, md5 = sha),
  file.path(audit_dir, "production_inputs_md5.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
print(sha)

cat("\n=== 2. Load Hep59.1.2.seu_20aug2025.rds (author-named local object) ===\n")
hep59 <- readRDS(files["hep59"])
cat("dims:", paste(dim(hep59), collapse = " x "), "\n")
cat("assays:", paste(Assays(hep59), collapse = ","), "\n")
cat("meta.data columns:\n")
print(colnames(hep59@meta.data))
cat("\nLiverForm table:\n")
print(table(hep59$LiverForm, useNA = "ifany"))

meta_cols_of_interest <- intersect(
  c("day", "replicate", "orig.ident", "LiverForm", "seurat_clusters"),
  colnames(hep59@meta.data)
)
cat("\nCross-tab day x LiverForm:\n")
if (all(c("day", "LiverForm") %in% colnames(hep59@meta.data))) {
  print(table(hep59$day, hep59$LiverForm, useNA = "ifany"))
}
cat("\nCross-tab replicate x LiverForm:\n")
if (all(c("replicate", "LiverForm") %in% colnames(hep59@meta.data))) {
  print(table(hep59$replicate, hep59$LiverForm, useNA = "ifany"))
}
cat("\norig.ident table:\n")
if ("orig.ident" %in% colnames(hep59@meta.data)) {
  print(table(hep59$orig.ident, useNA = "ifany"))
}

cat("\nSeurat object version:", as.character(hep59@version), "\n")
cat("commands stored:\n")
print(names(hep59@commands))

# Write full source metadata for cross-referencing
hep59_meta <- hep59@meta.data
hep59_meta$source_cell_id <- rownames(hep59_meta)
write.table(
  hep59_meta,
  file.path(audit_dir, "hep59_source_metadata.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

cat("\n=== 3. Reproduce script's cell-ID mapping/subset logic ===\n")
cells_to_subset <- data.frame(
  cells = rownames(hep59@meta.data),
  liverForm = hep59@meta.data$LiverForm
) |>
  mutate(
    cells_name = sub(".*_", "", cells),
    cells_num = sub("_.*", "", cells),
    new_cells = case_when(
      cells_num == "51inf" ~ "612",
      cells_num == "52inf" ~ "610",
      cells_num == "91inf" ~ "611",
      cells_num == "92inf" ~ "609"
    ),
    new_cells_name = paste0(new_cells, "_", cells_name)
  )

cat("Any cells_num values NOT in the 4 expected groups (unmapped new_cells)?\n")
print(table(cells_to_subset$cells_num, is.na(cells_to_subset$new_cells)))

cat("Total source Hep59 cells:", nrow(cells_to_subset), "\n")
cat("Unique mapped new_cells_name:", length(unique(cells_to_subset$new_cells_name)), "\n")
cat("Duplicated new_cells_name count:", sum(duplicated(cells_to_subset$new_cells_name)), "\n")

write.table(
  cells_to_subset |> select(cells, new_cells_name, liverForm),
  file.path(audit_dir, "hep59_to_study_cellid_map.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

cat("\n=== 4. Load ruberto2022_1.rds (PlaViSca study object) ===\n")
study <- readRDS(file.path(proj_root, "ruberto2022_1.rds"))
cat("dims:", paste(dim(study), collapse = " x "), "\n")
cat("n cells:", ncol(study), "\n")
cat("unique colnames:", length(unique(colnames(study))), "\n")
cat("meta.data columns:\n")
print(colnames(study@meta.data))

cat("\nrun_id table:\n")
print(table(study$run_id, useNA = "ifany"))
cat("\nday_post_infection table:\n")
print(table(study$day_post_infection, useNA = "ifany"))
cat("\ntreatment table:\n")
print(table(study$treatment, useNA = "ifany"))
cat("\nbiological_replicate table:\n")
print(table(study$biological_replicate, useNA = "ifany"))
cat("\nliver_form table (post-join):\n")
print(table(study$liver_form, useNA = "ifany"))
cat("\nNA liver_form count:", sum(is.na(study$liver_form)), "\n")

cat("\nCross tab run_id x day_post_infection x treatment x biological_replicate:\n")
print(
  study@meta.data |>
    count(run_id, day_post_infection, treatment, biological_replicate)
)

cat("\nCross tab run_id x liver_form:\n")
print(table(study$run_id, study$liver_form, useNA = "ifany"))

# Verify barcode set == cells_to_subset new_cells_name set
study_barcodes <- study$barcode
mapped_barcodes <- cells_to_subset$new_cells_name
cat("\nSet comparison, study object barcodes vs Hep59-derived mapped barcodes:\n")
cat("study cells:", length(study_barcodes), "\n")
cat("mapped (Hep59-derived) unique cells:", length(unique(mapped_barcodes)), "\n")
cat("study cells NOT in mapped set:", sum(!study_barcodes %in% mapped_barcodes), "\n")
cat("mapped cells NOT in study set:", sum(!unique(mapped_barcodes) %in% study_barcodes), "\n")

study_meta <- study@meta.data
study_meta$study_object_cell_id <- rownames(study_meta)
write.table(
  study_meta,
  file.path(audit_dir, "study_object_metadata.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

cat("\n=== 5. Exact-duplicate RNA count profiles within ruberto2022_1.rds ===\n")
counts <- GetAssayData(study, layer = "counts")
cat("counts dim:", paste(dim(counts), collapse = " x "), "\n")
cat("counts are integer/nonnegative:", all(counts@x == round(counts@x)), all(counts@x >= 0), "\n")

col_hashes <- apply(counts, 2, function(col) digest(col, algo = "sha256"))
dup_groups <- split(colnames(counts), col_hashes)
dup_groups <- dup_groups[sapply(dup_groups, length) > 1]
cat("Number of exact-duplicate RNA profile groups:", length(dup_groups), "\n")
if (length(dup_groups) > 0) {
  for (g in dup_groups) {
    tot <- sum(counts[, g[1]])
    cat("  group size", length(g), "; profile colSum =", tot, "; first 3:",
        paste(head(g, 3), collapse = ","), "\n")
  }
}
dup_df <- do.call(rbind, lapply(seq_along(dup_groups), function(i) {
  data.frame(group = i, cell = dup_groups[[i]])
}))
write.table(
  if (is.null(dup_df)) data.frame(group = integer(0), cell = character(0)) else dup_df,
  file.path(audit_dir, "exact_duplicate_profiles_ruberto1.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

cat("\n=== 6. Expression preservation vs raw STARsolo GeneFull/raw matrices ===\n")
# rtracklayer/GenomeInfoDb fails to load in this Pixi env (missing
# GenomeInfoDbData); the rRNA ID list was pre-extracted with awk from the same
# GFF file instead (see rRNA_ids.txt), reproducing exactly rRNA_pv$ID in the
# production script (50 rRNA features, same filter: type == "rRNA").
rRNA_pv <- data.frame(ID = readLines(file.path(audit_dir, "rRNA_ids.txt")))

run_map <- c("609" = "SRR19573609", "610" = "SRR19573610", "611" = "SRR19573611", "612" = "SRR19573612")
mismatch_total <- 0
compare_rows <- list()
for (suf in names(run_map)) {
  srr <- run_map[[suf]]
  path <- file.path(proj_root, "counts", "36093191", paste0(srr, "_solo_out/Solo.out/GeneFull/raw"))
  # Plain (non-gzipped) STARsolo output; Seurat::Read10X expects the .gz
  # naming convention, so read the MatrixMarket triplet directly instead.
  raw <- Matrix::readMM(file.path(path, "matrix.mtx"))
  raw_features <- read.delim(file.path(path, "features.tsv"), header = FALSE)
  raw_barcodes <- read.delim(file.path(path, "barcodes.tsv"), header = FALSE)
  rownames(raw) <- raw_features$V1
  colnames(raw) <- raw_barcodes$V1
  raw <- raw[!rownames(raw) %in% rRNA_pv$ID, ]
  study_cols <- colnames(study)[study$run_id == srr]
  study_barcode_suffix <- sub(paste0("^", suf, "_"), "", study_cols)
  present_in_raw <- study_barcode_suffix %in% colnames(raw)
  compare_rows[[suf]] <- data.frame(
    run = srr,
    n_study_cells = length(study_cols),
    n_present_in_raw = sum(present_in_raw)
  )
  if (!all(present_in_raw)) {
    cat("WARNING: missing barcodes in raw matrix for", srr, "\n")
    next
  }
  sub_raw <- raw[, study_barcode_suffix, drop = FALSE]
  sub_study <- counts[, study_cols, drop = FALSE]
  common_genes <- intersect(rownames(sub_raw), rownames(sub_study))
  # study object gene IDs are hyphenated versions of raw underscored IDs (Seurat auto-conversion)
  rn_raw_conv <- gsub("_", "-", rownames(sub_raw))
  m <- match(rn_raw_conv, rownames(sub_study))
  ok <- !is.na(m)
  a <- as.matrix(sub_raw[ok, , drop = FALSE])
  b <- as.matrix(sub_study[m[ok], , drop = FALSE])
  rownames(a) <- rownames(b) <- rn_raw_conv[ok]
  colnames(a) <- study_cols
  diff <- sum(a != b)
  cat(srr, ": genes matched", sum(ok), "of", nrow(sub_raw), "raw genes; unequal entries:", diff, "\n")
  extra_genes_raw <- rownames(sub_raw)[!ok]
  if (length(extra_genes_raw) > 0) {
    cat("  raw genes not found in study object (checking they are all-zero in these cells):",
        length(extra_genes_raw), "; nonzero among them:",
        sum(as.matrix(sub_raw[extra_genes_raw, , drop = FALSE]) != 0), "\n")
  }
  mismatch_total <- mismatch_total + diff
}
cat("\nTotal unequal count entries across all 4 runs:", mismatch_total, "\n")
write.table(
  do.call(rbind, compare_rows),
  file.path(audit_dir, "expression_lineage_summary.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

cat("\n=== 7. Gene-ID collision check (underscore -> hyphen truncation) ===\n")
raw_ids <- rownames(counts)
cat("n features in study object:", length(raw_ids), "\n")
cat("n unique:", length(unique(raw_ids)), "\n")

cat("\n=== 8. Near-empty parasite count profile investigation ===\n")
study_umi <- Matrix::colSums(counts)
cat("Distribution of per-cell total study-object RNA UMI:\n")
print(summary(study_umi))
cat("\nCells with 0 total UMI:", sum(study_umi == 0), "\n")
cat("Cells with 1 total UMI:", sum(study_umi == 1), "\n")
cat("Cells with <=1 total UMI:", sum(study_umi <= 1), "\n")
cat("Cells with <=10 total UMI:", sum(study_umi <= 10), "\n")

near_empty <- names(study_umi)[study_umi <= 1]
by_run <- table(study$run_id[near_empty])
cat("\nNear-empty (<=1 UMI) cells by run_id:\n")
print(by_run)

# Cross-reference against Hep59 (author object) nCount_RNA/nCount_SCT for the
# same cells, via the same cell-ID mapping used in section 3.
hep59_lookup <- cells_to_subset
rownames(hep59_lookup) <- hep59_lookup$new_cells_name
near_empty_hep59 <- hep59@meta.data[hep59_lookup[near_empty, "cells"], c("nCount_RNA", "nCount_SCT", "day", "replicate", "LiverForm")]
near_empty_hep59$study_object_umi <- study_umi[near_empty]
near_empty_hep59$study_cell_id <- near_empty
cat("\nAuthor-object (Hep59) nCount_RNA / nCount_SCT for these same near-empty cells:\n")
print(summary(near_empty_hep59$nCount_RNA))
print(summary(near_empty_hep59$nCount_SCT))

# Also compare for a set of NON-near-empty cells as a control
non_empty <- names(study_umi)[study_umi > 1]
non_empty_hep59_nCount <- hep59@meta.data[hep59_lookup[non_empty, "cells"], "nCount_RNA"]
cat("\nControl: author-object nCount_RNA for cells NOT near-empty in study object:\n")
print(summary(non_empty_hep59_nCount))

write.table(
  near_empty_hep59,
  file.path(audit_dir, "near_empty_cells_vs_hep59.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

cat("\n=== DONE inspect_ruberto2022_1.R ===\n")
