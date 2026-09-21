#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
  library(readxl)
  library(xml2)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
audit_dir <- normalizePath(dirname(sub("^--file=", "", script_arg[[1]])))
pre_dir <- normalizePath(file.path(audit_dir, "../.."))
counts_dir <- file.path(pre_dir, "counts", "32365102")
runs <- sprintf("SRR11008%d", 269:278)
prefix_num <- setNames(as.character(269:278), runs)

write_tsv <- function(x, name) {
  write.table(x, file.path(audit_dir, name), sep = "\t", quote = FALSE,
              row.names = FALSE, na = "NA")
}

source <- readRDS(file.path(pre_dir, "sa2020.rds"))
source_counts <- LayerData(source, assay = "RNA", layer = "counts")
source_ids <- colnames(source)

final <- readRDS(file.path(pre_dir, "../PlaViSca/data/cleaned_dataset.rds"))$mr_data
final <- final[final$study_label == "Sa2020", , drop = FALSE]

# Exact STARsolo -> sa2020.rds count and identifier lineage.
lineage_rows <- list()
all_expected_ids <- character()
all_matrix_columns <- list()
feature_totals <- NULL
for (run in runs) {
  d <- file.path(counts_dir, paste0(run, "_solo_out"), "Solo.out", "GeneFull", "filtered")
  barcodes <- readLines(file.path(d, "barcodes.tsv"))
  features <- read.delim(file.path(d, "features.tsv"), header = FALSE,
                         stringsAsFactors = FALSE, check.names = FALSE)
  mat <- readMM(file.path(d, "matrix.mtx"))
  rownames(mat) <- features[[1]]
  colnames(mat) <- barcodes
  if (is.null(feature_totals)) feature_totals <- setNames(numeric(nrow(mat)), rownames(mat))
  feature_totals[rownames(mat)] <- feature_totals[rownames(mat)] + Matrix::rowSums(mat)
  expected_ids <- paste0(prefix_num[[run]], "_", barcodes)
  all_expected_ids <- c(all_expected_ids, expected_ids)

  # Seurat replaces underscores in feature names with hyphens. Compare every
  # retained feature and cell after exact row/column alignment.
  transformed_features <- gsub("_", "-", rownames(mat), fixed = TRUE)
  keep <- transformed_features %in% rownames(source_counts)
  compare_mat <- mat[keep, , drop = FALSE]
  rownames(compare_mat) <- transformed_features[keep]
  compare_mat <- compare_mat[rownames(source_counts), , drop = FALSE]
  object_mat <- source_counts[, expected_ids, drop = FALSE]
  delta <- object_mat - compare_mat

  all_matrix_columns[[run]] <- compare_mat
  lineage_rows[[run]] <- data.frame(
    run_id = run,
    numeric_prefix = prefix_num[[run]],
    starsolo_cells = length(barcodes),
    starsolo_features = nrow(mat),
    retained_object_features = nrow(compare_mat),
    exact_id_order = identical(colnames(object_mat), expected_ids),
    source_object_cells = sum(source$run_id == run),
    missing_from_source_object = length(setdiff(expected_ids, source_ids)),
    extra_in_source_object = length(setdiff(source_ids[source$run_id == run], expected_ids)),
    unequal_count_entries = sum(delta@x != 0),
    source_umi = sum(mat),
    retained_feature_umi = sum(compare_mat),
    object_umi = sum(object_mat),
    stringsAsFactors = FALSE
  )
}
lineage <- do.call(rbind, lineage_rows)
write_tsv(lineage, "sa2020_count_lineage.tsv")

gff <- read.delim(file.path(pre_dir, "scripts", "ref", "PlasmoDB-68_PvivaxP01.gff"),
                  header = FALSE, comment.char = "#", quote = "", stringsAsFactors = FALSE)
rr_attr <- gff[gff[[3]] == "rRNA", 9]
rr_id <- sub(";.*$", "", sub("^ID=", "", rr_attr))
rr_id <- gsub("\\.1", "", rr_id)
feature_disposition <- data.frame(
  source_feature_id = names(feature_totals),
  seurat_feature_id = gsub("_", "-", names(feature_totals), fixed = TRUE),
  total_umi = as.numeric(feature_totals),
  annotated_rrna_removed = names(feature_totals) %in% rr_id,
  retained_in_source_object = gsub("_", "-", names(feature_totals), fixed = TRUE) %in% rownames(source_counts),
  stringsAsFactors = FALSE
)
write_tsv(feature_disposition, "sa2020_feature_disposition.tsv")

# Full cell-set identity and final metadata preservation.
shared_meta <- intersect(names(source@meta.data), names(final))
meta_rows <- lapply(shared_meta, function(field) {
  a <- as.character(source@meta.data[source_ids, field])
  b <- as.character(final[source_ids, field])
  data.frame(field = field, compared = length(a), disagreements = sum(a != b, na.rm = TRUE),
             source_na = sum(is.na(a)), final_na = sum(is.na(b)), stringsAsFactors = FALSE)
})
write_tsv(do.call(rbind, meta_rows), "sa2020_source_final_metadata.tsv")

# Validated author-prefix bridge. Table 1's sample identities and parasite
# totals, ENA aliases, and barcode overlaps are all retained in the output.
bridge <- data.frame(
  author_prefix = c("NIH_Sa", "NIH_CQ", "NIH_Ao", "PB_MACS", "Ches_Sa",
                    "Ches_CQ", "Ches_Ao", "AMRU_Sa", "AMRU_CQ", "AMRU_Ao"),
  run_id = runs,
  strain = c("NIH-1993", "Indonesia-I/CDC", "Indonesia-I/CDC", "NIH-1993",
             "Chesson", "Chesson", "Chesson", "AMRU-I", "AMRU-I", "AMRU-I"),
  host_species = c("Saimiri boliviensis", "Aotus nancymaae", "Aotus nancymaae",
                   "Saimiri boliviensis", "Saimiri boliviensis", "Aotus nancymaae",
                   "Aotus nancymaae", "Saimiri boliviensis", "Aotus nancymaae",
                   "Aotus nancymaae"),
  host_id = c("3879", "86574", "86574", "5541", "4215", "86436", "86436",
              "5107", "86416", "86416"),
  dpi = c(15, 16, 12, 21, 24, 16, 12, 24, 16, 12),
  treatment = c("untreated", "CQ 16h 5 mg/kg", "untreated", "untreated", "untreated",
                "CQ 16h 10 mg/kg", "untreated", "untreated", "CQ 16h 10 mg/kg", "untreated"),
  publication_parasites = c(928, 589, 2098, 267, 249, 521, 1037, 22, 1795, 1709),
  stringsAsFactors = FALSE
)

fig1 <- read_excel(file.path(audit_dir, "sa2020_s030.xlsx"), sheet = "data 1B")
fig1_ids <- as.character(fig1$GEM[!grepl("^Pberg_", fig1$GEM)])
fig1_prefix <- sub("_[^_]+$", "", fig1_ids)
fig1_barcode <- sub("^.*_", "", fig1_ids)
main <- data.frame(author_cell_id = fig1_ids, author_prefix = fig1_prefix,
                   barcode = fig1_barcode, stringsAsFactors = FALSE)
main$run_id <- bridge$run_id[match(main$author_prefix, bridge$author_prefix)]
main$plavisca_cell_id <- paste0(prefix_num[main$run_id], "_", main$barcode)
main$plavisca_present <- main$plavisca_cell_id %in% source_ids

for (i in seq_len(nrow(bridge))) {
  p <- bridge$author_prefix[i]
  run <- bridge$run_id[i]
  author_bc <- unique(main$barcode[main$author_prefix == p])
  star_bc <- sub("^[^_]+_", "", source_ids[source$run_id == run])
  bridge$figure1_cells[i] <- sum(main$author_prefix == p)
  bridge$starsolo_cells[i] <- length(star_bc)
  bridge$barcode_overlap[i] <- length(intersect(author_bc, star_bc))
  bridge$author_missing_from_plavisca[i] <- length(setdiff(author_bc, star_bc))
  bridge$plavisca_outside_author_main[i] <- length(setdiff(star_bc, author_bc))
}
write_tsv(bridge, "sa2020_authoritative_sample_crosswalk.tsv")
write_tsv(main, "sa2020_main_analysis_cell_lineage.tsv")

alt_raw <- read_excel(file.path(audit_dir, "sa2020_s036.xlsx"))
alt <- data.frame(author_cell_id = as.character(alt_raw$Cell), stringsAsFactors = FALSE)
alt$author_prefix <- sub("_[^_]+$", "", alt$author_cell_id)
alt$barcode <- sub("^.*_", "", alt$author_cell_id)
alt$run_id <- bridge$run_id[match(alt$author_prefix, bridge$author_prefix)]
alt$plavisca_cell_id <- paste0(prefix_num[alt$run_id], "_", alt$barcode)
alt$plavisca_present <- alt$plavisca_cell_id %in% source_ids
write_tsv(alt, "sa2020_lower_threshold_cell_lineage.tsv")

# Author Fig 5 annotations and the two independent code defects: recycled
# equality and incorrect/omitted prefix conversion.
gam <- read_excel(file.path(audit_dir, "sa2020_s034.xlsx"))
gam <- data.frame(sample = as.character(gam$Sample), source_type = as.character(gam$Type),
                  stringsAsFactors = FALSE)
gam$author_prefix <- sub("_[^_]+$", "", gam$sample)
gam$barcode <- sub("^.*_", "", gam$sample)
gam$source_label <- ifelse(gam$source_type == "Female", "Female gametocyte",
                           ifelse(gam$source_type == "Male", "Male gametocyte", "Asexual"))
gam$row_number <- seq_len(nrow(gam))
target <- c("Male gametocyte", "Female gametocyte")
gam$recycled_filter_retained <- gam$source_label == rep(target, length.out = nrow(gam))
gam$correct_filter_retained <- gam$source_label %in% target
gam$run_id <- bridge$run_id[match(gam$author_prefix, bridge$author_prefix)]
gam$correct_cell_id <- paste0(prefix_num[gam$run_id], "_", gam$barcode)
wrong_prefix <- c(AMRU_Ao="278", AMRU_CQ="277", AMRU_Sa="276", Ches_Ao="275",
                  Ches_CQ="274", Ches_Sa="273", NIH_Ao="269", NIH_CQ="269", NIH_Sa="272")
gam$production_cell_id <- ifelse(gam$author_prefix %in% names(wrong_prefix),
                                 paste0(wrong_prefix[gam$author_prefix], "_", gam$barcode), NA)
gam$correct_id_present <- gam$correct_cell_id %in% source_ids
gam$production_id_present <- gam$production_cell_id %in% source_ids
gam$final_stage <- as.character(final[gam$correct_cell_id, "life_cycle_stage"])
write_tsv(gam, "sa2020_source_annotation_lineage.tsv")

annotation_comparison <- as.data.frame(table(
  author_prefix = gam$author_prefix[gam$correct_filter_retained & gam$correct_id_present],
  source_label = gam$source_label[gam$correct_filter_retained & gam$correct_id_present],
  final_stage = gam$final_stage[gam$correct_filter_retained & gam$correct_id_present]
), stringsAsFactors = FALSE)
annotation_comparison <- annotation_comparison[annotation_comparison$Freq > 0, ]
write_tsv(annotation_comparison, "sa2020_source_vs_final_gametocytes.tsv")

production_effect <- data.frame(
  author_prefix = unique(gam$author_prefix), stringsAsFactors = FALSE
)
for (i in seq_len(nrow(production_effect))) {
  p <- production_effect$author_prefix[i]
  z <- gam$author_prefix == p
  production_effect$source_sex_rows[i] <- sum(z & gam$correct_filter_retained)
  production_effect$source_sex_present_correct_bridge[i] <- sum(z & gam$correct_filter_retained & gam$correct_id_present)
  production_effect$recycled_rows[i] <- sum(z & gam$recycled_filter_retained)
  production_effect$recycled_present_correct_bridge[i] <- sum(z & gam$recycled_filter_retained & gam$correct_id_present)
  production_effect$production_code_matches[i] <- sum(z & gam$recycled_filter_retained & gam$production_id_present)
  production_effect$production_mapping_is_correct[i] <- p %in% c("AMRU_Ao", "AMRU_CQ", "AMRU_Sa", "Ches_Ao", "Ches_CQ", "Ches_Sa")
}
write_tsv(production_effect, "sa2020_gametocyte_code_effect_by_prefix.tsv")

sex <- gam[gam$correct_filter_retained, ]
summary_rows <- data.frame(
  measure = c("workbook_rows", "workbook_unique_ids", "workbook_ids_in_main_list",
              "main_list_ids_absent_from_workbook", "workbook_ids_outside_main_list",
              "author_asexual_rows", "author_female_rows", "author_male_rows",
              "recycled_filter_retained", "recycled_filter_female", "recycled_filter_male",
              "correct_filter_retained", "correct_bridge_sex_matches",
              "correct_bridge_sex_absent", "correct_bridge_recycled_matches",
              "correct_bridge_omitted_by_recycling_among_present",
              "production_code_matches", "correct_bridge_matches_labeled_same_in_final",
              "correct_bridge_matches_labeled_differently_in_final"),
  value = c(nrow(gam), length(unique(gam$sample)), sum(gam$sample %in% main$author_cell_id),
            length(setdiff(main$author_cell_id, gam$sample)), length(setdiff(gam$sample, main$author_cell_id)),
            sum(gam$source_label == "Asexual"), sum(gam$source_label == "Female gametocyte"),
            sum(gam$source_label == "Male gametocyte"), sum(gam$recycled_filter_retained),
            sum(gam$recycled_filter_retained & gam$source_label == "Female gametocyte"),
            sum(gam$recycled_filter_retained & gam$source_label == "Male gametocyte"),
            sum(gam$correct_filter_retained), sum(sex$correct_id_present), sum(!sex$correct_id_present),
            sum(gam$recycled_filter_retained & gam$correct_id_present),
            sum(gam$correct_filter_retained & !gam$recycled_filter_retained & gam$correct_id_present),
            sum(gam$recycled_filter_retained & gam$production_id_present),
            sum(sex$correct_id_present & sex$final_stage == sex$source_label),
            sum(sex$correct_id_present & sex$final_stage != sex$source_label)),
  stringsAsFactors = FALSE
)
write_tsv(summary_rows, "sa2020_gametocyte_filter_summary.tsv")

# Detect exact duplicate expression columns in the source object. A collision is
# confirmed by a final sparse-column equality check, not by the signature alone.
signature <- vapply(seq_len(ncol(source_counts)), function(j) {
  col <- source_counts[, j, drop = FALSE]
  paste(paste(col@i, format(col@x, scientific = FALSE), sep = ":"), collapse = ",")
}, character(1))
groups <- split(seq_along(signature), signature)
groups <- groups[lengths(groups) > 1]
dup_rows <- list()
k <- 0L
for (g in groups) {
  for (a in seq_len(length(g) - 1L)) for (b in (a + 1L):length(g)) {
    if (identical(source_counts[, g[a]], source_counts[, g[b]])) {
      k <- k + 1L
      dup_rows[[k]] <- data.frame(cell_1 = source_ids[g[a]], run_1 = source$run_id[g[a]],
                                  cell_2 = source_ids[g[b]], run_2 = source$run_id[g[b]],
                                  same_barcode = sub("^[^_]+_", "", source_ids[g[a]]) == sub("^[^_]+_", "", source_ids[g[b]]),
                                  umi = sum(source_counts[, g[a]]), stringsAsFactors = FALSE)
    }
  }
}
dups <- if (length(dup_rows)) do.call(rbind, dup_rows) else
  data.frame(cell_1=character(),run_1=character(),cell_2=character(),run_2=character(),same_barcode=logical(),umi=numeric())
write_tsv(dups, "sa2020_exact_expression_duplicates.tsv")

inventory <- data.frame(
  measure = c("starsolo_total", "source_object_cells", "source_object_unique_ids",
              "final_cells", "final_unique_ids", "source_missing_from_final", "final_extra_vs_source",
              "publication_main_cells", "publication_main_unique_ids", "publication_main_present",
              "publication_main_missing", "plavisca_outside_publication_main",
              "lower_threshold_supplement_cells", "lower_threshold_supplement_unique_ids",
              "lower_threshold_supplement_present", "plavisca_outside_lower_threshold_supplement",
              "main_ids_absent_from_lower_threshold_supplement",
              "duplicate_expression_pairs", "barcode_tokens_reused_across_runs",
              "maximum_run_multiplicity_per_barcode", "input_features", "retained_features",
              "rrna_features_removed", "zero_count_features_removed"),
  value = c(length(all_expected_ids), ncol(source), length(unique(source_ids)), nrow(final),
            length(unique(rownames(final))), length(setdiff(source_ids, rownames(final))),
            length(setdiff(rownames(final), source_ids)), nrow(main), length(unique(main$author_cell_id)),
            sum(main$plavisca_present), sum(!main$plavisca_present),
            length(setdiff(source_ids, main$plavisca_cell_id)), nrow(alt), length(unique(alt$author_cell_id)),
            sum(alt$plavisca_present), length(setdiff(source_ids, alt$plavisca_cell_id)),
            length(setdiff(main$author_cell_id, alt$author_cell_id)), nrow(dups),
            sum(table(sub("^[^_]+_", "", source_ids)) > 1),
            max(table(sub("^[^_]+_", "", source_ids))), nrow(feature_disposition), nrow(source_counts),
            sum(feature_disposition$annotated_rrna_removed),
            sum(feature_disposition$total_umi == 0 & !feature_disposition$annotated_rrna_removed)),
  stringsAsFactors = FALSE
)
write_tsv(inventory, "sa2020_inventory_summary.tsv")

cat("Sa2020 audit completed\n")
print(inventory)
print(summary_rows)
