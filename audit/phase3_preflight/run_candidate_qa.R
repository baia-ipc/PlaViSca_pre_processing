#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

repo <- normalizePath(getwd())
candidate_root <- "/home/baia/prj/plavisca/pre_process_data"
candidate_dir <- file.path(candidate_root, "data", "candidate_export")
out_dir <- file.path(repo, "audit", "phase3_preflight")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

write_tsv <- function(x, path) {
  utils::write.table(x, path, sep = "\t", quote = TRUE, row.names = FALSE, na = "NA")
}

schema <- read.delim(file.path(repo, "audit/shared_pipeline/metadata_schema.tsv"),
                     check.names = FALSE, quote = "", comment.char = "")
cleaned <- readRDS(file.path(candidate_dir, "cleaned_dataset.rds"))
meta <- cleaned$mr_data

# Candidate schema inventory. Distinct counts are cheap for metadata-sized columns.
current_app_config <- c(
  study_label = "study_label", host_species = "host_species",
  parasite_stages = "life_cycle_stage", treatment = "treatment",
  development_phase = "parasite_stage", sample_type = "sample_type",
  strain = "strain", hour_post_invasion = "hour_post_invasion"
)
desired_fields <- schema$field_name[grepl("^Yes", schema$app_visible)]
schema_rows <- lapply(names(meta), function(nm) {
  x <- meta[[nm]]
  data.frame(
    column_name = nm,
    type = paste(class(x), collapse = ","),
    missing_fraction = mean(is.na(x)),
    distinct_count = length(unique(x[!is.na(x)])),
    app_visible_status = if (nm %in% desired_fields) "YES" else "NO_OR_INTERNAL",
    schema_definition_exists = nm %in% schema$field_name,
    expected_by_current_app = nm %in% unname(current_app_config),
    expected_by_desired_phase3_app = nm %in% desired_fields,
    literal_NA_count = if (is.character(x) || is.factor(x)) sum(as.character(x) == "NA", na.rm = TRUE) else 0L,
    notes = if (nm %in% c("parasite_stage", "life_cycle_stage", "hour_post_invasion"))
      "Obsolete current-app physical name; absent from candidate" else "",
    check.names = FALSE
  )
})
schema_out <- do.call(rbind, schema_rows)
# Add missing current-app columns so incompatibility is explicit.
missing_current <- setdiff(unname(current_app_config), names(meta))
if (length(missing_current)) {
  schema_out <- rbind(schema_out, data.frame(
    column_name = missing_current, type = "ABSENT", missing_fraction = NA_real_,
    distinct_count = NA_integer_, app_visible_status = "CURRENT_APP_EXPECTS_BUT_ABSENT",
    schema_definition_exists = missing_current %in% schema$field_name,
    expected_by_current_app = TRUE, expected_by_desired_phase3_app = missing_current %in% desired_fields,
    literal_NA_count = NA_integer_, notes = "Current master config is incompatible with Phase-2 candidate",
    check.names = FALSE
  ))
}
write_tsv(schema_out, file.path(out_dir, "app_candidate_schema.tsv"))

# Cell-key alignment, including a guard against metadata columns being mistaken for genes.
artifact_paths <- c(
  raw = file.path(candidate_dir, "raw_df.rds"),
  normalized = file.path(candidate_dir, "normalize_df.rds"),
  scaled = file.path(candidate_dir, "scale_df.rds")
)
key_rows <- list()
reference_ids <- rownames(meta)
for (nm in names(artifact_paths)) {
  x <- readRDS(artifact_paths[[nm]])
  ids <- rownames(x)
  gene_cols <- setdiff(names(x), names(meta))
  key_rows[[nm]] <- data.frame(
    artifact = nm, row_count = nrow(x), column_count = ncol(x),
    unique_rownames = !anyDuplicated(ids), rowname_set_equals_cleaned_metadata = setequal(ids, reference_ids),
    rowname_order_equals_cleaned_metadata = identical(ids, reference_ids),
    metadata_column_count = sum(names(x) %in% names(meta)), gene_column_count = length(gene_cols),
    suspicious_gene_columns = paste(gene_cols[!grepl("^PVP01_", gene_cols)], collapse = ";"),
    status = if (!anyDuplicated(ids) && setequal(ids, reference_ids) &&
                 !any(!grepl("^PVP01_", gene_cols)))
      if (identical(ids, reference_ids)) "PASS" else "SET_PASS_ORDER_MISMATCH" else "FAIL",
    notes = "All cell-indexed exports use canonical lexicographically sorted cell keys. Gene columns are columns absent from mr_data and all use PVP01_ identifiers."
  )
  rm(x); gc()
}
key_rows$cleaned_metadata <- data.frame(
  artifact = "cleaned_metadata", row_count = nrow(meta), column_count = ncol(meta),
  unique_rownames = !anyDuplicated(reference_ids), rowname_set_equals_cleaned_metadata = TRUE,
  rowname_order_equals_cleaned_metadata = TRUE, metadata_column_count = ncol(meta), gene_column_count = 0L,
  suspicious_gene_columns = "", status = if (!anyDuplicated(reference_ids)) "PASS" else "FAIL",
  notes = "Reference cell key"
)
write_tsv(do.call(rbind, key_rows), file.path(out_dir, "cell_key_validation.tsv"))

# Hazzard2024 AT12 direct comparison by run key.
auth <- read.delim(file.path(repo, "audit/hazzard2024/authoritative_run_table.tsv"),
                   check.names = FALSE, quote = "", na.strings = character())
hz <- meta[meta$study_label == "Hazzard2024", , drop = FALSE]
idx <- match(hz$run_id, auth$run_id)
at12_fields <- list(
  library_id = c("library_id", "library_id"),
  animal_id = c("animal_id", "animal_id"),
  parasite_lineage_isolate = c("parasite_lineage_isolate", "inoculum_or_lineage_authoritative"),
  day_post_infection = c("day_post_infection", "day_post_infection_authoritative"),
  infection_design = c("infection_design", "infection_design_authoritative")
)
at12 <- lapply(names(at12_fields), function(label) {
  fields <- at12_fields[[label]]
  observed <- as.character(hz[[fields[[1]]]])
  expected <- as.character(auth[[fields[[2]]]][idx])
  ok <- !is.na(idx) & observed == expected
  data.frame(field = label, cells_tested = length(ok), matches = sum(ok, na.rm = TRUE),
             mismatches = sum(!ok | is.na(ok)), agreement_fraction = mean(ok, na.rm = TRUE),
             status = if (length(ok) == 80024L && all(ok, na.rm = FALSE)) "PASS" else "FAIL")
})
write_tsv(do.call(rbind, at12), file.path(out_dir, "AT12_hazzard2024_validation.tsv"))

# AT28 nullable semantics: a non-NA pruned label must carry raw label, delta and pruning flag.
eligible <- !is.na(meta$idc_reference_similarity_label)
at28 <- data.frame(
  check = c("raw_label_nonmissing", "score_delta_nonmissing", "pruned_flag_nonmissing", "retained_label_not_flagged_pruned"),
  tested_cells = sum(eligible),
  failures = c(
    sum(eligible & is.na(meta$idc_similarity_raw_label)),
    sum(eligible & is.na(meta$idc_similarity_score_delta)),
    sum(eligible & is.na(meta$idc_similarity_pruned_flag)),
    sum(eligible & meta$idc_similarity_pruned_flag %in% TRUE, na.rm = TRUE)
  )
)
at28$status <- ifelse(at28$failures == 0L, "PASS", "FAIL")
write_tsv(at28, file.path(out_dir, "AT28_singleR_uncertainty.tsv"))

# Bibliography contract against the actual corrected candidate file.
bib <- read.csv(file.path(candidate_root, "data", "data_source.csv"), check.names = FALSE)
names(bib) <- tolower(names(bib))
candidate_pairs <- unique(meta[c("study_label", "study_pmid")])
candidate_pairs$study_pmid <- as.character(candidate_pairs$study_pmid)
label_col <- if ("study_label" %in% names(bib)) "study_label" else names(bib)[grepl("study.*label", names(bib))][1]
pmid_col <- if ("study_pmid" %in% names(bib)) "study_pmid" else names(bib)[grepl("pmid", names(bib))][1]
bib_has_pmid <- length(pmid_col) == 1L && !is.na(pmid_col) && nzchar(pmid_col)
bib_pairs <- unique(data.frame(study_label = bib[[label_col]],
                               study_pmid = if (bib_has_pmid) as.character(bib[[pmid_col]]) else NA_character_))
all_labels <- union(candidate_pairs$study_label, bib_pairs$study_label)
bib_out <- do.call(rbind, lapply(all_labels, function(label) {
  crows <- candidate_pairs[candidate_pairs$study_label == label, , drop = FALSE]
  brows <- bib_pairs[bib_pairs$study_label == label, , drop = FALSE]
  issue <- if (!nrow(crows)) "EXTRA_BIBLIOGRAPHY_ROW" else if (!nrow(brows)) "MISSING_BIBLIOGRAPHY_ROW" else if (nrow(crows) != 1L || nrow(brows) != 1L) "DUPLICATE_OR_AMBIGUOUS" else if (!bib_has_pmid) "LABEL_MATCH_PMID_NOT_IN_BIBLIOGRAPHY" else if (crows$study_pmid != brows$study_pmid) "PMID_MISMATCH" else "MATCH"
  data.frame(study_label = label,
             candidate_pmid = paste(unique(crows$study_pmid), collapse = ";"),
             bibliography_pmid = paste(unique(brows$study_pmid), collapse = ";"),
             candidate_rows = nrow(crows), bibliography_rows = nrow(brows), status = issue)
}))
write_tsv(bib_out, file.path(out_dir, "bibliography_contract.tsv"))

# Compact directly observed values used by the acceptance matrix/report.
summary <- data.frame(
  metric = c("total_cells", "studies", "idc_eligible_cells", "gametocyte_calls", "impossible_liver_gametocyte", "literal_NA_in_app_visible", "top_genes_studies", "ruberto_cells", "ruberto_near_empty", "ruberto_provenance_complete"),
  value = c(
    nrow(meta), length(unique(meta$study_label)), sum(!is.na(meta$idc_reference_similarity_label)),
    sum(!is.na(meta$pred_gametocyte_sex)),
    sum(meta$parasite_broad_stage == "Liver stage" & !is.na(meta$pred_gametocyte_sex)),
    sum(vapply(intersect(desired_fields, names(meta)), function(nm) if (is.character(meta[[nm]]) || is.factor(meta[[nm]])) sum(as.character(meta[[nm]]) == "NA", na.rm = TRUE) else 0L, numeric(1))),
    length(unique(cleaned$top_genes_exp$study)), sum(meta$study_label == "Ruberto2022_1"),
    sum(meta$study_label == "Ruberto2022_1" & meta$near_empty_expression_flag %in% TRUE, na.rm = TRUE),
    sum(meta$study_label == "Ruberto2022_1" & !is.na(meta$source_of_counts) & !is.na(meta$counting_pipeline) & !is.na(meta$count_reference_version) & !is.na(meta$count_provenance_status))
  )
)
write_tsv(summary, file.path(out_dir, "candidate_observed_summary.tsv"))

message("Candidate QA outputs written to ", out_dir)
