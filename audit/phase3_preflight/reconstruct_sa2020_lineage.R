#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(Seurat))

repo <- normalizePath(getwd())
candidate_root <- Sys.getenv("PLAVISCA_CANDIDATE_ROOT", repo)
out_path <- file.path(repo, "audit/mancio_silva2022/sa2020_authoritative_cell_lineage.tsv")

study <- readRDS(file.path(candidate_root, "sa2020.rds"))
final <- readRDS(file.path(candidate_root, "data/candidate_export/cleaned_dataset.rds"))$mr_data
final_ids <- rownames(final)[final$study_label == "Sa2020"]

bridge <- read.delim(file.path(repo, "audit/mancio_silva2022/sa2020_authoritative_sample_crosswalk.tsv"),
                     check.names = FALSE, quote = "", stringsAsFactors = FALSE)
main <- read.delim(file.path(repo, "audit/mancio_silva2022/sa2020_main_analysis_cell_lineage.tsv"),
                   check.names = FALSE, quote = "", stringsAsFactors = FALSE)
lower <- read.delim(file.path(repo, "audit/mancio_silva2022/sa2020_lower_threshold_cell_lineage.tsv"),
                    check.names = FALSE, quote = "", stringsAsFactors = FALSE)

study_ids <- colnames(study)
md <- study@meta.data[study_ids, , drop = FALSE]
raw_barcode <- sub("^[0-9]{3}_", "", study_ids)
run_suffix <- sub("_.*$", "", study_ids)
expected_run <- paste0("SRR11008", run_suffix)
stopifnot(ncol(study) == 9766L, !anyDuplicated(study_ids), identical(as.character(md$run_id), expected_run))

author_prefix <- bridge$author_prefix[match(md$run_id, bridge$run_id)]
stopifnot(!anyNA(author_prefix))
derived_author_id <- paste0(author_prefix, "_", raw_barcode)
main_ids <- unique(main$author_cell_id)
lower_ids <- unique(lower$author_cell_id)

lineage <- data.frame(
  source_cell_identifier = paste(md$run_id, raw_barcode, sep = ":"),
  source_run_id = as.character(md$run_id),
  source_library_prefix = author_prefix,
  source_raw_barcode = raw_barcode,
  author_cell_id = ifelse(derived_author_id %in% union(main_ids, lower_ids), derived_author_id, NA_character_),
  author_main_workbook_member = derived_author_id %in% main_ids,
  author_lower_threshold_workbook_member = derived_author_id %in% lower_ids,
  study_cell_identifier = study_ids,
  final_candidate_cell_identifier = study_ids,
  inclusion_status = "retained_phase2",
  main_analysis_member = as.logical(md$main_analysis_member),
  source_count_object = "STARsolo GeneFull/filtered matrix audited against sa2020.rds",
  reconstruction_evidence_method = "exact run-key plus raw-barcode identity; sa2020_count_lineage.tsv proves exact source-matrix column order/count equality per run",
  final_candidate_present = study_ids %in% final_ids,
  stringsAsFactors = FALSE
)

stopifnot(
  nrow(lineage) == 9766L,
  !anyDuplicated(lineage$source_cell_identifier),
  !anyDuplicated(lineage$study_cell_identifier),
  !anyDuplicated(lineage$final_candidate_cell_identifier),
  all(lineage$final_candidate_present),
  setequal(lineage$final_candidate_cell_identifier, final_ids),
  identical(sort(lineage$final_candidate_cell_identifier), sort(final_ids))
)

write.table(lineage, out_path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
message("Wrote ", out_path, ": 9766 retained cells; 0 missing; 0 duplicates; 0 extras")
