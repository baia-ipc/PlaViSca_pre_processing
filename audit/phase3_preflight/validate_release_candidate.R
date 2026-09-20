#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
repo <- normalizePath(getwd())
candidate_root <- Sys.getenv("PLAVISCA_CANDIDATE_ROOT", "/home/baia/prj/plavisca/pre_process_data")
candidate_dir <- file.path(candidate_root, "data/candidate_export")
out_dir <- file.path(repo, "audit/phase3_preflight")

write_tsv <- function(x, name) write.table(x, file.path(out_dir, name), sep = "\t", quote = FALSE,
                                            row.names = FALSE, na = "NA")
cleaned <- readRDS(file.path(candidate_dir, "cleaned_dataset.rds"))
meta <- cleaned$mr_data
canonical <- sort(rownames(meta))

paths <- c(raw = "raw_df.rds", normalized = "normalize_df.rds", scaled = "scale_df.rds")
order_rows <- lapply(names(paths), function(nm) {
  x <- readRDS(file.path(candidate_dir, paths[[nm]]))
  ids <- rownames(x)
  out <- data.frame(
    artifact = paths[[nm]], rows = nrow(x), unique_keys = !anyDuplicated(ids),
    identical_set_to_mr_data = setequal(ids, rownames(meta)),
    identical_order_to_mr_data = identical(ids, rownames(meta)),
    canonical_lexicographic_order = identical(ids, canonical),
    status = if (nrow(x) == 105963L && !anyDuplicated(ids) && identical(ids, rownames(meta)) && identical(ids, canonical)) "PASS" else "FAIL"
  )
  rm(x); gc(); out
})
order_rows[[length(order_rows) + 1L]] <- data.frame(
  artifact = "cleaned_dataset.rds$mr_data", rows = nrow(meta), unique_keys = !anyDuplicated(rownames(meta)),
  identical_set_to_mr_data = TRUE, identical_order_to_mr_data = TRUE,
  canonical_lexicographic_order = identical(rownames(meta), canonical),
  status = if (nrow(meta) == 105963L && !anyDuplicated(rownames(meta)) && identical(rownames(meta), canonical)) "PASS" else "FAIL"
)
order_out <- do.call(rbind, order_rows)
write_tsv(order_out, "AT17_cell_order_validation.tsv")

missing_broad <- is.na(meta$parasite_broad_stage)
missing_detail <- is.na(meta$harmonized_life_cycle_stage)
missing <- meta[missing_broad & missing_detail, , drop = FALSE]
stage_out <- data.frame(
  cell_id = rownames(missing), study_label = missing$study_label,
  total_umi_count = missing$total_umi_count,
  near_empty_expression_flag = missing$near_empty_expression_flag,
  dot1_duplicate_flag = missing$dot1_duplicate_flag,
  source_life_cycle_stage = missing$source_life_cycle_stage,
  idc_reference_similarity_label = missing$idc_reference_similarity_label,
  documented_class = ifelse(
    missing$study_label == "Hazzard2024" & missing$total_umi_count < 10,
    "below_10_UMI_inference_threshold",
    ifelse(missing$study_label == "Mancio-Silva2022" & missing$dot1_duplicate_flag %in% TRUE,
           "retained_flagged_dot1_duplicate_without_stage", "UNEXPECTED")
  )
)
write_tsv(stage_out, "stage_missing_validation.tsv")

stopifnot(
  all(order_out$status == "PASS"),
  sum(missing_broad != missing_detail) == 0L,
  nrow(missing) == 21L,
  sum(stage_out$documented_class == "below_10_UMI_inference_threshold") == 7L,
  sum(stage_out$documented_class == "retained_flagged_dot1_duplicate_without_stage") == 14L,
  !any(stage_out$documented_class == "UNEXPECTED")
)
message("Release candidate validation PASS: canonical order and 21 documented paired-stage missing cells")
