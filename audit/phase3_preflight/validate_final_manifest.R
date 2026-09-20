#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
repo <- normalizePath(getwd())
path <- file.path(repo, "audit/phase3_preflight/final_release_manifest.tsv")
m <- read.delim(path, check.names = FALSE, quote = "", stringsAsFactors = FALSE)
sha256 <- function(path) strsplit(system2("sha256sum", path, stdout = TRUE), " +")[[1]][1]

artifact <- m$path != "."
resolved <- ifelse(grepl("^/", m$path), m$path, file.path(repo, m$path))
exists <- rep(TRUE, nrow(m)); exists[artifact] <- file.exists(resolved[artifact])
hash_match <- rep(TRUE, nrow(m)); hash_match[artifact] <- vapply(resolved[artifact], sha256, character(1)) == m$sha256[artifact]
size_match <- rep(TRUE, nrow(m)); size_match[artifact] <- file.info(resolved[artifact])$size == m$byte_size[artifact]

required_keys <- c(
  "final_preprocessing_code_commit", "phase2_scientific_build_commit", "independent_QA_commit",
  "release_date", "R", "Seurat", "SeuratObject", "harmony", "SingleR", "pixi_lock",
  "integration_seed", "integration_feature_count", "total_cell_count",
  "Mancio-Silva2022_study_object", "Sa2020_study_object", "Hazzard2022_study_object",
  "Hazzard2024_study_object", "Ruberto2022_1_study_object", "Ruberto2022_2_study_object",
  "Ruberto2022_1_author_source", "pv_all_studies", "normalize_df", "raw_df", "scale_df",
  "cleaned_dataset", "data_source"
)

checks <- data.frame(
  check = c("manifest_unique_keys", "required_records_present", "artifact_files_exist",
            "artifact_sha256_match", "artifact_sizes_match", "verified_status_explicit",
            "six_study_counts_present", "all_final_outputs_have_cell_counts"),
  status = c(
    if (!anyDuplicated(m$key)) "PASS" else "FAIL",
    if (all(required_keys %in% m$key)) "PASS" else "FAIL",
    if (all(exists)) "PASS" else "FAIL",
    if (all(hash_match)) "PASS" else "FAIL",
    if (all(size_match)) "PASS" else "FAIL",
    if (all(m$verification_status == "AVAILABLE_AND_VERIFIED")) "PASS" else "FAIL",
    if (sum(m$record_type == "study_count") == 6L) "PASS" else "FAIL",
    if (all(!is.na(m$cell_count[m$record_type %in% c("scientific_output", "application_output")]))) "PASS" else "FAIL"
  ),
  detail = c(
    paste(nrow(m), "records"),
    paste(length(required_keys), "required records"),
    paste(sum(artifact), "artifact paths checked"),
    paste(sum(artifact), "SHA256 values recomputed"),
    paste(sum(artifact), "byte sizes checked"),
    "availability and verification state is explicit",
    "all six exact study counts recorded",
    "atlas and all application outputs carry row/cell counts"
  ), stringsAsFactors = FALSE
)
write.table(checks, file.path(repo, "audit/phase3_preflight/AT23_manifest_validation.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
if (any(checks$status != "PASS")) quit(status = 1L)
message("AT23 PASS: authoritative manifest complete and all hashes verified")
