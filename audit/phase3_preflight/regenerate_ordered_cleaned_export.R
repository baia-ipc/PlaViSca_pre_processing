#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

candidate_root <- Sys.getenv("PLAVISCA_CANDIDATE_ROOT", "/home/baia/prj/plavisca/pre_process_data")
candidate_dir <- file.path(candidate_root, "data/candidate_export")
cleaned_path <- file.path(candidate_dir, "cleaned_dataset.rds")
expression_paths <- c(
  raw = file.path(candidate_dir, "raw_df.rds"),
  normalized = file.path(candidate_dir, "normalize_df.rds"),
  scaled = file.path(candidate_dir, "scale_df.rds")
)

cleaned <- readRDS(cleaned_path)
old_mr <- cleaned$mr_data
canonical_ids <- sort(rownames(old_mr))
stopifnot(length(canonical_ids) == 105963L, !anyDuplicated(canonical_ids))

for (nm in names(expression_paths)) {
  x <- readRDS(expression_paths[[nm]])
  stopifnot(nrow(x) == 105963L, !anyDuplicated(rownames(x)), identical(rownames(x), canonical_ids))
  rm(x); gc()
}

non_cell_names <- setdiff(names(cleaned), "mr_data")
non_cell_before <- cleaned[non_cell_names]
cleaned$mr_data <- old_mr[canonical_ids, , drop = FALSE]
stopifnot(
  identical(rownames(cleaned$mr_data), canonical_ids),
  identical(cleaned$mr_data, old_mr[canonical_ids, , drop = FALSE]),
  identical(cleaned[non_cell_names], non_cell_before)
)

tmp <- tempfile("cleaned_dataset_ordered_", tmpdir = candidate_dir, fileext = ".rds")
on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
saveRDS(cleaned, tmp)
roundtrip <- readRDS(tmp)
stopifnot(identical(roundtrip, cleaned))
if (!file.rename(tmp, cleaned_path)) stop("Atomic replacement of cleaned_dataset.rds failed")
message("Rewrote only cleaned_dataset.rds with canonical lexicographic cell-key order; non-cell tables unchanged")
