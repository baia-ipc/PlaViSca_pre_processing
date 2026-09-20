#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
candidate_root <- Sys.getenv("PLAVISCA_CANDIDATE_ROOT", "/home/baia/prj/plavisca/pre_process_data")
csv_path <- file.path(candidate_root, "data/data_source.csv")
meta <- readRDS(file.path(candidate_root, "data/candidate_export/cleaned_dataset.rds"))$mr_data
bib <- read.csv(csv_path, check.names = FALSE)

pairs <- unique(meta[c("study_label", "study_pmid")])
stopifnot(nrow(pairs) == 6L, !anyDuplicated(pairs$study_label), !anyNA(pairs$study_pmid))
idx <- match(bib$Study_label, pairs$study_label)
stopifnot(nrow(bib) == 6L, !anyDuplicated(bib$Study_label), !anyNA(idx))
bib$study_pmid <- as.character(pairs$study_pmid[idx])

# Keep the production script's documented column order.
bib <- bib[c("Publication_date", "Study_label", "study_pmid", "Authors", "Journal", "Title", "DOI", "Number_of_cells")]
tmp <- tempfile("data_source_", tmpdir = dirname(csv_path), fileext = ".csv")
on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
write.csv(bib, tmp, row.names = FALSE, na = "")
check <- read.csv(tmp, check.names = FALSE)
stopifnot(identical(as.character(check$study_pmid), as.character(pairs$study_pmid[match(check$Study_label, pairs$study_label)])))
if (!file.rename(tmp, csv_path)) stop("Atomic replacement of data_source.csv failed")
message("Added verified study_pmid for all six bibliography rows")
