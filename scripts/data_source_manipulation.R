library(readr)
library(dplyr)
library(stringr)

.plavisca_root <- local({
  candidates <- unique(c(
    Sys.getenv("PLAVISCA_PREPROCESS_ROOT", unset = NA_character_),
    getwd(),
    "/home/sopheap/pvsca_b/pre_process_data"
  ))
  candidates <- candidates[!is.na(candidates) & nzchar(candidates)]
  hit <- candidates[file.exists(file.path(candidates, "scripts", "pipeline_lib.R"))]
  if (length(hit) == 0) {
    stop(
      "Could not locate the pre_process_data project root (looked for scripts/",
      "pipeline_lib.R under $PLAVISCA_PREPROCESS_ROOT, the current working ",
      "directory, and the historical hard-coded path). Set the ",
      "PLAVISCA_PREPROCESS_ROOT environment variable to the pre_process_data ",
      "directory's absolute path, or run this script from that directory.",
      call. = FALSE
    )
  }
  normalizePath(hit[[1]])
})
setwd(.plavisca_root)
source("scripts/pipeline_lib.R")

# Publication dates
# D057: Hazzard2022 (row 4) corrected from "16-Nov-22" to the actual
# publication date "16-Dec-22".
Publication_date <- c(
  "4-May-20",
  "4-Aug-22",
  "25-Aug-22",
  "16-Dec-22",
  "2-Sep-24",
  "19-Apr-22"
)

# Study labels: sourced from the single shared STUDY_LABELS constant
# (scripts/pipeline_lib.R) so this table can never again drift from the
# per-study scripts/singleR.R/app spelling (fixes D056: "Sar2020" typo and
# the "Mancio-Silva2022" vs "Mancio Silva2022" mismatch).
Study_label <- c(
  STUDY_LABELS[["sa2020"]],
  STUDY_LABELS[["ruberto2022_2"]],
  STUDY_LABELS[["ruberto2022_1"]],
  STUDY_LABELS[["hazzard2022"]],
  STUDY_LABELS[["hazzard2024"]],
  STUDY_LABELS[["mancio_silva2022"]]
)

# D067 release closure: authoritative PMIDs are the exact per-study constants
# already carried by every atlas cell. Keeping the key in the bibliography
# permits a direct startup foreign-key check instead of a label-only warning.
study_pmid <- c(
  "32365102", # Sa2020
  "35926062", # Ruberto2022_2
  "36093191", # Ruberto2022_1
  "36525464", # Hazzard2022
  "39223117", # Hazzard2024
  "35443155"  # Mancio-Silva2022
)

# Authors
Authors <- c(
  "Sà JM et al.",
  "Ruberto AA et al.",
  "Ruberto AA et al.",
  "Hazzard B et al.",
  "Hazzard B et al.",
  "Mancio-Silva L et al."
)

# Journals
Journal <- c(
  "PLoS Biol",
  "PLoS Negl Trop Dis",
  "Front. Cell. Infect. Microbiol.",
  "PLoS Negl Trop Dis",
  "Nat Commun",
  "Cell Host & Microbe"
)

# Titles
Title <- c(
  "Single-cell transcription analysis of <em>Plasmodium vivax</em> blood-stage parasites identifies stage- and species-specific profiles of expression",
  "Single-cell RNA sequencing of <em>Plasmodium vivax</em> sporozoites reveals stage- and species-specific transcriptomic signatures",
  "Single-cell RNA profiling of <em>Plasmodium vivax</em>-infected hepatocytes reveals parasite- and host- specific transcriptomic signatures and therapeutic targets",
  "Long read single cell RNA sequencing reveals the isoform diversity of <em>Plasmodium vivax</em> transcripts",
  "Single-cell analyses of polyclonal <em>Plasmodium vivax</em> infections and their consequences on parasite transmission",
  "A single-cell liver atlas of <em>Plasmodium vivax</em> infection"
)

# DOIs with links
DOI <- c(
  "<a href='https://doi.org/10.1371/journal.pbio.3000711' target='_blank'>https://doi.org/10.1371/journal.pbio.3000711</a>",
  "<a href='https://doi.org/10.1371/journal.pntd.0010633' target='_blank'>https://doi.org/10.1371/journal.pntd.0010633</a>",
  "<a href='https://doi.org/10.3389/fcimb.2022.986314' target='_blank'>https://doi.org/10.3389/fcimb.2022.986314</a>",
  "<a href='https://doi.org/10.1371/journal.pntd.0010991' target='_blank'>https://doi.org/10.1371/journal.pntd.0010991</a>",
  "<a href='https://doi.org/10.1038/s41467-024-51949-8' target='_blank'>https://doi.org/10.1038/s41467-024-51949-8</a>",
  "<a href='https://doi.org/10.1016/j.chom.2022.03.034' target='_blank'>https://doi.org/10.1016/j.chom.2022.03.034</a>"
)

# Number of cells
Number_of_cells <- c(
  9766,
  9947,
  1438,
  3294,
  80024,
  1494
)

# D067: keyed, fail-loud construction (never a bare cbind of independently
# maintained vectors) - all seven vectors must have exactly 6 entries, one
# per study, in the same order.
n_studies <- length(STUDY_LABELS)
for (v in list(Publication_date, Study_label, study_pmid, Authors, Journal, Title, DOI, Number_of_cells)) {
  assert_cardinality(v, n_studies, label = "data_source_manipulation.R study vector")
}
assert_unique_cell_ids(Study_label, label = "data_source.csv Study_label")
unexpected_labels <- setdiff(Study_label, STUDY_LABELS)
if (length(unexpected_labels) > 0) {
  pipeline_fail(paste0(
    "data_source.csv Study_label value(s) not in the canonical STUDY_LABELS constant: ",
    paste(unexpected_labels, collapse = ", ")
  ))
}

df <- data.frame(
  Publication_date = Publication_date,
  Study_label = Study_label,
  study_pmid = study_pmid,
  Authors = Authors,
  Journal = Journal,
  Title = Title,
  DOI = DOI,
  Number_of_cells = Number_of_cells,
  stringsAsFactors = FALSE
)

# Ensure the folder exists, then write to CSV

write_csv(df, "data/data_source.csv")
record_build_manifest(
  artifact_path = "data/data_source.csv",
  script_path = "scripts/data_source_manipulation.R",
  notes = "Fixes D056 (Study_label sourced from shared STUDY_LABELS constant),D057 (Hazzard2022 publication date)"
)

# Print the data frame
print(df)
