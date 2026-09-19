library(DropletUtils)
library(tidyverse)
library(rtracklayer)
library(Seurat)
library(scater)
library(SingleCellExperiment)

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

rds_files <- list.files(pattern = "^hazzard2024_pv.*\\.rds$")

hazzard2024_list <- lapply(rds_files, readRDS)

# Merge all
pv.combined.all <- Reduce(function(x, y) merge(x, y = y), hazzard2024_list)

# Join layers if needed (multi-layer data)
pv.combined.all <- JoinLayers(pv.combined.all)

# read preocess data from study hazzard2024
df <- read.delim("data/Proccessed_Data.txt") |>
  rownames_to_column("barcode") |>
  mutate(
    prefix = sub("(_[^_]+)$", "", barcode),
    cells = sub("^[^_]+_[^_]+_", "", barcode)
  ) |>
  select(prefix, cells) |>
  mutate(
    prefix = case_when(
      prefix == "4881_1" ~ "994",
      prefix == "5142_1" ~ "992",
      prefix == "5142_2" ~ "991",
      prefix == "5163_1" ~ "981",
      prefix == "5163_2" ~ "980",
      prefix == "5164_1" ~ "976",
      prefix == "5164_2" ~ "975",
      prefix == "5309_1" ~ "956",
      prefix == "5309_2" ~ "955",
      prefix == "5350_1" ~ "972",
      prefix == "5350_2" ~ "970",
      prefix == "5350_3" ~ "969",
      prefix == "5370_1" ~ "968",
      prefix == "5370_2" ~ "967",
      prefix == "5370_3" ~ "966",
      prefix == "5537_1" ~ "987",
      prefix == "5537_2" ~ "986",
      prefix == "5537_3" ~ "985",
      prefix == "5537_4" ~ "984",
      prefix == "5550_1" ~ "983",
      prefix == "5708_1" ~ "993",
      prefix == "5708_2" ~ "982",
      prefix == "5708_3" ~ "971",
      prefix == "5708_4" ~ "960",
    ),
    new_cells = paste0(prefix, "_", cells)
  )

assert_no_missing_mandatory(df$new_cells, label = "df$new_cells (retained-cell crosswalk)")
assert_unique_cell_ids(df$new_cells, label = "df$new_cells (retained-cell crosswalk)")

pv.combined.all <- subset(
  pv.combined.all,
  cells = df$new_cells
)

# Conservative cell-membership policy for Phase 1: preserve the exact
# 80,024-cell Hazzard2024 inclusion population.
assert_cell_count(ncol(pv.combined.all), 80024L, label = "hazzard2024.rds")
assert_unique_cell_ids(colnames(pv.combined.all), label = "hazzard2024.rds colnames")
assert_no_literal_na_string(pv.combined.all$library_id, label = "hazzard2024$library_id")
assert_no_literal_na_string(pv.combined.all$day_post_infection, label = "hazzard2024$day_post_infection")

# Save the final object
saveRDS(pv.combined.all, file = "hazzard2024.rds")
record_build_manifest(
  artifact_path = "hazzard2024.rds",
  script_path = "scripts/hazzard2024_merge_all.R",
  cell_count = ncol(pv.combined.all),
  notes = "Fixes D026-D030,D032,D033,D069 (keyed run table); preserves exact 80,024-cell population per Phase 1 conservative membership policy"
)

# Cleanup
rm(list = ls())
gc()
