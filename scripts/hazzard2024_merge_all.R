library(DropletUtils)
library(tidyverse)
library(rtracklayer)
library(Seurat)
library(scater)
library(SingleCellExperiment)

setwd("/home/sopheap/pvsca_b/pre_process_data")

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

pv.combined.all <- subset(
  pv.combined.all,
  cells = df$new_cells
)
# Save the final object
saveRDS(pv.combined.all, file = "hazzard2024.rds")

# Cleanup
rm(list = ls())
gc()
