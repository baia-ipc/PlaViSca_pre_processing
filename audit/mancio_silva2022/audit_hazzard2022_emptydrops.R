#!/usr/bin/env Rscript

# Exact audit-only reproduction of the active Hazzard2022 emptyDrops block.
suppressPackageStartupMessages({
  library(DropletUtils)
  library(SingleCellExperiment)
  library(Seurat)
})

audit_dir <- normalizePath("../audit/mancio_silva2022", mustWork = TRUE)
data_dir <- normalizePath("../counts/36525464", mustWork = TRUE)
runs <- c("SRR20710498", "SRR20710499", "SRR20710500", "SRR20710501")
write_tsv <- function(x, name) write.table(x, file.path(audit_dir, name), sep = "\t",
  quote = FALSE, row.names = FALSE, na = "NA")

gff <- read.delim("ref/PlasmoDB-68_PvivaxP01.gff", comment.char = "#", header = FALSE,
                  quote = "", stringsAsFactors = FALSE)
id <- sub(".*(?:^|;)ID=([^;]+).*", "\\1", gff[[9]], perl = TRUE)
rrna <- gsub("\\.1", "", id[gff[[3]] == "rRNA"])
study <- readRDS("../hazzard2022.rds")

set.seed(123456)
summary_rows <- list(); lineage_rows <- list()
for (i in seq_along(runs)) {
  run <- runs[i]
  message("Reading and testing ", run)
  path <- file.path(data_dir, paste0(run, "_solo_out/Solo.out/GeneFull/raw"))
  sce <- read10xCounts(path, col.names = TRUE)
  raw_features_before_rrna <- nrow(sce)
  sce <- sce[!rownames(sce) %in% rrna, ]
  eout <- emptyDrops(sce, lower = 1)
  keep <- which(eout$FDR <= 0.001)
  kept <- colnames(sce)[keep]
  expected <- sub("^[0-9]+_", "", colnames(study)[study$run_id == run])
  filtered_path <- file.path(data_dir, paste0(run, "_solo_out/Solo.out/GeneFull/filtered/barcodes.tsv"))
  star_filtered <- readLines(filtered_path)
  nz <- which(Matrix::colSums(counts(sce)) > 0)
  status <- data.frame(run_id = run, source_barcode = colnames(sce)[nz],
    total = as.numeric(Matrix::colSums(counts(sce)))[nz],
    emptydrops_fdr = eout$FDR[nz], emptydrops_retained = nz %in% keep,
    study_retained = colnames(sce)[nz] %in% expected,
    starsolo_filtered = colnames(sce)[nz] %in% star_filtered,
    stringsAsFactors = FALSE)
  lineage_rows[[i]] <- status
  summary_rows[[i]] <- data.frame(
    run_id = run, raw_features_before_rrna = raw_features_before_rrna,
    raw_features_after_rrna = nrow(sce), raw_barcodes = ncol(sce), nonzero_barcodes = length(nz),
    emptydrops_retained = length(kept), study_cells = length(expected),
    emptydrops_study_intersection = length(intersect(kept, expected)),
    emptydrops_only = length(setdiff(kept, expected)), study_only = length(setdiff(expected, kept)),
    starsolo_filtered = length(star_filtered), emptydrops_starsolo_intersection = length(intersect(kept, star_filtered)),
    emptydrops_not_starsolo_filtered = length(setdiff(kept, star_filtered)),
    starsolo_filtered_not_emptydrops = length(setdiff(star_filtered, kept)),
    stringsAsFactors = FALSE)
  rm(sce, eout, status); gc()
}
write_tsv(do.call(rbind, summary_rows), "hazzard2022_emptydrops_summary.tsv")
write_tsv(do.call(rbind, lineage_rows), "hazzard2022_emptydrops_lineage.tsv")
cat("Exact emptyDrops reproduction complete\n")
