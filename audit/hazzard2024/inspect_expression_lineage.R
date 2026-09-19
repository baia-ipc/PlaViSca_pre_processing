#!/usr/bin/env Rscript

# Exhaustive study-RDS to deployed raw/normalized expression comparison.
# This does not compare PlaViSca counts to the unavailable Zenodo count-table ZIP.

suppressPackageStartupMessages(library(Seurat))

audit_dir <- normalizePath("../audit/hazzard2024", mustWork = TRUE)
pre_dir <- normalizePath("..", mustWork = TRUE)
app_dir <- normalizePath("../../PlaViSca/data", mustWork = TRUE)

write_tsv <- function(x, filename) {
    write.table(
        x,
        file.path(audit_dir, filename),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE,
        na = "NA"
    )
}

study <- readRDS(file.path(pre_dir, "hazzard2024.rds"))
counts <- GetAssayData(study, assay = "RNA", layer = "counts")
counts_by_cell <- Matrix::t(counts)
cell_ids <- colnames(counts)
app_features <- gsub("-", "_", rownames(counts), fixed = TRUE)
stopifnot(!anyDuplicated(app_features))

compare_artifact <- function(filename, mode) {
    artifact <- readRDS(file.path(app_dir, filename))
    idx <- match(cell_ids, rownames(artifact))
    feature_idx <- match(app_features, colnames(artifact))
    stopifnot(!anyNA(idx), !anyNA(feature_idx))

    unequal <- 0
    max_abs_difference <- 0
    compared <- 0
    totals <- Matrix::rowSums(counts_by_cell)
    blocks <- split(seq_len(nrow(counts)),
        ceiling(seq_len(nrow(counts)) / 100L))
    for (block in blocks) {
        observed <- as.matrix(artifact[idx, feature_idx[block], drop = FALSE])
        expected <- as.matrix(counts_by_cell[, block, drop = FALSE])
        if (mode == "normalized") {
            expected <- log1p(sweep(expected, 1L, totals, "/") * 10000)
        }
        differences <- abs(observed - expected)
        unequal <- unequal + sum(differences > if (mode == "raw") 0 else 1e-12)
        max_abs_difference <- max(max_abs_difference, differences)
        compared <- compared + length(expected)
    }
    rm(artifact)
    gc()

    data.frame(
        artifact = filename,
        comparison = mode,
        cells = length(cell_ids),
        source_features = nrow(counts),
        compared_values = compared,
        unequal_values = unequal,
        max_abs_difference = format(max_abs_difference, scientific = TRUE),
        source_cells_missing = sum(is.na(idx)),
        source_features_missing = sum(is.na(feature_idx)),
        stringsAsFactors = FALSE
    )
}

raw_result <- compare_artifact("raw_df.rds", "raw")
normalized_result <- compare_artifact("normalize_df.rds", "normalized")
results <- rbind(raw_result, normalized_result)
stopifnot(all(results$unequal_values == 0))
write_tsv(results, "expression_concordance.tsv")

writeLines(
    c(
        "PASS all 526797992 study-RDS count values equal deployed raw_df values",
        "PASS all 526797992 expected log1p(count/cell_total*10000) values equal deployed normalize_df values within 1e-12",
        "LIMIT this establishes study-RDS to deployed expression lineage, not deposited author count-table to PlaViSca equivalence"
    ),
    file.path(audit_dir, "expression_validation.txt")
)
