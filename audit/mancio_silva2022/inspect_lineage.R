# Read-only production-data inspection; writes evidence only under this audit directory.
# Run from pre_process_data/scripts with env -u R_LIBS_USER pixi run Rscript --vanilla ../audit/mancio_silva2022/inspect_lineage.R
suppressPackageStartupMessages(library(Seurat))
p <- "../audit/mancio_silva2022"
sink(file.path(p, "focused_validation.txt"), split = TRUE)
write_tsv <- function(x, name) write.table(x, file.path(p, name), sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
s <- readRDS("../data/PvData_CHM_Final.RDS")
m <- s@meta.data
a <- GetAssayData(s, assay = "RNA", layer = "counts")
ids <- colnames(a)
gene <- vapply(strsplit(rownames(a), "-", fixed = TRUE), function(x) paste(head(x, 2), collapse = "_"), character(1))
stopifnot(!anyDuplicated(gene), !anyDuplicated(ids))
hash <- vapply(seq_len(ncol(a)), function(i) digest::digest(as.numeric(a[, i]), algo = "sha256"), character(1))
dup <- duplicated(hash) | duplicated(hash, fromLast = TRUE)
groups <- split(which(dup), hash[dup])
for (g in groups) stopifnot(all(vapply(g, function(i) all(a[, i] == a[, g[1]]), logical(1))))
cat("EXHAUSTIVE RNA COUNT DUPLICATES:", length(groups), "groups;", sum(dup), "cells\n")
cell <- read.delim(file.path(p, "cell_metadata_evidence.tsv"), check.names = FALSE)
cross <- read.delim(file.path(p, "source_group_accession_map.tsv"), check.names = FALSE)
out <- cell[match(ids[dup], cell$source_cell_id), ]
out$duplicate_group <- match(hash[dup], unique(hash[dup]))
out$count_sha256 <- hash[dup]
out$associated_runs <- cross$run_accessions[match(out$orig.ident, cross$orig.ident)]
out$source_column <- which(dup)
out$source_nCount_RNA <- m$nCount_RNA[dup]
out$source_nFeature_RNA <- m$nFeature_RNA[dup]
write_tsv(out, "all_exact_count_duplicates.tsv")
cat("Duplicate group sizes:\n"); print(table(lengths(groups)))
write_tsv(data.frame(source_cell_id = ids, source_position = seq_along(ids), count_sha256 = hash), "source_count_fingerprints.tsv")
cat("SOURCE COMMANDS\n"); print(names(s@commands))
for (n in names(s@commands)) {cat(n, "\n"); print(s@commands[[n]]@params)}

clean <- readRDS("../../PlaViSca/data/cleaned_dataset.rds")$mr_data
ci <- match(ids, rownames(clean)); stopifnot(!anyNA(ci), !anyDuplicated(rownames(clean)))
cm <- clean[ci, , drop = FALSE]
lineage <- data.frame(source_cell_id = ids, source_position = seq_along(ids), cleaned_position = ci)
metadata_results <- list(); expression_results <- list(); schemas <- list()
source_dense <- t(as.matrix(a)); colnames(source_dense) <- gene
expected_norm <- log1p(source_dense / rowSums(source_dense) * 10000)
source_data <- t(as.matrix(GetAssayData(s, assay = "RNA", layer = "data")))
colnames(source_data) <- gene
for (kind in c("raw", "normalize", "scale")) {
  cat("READING", kind, "artifact\n"); flush.console()
  f <- readRDS(paste0("../../PlaViSca/data/", kind, "_df.rds"))
  cat("DIMENSIONS", dim(f), "\n")
  stopifnot(!anyDuplicated(rownames(f)), !anyDuplicated(colnames(f)))
  ix <- match(ids, rownames(f)); stopifnot(!anyNA(ix))
  lineage[[paste0(kind, "_position")]] <- ix
  schemas[[kind]] <- data.frame(artifact = kind, column = colnames(f), position = seq_len(ncol(f)))
  cat("All row set matches cleaned", setequal(rownames(f), rownames(clean)), "; sorted", identical(rownames(f), sort(rownames(f))), "; source order", identical(ix, sort(ix)), "\n")
  small <- f[ix, , drop = FALSE]
  rm(f); gc()
  common_meta <- intersect(names(cm), names(small))
  metadata_results[[kind]] <- do.call(rbind, lapply(common_meta, function(k) {
    x <- cm[[k]]; y <- small[[k]]
    eq <- (is.na(x) & is.na(y)) | (!is.na(x) & !is.na(y) & as.character(x) == as.character(y))
    data.frame(artifact = kind, field = k, equal_cells = sum(eq), different_cells = sum(!eq))
  }))
  gg <- intersect(gene, names(small))
  mat <- as.matrix(small[, gg, drop = FALSE])
  ref <- switch(kind, raw = source_dense[, gg, drop = FALSE], normalize = expected_norm[, gg, drop = FALSE], scale = NULL)
  cat("SOURCE GENES PRESENT", length(gg), "missing", length(setdiff(gene, gg)), "\n")
  if (!is.null(ref)) {
    delta <- abs(mat - ref)
    lineage[[paste0(kind, "_mismatched_genes_exact")]] <- rowSums(delta != 0)
    lineage[[paste0(kind, "_max_absolute_error")]] <- apply(delta, 1, max)
    expression_results[[kind]] <- data.frame(artifact = kind, reference = ifelse(kind == "raw", "source_RNA_counts", "log1p_source_counts_per_10000"), genes = length(gg), unequal_exact = sum(delta != 0), unequal_above_1e_12 = sum(delta > 1e-12), max_absolute_error = max(delta))
    print(expression_results[[kind]])
    if (kind == "normalize") {
      ds <- abs(mat - source_data[, gg, drop = FALSE])
      cat("Normalized vs source RNA data: max", max(ds), "unequal >1e-12", sum(ds > 1e-12), "\n")
    }
  }
  extra <- setdiff(grep("^PVP?01_", names(small), value = TRUE), gene)
  cat("Extra parasite-gene columns", length(extra), "; nonzero entries", if(length(extra)) sum(as.matrix(small[,extra,drop=FALSE]) != 0) else 0, "\n")
  if (kind == "raw") {
    lineage$raw_nCount_matches <- rowSums(mat) == Matrix::colSums(a)
    lineage$raw_nFeature_matches <- rowSums(mat > 0) == Matrix::colSums(a > 0)
  }
  rm(small, mat); gc()
}
write_tsv(lineage, "source_to_final_lineage.tsv")
write_tsv(do.call(rbind, metadata_results), "artifact_metadata_concordance.tsv")
write_tsv(do.call(rbind, expression_results), "expression_concordance.tsv")
write_tsv(do.call(rbind, schemas), "artifact_column_schemas.tsv")
write_tsv(data.frame(source_cell_id = ids, cm, check.names = FALSE), "final_mancio_metadata.tsv")
cat("CLEANED ROW ORDER IS SOURCE ORDER", identical(ci, sort(ci)), "\n")
cat("METADATA DIFFERENCES BETWEEN ARTIFACTS\n"); print(subset(do.call(rbind, metadata_results), different_cells > 0))
print(sessionInfo())
sink()
