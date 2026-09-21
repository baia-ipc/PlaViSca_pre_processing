# Read-only upstream duplicate tracing. Writes evidence under the audit directory.
# Run in scripts Pixi environment, with one extracted GEO filename prefix argument.
suppressPackageStartupMessages(library(Seurat))
p <- "../audit/mancio_silva2022"
prefix <- commandArgs(trailingOnly = TRUE)[1]
stopifnot(!is.na(prefix), !grepl("/", prefix))
sink(file.path(p, paste0(prefix, "_trace.txt")), split = TRUE)
s <- readRDS("../data/PvData_CHM_Final.RDS")
a <- GetAssayData(s, assay = "RNA", layer = "counts")
d <- read.delim(file.path(p, "all_exact_count_duplicates.tsv"), check.names = FALSE)
d <- d[d$Sample == if(grepl("Infection1", prefix)) "Infection1" else "Infection2", ]
cell_file <- file.path(p, paste0(prefix, "_CellIDs.txt.gz"))
if (!file.exists(cell_file)) cell_file <- file.path(p, paste0(prefix, "_CellIDs.1.txt.gz"))
cells <- read.table(gzfile(cell_file), header = TRUE)[[1]]
genes <- read.table(gzfile(file.path(p, paste0(prefix, "_genes.txt.gz"))), header = TRUE)[[1]]
key <- sub("[.]1$", "", d$source_cell_id)
# Match exact IDs first. Report barcode-token candidates too; do not assume they belong to the same array.
token <- sub("^[^_]+_", "", key)
candidate <- which(cells %in% key | sub("^[^_]+_", "", cells) %in% token)
cat("UPSTREAM", prefix, "cells", length(cells), "unique", length(unique(cells)), "genes", length(genes), "\n")
cat("Candidate columns:"); print(data.frame(column = candidate, cell = cells[candidate]))
if (!length(candidate)) {
  absent <- data.frame(source_cell_id = d$source_cell_id,
    expected_upstream_cell_id = key, source_barcode_token = token,
    exact_identifier_matches = vapply(key, function(k) sum(cells == k), integer(1)),
    barcode_token_matches = vapply(token, function(k) sum(sub("^[^_]+_", "", cells) == k), integer(1)),
    upstream_total_columns = length(cells),
    status = "No corresponding called barcode column in this deposited matrix")
  print(absent)
  write.table(absent, file.path(p, paste0(prefix, "_comparison.tsv")), sep = "\t", row.names = FALSE, quote = FALSE)
  cat("No count-value comparison possible without an identified upstream column; matrix entries were not scanned.\n")
  sink()
  quit(save = "no", status = 0)
}
con <- gzfile(file.path(p, paste0(prefix, "_counts.mtx.gz")), "rt")
header <- readLines(con, n = 1)
while (startsWith(header, "%")) header <- readLines(con, n = 1)
dims <- scan(text = header, quiet = TRUE)
stopifnot(dims[1] == length(genes), dims[2] == length(cells))
pattern <- paste0("^[0-9]+[[:space:]]+(", paste(candidate, collapse = "|"), ")[[:space:]]")
entries <- list(); total <- 0L
repeat {
  lines <- readLines(con, n = 200000L)
  if (!length(lines)) break
  selected <- grep(pattern, lines, value = TRUE, perl = TRUE)
  if (length(selected)) entries[[length(entries) + 1L]] <- matrix(scan(text = paste(selected, collapse = "\n"), quiet = TRUE), ncol = 3, byrow = TRUE)
  total <- total + length(lines)
  if (total %% 20000000L == 0L) {cat("Scanned", total, "of", dims[3], "entries\n"); flush.console()}
}
close(con)
stopifnot(total == dims[3])
e <- do.call(rbind, entries)
z <- Matrix::sparseMatrix(i = e[,1], j = match(e[,2], candidate), x = e[,3], dims = c(length(genes), length(candidate)), dimnames = list(genes, cells[candidate]))
ri <- match(rownames(a), genes)
cat("COMPARED GENES", sum(!is.na(ri)), "; source counts in absent genes", sum(a[is.na(ri), d$source_cell_id]), "\n")
out <- do.call(rbind, lapply(seq_len(nrow(d)), function(i) {
  # Report all same-token candidates, including competing groups.
  hit <- which(sub("^[^_]+_", "", colnames(z)) == token[i])
  do.call(rbind, lapply(hit, function(j) data.frame(
    source_cell_id = d$source_cell_id[i], upstream_cell_id = colnames(z)[j],
    exact_group_key_match = colnames(z)[j] == key[i], upstream_column = candidate[j],
    upstream_columns_for_exact_key = sum(cells == key[i]),
    compared_genes = sum(!is.na(ri)),
    source_counts_in_missing_genes = sum(a[is.na(ri), d$source_cell_id[i]]),
    different_genes = sum(z[ri[!is.na(ri)], j] != a[!is.na(ri), d$source_cell_id[i]]),
    upstream_parasite_total = sum(z[grepl("^PVP01", rownames(z)), j]),
    upstream_all_species_total = sum(z[, j])
  )))
}))
print(out)
write.table(out, file.path(p, paste0(prefix, "_comparison.tsv")), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(data.frame(gene = genes, as.matrix(z), check.names = FALSE), file.path(p, paste0(prefix, "_selected_columns.tsv")), sep = "\t", row.names = FALSE, quote = FALSE)
cat("All-species candidate pair differences\n")
if(ncol(z) > 1) for(i in seq_len(ncol(z)-1)) for(j in seq.int(i+1,ncol(z))) cat(colnames(z)[i],colnames(z)[j],sum(z[,i]!=z[,j]),"\n")
sink()
