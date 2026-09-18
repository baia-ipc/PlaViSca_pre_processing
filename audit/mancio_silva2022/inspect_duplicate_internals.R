# Diagnostic only; run in the scripts Pixi environment. Writes audit tables only.
suppressPackageStartupMessages(library(Seurat))
p <- "../audit/mancio_silva2022"
s <- readRDS("../data/PvData_CHM_Final.RDS")
d <- read.delim(file.path(p, "all_exact_count_duplicates.tsv"), check.names = FALSE)
result <- list()
for (an in names(s@assays)) {
  for (sl in c("counts", "data", "scale.data")) {
    z <- tryCatch(GetAssayData(s, assay = an, layer = sl), error = function(e) NULL)
    if (is.null(z) || !nrow(z)) next
    for (g in unique(d$duplicate_group)) {
      cc <- d$source_cell_id[d$duplicate_group == g]
      if (all(cc %in% colnames(z))) result[[length(result) + 1L]] <- data.frame(group = g, assay = an, layer = sl, features = nrow(z), max_difference = max(abs(z[, cc[1]] - z[, cc[2]])))
    }
  }
}
write.table(do.call(rbind, result), file.path(p, "duplicate_assay_comparison.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
model <- lapply(names(s[["SCT"]]@SCTModel.list), function(n) {
  cc <- rownames(s[["SCT"]]@SCTModel.list[[n]]@cell.attributes)
  data.frame(source_cell_id = d$source_cell_id, model = n, present = d$source_cell_id %in% cc)
})
write.table(do.call(rbind, model), file.path(p, "duplicate_sct_model_membership.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
cat("Completed assay and model duplicate comparisons\n")
