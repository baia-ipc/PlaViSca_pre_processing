# Diagnostic only. Evaluates inspected study-construction expressions in memory.
# Run from scripts with the Pixi R environment. No production object is written.
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(dplyr))
p <- "../audit/mancio_silva2022"
final <- readRDS("../../PlaViSca/data/cleaned_dataset.rds")$mr_data
final <- final[grepl("Mancio", final$study_label), ]
results <- list()
for (version in c("historical_c80580f", "current")) {
  path <- if(version == "current") "silva.R" else file.path(p, "historical_c80580f_silva.R")
  txt <- gsub("data/", "../data/", readLines(path), fixed = TRUE)
  ee <- new.env()
  for (expr in parse(text = txt)) {
    if (is.call(expr) && as.character(expr[[1]]) %in% c("setwd", "saveRDS")) next
    eval(expr, envir = ee)
  }
  mm <- ee$silva@meta.data
  stopifnot(setequal(rownames(mm), rownames(final)))
  mm <- mm[match(rownames(final), rownames(mm)), ]
  names(mm)[names(mm) == "nCount_RNA"] <- "n_count_rna"
  names(mm)[names(mm) == "nFeature_RNA"] <- "n_feature_rna"
  cols <- intersect(names(mm), names(final))
  results[[version]] <- do.call(rbind, lapply(cols, function(k) {
    x <- as.character(mm[[k]]); y <- as.character(final[[k]])
    eq <- (is.na(x) & is.na(y)) | (!is.na(x) & !is.na(y) & x == y)
    data.frame(version = version, field = k, equal = sum(eq), different = sum(!eq))
  }))
  if (version == "historical_c80580f") {
    stopifnot(all(results[[version]]$different == 0))
    v <- as.numeric(final$hour_post_invasion)
    pred <- ifelse(v < 18, "Ring", ifelse(v < 30, "Trophozoite", ifelse(v < 46, "Schizont", "Merozoite")))
    gam <- final$pred_gametocyte %in% c("Male gametocyte", "Female gametocyte")
    pred[gam] <- final$pred_gametocyte[gam]
    r <- mm$refine_state
    pred[!is.na(r)] <- r[!is.na(r)]
    pred[pred == "Schizont"] <- "Schizont (Liver stage)"
    stopifnot(all(pred == final$life_cycle_stage))
    stage <- data.frame(source_cell_id = rownames(final), historical_refine_state = r,
      final_refine_state = final$refine_state, final_life_cycle_stage = final$life_cycle_stage,
      final_pred_gametocyte = final$pred_gametocyte, reconstructed_stage = pred,
      stage_origin = ifelse(is.na(r), "PlaViSca fallback (IDC or cluster prediction)",
        ifelse(r == "Hypnozoite", "source State", "source Table S2 topcell reference similarity")))
    write.table(stage, file.path(p, "historical_stage_comparison.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  }
}
write.table(do.call(rbind, results), file.path(p, "historical_script_metadata_comparison.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
cat("Historical metadata and stage-construction comparisons passed\n")
