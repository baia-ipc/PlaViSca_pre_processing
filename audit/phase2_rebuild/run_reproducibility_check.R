#!/usr/bin/env Rscript
# Phase 2 Part V: reproducibility and stability (AT29, AT30).
# Reads only pv_all_studies.rds (post-singleR.R). Writes only into
# audit/phase2_rebuild/. No production artifact modified.
# Run with: env -u R_LIBS_USER pixi run --manifest-path scripts/pixi.toml \
#   Rscript --vanilla audit/phase2_rebuild/run_reproducibility_check.R

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})

.plavisca_root <- local({
  candidates <- unique(c(
    Sys.getenv("PLAVISCA_PREPROCESS_ROOT", unset = NA_character_),
    getwd(),
    "/home/sopheap/pvsca_b/pre_process_data"
  ))
  candidates <- candidates[!is.na(candidates) & nzchar(candidates)]
  hit <- candidates[file.exists(file.path(candidates, "scripts", "pipeline_lib.R"))]
  if (length(hit) == 0) stop("Could not locate pre_process_data root.", call. = FALSE)
  normalizePath(hit[[1]])
})
setwd(.plavisca_root)
source("scripts/pipeline_lib.R")

so <- readRDS("pv_all_studies.rds")
out_dir <- "audit/phase3_preflight"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

results <- list()

## ------------------------------------------------------------------
## AT29: same-seed determinism. Re-run FindClusters at the SAME
## resolution/seed/graph already stored on the object and confirm the
## cluster assignment is identical to what was saved - this tests that
## the deterministic step of the pipeline (clustering on an already-
## fixed graph/PCA) reproduces exactly, not that the whole stochastic
## pipeline (Harmony/UMAP, which use their own internal seeds) is
## bit-identical end to end.
## ------------------------------------------------------------------
stored_clusters <- so$seurat_clusters
so_rerun <- FindClusters(so, resolution = 0.05, graph.name = "RNA_snn")
rerun_clusters <- so_rerun$seurat_clusters
identical_clusters <- identical(as.character(stored_clusters), as.character(rerun_clusters))
same_seed_ari <- if (requireNamespace("mclust", quietly = TRUE)) {
  mclust::adjustedRandIndex(stored_clusters, rerun_clusters)
} else {
  mean(as.character(stored_clusters) == as.character(rerun_clusters))
}
results$AT29 <- data.frame(
  test = "AT29_same_seed_determinism",
  identical_cluster_labels = identical_clusters,
  agreement_metric = same_seed_ari,
  metric_type = if (requireNamespace("mclust", quietly = TRUE)) "adjusted_rand_index" else "fraction_identical_labels",
  note = "Re-run FindClusters at the same resolution/graph already on the saved object"
)

## ------------------------------------------------------------------
## AT30: seed/resolution stability of biologically important
## populations (marker-validated gametocyte assignments). Re-cluster
## at a perturbed resolution and check whether cells originally called
## Male/Female gametocyte are still assigned to a cluster whose marker
## scores clear the same Female/Male threshold used in singleR.R - the
## question is whether the BIOLOGICAL CALL is stable, not whether every
## UMAP coordinate or numeric cluster ID matches.
## ------------------------------------------------------------------
# AddModuleScore() auto-converts underscore feature names to dashes
# internally before matching (a Seurat quirk, confirmed by the
# "replacing with dashes" warning it emits) - so the marker list passed
# to it must stay in the ORIGINAL dash format (as singleR.R itself uses
# them), not the post-save underscore format rownames(so) is actually in.
femaleGams <- c(
  "PVP01-1207200", "PVP01-0616100", "PVP01-1119300", "PVP01-1465500",
  "PVP01-1259400", "PVP01-1441000", "PVP01-0702600", "PVP01-1306800",
  "PVP01-0946800", "PVP01-0517400", "PVP01-1027600", "PVP01-1003000",
  "PVP01-1143200", "PVP01-1017500", "PVP01-1024300", "PVP01-0806000",
  "PVP01-1240300", "PVP01-0603600", "PVP01-0712800"
)
maleGams <- c(
  "PVP01-1262200", "PVP01-0530800", "PVP01-1412100", "PVP01-1025600",
  "PVP01-1266500", "PVP01-1229400"
)
rownames_dash <- gsub("_", "-", rownames(so))
femaleGams <- intersect(femaleGams, rownames_dash)
maleGams <- intersect(maleGams, rownames_dash)

perturbed_res <- 0.08 # a different, still-documented resolution (vs 0.05 used in production)
so_perturbed <- FindClusters(so, resolution = perturbed_res, graph.name = "RNA_snn")
so_perturbed <- AddModuleScore(so_perturbed, features = list(femaleGams), name = "femaleGams_perturbed")
so_perturbed <- AddModuleScore(so_perturbed, features = list(maleGams), name = "maleGams_perturbed")

original_gametocytes <- which(!is.na(so$pred_gametocyte_sex))
n_original_gametocytes <- length(original_gametocytes)

if (n_original_gametocytes > 0) {
  perturbed_cluster_summary <- so_perturbed@meta.data %>%
    group_by(seurat_clusters) %>%
    summarise(
      mean_female = mean(femaleGams_perturbed1, na.rm = TRUE),
      mean_male = mean(maleGams_perturbed1, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(call = case_when(
      mean_female >= 0.10 & (mean_female - mean_male) >= 0.05 ~ "Female gametocyte",
      mean_male >= 0.10 & (mean_male - mean_female) >= 0.05 ~ "Male gametocyte",
      TRUE ~ "Asexual"
    ))
  cluster_to_label_perturbed <- setNames(perturbed_cluster_summary$call, perturbed_cluster_summary$seurat_clusters)
  perturbed_call_for_original <- unname(cluster_to_label_perturbed[as.character(so_perturbed$seurat_clusters[original_gametocytes])])
  original_call <- so$pred_gametocyte_sex[original_gametocytes]
  n_stable <- sum(perturbed_call_for_original == original_call, na.rm = TRUE)
  stability_fraction <- n_stable / n_original_gametocytes
} else {
  stability_fraction <- NA_real_
  n_stable <- 0L
}

results$AT30 <- data.frame(
  test = "AT30_gametocyte_call_stability_under_perturbation",
  perturbed_resolution = perturbed_res,
  n_original_gametocyte_cells = n_original_gametocytes,
  n_stable_under_perturbation = n_stable,
  stability_fraction = stability_fraction,
  note = "Fraction of originally-called Male/Female gametocyte cells whose new (perturbed-resolution) cluster still clears the same marker-score threshold for the SAME sex call"
)

reproducibility_tbl <- bind_rows(results)
write.table(reproducibility_tbl, file.path(out_dir, "reproducibility_check.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

cat("=== Reproducibility/stability check ===\n")
print(reproducibility_tbl)
cat("\nWrote audit/phase3_preflight/reproducibility_check.tsv\n")
