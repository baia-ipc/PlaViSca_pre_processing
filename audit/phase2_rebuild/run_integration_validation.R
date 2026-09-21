#!/usr/bin/env Rscript
# Phase 2 Part III/VII: integration validation (AT24-AT27).
# Reads only pv_all_studies.rds (post-singleR.R, fully annotated).
# Writes only into audit/phase2_rebuild/. No production artifact modified.
# Run with: env -u R_LIBS_USER pixi run --manifest-path scripts/pixi.toml \
#   Rscript --vanilla audit/phase2_rebuild/run_integration_validation.R
# (from pre_process_data/)

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
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
meta <- so@meta.data

out_dir <- "audit/phase3_preflight"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## ------------------------------------------------------------------
## Shared neighbor-composition helper: for a given PCA embedding, find
## each cell's k nearest neighbors and summarize the study/stage
## composition of its neighborhood. A simple, transparent kNN-based
## mixing metric (no extra package required beyond what's already
## installed) - not kBET/LISI, but the same underlying idea: how
## homogeneous vs mixed is a cell's local neighborhood.
## ------------------------------------------------------------------
knn_composition <- function(embedding, k = 30, group_vec, sample_n = 3000, seed = 123) {
  set.seed(seed)
  n <- nrow(embedding)
  idx <- if (n > sample_n) sample(seq_len(n), sample_n) else seq_len(n)
  nn <- RANN::nn2(embedding, embedding[idx, , drop = FALSE], k = k + 1)$nn.idx[, -1, drop = FALSE]
  # for each sampled cell, fraction of its k neighbors that share its own group
  own_group <- group_vec[idx]
  frac_same_group <- vapply(seq_along(idx), function(i) {
    neighbor_groups <- group_vec[nn[i, ]]
    mean(neighbor_groups == own_group[i], na.rm = TRUE)
  }, numeric(1))
  data.frame(cell_idx = idx, group = own_group, frac_same_group = frac_same_group)
}

## ------------------------------------------------------------------
## AT24/AT26/AT27: per-study neighborhood composition, unintegrated vs
## integrated PCA. A well-mixed integration should show a LOWER
## frac_same_study in integrated space than unintegrated space, for
## biologically comparable populations - NOT for populations that are
## legitimately distinct (e.g. liver vs blood vs sporozoite).
## ------------------------------------------------------------------
pca_u <- Embeddings(so, "pca_unintegrated")[, 1:30]
pca_i <- Embeddings(so, "pca_integrated")[, 1:30]
study_vec <- meta$study_label
stage_vec <- meta$parasite_broad_stage

comp_u_study <- knn_composition(pca_u, group_vec = study_vec)
comp_i_study <- knn_composition(pca_i, group_vec = study_vec)

per_study_mixing <- comp_u_study %>%
  mutate(space = "unintegrated") %>%
  bind_rows(comp_i_study %>% mutate(space = "integrated")) %>%
  group_by(group, space) %>%
  summarise(mean_frac_same_study = mean(frac_same_group), n = n(), .groups = "drop") %>%
  pivot_wider(names_from = space, values_from = c(mean_frac_same_study, n)) %>%
  rename(study = group)

write.table(per_study_mixing, file.path(out_dir, "per_study_neighborhood_composition.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

## AT26: batch mixing evaluated ONLY within comparable biological
## populations (same broad stage) - mixing liver+sporozoite+blood
## together is not "batch effect", it is real biology, so this is
## computed separately, stratified by parasite_broad_stage.
mixing_by_stage <- lapply(unique(stage_vec[!is.na(stage_vec)]), function(st) {
  cells_in_stage <- which(stage_vec == st)
  if (length(cells_in_stage) < 50) {
    return(data.frame(broad_stage = st, note = "too few cells to evaluate (<50)"))
  }
  sub_pca_u <- pca_u[cells_in_stage, , drop = FALSE]
  sub_pca_i <- pca_i[cells_in_stage, , drop = FALSE]
  sub_study <- study_vec[cells_in_stage]
  n_studies_in_stage <- length(unique(sub_study))
  if (n_studies_in_stage < 2) {
    return(data.frame(
      broad_stage = st, n_cells = length(cells_in_stage),
      n_studies_present = n_studies_in_stage,
      note = "single study - mixing not applicable"
    ))
  }
  cu <- knn_composition(sub_pca_u, group_vec = sub_study, sample_n = min(2000, length(cells_in_stage)))
  ci <- knn_composition(sub_pca_i, group_vec = sub_study, sample_n = min(2000, length(cells_in_stage)))
  data.frame(
    broad_stage = st, n_cells = length(cells_in_stage), n_studies_present = n_studies_in_stage,
    mean_frac_same_study_unintegrated = mean(cu$frac_same_group),
    mean_frac_same_study_integrated = mean(ci$frac_same_group),
    note = ""
  )
})
mixing_by_stage_df <- bind_rows(mixing_by_stage)
write.table(mixing_by_stage_df, file.path(out_dir, "batch_mixing_within_comparable_populations.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

## AT27: minority-study retention - does the dominant study (Hazzard2024,
## ~75% of cells) erase minority-study signal? Check whether each
## minority study still forms its own identifiable neighborhood
## structure (i.e. is not maximally diluted to background composition)
## post-integration, by comparing each study's mean_frac_same_study
## in integrated space against the null expectation if cells mixed
## uniformly at random (proportional to that study's overall share).
study_share <- prop.table(table(study_vec))
per_study_mixing$expected_frac_if_random <- unname(study_share[per_study_mixing$study])
per_study_mixing$retention_ratio_integrated <- per_study_mixing$mean_frac_same_study_integrated / per_study_mixing$expected_frac_if_random
write.table(per_study_mixing, file.path(out_dir, "per_study_neighborhood_composition.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

## AT25: marker gene preservation - known lifecycle marker genes should
## still show broad_stage-specific expression patterns post-integration
## (checked on RNA expression directly, not on Harmony-space coordinates
## - AT24's contract that "RNA-expression integrity is kept explicitly
## separate from Harmony-space").
marker_genes <- c(
  "PVP01_1207200", "PVP01_0616100", # femaleGams (subset)
  "PVP01_1262200", "PVP01_0530800" # maleGams (subset)
)
marker_genes <- intersect(marker_genes, rownames(so))
if (length(marker_genes) > 0) {
  expr <- FetchData(so, vars = marker_genes, layer = "data")
  expr$broad_stage <- stage_vec
  marker_preservation <- expr %>%
    filter(!is.na(broad_stage)) %>%
    group_by(broad_stage) %>%
    summarise(across(all_of(marker_genes), ~ mean(.x, na.rm = TRUE)), n = n(), .groups = "drop")
  write.table(marker_preservation, file.path(out_dir, "marker_gene_preservation.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
} else {
  write.table(data.frame(note = "no marker genes found in this build's feature panel"),
    file.path(out_dir, "marker_gene_preservation.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
}

## ------------------------------------------------------------------
## Combined integration_validation.tsv summary (what the test suite
## checks for existence/non-emptiness).
## ------------------------------------------------------------------
summary_tbl <- per_study_mixing %>%
  transmute(
    metric = "per_study_neighborhood_composition", study,
    unintegrated = mean_frac_same_study_unintegrated,
    integrated = mean_frac_same_study_integrated,
    retention_ratio_integrated
  )
write.table(summary_tbl, file.path(out_dir, "integration_validation.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

cat("=== per-study neighborhood composition (mean fraction of kNN from same study) ===\n")
print(per_study_mixing)
cat("\n=== batch mixing within comparable (same broad-stage) populations ===\n")
print(mixing_by_stage_df)
cat("\nWrote integration_validation.tsv, per_study_neighborhood_composition.tsv, batch_mixing_within_comparable_populations.tsv, marker_gene_preservation.tsv\n")
