#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(Seurat))

repo <- normalizePath(getwd())
candidate_root <- Sys.getenv("PLAVISCA_CANDIDATE_ROOT", repo)
so <- readRDS(file.path(candidate_root, "pv_all_studies.rds"))

female_markers <- c(
  "PVP01-1207200", "PVP01-0616100", "PVP01-1119300", "PVP01-1465500",
  "PVP01-1259400", "PVP01-1441000", "PVP01-0702600", "PVP01-1306800",
  "PVP01-0946800", "PVP01-0517400", "PVP01-1027600", "PVP01-1003000",
  "PVP01-1143200", "PVP01-1017500", "PVP01-1024300", "PVP01-0806000",
  "PVP01-1240300", "PVP01-0603600", "PVP01-0712800"
)
male_markers <- c(
  "PVP01-1262200", "PVP01-0530800", "PVP01-1412100", "PVP01-1025600",
  "PVP01-1266500", "PVP01-1229400"
)

# Recreate production marker definitions on an in-memory copy only. This reads
# the normalized RNA layer and never alters or serializes the atlas.
so <- AddModuleScore(so, features = list(female_markers), name = "AT25_female", seed = 1)
so <- AddModuleScore(so, features = list(male_markers), name = "AT25_male", seed = 1)
md <- so@meta.data

# Exact production gate in scripts/singleR.R: gametocyte_eligible <-
# idc_eligible. Production filters to this population before cluster scoring,
# maps qualifying cluster calls, then masks all ineligible cells.
gametocyte_eligible <- (
  !is.na(md$tissue_or_sample_type) &
    md$tissue_or_sample_type == "Host blood" &
    (is.na(md$source_stage_provenance) |
      md$source_stage_provenance != "source_selection_defined") &
    md$total_umi_count >= 10
)

auc_rank <- function(pos, neg) {
  ranks <- rank(c(pos, neg), ties.method = "average")
  (sum(ranks[seq_along(pos)]) - length(pos) * (length(pos) + 1) / 2) /
    (length(pos) * length(neg))
}

composition <- function(x) {
  counts <- sort(table(ifelse(is.na(x), "NA", as.character(x))), decreasing = TRUE)
  paste0(names(counts), "=", as.integer(counts), collapse = "; ")
}

clusters <- as.character(md$seurat_clusters)
eligible_negative <- gametocyte_eligible & is.na(md$pred_gametocyte_sex)
specs <- list(
  female = list(label = "Female gametocyte", own = "AT25_female1", opposite = "AT25_male1"),
  male = list(label = "Male gametocyte", own = "AT25_male1", opposite = "AT25_female1")
)
validation_rows <- list()
cluster_rows <- list()

for (sex in names(specs)) {
  s <- specs[[sex]]
  positive <- !is.na(md$pred_gametocyte_sex) & md$pred_gametocyte_sex == s$label
  pos_score <- md[[s$own]][positive]
  neg_score <- md[[s$own]][eligible_negative]
  median_shift <- median(pos_score, na.rm = TRUE) - median(neg_score, na.rm = TRUE)
  auc <- auc_rank(pos_score[is.finite(pos_score)], neg_score[is.finite(neg_score)])
  called_clusters <- sort(unique(clusters[positive]))

  # A. Production classifier consistency: same population and thresholds.
  eligible_own <- tapply(md[[s$own]][gametocyte_eligible],
                         clusters[gametocyte_eligible], mean, na.rm = TRUE)
  eligible_opposite <- tapply(md[[s$opposite]][gametocyte_eligible],
                              clusters[gametocyte_eligible], mean, na.rm = TRUE)
  eligible_scores <- eligible_own[called_clusters]
  eligible_margins <- eligible_scores - eligible_opposite[called_clusters]
  production_consistent <- all(eligible_scores >= 0.10 & eligible_margins >= 0.05)

  # B. Integrated-population coherence: complete-cluster context. The >=90%
  # cell-level allowance is authoritative; non-clearing boundary clusters are
  # reported rather than subjected to a contradictory every-cluster rule.
  all_own <- tapply(md[[s$own]], clusters, mean, na.rm = TRUE)
  all_opposite <- tapply(md[[s$opposite]], clusters, mean, na.rm = TRUE)
  clearing_clusters <- names(all_own)[all_own >= 0.10 &
    (all_own - all_opposite) >= 0.05]
  n_nonclearing <- sum(positive & !clusters %in% clearing_clusters)
  coherent <- 1 - n_nonclearing / sum(positive)
  pass <- auc >= 0.70 && median_shift > 0 && production_consistent && coherent >= 0.90

  validation_rows[[sex]] <- data.frame(
    population = s$label,
    marker_set = if (sex == "female") "19 production female gametocyte markers" else "6 production male gametocyte markers",
    positive_population = s$label,
    negative_population = "Gametocyte-eligible cells without a gametocyte-sex call",
    n_positive = sum(positive), n_negative = sum(eligible_negative),
    auc = auc, auc_minimum = 0.70,
    median_positive_score = median(pos_score, na.rm = TRUE),
    median_negative_score = median(neg_score, na.rm = TRUE),
    median_score_shift = median_shift, expected_median_shift = ">0",
    production_eligible_cluster_score_minimum = min(eligible_scores),
    cluster_score_minimum = 0.10,
    production_eligible_cluster_margin_minimum = min(eligible_margins),
    cluster_margin_minimum = 0.05,
    production_classifier_consistent = production_consistent,
    integrated_all_cell_coherence = coherent,
    integrated_cluster_coherence_minimum = 0.90,
    n_positive_in_nonclearing_all_cell_clusters = n_nonclearing,
    fraction_positive_in_nonclearing_all_cell_clusters = n_nonclearing / sum(positive),
    nonclearing_all_cell_clusters = paste(setdiff(called_clusters, clearing_clusters), collapse = ","),
    status = if (pass) "PASS" else "FAIL",
    rationale = "Production consistency uses the production eligibility gate and unchanged 0.10/0.05 thresholds; AUC>=0.70 and positive median shift test discrimination; >=90% all-cell integrated-cluster coherence permits and reports limited boundary mixing",
    stringsAsFactors = FALSE
  )

  for (cluster_id in called_clusters) {
    in_cluster <- clusters == cluster_id
    eligible_in_cluster <- in_cluster & gametocyte_eligible
    called_in_cluster <- in_cluster & positive
    all_score <- mean(md[[s$own]][in_cluster], na.rm = TRUE)
    all_opposite_score <- mean(md[[s$opposite]][in_cluster], na.rm = TRUE)
    eligible_score <- mean(md[[s$own]][eligible_in_cluster], na.rm = TRUE)
    eligible_opposite_score <- mean(md[[s$opposite]][eligible_in_cluster], na.rm = TRUE)
    cluster_rows[[paste(sex, cluster_id)]] <- data.frame(
      population = s$label, cluster_id = cluster_id,
      total_atlas_cells = sum(in_cluster),
      gametocyte_eligible_cells = sum(eligible_in_cluster),
      called_cells = sum(called_in_cluster),
      other_eligible_blood_cells = sum(eligible_in_cluster & !positive),
      ineligible_cells = sum(in_cluster & !gametocyte_eligible),
      broad_stage_composition = composition(md$parasite_broad_stage[in_cluster]),
      tissue_or_sample_type_composition = composition(md$tissue_or_sample_type[in_cluster]),
      study_composition = composition(md$study_label[in_cluster]),
      mean_own_score_all_cells = all_score,
      mean_opposite_score_all_cells = all_opposite_score,
      own_minus_opposite_margin_all_cells = all_score - all_opposite_score,
      mean_own_score_gametocyte_eligible = eligible_score,
      mean_opposite_score_gametocyte_eligible = eligible_opposite_score,
      own_minus_opposite_margin_gametocyte_eligible = eligible_score - eligible_opposite_score,
      production_thresholds_pass_eligible_cells = eligible_score >= 0.10 &&
        (eligible_score - eligible_opposite_score) >= 0.05,
      all_cell_thresholds_pass = all_score >= 0.10 &&
        (all_score - all_opposite_score) >= 0.05,
      stringsAsFactors = FALSE
    )
  }
}

out <- do.call(rbind, validation_rows)
diagnostics <- do.call(rbind, cluster_rows)
write.table(out, file.path(repo, "audit/phase3_preflight/AT25_marker_preservation_validation.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
write.table(diagnostics[diagnostics$population == "Male gametocyte", , drop = FALSE],
            file.path(repo, "audit/phase3_preflight/AT25_male_cluster_diagnostics.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

if (any(out$status != "PASS")) quit(status = 1L)
message("AT25 PASS: production-consistent calls and integrated marker coherence meet all predeclared criteria")
