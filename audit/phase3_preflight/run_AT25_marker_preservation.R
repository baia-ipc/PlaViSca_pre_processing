#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(Seurat))

repo <- normalizePath(getwd())
candidate_root <- Sys.getenv("PLAVISCA_CANDIDATE_ROOT", "/home/baia/prj/plavisca/pre_process_data")
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

# Recreate the production module scores on an in-memory copy only. This reads
# the saved normalized RNA layer and does not alter or serialize the atlas.
so <- AddModuleScore(so, features = list(female_markers), name = "AT25_female", seed = 1)
so <- AddModuleScore(so, features = list(male_markers), name = "AT25_male", seed = 1)
md <- so@meta.data

# Predeclared criterion, derived from the production classifier rather than
# tuned to these results. A marker-defined sex population passes when:
# (1) its own module score separates it from eligible uncalled blood cells
#     with AUC >= 0.70 (useful discrimination) and a positive median shift,
# (2) every integrated cluster carrying that call retains the production
#     mean-score threshold >= 0.10 and own-minus-opposite margin >= 0.05, and
# (3) >= 90% of called cells belong to those threshold-clearing integrated
#     clusters (coherent representation; 90% allows limited boundary mixing).
auc_rank <- function(pos, neg) {
  values <- c(pos, neg)
  ranks <- rank(values, ties.method = "average")
  (sum(ranks[seq_along(pos)]) - length(pos) * (length(pos) + 1) / 2) /
    (length(pos) * length(neg))
}

eligible_negative <- md$parasite_broad_stage %in% "Blood stage" & is.na(md$pred_gametocyte_sex)
specs <- list(
  female = list(label = "Female gametocyte", own = "AT25_female1", opposite = "AT25_male1"),
  male = list(label = "Male gametocyte", own = "AT25_male1", opposite = "AT25_female1")
)

rows <- lapply(names(specs), function(sex) {
  s <- specs[[sex]]
  positive <- md$pred_gametocyte_sex == s$label & !is.na(md$pred_gametocyte_sex)
  pos_score <- md[[s$own]][positive]
  neg_score <- md[[s$own]][eligible_negative]
  median_shift <- median(pos_score, na.rm = TRUE) - median(neg_score, na.rm = TRUE)
  auc <- auc_rank(pos_score[is.finite(pos_score)], neg_score[is.finite(neg_score)])

  clusters <- as.character(md$seurat_clusters)
  called_clusters <- unique(clusters[positive])
  cluster_own <- tapply(md[[s$own]], clusters, mean, na.rm = TRUE)
  cluster_opposite <- tapply(md[[s$opposite]], clusters, mean, na.rm = TRUE)
  minimum_called_cluster_score <- min(cluster_own[called_clusters])
  minimum_called_cluster_margin <- min(cluster_own[called_clusters] - cluster_opposite[called_clusters])
  clearing_clusters <- names(cluster_own)[cluster_own >= 0.10 & (cluster_own - cluster_opposite) >= 0.05]
  coherent <- mean(clusters[positive] %in% clearing_clusters)
  pass <- auc >= 0.70 && median_shift > 0 && minimum_called_cluster_score >= 0.10 &&
    minimum_called_cluster_margin >= 0.05 && coherent >= 0.90
  data.frame(
    population = s$label,
    marker_set = if (sex == "female") "19 production female gametocyte markers" else "6 production male gametocyte markers",
    positive_population = s$label,
    negative_population = "Blood-stage eligible cells without a gametocyte-sex call",
    n_positive = sum(positive), n_negative = sum(eligible_negative),
    auc = auc, auc_minimum = 0.70,
    median_positive_score = median(pos_score, na.rm = TRUE),
    median_negative_score = median(neg_score, na.rm = TRUE),
    median_score_shift = median_shift, expected_median_shift = ">0",
    minimum_called_cluster_score = minimum_called_cluster_score, cluster_score_minimum = 0.10,
    minimum_called_cluster_margin = minimum_called_cluster_margin, cluster_margin_minimum = 0.05,
    integrated_cluster_coherence = coherent, integrated_cluster_coherence_minimum = 0.90,
    status = if (pass) "PASS" else "FAIL",
    rationale = "AUC>=0.70 is useful discrimination; positive median shift supplies direction; 0.10 activation and 0.05 sex margin are production classifier thresholds; >=90% integrated-cluster coherence permits limited boundary mixing",
    stringsAsFactors = FALSE
  )
})

out <- do.call(rbind, rows)
write.table(out, file.path(repo, "audit/phase3_preflight/AT25_marker_preservation_validation.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
if (any(out$status != "PASS")) quit(status = 1L)
message("AT25 PASS: both marker-defined sex populations meet all predeclared criteria")
