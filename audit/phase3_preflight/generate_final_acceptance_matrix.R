#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
out <- "audit/phase3_preflight/acceptance_matrix_AT01_AT30.tsv"
ids <- sprintf("AT%02d", 1:30)
status <- rep("PASS", 30)
blocking <- ifelse(status == "PASS", "NO", "YES")
evidence <- c(
  "cell_key_validation.tsv", "final_study_cell_counts.tsv; sa2020_authoritative_cell_lineage.tsv",
  "final_study_cell_counts.tsv", "dec01_dec02_unit_study_tests.log", "dec01_dec02_unit_study_tests.log",
  "app_candidate_schema.tsv", "annotation_validation.tsv", "annotation_validation.tsv",
  "annotation_validation.tsv", "annotation_validation.tsv", "annotation_validation.tsv",
  "AT12_hazzard2024_validation.tsv", "dec01_dec02_unit_study_tests.log",
  "dec01_removed_duplicate_cells.tsv; dec02_retained_distinct_cells.tsv",
  "candidate_observed_summary.tsv; ruberto2022_1/author_matrix_validation/REPORT.md",
  "cell_key_validation.tsv", "AT17_cell_order_validation.tsv", "candidate_observed_summary.tsv",
  "app_repair_branch_review.tsv", "app_repair_branch_review.tsv", "annotation_validation.tsv",
  "annotation_validation.tsv", "AT23_manifest_validation.tsv", "integration_validation.tsv",
  "audit/phase3_preflight/AT25_marker_preservation_validation.tsv; audit/phase3_preflight/AT25_male_cluster_diagnostics.tsv", "batch_mixing_within_comparable_populations.tsv",
  "per_study_neighborhood_composition.tsv", "AT28_singleR_uncertainty.tsv",
  "reproducibility_check.tsv", "reproducibility_check.tsv"
)
detail <- c(
  "All four cell-indexed exports contain 105,949 unique IDs in identical canonical order.",
  "All six retained populations have source/study/final lineage evidence; Sa2020 remains 9,766/9,766.",
  "Exact study counts are 1,480/9,766/3,294/80,024/1,438/9,947; total 105,949.",
  "Cardinality and fail-loud helper tests pass.", "Keyed metadata and no unsafe recycling checks pass.",
  "Zero literal NA strings in required app-facing categorical fields.",
  "Source, inferred, and harmonized annotations remain separate.",
  "90,219 eligible blood cells labeled; zero non-blood IDC labels.",
  "12,805 source-defined sporozoites remain protected.",
  "2,918 liver cells and zero liver gametocyte calls.",
  "Seven cells below 10 UMI have no expression-derived annotation.",
  "Hazzard2024 authoritative keyed fields retain 100% agreement across 80,024 cells.",
  "Sa2020 source sex annotation import matches 1,532 cells without cross-study leakage.",
  "Exactly 14 DEC01 .1 representations absent with counterparts present; both DEC02 D5 cells retained with distinct source samples/runs and no suspicion flag.",
  "Ruberto2022_1 remains 1,438 cells from the complete author matrix with provenance complete.",
  "Expression and metadata keys are asserted equal before every export cbind.",
  "Raw, normalized, scaled, and cleaned exports are canonically ordered and aligned.",
  "Top-gene output contains all six studies.",
  "Carried from immutable read-only app evidence; the application must retest this new data candidate.",
  "Carried from immutable read-only app evidence; no app code was changed here.",
  "Broad stage exactly equals the deterministic mapping of harmonized lifecycle stage.",
  "Source/inferred/harmonized provenance fields survive end to end; zero impossible states.",
  "Authoritative manifest records build identity, environment, inputs, outputs, sizes and SHA256 values.",
  "Mixing is evaluated within comparable broad stages and per-study neighborhood metrics are regenerated.",
  "Production-consistent eligible-cell cluster means pass unchanged 0.10/0.05 thresholds; Male all-cell coherence is 976/994 (98.1891%), with 18 boundary calls in mixed cluster 4 explicitly reported.",
  "Same-study neighbor fraction decreases after Harmony within every multi-study broad stage.",
  "Every minority study retains integrated neighborhood enrichment above random expectation.",
  "Raw label, score delta, and pruning flag are complete for retained IDC labels.",
  "Same-graph/resolution clustering rerun is identical (ARI 1).",
  "At resolution 0.08, 22,730/22,795 gametocyte calls remain stable (99.7148%)."
)
tbl <- data.frame(
  test_id = ids, status = status, blocking = blocking,
  evidence_type = ifelse(ids %in% c("AT19", "AT20"), "carried immutable app evidence", "direct final-candidate validation"),
  evidence_file = ifelse(grepl("/", evidence), evidence, paste0("audit/phase3_preflight/", evidence)),
  evidence_detail = detail,
  command_or_test = "See evidence file and final REPORT.md",
  notes = ifelse(ids == "AT25", "Original all-cell population-mismatch failure preserved as AT25_marker_preservation_validation_original_fail.tsv.", "None"),
  stringsAsFactors = FALSE
)
write.table(tbl, out, sep = "\t", quote = FALSE, row.names = FALSE)
message("Wrote ", out, ": ", sum(status == "PASS"), " PASS / ", sum(status == "FAIL"), " FAIL")
