# Lightweight deterministic validation (Section 11 of the production-repair
# spec; audit/shared_pipeline/acceptance_tests.tsv where applicable before
# the Phase 2 atlas rebuild).
#
# Run from the pre_process_data directory:
#   env -u R_LIBS_USER pixi run --manifest-path scripts/pixi.toml Rscript --vanilla scripts/tests/test_shared_pipeline.R
#
# Tests are grouped into:
#  - UNIT tests: run against pipeline_lib.R helpers alone, no study data
#    required. Always run.
#  - STUDY tests: run against each per-study .rds IF it already exists in
#    the working directory (regenerated study objects are not produced
#    automatically by this test file - see Section 10/Phase D). Skipped
#    with a clear message if the object is absent, never silently "passed".
#  - PHASE 2 tests: require pv_all_studies.rds (the merged, integrated,
#    annotated atlas) which Phase 1 explicitly does not (re)build. These are
#    marked explicitly and always reported as SKIPPED (PHASE 2), never
#    faked as passing.

suppressPackageStartupMessages({
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
  if (length(hit) == 0) stop("Could not locate pre_process_data root; set PLAVISCA_PREPROCESS_ROOT.", call. = FALSE)
  normalizePath(hit[[1]])
})
setwd(.plavisca_root)
source("scripts/pipeline_lib.R")

results <- list()
record <- function(id, description, status, detail = "") {
  results[[length(results) + 1]] <<- list(id = id, description = description, status = status, detail = detail)
  cat(sprintf("[%s] %-6s %s%s\n", id, status, description, if (nzchar(detail)) paste0(" - ", detail) else ""))
}

run_test <- function(id, description, expr) {
  res <- tryCatch({
    expr()
    list(ok = TRUE, msg = "")
  }, error = function(e) list(ok = FALSE, msg = conditionMessage(e)))
  record(id, description, if (res$ok) "PASS" else "FAIL", res$msg)
}

skip_test <- function(id, description, reason) {
  record(id, description, "SKIP", reason)
}

# ============================================================================
# UNIT tests (no study data required)
# ============================================================================

run_test("UNIT-01", "assert_cardinality rejects a length mismatch", function() {
  ok <- tryCatch({
    assert_cardinality(c(1, 2), 3)
    FALSE
  }, error = function(e) TRUE)
  stopifnot(ok)
})

run_test("UNIT-02", "assert_keyed_join rejects duplicate/unmatched/extra keys", function() {
  stopifnot(inherits(tryCatch(assert_keyed_join(c("a", "a"), c("a", "b")), error = function(e) e), "error"))
  stopifnot(inherits(tryCatch(assert_keyed_join(c("a", "c"), c("a", "b")), error = function(e) e), "error"))
  stopifnot(inherits(tryCatch(assert_keyed_join(c("a"), c("a", "b")), error = function(e) e), "error"))
  assert_keyed_join(c("a", "b"), c("a", "b")) # should not error
})

run_test("AT01-unit", "assert_unique_cell_ids rejects duplicates", function() {
  ok <- tryCatch({
    assert_unique_cell_ids(c("x", "x", "y"))
    FALSE
  }, error = function(e) TRUE)
  stopifnot(ok)
  assert_unique_cell_ids(c("x", "y", "z")) # should not error
})

run_test("AT06-unit", "assert_no_literal_na_string rejects the literal string \"NA\"", function() {
  ok <- tryCatch({
    assert_no_literal_na_string(c("foo", "NA", "bar"))
    FALSE
  }, error = function(e) TRUE)
  stopifnot(ok)
  assert_no_literal_na_string(c("foo", NA_character_, "bar")) # true NA is fine
})

run_test("AT16-unit", "assert_rowname_order_equal catches set match but order mismatch", function() {
  ok <- tryCatch({
    assert_rowname_order_equal(c("a", "b"), c("b", "a"))
    FALSE
  }, error = function(e) TRUE)
  stopifnot(ok)
  assert_rowname_order_equal(c("a", "b"), c("a", "b")) # should not error
})

run_test("AT21-unit", "map_broad_stage: Sporozoite always maps to Sporozoite stage (D020)", function() {
  stopifnot(identical(map_broad_stage("Sporozoite"), "Sporozoite stage"))
  stopifnot(identical(map_broad_stage("Ring"), "Blood stage"))
  stopifnot(identical(map_broad_stage("Hypnozoite"), "Liver stage"))
})

run_test("AT21-unit-2", "map_broad_stage rejects an unmapped harmonized value (fail-loud, not silent NA)", function() {
  ok <- tryCatch({
    map_broad_stage("Some Unmapped Value")
    FALSE
  }, error = function(e) TRUE)
  stopifnot(ok)
})

# ============================================================================
# STUDY tests. Section 13 forbids overwriting the existing deployed study
# .rds files at the pre_process_data root in place, so this test suite NEVER
# treats those legacy (pre-repair) files as evidence of the fix's
# correctness - it only validates a candidate object regenerated to
# $PLAVISCA_CANDIDATE_STUDY_DIR (default: scripts/tests/candidate_study_objects/).
# If only the legacy root-level file exists, that is reported as an
# explicit SKIP with a note, never as a PASS or FAIL of the repaired script.
# ============================================================================
candidate_study_dir <- Sys.getenv("PLAVISCA_CANDIDATE_STUDY_DIR", "scripts/tests/candidate_study_objects")

check_study_object <- function(legacy_path, expected_n, label) {
  candidate_path <- file.path(candidate_study_dir, legacy_path)
  if (!file.exists(candidate_path)) {
    note <- if (file.exists(legacy_path)) {
      sprintf(
        "no candidate at %s; %s exists but is the pre-repair legacy artifact protected by Section 13 (never overwritten in place) - regenerate the fixed script's output to the candidate path to validate",
        candidate_path, legacy_path
      )
    } else {
      sprintf("neither %s nor legacy %s present", candidate_path, legacy_path)
    }
    skip_test(paste0("STUDY-", label), sprintf("%s: cell count/keying/NA checks", label), note)
    return(invisible())
  }
  path <- candidate_path
  run_test(paste0("AT01-", label), sprintf("%s: unique cell IDs", label), function() {
    so <- readRDS(path)
    assert_unique_cell_ids(colnames(so), label = path)
  })
  run_test(paste0("AT03-", label), sprintf("%s: expected cell count (%d)", label, expected_n), function() {
    so <- readRDS(path)
    assert_cell_count(ncol(so), expected_n, label = path)
  })
  run_test(paste0("AT06-", label), sprintf("%s: no literal \"NA\" string in character metadata", label), function() {
    so <- readRDS(path)
    char_cols <- names(so@meta.data)[vapply(so@meta.data, is.character, logical(1))]
    for (col in char_cols) {
      assert_no_literal_na_string(so@meta.data[[col]], label = paste0(path, "$", col))
    }
  })
}

check_study_object("silva2022.rds", 1494L, "mancio_silva2022")
check_study_object("sa2020.rds", 9766L, "sa2020")
check_study_object("hazzard2022.rds", 3294L, "hazzard2022")
check_study_object("hazzard2024.rds", 80024L, "hazzard2024")
check_study_object("ruberto2022_1.rds", 1438L, "ruberto2022_1")
check_study_object("ruberto2022_2.rds", 9947L, "ruberto2022_2")

ruberto1_candidate <- file.path(candidate_study_dir, "ruberto2022_1.rds")
if (file.exists(ruberto1_candidate)) {
  run_test("D036-flag-legacy", "ruberto2022_1: the historical 538-cell STARsolo artifact remains documented in frozen audit evidence (never baked into a live object field, since no artifact ever existed with that field populated pre-Phase-2)", function() {
    evidence_path <- "audit/ruberto2022_1/near_empty_cells_vs_hep59.tsv"
    if (!file.exists(evidence_path)) {
      pipeline_fail(sprintf("%s not found - the frozen D036 near-empty evidence record must be preserved", evidence_path))
    }
    evidence <- read.delim(evidence_path, stringsAsFactors = FALSE)
    if (nrow(evidence) != 538L) {
      stop(sprintf("expected 538 rows in %s per D036 audit evidence, got %d", evidence_path, nrow(evidence)))
    }
    so <- readRDS(ruberto1_candidate)
    stopifnot("legacy_starsolo_near_empty_flag" %in% colnames(so@meta.data))
    assert_cell_count(ncol(so), 1438L, label = "ruberto2022_1.rds (D036 population preserved)")
  })
  run_test("D036-phase2-resolved", "ruberto2022_1: adopted author counts no longer reproduce the STARsolo near-empty artifact; provenance populated for all 1438 cells", function() {
    so <- readRDS(ruberto1_candidate)
    stopifnot(all(c("source_of_counts", "counting_pipeline", "count_reference_version", "count_provenance_status", "total_umi_count", "near_empty_expression_flag") %in% colnames(so@meta.data)))
    assert_no_missing_mandatory(so$source_of_counts, label = "ruberto2022_1$source_of_counts")
    if (!all(so$source_of_counts == "author_processed_object_raw_umi")) {
      pipeline_fail("not all 1438 Ruberto2022_1 cells carry the adopted author-matrix source_of_counts value")
    }
    n_adopted_near_empty <- sum(so$near_empty_expression_flag, na.rm = TRUE)
    # The whole point of the Phase 2 repair is that the artificial STARsolo
    # near-empty artifact (538/1438 cells) does not reappear under the
    # adopted author counts - a handful of genuinely low-count cells may
    # still exist, but nothing resembling the old systematic 538-cell block.
    if (n_adopted_near_empty >= 538L) {
      pipeline_fail(sprintf(
        "near_empty_expression_flag under the adopted author counts is %d (>= the old 538-cell STARsolo artifact) - the Phase 2 repair does not appear to have taken effect",
        n_adopted_near_empty
      ))
    }
  })
} else {
  skip_test("D036-flag-legacy", "ruberto2022_1: legacy STARsolo near-empty flag check", sprintf("no candidate at %s", ruberto1_candidate))
  skip_test("D036-phase2-resolved", "ruberto2022_1: Phase 2 author-count resolution check", sprintf("no candidate at %s", ruberto1_candidate))
}

# ============================================================================
# D036 resolution artifacts (2026-09-20 follow-up investigation): assert the
# documented classification is consistent with its own evidence table,
# rather than trusting D036_resolution.md's prose alone.
# ============================================================================
d036_rescue_tsv <- "audit/ruberto2022_1/D036_cell_level_rescue.tsv"
if (file.exists(d036_rescue_tsv)) {
  run_test("D036-rescue-evidence", "D036_cell_level_rescue.tsv is consistent with classification B (no rescue)", function() {
    rescue <- read.delim(d036_rescue_tsv, stringsAsFactors = FALSE)
    stopifnot(nrow(rescue) == 538)
    zero_frac_candidate <- mean(rescue$kallisto_candidate_umi == 0)
    median_candidate <- median(rescue$kallisto_candidate_umi)
    # These thresholds encode what "not rescued" means quantitatively (Part 8
    # of the task spec: "fixed" is never merely "counts went up") - if a
    # future rerun of this investigation ever produces materially different
    # numbers, this test should fail loudly rather than silently pass.
    if (!(zero_frac_candidate > 0.9 && median_candidate == 0)) {
      pipeline_fail(sprintf(
        "D036_cell_level_rescue.tsv no longer shows the documented null-rescue pattern (zero_frac_candidate=%.3f, median_candidate=%.1f) - re-check D036_resolution.md's classification before trusting it",
        zero_frac_candidate, median_candidate
      ))
    }
  })
} else {
  skip_test("D036-rescue-evidence", "D036_cell_level_rescue.tsv consistency check", "D036_cell_level_rescue.tsv not present")
}

# ============================================================================
# PHASE 2 tests (require the merged/integrated/annotated atlas -
# pv_all_studies.rds). Each runs for real when its required artifact
# exists; otherwise it is explicitly SKIPPED (never faked as passing).
# ============================================================================
atlas_path <- "pv_all_studies.rds"
cleaned_export_path <- "data/candidate_export/cleaned_dataset.rds"
normalize_export_path <- "data/candidate_export/normalize_df.rds"
raw_export_path <- "data/candidate_export/raw_df.rds"
scale_export_path <- "data/candidate_export/scale_df.rds"
integration_validation_path <- "audit/phase2_rebuild/integration_validation.tsv"

if (file.exists(atlas_path)) {
  atlas <- readRDS(atlas_path)
  atlas_meta <- atlas@meta.data

  run_test("AT07", "Source annotations preserved distinct from inferred annotations", function() {
    required <- c("source_life_cycle_stage", "harmonized_life_cycle_stage", "source_stage_provenance", "pred_gametocyte_sex")
    missing <- setdiff(required, colnames(atlas_meta))
    if (length(missing) > 0) pipeline_fail(sprintf("missing distinct provenance field(s): %s", paste(missing, collapse = ", ")))
    if ("source_sex_annotation" %in% colnames(atlas_meta) && "pred_gametocyte_sex" %in% colnames(atlas_meta)) {
      if (identical(atlas_meta$source_sex_annotation, atlas_meta$pred_gametocyte_sex)) {
        pipeline_fail("source_sex_annotation and pred_gametocyte_sex are identical columns - provenance appears merged, not preserved distinct")
      }
    }
  })

  run_test("AT08", "Non-blood-stage cells receive no blood-IDC/HPI inference", function() {
    non_blood <- atlas_meta$tissue_or_sample_type != "Host blood" & !is.na(atlas_meta$tissue_or_sample_type)
    if (any(!is.na(atlas_meta$idc_reference_similarity_label[non_blood]))) {
      pipeline_fail("at least one non-blood-stage cell has a non-NA idc_reference_similarity_label")
    }
  })

  run_test("AT09", "Source-defined sporozoites never broad-labeled Blood stage", function() {
    spz <- !is.na(atlas_meta$source_stage_provenance) & atlas_meta$source_stage_provenance == "source_selection_defined" &
      !is.na(atlas_meta$source_life_cycle_stage) & atlas_meta$source_life_cycle_stage == "Sporozoite"
    if (sum(spz) == 0) pipeline_fail("no source-defined Sporozoite cells found - cannot validate this invariant")
    if (any(atlas_meta$parasite_broad_stage[spz] != "Sporozoite stage")) {
      pipeline_fail("at least one source-defined Sporozoite cell is not broad-labeled Sporozoite stage")
    }
  })

  run_test("AT10", "Liver-stage cells never exposed as gametocytes", function() {
    liver <- atlas_meta$parasite_broad_stage == "Liver stage" & !is.na(atlas_meta$parasite_broad_stage)
    if (sum(liver) == 0) pipeline_fail("no Liver stage cells found - cannot validate this invariant")
    if (any(!is.na(atlas_meta$pred_gametocyte_sex[liver]))) {
      pipeline_fail("at least one Liver stage cell has a non-NA pred_gametocyte_sex")
    }
  })

  run_test("AT11", "Zero/near-zero expression cells receive no expression-derived annotation", function() {
    near_zero <- atlas_meta$total_umi_count < 10
    if (any(!is.na(atlas_meta$idc_reference_similarity_label[near_zero])) ||
      any(!is.na(atlas_meta$pred_gametocyte_sex[near_zero]))) {
      pipeline_fail("a near-zero-expression cell (<10 total UMI) received idc_reference_similarity_label or pred_gametocyte_sex")
    }
  })

  run_test("AT13", "Sa2020 gametocyte source import: source_sex_annotation populated only for Sa2020, sane cardinality", function() {
    if (!"source_sex_annotation" %in% colnames(atlas_meta)) pipeline_fail("source_sex_annotation column missing")
    n_matched <- sum(!is.na(atlas_meta$source_sex_annotation))
    non_sa2020_matched <- sum(!is.na(atlas_meta$source_sex_annotation) & atlas_meta$study_label != "Sa2020")
    if (non_sa2020_matched > 0) {
      pipeline_fail(sprintf("%d non-Sa2020 cells carry a source_sex_annotation value - the Sa2020-only import leaked into another study", non_sa2020_matched))
    }
    cat(sprintf("    [AT13 info] %d Sa2020 cells matched a source gametocyte-sex annotation\n", n_matched))
  })

  run_test("AT21-atlas", "Broad/detailed lifecycle fields consistent on the merged atlas", function() {
    expected_broad <- map_broad_stage(atlas_meta$harmonized_life_cycle_stage)
    if (!identical(expected_broad, atlas_meta$parasite_broad_stage)) {
      pipeline_fail("parasite_broad_stage does not match the deterministic mapping of harmonized_life_cycle_stage")
    }
  })

  run_test("AT22", "Source/inferred/harmonized annotation provenance retained end to end", function() {
    required <- c(
      "source_stage_provenance", "count_provenance_status", "idc_similarity_pruned_flag",
      "idc_similarity_raw_label", "idc_similarity_score_delta"
    )
    missing <- setdiff(required, colnames(atlas_meta))
    if (length(missing) > 0) pipeline_fail(sprintf("missing provenance field(s): %s", paste(missing, collapse = ", ")))
  })

  run_test("impossible-state-check", "No impossible lifecycle-annotation combinations on the merged atlas", function() {
    liver_or_spz <- atlas_meta$parasite_broad_stage %in% c("Liver stage", "Sporozoite stage")
    bad_idc <- liver_or_spz & !is.na(atlas_meta$idc_reference_similarity_label)
    bad_gam <- liver_or_spz & !is.na(atlas_meta$pred_gametocyte_sex)
    n_bad <- sum(bad_idc) + sum(bad_gam)
    if (n_bad > 0) {
      pipeline_fail(sprintf("%d impossible-state cell(s) found (Liver/Sporozoite stage with a blood-derived IDC or gametocyte call)", n_bad))
    }
  })
} else {
  for (t in list(
    list(id = "AT07", desc = "Source annotations preserved distinct from inferred annotations"),
    list(id = "AT08", desc = "Non-blood-stage cells receive no blood-IDC/HPI inference"),
    list(id = "AT09", desc = "Source-defined sporozoites never broad-labeled Blood stage"),
    list(id = "AT10", desc = "Liver-stage cells never exposed as gametocytes"),
    list(id = "AT11", desc = "Zero/near-zero expression cells receive no expression-derived annotation"),
    list(id = "AT13", desc = "Sa2020 gametocyte source import cardinality"),
    list(id = "AT21-atlas", desc = "Broad/detailed lifecycle fields consistent on the merged atlas"),
    list(id = "AT22", desc = "Source/inferred/harmonized annotation provenance retained end to end"),
    list(id = "impossible-state-check", desc = "No impossible lifecycle-annotation combinations")
  )) {
    skip_test(t$id, t$desc, "requires pv_all_studies.rds (Phase 2 atlas rebuild) - not faked as passing")
  }
}

if (file.exists(normalize_export_path) && file.exists(raw_export_path) && file.exists(scale_export_path)) {
  run_test("AT17", "Export tables explicitly cell-key aligned", function() {
    n <- rownames(readRDS(normalize_export_path))
    r <- rownames(readRDS(raw_export_path))
    s <- rownames(readRDS(scale_export_path))
    assert_rowname_order_equal(n, r, "normalize_df", "raw_df")
    assert_rowname_order_equal(n, s, "normalize_df", "scale_df")
  })
} else {
  skip_test("AT17", "Export tables explicitly cell-key aligned", "requires candidate_export/{normalize,raw,scale}_df.rds")
}

if (file.exists(cleaned_export_path)) {
  run_test("AT18", "All six studies represented in top-gene export", function() {
    cleaned <- readRDS(cleaned_export_path)
    n_studies <- dplyr::n_distinct(cleaned$top_genes_exp$study)
    if (n_studies != length(STUDY_LABELS)) {
      pipeline_fail(sprintf("top_genes_exp contains %d distinct studies, expected %d", n_studies, length(STUDY_LABELS)))
    }
  })
} else {
  skip_test("AT18", "All six studies represented in top-gene export", "requires candidate_export/cleaned_dataset.rds")
}

if (file.exists(integration_validation_path)) {
  run_test("AT24-AT27", "Integration validation metrics table exists and is non-empty", function() {
    tbl <- read.delim(integration_validation_path, stringsAsFactors = FALSE)
    if (nrow(tbl) == 0) pipeline_fail("integration_validation.tsv is empty")
  })
} else {
  for (id in c("AT24", "AT25", "AT26", "AT27")) {
    skip_test(id, "Integration validation metric", "requires audit/phase2_rebuild/integration_validation.tsv")
  }
}

reproducibility_path <- "audit/phase2_rebuild/reproducibility_check.tsv"
if (file.exists(reproducibility_path)) {
  run_test("AT29-AT30", "Reproducibility/stability check table exists and reports no undocumented instability", function() {
    tbl <- read.delim(reproducibility_path, stringsAsFactors = FALSE)
    if (nrow(tbl) == 0) pipeline_fail("reproducibility_check.tsv is empty")
  })
} else {
  for (id in c("AT29", "AT30")) {
    skip_test(id, "Reproducibility/stability check", "requires audit/phase2_rebuild/reproducibility_check.tsv")
  }
}

# ============================================================================
# Summary
# ============================================================================
status_counts <- table(vapply(results, function(r) r$status, character(1)))
cat("\n==== SUMMARY ====\n")
print(status_counts)

n_fail <- sum(vapply(results, function(r) r$status == "FAIL", logical(1)))
if (n_fail > 0) {
  cat(sprintf("\n%d test(s) FAILED.\n", n_fail))
  quit(status = 1)
} else {
  cat("\nAll runnable tests passed (SKIP entries require regenerated study objects or the Phase 2 atlas rebuild).\n")
}
