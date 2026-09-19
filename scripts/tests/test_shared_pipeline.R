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
  run_test("D036-flag", "ruberto2022_1: 538 cells explicitly flagged near_empty_expression_flag, none excluded", function() {
    so <- readRDS(ruberto1_candidate)
    stopifnot("near_empty_expression_flag" %in% colnames(so@meta.data))
    n_flagged <- sum(so$near_empty_expression_flag, na.rm = TRUE)
    if (n_flagged != 538L) {
      stop(sprintf("expected 538 flagged cells per D036 audit evidence, got %d", n_flagged))
    }
    assert_cell_count(ncol(so), 1438L, label = "ruberto2022_1.rds (D036 population preserved)")
  })
} else {
  skip_test("D036-flag", "ruberto2022_1: D036 near-empty flag check", sprintf("no candidate at %s", ruberto1_candidate))
}

# ============================================================================
# PHASE 2 tests (require the merged/integrated/annotated atlas -
# pv_all_studies.rds - which Phase 1 does not regenerate)
# ============================================================================

phase2_tests <- list(
  list(id = "AT07", desc = "Source annotations preserved distinct from inferred annotations"),
  list(id = "AT08", desc = "Non-blood-stage cells receive no blood-IDC/HPI inference"),
  list(id = "AT09", desc = "Source-defined sporozoites never broad-labeled Blood stage"),
  list(id = "AT10", desc = "Liver-stage cells never exposed as gametocytes"),
  list(id = "AT11", desc = "Zero/near-zero expression cells receive no expression-derived annotation"),
  list(id = "AT13", desc = "Sa2020 gametocyte source import reproduces exactly 1532 matched rows"),
  list(id = "AT17", desc = "Export tables explicitly cell-key aligned"),
  list(id = "AT18", desc = "All six studies represented in top-gene export"),
  list(id = "AT22", desc = "Source/inferred/harmonized annotation provenance retained end to end"),
  list(id = "AT24", desc = "Integration validation: batch mixing metrics"),
  list(id = "AT25", desc = "Integration validation: marker gene preservation"),
  list(id = "AT26", desc = "Integration validation: batch mixing only within comparable populations"),
  list(id = "AT27", desc = "Integration validation: study dominance does not erase minority signal"),
  list(id = "AT29", desc = "Same-seed clustering determinism"),
  list(id = "AT30", desc = "Cluster identity corresponds to marker-defined populations under perturbation")
)
for (t in phase2_tests) {
  skip_test(t$id, t$desc, "requires pv_all_studies.rds (Phase 2 atlas rebuild) - not faked as passing in Phase 1")
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
