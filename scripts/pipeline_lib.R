# pipeline_lib.R
#
# Shared production-repair helpers for the PlaViSca pre-processing pipeline.
# Implements the fail-loud assertion discipline, controlled vocabularies, and
# deterministic mapping tables specified in
# audit/shared_pipeline/{metadata_schema,annotation_precedence,acceptance_tests}.tsv.
#
# Sourced by every per-study preprocessing script and by singleR.R/
# integration.R/flatten_data.R. Never itself reads or writes study data.

suppressPackageStartupMessages(library(dplyr))

# ---------------------------------------------------------------------------
# Controlled vocabularies (fixes D007/D056: one canonical string per study,
# never re-typed per script)
# ---------------------------------------------------------------------------

STUDY_LABELS <- c(
  mancio_silva2022 = "Mancio-Silva2022",
  sa2020 = "Sa2020",
  ruberto2022_1 = "Ruberto2022_1",
  ruberto2022_2 = "Ruberto2022_2",
  hazzard2022 = "Hazzard2022",
  hazzard2024 = "Hazzard2024"
)

# parasite_broad_stage is a deterministic, version-controlled image of
# harmonized_life_cycle_stage (metadata_schema.tsv row 41). It must never be
# independently overwritten - always derive it with map_broad_stage().
HARMONIZED_TO_BROAD_STAGE <- c(
  "EEF" = "Liver stage",
  "Hypnozoite" = "Liver stage",
  "Schizont (Liver stage)" = "Liver stage",
  "Ring" = "Blood stage",
  "Trophozoite" = "Blood stage",
  "Schizont (Blood stage)" = "Blood stage",
  "Merozoite" = "Blood stage",
  "Male gametocyte" = "Blood stage",
  "Female gametocyte" = "Blood stage",
  "Asexual (unresolved stage)" = "Blood stage",
  "Sporozoite" = "Sporozoite stage"
)

# tissue_or_sample_type canonical vocabulary (fixes the D058 vocabulary-drift
# root cause: "Host blood" vs "Mammalian host: blood" etc. must never diverge
# again). Every per-study script must emit exactly one of these values.
TISSUE_OR_SAMPLE_TYPE <- c(
  "Host blood",
  "Host liver",
  "Vector salivary gland",
  "Vector midgut"
)

# Populations that are never eligible for blood-stage IDC reference
# similarity classification (D058), regardless of expression magnitude.
NON_BLOOD_TISSUE_TYPES <- c("Host liver", "Vector salivary gland", "Vector midgut")

# ---------------------------------------------------------------------------
# Fail-loud assertion helpers (AT04, AT05: no vector-length mismatch, no
# positional recycling ever silently accepted)
# ---------------------------------------------------------------------------

#' Stop with a clearly labeled pipeline error (never a silent warning).
pipeline_fail <- function(msg) {
  stop(paste0("[PIPELINE ASSERTION FAILED] ", msg), call. = FALSE)
}

#' Assert that `x` has exactly `n` elements before it is used as a per-cell
#' metadata vector. Fixes D002,D027,D028,D033 (AT04).
assert_cardinality <- function(x, n, label = deparse(substitute(x))) {
  if (length(x) != n) {
    pipeline_fail(sprintf(
      "cardinality mismatch for '%s': length %d, expected %d (per-cell assignment must never rely on recycling)",
      label, length(x), n
    ))
  }
  invisible(TRUE)
}

#' Assert that every id in `keys` is present exactly once in `table_keys`
#' (no duplicate run keys, no unmatched runs, no unexpected extra runs).
#' Fixes AT05 and the "fail loudly" requirements in section 1 of the repair
#' spec.
assert_keyed_join <- function(keys, table_keys, label = "keys") {
  dup_keys <- table_keys[duplicated(table_keys)]
  if (length(dup_keys) > 0) {
    pipeline_fail(sprintf(
      "duplicate run keys in reference table for '%s': %s",
      label, paste(unique(dup_keys), collapse = ", ")
    ))
  }
  unmatched <- setdiff(keys, table_keys)
  if (length(unmatched) > 0) {
    pipeline_fail(sprintf(
      "unmatched runs for '%s': %s (present in data, absent from keyed reference table)",
      label, paste(unmatched, collapse = ", ")
    ))
  }
  extra <- setdiff(table_keys, keys)
  if (length(extra) > 0) {
    pipeline_fail(sprintf(
      "unexpected extra runs in keyed reference table for '%s': %s (present in table, absent from data)",
      label, paste(extra, collapse = ", ")
    ))
  }
  invisible(TRUE)
}

#' Assert no duplicated final cell IDs (AT01).
assert_unique_cell_ids <- function(ids, label = "cell IDs") {
  dups <- ids[duplicated(ids)]
  if (length(dups) > 0) {
    pipeline_fail(sprintf(
      "duplicated final cell IDs for '%s': %d duplicates (e.g. %s)",
      label, length(dups), paste(utils::head(unique(dups), 5), collapse = ", ")
    ))
  }
  invisible(TRUE)
}

#' Assert an object has exactly `expected_n` cells (unexpected cell count).
assert_cell_count <- function(object_ncol, expected_n, label = "object") {
  if (object_ncol != expected_n) {
    pipeline_fail(sprintf(
      "unexpected cell count for '%s': got %d, expected %d",
      label, object_ncol, expected_n
    ))
  }
  invisible(TRUE)
}

#' Assert no mandatory source ID is missing (true NA disallowed for
#' mandatory identity fields).
assert_no_missing_mandatory <- function(x, label = deparse(substitute(x))) {
  if (anyNA(x) || any(x == "" & !is.na(x))) {
    pipeline_fail(sprintf(
      "missing mandatory source ID(s) in '%s': %d NA/empty of %d",
      label, sum(is.na(x) | x == ""), length(x)
    ))
  }
  invisible(TRUE)
}

#' Assert a character vector never contains the literal string "NA" (true NA
#' only). Fixes D024,D030,D062 (AT06).
assert_no_literal_na_string <- function(x, label = deparse(substitute(x))) {
  bad <- !is.na(x) & x == "NA"
  if (any(bad)) {
    pipeline_fail(sprintf(
      "literal string \"NA\" found in '%s' (%d cells) - use true NA (NA_character_), never the string \"NA\"",
      label, sum(bad)
    ))
  }
  invisible(TRUE)
}

#' Assert exact rowname-set equality between a metadata table and an
#' expression matrix before any join/cbind (AT16,AT17; fixes D053,D055).
assert_rowname_set_equal <- function(a_names, b_names, label_a = "a", label_b = "b") {
  if (!identical(sort(a_names), sort(b_names))) {
    only_a <- setdiff(a_names, b_names)
    only_b <- setdiff(b_names, a_names)
    pipeline_fail(sprintf(
      "rowname set mismatch between '%s' and '%s': %d only in %s, %d only in %s",
      label_a, label_b, length(only_a), label_a, length(only_b), label_b
    ))
  }
  invisible(TRUE)
}

#' Assert exact rowname order equality (stronger than set equality; used
#' immediately before a positional cbind, AT16).
assert_rowname_order_equal <- function(a_names, b_names, label_a = "a", label_b = "b") {
  assert_rowname_set_equal(a_names, b_names, label_a, label_b)
  if (!identical(a_names, b_names)) {
    pipeline_fail(sprintf(
      "rowname ORDER mismatch between '%s' and '%s' (sets match but order differs) - reorder explicitly before cbind",
      label_a, label_b
    ))
  }
  invisible(TRUE)
}

# ---------------------------------------------------------------------------
# Deterministic broad-stage derivation (fixes D020: makes "Sporozoite
# detailed stage + Blood broad stage" impossible by construction; AT09,AT21)
# ---------------------------------------------------------------------------

#' Deterministically derive parasite_broad_stage from
#' harmonized_life_cycle_stage. Never call this ad hoc with a hand-written
#' case_when - always route through this single mapping table so broad/
#' detailed can never silently diverge again.
map_broad_stage <- function(harmonized_life_cycle_stage) {
  unmapped <- setdiff(
    unique(harmonized_life_cycle_stage[!is.na(harmonized_life_cycle_stage)]),
    names(HARMONIZED_TO_BROAD_STAGE)
  )
  if (length(unmapped) > 0) {
    pipeline_fail(sprintf(
      "harmonized_life_cycle_stage value(s) not present in HARMONIZED_TO_BROAD_STAGE mapping table: %s",
      paste(unmapped, collapse = ", ")
    ))
  }
  unname(HARMONIZED_TO_BROAD_STAGE[harmonized_life_cycle_stage])
}

# ---------------------------------------------------------------------------
# Keyed run-table metadata application (fixes D002,D026-D028,D033: replaces
# `rep()`/positional-vector per-cell metadata assignment with an explicit
# keyed join against a committed run table, with cardinality/key assertions)
# ---------------------------------------------------------------------------

#' Split a list of per-run Seurat objects into those that ARE keyed in
#' run_table (survive to the final retained population) and those that are
#' NOT (e.g. Hazzard2024's imported-but-not-QC-retained runs, which a later
#' explicit subset() step in the merge script removes before saveRDS). This
#' makes explicit, at the earliest responsible step, exactly which runs are
#' authoritative vs transient/pre-filter - never a silent positional
#' assumption.
split_run_membership <- function(seurat_list, run_table) {
  run_ids <- names(seurat_list)
  list(
    retained = seurat_list[run_ids %in% run_table$run_id],
    pending_filter = seurat_list[!run_ids %in% run_table$run_id]
  )
}

#' Apply per-run metadata from a keyed run table to a list of per-run Seurat
#' objects, matching strictly on run_id. Fails loudly on any duplicate,
#' unmatched, or unexpected-extra run key (AT05).
#'
#' @param seurat_list named list of per-run Seurat objects, names = run_id
#' @param run_table data.frame with one row per run_id and the metadata
#'   columns to assign (must include column `run_id`)
#' @param fields character vector of column names in run_table to copy onto
#'   each Seurat object's metadata
apply_keyed_run_metadata <- function(seurat_list, run_table, fields) {
  run_ids <- names(seurat_list)
  assert_keyed_join(run_ids, run_table$run_id, label = "run_id")
  missing_fields <- setdiff(fields, colnames(run_table))
  if (length(missing_fields) > 0) {
    pipeline_fail(sprintf(
      "run_table is missing required field(s): %s",
      paste(missing_fields, collapse = ", ")
    ))
  }
  for (rid in run_ids) {
    row <- run_table[run_table$run_id == rid, , drop = FALSE]
    if (nrow(row) != 1) {
      pipeline_fail(sprintf("run_table has %d rows for run_id '%s', expected exactly 1", nrow(row), rid))
    }
    n_cells <- ncol(seurat_list[[rid]])
    for (f in fields) {
      seurat_list[[rid]][[f]] <- rep(row[[f]], n_cells)
    }
  }
  seurat_list
}

# ---------------------------------------------------------------------------
# Build manifest (Part 9 / AT23)
# ---------------------------------------------------------------------------

#' Record a build-manifest row for a produced artifact. Appends to
#' build_manifest.tsv at the repository root. Never overwrites prior rows.
record_build_manifest <- function(
  artifact_path,
  script_path,
  cell_count = NA_integer_,
  seed = NA_integer_,
  notes = "",
  manifest_path = "build_manifest.tsv"
) {
  sha256 <- function(path) {
    if (!file.exists(path)) {
      return(NA_character_)
    }
    tryCatch(
      as.character(openssl::sha256(file(path, raw = TRUE))),
      error = function(e) {
        tryCatch(
          digest::digest(file = path, algo = "sha256"),
          error = function(e2) NA_character_
        )
      }
    )
  }

  git_commit <- tryCatch(
    system2("git", c("rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE),
    error = function(e) NA_character_
  )
  if (length(git_commit) == 0) git_commit <- NA_character_

  pkgs <- c("Seurat", "SeuratObject", "harmony", "SingleR")
  pkg_versions <- vapply(pkgs, function(p) {
    tryCatch(as.character(utils::packageVersion(p)), error = function(e) NA_character_)
  }, character(1))

  row <- data.frame(
    date = as.character(Sys.time()),
    artifact_path = artifact_path,
    output_sha256 = sha256(artifact_path),
    script_path = script_path,
    git_commit = if (length(git_commit) > 0) git_commit[1] else NA_character_,
    r_version = paste(R.version$major, R.version$minor, sep = "."),
    Seurat_version = pkg_versions[["Seurat"]],
    SeuratObject_version = pkg_versions[["SeuratObject"]],
    harmony_version = pkg_versions[["harmony"]],
    SingleR_version = pkg_versions[["SingleR"]],
    cell_count = cell_count,
    seed = seed,
    notes = notes,
    stringsAsFactors = FALSE
  )

  if (file.exists(manifest_path)) {
    utils::write.table(
      row,
      manifest_path,
      append = TRUE,
      sep = "\t",
      row.names = FALSE,
      col.names = FALSE,
      quote = FALSE
    )
  } else {
    utils::write.table(
      row,
      manifest_path,
      append = FALSE,
      sep = "\t",
      row.names = FALSE,
      col.names = TRUE,
      quote = FALSE
    )
  }
  invisible(row)
}
