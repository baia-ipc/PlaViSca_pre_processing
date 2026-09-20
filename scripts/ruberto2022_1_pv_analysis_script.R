##### --------------------------------------------------------------------------
##### -Ruberto, CTEGD, Institute of Bioinformatics, UGA
##### --------------------------------------------------------------------------
#####
##### PHASE 2 (2026-09-20): this study's expression source is now the
##### author-derived raw UMI matrix (data/Hep59.1.2.seu_20aug2025.rds), not
##### PlaViSca's own STARsolo reconstruction - see
##### audit/ruberto2022_1/D036_resolution.md "PHASE 2 UPDATE" and
##### audit/shared_pipeline/decision_register.tsv DEC09 2026-09-20 UPDATE.
#####
##### This script therefore no longer reads counts/36093191/*_solo_out/ at
##### all: cell membership (the exact 1,438 cells), per-run metadata, and
##### expression counts are all derivable without touching the STARsolo raw
##### barcode space (millions of candidate barcodes) that the pre-Phase-2
##### version of this script loaded solely to build an object it then threw
##### 99.9% of away. The only STARsolo-derived information retained is a
##### historical/diagnostic comparison (legacy_starsolo_total_umi_count/
##### legacy_starsolo_near_empty_flag), reused directly from this study's
##### existing Phase-1 output object rather than recomputed from raw counts.
##### --------------------------------------------------------------------------

##### --------------------------------------------------------------------------
##### 1. Load libraries
##### --------------------------------------------------------------------------
library(tidyverse)
library(Seurat)

# set current working directory
.plavisca_root <- local({
  candidates <- unique(c(
    Sys.getenv("PLAVISCA_PREPROCESS_ROOT", unset = NA_character_),
    getwd(),
    "/home/sopheap/pvsca_b/pre_process_data"
  ))
  candidates <- candidates[!is.na(candidates) & nzchar(candidates)]
  hit <- candidates[file.exists(file.path(candidates, "scripts", "pipeline_lib.R"))]
  if (length(hit) == 0) {
    stop(
      "Could not locate the pre_process_data project root (looked for scripts/",
      "pipeline_lib.R under $PLAVISCA_PREPROCESS_ROOT, the current working ",
      "directory, and the historical hard-coded path). Set the ",
      "PLAVISCA_PREPROCESS_ROOT environment variable to the pre_process_data ",
      "directory's absolute path, or run this script from that directory.",
      call. = FALSE
    )
  }
  normalizePath(hit[[1]])
})
setwd(.plavisca_root)
source("scripts/pipeline_lib.R")

##### --------------------------------------------------------------------------
##### 2. Load the author-derived expression source and build the validated
##### author-cell-id -> PlaViSca-cell-id crosswalk (51inf->612, 52inf->610,
##### 91inf->611, 92inf->609 - see audit/ruberto2022_1/AUDIT.md and
##### author_matrix_validation/cell_crosswalk_validation.tsv).
##### --------------------------------------------------------------------------
study_num <- "36093191"

Hep59.1.2.seu <- readRDS("data/Hep59.1.2.seu_20aug2025.rds")

author_crosswalk <- data.frame(
  author_cell_id = rownames(Hep59.1.2.seu@meta.data),
  liver_form = Hep59.1.2.seu@meta.data$LiverForm
) |>
  mutate(
    cells_name = sub(".*_", "", author_cell_id),
    cells_num = sub("_.*", "", author_cell_id),
    run_suffix = case_when(
      cells_num == "51inf" ~ "612",
      cells_num == "52inf" ~ "610",
      cells_num == "91inf" ~ "611",
      cells_num == "92inf" ~ "609"
    ),
    plavisca_cell_id = paste0(run_suffix, "_", cells_name)
  )

if (nrow(author_crosswalk) != 1438L) {
  pipeline_fail(sprintf("Hep59.1.2.seu has %d cells, expected exactly 1438", nrow(author_crosswalk)))
}
assert_unique_cell_ids(author_crosswalk$author_cell_id, label = "author_crosswalk$author_cell_id")
assert_unique_cell_ids(author_crosswalk$plavisca_cell_id, label = "author_crosswalk$plavisca_cell_id")
if (anyNA(author_crosswalk$run_suffix)) {
  pipeline_fail("author_crosswalk contains an author cell-ID prefix that did not map to a known run (51inf/52inf/91inf/92inf)")
}
rownames(author_crosswalk) <- author_crosswalk$plavisca_cell_id

##### --------------------------------------------------------------------------
##### 3. Per-run metadata (D038,D039,D069 fixes, unchanged from Phase 1): four
##### source libraries/runs (SRR19573609-612), two biological replicates x two
##### treatment arms. These are per-run CONSTANTS, not derived from any count
##### matrix, so they are assigned directly to the final 1,438-cell population
##### via a keyed join on run_suffix - no STARsolo object required to compute
##### them.
##### --------------------------------------------------------------------------
run_table <- data.frame(
  run_suffix = c("609", "610", "611", "612"),
  run_id = c("SRR19573609", "SRR19573610", "SRR19573611", "SRR19573612"),
  study_pmid = study_num,
  study_label = STUDY_LABELS[["ruberto2022_1"]],
  pub_year = 2022,
  parasite_origin_location = "Cambodia_Mondulkiri",
  # D038,D069: corrected spelling; V3 chemistry not independently confirmed
  # for this study, so left unqualified.
  sc_technology = "10x_Chromium",
  sequencer_registered = "HiSeq X Ten",
  host_species = "Homo sapiens",
  host_taxid = 9606L,
  donor_id = "BioIVT:BGW",
  sample_type = "Mammalian host: hepatocyte",
  tissue_or_sample_type = "Host liver",
  strain = "Cambodia field isolate",
  day_post_infection = c(9, 5, 9, 5),
  # D039: "No_Treatment" is already the correct current-script value (the
  # deployed "None" string is stale/historical drift, not reproduced here).
  source_treatment = c("MMV390048", "None", "MMV390048", "None"),
  treatment = c("MMV390048", "No_Treatment", "MMV390048", "No_Treatment"),
  # Four source libraries = two biological replicates x two treatment arms,
  # explicitly identified (not just implied by run order).
  biological_replicate_id = c("2", "2", "1", "1"),
  biological_replicate_type = "infection",
  stringsAsFactors = FALSE
)
assert_keyed_join(author_crosswalk$run_suffix, run_table$run_suffix, label = "run_suffix")

metadata_df <- author_crosswalk |>
  left_join(run_table, by = "run_suffix") |>
  mutate(
    barcode = plavisca_cell_id,
    host_id = donor_id, # retained for backward compatibility with existing app schema
    source_life_cycle_stage = liver_form,
    source_stage_provenance = "source_selection_defined"
  )
rownames(metadata_df) <- metadata_df$plavisca_cell_id

##### --------------------------------------------------------------------------
##### 4. Adopt the complete author-derived raw UMI matrix for ALL 1,438 cells
##### (D036/DEC09 project-lead sign-off). Reindex (not just label) the author
##### RNA counts matrix onto PlaViSca cell IDs, keyed and asserted - never
##### positional.
##### --------------------------------------------------------------------------
author_counts_full <- GetAssayData(Hep59.1.2.seu, assay = "RNA", layer = "counts")
author_counts <- author_counts_full[, author_crosswalk[author_crosswalk$plavisca_cell_id, "author_cell_id"], drop = FALSE]
colnames(author_counts) <- author_crosswalk$plavisca_cell_id

if (any(author_counts@x < 0) || any(author_counts@x != round(author_counts@x))) {
  pipeline_fail("author RNA counts layer contains non-integer or negative values - refusing to adopt as raw UMI counts")
}

# Feature-space handling (Phase 2 spec Part 5/6): zero-fill ONLY the PlaViSca
# features explicitly classified not_part_of_author_reference in
# audit/ruberto2022_1/phase2_feature_coverage.tsv - never a blanket
# zero-pad. That table documents, per feature, why it is absent from the
# author matrix (author-side protein_coding_gene-only reference filter,
# PlasmoDB-51 vs PlasmoDB-68 annotation growth); see the table for evidence.
# The PlaViSca feature panel is taken directly from that coverage table
# (built once, from the pre-Phase-2 ruberto2022_1.rds feature set, which is
# identical across all STARsolo-derived studies' shared reference).
feature_coverage <- read.delim("audit/ruberto2022_1/phase2_feature_coverage.tsv", stringsAsFactors = FALSE)
plavisca_features <- feature_coverage$plavisca_feature
if (!all(rownames(author_counts) %in% plavisca_features)) {
  pipeline_fail("author RNA counts matrix contains feature IDs absent from the PlaViSca feature panel in phase2_feature_coverage.tsv")
}
zero_fill_features <- feature_coverage$plavisca_feature[feature_coverage$classification == "not_part_of_author_reference"]
non_author_features <- setdiff(plavisca_features, rownames(author_counts))
if (!setequal(non_author_features, zero_fill_features)) {
  pipeline_fail("zero-filled feature set does not exactly match the audited not_part_of_author_reference feature list in phase2_feature_coverage.tsv")
}

full_counts <- Matrix::Matrix(
  0,
  nrow = length(plavisca_features), ncol = nrow(author_crosswalk),
  dimnames = list(plavisca_features, author_crosswalk$plavisca_cell_id),
  sparse = TRUE
)
full_counts[rownames(author_counts), ] <- author_counts

##### --------------------------------------------------------------------------
##### 5. Assemble the final Seurat object directly from the adopted counts
##### and the keyed metadata table (no STARsolo-derived shell object).
##### --------------------------------------------------------------------------
pv.combined.all <- CreateSeuratObject(
  counts = full_counts,
  meta.data = metadata_df[colnames(full_counts), setdiff(colnames(metadata_df), c("author_cell_id", "cells_name", "cells_num", "run_suffix", "plavisca_cell_id")), drop = FALSE] # Phase 2 fix: liver_form must be RETAINED - singleR.R Part 5c reads it to set harmonized_life_cycle_stage precedence for this study
)
# The adopted counts are raw UMI counts, not author-normalized values (the
# author RNA$data/SCT slots are explicitly NOT used - see Phase 2 spec Part
# 3). 'data' starts identical to 'counts'; integration.R performs its own
# NormalizeData() from 'counts' regardless.
pv.combined.all <- SetAssayData(pv.combined.all, assay = "RNA", layer = "data", new.data = full_counts)

# Explicit count provenance (Phase 2 spec Part I.4), populated for every one
# of the 1,438 cells (all now share the same adopted source).
pv.combined.all <- set_count_provenance(
  pv.combined.all,
  source_of_counts = "author_processed_object_raw_umi",
  counting_pipeline = "kallisto_bustools",
  count_reference_version = "PlasmoDB-51_PvivaxP01_AnnotatedTranscripts",
  count_provenance_status = "computationally_validated_zenodo_byte_identity_unconfirmed"
)
pv.combined.all$count_source_rds_filename <- "Hep59.1.2.seu_20aug2025.rds"
pv.combined.all$count_source_rds_sha256 <- tryCatch(
  as.character(openssl::sha256(file("data/Hep59.1.2.seu_20aug2025.rds", raw = TRUE))),
  error = function(e) digest::digest(file = "data/Hep59.1.2.seu_20aug2025.rds", algo = "sha256")
)

# total_umi_count/near_empty_expression_flag describe the ADOPTED (author)
# counts - this is the field singleR.R and other downstream QC consult, so
# it must reflect the counts actually used.
adopted_total_umi <- Matrix::colSums(full_counts)
pv.combined.all$total_umi_count <- adopted_total_umi
pv.combined.all$near_empty_expression_flag <- adopted_total_umi <= 1

##### --------------------------------------------------------------------------
##### 6. Historical/diagnostic-only comparison (D036 characterization): reuse
##### the pre-Phase-2 STARsolo-derived per-cell UMI totals for these exact
##### 1,438 cells directly from this study's existing (about-to-be-replaced)
##### output object, rather than re-deriving them from the raw STARsolo
##### barcode space. These fields are NEVER used for QC/eligibility gating
##### downstream (see legacy_starsolo_* notes in metadata_schema.tsv) - only
##### for the Part VI old-vs-new audit comparison. If no prior build exists
##### (e.g. a from-scratch checkout that has never produced ruberto2022_1.rds
##### before), the comparison is explicitly marked unavailable rather than
##### silently recomputed from a different, unaudited source.
##### --------------------------------------------------------------------------
legacy_path <- "ruberto2022_1.rds"
if (file.exists(legacy_path)) {
  legacy <- tryCatch(readRDS(legacy_path), error = function(e) NULL)
  has_legacy_fields <- !is.null(legacy) &&
    all(c("total_umi_count", "near_empty_expression_flag") %in% colnames(legacy@meta.data)) &&
    setequal(colnames(legacy), colnames(pv.combined.all))
  if (has_legacy_fields) {
    legacy_meta <- legacy@meta.data[colnames(pv.combined.all), , drop = FALSE]
    pv.combined.all$legacy_starsolo_total_umi_count <- legacy_meta$total_umi_count
    pv.combined.all$legacy_starsolo_near_empty_flag <- legacy_meta$near_empty_expression_flag
  } else {
    warning("ruberto2022_1.rds exists but is not a recognizable pre-Phase-2 STARsolo build (missing fields or cell-set mismatch) - legacy_starsolo_* fields set to NA")
    pv.combined.all$legacy_starsolo_total_umi_count <- NA_real_
    pv.combined.all$legacy_starsolo_near_empty_flag <- NA
  }
} else {
  warning("No pre-existing ruberto2022_1.rds found - legacy_starsolo_* historical comparison fields set to NA (the frozen D036 record remains in audit/ruberto2022_1/, independent of this field)")
  pv.combined.all$legacy_starsolo_total_umi_count <- NA_real_
  pv.combined.all$legacy_starsolo_near_empty_flag <- NA
}

##### --------------------------------------------------------------------------
##### 7. Final assertions and save
##### --------------------------------------------------------------------------
# Conservative membership policy (unchanged by Phase 2): preserve the exact
# 1,438-cell population, including all 538 previously-flagged
# replicate-2 cells - do not exclude them; they now carry real
# author-derived expression instead of a near-empty artifact.
assert_cell_count(ncol(pv.combined.all), 1438L, label = "ruberto2022_1.rds")
assert_unique_cell_ids(colnames(pv.combined.all), label = "ruberto2022_1.rds colnames")

n_legacy_near_empty <- sum(pv.combined.all$legacy_starsolo_near_empty_flag, na.rm = TRUE)
n_adopted_near_empty <- sum(pv.combined.all$near_empty_expression_flag)
if (!anyNA(pv.combined.all$legacy_starsolo_near_empty_flag) && n_legacy_near_empty != 538L) {
  warning(sprintf(
    "legacy_starsolo_near_empty_flag count is %d, expected 538 per D036 audit evidence - re-check the historical STARsolo record before trusting this build",
    n_legacy_near_empty
  ))
}
cat(sprintf(
  "Phase 2 Ruberto2022_1: %d/%d cells were STARsolo-near-empty (legacy, superseded); %d/%d cells are near-empty under the newly adopted author counts.\n",
  n_legacy_near_empty, ncol(pv.combined.all), n_adopted_near_empty, ncol(pv.combined.all)
))

# Save Seurat object
saveRDS(pv.combined.all, file = "ruberto2022_1.rds")
record_build_manifest(
  artifact_path = "ruberto2022_1.rds",
  script_path = "scripts/ruberto2022_1_pv_analysis_script.R",
  cell_count = ncol(pv.combined.all),
  notes = sprintf(
    "PHASE 2: D036/DEC09 resolved - adopted author kallisto/bustools raw UMI counts (Hep59.1.2.seu_20aug2025.rds) for all 1438 cells, replacing the STARsolo reconstruction; %d/2089 zero-filled reference-absent features; %d cells near-empty under legacy STARsolo (historical), %d cells near-empty under adopted author counts. STARsolo raw barcode space was not re-read (cell membership and metadata are fully determined by the author crosswalk and per-run constants); legacy comparison fields sourced from the prior build's own output.",
    length(zero_fill_features), n_legacy_near_empty, n_adopted_near_empty
  )
)

# Remove all object
rm(list = ls())
gc()
