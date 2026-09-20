library(tidyverse)
library(Seurat)
# DropletUtils/rtracklayer/scater/SingleCellExperiment are only needed by the
# slow per-chunk STARsolo-reprocessing fallback path below - loaded lazily.

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
##### PHASE 2 efficiency fix (2026-09-20, at user request): expression
##### counts/cell-membership are unchanged by Phase 1/2 for this study (only
##### metadata was buggy - D026-D030,D032,D033,D069, all fixed via the keyed
##### run table below); Part II of the Phase 2 spec requires non-Ruberto2022_1
##### studies' expression matrices to remain consistent with their audited
##### production source. When a prior hazzard2024.rds with the expected
##### 80,024-cell population already exists, its RNA counts and cell
##### membership are reused directly and ONLY the metadata is recomputed
##### fresh from scripts/ref/hazzard2024_run_table.tsv (a keyed join on
##### run_id, derived from each cell's `<3-digit-SRR-suffix>_<barcode>`
##### colname - run_id = paste0("SRR27021", suffix), confirmed against every
##### suffix actually present). This is bit-identical to a full rebuild for
##### counts/membership and skips re-running all 6 STARsolo chunk scripts +
##### their re-merge entirely. A from-scratch checkout with no prior artifact
##### (and no chunk .rds files either) fails loudly, instructing the operator
##### to run the chunk scripts first.
##### --------------------------------------------------------------------------
run_table <- read.delim(
  "scripts/ref/hazzard2024_run_table.tsv",
  stringsAsFactors = FALSE
)

existing_path <- "hazzard2024.rds"
reused <- NULL
if (file.exists(existing_path)) {
  obj <- tryCatch(readRDS(existing_path), error = function(e) NULL)
  if (!is.null(obj) && inherits(obj, "Seurat") && ncol(obj) == 80024L) {
    suf <- sub("^([0-9]{3})_.*$", "\\1", colnames(obj))
    run_id_per_cell <- paste0("SRR27021", suf)
    if (!anyNA(suf) && all(run_id_per_cell %in% run_table$run_id)) {
      reused <- obj
    }
  }
}

if (!is.null(reused)) {
  message("hazzard2024: reused expression counts and cell membership from existing hazzard2024.rds (skipped all 6 STARsolo chunk scripts); recomputing metadata from the current keyed run-table logic only.")

  # Build the ENTIRE metadata data.frame first and pass it to
  # CreateSeuratObject() in one call, rather than looping ~20+ separate
  # `pv.combined.all[[f]] <-` assignments after construction - at this
  # study's 80,024-cell scale, Seurat's per-column metadata setter was
  # observed to be prohibitively slow (>15 min, still not finished); a
  # single meta.data= construction is the same pattern already used
  # successfully for Ruberto2022_1 and completes in seconds.
  reused_counts <- GetAssayData(reused, assay = "RNA", layer = "counts")
  cell_ids <- colnames(reused_counts)
  suf <- sub("^([0-9]{3})_.*$", "\\1", cell_ids)
  run_id_per_cell <- paste0("SRR27021", suf)
  assert_keyed_join(unique(run_id_per_cell), run_table$run_id, label = "hazzard2024 fast-path run_id")

  meta_df <- data.frame(run_id = run_id_per_cell, stringsAsFactors = FALSE) |>
    left_join(run_table, by = "run_id") |>
    mutate(
      study_pmid = "39223117",
      study_label = STUDY_LABELS[["hazzard2024"]],
      pub_year = 2024L,
      sc_technology = technology,
      sequencer = sequencer_registered,
      # No drug-treatment arm exists in this study (infection design - mono/
      # consecutive/simultaneous/sporozoite - is a separate covariate
      # carried in infection_design, not a treatment arm); "No_Treatment"
      # is genuinely correct, never a literal "NA" string.
      treatment = "No_Treatment",
      source_treatment = NA_character_,
      host_id = library_id,
      sample_type = "Host blood",
      # Pre-existing Phase 1 gap (present in every prior version of this
      # script and the per-chunk STARsolo scripts, not introduced in Phase
      # 2): only `sample_type` was ever set, never the canonical
      # `tissue_or_sample_type` field singleR.R's IDC eligibility gate
      # (Part 1) actually reads. Without it, all 80,024 Hazzard2024 cells
      # silently fail `tissue_or_sample_type == "Host blood"` (NA, not
      # TRUE) and are wrongly excluded from blood-stage IDC/gametocyte
      # classification - discovered when the Phase 2 atlas rebuild's
      # eligible-cell count matched Sa2020+Hazzard2022-blood exactly (10202)
      # with zero Hazzard2024 cells included. Fixed here for all cells;
      # this is the correct TISSUE_OR_SAMPLE_TYPE-vocabulary value
      # (pipeline_lib.R) since every Hazzard2024 library is a Saimiri
      # boliviensis peripheral blood draw (see host_provenance/host_species
      # in scripts/ref/hazzard2024_run_table.tsv).
      tissue_or_sample_type = "Host blood",
      retained_by_plavisca_qc = TRUE,
      barcode = cell_ids
    )
  rownames(meta_df) <- cell_ids

  pv.combined.all <- CreateSeuratObject(
    counts = reused_counts,
    meta.data = meta_df,
    min.cells = 0, min.features = 0
  )
} else {
  message("hazzard2024: no reusable prior hazzard2024.rds found - falling back to the full per-chunk STARsolo reprocessing path (requires hazzard2024_pv_combined_*.rds chunk files, produced by running each hazzard2024_pv_combined_*.R script first).")

  rds_files <- list.files(pattern = "^hazzard2024_pv.*\\.rds$")
  if (length(rds_files) == 0) {
    pipeline_fail("no hazzard2024_pv_combined_*.rds chunk files found and no reusable hazzard2024.rds exists - run the 6 STARsolo chunk scripts (hazzard2024_pv_combined_*.R) first")
  }
  hazzard2024_list <- lapply(rds_files, readRDS)
  pv.combined.all <- Reduce(function(x, y) merge(x, y = y), hazzard2024_list)
  pv.combined.all <- JoinLayers(pv.combined.all)

  # read process data from study hazzard2024
  df <- read.delim("data/Proccessed_Data.txt") |>
    rownames_to_column("barcode") |>
    mutate(
      prefix = sub("(_[^_]+)$", "", barcode),
      cells = sub("^[^_]+_[^_]+_", "", barcode)
    ) |>
    select(prefix, cells) |>
    mutate(
      prefix = case_when(
        prefix == "4881_1" ~ "994",
        prefix == "5142_1" ~ "992",
        prefix == "5142_2" ~ "991",
        prefix == "5163_1" ~ "981",
        prefix == "5163_2" ~ "980",
        prefix == "5164_1" ~ "976",
        prefix == "5164_2" ~ "975",
        prefix == "5309_1" ~ "956",
        prefix == "5309_2" ~ "955",
        prefix == "5350_1" ~ "972",
        prefix == "5350_2" ~ "970",
        prefix == "5350_3" ~ "969",
        prefix == "5370_1" ~ "968",
        prefix == "5370_2" ~ "967",
        prefix == "5370_3" ~ "966",
        prefix == "5537_1" ~ "987",
        prefix == "5537_2" ~ "986",
        prefix == "5537_3" ~ "985",
        prefix == "5537_4" ~ "984",
        prefix == "5550_1" ~ "983",
        prefix == "5708_1" ~ "993",
        prefix == "5708_2" ~ "982",
        prefix == "5708_3" ~ "971",
        prefix == "5708_4" ~ "960",
      ),
      new_cells = paste0(prefix, "_", cells)
    )

  assert_no_missing_mandatory(df$new_cells, label = "df$new_cells (retained-cell crosswalk)")
  assert_unique_cell_ids(df$new_cells, label = "df$new_cells (retained-cell crosswalk)")

  pv.combined.all <- subset(
    pv.combined.all,
    cells = df$new_cells
  )
}

# ============================================================================
# DEC07 (2026-09-19 UPDATE, "interim state" recommendation): the authors'
# own PseudoGroup classification (data/Proccessed_Data.txt) must be stored
# verbatim as source_pseudogroup, distinct from and never overwritten by
# PlaViSca's own cluster/marker-derived pred_gametocyte_sex - this field was
# previously never imported at all (a genuine Phase 1 gap, not merely a
# rename). "Sexual" cells (21,009) are the authors' own sexual-stage call;
# GroupA/B/C are their asexual-stage sub-groupings. Neither label is
# preferred by default; reconciliation with pred_gametocyte_sex is deferred
# to adjudicated_display_stage (DEC07 (c), unchanged).
# ============================================================================
pseudogroup_crosswalk <- read.delim("data/Proccessed_Data.txt") |>
  rownames_to_column("barcode") |>
  transmute(
    prefix = sub("(_[^_]+)$", "", barcode),
    cells = sub("^[^_]+_[^_]+_", "", barcode),
    PseudoGroup
  ) |>
  mutate(
    prefix = case_when(
      prefix == "4881_1" ~ "994",
      prefix == "5142_1" ~ "992",
      prefix == "5142_2" ~ "991",
      prefix == "5163_1" ~ "981",
      prefix == "5163_2" ~ "980",
      prefix == "5164_1" ~ "976",
      prefix == "5164_2" ~ "975",
      prefix == "5309_1" ~ "956",
      prefix == "5309_2" ~ "955",
      prefix == "5350_1" ~ "972",
      prefix == "5350_2" ~ "970",
      prefix == "5350_3" ~ "969",
      prefix == "5370_1" ~ "968",
      prefix == "5370_2" ~ "967",
      prefix == "5370_3" ~ "966",
      prefix == "5537_1" ~ "987",
      prefix == "5537_2" ~ "986",
      prefix == "5537_3" ~ "985",
      prefix == "5537_4" ~ "984",
      prefix == "5550_1" ~ "983",
      prefix == "5708_1" ~ "993",
      prefix == "5708_2" ~ "982",
      prefix == "5708_3" ~ "971",
      prefix == "5708_4" ~ "960",
    ),
    new_cells = paste0(prefix, "_", cells)
  ) |>
  select(new_cells, PseudoGroup)

pseudogroup_match <- match(colnames(pv.combined.all), pseudogroup_crosswalk$new_cells)
pv.combined.all$source_pseudogroup <- pseudogroup_crosswalk$PseudoGroup[pseudogroup_match]
n_pseudogroup_matched <- sum(!is.na(pv.combined.all$source_pseudogroup))
cat(sprintf(
  "[hazzard2024_merge_all.R] source_pseudogroup imported for %d/%d cells (%s)\n",
  n_pseudogroup_matched, ncol(pv.combined.all),
  paste(names(table(pv.combined.all$source_pseudogroup)), table(pv.combined.all$source_pseudogroup), sep = "=", collapse = ", ")
))

# Explicit count provenance (Phase 2 spec Part I.4).
# QC UMI fields required by singleR.R's IDC eligibility gate (Part IV) -
# must be populated for every study, not only Ruberto2022_1.
pv.combined.all <- set_qc_umi_fields(pv.combined.all)

pv.combined.all <- set_count_provenance(
  pv.combined.all,
  source_of_counts = "plavisca_starsolo_raw_umi",
  counting_pipeline = "STARsolo",
  count_reference_version = "PlasmoDB-68_PvivaxP01",
  count_provenance_status = "plavisca_production_pipeline"
)

# Conservative cell-membership policy for Phase 1: preserve the exact
# 80,024-cell Hazzard2024 inclusion population.
assert_cell_count(ncol(pv.combined.all), 80024L, label = "hazzard2024.rds")
assert_unique_cell_ids(colnames(pv.combined.all), label = "hazzard2024.rds colnames")
assert_no_literal_na_string(pv.combined.all$library_id, label = "hazzard2024$library_id")
assert_no_literal_na_string(pv.combined.all$day_post_infection, label = "hazzard2024$day_post_infection")

# Save the final object
saveRDS(pv.combined.all, file = "hazzard2024.rds")
record_build_manifest(
  artifact_path = "hazzard2024.rds",
  script_path = "scripts/hazzard2024_merge_all.R",
  cell_count = ncol(pv.combined.all),
  notes = "Fixes D026-D030,D032,D033,D069 (keyed run table); preserves exact 80,024-cell population per Phase 1 conservative membership policy"
)

# Cleanup
rm(list = ls())
gc()
