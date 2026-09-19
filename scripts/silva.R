library(Seurat)
library(dplyr)

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

# ---------------------------------------------------------------------------
# Source metadata (Table S2; preserved verbatim as source_* fields - D008,D009)
# ---------------------------------------------------------------------------
table_s2 <- readxl::read_xls("data/TableS2_Cell_Metadata.xls", skip = 1) |>
  mutate(
    refine_state = case_when(
      State == "Replicative" ~ "Schizont (liver stage)",
      TRUE ~ State
    )
  )

so <- readRDS("data/PvData_CHM_Final.RDS")

raw <- GetAssayData(so, assay = "RNA", layer = "counts")
rownames(raw) <- sub("^([^ -]*-[^ -]*).*", "\\1", rownames(raw))

silva <- CreateSeuratObject(counts = raw)
n_cells <- ncol(silva)
assert_cell_count(n_cells, 1494L, label = "Mancio-Silva2022 source object")

# ---------------------------------------------------------------------------
# Keyed source-group/accession crosswalk (fixes D002,D005: no more
# `rep(paste0("SRR18134", 227:284))` positional recycling; every run/
# sequencer value now traces to a keyed GEO/ENA crosswalk row)
# ---------------------------------------------------------------------------
run_table <- read.delim(
  "scripts/ref/mancio_silva2022_run_table.tsv",
  stringsAsFactors = FALSE
)

raw_colnames <- colnames(silva)
is_dot1_duplicate <- grepl("\\.1$", raw_colnames)
base_colname <- sub("\\.1$", "", raw_colnames)
prefix_key <- sub("_[ACGTN]+$", "", base_colname)

assert_keyed_join(unique(prefix_key), run_table$prefix_key, label = "Mancio-Silva2022 source-group prefix_key")

crosswalk <- run_table[match(prefix_key, run_table$prefix_key), ]

# D001: correct PMID (was silently reusing Hazzard2024's 39223117)
silva$study_pmid <- rep("35443155", n_cells)

# D007: one canonical study_label string
silva$study_label <- rep(STUDY_LABELS[["mancio_silva2022"]], n_cells)
silva$pub_year <- rep(2022L, n_cells)

# D002,D005: run/sequencer provenance from the keyed crosswalk, never
# positional recycling; multi-run groups keep their full run set as a
# documented semicolon-joined string rather than an arbitrary single value.
silva$run_id <- crosswalk$run_accession
silva$sequencer_registered <- crosswalk$sequencer_registered
# D005: publication-protocol conflict for Infection2/SEQ groups is not
# silently resolved - no independent publication-protocol sequencer text was
# recoverable for this study during the audit, so this stays NA (true
# missing, never a fabricated value) pending further provenance review.
silva$sequencer_protocol_publication <- rep(NA_character_, n_cells)

# D003: source day parsed from the validated keyed crosswalk (was: 1222/1494
# cell-level disagreements from a hand-rolled day vector)
silva$day_post_infection <- as.integer(crosswalk$day_post_infection)
assert_no_missing_mandatory(silva$day_post_infection, label = "silva$day_post_infection")

# D004: treatment harmonization distinguishes study_control from real
# No_Treatment; source_treatment (CTRL/PI4K) is preserved verbatim.
silva$source_treatment <- crosswalk$source_treatment
silva$treatment <- case_when(
  crosswalk$source_treatment == "CTRL" ~ "study_control",
  crosswalk$source_treatment == "PI4K" ~ "PI4K_inhibitor",
  TRUE ~ NA_character_
)
assert_no_literal_na_string(silva$treatment, label = "silva$treatment")

# D006: parasite vs host geography kept separate, never concatenated
silva$parasite_origin_location <- rep("Thailand_Ubon-Ratchathani", n_cells)
silva$host_provenance <- rep("US-procured primary human hepatocytes (BioIVT)", n_cells)
silva$experimental_site <- rep(NA_character_, n_cells)
silva$sequencing_site <- rep(NA_character_, n_cells)
silva$registry_location <- rep(NA_character_, n_cells)

silva$sc_technology <- rep("Seq_Well", n_cells)
silva$host_species <- rep("Homo sapiens", n_cells)
silva$host_taxid <- rep(9606L, n_cells)
silva$donor_id <- rep(NA_character_, n_cells) # exact host donor identity unresolved (D006 open provenance question)
silva$biological_replicate_id <- crosswalk$source_group
silva$biological_replicate_type <- rep("infection", n_cells)
silva$tissue_or_sample_type <- rep("Host liver", n_cells)
silva$sample_type <- rep("Mammalian host: hepatocyte", n_cells) # retained for backward compatibility with existing app schema
silva$strain <- rep("Patient isolate", n_cells)

# D013: num_srr removed from cell metadata; run counts belong in study-level
# provenance, not an ambiguous per-cell scalar.
silva$source_run_count_deposited <- rep(58L, n_cells)
silva$source_run_count_imported <- rep(58L, n_cells)
silva$source_run_count_retained <- rep(58L, n_cells)

# ---------------------------------------------------------------------------
# D008,D009,D010: preserve source cell/barcode provenance exactly, and join
# State/AP2G/sex-state by an EXACT keyed match on Table S2's Updated_names
# (all 1494 raw colnames match exactly 1:1; the 14 proven ".1" duplicate
# records have no Table S2 match by design - see DEC01) rather than the
# previous fragile grep-prefix matching.
# ---------------------------------------------------------------------------
silva$source_cell_id <- raw_colnames
silva$reconstructed_barcode <- sub("^[^_]+_", "", base_colname)

s2_match <- match(raw_colnames, table_s2$Updated_names)
silva$source_state <- table_s2$State[s2_match]
silva$source_ap2g <- table_s2$AP2G[s2_match]
silva$source_sex_state <- table_s2$sex.state[s2_match]
silva$source_orig_ident <- table_s2$orig.ident[s2_match]
silva$refine_state <- table_s2$refine_state[s2_match]

# ---------------------------------------------------------------------------
# D011,D012: explicit duplicate/provenance flags (Phase 1 conservative
# policy - do NOT remove any cells; preserve flags for the approved
# Phase-2 policy per DEC01/DEC02)
# ---------------------------------------------------------------------------
silva$dot1_duplicate_flag <- is_dot1_duplicate
silva$dot1_duplicate_of <- ifelse(is_dot1_duplicate, base_colname, NA_character_)

d5_cross_array_pair <- prefix_key %in% c("D5Seq1", "D5Seq2")
silva$d5_cross_array_duplicate_flag <- d5_cross_array_pair

silva$source_stage_provenance <- ifelse(
  is.na(silva$refine_state),
  "singleR_inferred",
  "source_state_defined"
)

# Add barcode column for consistency with other studies (kept distinct from
# the raw reconstructed_barcode - see metadata_schema.tsv)
silva$barcode <- colnames(silva)

# Conservative Phase 1 membership policy: preserve all 1494 cells, including
# the 14 proven .1 duplicates and both members of the unresolved D5 pair.
assert_cell_count(ncol(silva), 1494L, label = "silva2022.rds")
assert_unique_cell_ids(colnames(silva), label = "silva2022.rds colnames")
assert_no_literal_na_string(silva$study_pmid, label = "silva$study_pmid")

saveRDS(silva, file = "silva2022.rds")
record_build_manifest(
  artifact_path = "silva2022.rds",
  script_path = "scripts/silva.R",
  cell_count = ncol(silva),
  notes = "Fixes D001-D010,D013; D011/D012 flagged not removed per Phase 1 conservative membership policy (DEC01/DEC02 pending)"
)
