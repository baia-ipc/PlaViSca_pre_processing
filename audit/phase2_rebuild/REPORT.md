# PlaViSca Phase-2 Rebuild — Candidate Atlas Report

**Date:** 2026-09-20
**Build-start git commit:** `c47791342bf47a035630c7606c4dfa1102fb4385` ("Validate Ruberto2022_1 author count matrix")
**Branch:** `plavisca-production-repair`
**Repository:** `pre_process_data` (this task does not touch the `PlaViSca` application repository)

This report documents the Phase-2 rebuild: adoption of the validated author-derived
raw UMI matrix for Ruberto2022_1 (all 1,438 cells), a full six-study atlas
regeneration, Harmony integration, SingleR/gametocyte re-annotation, candidate
app-facing exports, and full acceptance-test validation.

---

## 1. Part I — Ruberto2022_1 author-matrix adoption

**Decision (DEC09, 2026-09-20 update):** use the complete author-derived raw UMI
matrix (`data/Hep59.1.2.seu_20aug2025.rds`, `RNA` counts layer) for **all 1,438**
Ruberto2022_1 cells, not just the 538 previously near-empty (D036) cells. Recorded
in `audit/shared_pipeline/decision_register.tsv` (DEC09) and
`audit/ruberto2022_1/D036_resolution.md` ("PHASE 2 UPDATE" section), without
deleting or rewriting the prior forensic history.

**Implementation:** `scripts/ruberto2022_1_pv_analysis_script.R` was rewritten to
build the crosswalk directly from the author object (no STARsolo raw-matrix I/O
needed — see Part 1.5 efficiency note below), substitute the author `RNA` counts
for all 1,438 cells, and zero-fill only the 2,089 PlaViSca features classified
`not_part_of_author_reference` (see Section 2 below). Explicit provenance fields
(`source_of_counts`, `counting_pipeline`, `count_reference_version`,
`count_provenance_status`, `count_source_rds_filename`, `count_source_rds_sha256`)
are populated for all 1,438 cells.

**Raw-count equality validation:** exact equality confirmed between the adopted
counts and the author matrix over the shared 4,722 features x 1,438 cells
(6,790,236 entries, **0 mismatches**). Zero-fill confirmed limited exactly to the
2,089 audited reference-absent features (100% zero, no unexpected nonzero rows
outside the author-provided feature set).

**Result:** `near_empty_expression_flag` (computed from the adopted counts) is
**0/1438** — the D036 artifact is fully resolved. `legacy_starsolo_near_empty_flag`
(538/1438) is retained as a frozen historical field, and the underlying 538-row
STARsolo evidence remains independently documented in
`audit/ruberto2022_1/near_empty_cells_vs_hep59.tsv` (539 lines incl. header) —
verified by the `D036-flag-legacy` test.

### 1.5 Efficiency note (implemented at user's explicit request mid-task)

The original Phase-1 Ruberto2022_1 script (and the other four STARsolo-derived
study scripts) re-read STARsolo's *unfiltered* raw barcode matrices (millions of
candidate barcodes per run) purely to reconstruct an intermediate object that was
then subset down to a few thousand real cells — expensive and, for
Ruberto2022_1/Sa2020/Hazzard2022/Ruberto2022_2, entirely unnecessary once a prior
correctly-membership build already exists. All five scripts (plus
`hazzard2024_merge_all.R`) were given an explicit "reuse existing counts, only
recompute metadata" fast path: when a prior study `.rds` with the exact expected
cell/feature structure exists, its RNA counts are reused directly (bit-identical
to a fresh reprocessing, since Phase 1/2 do not change these studies' count-
generation logic) and only metadata is recomputed from the current (corrected)
script logic. Falls back to the original full STARsolo path if no reusable
artifact exists (e.g. a from-scratch checkout). This reduced per-study
regeneration time from tens of minutes (Ruberto2022_1's first, killed run: ~40
min still incomplete) to 18-63 seconds each.

---

## 2. Part I.5 — Feature-space handling

`audit/ruberto2022_1/phase2_feature_coverage.tsv` (6,811 rows) classifies every
PlaViSca Ruberto2022_1 feature. Evidence: the authors' own workflow
(`audit/ruberto2022_1/ruberto_github_additional_file_1.Rmd`/`_2.Rmd`) downloads
`PlasmoDB-51_PvivaxP01.gff` and builds `gene.info` via
`filter(type == "protein_coding_gene")` — the identical filter logic PlaViSca's
own script uses, just against an older annotation release (PlasmoDB-51 vs
PlaViSca's PlasmoDB-68).

| classification | n features | evidence |
|---|---|---|
| `measured_and_present` | 4,722 | present in the author RNA counts assay, confirmed exact match |
| `not_part_of_author_reference` (biotype-excluded) | 288 (103 ncRNA_gene + 185 pseudogene) | author's `protein_coding_gene`-only filter categorically excludes these biotypes, independent of GFF release version — confirmed 0/288 present in author matrix |
| `not_part_of_author_reference` (release-version gap) | 1,801 (protein_coding_gene) | present in PlasmoDB-68 but absent from the author's PlasmoDB-51-derived gene set; author applied no expression-based gene filtering (`CreateSeuratObject()` calls use no `min.cells`) |

No feature was zero-filled without this evidence trail (Part 5 requirement: never
treat reference-absent genes as unconditional biological zero).

---

## 3. Part II — Six-study regeneration and exact cell counts

| Study | Expected | Actual | Match |
|---|---|---|---|
| Mancio-Silva2022 | 1,494 | 1,494 | yes |
| Sa2020 | 9,766 | 9,766 | yes |
| Hazzard2022 | 3,294 | 3,294 | yes |
| Hazzard2024 | 80,024 | 80,024 | yes |
| Ruberto2022_1 | 1,438 | 1,438 | yes |
| Ruberto2022_2 | 9,947 | 9,947 | yes |
| **Total** | **105,963** | **105,963** | yes |

No unexpected cell membership changes. See `study_cell_counts.tsv` and
`study_expression_qc.tsv` for per-study QC (median/mean/min/max UMI, feature
count, count provenance, integer/non-negative confirmation).

### Additional gaps found and fixed during Phase 2

All pre-existing Phase-1 defects, never previously exercised because this is the
first full end-to-end atlas run. Documented here for full transparency, fixed in
the relevant production scripts, and (where the already-completed Harmony/SingleR
compute made a full re-run costly) additionally verified via a metadata-only
patch onto the already-integrated/annotated atlas, using the exact same
deterministic logic the corrected script would have applied:

1. **`total_umi_count` never populated for Sa2020/Hazzard2022/Ruberto2022_2**
   (the shared `set_qc_umi_fields()` helper was added to `pipeline_lib.R` after
   these three studies were first regenerated, and they were not re-run before
   this was discovered via `singleR.R`'s IDC eligibility gate returning `NA`).
   Fixed: all three re-run; atlas patched.
2. **`tissue_or_sample_type` never set for Hazzard2024** (only the non-canonical
   `sample_type` field was ever populated) — silently excluded all 80,024
   Hazzard2024 cells from blood-stage IDC/gametocyte eligibility. Fixed in
   `hazzard2024_merge_all.R` and all 6 STARsolo chunk scripts (for
   from-scratch-checkout correctness); atlas patched. Eligible-cell count went
   from 10,202 (Sa2020 + Hazzard2022-blood only) to the correct 90,219.
3. **`liver_form` column accidentally dropped** from the rewritten
   Ruberto2022_1 script's final object (an artifact of the STARsolo-avoidance
   rewrite in this task, not a pre-existing Phase-1 defect) — broke
   `singleR.R`'s Part 5c liver-stage annotation precedence. Fixed; atlas patched.
4. **`"Schizont (liver stage)"` (lowercase) vs `"Schizont (Liver stage)"`**
   case mismatch between `silva.R`'s `refine_state` values and
   `pipeline_lib.R`'s `HARMONIZED_TO_BROAD_STAGE` mapping table — crashed
   `map_broad_stage()` the first time a Mancio-Silva2022 cell actually reached
   it. Fixed; atlas patched.
5. **Bare `"Schizont"` heuristic-bin string** in `singleR.R`'s IDC hour->bin
   `case_when` (should have been `"Schizont (Blood stage)"` to match the
   mapping table's key) — crashed `map_broad_stage()` for any cell landing in
   the 30-46h bin. Fixed.
6. **Hazzard2022's 2,858 Anopheles sporozoite cells (runs 498/499) never
   received `source_life_cycle_stage`/`source_stage_provenance`**, despite
   `annotation_precedence.tsv` explicitly documenting that they should — these
   cells fell through the annotation precedence chain with `NA`
   `harmonized_life_cycle_stage`. Fixed in `hazzard2022_pv_analysis_script.R`;
   atlas patched. `parasite_broad_stage` NA count dropped from 2,879 to 21
   (the residual 21 are 14 Mancio-Silva `.1`-duplicate cells with no Table S2
   match, by design, plus 7 near-zero-expression cells correctly excluded).
7. **DEC07 (`source_pseudogroup`) never implemented**: the authors' own
   Hazzard2024 `PseudoGroup` classification (`data/Proccessed_Data.txt`) was
   never imported at all, despite `decision_register.tsv` DEC07's explicit
   "interim state" recommendation to store it verbatim, separate from PlaViSca's
   inferred `pred_gametocyte_sex`. Implemented in `hazzard2024_merge_all.R`;
   100% of 80,024 cells matched (GroupA 39,444 / GroupB 12,422 / GroupC 7,149 /
   Sexual 21,009).
8. **`singleR.R`'s own AT21 self-check used `identical()`** on
   `parasite_broad_stage` immediately after a `subset()` call, which can change
   vector attributes (names) without changing values — a false-positive-prone
   check that never actually ran successfully before. Rewritten to compare
   values position-by-position, correctly ignoring cosmetic attribute changes.

None of these fixes altered cell membership or raw expression counts; all are
strictly metadata/annotation-logic corrections, verified against the acceptance
test suite (Section 7) both before and after.

---

## 4. Part III — Integration feature universe (Part 6 of the spec)

| Quantity | Value |
|---|---|
| Atlas feature universe (union, 6 studies) | 6,860 |
| Integration-eligible common genes (intersection, all 6 studies) | 5,203 |
| Genes lost to reference incompatibility (union - intersection) | 1,657 |
| Highly variable genes selected (from the 5,203-gene intersection) | 2,000 |
| Genes contributed to PCA | 2,000 (the above HVGs) |

HVG selection and PCA/Harmony were computed **only** from the 5,203-gene
intersection (a separate subsetted copy), never from the full union with
structural zero-padding — avoiding exactly the "encode study/reference identity
into the matrix" failure mode Part 6 warns against. The main merged object
retains its full 6,860-feature union for per-study raw/normalized exports.
`feature_universe.tsv` (per-feature, per-study measured/not-measured) and
`feature_availability_by_study.tsv` (long format) preserve the
zero-expression-vs-not-measured distinction for downstream (future Phase 3 app)
use.

---

## 5. Part III — Harmony integration validation

Harmony integration succeeded (`harmony` R package was missing from the pixi
environment — added to `scripts/pixi.toml` and installed). Layer->study
diagnostics confirmed exactly 6 layers, each mapping to exactly one study, all 6
studies represented, before every Harmony run.

**Batch mixing within comparable (same broad-stage) populations**
(`batch_mixing_within_comparable_populations.tsv`) — mean fraction of each
cell's 30 nearest neighbors sharing its own study, unintegrated vs integrated:

| Broad stage | n cells | n studies | Unintegrated | Integrated |
|---|---|---|---|---|
| Sporozoite stage | 12,805 | 2 | 0.844 | 0.659 |
| Blood stage | 90,219 | 3 | 0.996 | 0.919 |
| Liver stage | 2,918 | 2 | 0.972 | 0.845 |

Same-study neighbor fraction **decreased** after integration in every
biologically comparable population — the expected signature of a working batch
correction (more cross-study mixing among cells of the same real biological
type). Liver/blood/sporozoite populations were evaluated **separately**, never
pooled together as if their separation were itself a batch effect (Part 13's
explicit instruction).

**Marker gene preservation** (`marker_gene_preservation.tsv`): gametocyte marker
genes (2 female + 2 male, subset shown) remain near-zero in Sporozoite-stage
cells post-integration (e.g. `PVP01_1207200` mean 0.0049) while showing real
signal in Blood/Liver stage cells (0.06-1.0) — expression biology preserved
through integration, not homogenized away.

**Minority-study retention** (`per_study_neighborhood_composition.tsv`):
per-study integrated-space same-study neighbor fraction ranges 0.636 (Hazzard2022,
the smallest study) to 0.998 (Hazzard2024, ~75% of the atlas) — every minority
study retains a distinguishable neighborhood signature; none is diluted to
background/random composition (`retention_ratio_integrated` column quantifies
this against the null expectation for each study's overall atlas share).

Full per-cell reduction data underlying these summaries: `pca_unintegrated`,
`umap_unintegrated`, `pca_integrated`, `umap_integrated`, `tsne_integrated` on
`pv_all_studies.rds`.

---

## 6. Part IV — SingleR / IDC eligibility and gametocyte classification

| Quantity | Value |
|---|---|
| IDC-eligible cells (Host blood, not source-selection-defined stage, >=10 UMI) | 90,219 / 105,963 |
| Cells receiving `idc_reference_similarity_label` | 90,219 (all and only the eligible set) |
| Non-blood cells with non-NA IDC label | 0 (AT08) |
| Cells called Male/Female gametocyte | 23,689 (22,705 Female + 984 Male) |
| Liver-stage cells with a gametocyte-sex call | 0 (AT10) |
| Near-zero-expression cells (<10 UMI) with any expression-derived annotation | 0 (AT11) |
| Sa2020 source gametocyte-sex import | **1,532** matched cells (exact match to the spec's AT13 requirement) |

All 7 checks in `annotation_validation.tsv` PASS; `impossible_state_check.tsv`
contains **0 rows** (0 impossible lifecycle/annotation combinations found).

`idc_reference_similarity_label` (renamed from `hour_post_invasion`),
`idc_similarity_raw_label`, `idc_similarity_score_delta`,
`idc_similarity_pruned_flag`, and `idc_heuristic_stage_bin` (explicitly labeled
PlaViSca-heuristic, not literature-sourced) are all persisted per cell.
`source_sex_annotation` (Sa2020) and `source_pseudogroup` (Hazzard2024) are
retained as source annotations, kept structurally separate from
`pred_gametocyte_sex` (PlaViSca-inferred); `adjudicated_display_stage` remains
NA pending a future reviewed reconciliation decision (DEC07).

Post-integration removal filter (D045/D049 policy): **0 cells removed** — no
Anopheles-host cell was ever classified as a gametocyte by the downstream
classifier without also being source-selection-defined.

---

## 7. Part VI — Ruberto2022_1: old vs. rebuilt positions

`ruberto_old_vs_rebuilt.tsv` compares the deployed app's integrated-UMAP
coordinates against the Phase-2 rebuilt coordinates for all 1,438 Ruberto2022_1
cells, split into the three groups the spec requires (517 formerly-zero-count +
21 formerly-one-count + 900 unaffected).

**Old (deployed) integrated UMAP** — standard deviation of coordinates:

| Group | sd(dim1) | sd(dim2) | sd(dim3) |
|---|---|---|---|
| 538 D036 cells | 0.090 | 0.122 | 0.184 |
| 900 unaffected cells | 0.894 | 1.403 | 2.988 |

The 538 D036 cells are ~5-16x tighter than the unaffected population — the
"artificial island" the task describes, consistent with near-identical
zero-count transcriptomes clustering together artificially.

**New (Phase 2 rebuilt) integrated UMAP:**

| Group | sd(dim1) | sd(dim2) | sd(dim3) |
|---|---|---|---|
| 538 D036 cells | 1.395 | 4.206 | 2.704 |
| 900 unaffected cells | 1.396 | 4.576 | 2.762 |

The 538 cells' spread is now essentially indistinguishable from the unaffected
population's spread — the artificial island has dissolved; these cells now
distribute according to real author-derived expression, exactly as Part 20
requires ("The rebuilt positions must now arise from real author-derived
expression... do not require them to occupy any predetermined location").

---

## 8. Part VII — Acceptance tests

`acceptance_test_results.tsv` (41 rows, generated directly from the test suite's
final run). **41 PASS, 0 FAIL, 0 SKIP.**

Includes: all UNIT tests (7), all STUDY tests (18, 3 per study x 6 studies), the
D036 legacy/resolution/rescue-evidence tests (3), and all atlas-dependent tests
that were previously SKIPPED under Phase 1 and now run for real: AT07, AT08,
AT09, AT10, AT11, AT13, AT17, AT18, AT21 (atlas), AT22, AT24-AT27 (via the
integration-validation table), AT29-AT30 (via the reproducibility-check table),
plus the impossible-state check (13 tests newly activated).

No test was converted to a pass without actually exercising its intended check —
each atlas-dependent test reads real fields from the built `pv_all_studies.rds`
or candidate export files and asserts the specific invariant it names.

---

## 9. Part IX — Build manifest

`build_manifest.tsv` (repo root, 19 rows) records, for every artifact produced
in this task: date, artifact path, output SHA256, script path, git commit
(`c477913...`), R version (4.3.3), Seurat (5.3.0), SeuratObject (5.2.0), harmony
(1.2.3), SingleR (2.4.0), cell count, and notes. `integration.R`'s row records
`seed = 123` (the `RunTSNE` seed; PCA/UMAP use Seurat/UWOT defaults, not
independently seeded in the production scripts).

---

## 10. Part X — Deployed vs. candidate comparison

`deployed_vs_candidate.tsv` — key differences, all classified:

- Total cells: 105,963 both (deployed already had this total; the difference is
  entirely in which 538 Ruberto2022_1 cells' *expression* is real vs.
  near-empty, not in overall membership) — **EXPECTED RECOMPUTATION**.
- Per-study cell counts: unchanged for all 6 studies — **EXPECTED RECOMPUTATION
  (membership unchanged)**, confirming no unapproved membership drift.
- Metadata field count: 52 (deployed) -> 88 (candidate) — **EXPECTED REPAIR**
  (schema renames + added provenance fields per Parts I.4/IV).
- Ruberto2022_1 near-empty cells: 538 -> 0 — **EXPECTED REPAIR** (the core
  purpose of this rebuild).
- Cells with an IDC/HPI label: 105,963 (deployed, assigned to *every* cell
  including liver/sporozoite — the D058 defect) -> 90,219 (candidate, correctly
  gated to blood-eligible cells) — **EXPECTED REPAIR**.
- Cells with a gametocyte call: 105,963 (deployed) -> 23,689 (candidate,
  correctly gated and marker-validated) — **EXPECTED REPAIR/RECOMPUTATION**.

No **UNEXPECTED - INVESTIGATE** or **REGRESSION** classifications were assigned.
Ordinary UMAP/PCA coordinate changes are not flagged as defects (per Part 25's
explicit instruction).

---

## 11. Candidate output paths and hashes

All large candidate artifacts live inside `pre_process_data/` (git-ignored,
`*.rds` never committed) — never inside the `PlaViSca` application repository,
and the deployed `PlaViSca/data/*.rds` files were never written to:

| Artifact | Path | Size (bytes) |
|---|---|---|
| Merged/integrated/annotated atlas | `pre_process_data/pv_all_studies.rds` | 639,532,902 |
| Candidate normalized export | `pre_process_data/data/candidate_export/normalize_df.rds` | 489,254,038 |
| Candidate raw export | `pre_process_data/data/candidate_export/raw_df.rds` | 109,371,235 |
| Candidate scaled export | `pre_process_data/data/candidate_export/scale_df.rds` | 205,718,696 |
| Candidate cleaned dataset | `pre_process_data/data/candidate_export/cleaned_dataset.rds` | 19,357,522 |
| Per-study objects | `pre_process_data/{silva2022,sa2020,hazzard2022,hazzard2024,ruberto2022_1,ruberto2022_2}.rds` | see `study_cell_counts.tsv` |

SHA256 checksums for every artifact: `build_manifest.tsv` (`output_sha256`
column).

---

## 12. Remaining scientific/documentary limitations

1. **Byte-level Zenodo provenance for `Hep59.1.2.seu_20aug2025.rds` remains
   unconfirmed** — recorded explicitly via `count_provenance_status =
   "computationally_validated_zenodo_byte_identity_unconfirmed"` on every
   Ruberto2022_1 cell, not silently treated as fully verified.
2. **The PlasmoDB-51 vs PlasmoDB-68 gene-set difference (1,801 genes) is
   evidence-based inference, not gene-by-gene confirmed** against the actual
   historical PlasmoDB-51 GFF (not available locally or via network in this
   environment) — flagged `not_part_of_author_reference` with an explicit
   "inferred" note in `phase2_feature_coverage.tsv`, not silently treated as
   confirmed.
3. **`adjudicated_display_stage` remains unpopulated** (DEC07 deferred) — the
   Sa2020 (188-cell) and Hazzard2024 (~2,387-cell, via `source_pseudogroup`)
   source-vs-inferred sex/stage disagreements are now both fully retrievable
   (source and inferred fields both present and distinct) but not yet
   reconciled into a single adjudicated display field; this is explicitly a
   future decision, not a Phase-2 blocker.
4. **App data contract changes (Part IX of `app_data_contract.tsv`)** — field
   renames (`development_phase`->`broad_stage`, `parasite_stages`->
   `detailed_stage`, HPI relabeling, gametocyte color-symmetry) are
   Phase-3/app-repair items, out of scope for this preprocessing-only task and
   not implemented here (per explicit task instruction not to touch the
   `PlaViSca` app repository).
5. **Mancio-Silva2022's 14 `.1`-duplicate cells and the D5 cross-array pair**
   remain included, flagged but unresolved (DEC01/DEC02 still open) — no
   membership change was made, per the conservative Phase-1/2 policy.

---

## 13. PHASE-2 RELEASE CANDIDATE determination

All promotion conditions (Part XI) are met:

- All acceptance tests pass (41/41 runnable; 0 FAIL; 0 SKIP).
- Ruberto2022_1 author counts correctly installed for all 1,438 cells, with
  exact raw-count equality confirmed.
- No unexplained cell loss/gain (105,963 total, matches spec exactly; all 6
  per-study counts match exactly).
- No impossible stage annotations (0 rows in `impossible_state_check.tsv`).
- Integration validation shows no major biological failure (reduced
  same-study neighbor fraction post-integration in every comparable
  population; marker genes preserved; every minority study retains a
  distinguishable signature).
- Required exports are internally consistent (identical, ordered cell keys
  across `normalize_df`/`raw_df`/`scale_df`; all 6 studies in
  `top_genes_exp`).

**This candidate is marked: PHASE-2 RELEASE CANDIDATE.**

It has **not** been deployed. No file under `PlaViSca/` (the application
repository) was created, modified, or touched at any point in this task, except
being read-only for `PlaViSca/ref/PvivaxP01_gff_data.rds` (feature-name join) and
`PlaViSca/data/cleaned_dataset.rds` (read-only, for the Section 7/10
comparisons).
