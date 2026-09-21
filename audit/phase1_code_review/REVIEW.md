# Independent review of PlaViSca Production Repair Phase 1

Review date: 2026-09-20. Reviewed baseline:
`1f5a42a44a83790e33ff40c15b3f6dcf6f1c91e5` (`34ac20c..1f5a42a`).

Scope was read-only inspection of Phase-1 production changes. No production
script, production RDS, count matrix, application file, or D036 evidence was
modified. D036 was not reinvestigated; its unresolved status is accepted from
the shared audit while a separate agent handles it. Review outputs are confined
to this directory.

## Verdict

Phase 1 is **not ready for a Phase-2 atlas rebuild**. Of 62 claimed defect
fixes, 32 are independently verified, 14 are partial, and 16 are incorrect.
The implementation contains multiple blocking execution defects in the shared
annotation and Harmony preparation paths, plus scientifically incorrect Mancio
stage/geography changes and an ineffective Hazzard2022 sporozoite repair.

The verified work is substantial: the Hazzard2024 24-run core mapping is exact,
Sa2020's 9,215-row main-analysis relation is correct, key metadata corrections
are present for Sa/Ruberto/Hazzard2024, the export row-alignment and all-study
top-gene fixes are structurally sound, and candidate export paths do not point
at `PlaViSca/data/`. Those successes do not offset the blockers below.

Full per-defect results are in `defect_fix_verification.tsv`; schema and
reference-table results are in `schema_implementation_check.tsv` and
`ref_table_validation.tsv`.

## BLOCKING before Phase 2

### B1. SingleR eligibility cannot operate on five studies

- **Location:** `scripts/singleR.R:46-53,81-93,153-160`; only
  `scripts/ruberto2022_1_pv_analysis_script.R:368-370` creates
  `total_umi_count`.
- **Problem:** `singleR.R` requires `total_umi_count`, but Sa2020,
  Hazzard2022, Hazzard2024, Mancio-Silva2022, and Ruberto2022_2 never populate
  it. After Seurat merge, the column can exist because Ruberto2022_1 supplied
  it, while all other studies receive missing values. Blood-cell eligibility
  therefore becomes `NA`, not `TRUE/FALSE`; `sum(idc_eligible)` is missing and
  `colnames(...)[idc_eligible]` contains missing cell names.
- **Consequence:** the atlas annotation path fails or passes invalid cells to
  `subset()` before SingleR runs. D021, D034 and D058 are not fixed.
- **Minimal correction:** create an integer/numeric `total_umi_count` for every
  study from the retained RNA counts at preprocessing time, assert it is
  nonmissing/nonnegative for every cell, and assert `idc_eligible` contains no
  missing values before subsetting.

### B2. Hazzard2022 source sporozoites are never supplied to precedence logic

- **Location:** `scripts/hazzard2022_pv_analysis_script.R:225-263` versus
  `scripts/singleR.R:263-276`.
- **Problem:** the preprocessing script sets tissue/sample strings but never
  sets `source_life_cycle_stage="Sporozoite"` or
  `source_stage_provenance="source_selection_defined"` for runs 498/499.
  `singleR.R` only performs its sporozoite override when both fields already
  exist.
- **Consequence:** all 2,858 source-defined Hazzard2022 sporozoites miss the
  repaired override; the confirmed broad/detailed contradiction is not
  structurally repaired.
- **Minimal correction:** populate and assert those immutable source fields in
  Hazzard2022 preprocessing using the two audited run IDs, before merge.

### B3. Marker-driven gametocyte code has two immediate execution defects

- **Location:** `scripts/singleR.R:333-350,525`.
- **Problem:** `gametocyte_eligible` is an unnamed logical vector, but
  `filter(gametocyte_eligible[cell_id])` indexes it with character cell IDs;
  this yields missing values and an empty summary. Separately, `write_tsv()` is
  called twice but is defined nowhere in `singleR.R` or `pipeline_lib.R`.
- **Consequence:** marker classification cannot complete. The first writer
  call fails even if the indexing bug is corrected. D048/D049 are not fixed.
- **Minimal correction:** use an explicitly named eligibility vector or join a
  keyed eligibility column into metadata, assert nonempty eligible rows, and
  define one tested audit/build-output TSV writer with explicit output paths.

### B4. Gametocyte calls do not affect most harmonized stages

- **Location:** `scripts/singleR.R:226-277,399-408`.
- **Problem:** `harmonized_life_cycle_stage` is first filled from IDC bins for
  eligible blood cells. A gametocyte call later overwrites only cells whose
  harmonized stage is `NA`. Thus a marker-called gametocyte with any IDC bin
  remains Ring/Trophozoite/Schizont/Merozoite in the harmonized field.
- **Consequence:** inferred sex and harmonized lifecycle disagree by design;
  the stated precedence/schema is not implemented and removal logic will not
  see most inferred gametocytes.
- **Minimal correction:** encode an explicit precedence function by population
  and provenance, with source locks first and a reviewed rule for marker sex
  versus IDC bins; test every branch on synthetic metadata.

### B5. Eligibility is circular/late for source sexual annotations

- **Location:** `scripts/singleR.R:81-89,422-479`; Hazzard2024 preprocessing.
- **Problem:** Sa2020 source sex is imported only after IDC and marker
  inference. Hazzard2024 author `PseudoGroup` is not imported into the repaired
  study metadata. Eligibility therefore cannot exclude the source-defined
  sexual populations that the audit explicitly required to be excluded from
  blood-IDC inference.
- **Consequence:** even after B1 is fixed, source sexual cells can receive IDC
  labels before their source status is known. This is the circular/late logic
  the repair was meant to eliminate.
- **Minimal correction:** import immutable source sex/PseudoGroup upstream of
  eligibility, then construct eligibility solely from source/tissue/count
  fields; inferred fields must never feed their own eligibility.

### B6. Harmony layer-to-study validation uses the same unsafe indexing pattern

- **Location:** `scripts/integration.R:89-116`.
- **Problem:** `pv.combined.all$study_label[cells_in_layer]` indexes a metadata
  column vector with character cell IDs. Data-frame/Seurat metadata columns are
  not guaranteed to carry cell names. The result can be all missing, causing
  the one-study-per-layer check to fail or report meaningless membership.
- **Consequence:** Phase 2 stops before Harmony; the new assertion does not
  establish its claimed invariant.
- **Minimal correction:** index the metadata data frame explicitly by row name,
  e.g. `pv.combined.all@meta.data[cells_in_layer, "study_label"]`, after
  asserting every layer cell exists exactly once in metadata. Then assert one
  layer per study and one study per layer.

### B7. Mancio stage handling introduces unsupported biology

- **Location:** `scripts/silva.R:31-37,133-155`.
- **Problem:** every source `State=="Replicative"` cell is converted to
  `Schizont (liver stage)`. The completed audit explicitly says Replicative is
  a broad source state and does not justify blanket mature-schizont conversion;
  historical detailed labels came from the separate Table S2 `topcell`
  reference-similarity field, which Phase 1 does not preserve. The 14 `.1`
  records also fail the exact `Updated_names` match and lose corresponding
  source metadata, despite their audited normalized-key multiplicity mapping.
- **Consequence:** 1,221 cells can be assigned unsupported liver schizont labels
  and source/inferred provenance is corrupted. D008/D009 are not fixed.
- **Minimal correction:** preserve source `State`, AP2G/sex state, and `topcell`
  separately; implement the audited normalized-key/multiplicity handling for
  the 14 records without overwriting unique IDs; do not derive Schizont from
  Replicative alone.

### B8. D036 and pending membership decisions remain external prerequisites

D036 is deliberately outside this review and is being handled independently.
It remains a shared-audit blocker to a scientifically interpretable full atlas
rebuild until that work is resolved. Likewise, approved policies for the
Mancio duplicate pairs and Hazzard2022 QC authority must be recorded before a
production rebuild that changes membership. A disposable technical test rebuild
could use the explicit conservative populations, but it must not be promoted or
called the repaired atlas.

## IMPORTANT (non-blocking for a corrected test rebuild)

### Marker thresholds are a new undocumented heuristic

`singleR.R:350-382` hard-codes female/male mean module-score thresholds of
0.10 and a 0.05 margin. These values do not come from the forensic evidence or
a cited biological validation. The two-margin rules prevent a cluster from
satisfying both calls simultaneously; ties/ambiguous clusters become Asexual,
and the script fails only when neither sex has any called cluster. It does not
require both expected sexes, nor does it distinguish genuine absence from a
bad threshold. The summary is written before the `call` column is added and
does not record thresholds. After execution bugs are fixed, this remains a
Phase-2 validation method, not a verified biological classifier. Thresholds,
marker coverage, null behavior, and perturbation stability require AT25/AT30.

### Mancio geography and run relation are contrary to authority

`silva.R:101` asserts `Thailand_Ubon-Ratchathani` for every cell although the
audit supports only Thailand at country level and leaves per-infection province
unresolved. `silva.R:77` writes semicolon-delimited multi-run sets into per-cell
`run_id`; the schema requires a separate group/run relation and a run accession
must be scalar. The committed 15-row crosswalk itself accurately reproduces the
audited group/GEO/run relation; its use is the problem.

### Study-level run counts remain per-cell and one value is wrong

Mancio and all Hazzard2024 chunk scripts replicate deposited/imported/retained
run counts onto every cell, contrary to the study-level schema. Sa2020 assigns
`source_run_count_retained <- ncol(pv.combined.all)`, i.e. 9,766 cells rather
than 10 retained runs. D013/D032 are therefore not implemented faithfully.

### Metadata schema is only partially adopted

The implementation mixes canonical fields with legacy aliases (`study_label`
instead of `study_id`, `run_id` instead of the run relation,
`sc_technology` instead of `technology`, `sample_type`, `strain`, `host_id`).
`source_stage_annotation` is absent, and source stage/provenance fields are not
populated consistently. See `schema_implementation_check.tsv`.

### Export is key-safe but memory-expensive and incompletely guarded

`flatten_data.R:133-165` correctly sorts and asserts cell-key equality before
each `cbind`, and lines 232-283 build all-study `top_genes_df`. The candidate
directory cannot overwrite `PlaViSca/data/` directly. However, three sparse
matrices are converted to dense matrices (`as.matrix`) and retained
simultaneously. At roughly 105,963 cells × ~6,860 genes, one double matrix is
about 5.8 GB before data-frame/copy overhead; three can exceed 17 GB before
intermediate copies. This can make Phase 2 impractical or fail by memory.
Duplicate cleaned metadata names and metadata/gene name collisions are not
explicitly rejected. Candidate files can also be overwritten without a build
ID, although deployed files remain protected.

Top-gene extraction calculates mean log-normalized expression per study after
normalizing from counts. That is reproducible as a descriptive “highest mean
expression” statistic, but it is not a differential marker analysis and should
not be presented as one.

### Build manifest does not meet AT23

`record_build_manifest()` records an output checksum, current Git commit,
selected package versions and one seed field. It records no input checksums,
script checksum, dirty-tree state, or checksums for all three expression
exports. `integration.R` records seed 123 even though 123 is explicitly passed
only to t-SNE; it does not establish that Harmony/PCA/UMAP/clustering all used
that seed. The primary `openssl` path also calls `file(path, raw=TRUE)`, an
invalid base-R signature, and silently relies on the digest fallback.

### Application defects D063-D068 were not implemented

The Phase-1 diff contains no PlaViSca application/config changes. Dead stage
columns, hard-coded label/color maps, unreachable NA context, asymmetric
gametocyte colors, missing app-startup bibliography join, and biological legend
ordering therefore remain. D067 is only partially addressed by validating the
bibliography build vectors against preprocessing constants; that is not the
required app-startup join to `mr_data`.

### Source-table results

- Hazzard2024: exact 24-run key set; library, animal, day, infection design,
  host, registered/protocol sequencer, host sex, registry location, and host
  provenance all match `audit/hazzard2024/authoritative_run_table.tsv`.
- Mancio: all 15 group/GEO/run/day/treatment/instrument relations match the
  audited table; eight groups correctly retain two associated runs.
- Sa2020: 9,215 unique author IDs; exactly 9,018 are marked present; run totals
  reproduce the audit.

## `pipeline_lib.R` helper review

The direct-assignment loop in `apply_keyed_run_metadata()` preserves each
Seurat object's cell order. It also rejects duplicate reference-table keys,
unmatched keys, extra reference keys, and missing requested fields. Those are
useful safeguards.

Hidden assumptions remain:

- `assert_keyed_join()` does not reject duplicated **data-side** keys despite
  its documentation, and accepts a single missing key when both sides contain
  `NA`. The production test's duplicate case also contains an extra reference
  key, so it passes for the wrong reason.
- `apply_keyed_run_metadata()` does not validate named-list uniqueness, missing
  keys, field types, controlled vocabularies, or typed optional missingness.
- `split_run_membership()` performs no validation before partitioning.
- `map_broad_stage()` silently mis-maps factor input because a factor indexes
  the named vector by integer level codes. The review test maps factor
  `Sporozoite;Ring` to `Liver stage;Liver stage`.
- `assert_rowname_set_equal()` does not itself reject duplicate keys; callers
  must remember a separate uniqueness assertion.
- no helper enforces treatment vocabulary, day range, stage provenance
  vocabulary, or required schema types across all studies.

Deterministic evidence is in `helper_edge_case_results.tsv`.

## Harmony preparation assessment

If B6 is corrected, the intended mechanism is clear: six joined single-study
objects are merged, yielding six count/data layers; `NormalizeData`, HVF,
scaling and PCA run on the multi-layer object; Seurat's layer-derived
`CreateIntegrationGroups` supplies Harmony with one group per verified study;
`JoinLayers` occurs only after `IntegrateLayers`. With six layers, each mapping
to one of six unique expected studies, a study cannot be silently split across
multiple layers without another study being absent, which the set assertion
would catch.

The change mainly makes the previously layer-derived behavior explicit and
testable; it does not demonstrate that the old integration grouped cells
incorrectly or that global study-level Harmony is biologically appropriate.
AT24-AT27 remain necessary and are correctly deferred.

## Study-script preservation assessment

- **Mancio-Silva2022:** count extraction and 1,494-cell retention are preserved;
  14 `.1` records and the D5 pair are flagged, not deleted. Metadata biology
  is not faithful because of B7, unsupported province assignment, incomplete
  source metadata, and scalar misuse of multi-run relations.
- **Sa2020:** the 9,766-cell STARsolo population is retained. The committed
  relation correctly gives 9,018 main-analysis members; it does not switch to
  the publication population. Animal pairing is explicit.
- **Hazzard2022:** the original emptyDrops code and exact 3,294 population are
  unchanged, appropriately leaving QC authority unresolved. Source-stage
  wiring is missing (B2).
- **Hazzard2024:** selection remains keyed to the 80,024 processed IDs; no old
  positional metadata vectors remain. The authoritative 24-run core fields
  match. Study-level run counts are at the wrong granularity.
- **Ruberto2022_1:** 1,438 cells are retained and 538 are flagged, not silently
  corrected or removed. No further D036 investigation was performed here.
- **Ruberto2022_2:** source selection remains 9,947 sporozoites with no intended
  count/membership change; source-stage lock is correctly created.

## Treatment semantics

Mancio source `CTRL/PI4K` is preserved and harmonized to
`study_control/PI4K_inhibitor`, correctly avoiding a false untreated claim.
Sa2020 preserves `None` as source text and uses `No_Treatment` only for the
documented untreated arm; chloroquine arms remain named. Ruberto2022_1
preserves source `None/MMV390048` and harmonizes appropriately. Hazzard2022,
Hazzard2024 and Ruberto2022_2 use `No_Treatment` while retaining infection
design separately where relevant. No production assignment deliberately uses
literal `"NA"` as treatment. A shared controlled-vocabulary validator is still
missing, and the export check ignores factor columns.

## Test-suite review

The checked-in test suite could not run in the baseline Pixi environment:
`env -u R_LIBS_USER pixi run Rscript --vanilla tests/test_shared_pipeline.R`
failed while loading `dplyr`. The baseline manifest lists the dependency, but
the realized environment is incomplete; unrelated concurrent worktree changes
attempting environment setup were not reviewed or committed here.

Independently of that environment failure, the suite has major coverage gaps:

- UNIT-02's duplicate-data-key case fails because of an extra table key, so it
  does not test the behavior it claims.
- No negative tests cover missing keys, factor mapping, factor/literal-NA
  handling, duplicate cleaned names, list-name duplicates, or optional typed NA.
- Study tests check only cell count, unique IDs and literal NA; they do not
  compare run metadata with authoritative tables, expression/count identity,
  source fields, treatment semantics, or approved populations.
- AT13 (the exact Sa2020 count of 1,532 matched source-sex rows) is deferred to
  Phase 2 even though its source-import transformation can be tested now.
- No test statically checks that every study emits fields required by
  `singleR.R`; this would have caught missing `total_umi_count` and Hazzard2022
  source-stage fields.
- No test executes the marker-summary code on a tiny synthetic metadata frame;
  this would catch the unnamed-vector indexing and undefined writer.
- No synthetic multi-layer integration test validates layer-to-study lookup.
- AT21 tests constants copied from the production mapping, not an independent
  authoritative mapping fixture; it also tests character input only.
- No export test covers duplicate gene names, metadata/gene collisions,
  factor preservation, zero-column selection, or memory-safe sparse handling.
- Build-manifest contents are not tested against AT23.

The review adds only small base-R/Python checks; it does not rerun the atlas.

## MINOR

- Several scripts source `pipeline_lib.R` twice and retain legacy alias fields,
  increasing schema ambiguity.
- Root discovery still includes a personal hard-coded fallback. Running from
  the `scripts/` directory without environment variables is not consistently
  supported even though historical audit commands used that directory.
- `record_build_manifest()` appends without locking or a build identifier.
- Hazzard2024 chunk logic is duplicated six times, increasing drift risk even
  though the current core table is correct.

## VERIFIED components

- 32 claimed defect fixes are independently verified; see the TSV for each.
- Hazzard2024 authoritative core run table: zero mismatches in directly
  comparable 24-run fields.
- Mancio and Sa relation tables reproduce their audited source mappings/counts.
- Study cell-count assertions preserve 1,494 / 9,766 / 3,294 / 80,024 / 1,438 /
  9,947 populations in the intended scripts.
- Correct `%in%` and NIH/PB_MACS mapping are present for Sa source sex import.
- `graph.name="RNA_snn"` and graph existence assertion replace the invalid API.
- Export row-key assertions and all-study top-gene accumulation are present.
- Candidate expression/cleaned exports cannot directly overwrite deployed
  `PlaViSca/data/*.rds` paths.

## Phase-2 readiness assessment

Not ready for an atlas rebuild. Correct B1-B7 and add focused synthetic tests
before even a test integration/annotation run. Resolve the externally owned
D036 blocker and record the outstanding membership policies before producing
or promoting a scientific candidate. Once those are addressed, a test rebuild
can proceed, but publication readiness additionally requires the IMPORTANT app,
schema, manifest, marker-validation, integration-validation, and export-memory
items above.

## Evidence manifest

| File | Purpose |
|---|---|
| `defect_fix_verification.tsv` | All 62 claimed fixes, location, summary, regression risk, verdict |
| `schema_implementation_check.tsv` | Field-by-field implementation review against the authoritative schema |
| `ref_table_validation.tsv` | Deterministic comparison of the three committed keyed tables to audit evidence |
| `helper_edge_cases.R` | Small base-R helper/character-index edge tests |
| `helper_edge_case_results.tsv` | Results of those edge tests |
| `review_checks.py` | Reproducible generator for the three review TSVs |

No large evidence was duplicated.
