# Mancio-Silva2022: metadata and accession audit

Audit date: 2026-09-19. Status: inspection complete; repairs proposed, not applied.

Preprocessing branch: `plavisca-audit-fixes`. Production-code baseline:
`e5b3de5cfaa93d4e1a6ed8c8c4722dbcf1ec9144`; this is not the current audit-branch HEAD.
The branch contains audit-only commits `b1b9962` and `4e34a94`, recording the report
and supporting evidence. No production pipeline scripts, source data or application
objects were changed. The present documentation revision is not committed.
Changes are confined to this audit directory and its evidence files. R inspections used
`env -u R_LIBS_USER pixi run Rscript --vanilla` from `pre_process_data/scripts`.
Runtime: R 4.3.3, Seurat 5.3.0, SeuratObject 5.2.0, readxl 1.4.5.

## Evidence and identity

The authors' [Zenodo version 2.0 deposit](https://zenodo.org/records/6280956)
links to GSE197409 and describes both source files. Their advertised MD5 checksums
exactly match the local files:

| File under `../../data` | MD5 |
|---|---|
| PvData_CHM_Final.RDS | b58ba79595bee10121dc374aa160bc30 |
| TableS2_Cell_Metadata.xls | 5450e12b6f2a51c0e5197c4f30fd91c8 |

This establishes identity with the deposited files, not the absence of scientific
errors in those files. The object has 1,494 cell columns, 78 source metadata fields,
and an RNA count matrix with 5,934 features. Counts inspected are nonnegative integers.
Source object version is 4.0.3. RNA counts, rather than SCT/integrated values, are
the appropriate existing input to preserve for the proposed metadata repair.
Gene-name truncation to the first two hyphen-delimited components creates no duplicate
identifiers in this object; annotation/reference equivalence is not established by that check.

Authoritative sources downloaded into this directory:

- `zenodo_6280956.json`: https://zenodo.org/api/records/6280956
- `GSE197409.txt`: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE197409&targ=self&form=text&view=full
- `GSE197409_samples.txt`: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE197409&targ=gsm&form=text&view=full
- `ena_SRP361380_runs.tsv`: ENA file-report API, accession SRP361380, result read_run;
  columns include run, experiment, sample, study, sample title, instrument, and library descriptors.
- `ena_SRP361380_files.tsv`: ENA file-report API for the same study, including FASTQ URLs.
- `SRX14282321.xml`: https://www.ebi.ac.uk/ena/browser/api/xml/SRX14282321
- `SRR18134252.xml`: https://www.ebi.ac.uk/ena/browser/api/xml/SRR18134252
- `article.xml`: https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12285726/fullTextXML
- Publication identity independently checked at https://pubmed.ncbi.nlm.nih.gov/35443155/.

Derived evidence: `geo_sample_summary.tsv`, `cell_metadata_evidence.tsv`,
`source_group_accession_map.tsv`, `represented_group_runs.tsv`,
`duplicate_count_pairs.tsv`, and `validation.txt`. These are audit outputs, not
replacement metadata or production pipeline inputs.

## Cell identifiers and source fields

The full source Seurat cell names are unique. Table S2 has 1,494 rows but only
1,480 distinct identifiers in BOTH `Cell` and `Updated_names`.
There are 14 duplicated identifier pairs. In the Seurat object, one member of
each pair has a terminal `.1`; 1,480 names match `Updated_names` exactly and the
remaining 14 match after removing precisely that suffix for comparison.
The normalized-key multiplicities match between the two files.
This normalization must not replace the unique original cell IDs.

`Cell` is recovered from the source `orig.ident`, underscore, and the barcode
part of `Updated_names`. This matches all 1,494 records. Prefix renaming, for
example `CAPD5r50` to `D5SEQ50`, is consistent within each group.

For duplicate keys, `orig.ident`, `Sample`, `Day`, `Treatment`, `State`, and
`sex.state` agree within each pair. Across all 1,494 cells, the Seurat and Table S2
values of the first five fields agree exactly. The sexual-state labels correspond
exactly as follows:

| Seurat source | Table S2 source | Cells |
|---|---|---:|
| Non.Committed_Hypnozoite | AP2Gneg_Hypnozoite | 231 |
| Non.Committed_Replicative | AP2Gneg_Replicative | 1,114 |
| Sexually.Committed_Hypnozoite | AP2Gpos_Hypnozoite | 42 |
| Sexually.Committed_Replicative | AP2Gpos_Replicative | 107 |

These are corresponding source annotations, not independently demonstrated mature
male/female gametocytes. Preserve both vocabularies and their provenance.

`orig.ident` identifies 15 complete source groups. `Sample` identifies two
infection experiments: Infection1 (1,270 cells) and Infection2 (224 cells).
Neither is a patient/donor identifier. The author deposit explains `SEQ50` and
`SEQ100` as inoculum groups, and `SEQ1`/`SEQ2` as replicate arrays in Infection2.
GEO specifies 50K/100K sporozoites for Infection1; 50/100 are not read lengths.

Source-field defects:

- `orig.ident1`: 21 actual missing values and four nonmissing `0.6` values that
  disagree with `orig.ident`. Other nonmissing entries agree. Do not use it as a join key.
- `sample.day`: 21 actual missing values; all other entries agree with
  `paste(Sample, Day, sep = "_")`. Preserve source values; derive a separate complete field.
- `treatment.day.sample`: complete and agrees with
  `paste(Treatment, Day, Sample, sep = "_")` for all cells.

## Accession membership and group mapping

GEO explicitly links GSE197409 to PMID 35443155 and PRJNA810174. ENA returns exactly
the 58 unique runs SRR18134227–SRR18134284 for SRP361380 / PRJNA810174.
Thus the accession range belongs to the correct single-cell study.

There are 44 deposited GEO samples / SRA experiments: 16 Infection1 infected
samples with one run each, and 28 Infection2 samples. The latter include 14
infected samples with two runs each and 14 mock/uninfected samples with one run each.
The 58 are study-level deposited runs, not 58 parasite-cell libraries retained in
the final object. The retained source groups map to 15 GEO samples and 23 associated runs.

The mapping below is reconstructed from the authors' documented naming scheme,
the source cell-level infection/day/treatment fields, and GEO sample titles and
characteristics. Every group has a unique matching GEO sample; GEO day and treatment
agree for every cell. Run membership follows the GEO sample accession in ENA.
It establishes sample-associated run sets, not direct observation of every cell's
reads in every run.

| orig.ident | Cells | GEO sample | Associated SRR run(s) |
|---|---:|---|---|
| CAPD1r100 | 207 | GSM5916705 | SRR18134228 |
| CAPD5r100 | 127 | GSM5916707 | SRR18134230 |
| CAPD5r50 | 177 | GSM5916706 | SRR18134229 |
| CAPD11r100 | 313 | GSM5916711 | SRR18134234 |
| CAPD11r50 | 213 | GSM5916710 | SRR18134233 |
| CAPD14pr100 | 71 | GSM5916719 | SRR18134242 |
| CAPD14r100 | 162 | GSM5916713 | SRR18134236 |
| SEQ4r1 | 63 | GSM5916726 | SRR18134251;SRR18134252 |
| SEQ4r2 | 17 | GSM5916727 | SRR18134253;SRR18134254 |
| SEQ5r1 | 1 | GSM5916730 | SRR18134257;SRR18134258 |
| SEQ5r2 | 1 | GSM5916731 | SRR18134259;SRR18134260 |
| SEQ8pr1 | 28 | GSM5916742 | SRR18134275;SRR18134276 |
| SEQ8r1 | 37 | GSM5916734 | SRR18134263;SRR18134264 |
| SEQ11p2 | 54 | GSM5916747 | SRR18134283;SRR18134284 |
| SEQ11r1 | 23 | GSM5916738 | SRR18134269;SRR18134270 |

For full sample titles, day, treatment, experiment accessions, and registered
instruments, see `source_group_accession_map.tsv`.

The sample/library structure is cells within Seq-Well arrays, nested in infection,
collection day, treatment, and inoculum/array group. GEO describes pooled culture
wells before array loading. Infection2 processed parasite counts combine matching
barcodes from pre- and post-capture data. Consequently one scalar SRR per cell
cannot adequately express its provenance. Do not reinterpret multiple sequencing
runs as biological replicates.

## Platform, treatment and geographic evidence

GEO/ENA register Infection1 as NextSeq 500 and Infection2 as Illumina NovaSeq 6000.
The currently retained groups therefore have registered labels for 1,270 and 224
cells respectively, rather than the recycled 416/1,078 split in the app.
However, the article describes NovaSeq for uncaptured material and NextSeq 500
for captured material from both infections. The Infection2 ENA experiment XML
itself includes a NextSeq capture protocol under a NovaSeq platform label.
This is an unresolved metadata discrepancy. Neither accession order nor read
counts establish which run is pre/post capture. Retain registered instrument and
publication protocol separately; do not claim verified run-specific instruments
or assign capture roles to SRRs without further evidence.

The article identifies the treatment as MMV390048 (PI4K inhibitor), 1 micromolar,
days 5–8. GEO identifies the relevant samples as Pi4K-treated or non-treated.
Thus source `Treatment` supports CTRL -> study control and PI4K ->
PI4K inhibitor treatment. Preserve `source_treatment = CTRL/PI4K`. Control status is relative to this experimental treatment;
it does not mean that the cultures underwent no other preparation or manipulation.

The article describes patient parasite collection in Tak, Songkla and
Ubon-Ratchathani, Thailand, and US-procured primary human hepatocytes. It does not
establish a per-infection province or support assigning Rockville to the parasites.
Use separate parasite-origin and host-donor-origin concepts. Thailand is supported
at country level for parasite provenance; a province for each infection remains
unresolved. `host_species = Homo sapiens` describes the hepatocyte host, not the
species of the parasite RNA features. Seq-Well is supported, with S3 protocol detail.

## Confirmed defects in PlaViSca construction/export

Checked against `PlaViSca/data/cleaned_dataset.rds` (`mr_data`, 1,494 Mancio rows),
using exact source cell names in row names. The final `barcode` column is missing
for all 1,494, despite its assignment in the current script. This is additional
evidence that the existing app artifact cannot simply be assumed to reproduce
the current checked-out scripts.

| Problem | Evidence | Proposed correction | Responsible files |
|---|---|---|---|
| Wrong PMID | 39223117 on all cells; GEO/PubMed support 35443155 | Set 35443155 | silva.R; later README/data-source consistency |
| Recycled run IDs | All 1,494 exactly match repeating the 58-run vector in source-cell order; reproduced in memory in Seurat | Join source group to GEO, keep run sets in a relation table | silva.R; export/app schema if needed |
| Wrong days | 1,222/1,494 cell-level disagreements | Parse validated source Day | silva.R |
| Lost PI4K treatment | All final rows No_Treatment; 153 source PI4K cells | Derive from source Treatment | silva.R |
| Platform vector and typo | Canonical sequencer missing on all; typo field has 416 NextSeq / 1,078 NovaSeq | Canonical field plus explicit registry/protocol provenance and unresolved status | silva.R |
| Geographic typo and unsupported specificity | Canonical field missing; typo value is USA_RockvilleThailand_Ubon-Ratchathani | Separate origin concepts; only supported country-level parasite location | silva.R; investigate historical transformation |
| Inconsistent study label | Script Mancio Silva2022; app Mancio_silva2022 | Mancio-Silva2022 | silva.R and downstream study selectors |
| Discarded source metadata | CreateSeuratObject(counts=raw) without source metadata | Preserve keyed original metadata and source fields | silva.R |
| Fragile state join | grep prefix matching rather than cell-keyed metadata | Copy validated source State by exact cell ID | silva.R |
| Missing final barcode | All 1,494 NA | Preserve source_cell_id; validate exact export correspondence | silva.R and export provenance |

Source day counts: d1=207, d4=80, d5=306, d8=65, d11=603, d14=233.
Source treatment counts: CTRL=1,341; PI4K=153.

Stage issue found during comparison: source State is Hypnozoite (273) or Replicative
(1,221). Final app stages split these into multiple inferred labels; one source
hypnozoite is called Female gametocyte. The other 272 remain Hypnozoite. Replicative
cells include Ring, Trophozoite, Merozoite, EEF, male/female gametocytes and liver
schizonts. This does not validate those inferences. Source State, source reference
similarity labels and PlaViSca predictions must remain separate. A blanket conversion
of Replicative to mature Schizont is not justified solely by source State.

## Repeated count profiles: conclusions from completed upstream tracing

All 14 same-group `.1` pairs are **proven duplicate representations of single
upstream GEO count columns**. Both source Seurat records in each pair match the
same upstream column across all 5,934 RNA genes. The mechanism that created the
second representation between the deposited upstream matrices and final authors'
object remains unknown; intent and the exact processing operation are not proven.
PlaViSca preserved these duplicate representations unchanged.

The cross-array pair `D5Seq1_CCCCGATTGACG` / `D5Seq2_CCCCGATTGACG` remains
**unresolved**. Two separately named upstream GEO columns exist, but their full
85,578-feature host-plus-parasite count profiles are identical. Biological
independence is not demonstrated. No survivor is selected for this pair.

The detailed upstream evidence is recorded in the focused duplicate-tracing
section and `exact_duplicate_trace.tsv`. All 1,494 deposited records remain
preserved for provenance and source-to-final comparison. The proposed future
analytic-inclusion policy below is separate from this unchanged source inventory.

## Exact proposed changes to silva.R (not applied)

1. Replace the hard-coded Sopheap working directory with explicit project-relative
   input/output resolution. Verify both input checksums and source cell uniqueness.
2. Retain extraction of the RNA counts. Retain the gene-ID transformation only with
   an explicit source-feature-to-output-feature mapping and uniqueness assertion;
   do not change values, cell order, or cell inclusion in this metadata repair.
3. Preserve all original metadata in a keyed source metadata artifact, retaining
   the original field names and values. In the new object expose `source_cell_id`,
   `source_orig_ident`, `source_orig_ident1`, `source_sample`, `source_day`,
   `source_treatment`, `source_sample_day`, `source_treatment_day_sample`,
   `source_state`, and `source_sex_state`. Label imported source computational
   annotations as author-provided rather than new PlaViSca computations.
4. Set `study_pmid = "35443155"`, `study_label = "Mancio-Silva2022"`,
   `pub_year = 2022`, and preserve DOI/GEO/BioProject/SRP provenance explicitly.
5. Delete all 58-element metadata vectors. After exact cell-ID alignment, derive
   `day_post_infection = as.integer(sub("^d", "", as.character(source_day)))`.
   Require source_day in d1/d4/d5/d8/d11/d14, with no missing or unexpected values.
6. Preserve `source_treatment = CTRL/PI4K` exactly. Derive a scoped
   `plavisca_treatment_arm`: CTRL -> `study_control`, PI4K -> `PI4K_inhibitor`.
   Display these as "Study control (CTRL)" and "PI4K inhibitor (MMV390048)".
   Never translate CTRL to No_Treatment. Keep compound, concentration, exposure
   timing, and any verified vehicle information separate. Fail on unknown source labels.
7. Set `infection_id` from source Sample; derive a separate complete
   `sample_day = paste(source_sample, source_day, sep = "_")` and
   `treatment_day_sample` similarly. Keep source composites unchanged, including
   their missing values. Do not use orig.ident1 to populate library identifiers.
8. Add a reviewed source-group-to-GEO crosswalk backed by the deposited naming
   notes and GEO characteristics. Join on source_orig_ident, validate source
   infection/day/treatment agreement, and attach geo_sample_id and experiment_id.
   Maintain one row per source group/run in a separate accession relation table.
   Use `run_ids` as the complete sorted associated-run set if cell-level export
   needs it. Do not select an arbitrary one of two runs for Infection2. If legacy
   scalar `run_id` must temporarily remain, populate it only for single-run groups
   and use NA plus explicit multi-run status for other groups; update consumers
   before regenerating production exports.
9. Remove `num_srr` from ordinary cell metadata. Store deposited study runs=58
   and retained-group associated runs=23 in study provenance, with accession
   snapshot date and scope. Store associated runs per source group=1 or 2 in
   the group/run relation. None of these is a biological replicate count.
10. Remove typo assignments. Use explicit `parasite_origin_country = Thailand`
    and a separate `host_tissue_provenance` record for US-procured primary human
    hepatocytes, with evidence and resolution. Province per infection remains
    unresolved. Do not conflate procurement country with donor residence or use
    a combined parasite/host geographic_location field. Keep host_id unknown;
    do not fill it with Infection1/2. Preserve patient-derived isolate status.
11. Keep a harmonized Seq-Well technology label with source protocol Seq-Well S3.
    Store GEO/ENA instrument as `sequencer_registered`; Infection1 sequencer can be
    NextSeq 500, whereas Infection2 physical sequencing provenance remains
    unresolved/multiple-protocol pending clarification. Record the discrepancy,
    rather than filling canonical sequencer with an unqualified NovaSeq claim.
12. Replace the grep loop with exact keyed source-State import. Preserve State and
    AP2G/sexual-state labels. Prefer a broad harmonized Replicative (liver stage)
    category to unsupported automatic conversion of every replicative cell to
    Schizont; review its downstream vocabulary explicitly before implementation.
    Keep stage predictions in separate fields. Set broad phase to Liver stage.
13. Preserve `source_cell_id` exactly, including `.1`. Derive a separate
    `source_barcode_token` from the documented identifier structure, recording
    removal of the uniqueness suffix only in that derived token. It is a processed
    barcode token; it is not verified as the uncorrected sequenced barcode. A truly
    raw barcode requires read-level evidence. Do not call source_cell_id a raw
    barcode, and do not overwrite production barcode fields before reviewing app
    compatibility. Assert metadata/count
    cell ordering, no unmatched joins, expected distributions, and zero controls
    in the parasite-group accession mapping before writing silva2022.rds.
14. Emit provenance (source checksums, metadata authority, crosswalk evidence and
    unresolved flags). Fail rather than silently recycling or introducing unknown
    source categories. Preserve duplicate evidence and flags; apply any future analytic
    inclusion change only through the separately reviewed policy below.

## Downstream regeneration and remaining decisions

After review and upstream repair, regenerate `silva2022.rds`, then propagate metadata
through the workflow producing `pv_all_studies.rds` and `flatten_data.R` outputs:
`normalize_df.rds`, `raw_df.rds`, `scale_df.rds`, and `cleaned_dataset.rds`, including
study-keyed summaries. `singleR.R` contains a legacy Mancio study selector and
downstream stage overrides that need a separate reviewed repair. Audit exact stage
provenance before changing biological labels. Do not patch final RDS files.

The central integrated object is absent. Reconstructing it is a later task, after
dataset and cell-inclusion decisions. Metadata-only fixes do not themselves prove
that expression or embeddings require recomputation; cell removal or count changes
would. Current reports of metadata correction must not be confused with integration
validation. Nothing has been regenerated in this audit.

Open questions: the unknown mechanism creating the 14 proven `.1` duplicate
representations; the origin and biological independence of the unresolved day-5
cross-array pair; per-run capture role and true physical instrument for Infection2;
per-infection parasite province/donor provenance; precise cell/read contribution
from each associated run; and the validity and applicability of inferred lifecycle
labels. Source-to-final cell identity, raw counts and comparable normalized
expression are resolved by the completed lineage checks. The historical script
state reproducing the observed metadata is identified; the exact executed build
and date remain unverified without a build manifest or the central object. These
remaining uncertainties do not prevent source-keyed day/treatment correction,
but do prevent declaring the whole study publication-ready.

### Proposed analytic-inclusion policy (not applied)

- Preserve all 1,494 deposited records and their original identifiers in the source
  inventory and audit crosswalk. Do not alter the deposited object or application
  artifacts in this documentation update.
- In a future analytic object, remove the 14 terminal `.1` copies identified in
  `exact_duplicate_trace.tsv`, retaining their corresponding unsuffixed records.
  This rule applies to these 14 proven duplicate representations, not to suffixes
  in general. Record each exclusion, its retained counterpart and the shared
  upstream GEO column in an explicit inclusion/exclusion manifest.
- Flag both `D5Seq1_CCCCGATTGACG` and `D5Seq2_CCCCGATTGACG` for a separate reviewed
  decision. Do not select a survivor or claim that they are independent biological
  cells. Removing only the 14 proven copies would leave 1,480 records, including
  this unresolved pair; 1,480 is not a validated independent-cell count.
- If this policy is approved and implemented later, regenerate affected analytic
  objects and downstream analyses from the responsible upstream step, validate
  the exclusion manifest and resulting cell inventory, and document the separate
  decision for the day-5 pair. No exclusion or regeneration has been applied.

Validation after correction: not applicable because no correction was applied.
Inspection validation is recorded in `validation.txt`. One initial inspection
assertion compared entire R table objects and failed on differing dimension labels;
comparison of explicit identifier names/multiplicities then passed. An exploratory
gene replacement used an incorrectly escaped replacement string; it was discarded
and the split-based uniqueness check passed. Neither diagnostic changed any object
on disk. Final tracked-file diff was empty before writing this audit document.

## Focused follow-up: source-to-final lineage (2026-09-19)

This section supersedes the earlier uncertainty about the origin of the current
Mancio metadata. It identifies an exact historical computational signature, while
not claiming a recovered execution log or definitive build timestamp.

### Cell and expression identity

`inspect_lineage.R` performs identifier-based comparisons and writes
`source_to_final_lineage.tsv` (one row for every source cell),
`expression_concordance.tsv`, `artifact_metadata_concordance.tsv`, and
`artifact_column_schemas.tsv`. All are audit evidence; no production object was saved.

- All 1,494 full source cell IDs occur exactly once in each of raw_df,
  normalize_df, scale_df, and cleaned_dataset$mr_data. Each artifact has unique
  row names and the same complete 105,963-row cell set.
- All 5,934 source RNA features map uniquely to app feature IDs by the existing
  gene-name truncation/hyphen-to-underscore rule. All 8,865,396 Mancio raw values
  agree exactly. The additional 926 parasite-gene columns in the app raw and
  normalized files are zero for these cells.
- App normalized expression exactly equals `log1p(count / cell_total * 10000)`
  for all 8,865,396 source-gene values. It also agrees with the authors' RNA data
  layer to floating-point precision (maximum difference 8.881784e-16).
  This is a numerical concordance check, not a normalization-method endorsement.
- scale_df contains 2,000 source genes; cell identity and shared metadata agree.
  Its integrated-workflow scaling was not equated to the authors' SCT residuals.
  Reconstructing scale parameters is outside these two focused investigations.
- The expression files are sorted by cell ID. cleaned_dataset$mr_data retains
  the source Mancio order. `source_to_final_lineage.tsv` records all positions;
  positional cbind between these orders would be incorrect.
- All 37 shared metadata/embedding columns agree across the four artifacts for
  all Mancio cells. Raw/normalized/scaled files omit 15 fields present in cleaned:
  n_count_rna, n_feature_rna, ten rna_snn_res columns, and three t_sne columns.

There is therefore a one-to-one source-to-final *record* lineage with unchanged
RNA counts. This does not establish 1,494 independent biological cells: duplicate
profiles are inherited from the deposited source object, and expression alone
cannot distinguish the two members of an exact-duplicate pair.

### Exact historical metadata construction recovered

Although the standalone preprocessing repository has only its initial commit,
the local PlaViSca repository retains earlier preprocessing scripts in Git history.
No branch was checked out. Archived script snapshots are saved under this audit
directory for evidence and must not be run as production scripts.

Commit `c80580f527cbea9641ddc2d6592c85609095b43c` (2025-09-18),
`pre_process_data/scripts/silva.R`, reproduces **all 19 comparable metadata fields
for all 1,494 cells**. The diagnostic evaluated the study construction in memory,
skipped setwd/saveRDS, and read only existing source files. See
`historical_script_metadata_comparison.tsv`.

Its distinctive assignments include Mancio_silva2022,
USA_RockvilleThailand_Ubon-Ratchathani, sc_technology=NA,
Infected host hepatocytes, parasite_stage=Liver stage, lower-case patient isolate,
and no barcode field. The same invalid positional run/day/platform vectors are
already present there. These are historical construction errors, not evidence
of corruption of expression during the final flattening step.

The old refine_state rule retained Hypnozoite from source State, otherwise used
Table S2 `topcell` (a source reference-similarity annotation), translated
Male/Female labels, and joined only exact Updated_names. This left precisely the
14 `.1` source records without refine_state.

Combining that historical rule with the old singleR.R override order and the
recorded final HPI/gametocyte predictions reproduces **all 1,494 final lifecycle
labels**. This reconstructs label propagation; it does not rerun or independently
validate SingleR or clustering. The resulting provenance breakdown is:

| Effective label origin | Cells |
|---|---:|
| Source State: Hypnozoite | 272 |
| Source Table S2 topcell reference similarity | 1,208 |
| PlaViSca IDC/cluster fallback for unmatched `.1` records | 14 |

The last group comprises four Female gametocyte, seven Schizont (Liver stage),
two Trophozoite and one Merozoite labels. One source Hypnozoite `.1` copy acquires
a Female gametocyte label through that fallback. See `historical_stage_comparison.tsv`.

The historically identifiable divergence points are:

| Commit/date | Checked-in change absent from the observed app signature |
|---|---|
| 03020ee, 2025-10-06 | Changed study label, sample terminology, strain case and State-based refine_state logic |
| 4875896, 2025-10-07 | Changed location/technology, expanded state matching with grep, added non-blood HPI exclusion and moved stage overrides |
| 5e8e144, 2025-10-09 | Added barcode assignment |
| Later 2026 terminology/export changes | Changed stage field vocabulary and output schema; current flatten_data.R removes fields retained by the app |

The historical export at c80580f also explains the 15 omitted metadata/embedding
columns: it selects original metadata names after clean_names has changed them,
and searches for `tsne` although the cleaned names contain `t_sne`. Its construction
of `Schizont (Liver stage)` from old life_cycle_stage=Schizont and broad
parasite_stage=Liver stage matches the app capitalization. The current app-side
export changed that lookup from parasite_stage to parasite_stages in commit
9739251 (2026-02-22); the standalone current flattening script instead comments
out that mutation and removes other observed columns.

The app artifacts are thus numerically consistent with the old study/export
logic and inconsistent with a fresh execution of the current scripts. It is
not possible to prove which exact historical commit was executed, whether mixed
versions were used, or when the artifacts were built: no run manifest or central
Seurat object is available. File modification times are not build provenance.

### Revised metadata design, authoritative for the next repair proposal

- Preserve source_treatment as CTRL/PI4K. Use a study-scoped harmonized arm,
  `study_control` versus `PI4K_inhibitor`, with explicit labels. CTRL is never
  recoded to No_Treatment; vehicle details remain unspecified unless verified.
- Separate parasite geographic origin from host-tissue procurement/provenance.
  Do not use a concatenated geographic_location as a source of either concept.
- Keep registered instrument, publication protocol evidence and conflict status
  separately. Do not silently resolve the Infection2 platform discrepancy.
- Keep source_cell_id, source group and processed barcode token separate.
  A uniqueness suffix is not part of the molecular barcode, and the processed
  token is not automatically the raw sequenced barcode.
- Put study-level run counts in a study provenance table, with snapshot and scope;
  maintain group-to-run relations separately from cells.
- Preserve source State, source sex.state/AP2G, source topcell reference similarity,
  and PlaViSca-inferred lifecycle labels in distinct fields with explicit authority.
  No inferred label may silently overwrite source State.

These are conceptual proposals only. Production code and data remain unchanged.

## Focused follow-up: exhaustive duplicate tracing (2026-09-19)

### Deliverables and scope

`exact_duplicate_trace.tsv` contains one row for every affected source record
(30 rows, 15 pairs), including full source cell ID, original Table S2 Cell and
Updated_names, orig.ident and orig.ident1, source group, Sample, Day, Treatment,
State, sexual-state annotation, GEO sample, associated runs, source/upstream column
positions, comparison metrics, paired ID and upstream identifiability conclusion.
`DUPLICATE_TABLE.md` is a readable rendering. No cells or source files were deleted.

An exhaustive SHA-256 grouping of all 1,494 full RNA count columns, followed by
elementwise verification of each candidate group, found exactly the previously
identified 15 pairs and no additional exact RNA duplicates. Each pair has two
members. The original object therefore contains 1,479 distinct RNA count profiles;
this is not, by itself, an estimate of the number of independent biological cells.

The pairs are identical not only in RNA counts but in RNA data, SCT counts/data/
scale.data, and integrated data/scale.data, wherever those layers have features.
All 30 records are present in the appropriate source SCT model's cell attributes.
See `duplicate_assay_comparison.tsv` and `duplicate_sct_model_membership.tsv`.
This proves the duplication is embedded in the deposited processed object; it
does not prove which exact author-side operation first introduced it.

### Upstream GEO matrices actually inspected

All three archives were downloaded from the official
[GEO supplementary directory](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE197nnn/GSE197409/suppl/),
retained under this audit directory, and extracted only here:

| Archive | Matrix dimensions (features x columns) | Unique column IDs | Inspection |
|---|---:|---:|---|
| GSE197409_Postcapture_Infection1.tar.gz | 48,509 x 159,984 | 159,984 | Complete coordinate stream scanned; affected columns and same-token alternatives compared |
| GSE197409_Postcapture_merged_Infection2_corrected.tar.gz | 85,578 x 209,899 | 209,899 | Complete coordinate stream scanned; three affected upstream columns compared |
| GSE197409_Precapture_Infection2.tar.gz | 90,659 x 280,000 | 280,000 | Complete column-ID list checked; no affected ID/token exists, so no count comparison is possible |

The combined Infection2 archive's cell-name file is explicitly named
`Postcapture_merged_Infection2_CellIDs.1.txt.gz`; that filename is unrelated to
the per-cell `.1` suffix. It contains unique cell IDs and no extra `.1` column
corresponding to the affected day-4 duplicate.

The streaming inspector verifies the MatrixMarket entry count, compares named
features and never treats column positions across objects as identities.
Infection1 has 5,868 matching source features; all remaining 66 source features
are zero in the affected cells. Infection2 has 5,913 matching features; all
remaining 21 source features are zero in the affected cells. Thus the zero
differences reported below account for the complete 5,934-feature source vectors.

### Fourteen same-group `.1` pairs

For every one of these pairs, both deposited Seurat columns match the same
single upstream GEO column exactly. The 13 Infection1 pairs correspond to 13
distinct post-capture columns. The Infection2 day-4 pair corresponds to combined
matrix column 33,197 (`D4Seq1_GCTGAATTTCGG`). There is no second independently
identified column for any `.1` copy in these deposited matrices.

For Infection1, the search also found eight other same-barcode-token columns in
different groups. Their count profiles differ from the affected source profile;
none offers an alternative exact match. Therefore the affected `.1` copies cannot
be explained by simply confusing a legitimate same-token column in a different
deposited array. Matching must use group plus barcode, not barcode alone.

Additional source-object evidence: the 13 Infection1 `.1` copies lie near the
end of the Infection1 block (positions 1,253–1,270, interspersed with other cells).
The day-4 `.1` copy is position 1,492. These copies have missing older metadata
such as orig.ident1 and sample.day, while their earlier counterparts generally
have those fields. Many older computational metadata fields also differ or are
missing; the audited experimental annotations and count profiles agree. In Table
S2 the paired parasite rows differ only in UMAP coordinates. Neither member is
independently represented in Table S2's host sheet under the corresponding key.

**Proven:** all 14 `.1` pairs are duplicate representations of single upstream
GEO count columns. The multiplicity changes from one named upstream count column to two
records in the authors' deposited final Seurat object and parasite metadata table;
PlaViSca preserves both count columns unchanged. Distinct-array identity is not
supported for these 14 pairs. Infection1 is post-capture only, excluding a simple
pre/post pair interpretation there. No matching called pre-capture barcode exists
for the day-4 pair either.

**Inference:** author-side concatenation, selection, reinsertion or renaming that
repeated an existing count column is the best-supported class of explanation.
The metadata pattern is consistent with addition of records lacking older
metadata, followed by downstream processing. A `.1` suffix is consistent with
name disambiguation, but neither the suffix nor this pattern identifies the exact
R command or establishes intent. No evidence of deliberately designed duplicated
biological observations was found.

### Cross-array day-5 pair

`D5Seq1_CCCCGATTGACG` and `D5Seq2_CCCCGATTGACG` have distinct source group names,
GEO samples and associated run sets. They are independently addressable as
columns 56,896 and 74,045, respectively, of the combined Infection2 matrix.

Those two upstream columns are exactly identical over **all 85,578 features**,
not only parasite genes. Each totals 612 counts across species, including 581
parasite counts. Their 5,934-feature parasite projections agree with both final
Seurat columns. They are final source positions 1,493 and 1,494 and differ in
only orig.ident among the 78 source metadata fields inspected. Both have missing
orig.ident1 and sample.day. Table S2 lists them separately as SEQ5r1 and SEQ5r2.

The pre-capture matrix has neither the full identifiers nor the CCCCGATTGACG
barcode token anywhere among 280,000 called columns. The day-4 GCTGAATTTCGG token
is likewise absent. This is absence from the deposited called-cell matrix, not
proof of no corresponding raw sequencing reads below its calling threshold.

**Proven:** the cross-array exact-profile duplication is already present at the
GEO combined-count stage, before the final parasite-only Seurat object. The
deposited files have two distinct column labels but supply no independent count
profile or pre-capture counterpart for either member.

**Inference/uncertainty:** a repeated profile or misassignment at or before
combined-matrix construction is more plausible than independent cells producing
identical full host-plus-parasite profiles. Barcode reuse alone does not explain
the identical full profile. The available combined matrix cannot distinguish
whether copying arose during pre/post merging, earlier library processing,
renaming, contamination or another error. Standalone post-capture-only Infection2
matrices and the authors' exact construction code are not in the inspected deposit.
Raw-read/UMI-level analysis or author clarification is needed to resolve this.

### What was and was not resolved

| Candidate explanation | `.1` pairs | Cross-array day-5 pair |
|---|---|---|
| Intentionally duplicated biological observations | No supporting evidence; intent unknown | No supporting evidence; intent unknown |
| Same barcode in genuinely distinct arrays | Contradicted by the shared matching source group/column for the pair | Distinct array labels exist, but count independence is not demonstrated |
| Repeated columns during processing/renaming | Supported class of explanation; localized between upstream matrix and final author object | Already present in GEO combined matrix; earlier operation unresolved |
| Separate pre/post representations of one cell | Not supported by deposited column identities; Infection1 has only post data | No matching called pre-capture column; simple two-representation explanation unsupported |
| Introduced by PlaViSca count export | Excluded by exact source/final count concordance | Excluded by exact source/final count concordance |

No deduplication or correction has been performed. The proposed analytic-inclusion
policy removes the 14 proven `.1` copies from a future analytic object while
preserving the full deposited inventory. The mechanism producing those copies
remains unknown. The day-5 cross-array pair remains flagged for a separate
decision: no survivor is chosen and biological independence is not demonstrated.

### Reproducibility and remaining limits

Run the diagnostic scripts from `pre_process_data/scripts`, always through:

```text
env -u R_LIBS_USER pixi run Rscript --vanilla ../audit/mancio_silva2022/inspect_lineage.R
env -u R_LIBS_USER pixi run Rscript --vanilla ../audit/mancio_silva2022/inspect_duplicate_internals.R
env -u R_LIBS_USER pixi run Rscript --vanilla ../audit/mancio_silva2022/inspect_upstream.R Postcapture_data_Infection1
env -u R_LIBS_USER pixi run Rscript --vanilla ../audit/mancio_silva2022/inspect_upstream.R Postcapture_merged_Infection2
env -u R_LIBS_USER pixi run Rscript --vanilla ../audit/mancio_silva2022/inspect_upstream.R Precapture_Infection2
env -u R_LIBS_USER pixi run Rscript --vanilla ../audit/mancio_silva2022/inspect_historical_lineage.R
```

Production input checksums are in `production_inputs_SHA256SUMS.txt`; audit
evidence checksums are in `SHA256SUMS.txt`. The historical source snapshots are
for inspection only: do not execute them directly because they contain saveRDS.
The historical diagnostic skips setwd/saveRDS and writes only audit TSV files.

The extra Zenodo CSV is explicitly advertised as SCT-transformed counts, not raw
RNA counts. Retrieval attempts through its content/file URLs returned HTTP 504;
it was not used as evidence. The local deposited RDS exposes the corresponding
SCT layers, whose duplicate concordance was inspected directly. Searches of the
deposit, publication links and exact source identifier/file names did not locate
the authors' study-specific cell-selection/renaming construction script. Generic
Drop-seq code is not a substitute for that missing provenance.

Diagnostic limitations/failures are recorded rather than hidden: one exploratory
inline source-internal assay check ended with an incomplete-expression parse
error after writing preliminary metadata evidence; the complete
inspect_duplicate_internals.R subsequently ran successfully. The initial
pre-capture inspector expected at least one candidate and stopped when none was
found; it was revised to emit an explicit absence result and completed. Neither
event modified production data. All conclusions above use completed comparisons.

Validation after production correction remains not applicable: no production
correction or regeneration was authorized or performed in this follow-up.


## Preliminary audit of the other PlaViSca studies and shared pipeline

**Cross-study addition, 2026-09-19; inspection and documentation only.** The
preceding detailed Mancio-Silva audit is unchanged. This addition checks the five
other study scripts, all six Hazzard2024 chunks and merger, shared integration,
SingleR, source-table and export code, existing study objects, and the deployed
`cleaned_dataset.rds`. It uses targeted publication/accession/deposit metadata.
It does not validate every input, alignment, count matrix, reference, orthology,
QC decision, embedding or biological annotation. No production script/data was
modified, no cells were removed, and no SingleR/Harmony/export was rerun.

The inspected preprocessing branch is `plavisca-audit-fixes`, HEAD `c2919da` at
inspection; production-code baseline remains `e5b3de5cfaa93d4e1a6ed8c8c4722dbcf1ec9144`.
Sol performed the bounded script review. R inspection used the prescribed Pixi
command from `scripts`, with Seurat 5.3.0. Historical comparisons read app Git
commit `c80580f` without checking out a branch. That commit is a matching historical
candidate, not a proven executed build. The central integrated object is absent.

Every finding below uses these categories:

- **CONFIRMED DEFECT:** a demonstrated code error, conflicting authoritative
  assignment, or deployed-data contradiction; code and deployed scope are stated.
- **STRONG SUSPICION / REQUIRES VALIDATION:** plausible scientific or lineage risk
  whose mechanism, biological validity or artifact impact is not established.
- **APPARENTLY CONSISTENT:** the named checks agree; this is not study validation.
- **NOT YET AUDITED:** no sufficient study-specific check was performed.

### Bibliographic/data-source table (separate from analytical metadata)

PubMed electronic dates, titles, DOI and PMID were checked for all five studies
(`cross_study_publications.tsv`). Publisher pages independently confirm the dates
and titles. The five script titles and DOIs agree with these publications, allowing
HTML formatting, terminal punctuation and whitespace. Years and PMIDs in the
study objects agree with the corresponding study; the source-table script itself
has no PMID column. Its `Number_of_cells` values are PlaViSca output counts, not
an independently established publication population.

| Study | PMID | DOI | Verified publication date | Publication/source cell-count meaning |
|---|---|---|---|---|
| Sa2020 | 32365102 | 10.1371/journal.pbio.3000711 | 2020-05-04 | Main analysis: 9,215 parasite transcriptomes after read-based QC; alternative lower-threshold analysis: 13,503 |
| Ruberto2022_1 | 36093191 | 10.3389/fcimb.2022.986314 | 2022-08-25 | 1,438 parasite transcriptomes across two infected-hepatocyte replicates |
| Ruberto2022_2 | 35926062 | 10.1371/journal.pntd.0010633 | 2022-08-04 | 9,947 salivary-gland sporozoites |
| Hazzard2022 | 36525464 | 10.1371/journal.pntd.0010991 | 2022-12-16 | Published sporozoite analysis: 2,609 + 2,363; exact total blood/source-list population not established here |
| Hazzard2024 | 39223117 | 10.1038/s41467-024-51949-8 | 2024-09-02 | 80,024 profiled parasites; deposited processed table also has 80,024 unique rows |

**CONFIRMED DEFECT, bibliography:** `data_source_manipulation.R:19` says `Sar2020`
instead of `Sa2020`; line 12 gives Hazzard2022 `16-Nov-22`, but the publication is
16 December 2022. The deployed `PlaViSca/data/data_source.csv` also has the wrong
November date and additionally lists that paper's journal as `PLoS Biol`, rather
than `PLoS Negl Trop Dis`. Its Ruberto author strings are `Anthony et al.` rather
than a surname-based citation, unlike the checked-out table script. Its schema
omits study/count columns, further showing that current source code and deployed
bibliography are not identical. These are bibliographic defects, not evidence of
incorrect counts or biological cell annotations.

Sources: [Sa2020](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3000711),
[Ruberto liver](https://doi.org/10.3389/fcimb.2022.986314),
[Ruberto sporozoites](https://journals.plos.org/plosntds/article?id=10.1371/journal.pntd.0010633),
[Hazzard2022](https://journals.plos.org/plosntds/article?id=10.1371/journal.pntd.0010991),
[Hazzard2024](https://doi.org/10.1038/s41467-024-51949-8).

### Cell inventories and scope of the boundary check

| Study | PMID | Study RDS / final cells | Established publication/source count | Inclusion authority in PlaViSca | Metadata risk | Annotation risk | Preliminary status |
|---|---|---|---|---|---|---|---|
| Sa2020 | 32365102 | 9,766 / 9,766 | 9,215 main analysis; 13,503 alternative | STARsolo GeneFull filtered cells, not an author-list subset | Handwritten library bridge, days/doses/animal IDs; geographic concept | Confirmed source-filter error; source-ID conversion suspect | CONFIRMED DEFECT plus REQUIRES VALIDATION |
| Ruberto2022_1 | 36093191 | 1,438 / 1,438 | 1,438 | Local Hep59 source-object IDs and LiverForm | Source-object chain of custody unresolved; parasite/host provenance separation | Source forms internally agree; blood-reference HPI inappropriate | APPARENTLY CONSISTENT selection, not fully audited |
| Ruberto2022_2 | 35926062 | 9,947 / 9,947 | 9,947 | Local combined object's Stage==Sporozoite | Local-object provenance and barcode bridge authority | Source-defined sporozoites; inappropriate HPI | APPARENTLY CONSISTENT selection, not fully audited |
| Hazzard2022 | 36525464 | 3,294 / 3,294 | 4,972 published sporozoites; blood total not established | Reprocessed raw GeneFull with emptyDrops | Registered instrument/protocol conflict; geography concept | Confirmed broad/detailed contradiction and HPI applicability | CONFIRMED DEFECT plus REQUIRES VALIDATION |
| Hazzard2024 | 39223117 | 80,024 / 80,024 | 80,024 processed rows, checksum validated | Deposited processed-cell list | Confirmed host-ID and strain errors; positional days; run-count semantics | Hard-coded cluster interpretation; study dominates integration | HIGH PRIORITY; CONFIRMED DEFECT plus REQUIRES VALIDATION |

For every study, pre-study IDs and final IDs are unique and identical in set and
within-study order: zero pre cells absent from final, zero final cells absent from
pre. This is an identity/inventory check, not a full raw-expression lineage check.
The objects have RNA `counts` only and no stored graphs/commands. Shared metadata
were compared by cell ID; changes in sample-type vocabulary, liver-form spelling,
Ruberto liver host wording and treatment terminology are recorded in
`cross_study_metadata_concordance.tsv`, not automatically classified as errors.
No assumption of positional equivalence was used to join separate objects.

### Sa2020: inclusion and source gametocyte annotations

**APPARENTLY CONSISTENT:** ENA PRJNA603327 contains exactly the ten used runs
SRR11008269–SRR11008278, all registered HiSeq 4000. Biosample strain/host and
sample aliases support the script's Saimiri/Aotus, NIH-1993/Indonesia-I/Chesson/AMRU-I
and CQ-bearing-library assignments. This does not validate exact host ID, infection
day, CQ dose/time or untreated/control semantics. Those finer assignments are
**NOT YET AUDITED** against original sample records/table.

`sa2020_pv_analysis_script.R` reads STARsolo `GeneFull/filtered`, removes rRNA,
uses `min.cells=1,min.features=1`, and saves the merged counts. Knee plots are
exploratory; `emptyDrops` is commented out and no author-cell-list subset is
active. **CONFIRMED DEFECT, terminology:** `10x_Chomium_V2` is misspelled; it does
not establish capture chemistry by itself. **REQUIRES VALIDATION, geography:**
`USA_Rockville` faithfully reflects registered `geo_loc_name`, but should not be
presented as original parasite geographic origin. It describes the deposited
experimental sample location; strain origins require a different provenance field.

The publication's 9,215 main-analysis transcriptomes follow PCR/read-based QC
with the higher 5,000-read cutoff and exclusion of very high-read droplets; its
13,503 alternative uses the lower 1,000-read cutoff. PlaViSca's 9,766 STARsolo
filtered cells use a different counting/QC definition. The net difference of
551 is **not by itself an analytical error** or proof of 551 additional cells.

Publication Fig 1B supplement `pbio.3000711.s030`, read in memory, contains
9,215 non-`Pberg` identifiers, including 267 `PB_MACS` records; it also contains
4,884 `Pberg` comparison cells. Candidate library correspondence inferred from
barcode overlaps is AMRU_Ao→278, AMRU_CQ→277, AMRU_Sa→276, Ches_Ao→275,
Ches_CQ→274, Ches_Sa→273, NIH_Ao→271, NIH_CQ→270, NIH_Sa→269, PB_MACS→272.
Under this explicitly **computational candidate bridge**, 9,018 author IDs map to
study cells, 197 author IDs do not, and 748 PlaViSca IDs are outside the mapped
main-analysis list. The mapped IDs are unique. The actual STARsolo filtered
barcode-row counts equal the study counts by run. **STRONG SUSPICION / REQUIRES
VALIDATION:** PlaViSca includes a different population; the source prefix/run bridge
must be independently validated before naming these as definitively extra/missing
biological cells. Barcode overlap alone is insufficient across separate arrays.
No inclusion change is proposed or applied in this preliminary check.

**CONFIRMED DEFECT, SingleR source-import code:**
`filter(type == c("Male gametocyte","Female gametocyte"))` is parity-dependent
recycled equality. It retains 789 of the spreadsheet's 1,571 sex-labeled rows
(749/1,477 female; 40/94 male), omitting 782. The local workbook matches the publisher download byte-for-byte
(checksum evidence saved). It has 9,158 rows, including
unmapped `PB_MACS`; this Fig 5 data table is not a complete 9,215-cell inclusion list.
Holding the production prefix conversion fixed, current filtering yields 355 final
ID matches versus 723 for `%in%`: 368 additional matches, of which 325 already agree
with final stage and 43 disagree. These are conditional code-impact measurements,
not a validated 43-cell biological correction list. All existing 355 matches agree
with final source labels (344 female, 11 male).

**STRONG SUSPICION / REQUIRES VALIDATION, source-ID conversion:** current code maps
both NIH_Ao and NIH_CQ to 269, maps NIH_Sa to 272, omits PB_MACS and has no 270/271
mapping. This disagrees strongly with the observed candidate bridge and with
269 being a registered Saimiri NIH sample rather than an Aotus/CQ library. Validate
author sample aliases and counts/read evidence before assigning corrected IDs;
otherwise same-token coincidences can create false annotation matches.

### Ruberto2022_1: liver source selection

**APPARENTLY CONSISTENT:** ENA PRJNA843856 has six runs; the four used infected
libraries are 609=d9/rep2, 610=d5/rep2, 611=d9/rep1, 612=d5/rep1, with two additional
uninfected controls 613/614. All register HiSeq X Ten. The publication supports
BGW BioIVT hepatocytes, Mondulkiri-derived parasites, v3 capture and d9 cultures
following MMV390048 treatment; these are distinct host-tissue and parasite-origin
concepts. Exact source-object rows agree with day/replicate mapping. The source
`condition` says `Infected`, not a detailed drug category: drug assignment draws
on the experimental protocol rather than that field alone.

The local `Hep59.1.2.seu_20aug2025.rds` has 1,438 unique IDs, with 1,147 Hypnozoites
and 291 Schizonts. Conversions 51inf→612,52inf→610,91inf→611,92inf→609 produce
1,438 unique IDs: zero unmapped/missing/excess IDs against study and final objects.
Final liver forms and final lifecycle labels agree for every cell, with 291
`Schizont (Liver stage)` and 1,147 Hypnozoite. This is **APPARENTLY CONSISTENT**
source-label propagation, not independent PlaViSca validation of liver identity.

**REQUIRES VALIDATION:** provenance of the dated local object. It contains prior
RNA/SCT processing commands, but object structure does not establish authorship
or identity with a publication deposit. The author repository and Zenodo 6463338
are authoritative starting points; the archive contains large output packages,
which were not downloaded. The inspected GitHub tree does not identify the dated
local filename. Preserve source LiverForm separately from SingleR/PlaViSca stage;
do not equate a liver Schizont with a blood IDC Schizont.

### Ruberto2022_2: source-defined sporozoites

**APPARENTLY CONSISTENT:** ENA PRJEB42435 contains ERR5087438–440 with rep1/2/3 aliases
and HiSeq X Ten. The local combined object has 9,947 Stage=Sporozoite and 8,304
BloodStage cells. Every selected sporozoite maps uniquely through Case_909→438,
Case_922→439,Case_923→440 to the study and final object, with zero missing/excess
cells. Source days are respectively 16,18,17 and counts 3,375/4,002/2,570, agreeing
with script assignments. The publication supports Cambodian Mondulkiri isolates,
Anopheles dirus salivary-gland collection and the 16–18-day window.

The all-Sporozoite override is **source/selection-defined**, not an independent
PlaViSca prediction. **REQUIRES VALIDATION:** chain of custody of
`6_PvSPZ.BS.combined.rds`, exact capture chemistry against original library protocol,
source filtering/count-reference equivalence and raw-expression reconstruction.
The authors' repository/Zenodo 6474355 are appropriate authorities; deposit manifests
were inspected, not downloaded in bulk. No full provenance proof follows from
9,947 matching the publication.

### Hazzard2022: confirmed stage contradiction, different QC

**APPARENTLY CONSISTENT:** ENA PRJNA863611 has eight runs: four PacBio and the four
used Illumina runs. Registered samples support 498=freeborni sporozoites,
499=stephensi sporozoites, 500/501=Saimiri boliviensis blood, and source strains
Sal1/Chesson for the two mosquito libraries, Sal1 for blood. Retained counts are
1,433/1,425/362/74. **CONFIRMED DEFECT, deployed annotation:** all 2,858 cells of
498/499 have detailed `life_cycle_stage=Sporozoite` but broad `parasite_stage=Blood
stage`. The app maps development_phase to parasite_stage, making this contradiction
user-visible. The source/run override corrects the detailed field without updating
the broad field; any future correction must occur upstream and propagate to exports.

The script actively runs seeded emptyDrops(lower=1,FDR<=0.001) on raw GeneFull
counts before object construction. The publication instead uses read/mapping-based
QC and reports 2,609 + 2,363 sporozoite transcriptomes above its low-read threshold.
**REQUIRES VALIDATION:** their relationship to the 2,858 retained PlaViSca
sporozoites; exact source-cell IDs, read/count definitions and blood population
remain unaudited. Different QC/alignment is established, not automatically invalid.

**REQUIRES VALIDATION, evidence conflicts:** script/study/final instrument is
NovaSeq6000 for every run, but ENA registers HiSeq2500 for 500/501 and NovaSeq6000
for 498/499; publication prose describes NovaSeq. Preserve registered instrument
and publication protocol separately pending reconciliation. Registered strain
labels also require reconciliation with publication prose describing NIH/Chesson
mosquito feeds; do not silently overwrite either authority. `USA_District of
Columbia` matches deposited sample location, not demonstrated original parasite
origin. Script days/host IDs are literal "NA" strings. Publication states sporozoite
collection 21 days post-feed, but precise run-level/time-origin mapping should be
verified before replacing missing day_post_infection; missing blood infection day
is not itself an error.

### Hazzard2024: high priority, inclusion authority does not validate metadata

**APPARENTLY CONSISTENT:** PMID39223117 and used accessions belong to the publication's
PRJNA1047651. ENA lists 43 deposited runs, all NovaSeq6000. The six chunk scripts
import 30 distinct runs without overlap: 955–956,959–965,966–972,975–976,980–987,
991–994 (prefix SRR27021). The merger's 24 source prefixes are exactly the 24 runs
represented in the 80,024-cell study/final object. Six imported runs contribute
no retained source IDs (959,961–965); 13 deposited runs are not imported by these
chunks. Their experimental roles/exclusion rationale are **NOT YET AUDITED**.
`num_srr=24` describes retained-library count, not all 43 deposited runs; record
these distinct study-provenance counts explicitly rather than recycling an
ambiguous ordinary cell metadata field.

**CONFIRMED DEFECT, vector construction; impact bounded:** every chunk's run_id
vector has n² entries because rep(vector,n) repeats the full vector. Loops use only
the first n entries, so this is not Mancio-style per-cell positional recycling.
66_72 has seven samples and eight days: loop assigns 31,30,29,31,31,30,24 and ignores
last29; no length-eight assignment occurs, so no execution error is implied.
91_94 has four samples but five strain entries, likewise ignoring its fifth entry.
These unchecked surplus vectors can conceal mapping errors. The duplicated
`run_id <- run_id <-` assignment is cosmetic. Source values are not validated
merely because execution succeeds.

**CONFIRMED DEFECT, reached study/final data and corroborated upstream:**

- SRR27021986 has host_id5537_3 for 3,489 cells; both accession alias and deposited
  processed source identify 5537_2.
- SRR27021993 has host_id5709_1 for 16,321 cells; both accession alias and deposited
  processed source identify 5708_1.
- All 16,321 cells of SRR27021993 carry strain NIH1993-F3, but registered strain is
  Chesson and deposited processed rows uniformly have Genotype=Chesson,
  Condition=mono_CH. The extra strain entry in 91_94 and applied fourth entry
  confirm the positional code defect's deployed effect.

`fed banimal 5164` is a confirmed source-code typo for imported run962, but run962
is excluded from the retained source list: it did **not** reach the checked
80,024-cell study/final metadata. Mixed NIH/Chesson strain and PNG/El Salvador
strings describe experimental infection/strain combinations rather than each
cell's genotype or a single collection place; exact authority and semantics of
these combinations, animal IDs and days require detailed audit. For example,
run972's processed cells include both NIH and Outcross genotypes: preserve
source cell genotype separately from library infection composition. Registered
Biosample host=Homo sapiens and USA:Bethesda for the three targeted libraries
conflicts with monkey experimental provenance/strain-origin metadata; retain
registry evidence separately rather than using it uncritically to change hosts.

**APPARENTLY CONSISTENT, source inclusion:** local `Proccessed_Data.txt` matches
Zenodo12775216 md5 `6ede1f7da4e5d57841c9f90883b2b2c1`. It has 80,024 rows and unique
source IDs; conversion produces 80,024 unique reconstructed IDs with zero
unmatched source or extra retained cells against study/final. All 24 retained
source groups are blood libraries; mosquito chunks do not contribute to this
retained processed population. Five chunk RDS objects remain: their union has
8,736,651 unique IDs, with no duplicates. They cover 75,996 final IDs; the missing
4,028 are run960, whose 59_65 chunk RDS is absent. Their 8,660,655 additional IDs
are pre-subset reconstructed candidates, not extra final analytical cells.
All 24 raw STARsolo barcode files are present; inventory is recorded. No full
raw-expression comparison or publication count-table/read provenance check was
performed. **NOT YET AUDITED:** expression equivalence, skipped-run rationale,
all days/animal provenance, and inclusion/QC derivation within the author analysis.

### Shared SingleR, HPI, lifecycle and gametocyte annotations

**CONFIRMED DEFECT, deployed HPI applicability:** all 105,963 final cells have a
nonmissing HPI label, and every study contains all seven labels 6,12,24,32,36,42,48.
This includes all 1,438 Ruberto liver cells, all 9,947 Ruberto sporozoites and the
2,858 Hazzard2022 sporozoites; the existing Mancio finding is independently
corroborated. The final vocabulary is harmonized (e.g. Host blood), while study
objects still use Mammalian host: blood. Current singleR.R masks with the latter
exact string, so vocabulary and execution order matter.

**Historical evidence refines the cause:** c80580f assigns HPI and detailed IDC
bins without any non-blood mask, and its late source overrides retain the HPI
label. Thus the artifact is consistent with an older unmasked pipeline; vocabulary
mismatch is a possible current/historical lineage risk, not the only proven cause.
The exact executed build is unknown. Current masks cannot be credited with having
protected the deployed artifact. Quantification by study/sample/phase/detailed
stage/HPI is in cross_study_hpi_stage.tsv and cross_study_final_stage_metadata.tsv.

The Zhu SMRU1 reference is a blood-stage time course. Its labels indicate reference
similarity, not measured infection age, and applicability to liver/sporozoite/sexual
forms is not established. Future meaningful HPI should be restricted to validated
appropriate asexual blood cells with explicit inference provenance; gametocyte
HPI also requires exclusion/applicability review. Preserve source labels and the
SingleR scores/delta/pruning confidence, not just hard labels. The prediction object
is transient; no explicit confidence persistence is present and no such fields were
found in the final metadata. **NOT YET AUDITED:** reference normalization,
orthology join cardinality/feature compatibility and inference performance.
0–<18 Ring,18–<30 Trophozoite,30–<46 Schizont,46–48 Merozoite are **PlaViSca heuristic
bins**, not directly validated universal lifecycle annotations.

**CONFIRMED DEFECT, graph argument name; actual behavior clarified:** current and
c80580f SingleR call FindClusters(graph_name="new_clustering"), but installed
Seurat uses graph.name. The underscore argument does not select that graph;
RNA_snn is the default. Repository code has no new_clustering graph creation;
integration creates RNA_snn from integrated PCA. A missing new_clustering graph
therefore does not by itself establish an execution failure. The central object
is absent, so historical graph identity/contents and exact versions are not
inspectable. **REQUIRES VALIDATION:** reproducibility of cluster numbering and
hard-coded cluster 2=female/cluster 4=male. Nineteen female/six male markers feed
AddModuleScore, but no explicit score threshold determines final cluster labeling.
Marker specificity and clustering biology remain unaudited.

Stored cluster-label counts and final gametocyte labels:

| Study | Cluster female / male | Final female / male | Source-defined propagation evidence |
|---|---|---|---|
| Sa2020 | 1,454 / 200 | 1,470 / 211 | 355 matched sex overrides (344F/11M); 328 already cluster female, 27 cluster asexual |
| Ruberto2022_1 | 20 / 0 | 0 / 0 | 1,438 source liver overrides; all 20 cluster-female cells end with liver forms |
| Ruberto2022_2 | 0 / 0 | 0 / 0 | 9,947 selection-defined sporozoite overrides |
| Hazzard2022 | 48 / 0 | 47 / 0 | One stored cluster-female mosquito cell ends as source/run-defined Sporozoite |
| Hazzard2024 | 20,292 / 0 | 20,292 / 0 | No independent sex-annotation authority demonstrated in this inspection |

Counts overlap across stages of assignment; they are not additive independent
inferences. Sa overrides are reconstructed through the suspect conversion, not
proof of a valid biological bridge. Liver labels are source forms, not a source
sex annotation. c80580f puts source liver/refine/sporozoite overrides after cluster
labeling; current code puts them before the cluster gametocyte overwrite, then
removes detailed male/female gametocytes with Anopheles hosts. This changed order
can lose source-defined labels and affects inclusion. **CONFIRMED code behavior:**
that subset is a real post-integration removal operation, not metadata cleanup.
**REQUIRES VALIDATION, executed impact:** zero cells are absent between any existing
study RDS and final rows, and no final mosquito gametocyte remains. One surviving
mosquito cell has stored cluster-female prediction and final Sporozoite; its exact
ID is saved. No historical removed-cell inventory is available, so deleted cells
cannot be reconstructed or an executed removal count claimed. No clustering or
annotation/removal was rerun for this audit.

Current new development_phase defaults unclassified detailed labels to Blood
stages and does not repair the old parasite_stage field consumed by the app.
Its source liver/sporozoite logic and obsolete Mancio label selector also differ
from deployed lineage. Preserve source LiverForm/State/sex labels, initial IDC
inference, cluster labels and final adjudication separately in a future repair.
This is a cross-study addition to the existing Mancio provenance recommendation,
not a change to the earlier Mancio conclusions or proposed inclusion policy.

### Integration: major scientific validation requirement, effective grouping caveat

The script merges six biologically disparate studies, normalizes RNA, selects 2,000
variable genes, scales/PCA, retains unintegrated embeddings and invokes
IntegrateLayers(method=HarmonyIntegration,group_by="study_label").
**Important verified API qualification:** installed HarmonyIntegration accepts
extra arguments but does not consume group_by; it constructs groups from assay
layers via CreateIntegrationGroups, then calls RunHarmony(vars_use="group").
The code expresses an intention to group by study, but that keyword does not prove
the actual grouping. Study objects each have one joined counts layer, so merged
layer membership may correspond to study; effective group-to-study correspondence
must be checked in the historical central object/build before asserting it.
Installed method source is captured; no integration was run.

**STRONG SUSPICION / REQUIRES VALIDATION, major scientific risk:** study/layer is
confounded with lifecycle, tissue, host, technology and experimental conditions.
This can remove true biology. Code inspection alone does not prove over-correction.
Future validation must compare integrated/unintegrated representations, preservation
of lifecycle/tissue marker patterns and trajectories, and batch mixing only within
biologically comparable subsets. Attractive mixing is insufficient. Use appropriate
RNA expression/counts, not Harmony/PCA/UMAP embeddings, as expression measurements
for differential expression. QC/source inclusion must be validated first.

### flatten_data.R: separate code defects from deployed export

**CONFIRMED DEFECT, current code and deployed output:** the loop constructs
combined top_genes_df, but save_data uses last-loop top_genes_exp instead. Deployed
cleaned_dataset$top_genes_exp has 14,940 rows/four columns, no study column, and its
cell IDs all belong to the 1,494 Mancio cells (10 gene rows/cell). Other studies'
top-gene summaries are missing from that saved component. Fix later at export,
not by manually patching the final RDS. hiv_data is a cosmetic HVF naming issue;
no additional analytical defect from that variable name was demonstrated.

**APPARENTLY CONSISTENT, scale layer spelling under installed API:** current
GetAssayData(layer="scale") resolves scale.data in a small in-memory Assay5
inspection under the Pixi versions. It is not automatically a v5 execution error.
The historical central assay is absent: actual historical scale layer/default
assay and scaled-value provenance remain **NOT YET AUDITED**. No production assay
was scaled or changed for this API check.

**CONFIRMED code/artifact divergence:** current export removes pred_gametocyte,
parasite_stages,liver_form,refine_state, but deployed mr_data still includes
pred_gametocyte,liver_form,refine_state and life_cycle_stage, not the current
parasite_stages/development_phase schema. Original metadata names are collected
before clean_names and then selected with any_of(metadata), potentially omitting
renamed columns, as previously traced for Mancio. Hard-coded old working paths
are present across scripts. These reproducibility/code concerns are separate from
valid expression in existing artifacts; non-Mancio expression values were not
fully compared in this preliminary audit.

### Correction plan, regeneration dependencies and next forensic audit

No correction is applied. Future small reviewed repairs should address source-table
bibliography, source-keyed host/strain fields, Sa library/sex annotation joins,
HPI applicability/provenance, broad/detailed stage consistency, graph/group API
arguments and combined top-gene export. First establish the correct executed
lineage and authoritative source bridges. Repair at the earliest responsible step:
study metadata constructions → study objects → central object/annotation → exports.
Cell/QC changes require recomputation of affected normalization/integration and
annotations; metadata-only changes should not be portrayed as integration validation.
Export-only summary errors belong in export. Preserve originals/source evidence and
separate harmonized/inferred fields. Validation after correction is not applicable;
validation performed here is documented by the evidence below.

Recommended full-audit order remains **Sa2020 first**, then **Hazzard2024**,
**Hazzard2022**, **Ruberto2022_1**, **Ruberto2022_2**. Sa has a demonstrated sex-filter
error plus a source-library bridge/inclusion discrepancy affecting biological
conclusions; validate that bridge and author QC before repairing annotations.
Hazzard2024 is equally urgent operationally because it dominates the atlas and
already has authoritative, high-impact metadata defects. Hazzard2022 follows for
its confirmed user-visible stage contradiction and QC/platform conflicts. The
Ruberto selections appear internally coherent but still need deposit lineage and
source-label validation. No study is declared publication-ready.

### New evidence and verification limits

All newly created files are named `cross_study_*` in this audit directory; exact
filenames, byte sizes and SHA256 are listed in `cross_study_evidence_manifest.tsv`.
`cross_study_documentation_validation.txt` records append-only preservation of the
previous audit and tracked production-code status. The pre-existing SHA256SUMS.txt
retains hashes from the prior audit package; its AUDIT.md hash predates this clearly
marked addition. Use the new cross-study manifest for this revision, without
rewriting old Mancio evidence.

Key evidence: cell_summary; run_metadata/final_run_metadata; metadata_concordance;
hpi_stage/final_stage_metadata; gametocyte/mosquito_cluster_gametocyte_candidates;
liver_mapping/liver_source_groups; spz_mapping/spz_source_groups;
hazzard2024_mapping/run_membership/source_conditions/existing_chunk_summary;
sa_source_filter/filter_summary/author_fig1B_ids/author_barcode_overlap/
candidate_lineage/candidate_additional_cells; publications and five ENA tables;
accession_sample_attributes; deposit_manifest; source_checksum_validation;
retrievals; two historical script snapshots; inspection source/log text files.
Author repositories were queried for file provenance without bulk downloads;
missing local dated-object identities remain unresolved. No deposited count archive
was downloaded, and no production object was regenerated. Evidence scripts only
read production files and write these audit tables/logs; the sole synthetic assay
was for API resolution, not a scientific reanalysis.
