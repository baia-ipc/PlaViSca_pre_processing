# Shared-pipeline and cross-study consolidation audit

Audit date: 2026-09-19. Status: consolidation complete; no production repairs applied.

This document is the bridge between the completed study-specific forensic-audit
phase and the future implementation/repair phase. It does not repeat the six
study audits; it establishes the authoritative defect register, the shared
mechanisms that produced cross-study errors, a provenance model for metadata,
the exact repair required at each pipeline layer, what must be regenerated,
what needs a human scientific decision, and objective acceptance tests.

Repository: `pre_process_data` (git, branch `plavisca-audit-fixes`), a nested
repository inside `/home/baia/prj/plavisca` (which is itself not a git
repository — only `pre_process_data` is). At the time of this audit,
`git status` showed a clean tree except for one untracked `data/` directory
(present before this audit started, left untouched), and `git log` showed ten
commits from `e5b3de5` (initial commit) through `9e6f947` ("Complete
Ruberto2022_2 forensic audit"), the last of the six study audits.

Runtime for all inspections in this consolidation:
`env -u R_LIBS_USER pixi run Rscript --vanilla ...` from `pre_process_data/scripts`,
R 4.3.3 / Seurat 5.3.0 / SeuratObject 5.2.0, matching the convention used by the
prior study audits. No production script, source data, study object, or
application code was modified. No cells were removed, no SingleR/Harmony/export
step was rerun against real data beyond small, read-only, dimension-only
diagnostics.

## Important caveat: Hazzard2022 is partially, not fully, audited

Unlike the other five studies, no dedicated `audit/hazzard2022/AUDIT.md` or
equivalent full publication/accession/QC/provenance audit file exists anywhere
in this repository. However, commit `ea006124` ("audit Hazzard2022 cell and
expression lineage") **already completed a deterministic forensic
cell/expression lineage audit**, adding `hazzard2022_cell_lineage.tsv`,
`hazzard2022_lineage_summary.tsv`, `hazzard2022_expression_concordance.tsv`,
`hazzard2022_raw_inventory.tsv`, `hazzard2022_exact_duplicate_profiles.tsv`,
annotation/metadata evidence, and the deterministic audit scripts
`audit_hazzard2022.R`/`audit_hazzard2022_emptydrops.R`. This evidence proves,
for all 3,294 cells: exact source→study-RDS→deployed cell lineage with zero
missing/extra cells at any boundary; zero study-RDS-vs-deployed raw and
normalized expression mismatches; zero exact-duplicate expression profiles;
and full raw-STARsolo-barcode membership for every retained cell. It is
therefore false to say that only a preliminary paragraph exists, or that
exact lineage/expression/duplicate checks were never performed.

What remains genuinely incomplete is the **publication/accession/QC/provenance
audit**: the only evidence for the sequencer conflict (`D023`), the literal
`"NA"` day/host strings (`D024`), and the QC/inclusion-authority
correspondence between PlaViSca's 2,858 retained cells and the publication's
2,609+2,363 (`D022`) is still the short preliminary/cross-study paragraph
inside `audit/mancio_silva2022/AUDIT.md` (lines 1017–1046), predating the
full-forensic-audit methodology applied to Sa2020, Hazzard2024, Ruberto2022_1,
and Ruberto2022_2.

The correct classification is therefore: **forensic cell/expression lineage
audit completed; authoritative publication/accession/QC/provenance audit
remains incomplete.** This consolidation treats Hazzard2022's lineage and
expression-concordance figures (3,294 cells, zero mismatches) as final, and
treats only the publication/QC/provenance figures (2,858 vs. 2,609+2,363, the
sequencer conflict, day/host provenance) as preliminary, and registers
completing that remaining reconciliation as a blocking decision (`DEC08`,
defect `D025`) before any Hazzard2022 QC/inclusion-driven repair is
finalized. This is a genuine, narrower gap in the forensic-audit phase, not
an oversight of this consolidation task, and it should be surfaced to the
user/project lead explicitly.

## 1. What is wrong?

Sixty-nine confirmed defects, source/provenance issues, reproducibility
defects, and export/UI defects are registered in `confirmed_defects.tsv`
(D001–D069), spanning:

- **Per-study metadata construction bugs** (Mancio-Silva2022 D001–D010,D013;
  Sa2020 D014–D017,D019; Hazzard2022 D020,D023,D024; Hazzard2024
  D026–D030,D032,D033; Ruberto2022_1 D038,D039; Ruberto2022_2 D040–D043):
  recycled/positional vectors standing in for keyed joins, discarded source
  metadata, conflated geography, misspelled/underspecified technology fields,
  ambiguous per-cell run counts, and literal `"NA"` strings.
- **A systemic, copy-pasted `"10x_Chomium_V2/V3"` misspelling** appearing
  independently in Sa2020, Hazzard2024, and Ruberto2022_1 (D069) — the same
  shared-constant bug, not three unrelated typos.
- **Shared annotation-pipeline defects in `singleR.R`** (D047,D048,D049,
  D058,D059,D060,D061; Parts 3, 5, 6 below): an invalid `FindClusters`
  argument name that silently falls back to the wrong graph; hard-coded,
  demonstrably unstable numeric cluster-to-sex-label mapping; a blood-stage
  reference applied blindly before any population gating; unsupported
  heuristic stage bins presented without qualification; no persisted
  inference uncertainty; and a post-integration cell-removal filter that is
  a latent, unlogged data-loss risk.
- **A shared integration defect in `integration.R`** (D050,D051; Part 7
  below): `IntegrateLayers(..., group_by = "study_label")` is silently
  ignored by both `IntegrateLayers` and the installed `HarmonyIntegration` —
  neither has a `group_by` formal parameter. The actual Harmony batch
  grouping is derived from per-merged-object layer structure via
  `CreateIntegrationGroups`, which today happens to have the same cardinality
  as the six studies, by coincidence of how the objects were merged, not by
  design.
- **Shared export defects in `flatten_data.R`** (D052–D055; Part 8 below): a
  variable-name bug that discards five of six studies' top-gene data; an
  `any_of(metadata)` selection evaluated against pre-`clean_names()` column
  names that silently drops Seurat-default metadata columns; proven
  divergence between current export code and the deployed `cleaned_dataset.rds`
  schema; and positional (non-keyed) `cbind` joins between metadata and
  expression matrices with no equality assertion.
- **Two major expression-layer defects**: Mancio-Silva2022's 14 proven `.1`
  duplicate representations of single upstream GEO count columns plus one
  unresolved cross-array duplicate pair (D011,D012 — scientific decisions,
  not code bugs); and Ruberto2022_1's near-total expression loss for 538/1,438
  biological-replicate-2 cells (D036 — the single most severe defect in the
  project, requiring FASTQ/alignment-level investigation).
- **Application-layer defects** (D062–D069; Part 9 below): literal `"NA"`
  strings surfaced as ordinary categories; dead code hardcoding non-existent
  column names that directly encodes the same plural/singular naming trap
  found in the export layer; hardcoded label/color maps that bypass the
  app's own config-indirection safety mechanism; and an unvalidated
  bibliography table.

## 2. What is merely uncertain?

- The **exact executed build** of every deployed artifact. No build manifest
  or central `pv_all_studies.rds` exists anywhere in the repository — this is
  the single biggest reproducibility gap in the whole pipeline (see
  `pipeline_dependency.tsv`, "final metadata harmonization" row). Multiple
  independent lines of evidence (historical script diffs, schema divergence,
  vocabulary drift) prove the deployed `cleaned_dataset.rds` was **not**
  produced by a fresh execution of the current checked-out scripts, but the
  exact historical commit/version mix cannot be recovered.
- **Mechanism** (not consequence) of several defects: the author-side process
  that created Mancio-Silva2022's 14 `.1` duplicates and the D5 cross-array
  pair (D011,D012); the root technical cause of Ruberto2022_1's replicate-2
  expression loss (D036 — consequence proven exactly, cause not yet
  investigated at the FASTQ/alignment level).
- **Scientific validity of the current global Harmony integration** given
  severe study/biology confounding (Hazzard2024 alone ≈75.5% of cells) — not
  provable or disprovable from code alone; requires the quantitative
  validation metrics in `acceptance_tests.tsv` (AT24–AT27), not visual UMAP
  inspection.
- **Hazzard2022's true retained-cell population and its correspondence to
  the publication's QC** (D022, DEC04), blocked on the remaining
  publication/accession/QC/provenance audit (D025, DEC08) — the forensic
  cell/expression lineage figures themselves are no longer uncertain (see
  §3 below).
- Several **source-provenance questions flagged but not resolved** by the
  original study audits and preserved here unchanged: Sa2020's exact
  NIH-1993-F3-equivalent strain origin question surfaces again in
  Hazzard2024 (NIH-1993-F3's relationship to Salvador-I is "closely related
  to but distinct from," not identical); Mancio-Silva2022's per-infection
  parasite province and Infection2 sequencer capture-role.

## 3. What is already trustworthy?

- **Cell lineage/identity** (source → study RDS → deployed export) is proven
  exact for all six studies: zero missing/extra cells at any boundary, for
  Mancio-Silva2022 (1,494), Sa2020 (9,766), Ruberto2022_1 (1,438),
  Ruberto2022_2 (9,947, the cleanest result in the project — zero mismatches
  at every boundary), Hazzard2022 (3,294, forensic lineage audit complete —
  see the caveat above; only its publication/QC correspondence remains
  preliminary), and Hazzard2024 (80,024).
- **Raw/normalized expression concordance boundaries must not be
  overstated.** What is proven, precisely, differs by study:
  - Mancio-Silva2022 and Ruberto2022_2: normalized-expression concordance
    with the stated formula (`log1p(count/cell_total*10000)`) is proven
    exact, and raw counts are proven byte-identical to source.
  - **Hazzard2024** proves study-RDS → deployed raw/normalized expression
    concordance (zero mismatches) and processed-list/study/deployed cell
    lineage identity, but does **not** prove author-deposited count-table →
    PlaViSca concordance: the author count-table ZIP (advertised MD5
    `276aa89b96636e9091a05cd43166465e`) could not be downloaded (repeated
    Zenodo requests returned HTTP 504), so study-RDS counts have never been
    checked against the authors' own deposited values.
  - **Hazzard2022** proves study-RDS → deployed raw/normalized expression
    concordance (zero mismatches, `hazzard2022_expression_concordance.tsv`)
    and raw-STARsolo-barcode membership for every retained cell
    (`hazzard2022_raw_inventory.tsv`), but does **not** currently prove
    STARsolo raw count *values* → study-RDS equality (only barcode presence
    was checked against the raw matrix, not elementwise count values), and
    comparable author-deposited raw/processed matrices were not available to
    check PlaViSca's independent emptyDrops reprocessing against the
    authors' own processing.
  - Ruberto2022_1's expression is proven exact **except** for the 538
    replicate-2 cells (D036).
- The **Zhu SMRU1 → PvP01 orthology join itself** is a clean 1:1 mapping (no
  duplicate `InputOrtholog`/`GeneID`, no one-to-many/many-to-one artifacts) —
  the reference-construction problem is about labeling/gating and an
  unlogged `drop_na()` count, not about join-cardinality corruption.
- The **`development_phase`/`parasite_stages` two-dimension filter design**
  in the live (non-dead) app code is a genuinely good separation of concerns,
  and the app's `get_data_column()`/`get_data_label()` config-indirection
  layer, where actually used, is a good pattern worth keeping and extending.

## 4. What must be corrected, at what earliest pipeline step, and what must be regenerated?

See `confirmed_defects.tsv` (per-defect detail), `pipeline_dependency.tsv`
(the effective dependency graph, stage by stage, with reproducibility
status), and `regeneration_matrix.tsv` (the exact regeneration footprint of
every defect/repair). In summary, ordered by dependency:

1. **Study-preprocessing metadata fixes** (most of D001–D046) are
   metadata-only: they require regenerating the affected study RDS and every
   downstream export, but **not** re-normalization, re-integration, or
   re-clustering, per the stated principle that a metadata-only change does
   not automatically require recomputing expression embeddings, *provided*
   cell IDs and counts remain unchanged. That said, a metadata-only fix is
   not free of merged-object impact: because the deployed exports are built
   from the merged atlas object's metadata, every corrected study-RDS field
   must still propagate into the merged/final atlas object's metadata before
   a corrected export can be produced — this is a metadata reconstruction/
   update of the merged object, not a normalization/HVF/PCA/Harmony rerun,
   and `regeneration_matrix.tsv` marks it explicitly as "Yes — metadata
   reconstruction/update only" rather than "No", so that a corrected export
   is never assumed to be producible from a stale merged-object metadata
   snapshot.
2. **Cell-membership-changing scientific decisions** (D011,D012,D018,D022,
   D036) require full downstream regeneration for the affected study —
   and, because the atlas is globally integrated, potentially for the whole
   merged object — once approved.
3. **Ruberto2022_1's expression repair (D036)** is the most consequential
   single item: because it changes actual count values for 538 cells, it
   requires FASTQ/alignment-level work before any downstream regeneration,
   and blocks trusting any current Ruberto2022_1 clustering/integration
   result in the meantime.
4. **Shared `singleR.R` fixes** (D047,D048,D049,D058,D059,D060,D061) require
   one atlas-wide rerun of clustering/SingleR/gametocyte-inference — this
   should happen exactly once, after all per-study cell-membership and
   expression decisions are frozen (see Part 13 staged strategy below), not
   once per defect.
5. **The `integration.R` Harmony-grouping fix** (D050,D051) requires a full
   atlas-wide Harmony rerun, and should happen before the shared `singleR.R`
   rerun (steps depend on the corrected integrated reduction).
6. **`flatten_data.R`/`data_source_manipulation.R` export fixes**
   (D052–D057) are cheap, export-only regenerations, but only meaningful
   once `pv_all_studies.rds` exists again (it currently does not) and the
   upstream repairs above are complete.
7. **Application-layer fixes** (D062–D069) require no preprocessing rebuild
   at all, except where they depend on the upstream "never emit literal NA"
   rule (D062) or the corrected config-key naming (D063,D064).

## 5. What requires a human scientific decision?

Twelve decisions are registered in `decision_register.tsv` (DEC01–DEC12),
covering: the two Mancio-Silva2022 duplicate-inclusion questions (DEC01,
DEC02); Sa2020's inclusion-population authority (DEC03); Hazzard2022's
QC/inclusion authority and sequencer-conflict resolution, both blocked on
completing its remaining publication/accession/QC/provenance audit (DEC04,
DEC05, DEC08 — the forensic cell/expression lineage audit is already
complete); Hazzard2024's
biological-replicate model and unresolved strain-origin geography (DEC06);
the sex-annotation authority question shared by Hazzard2024 and Sa2020
(DEC07); Ruberto2022_1's expression-repair strategy (DEC09); the
post-integration removal-filter policy (DEC10); the global-vs-stratified
integration question (DEC11); and the cluster-to-sex hard-coding redesign
(DEC12). None of these has been resolved by this consolidation; each carries
alternatives, evidence, downstream consequence, and (where evidence
supports one) a recommendation, per the task's instruction never to convert
an unresolved scientific question into an automatic correction.

## 6. What exact tests will define a successful repaired build?

Thirty machine-testable acceptance criteria are registered in
`acceptance_tests.tsv` (AT01–AT30), covering identity/lineage integrity,
metadata-vector-length and keying discipline, the "never literal NA" rule,
source/inferred provenance separation, biological-validity gates (non-blood
HPI, sporozoite broad-stage labeling, liver-cell gametocyte exposure,
zero-expression annotation), study-specific corrected-value checks, export
key-alignment discipline, application-contract documentation, and — for
integration specifically — quantitative validation metrics (batch mixing
within biologically comparable populations only, marker preservation, study
dominance checks) rather than visual UMAP inspection alone. Cluster
reproducibility is deliberately split into two tests: same-seed determinism
(AT29) and marker/population correspondence under seed/resolution
perturbation (AT30), so that a numeric cluster ID is never treated as a
stable biological identity by either test alone.

## Parts 3–9: detailed technical audits

The following sections summarize the deep technical audits performed for
this consolidation. Full code citations (file:line) are preserved in
`confirmed_defects.tsv`'s `evidence` column and in the underlying agent
findings; this section gives the narrative conclusions.

### Part 3 — `singleR.R` (IDC/HPI reference, uncertainty, heuristic bins)

The Zhu SMRU1 → PvP01 reference construction (`singleR.R:20-72`) is a clean
1:1 orthology join (4,912-row table, 0 duplicated keys) that silently drops
318 of 5,226 Zhu SMRU1 genes via `drop_na()` with no logged count, yielding a
**4,908-gene × 7-timepoint reference** whose values (-4.0 to 15.2) are
single-replicate, log2-scale bulk microarray intensities stored under the
misleading `logcounts` assay name — not scRNA-seq-comparable counts (D059).
Test/reference intersection is ≈4,848/4,908 genes (91–99%) for one
representative study.

The values `6,12,24,32,36,42,48` are SingleR's classification labels — "which
reference timepoint profile does this cell correlate best with" — copied
verbatim into a field named `hour_post_invasion` with no interpolation or
confidence weighting (`singleR.R:67-80`). This is not a measured biological
timestamp (D058). Classification is run on the **entire combined dataset**
with no pre-filter by tissue/population type; gating happens only afterward
via post-hoc `NA`-masking on a `sample_type` string ("Mammalian host: blood")
that does not match every study object's actual vocabulary, so historically
built objects can retain HPI for non-blood cells regardless of the current
masking code's intent.

Only the final pruned/binned label is ever persisted (D061); SingleR's raw
label, score matrix, and `delta.next` exist only transiently and are never
written to metadata or the saved `.rds`. The heuristic stage bins
(0–<18 Ring, 18–<30 Trophozoite, 30–<46 Schizont, 46–48 Merozoite,
`singleR.R:84-93`) carry no citation anywhere in the script, the audit
documentation, or the manuscript, and are confirmed PlaViSca-invented
heuristics, not Zhu et al. 2016-sourced boundaries (D060).

### Part 4 — Lifecycle-annotation precedence

`annotation_precedence.tsv` gives the deterministic, per-study/population
precedence rule the repaired pipeline must follow, distinguishing
`source_life_cycle_stage`, `idc_reference_similarity_label`,
`pred_gametocyte_sex`, `harmonized_life_cycle_stage`, and an optional
`adjudicated_display_stage` (Part 10's four-tier schema). The central
finding is that **precedence is currently order-dependent and undocumented**:
the historical script applied source/liver/sporozoite overrides *after*
cluster-based gametocyte labeling, while the current script applies them
*before* it, then removes detailed labels for Anopheles-host cells — a
silent change in which source ever wins over an inference (shared root cause
of D045, D048, D049).

### Part 5 — Gametocyte annotation

`FindClusters(graph_name = "new_clustering")` (`singleR.R:183-187`) uses an
argument name (`graph_name`) that does not exist on the installed Seurat
5.3.0 `FindClusters.Seurat` (the real argument is `graph.name`); the bogus
argument is silently absorbed via `CheckDots()` (a warning, not an error),
so the call falls back to the default `"RNA_snn"` graph already built by
`integration.R`. **No graph literally named `"new_clustering"` is ever
created anywhere in the repository** (D047). This is a confirmed
code/reproducibility defect — the intended graph-selection argument does not
work as written — but it does **not** by itself establish that the fallback
`"RNA_snn"` graph was biologically the wrong graph to cluster on; whether
clustering fidelity was actually compromised requires the marker/population
validation in `AT29`/`AT30`, not an assumption drawn from the API mismatch
alone. Cluster `"2"` = Female
gametocyte and cluster `"4"` = Male gametocyte are hard-coded in
`RenameIdents` (`singleR.R:296-320`); the historical script version instead
had cluster `"5"` = Male, cluster `"4"` = Asexual — **proving cluster
indices already shifted between two versions of the same repository** on
comparable input (D048). Given Seurat/Louvain cluster numbering is inherently
order/seed/resolution-dependent, this hard-coding is not reproducible after
any pipeline rebuild without manual re-inspection.

Cross-study source-vs-cluster gametocyte disagreement is substantial and
quantified: Sa2020 188/1,532 present-cell disagreements after correcting the
import bugs; Hazzard2024 ≈2,387 discordant cells between author
`PseudoGroup=Sexual` (21,009) and deployed Female-gametocyte calls (20,292).
This consolidation does not rerun clustering to validate the existing
labels; it specifies (DEC07, DEC12) what must be redesigned/revalidated.

### Part 6 — Post-integration cell removal

`singleR.R:460-465` removes any cell with `parasite_stages` in
(Male/Female gametocyte) **and** `host_species` matching `"Anopheles"`. Today
this removes **zero** cells, and an earlier, separate override block
(`singleR.R:142-146`) does force-relabel the two known Anopheles-host
populations (Ruberto2022_2, Hazzard2022 runs 498/499) to "Sporozoite" before
this filter runs. However, that earlier override does **not** structurally
protect these cells: a later, unconditional gametocyte-cluster override
(`singleR.R:329-331`, D045) runs *after* the Sporozoite override and *before*
the deletion filter, and can silently re-overwrite any Sporozoite-labeled
cell back to Male/Female gametocyte purely from its expression-cluster
membership — with no exclusion condition checking whether the cell already
carries a source-selection-defined stage. Zero cells are removed in the
deployed build today only because none of these source-defined mosquito-host
cells currently lands in the female/male gametocyte cluster; this is a
property of today's data, not a property of the code's ordering. Ruberto2022_2's
own audit (`D045`) already describes this mechanism correctly; the ordering
here reconciles `D049` to the same finding. Any future or currently-known
Anopheles-host population remains fully exposed to silent, unlogged deletion
the moment it (or a rerun's clustering) lands in the gametocyte cluster
(D049). The repaired policy (DEC10) must therefore fix `D045` and `D049`
together: gate the removal filter on source-defined population type, not
solely on downstream classifier output; scope the gametocyte-cluster
override to skip any cell already carrying a source-selection-defined stage;
and log every removed cell. Because zero cells are removed in the currently
deployed build, fixing or redesigning this filter does not itself require
any cell-membership regeneration today; the requirement is CONDITIONAL/
POSSIBLE, not certain — it would only materialize on a future rerun where
clustering assigns a source-defined mosquito-host cell to the gametocyte
cluster while `D045` remains unfixed.

### Part 7 — Harmony/integration

Confirmed by direct inspection of installed package source:
`IntegrateLayers` and `HarmonyIntegration` have **no `group_by` formal
parameter** — `group_by="study_label"` in `integration.R:43-50` is silently
absorbed into `...` and never read. The actual Harmony batch variable comes
from `CreateIntegrationGroups`, which derives one group per data layer
present on the assay (from `merge()`'s internal `cells` slot), **not from
any metadata column** (D050). Because `JoinLayers()` runs *after*
`IntegrateLayers()` (`integration.R:52` vs `:43-50`), the six per-study
objects remain as six unjoined layers at integration time, so
`CreateIntegrationGroups` happens to return six groups — coincidentally
matching the six studies today, but not by design, and would behave
differently (or error) if `JoinLayers` timing, merge order, or the number of
pre-existing layers per input ever changed (D051). This is a confirmed
API/intent mismatch and fragility, not by itself confirmed proof that the
atlas was integrated with the wrong grouping: because the six layers present
at `IntegrateLayers()` time currently map one-to-one to the six studies, the
effective layer-derived grouping may already equal the intended per-study
grouping today. The requirement to make this grouping explicit and testable
stands regardless of whether current output happens to be correct. Given Hazzard2024's
~75.5% atlas share and simultaneous confounding across lifecycle/tissue/
technology/condition, one global integration's validity cannot be asserted
from code alone — it requires the quantitative tests in AT24–AT27
(DEC11 registers the global-vs-stratified question as an open decision).

### Part 8 — `flatten_data.R` and export

Four confirmed defects: (1) `save_data` stores the stale last-loop-iteration
`top_genes_exp` instead of the correctly-accumulated `top_genes_df`, so
deployed `cleaned_dataset$top_genes_exp` contains only Mancio-Silva2022's
14,940 rows, missing all five other studies (D052); (2) the `metadata`
column-name vector is captured *before* `clean_names()` but `any_of(metadata)`
is evaluated *after* it, silently dropping Seurat-default columns like
`orig.ident`/`nCount_RNA`/`nFeature_RNA` from every expression export
(D053); (3) current export code removes `pred_gametocyte`/`parasite_stages`/
`liver_form`/`refine_state` columns that the deployed `cleaned_dataset.rds`
still contains under the old 52-column schema, proving the deployed artifact
predates the current export code (D054); (4) `red_df`/`norm_exp`/`raw_exp`/
`scale_exp` are each independently sorted by rowname then combined via
positional `cbind` with no explicit key-equality assertion — safe only if
cell-ID sets are identical, which is assumed, never verified (D055).
`study_label` spelling is also confirmed inconsistent across scripts
(`"Sar2020"`, `"Mancio Silva2022"` vs `"Mancio-Silva2022"`, D056), and
`data_source_manipulation.R`/`data_source.csv` carry independently confirmed
bibliographic errors (D057).

### Part 9 — Application data contract

See `app_data_contract.tsv` for the full field-by-field mapping. Highest-signal
findings: literal `"NA"` strings cover 88% of `treatment` and ~3% of
`day_post_infection` and are never distinguished from real missingness
anywhere in the app (D062); a sourced-but-dead module
(`metadata_filter_module.R`) hardcodes references to `mr_data$development_phase`
and `mr_data$parasite_stages`, **neither of which exists** in the shipped
schema — a live instance of the same plural/singular naming trap found in the
export layer, currently inert only because nothing calls it (D063); two
hardcoded label/color maps in `server.R` duplicate the config-driven
`get_data_label()` logic by key on physical column name, and would silently
desync if `config/data_columns.yml` were ever repointed (D064); an
unreachable `"blood_stage"` NA-context branch references a column that does
not exist (D065); the "Female gametocyte" color override has no "Male
gametocyte" counterpart (D066); and the bibliography table
(`data/data_source.csv`) has no programmatic join or validation against
`mr_data$study_label`/`study_pmid`/`pub_year`, so the two already-confirmed
string mismatches (D056,D057) are currently undetectable by the app itself
(D067). This confirms the task's premise directly: the Hazzard2022
broad/detailed-stage contradiction (D020) and these app-layer findings both
show that preprocessing-schema changes cannot be designed without checking
the application consumer at the same time — hence `metadata_schema.tsv` and
`app_data_contract.tsv` are designed together in this consolidation.

## Parts 10–14 pointers

- **Part 10 (target metadata schema):** `metadata_schema.tsv`.
- **Part 11 (automatic fixes vs. scientific decisions):** the
  `human_decision_required` column in `confirmed_defects.tsv` plus the full
  `decision_register.tsv`.
- **Part 12 (regeneration dependency matrix):** `regeneration_matrix.tsv`.
- **Part 13 (staged repair strategy):** see below.
- **Part 14 (acceptance tests):** `acceptance_tests.tsv`.

### Staged repair strategy (Part 13)

Refining the suggested phases against the actual dependency graph established
above:

- **Phase A — Resolve blocking scientific decisions.** At minimum DEC01,
  DEC02, DEC03, DEC08 (which gates DEC04/DEC05), DEC09. No production
  rebuild starts before these are settled, because each can change cell
  membership or expression values.
- **Phase B — Fix source tables and study preprocessing.** All per-study
  metadata-only defects (D001–D010,D013–D017,D019,D020,D023,D024,D026–D030,
  D032,D033,D038–D043), using keyed metadata tables, never positional
  vectors. Independent per study; can proceed in parallel with Phase A for
  studies whose decisions are already resolved (Sa2020's automatic fixes,
  Hazzard2024's automatic fixes).
- **Phase C — Repair/regenerate defective expression datasets.**
  Ruberto2022_1's replicate-2 expression (D036), following DEC09's chosen
  strategy. This is the highest-risk, highest-effort single item and should
  be resourced accordingly; nothing atlas-wide should proceed until this is
  resolved (or explicitly deferred with the affected cells flagged).
- **Phase D — Regenerate validated study objects.** Re-run each of the six
  per-study preprocessing scripts with assertions for cell counts, keys,
  metadata vector lengths, and count integrity (AT02–AT05). Freeze cell
  membership and counts before Phase E.
- **Phase E — Rebuild merged atlas/integration.** Fix D050/D051 (Harmony
  grouping) and re-run merge → normalization → HVF → scaling → PCA → Harmony
  → neighbors → clustering exactly once, only after Phase D's populations
  are frozen. Evaluate DEC11 (global vs. stratified) using AT24–AT27 before
  finalizing.
- **Phase F — Recompute shared annotations.** Fix D047 (graph.name), D048
  (cluster hard-coding redesign per DEC12), D049 (removal policy per DEC10),
  D058/D059/D060/D061 (HPI naming/gating/uncertainty), then re-run SingleR
  and gametocyte inference once, using the corrected applicability gates and
  provenance model from `annotation_precedence.tsv`.
- **Phase G — Regenerate exports.** Fix D052–D057 in `flatten_data.R`/
  `data_source_manipulation.R`, with explicit key validation (AT16,AT17) and
  the full corrected `metadata_schema.tsv` field set.
- **Phase H — Update application consumer.** Apply the `app_data_contract.tsv`
  changes (D062–D069) against the newly versioned data contract; add the
  config-key renames so `development_phase`/`parasite_stages` naming can
  never again silently reference a non-existent column.
- **Phase I — Final validation and publication freeze.** Run the full
  `acceptance_tests.tsv` suite; produce the build manifest (AT23: input
  checksums, git commit, package versions, seeds, output checksums) that has
  never existed for this pipeline before.

## Pointer from the existing overall audit documentation

A one-line cross-reference should be added to `pre_process_data/audit/README.md`
(or equivalent top-level audit index, if the user maintains one) pointing to
this directory as the consolidation/repair-planning phase that follows the
six completed study audits. This document does not duplicate that content
here to avoid maintaining two copies.
