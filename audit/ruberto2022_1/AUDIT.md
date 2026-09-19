# Ruberto2022_1 forensic audit

Audit date: 2026-09-19. Status: inspection complete; repairs proposed, not
applied. Scope is Ruberto et al. 2022 (liver-stage study, PMID 36093191) only.
This audit continues from, and does not repeat, the completed Mancio-Silva2022,
Sa2020 and Hazzard2024 audits recorded under
`pre_process_data/audit/mancio_silva2022/AUDIT.md` and
`pre_process_data/audit/hazzard2024/AUDIT.md`. Its starting point is the
existing preliminary "Ruberto2022_1: liver source selection" finding in the
Mancio-Silva2022 cross-study section, which classified this study
**APPARENTLY CONSISTENT** for cell-lineage/label propagation and **REQUIRES
VALIDATION** for object provenance. This audit confirms lineage, resolves most
of the provenance question, and **overturns the preliminary "APPARENTLY
CONSISTENT" classification with a new confirmed expression-data defect** that
the preliminary pass did not check.

No production R script, source count matrix, production RDS, integrated
object, metadata, or application code was modified. SingleR/Harmony were not
rerun. All new files are audit documentation, authoritative-source snapshots,
or deterministic audit-only evidence under this directory. R diagnostics used
`env -u R_LIBS_USER pixi run Rscript --vanilla` from `pre_process_data/scripts`.
`rtracklayer`/`GenomeInfoDb` fail to load in this Pixi environment (missing
`GenomeInfoDbData`); the required rRNA feature-ID list was instead extracted
with `awk` directly from the same GFF file used by the production script
(`rRNA_ids.txt`, 50 IDs, `type == "rRNA"`, matching the production filter
exactly). This is an environment limitation, not a change to any script logic.

## Executive finding

Cell inventory and identifier lineage are exact and fully resolved: the
authors' own local object, the PlaViSca study object, and the deployed
`cleaned_dataset.rds` all agree on precisely **1,438 cells**, with zero missing
or extra IDs in any direction. Day, replicate, treatment and run assignment are
now independently confirmed correct against ENA sample metadata and the
published Methods text (previously only "consistent," not source-verified).
Final `life_cycle_stage` labels agree with source `liver_form` for all 1,438
cells (1,147 Hypnozoite, 291 Schizont (Liver stage)), confirming the
preliminary label-propagation finding.

However, **direct inspection of the retained RNA count matrix, not performed
in the preliminary pass, finds that 538 of the 1,438 retained cells (37.4%)
carry an essentially empty P. vivax expression profile (517 with exactly zero
total UMI, 21 with exactly one) in `ruberto2022_1.rds`, and this is unchanged
in the deployed `cleaned_dataset$mr_data`/`raw_df.rds`.** These 538 cells are
concentrated almost entirely in one biological replicate: 89/90 cells of run
SRR19573609 and 449/457 cells of run SRR19573610 (both replicate 2), versus 0
of 891 replicate-1 cells (runs SRR19573611/612). The same 538 cells have
substantial real signal in the authors' own local object (median 718.5, mean
1,523.6, max 18,273 `nCount_RNA`), proving they are genuine non-empty cells in
the original data and that the emptiness is a defect of PlaViSca's own
count-matrix generation for these two runs, not of the underlying biology or
of the day/treatment/run metadata join (independently verified correct). All
1,438 cells, including these 538, still receive full downstream biological
annotations in the deployed app: a `pred_gametocyte` sex/stage call (20 of
them "Female gametocyte") and a `hour_post_invasion` label, both computed from
expression that is at or near zero for over a third of the study.

| Confirmed problem | Affected deployed cells | Exact scope |
|---|---:|---|
| Empty/near-empty retained expression | 538 | 517 exactly 0 UMI, 21 exactly 1 UMI; runs SRR19573609 (89/90) and SRR19573610 (449/457), i.e. all of biological_replicate 2 |
| Gametocyte-sex prediction on liver-stage cells | 1,438 (20 labeled Female gametocyte) | Liver hypnozoites/schizonts cannot be gametocytes; predictions include the 538 near-empty cells |
| HPI assigned despite no measured-HPI meaning | 1,438 | seven blood-reference similarity labels on every Ruberto2022_1 cell, corroborating the general cross-study HPI finding for this specific study |
| Misspelled/unverified technology | 1,438 | `10x_Chomium_V3` on every cell; V3 chemistry not independently confirmed in the inspected publication text |
| Deployed treatment vocabulary drift | 1,438 | deployed `None`/`MMV390048` vs current script's `No_Treatment`/`MMV390048`; consistent with the already-documented historical-vs-current-script divergence, not a new script defect |

## Authoritative publication and accession identity

| Concept | Verified value | Authority |
|---|---|---|
| PMID | 36093191 | PubMed / article XML |
| DOI | 10.3389/fcimb.2022.986314 | article XML |
| Title | Single-cell RNA profiling of *Plasmodium vivax*-infected hepatocytes reveals parasite- and host- specific transcriptomic signatures and therapeutic targets | article XML |
| Publication date | 2022-08-25 | Europe PMC record |
| PMC | PMC9453201 | Europe PMC record |
| BioProject | PRJNA843856 | data-availability statement and ENA |
| Deposited runs used by PlaViSca | SRR19573609–SRR19573612 (infected hepatocytes) | ENA run report |
| Deposited runs NOT used | SRR19573613, SRR19573614 (uninfected controls) | ENA run report |
| Author code | GitHub `AnthonyRuberto/Pv_LS_singleCell` | publication code-availability statement |
| Author archived outputs | Zenodo 10.5281/zenodo.6463338 | publication data-availability statement; not bulk-downloaded (large archive) |

`article.xml`, `ena_PRJNA843856_runs.tsv`, and the six author R Markdown
workflow files plus README (`ruberto_github_*`) are retained here. The GitHub
default-branch tree listing is `github_tree.json`.

## Exact day/replicate/treatment/run identity, now source-verified

The preliminary audit inferred day/replicate/run correspondence from the
script's own naming convention. This audit instead verifies it independently
against ENA sample titles and the publication Methods:

| Run | ENA sample alias / title | Script day | Script treatment | Script replicate |
|---|---|---:|---|---:|
| SRR19573609 | Pv_infected_hepatocyte_9dpi_rep2 | 9 | MMV390048 | 2 |
| SRR19573610 | Pv_infected_hepatocyte_5dpi_rep2 | 5 | No_Treatment | 2 |
| SRR19573611 | Pv_infected_hepatocyte_9dpi_rep1 | 9 | MMV390048 | 1 |
| SRR19573612 | Pv_infected_hepatocyte_5dpi_rep1 | 5 | No_Treatment | 1 |

All four agree exactly. The publication Methods independently confirm the
treatment/day logic itself, not only the accession-to-day mapping: "Hepatocytes
were treated with 1 µM MMV390048 on days 5, 6, and 7 post-infection to kill
replicating schizonts, resulting in cultures enriched with hypnozoites at 9
days post-infection" and "Hepatocyte cultures infected with P. vivax were
processed prior to- and post- MMV390048 treatment (5- and 9-days
post-infection, respectively)." This **resolves** the preliminary audit's open
qualification that "drug assignment draws on the experimental protocol rather
than that field alone": the protocol itself is now directly confirmed by the
article text, not only inferred. All four runs are independently confirmed
Illumina HiSeq X Ten by ENA, matching the publication and the script's
`sequencer` field; there is no platform conflict here (unlike Hazzard2022).

The publication's stated total, "we characterized 1,438 parasite transcriptomes
across two independent experimental replicates," matches the retained cell
count exactly.

## Author-object provenance: resolved by the author's own deposited code

The preliminary audit flagged `Hep59.1.2.seu_20aug2025.rds` chain of custody as
**REQUIRES VALIDATION** because the inspected GitHub tree did not obviously
identify that dated filename. This audit fetched and read the six author R
Markdown workflow files directly (`ruberto_github_additional_file_1.Rmd`
through `_6.Rmd`). `additional_file_3.Rmd` constructs and saves an object named
exactly `Hep59.1.2.seu` (`saveRDS(Hep59.1.2.seu, "outputs/Hep59.1.2.seu.rds")`,
line 205), built as `merge(Hep59.1.seu, Hep59.2.seu)`, itself
`merge(Hep51a.seu.viv.f1, Hep91a.seu.viv.f1)` (replicate 1, day 5 + day 9) and
`merge(Hep54a.seu.viv.f1, Hep94a.seu.viv.f1)` (replicate 2, day 5 + day 9).
`additional_file_4.Rmd` reads exactly `outputs/Hep59.1.2.seu.rds` and assigns
`RenameIdents(Hep59.1.2.seu, `1` = "Hypnozoites", `2` = "Schizonts")` on
re-clustered `seurat_clusters` of the full merged 1,438-cell object — the exact
two-level `LiverForm` vocabulary found in the local file. The local filename
differs only by an appended `_20aug2025` date suffix.

This is a strong name/structure/construction-logic match to the authors' own
documented pipeline, not a PlaViSca-invented label. It is **not** a byte-level
provenance proof: the Zenodo archive (10.5281/zenodo.6463338) containing the
canonical `outputs/Hep59.1.2.seu.rds` was not downloaded (large output
package), so an MD5/checksum match to the deposited file is not established.
Local object metadata is otherwise informative and consistent: `day` levels
are `Five`/`Nine`, `replicate` is `1`/`2`, Seurat object version is 4.0.3, and
stored `commands` show `NormalizeData`, `SCTransform`, `RunPCA`, `RunUMAP`,
`RunTSNE`, `FindNeighbors`, `FindClusters` — a complete, plausible single
build history for a merged, re-clustered object, with no signs of manual
metadata patching. `additional_file_2.Rmd` independently corroborates the
day-9 "solely hypnozoites" biology: only day-9 objects are saved directly as
the hypnozoite-only ("small") subset with no schizont ("large") counterpart
extracted, while day-5 objects are split into both `.large` (Schizonts) and
`.small` (Hypnozoites) subsets by the same `1`/`2` cluster-rename rule. This
is corroborated locally: **zero of the 345 day-9 cells carry a Schizont
LiverForm label** (day x LiverForm cross-tab: Nine → 345 Hypnozoites, 0
Schizonts), exactly matching this documented construction.

## Source-to-final cell lineage: exact, three-way

| Object | Cells | Unique IDs |
|---|---:|---:|
| `Hep59.1.2.seu_20aug2025.rds` (author-named local object) | 1,438 | 1,438 |
| `ruberto2022_1.rds` (PlaViSca study object) | 1,438 | 1,438 |
| `cleaned_dataset$mr_data`, Ruberto2022_1 rows (deployed app) | 1,438 | 1,438 |

The script's `case_when` conversion (`51inf`→`612`, `52inf`→`610`, `91inf`→
`611`, `92inf`→`609`) maps every one of the 1,438 Hep59 cell IDs to a unique
`<SRR-suffix>_<barcode>` study identifier with **zero unmapped, duplicated, or
excess IDs** in either direction (`hep59_to_study_cellid_map.tsv`). The
resulting `liver_form` join into the study object has **zero missing values**
(`study_object_metadata.tsv`). The deployed `cleaned_dataset.rds` rows for
`study_label == "Ruberto2022_1"` are exactly the same 1,438 identifiers as the
study object, with zero missing on either side
(`deployed_ruberto1_metadata.tsv`). `life_cycle_stage` in the deployed data
agrees with study-object `liver_form` for all 1,438 cells with **zero
disagreements** (1,147 Hypnozoite / 291 Schizont (Liver stage) both sides).
This directly confirms and sharpens the preliminary "final liver forms and
final lifecycle labels agree for every cell" finding with an explicit
cell-by-cell cross-tabulation rather than aggregate counts alone.

Note the app's deployed `parasite_stage` field is a distinct, coarser field
(`"Liver stage"` for all 1,438 cells, a broad-phase label), and `treatment` is
deployed as `None`/`MMV390048` rather than the current script's
`No_Treatment`/`MMV390048` — consistent with the general historical-vs-current
production-script divergence already documented for other studies; preserve
`No_Treatment` going forward and do not read the deployed `None` string as
evidence of a distinct code path.

## Confirmed defect: near-total expression loss in one full replicate

Gene-level count preservation was checked between `ruberto2022_1.rds` and
PlaViSca's own STARsolo `GeneFull/raw` matrices for the same four runs
(`inspect_ruberto2022_1.R`, section 6). After removing the same 50 GFF-annotated
rRNA features the production script removes, all 6,811 study-object features
matched raw features exactly, and there were **zero unequal count entries**
across all 1,438 retained cells and all four runs. `ruberto2022_1.rds` is
therefore an exact, unmodified copy of PlaViSca's own realigned raw counts for
its retained barcodes; this rules out an import/subsetting bug as the cause of
what follows. Unlike Mancio-Silva2022 (where raw counts are read directly from
the authors' own deposited processed object), Ruberto2022_1's counts are
independently regenerated by PlaViSca's own STARsolo alignment of the raw
FASTQ files, so this check establishes internal reproducibility of the
production pipeline, not preservation of the authors' original processed
values; the authors' own count matrix was not deposited as a directly
comparable raw file in the inspected evidence.

That same count matrix, however, is empty or near-empty for a very specific
and large subset of cells:

| Metric | Value |
|---|---:|
| Cells with exactly 0 total UMI | 517 |
| Cells with exactly 1 total UMI | 21 |
| Cells with ≤1 total UMI | 538 (37.4% of 1,438) |
| Median total UMI, all 1,438 cells | 117 |
| Mean total UMI, all 1,438 cells | 801.8 |

These 538 cells are not spread evenly across the study: **89 of 90 SRR19573609
cells and 449 of 457 SRR19573610 cells** are near-empty, versus **0 of 255
SRR19573611 cells and 0 of 636 SRR19573612 cells**. SRR19573609 and
SRR19573610 are exactly the two runs belonging to biological replicate 2 (day
9 and day 5 respectively); SRR19573611/612 are replicate 1. In other words,
**98.4% of all replicate-2 cells (538 of 547) are near-empty in PlaViSca's own
count matrix, and 0% of replicate-1 cells are.**

The same 538 cells have clearly real signal in the authors' own local object
for the identical cell IDs (`near_empty_cells_vs_hep59.tsv`): `nCount_RNA`
median 718.5, mean 1,523.6, max 18,273; `nCount_SCT` median 719, mean 820.3.
Non-near-empty study cells have author-object `nCount_RNA` median 983, mean
4,117 — of the same general order, confirming these are ordinary retained
cells in the authors' data, not cells the authors themselves flagged as
low-quality. **This demonstrates the emptiness is specific to PlaViSca's own
replicate-2 count regeneration, not to the cells' underlying biology, and not
to the day/treatment/replicate/run metadata assignment**, which was
independently verified correct against ENA and the publication above.

STARsolo's own per-run `Summary.csv` is consistent with a technical fault
localized to these two FASTQ inputs/alignments, though the audit does not
establish its exact mechanism:

| Run | Estimated Number of Cells | Median UMI per Cell (STARsolo) | Reads Mapped to Genome (Unique+Multi) |
|---|---:|---:|---:|
| SRR19573609 (rep2, d9) | 61,444 | 84 | 6.69% |
| SRR19573610 (rep2, d5) | 42,071 | 73 | 6.23% |
| SRR19573611 (rep1, d9) | 276 | 179 | 0.74% |
| SRR19573612 (rep1, d5) | 156 | 4,706 | 2.54% |

The rep2 runs report tens of thousands of "cells" with very low per-cell UMI —
consistent with a diffuse, ambient-like signal spread across a much larger
fraction of the barcode space than a true single-cell library, rather than the
concentrated signal STARsolo finds for rep1. This pattern (a higher overall
genome-mapping fraction together with an implausibly large "cell" count and
very low per-cell depth) is **consistent with, but does not prove**, a
technical fault specific to the rep2 FASTQ files or their alignment/barcode
extraction; no FASTQ-level or alignment-log-level root-cause investigation was
performed in this audit, and no claim is made about which specific step
introduced it. What is proven is the downstream consequence: PlaViSca's own
regenerated counts do not recover the authors' real replicate-2 parasite
signal for the retained cell barcodes.

Among the 900 non-near-empty cells, three small exact-duplicate profile groups
were found by exhaustive SHA-256 grouping of all 1,438 full 6,811-feature
count columns (`exact_duplicate_profiles_ruberto1.tsv`): two groups of 10
cells each sharing a single incidental 1-UMI count, and the 517-cell all-zero
group itself, which trivially collapses to one hash. **These are degenerate
artifacts of the near-empty-count defect above, not evidence of Mancio-style
biological/technical profile duplication**, and must not be reported as
duplicate cells in the Mancio `.1`-pair sense. Excluding these three
degenerate groups, there are no other exact-duplicate RNA profiles among the
1,438 cells.

Both `raw_df.rds` and `normalize_df.rds` in the deployed
`PlaViSca/data` directory inherit this defect unchanged: the 538 near-empty
study cells have deployed `raw_df` P. vivax gene-column row-sums with median 0,
mean 0.039, and max 1 — i.e. the emptiness reaches the live application
exactly as found in the study object, with no intervening masking, imputation
or exclusion.

## Confirmed defect: biologically inapplicable downstream annotations reach all cells, including the empty ones

The deployed `cleaned_dataset$mr_data` assigns a `pred_gametocyte` SingleR call
to all 1,438 Ruberto2022_1 cells, including 20 labeled `Female gametocyte`, and
the remainder distributed across 19 numbered `Asexual` blood-cycle clusters.
Liver-stage hypnozoites and schizonts are pre-gametocyte forms; gametocyte
commitment occurs in the blood stage, so a liver-stage cell cannot correctly
be a measured gametocyte, and the blood-stage `Asexual N` vocabulary itself
does not apply to liver forms. This is not a new mechanism — it is the same
blood-stage-reference SingleR call already flagged generally in the
Mancio-Silva2022 cross-study section — but this audit newly confirms its exact
per-study scope for Ruberto2022_1: **all 1,438 cells**, not a subset, and
critically, the call is made **without regard to the 538 near-empty cells
above**: a SingleR/cluster call cannot be biologically informative when
computed from an all-zero or near-zero expression vector, so at minimum those
538 calls (including any of the 20 Female gametocyte calls that fall among
them) carry no supporting signal at all.

`hour_post_invasion` is likewise nonmissing for all 1,438 cells (6→32, 12→51,
24→262, 32→128, 36→925, 42→24, 48→16), independently corroborating and
quantifying the general cross-study HPI-applicability finding specifically for
this study. The Zhu blood-stage IDC reference has no established applicability
to liver-stage hypnozoites/schizonts.

## Other metadata checks

- `sc_technology = "10x_Chomium_V3"` is a confirmed misspelling ("Chomium");
  the inspected publication text and GitHub workflow files describe 10x
  Genomics single-cell partitioning but the exact V3 chemistry claim was not
  independently located in the inspected publication passages. Treat V3 as
  unverified, analogous to the Sa2020 `10x_Chomium_V2` finding.
- `geographic_location = "Cambodia_Mondulkiri"` is confirmed by the
  publication's description of blood collection in Mondulkiri province,
  Cambodia, and describes parasite origin specifically; it is not conflated
  with host-hepatocyte provenance in this script (unlike Mancio-Silva2022 and
  Sa2020, Ruberto2022_1 already keeps `geographic_location` and `host_id`
  separate).
- `host_id = "BioIVT:BGW"` is confirmed by the publication ("Lot BGW,
  BioIVT" / "donor BGW, BioIVT"), correctly describing the hepatocyte donor,
  not the parasite.
- `strain = "Cambodia field isolate"` matches the publication's description of
  Cambodian field isolates; no distinct named strain identifier was found.
- `host_species = "Homo sapiens"` correctly describes the hepatocyte host, with
  no parasite/host species conflation.
- No gene-ID collisions: all 6,811 study-object feature IDs are unique after
  Seurat's underscore-to-hyphen conversion.
- `num_srr = 4` matches the count of runs actually imported and retained
  (609–612); the two deposited uninfected-control runs (613, 614) are
  correctly excluded from this study script. As in other studies, this is a
  retained-run count, not a biological-replicate count, and should be moved to
  study-level provenance rather than ordinary per-cell metadata in any future
  repair, consistent with the general policy already proposed for other
  studies.

## What was and was not resolved

| Question | Status |
|---|---|
| Exact 1,438-cell inventory and identifier lineage (author object → study object → deployed app) | **Resolved**: exact three-way match, zero missing/extra either direction |
| Day/replicate/treatment/run assignment | **Resolved**: independently confirmed against ENA and the published Methods text, not only internal script consistency |
| `Hep59.1.2.seu_20aug2025.rds` name/construction-logic provenance | **Resolved** against the authors' own GitHub workflow code (object name, merge structure, and cluster-rename logic match exactly) |
| `Hep59.1.2.seu_20aug2025.rds` byte-level identity to the Zenodo deposit | **Not resolved**: the Zenodo archive was not downloaded; no checksum comparison was made |
| RNA count concordance, `ruberto2022_1.rds` vs PlaViSca's own STARsolo output | **Resolved**: exact, zero unequal entries |
| RNA count concordance vs the authors' own original processed counts | **Not established**: no directly comparable author raw/processed count file was located in the inspected evidence; PlaViSca's counts are an independent realignment |
| Cause of near-total replicate-2 signal loss | **Not resolved**: the consequence (538 near-empty retained cells, concentrated exactly in replicate 2) is proven; the technical root cause in FASTQ/alignment is not established by this audit |
| Biological validity of gametocyte-sex and HPI calls on liver-stage cells | **Not resolved as biologically valid**; confirmed as methodologically inapplicable, consistent with the general cross-study finding, now quantified per-study |
| Exact-duplicate biological profiles (Mancio-`.1`-style) | **Resolved**: none found beyond the degenerate near-empty-count artifacts, which are explicitly not biological duplicates |

## Proposed repair policy (not applied)

1. Preserve `source_cell_id` (the Hep59 identifiers) and the validated
   `51inf/52inf/91inf/92inf → 612/610/611/609` crosswalk as an explicit,
   reviewed relation, exactly as already recommended in general for other
   studies' source-group-to-accession maps. Fail on any future unmapped or
   duplicated ID.
2. Before any further use of Ruberto2022_1 expression data (differential
   expression, integration, marker scoring), explicitly flag and consider
   excluding the 538 near-empty cells identified in
   `near_empty_cells_vs_hep59.tsv`, or regenerate their counts from a
   corrected realignment of SRR19573609/610. Do not silently drop them without
   review: they are genuine author-retained cells with real biology in the
   original data, not authors' own low-quality exclusions.
3. Investigate the SRR19573609/610 alignment at the FASTQ/STARsolo-log level
   before deciding on a repair (reprocess vs. exclude vs. flag-and-retain).
   This audit establishes the defect and its scope but not its technical cause.
4. Do not compute or display `pred_gametocyte`/`hour_post_invasion` for
   Ruberto2022_1 (or other confirmed non-blood-stage studies) without an
   explicit sample_type/development_phase gate, consistent with the general
   cross-study HPI/gametocyte recommendation; additionally gate on nonzero
   expression so that near-empty cells do not receive any stage/sex call.
   Preserve source `liver_form` and any future replacement prediction as
   clearly separate, provenance-labeled fields, as already established for
   Mancio-Silva2022's `liver_form`/SingleR separation.
5. Correct `10x_Chomium_V3` to `10x_Chromium`; do not assert V3 chemistry
   without further evidence, consistent with the Sa2020 precedent.
6. Preserve the deployed `treatment` vocabulary difference (`None` vs.
   `No_Treatment`) as a known historical-build artifact; standardize on
   `No_Treatment` in any future regeneration, never silently reintroducing
   `None`.
7. Move `num_srr` out of ordinary per-cell metadata into study-level
   provenance, consistent with the general policy already proposed for other
   studies.

**Final Ruberto2022_1 classification:** confirmed near-total expression loss
for biological replicate 2 (538/1,438 cells, previously undetected);
confirmed inapplicable gametocyte-sex and HPI predictions reaching all 1,438
cells, including the near-empty ones; resolved cell/day/replicate/treatment/run
lineage and author-object construction-logic provenance; resolved absence of
Mancio-style biological duplicate profiles; unresolved technical root cause of
the replicate-2 signal loss and byte-level Zenodo provenance of the local
author object. The preliminary "APPARENTLY CONSISTENT" classification from the
Mancio-Silva2022 cross-study section is superseded for expression data by this
audit's confirmed defect; the preliminary label-propagation and accession
findings are confirmed and sharpened, not overturned.

## Reproducibility

```text
env -u R_LIBS_USER pixi run Rscript --vanilla ../audit/ruberto2022_1/inspect_ruberto2022_1.R
env -u R_LIBS_USER pixi run Rscript --vanilla ../audit/ruberto2022_1/inspect_deployed.R
```

Both scripts read only existing production/app files and write only into this
audit directory. Production input checksums are in `production_inputs_md5.tsv`.
Run logs are saved as `inspect_ruberto2022_1.log` and `inspect_deployed.log`.
No production script, source data, study object, integrated object or deployed
app artifact was changed. No cells were removed and no correction was applied;
validation after correction is not applicable.
