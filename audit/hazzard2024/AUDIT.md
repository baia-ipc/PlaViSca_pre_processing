# Hazzard2024 forensic audit

Audit date: 2026-09-19. Status: inspection complete; repairs proposed, not
applied. Scope is Hazzard et al. 2024 / PMID 39223117 only, except for the shared
downstream code needed to trace Hazzard lifecycle and HPI annotations.

No production R script, source count matrix, production RDS, integrated object,
metadata, or application code was modified. Harmony and SingleR were not rerun.
All new files are audit documentation, authoritative-source snapshots, or
deterministic audit-only evidence under this directory. R diagnostics used
`env -u R_LIBS_USER pixi run Rscript --vanilla` from `pre_process_data/scripts`.

## Executive finding

The deployed Hazzard2024 inventory is **exactly 80,024 cells**, not merely
approximately 80,024. The publication's 24 blood-library counts, the 80,024-row
Zenodo processed table, the Hazzard study RDS, and the deployed app rows agree
exactly in inventory. All source identifiers are unique and transform one-to-one
to unique PlaViSca identifiers. No exact RNA count-profile duplicates were found
within the study RDS.

That strong inventory result does **not** validate the attached metadata. The
following errors are confirmed in both `hazzard2024.rds` and the deployed
`cleaned_dataset.rds`:

| Confirmed problem | Affected deployed cells | Exact scope |
|---|---:|---|
| Wrong library/host identifier | 19,810 | SRR27021986: 3,489; SRR27021993: 16,321 |
| Wrong parasite strain | 16,321 | all SRR27021993 cells are deposited/published Chesson but labeled NIH1993-F3 |
| Wrong day post infection | 9,644 | SRR27021970: 3,177; SRR27021971: 3,517; SRR27021972: 2,950 |
| Misspelled/underspecified technology | 80,024 | `10x_Chomium_V3` on every cell |
| Missing biological-replicate field | 80,024 | literal `"NA"`; the ten animal identities are known |
| HPI assigned despite no measured-HPI meaning | 80,024 | seven blood-reference similarity labels on every Hazzard cell |
| HPI assigned to author `Sexual` group | 21,009 | especially inappropriate for asexual blood-cycle time interpretation |

The source/deployed stage comparison also shows a substantive authority conflict.
The deposited author table has 21,009 `PseudoGroup=Sexual` cells. PlaViSca calls
20,292 cells `Female gametocyte`, but these are not the same set: 1,552 author
Sexual cells receive an asexual final label, while 835 author GroupA/B/C cells
become Female gametocyte. There is no final male-gametocyte call. The final
20,292 female labels come from a hard-coded integrated-cluster interpretation,
not a keyed import of the author `PseudoGroup` field.

## Authoritative publication and accession identity

| Concept | Verified value | Authority |
|---|---|---|
| PMID | 39223117 | PubMed / article XML |
| DOI | 10.1038/s41467-024-51949-8 | article XML |
| Title | Single-cell analyses of polyclonal *Plasmodium vivax* infections and their consequences on parasite transmission | article XML |
| Publication date | 2024-09-02 | article XML |
| PMC | PMC11369214 | article XML |
| BioProject | PRJNA1047651 | publication data-availability statement and ENA |
| ENA/SRA study | SRP475376 | ENA run report |
| GEO | none identified or cited for this study | publication/accession records |
| Sequencing runs deposited | 43 unique SRR accessions | ENA snapshot |
| Deposited processed data | Zenodo 12775216 | publication data-availability statement |
| Deposited analysis code | GitHub `bhazzard11/Coinfection_Analysis`, archived as Zenodo 12775129 | publication code-availability statement |

The article XML, ENA project/run metadata, NCBI BioSample XML, publication
Supplemental Data 1, author-code release, and checksummed processed table are
retained here. `Proccessed_Data.txt` has MD5
`6ede1f7da4e5d57841c9f90883b2b2c1`, exactly matching the Zenodo-advertised
checksum. Supplemental Data 1 has MD5
`702418299ebb4d4e3af441fec9966abc` and was downloaded from the publisher.

## Deposited-run membership

ENA reports 43 runs, all registered as Illumina NovaSeq 6000. Their accession
aliases divide into 24 blood scRNA-seq runs, six salivary-gland sporozoite runs,
eleven `AS` auxiliary runs, and two WGS runs. The publication describes 24 blood
libraries, eight sporozoite libraries, and two WGS libraries. The relation of the
eight publication sporozoite libraries to the six deposited SPZ run accessions is
not fully resolved by the run aliases alone and is irrelevant to the deployed
blood-cell population, but it must not be silently called one-to-one.

The six PlaViSca chunk scripts import 30 distinct runs: all 24 retained blood
runs plus the six SPZ runs. The deposited processed-cell-list subset retains only
the 24 blood runs. The eleven auxiliary and two WGS runs are not imported. The
six imported SPZ runs contribute zero deployed Hazzard2024 cells. Exact membership,
sample/experiment accessions, aliases and retained counts are in
`all_deposited_run_membership.tsv`.

`num_srr=24` on every deployed cell is therefore a retained blood-library count.
It is neither the deposited run count (43), the imported run count (30), nor a
biological-replicate count. It should be removed from ordinary cell metadata and
replaced with explicitly scoped study/run provenance.

## Study design and authoritative retained-library structure

The publication used ten female *Saimiri boliviensis* monkeys. It generated
10X Genomics 3-prime scRNA-seq libraries (catalog 1000268), loading an estimated
5,000 enriched infected cells per sample, and sequenced the libraries on an
Illumina NovaSeq. Blood-stage cells were retained by the authors at least
1,000 unique reads mapped within annotated genes. This produced:

| Infection design | Animals | Libraries per animal | Published retained cells |
|---|---:|---:|---:|
| Chesson mono-infection | 3 | 1–2 | 22,847 |
| NIH-1993-F3 mono-infection | 3 | 1–2 | 16,292 |
| Consecutive infection | 2 | 3 | 18,286 |
| Simultaneous infection | 2 | 2 | 11,715 |
| Infection initiated by sporozoites | 2 | 3 | 10,884 |
| **Total retained blood cells** | **10 distinct animals** | **24 libraries** | **80,024** |

Animal is the biological-replicate unit; serial libraries from one animal are
repeated measures, not independent biological replicates. The complete run-level
source-of-truth table is `authoritative_run_table.tsv`. It records exact library
ID, animal ID, day, infection design, inoculum/lineage description, publication
read and cell counts, SRR/SRX/SAMN/SRP accessions, registered platform, deposited
processed count, deployed values, and provenance notes.

The processed author table contains per-cell `Genotype`, `Condition`, allele
counts, pseudotime and `PseudoGroup`. Its genotype inventory is Chesson 36,844,
NIH 38,650, Outcross 4,529 and one missing value. Its condition inventory is
mono_CH 22,847; mono_NIH 16,292; Coinfection 11,715; NIH_CH 10,113; CH_NIH 8,173;
and Outcross 10,884. These cell genotypes must remain separate from library
inoculum composition and infection design.

## Provenance concepts that must remain separate

The existing `geographic_location` field conflates several concepts and is not a
valid host or experimental-location field.

- **Parasite provenance.** Chesson is a historical Papua New Guinea lineage.
  NIH-1993-F3 is described by the cited lineage work as closely related to, but
  distinct from, Salvador-I. The exact geographic origin of NIH-1993 itself is
  not established by this publication. Salvador-I's El Salvador origin must not
  automatically be transferred to NIH-1993-F3.
- **Host provenance.** The publication identifies female *Saimiri boliviensis*
  obtained from NIH-approved sources. It does not report donor/birth geography.
- **Experimental site.** Animal work is under NIAID/NIH protocols. An institutional
  location must not be interpreted as animal or parasite origin.
- **Sequencing site.** The publication and Supplemental Data 1 specify the NovaSeq
  platform but do not state a physical sequencing site.
- **Registry location.** NCBI BioSample records use `USA: Bethesda`, `Homo sapiens`
  and generic titles such as “Invertebrate sample from Plasmodium vivax” even for
  these monkey-blood libraries. These registry fields conflict with the publication
  biology and must be retained as registry evidence, not used uncritically as
  canonical host or origin metadata.

Current deployed geography counts are Papua New Guinea 22,847; El Salvador 8,256;
and Papua New Guinea/El Salvador 48,921. The last category incorrectly includes
8,036 NIH-only mono-infection cells merely because the chunk script assigns one
value to all eight runs. All 80,024 geography values are strain-origin proxies in
a generically named field; they do not describe host provenance, experimental
site or sequencing site.

## Cell inclusion, identifiers and count lineage

`hazzard2024_merge_all.R` reads the deposited `Proccessed_Data.txt`, converts each
author ID from `animal_timepoint_barcode` to `runSuffix_barcode` using a 24-entry
handwritten crosswalk, and subsets the merged candidate object to those IDs.
The audit reproduced this transformation:

- 80,024 source rows and 80,024 distinct source IDs;
- 80,024 distinct reconstructed PlaViSca IDs;
- zero reconstructed IDs missing from `hazzard2024.rds`;
- zero extra study-RDS cells;
- zero study-RDS cells absent from the deployed Hazzard rows;
- zero deployed Hazzard rows absent from the study RDS;
- 6,583 RNA features and 47,311,903 nonzero study-RDS counts;
- zero exact duplicate RNA count profiles among the 80,024 columns.

`cell_lineage.tsv` gives the complete source-to-study-to-deployed record lineage.
The five available chunk RDS files contain millions of pre-subset barcode
candidates; those are not evidence of millions of additional biological cells.
The missing 59–65 chunk RDS prevents repeating that preliminary candidate-union
check for run 960, but the final inclusion identity is complete.

The author count-table ZIP was advertised with MD5
`276aa89b96636e9091a05cd43166465e` but repeated Zenodo content requests returned
HTTP 504 during this audit. Consequently, the study-RDS counts have not been
proved equal to the author's deposited count-table values. The available evidence
does prove processed-list inclusion and study-to-deployed identity; it must not be
overstated as raw/source-count correctness.

## Exact metadata defects and mechanism

### Host/library identifiers

The deployed `host_id` actually holds library/timepoint identifiers for nearly
all rows, not the animal identifier. Two handwritten values are additionally
wrong:

| Run | Published/source library | Deployed value | Cells |
|---|---|---|---:|
| SRR27021986 | 5537_2 | 5537_3 | 3,489 |
| SRR27021993 | 5708_1 | 5709_1 | 16,321 |

Thus 19,810 cells have a wrong library identifier. A corrected schema needs both
`animal_id` (e.g. 5537) and `library_id` (e.g. 5537_2); neither should be called
the other.

### Parasite strain/genotype

All 16,321 SRR27021993 cells are a Chesson mono-infection in publication
Supplemental Data 1. All 16,321 deposited processed rows have `Genotype=Chesson`
and `Condition=mono_CH`. PlaViSca labels every one `NIH1993-F3`. The four-run
91–94 script has five strain values and uses the fourth (`NIH1993-F3`) for run
993, leaving the fifth (`Chesson`) unused. This is a direct positional construction
error, not an interpretive disagreement.

For mixed libraries, the deployed coarse strain strings describe the inoculum or
infection history, while deposited per-cell genotype varies. Both are useful but
must be separate fields. For example, run 972 contains NIH and Outcross cells;
one scalar library-composition string cannot replace cell genotype.

### Days post infection

The seven-run 66–72 script supplies eight days and loops over only seven entries.
The surplus final value is unused, but three used positions are also wrong:

| Run / library | Deployed day | Published day | Cells |
|---|---:|---:|---:|
| SRR27021970 / 5350_2 | 31 | 30 | 3,177 |
| SRR27021971 / 5708_3 | 30 | 24 | 3,517 |
| SRR27021972 / 5350_1 | 24 | 29 | 2,950 |

This affects 9,644 cells. The unused eighth value does not cause an R length error,
but it concealed a shifted/manual mapping defect.

### Other metadata defects and ambiguities

- `run_id <- rep(run_vector, n)` creates n-squared values in every chunk. Only
  the first n are indexed, so it does not recycle run IDs across cells, but it is
  an unchecked construction defect.
- `10x_Chomium_V3` is misspelled for all 80,024 cells and overstates chemistry
  precision relative to the safest publication-supported label. Preserve the
  10X 3-prime protocol and catalog number; normalize chemistry only with explicit
  catalog authority.
- `biological_replicate` is the literal string `"NA"` for all 80,024 cells even
  though animal IDs and repeated-measure structure are known.
- `No_Treatment` is not an adequate design label. There is no drug-treatment arm,
  but inoculum design, consecutive/simultaneous infection, route, timepoint and
  splenectomy are separate experimental covariates.
- The study scripts call the material `Mammalian host: blood`; deployed export
  changes all 80,024 to `Host blood`. This is a vocabulary transformation, not a
  demonstrated biological error. Treatment similarly becomes `NA` downstream.
- `num_srr=24` is ambiguous as described above.

## Shared HPI and lifecycle impact, restricted to Hazzard2024

All 80,024 deployed Hazzard cells have one of seven `hour_post_invasion` values:
6 (3,313), 12 (4,686), 24 (15,723), 32 (7,055), 36 (33,068), 42 (10,215), and
48 (5,964). These are SingleR similarity labels against the Zhu SMRU1 asexual
blood-stage reference. They are not measured time since invasion and no inference
confidence, score, delta or pruning information is persisted in the deployed
metadata.

The deployed lifecycle inventory is Female gametocyte 20,292; Ring 6,732;
Trophozoite 14,766; Schizont (Blood stage) 35,984; and Merozoite 2,250. Every one
of the 20,292 final Female gametocytes also retains an HPI label. The author source
has `PseudoGroup` totals GroupA 39,444, GroupB 12,422, GroupC 7,149 and Sexual
21,009. Cross-tabulation is in `source_vs_deployed_lifecycle.tsv`.

The final female label is a hard-coded interpretation of integrated cluster 2.
The scripts calculate female/male marker module scores but do not use an explicit
score threshold to assign sex; they hard-code cluster 2 as female and cluster 4
as male. Hazzard contributes 20,292 cells to cluster 2 and zero to cluster 4.
The `FindClusters` call uses `graph_name`, not Seurat's `graph.name`, so it does
not select the claimed `new_clustering` graph; no such graph is created in the
repository. Cluster numbering and the 20,292-cell sex assignment are therefore
not independently validated.

The deposited author `PseudoGroup` is itself a computational annotation based on
author PCA/UMAP and marker interpretation, not direct morphology. It should be
preserved with its authority, alongside PlaViSca inference fields, rather than
silently overwritten. Future HPI display should say reference similarity and be
restricted to a validated applicable asexual population. Author Sexual cells and
PlaViSca gametocyte calls must not carry an asexual HPI as though it were measured
infection age.

Hazzard2024 comprises 80,024 of the deployed atlas's 105,963 cells (75.52%). It
therefore dominates the shared Harmony geometry, neighbor graph, clustering and
any cluster-based annotation. This is an impact statement, not evidence that
Harmony is invalid. The central integrated object is absent and Harmony was not
rerun, so embeddings/clusters cannot be reproduced or repaired within this audit.

## Proposed repair design (not applied)

1. Replace all parallel handwritten metadata vectors with a reviewed 24-row
   library table keyed by full `run_id`; assert one-to-one run mapping, expected
   row count, no unknown runs and publication cell totals before construction.
2. Preserve `source_cell_id`, `source_group`, processed barcode token, source
   `Genotype`, `Condition`, allele counts, pseudotime and `PseudoGroup` from the
   deposited table by exact cell-key join.
3. Store `animal_id` and `library_id` separately. Set animal as the biological
   replicate and model its multiple libraries as repeated measures.
4. Correct run 986 to library 5537_2 and run 993 to library 5708_1. Correct run
   993's library inoculum and every deposited cell genotype to Chesson without
   collapsing the two concepts.
5. Correct days for runs 970, 971 and 972 to 30, 24 and 29. Delete surplus vector
   entries and fail on any vector/table cardinality mismatch.
6. Replace generic geography with separate fields for parasite lineage provenance,
   host provenance, experimental site, sequencing site and registry location,
   each carrying an authority/status note. Do not assign Salvador-I's origin to
   NIH-1993-F3 without direct evidence.
7. Use explicit study-level counts: deposited runs 43, imported runs 30, retained
   blood libraries 24. Do not store an ambiguous `num_srr` on every cell.
8. Preserve author stages and PlaViSca predictions in distinct fields. Do not
   expose SingleR reference labels as measured HPI. Revalidate Hazzard-only
   lifecycle/sex mapping before applying any cluster overwrite.
9. If only metadata is corrected with identical cells/counts, regenerate the
   Hazzard study object and downstream exports through the normal reviewed
   pipeline. If cell inclusion or counts change, rerun the affected integration
   and downstream analyses. Never patch production RDS files in place.
10. Add a build manifest containing source checksums, script commit, package
    versions, run table version, cell-inclusion manifest, random seeds and output
    checksums.

## Reproducibility and evidence index

Run the completed diagnostics from `pre_process_data/scripts`:

```text
env -u R_LIBS_USER pixi run Rscript --vanilla ../audit/hazzard2024/inspect_hazzard2024.R
env -u R_LIBS_USER pixi run Rscript --vanilla ../audit/hazzard2024/inspect_expression_lineage.R
```

Key evidence files:

- `authoritative_run_table.tsv`: primary run-level source of truth.
- `all_deposited_run_membership.tsv`: all 43 runs and import/retention status.
- `cell_lineage.tsv`: all 80,024 source-to-study-to-deployed records.
- `metadata_discrepancies.tsv`: exact issue/run/value/cell-count evidence.
- `source_genotype_condition_counts.tsv`: deposited per-cell genotype/condition.
- `source_vs_deployed_lifecycle.tsv`: author pseudotime group versus PlaViSca HPI,
  lifecycle and cluster interpretation.
- `hpi_lifecycle_impact.tsv`: exact Hazzard-only downstream annotation counts.
- `exact_count_duplicates.tsv`: exhaustive within-study exact-profile result.
- `study_to_deployed_metadata_concordance.tsv`: fieldwise propagation evidence.
- `expression_concordance.tsv`: study-to-deployed raw/normalized comparison.
- `inventory_summary.tsv`, `validation.txt`, `expression_validation.txt`: compact
  machine-readable conclusions and explicit limits.

Outstanding limits are the unavailable author count-table download, the absence
of the central integrated Seurat object/build manifest, unresolved one-to-one
mapping of eight publication sporozoite libraries to six SPZ runs, physical
sequencing site, NIH-1993-F3's exact geographic origin, and independent biological
validation of source and PlaViSca stage/sex annotations. None prevents the exact
run-keyed correction of the confirmed host/library, strain and day defects.
