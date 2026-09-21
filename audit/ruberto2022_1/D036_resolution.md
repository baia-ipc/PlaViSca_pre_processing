# D036 resolution: Ruberto2022_1 replicate-2 near-total expression loss

Status: **B - public-deposit/provenance limitation, confirmed not fixable by
reprocessing.** Investigated to a definitive conclusion in this task; NOT
silently repaired; the 538 affected cells remain flagged and retained.

Date: 2026-09-20.

## Summary of what changed since the Phase 1 blocker report

Phase 1 (`D036_blocker_report.md`) ruled out truncated download, chemistry/
parameter mismatch, and threshold artifacts, but could not distinguish
between "genuinely low parasite signal in the public FASTQ" and "PlaViSca's
alignment/reference configuration fails to recover signal that is present" -
it recommended a controlled kallisto\|bustools re-alignment (matching the
authors' own successful method) as the next decisive step, contingent on a
working Bioconductor environment and either cached or freshly obtained raw
FASTQ.

This task:

1. **Fixed the pixi/R environment** (Part 0 below) - `DropletUtils`,
   `rtracklayer`, `GenomeInfoDb`/`GenomeInfoDbData` all load reproducibly now.
2. **Located the exact raw FASTQ files already on local storage**
   (`/pool/baia/refdata/PlaViSca/fq_data/36093191/`) - the same files
   STARsolo originally aligned. No new download was needed; content was
   independently re-verified against ENA (exact read-count match via full
   decompression for all 4 runs - see `D036_fastq_provenance.tsv`).
3. **Recovered the authors' exact alignment/counting protocol** from their
   Zenodo deposit (10.5281/zenodo.6463338): `kallisto bus -x 10xv3` against a
   **combined human GRCh38 + P. vivax PlasmoDB-51 kallisto index** (k=31),
   followed by `bustools correct|sort|count` - fundamentally different from
   PlaViSca's STARsolo full-genome alignment against a **parasite-only**
   PlasmoDB-68 reference (see `D036_reference_comparison.tsv`).
4. **Ran a controlled kallisto\|bustools re-alignment** of all four local,
   content-verified FASTQ pairs (609, 610 affected; 611, 612 unaffected
   controls) against a parasite-only PlasmoDB-51 transcript index (the
   parasite side of the authors' combined index, obtained from the same
   Zenodo deposit), using the authors' own 10x v3 whitelist and the standard
   `bustools correct|sort|count` workflow (`D036_candidate_processing.tsv`).
5. **Compared the result at the exact 538 previously-affected cell-barcode
   level**, and against the 900 unaffected cells in 611/612 as a positive
   control for the method itself (`D036_cell_level_rescue.tsv`).

## Decisive finding

**Total P. vivax-mappable signal is abundant in the affected runs - by both
methods, and even proportionally *higher* than in the successful control
runs:**

| Run | STARsolo genome-mapped (Unique+Multi) | kallisto pseudoaligned (parasite-only index) |
|---|---:|---:|
| SRR19573609 (rep2, d9, **affected**) | 6.69% | **2.6%** (11,613,567 / 440,660,768 reads) |
| SRR19573610 (rep2, d5, **affected**) | 6.23% | **2.6%** (11,206,148 / 433,820,693 reads) |
| SRR19573611 (rep1, d9, control) | 0.74% | 0.6% (2,246,795 / 402,658,977 reads) |
| SRR19573612 (rep1, d5, control) | 2.54% | 2.0% (8,529,314 / 424,172,132 reads) |

Both the original STARsolo alignment and this investigation's independent
kallisto pseudoalignment agree: **the affected runs have as much or more raw
parasite-mappable sequence than the unaffected runs.** This rules out "low
intrinsic parasite content in the public deposit" as a sufficient
explanation by itself.

**But that signal is not attributable to the correct cell barcodes, under
either method:**

| | 609 known cells (n=89) | 610 known cells (n=449) | 611 known cells (n=255, control) | 612 known cells (n=636, control) |
|---|---:|---:|---:|---:|
| kallisto median UMI | 0 | (see D036_cell_level_rescue.tsv) | 604.0 | 1236.5 |
| kallisto zero-count fraction | 0.978 | (see D036_cell_level_rescue.tsv) | 0.000 | 0.000 |
| kallisto rank among barcodes-with-signal (median; 1=highest UMI) | 108,847 of 195,430 | (see table) | 128 of 120,599 | 318 of 166,988 |
| kallisto max single-cell UMI observed | 1 | (see table) | 14,115 | 101,342 |

The known real-cell barcodes for the **control** runs rank at or near the
very top of the kallisto UMI distribution (median rank in the low hundreds,
best rank = 1, zero cells with 0 UMI) - direct proof this investigation's
pipeline (index, whitelist, barcode correction, counting) correctly
recovers real single-cell signal when the barcode-to-cell mapping is valid.

The known real-cell barcodes for the **affected** runs are statistically
indistinguishable from background: median rank in the hundred-thousands (out
of ~170,000-195,000 barcodes carrying any signal at all), 98% exactly zero
UMI, maximum observed UMI of 1. The abundant parasite-mappable reads in
609/610 are instead spread diffusely across ~170,000-195,000 distinct
barcodes (mean ~12-14 UMI/barcode) with no barcode standing out as a
plausible real cell - and the top-UMI barcodes in each affected run do
**not** correspond to any of the 1,438 known real cell barcodes from *any*
of the four runs (ruling out a simple barcode/run swap between the four
Ruberto2022_1 FASTQ files themselves).

## Interpretation

This is not an alignment-sensitivity problem (kallisto's much more
permissive k-mer pseudoalignment against a curated transcriptome, with no
requirement to discover splice junctions de novo, gives the same null result
as STAR's stricter spliced genome alignment). It is not a truncated/corrupt
download (content-verified, exact read-count match to ENA). It is not a
chemistry/parameter misconfiguration (identical command line/whitelist/read
structure across all four runs, and the exact same `soloBarcodeMate`/
`clip5pNbases` convention that works correctly for every other PlaViSca
study). It is not a simple barcode/run swap among the four local FASTQ files
(the affected runs' top-signal barcodes match no known real cell from any
run).

The evidence is consistent with a **library-level or public-deposit-level
problem specific to biological replicate 2** that occurred upstream of any
alignment step - most plausibly, either:

- severe ambient/diffuse RNA contamination in the replicate-2 library
  preparation (parasite transcript fragments not properly compartmentalized
  into droplets, so they appear at low level across a huge number of
  barcodes rather than concentrated in the true cells), or
- a barcode/UMI-corrupting event specific to how these two runs were
  captured, multiplexed, or deposited (e.g. an index-hopping-like effect
  that decouples the true cell barcode from its correct transcript pool, at
  a scale far exceeding typical background hopping rates).

Distinguishing these two sub-mechanisms would require material this
investigation does not have access to (e.g. the authors' internal
pre-deposit intermediate files, or a second, non-SRA copy of the raw
replicate-2 data) and does not change the practical conclusion below.

## Why this rules out Option A (repaired from public raw data)

No credible reprocessing strategy on the public FASTQ can concentrate signal
that is not there to be concentrated: the pseudoalignment is not failing to
*find* parasite reads (it finds proportionally more of them than in the
successful control runs); it fails to find them *at the correct barcodes*.
Re-running STARsolo with different parameters, a different reference
version, a combined host+parasite index, or a different aligner entirely
would not change which barcode each read carries - that assignment is fixed
by the FASTQ content itself. This investigation's kallisto experiment is
effectively already the "try a fundamentally different, more sensitive
processing approach" step DEC09 asked for, and it returns the same negative
result at the cell-barcode level that STARsolo did.

**D036 is therefore classified B, not A: this is a source/provenance
limitation of the public deposit relative to the authors' processed object,
not a PlaViSca pipeline defect.** No further reprocessing of these two runs
is recommended.

## Author-count fallback evaluation (Part 10 - evaluation only, not applied)

`data/Hep59.1.2.seu_20aug2025.rds` (`RNA` assay, `counts` layer) was
inspected directly:

- Built with Seurat 4.0.3 (historical object).
- 4,722 *P. vivax*-only features (`PVP01-` prefixed, hyphenated - i.e.
  already restricted to the parasite side of the authors' combined
  human+parasite kallisto index), 1,438 cells - the exact retained
  population.
- `counts` values are integers (confirmed: `all(counts == round(counts))`
  on a sample), consistent with genuine raw UMI counts, not a normalized
  representation. `data` differs from `counts` (not identical), consistent
  with `data` holding a separately normalized layer alongside raw `counts`.
- Provenance is now well-established (not "unresolved" as Phase 1 reported):
  traced to `kallisto bus -x 10xv3` (kallisto 0.46.0) against
  `hs38_103ens_pvp01_51PlamsoDB_tr_index.idx` (Homo_sapiens.GRCh38.cdna.all
  + PlasmoDB-51_PvivaxP01_AnnotatedTranscripts, k=31), followed by
  `bustools correct|sort|count` (bustools 0.39.3) with `emptyDrops(lower=1000)`
  cell calling - see `Hep_KallistoBUS_10xV3_aln.sh`,
  `Homo38_v103ens_PvivaxP01_PlasmoDBv51_Index.sh`,
  `Hep_KallistoBUS_10xV3_count.sh`, and
  `ruberto_github_additional_file_1.Rmd` (Zenodo 10.5281/zenodo.6463338).
- This investigation's own kallisto re-alignment used the *same*
  methodology family (kallisto\|bustools, parasite transcriptome
  pseudoalignment) as the authors and reproduces the exact same null result
  for the 538 cells as STARsolo. This means the authors' own recovered
  signal for these cells did **not** come from applying the same method to
  the same public SRA deposit - it must derive from their internal working
  copy of the data (pre-deposit), a different FASTQ generation/version, or
  additional per-sample QC/processing not fully documented in the
  reproducible shell scripts. This is itself new evidence supporting a
  deposit-vs-internal-data provenance gap, not a PlaViSca-side error.
- Feature coverage: `Hep59.1.2.seu`'s 4,722 P. vivax genes are a subset of
  PlaViSca's 6,811/6,861 (post-rRNA-removal) gene set - consistent with the
  authors dropping never-expressed genes after their own cell-calling, not
  evidence of a fundamentally different gene model.
- Mixing these author-derived counts for 538/1,438 Ruberto2022_1 cells with
  independently STARsolo-derived counts for the other 900 Ruberto2022_1
  cells and all five other PlaViSca studies would introduce a real,
  documented pipeline inconsistency (different aligner, different reference
  version/scope, different cell-calling method, different software era) for
  a defined, flagged subset of one study's cells.

**Recommendation ranking** (per Part 10 - evaluation only; none of these are
applied in this task):

1. **Recommended default: option 3** - retain and explicitly flag the 538
   affected cells (already implemented: `near_empty_expression_flag`,
   `total_umi_count` in `ruberto2022_1_pv_analysis_script.R`), and exclude
   them from expression-dependent integrated analyses (differential
   expression, marker scoring, clustering-driven annotation) until a
   project-lead scientific decision is made. This preserves the Phase 1
   conservative membership policy and introduces no cross-pipeline
   inconsistency.
2. **Possible alternative, requires project-lead approval: option 2** - use
   the author-deposited `Hep59.1.2.seu` raw counts for specifically these
   538 cell IDs, with an explicit `source_of_counts` provenance field
   distinguishing them from the other 900 STARsolo-derived
   Ruberto2022_1 cells. Now technically well-documented (see provenance
   trace above) and feature-compatible (both indexed by the same
   PVP01-prefixed gene IDs), but knowingly mixes two different alignment
   methodologies within one study's retained population and should not be
   adopted without explicit sign-off.
3. **Not recommended: option 1** (further public-FASTQ reprocessing) - this
   investigation's controlled experiment already used a fundamentally
   different, more sensitive method on the same verified-intact FASTQ and
   found the same null result; further reprocessing attempts are very
   unlikely to succeed without new source material.
4. **Not recommended as a default: option 4** (exclude from atlas) - a
   membership-changing decision Phase 1's conservative policy explicitly
   defers; nothing in this investigation newly justifies overriding that
   policy.
5. **Superseded**: option 5 (keep D036 unresolved) - this investigation
   reaches a decisive characterization (classification B), so "unresolved"
   is no longer accurate; the remaining question is a scientific
   inclusion/data-source policy decision (1 vs 2 above), not further
   technical investigation.

## Part 0: environment fix (for reproducibility of this investigation)

`scripts/pixi.toml` gained an `[activation]` script
(`scripts/setup/ensure_bioc_data_packages.sh`) that idempotently installs
`GenomeInfoDbData` into the pixi environment's own R library. Root cause:
bioconda's `bioconductor-genomeinfodbdata` recipe fetches its actual (large)
data payload via a conda `post-link.sh` hook at install time; pixi does not
execute conda post-link/pre-unlink hooks (a deliberate pixi/rattler design
choice), so the payload was never downloaded even though `conda-meta`
reported the package installed. The fix reproduces exactly what that hook
would have done (download the official Bioconductor data tarball, verify
its MD5, `R CMD INSTALL` it into `$CONDA_PREFIX/lib/R/library`), scoped
entirely to the pixi-managed library - never a personal/user library - and
is a no-op if the package is already present (verified: first activation
after removal took ~8s to fetch+install; a second activation completed in
under 2s with the package already present). `scripts/tests/test_environment.R`
(`pixi run test-env`) is a smoke test confirming all 12 required packages
(including `DropletUtils`, `rtracklayer`, `GenomeInfoDb`, `GenomeInfoDbData`)
load. `pixi.lock` was not modified - no new dependencies were added, only
the activation-time completion of an existing declared dependency.

## Cell-membership and count decisions in this task

- **No cell-membership change.** Ruberto2022_1 remains exactly 1,438 cells.
- **No count value change.** `ruberto2022_1.rds` (and its production script)
  are unchanged by this investigation; `near_empty_expression_flag`/
  `total_umi_count` (added in Phase 1) remain accurate and are NOT removed,
  since the underlying near-empty condition is confirmed genuine, not a
  PlaViSca-side artifact.
- **No candidate repaired study object was created**, because no repair was
  found. This is consistent with the Section 9 decision tree's Outcome B:
  document the limitation, do not fabricate a computational rescue.

## Remaining Phase-2 blocker status

D036 (or more precisely, the now-better-characterized public-deposit
provenance gap it represents) **remains a blocker for the full Phase 2
atlas rebuild**, in the sense that the scientific inclusion decision
(recommendation ranking above) still requires project-lead sign-off before
the affected 538 cells' expression values can be trusted in any
expression-dependent downstream analysis. It is no longer an open technical
question; it is a scientific/policy question.

---

## PHASE 2 UPDATE (2026-09-20): project-lead sign-off and formal decision

The Phase-2 project lead has now made the scientific decision required by
the "Remaining Phase-2 blocker status" section above, on the basis of the
independent `audit/ruberto2022_1/author_matrix_validation/` study (commit
`c477913`), which supersedes the weaker "possible alternative" ranking of
option (b) given in Part 10 above.

**Decision:** Use the complete author-derived raw UMI matrix
(`data/Hep59.1.2.seu_20aug2025.rds`, `RNA` assay `counts` layer) as the
canonical expression source for **all 1,438** Ruberto2022_1 cells - not only
the 538 previously near-empty (D036-affected) cells. Do not mix author
counts with STARsolo counts within the study, and do not exclude the 538
D036 cells.

**Why (consolidated):**

1. Public FASTQ/STARsolo reconstruction is not reproducible for the 538
   cells (this document, Sections above).
2. Independent kallisto/bustools reprocessing, using the authors' own
   successful methodology on the same local FASTQ, did not rescue them
   either - the defect is upstream of alignment method (library/deposit
   level), so no further reprocessing is expected to help.
3. The author `RNA` counts layer is confirmed non-negative integer raw UMI
   counts (`author_matrix_validation/rna_counts_validation.tsv`:
   `noninteger_nonzero_entries: 0`, `SCT_assay_present: FALSE`).
4. The 1,438-cell crosswalk (`51inf->612, 52inf->610, 91inf->611,
   92inf->609`) is exact and one-to-one, with zero unmapped or duplicated
   IDs (`author_matrix_validation/cell_crosswalk_validation.tsv`).
5. The 900 unaffected cells show strong positive expression/cell-identity
   concordance between the author matrix and PlaViSca's own STARsolo
   reconstruction (`cell_identity_test.tsv`: true-counterpart rank-1 match
   85-100% depending on run, versus a 0.16-0.39% chance rate) - this
   supports that, for the unaffected population, the author matrix
   measures the same biological cells PlaViSca's pipeline measures, not a
   relabeled or different population, which in turn supports trusting the
   author matrix for the 538 cells where PlaViSca has no signal to compare
   against.
6. All 4,722 author PVP01 gene IDs map exactly and unambiguously into
   PlaViSca's 6,811-feature panel (`author_gene_mapping.tsv`: 4722/4722
   `exact_string_match`, 0 ambiguous/unmapped/duplicate).

**What this decision does NOT resolve** (explicit residual caveat, carried
forward unchanged): byte-level identity of the locally dated object
(`Hep59.1.2.seu_20aug2025.rds`) to the canonical Zenodo-deposited RDS
remains unconfirmed. This is a provenance-completeness limitation, not a
reason to distrust the validated counts - it is recorded per-cell via the
new `count_provenance_status` field rather than silently dropped. The
public-data discrepancy characterized in `audit/ruberto2022_1/AUDIT.md` is
**not** being called "fixed" by reprocessing; it remains an unreproducible
public-deposit defect. It is resolved *operationally*, for the purpose of
building the Phase-2 atlas, by adopting the independently validated author
matrix as this study's expression source, with explicit, auditable
provenance carried at cell/study metadata level.

**Implementation:** `scripts/ruberto2022_1_pv_analysis_script.R` is updated
to read `RNA` counts from `data/Hep59.1.2.seu_20aug2025.rds` for all 1,438
cells (after the existing crosswalk-based subset/rename), zero-fill only
the PlaViSca features classified `not_part_of_author_reference` in
`audit/ruberto2022_1/phase2_feature_coverage.tsv` (never interpreting
reference-absent genes as unconditional biological zero without that
evidence trail), and populate `source_of_counts`, `counting_pipeline`,
`count_reference_version`, `count_provenance_status` for every cell. D036
is closed as a Phase-2 blocker; `near_empty_expression_flag`/
`total_umi_count` are retained as historical/diagnostic fields describing
the now-superseded STARsolo reconstruction, not the adopted counts.
