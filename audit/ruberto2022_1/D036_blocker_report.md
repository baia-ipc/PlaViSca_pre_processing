# D036 blocker report: Ruberto2022_1 replicate-2 near-total expression loss

Status: **investigated, not resolved. Remains a hard Phase 2 blocker.**
Date: 2026-09-19 (Phase 1 production repair).

## Summary

538 of 1,438 Ruberto2022_1 cells (overwhelmingly biological replicate 2:
runs `SRR19573609` and `SRR19573610`) have near-total expression loss (517
cells with exactly 0 UMI, 21 cells with 1 UMI). PlaViSca's own
STARsolo-output-vs-study-object concordance is exact for these cells (ruling
out an R-level import/merge bug - already established in
`audit/ruberto2022_1/AUDIT.md`). The authors' own local object
(`data/Hep59.1.2.seu_20aug2025.rds`) shows real signal for the identical
cell IDs (median nCount_RNA 718.5).

This Phase 1 investigation went one level deeper: comparing FASTQ/alignment-
level evidence for the two affected runs (609/610) against the two
unaffected runs (611/612) from the same study, using STARsolo's own logs
and output matrices (no raw FASTQ is stored in this repository).

## Root-cause investigation findings

### Data available vs unavailable

- Raw FASTQ: **not present locally** anywhere in the repository. All checks
  below use registered ENA metadata (`audit/ruberto2022_1/ena_PRJNA843856_runs.tsv`)
  compared against STARsolo's own "Number of input reads" - no direct
  read/quality inspection of FASTQ content was possible.
- STARsolo logs/outputs: **present and complete for all 4 runs** -
  `Log.final.out`, `Log.out` (full command line + per-barcode filtering
  diagnostics), `Solo.out/GeneFull/{raw,filtered}/{barcodes.tsv,features.tsv,matrix.mtx}`,
  `Summary.csv`, `UMIperCellSorted.txt`, `Barcodes.stats` all exist under
  `counts/36093191/SRR19573609{,10,11,12}_solo_out/`.
- No alignment wrapper/launcher script exists in this repository (only
  `counts/` outputs); the STARsolo command line was recovered from
  `Log.out` itself.
- Neither `production_inputs_md5.tsv` nor `SHA256SUMS.txt` contains FASTQ
  checksums - FASTQ provenance/integrity is not independently verifiable
  from repository evidence alone.

### Quantitative comparison (609/610 = replicate 2/affected vs 611/612 = replicate 1/unaffected)

| Metric | 609 (rep2,d9) | 610 (rep2,d5) | 611 (rep1,d9) | 612 (rep1,d5) |
|---|---:|---:|---:|---:|
| ENA registered read_count | 440,660,768 | 433,820,693 | 402,658,977 | 424,172,132 |
| STARsolo "Number of input reads" | 440,660,768 | 433,820,693 | 402,658,977 | 424,172,132 |
| -> exact match to ENA, both runs | yes | yes | yes | yes |
| Uniquely mapped reads % (genome) | 3.47% | 3.15% | 0.55% | 2.28% |
| Reads mapped Unique+Multi % | 6.69% | 6.23% | 0.74% | 2.54% |
| Reads with valid barcodes | 97.72% | 97.88% | 98.89% | 98.77% |
| Sequencing saturation | 0.576 | 0.677 | 0.180 | 0.307 |
| STARsolo auto "Estimated Number of Cells" | 61,444 | 42,071 | 276 | 156 |
| Mean/Median UMI per (auto) cell | 88 / 84 | 82 / 73 | 305 / 179 | 5,852 / 4,706 |
| Mean/Median GeneFull per cell | 2 / 2 | 3 / 3 | 221 / 154 | 1,667 / 1,680 |
| Total GeneFull genes detected | 2,900 | 3,930 | 4,382 | 4,886 |
| **Max single-barcode UMI, entire raw matrix (6,794,880 barcodes)** | **1,593** | **2,179** | 4,402 | 32,398 |
| Raw matrix dims (genes x barcodes) | identical (6,861 x 6,794,880) across all 4 runs | | | |
| STARsolo command line (chemistry, whitelist, clip params) | identical across all 4 runs except run ID | | | |

Author's local object, for the 538 cells PlaViSca shows as near-empty (89
from 609, 449 from 610): nCount_RNA median 718.5, mean 1,523.6, **max
18,273**.

### Assessment

- **Ruled out: FASTQ truncation/incomplete download.** STARsolo's own input
  read count matches the ENA-registered `read_count` exactly, for all 4
  runs, to the read.
- **Ruled out: chemistry/barcode-length/whitelist parameter mismatch.** The
  four `Log.out` "Final effective command line" blocks are byte-identical
  except for the run ID/output prefix. `Barcodes.stats` valid-barcode
  fractions (97.7-98.9%) are essentially equal across all 4 runs.
- **Ruled out as the sole/primary mechanism: a downstream cell-selection/
  knee-threshold artifact.** STARsolo's automatic knee filter is badly
  miscalibrated for 609/610 (tens of thousands of "cells" at mean UMI/cell
  of only ~85), but this cannot be the mechanism behind D036 itself: the
  study's 1,438 retained cells come from a fixed author-supplied barcode
  list, not STARsolo's auto-filter, and decisively, **the maximum UMI value
  for any single barcode anywhere in the entire 6,794,880-barcode raw
  matrix is only 1,593 (609) / 2,179 (610)**, while the authors' local
  object attributes up to 18,273 UMI to specific cells with those same
  barcode IDs. No re-thresholding of PlaViSca's own STARsolo output can
  produce a value that never appears in the matrix at all - this is a
  magnitude/ceiling discrepancy, not a selection-boundary discrepancy.
- **Most consistent with:** a genuine alignment-output deficiency specific
  to these two runs, of unresolved deeper origin. Two sub-possibilities
  cannot be distinguished with evidence available in this repository:
  1. The FASTQ data registered/deposited under SRR19573609/610 in SRA/ENA
     itself carries far less recoverable *P. vivax* parasite signal than
     whatever the authors used to build their local object (a
     provenance mismatch between the SRA deposit and the authors' own
     working data) - all technical metrics (read count, valid-barcode
     rate) look normal, i.e. these are well-formed reads that simply do
     not align to the parasite reference at anywhere near the rate implied
     by the authors' counts.
  2. PlaViSca's alignment reference/parameters recover substantially less
     signal than whatever pipeline/reference the authors used, specific to
     replicate 2 - considered less likely as the *sole* explanation, since
     the identical reference/parameters were used for replicate 1
     (611/612), which shows STARsolo signal of the same broad order of
     magnitude as the authors report.
  3. No FASTQ swap between runs is evident: barcode format, whitelist match
     rate, and cell-ID mapping (`hep59_to_study_cellid_map.tsv`, verified
     exact 1:1 in the prior audit) are all consistent per-run.
- **Confidence:** High confidence the cause is upstream of R/PlaViSca's
  import and merge logic (already proven exact) and specific to the
  STARsolo count-generation step for these two runs. Moderate confidence it
  reflects either an SRA-deposit-level data problem or an
  alignment-reference/pipeline difference versus the authors' own
  (undocumented) processing - this specific sub-mechanism is **inconclusive
  from available evidence** and cannot be resolved without either the raw
  FASTQ content itself or the authors' exact alignment protocol.

### Prioritized next steps (Phase 2)

1. Fresh-download FASTQ for SRR19573609/610 directly from ENA
   (`ena_PRJNA843856_runs.tsv` fastq_ftp URLs) and spot-check read content
   (not just count) - confirm the copy used for this STARsolo run was not
   itself corrupted despite matching read counts.
2. Run FastQC/read-composition estimation (e.g. kraken2, or a BLAST sample)
   on 609/610 vs 611/612 FASTQ: is replicate 2 predominantly host (human)
   reads with genuinely minimal parasite content, or does it have
   recoverable parasite signal that PlaViSca's index/parameters are failing
   to capture?
3. Contact study authors (or search PMID 36093191 supplementary methods)
   for their exact alignment reference/version and per-accession cell
   counts, to determine whether their local object was built from these
   same SRA accessions or from an internal, non-deposited FASTQ/processing.
4. If (2) confirms genuinely low intrinsic parasite signal in the deposited
   replicate-2 FASTQ, this is very likely **not fixable by reprocessing**
   with different STARsolo parameters - the raw data itself lacks the
   signal, and DEC09 option (c) (flag/retain with caveat, or exclude per a
   separately approved decision) becomes the only remaining path.

## Phase 1 disposition

Per DEC09 ("pursue investigation/realignment first; fall back to explicit
flagging if reprocessing does not resolve the issue within a reasonable
effort; do not pursue substituting the authors' counts without first
resolving the author object's Zenodo provenance"):

- Root-cause investigation was performed (this report).
- Read counts already match the registered SRA deposit exactly (rules out a
  simple truncated-download fix).
- Chemistry/parameters are already correct and identical to the unaffected
  runs (rules out a parameter-tuning fix).
- The defect is a magnitude/ceiling-level discrepancy versus the authors'
  data, not a threshold/selection issue - none of these are resolvable by
  adjusting STARsolo flags on the existing FASTQ.
- **D036 is NOT repaired in Phase 1.** No FASTQ realignment was attempted
  (Phase 1 explicitly forbids the atlas-wide rebuild this would gate, and
  the investigation above shows realignment alone is unlikely to resolve
  it without first confirming FASTQ-content composition - step 2 above).
- `scripts/ruberto2022_1_pv_analysis_script.R` now computes `total_umi_count`
  and `near_empty_expression_flag` per cell and asserts the flagged count
  equals 538 (warning, not silent, if it drifts), per the conservative
  Phase 1 membership policy: **the 538 cells remain in the 1,438-cell
  population, explicitly flagged, not excluded.**
- D036 remains the single hard Phase 2 blocker preventing the full
  atlas-wide Harmony/SingleR/UMAP rebuild, exactly as this phase's task
  specification anticipated.
