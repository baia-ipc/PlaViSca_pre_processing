# Phase-3 preflight independent QA report

## Verdict

The Phase-2 candidate is scientifically plausible and suitable as the **data basis** for an app repair, but it is **not an acceptance-complete release candidate**. Independent closure yields 23 PASS, 4 FAIL, 1 PENDING_DECISION, and 2 PENDING_APP. The blocking evidence gaps are AT02 (Sa2020 lineage), AT17 (cleaned-metadata row order), AT23 (manifest provenance), and AT25 (no quantitative marker-preservation pass criterion). DEC01/DEC02 keep AT14 pending. The application itself is not declared ready because no immutable remote repair branch was available.

## Scope and immutable inputs

- QA base: `54c6b9581ff5268422ddffda42ce457bf09c5b55`
- QA branch: `plavisca-phase3-qa`
- Isolated worktree: `/home/baia/prj/plavisca/pre_process_data_phase3_qa`
- Candidate binaries were read only from `/home/baia/prj/plavisca/pre_process_data`.
- App master baseline SHA: `020f9cbd651750ba92c037399fb84d3b17d389e8`, inspected read-only before the app agent switched its worktree.

Initial ETA was 3–4 hours. The evidence review completed within that window; the remaining time was used for immutable-branch recheck, artifact validation, and commit hygiene.

## Acceptance closure

The canonical result is `acceptance_matrix_AT01_AT30.tsv`. New direct evidence was created for AT02, AT12, AT17, AT23, and AT28; static/semantic reviews directly closed AT04, AT05, and AT24–AT30. AT15's criterion was minimally amended to recognize the already approved DEC09 complete-author-matrix repair path.

- **AT02 — FAIL:** five studies have complete committed source→study→final evidence. Sa2020's broadest committed cell-level lineage has 9,432 present unique IDs for 9,766 retained cells, leaving 334 without an explicit source→study row.
- **AT04 — PASS:** production assignments use `ncol`-derived exact lengths or keyed helpers, with fail-loud assertion primitives and committed unit evidence. The current host's R 4.6.1 cannot load older compiled `rlang`/Seurat packages, so the unit suite could not be replayed in this host environment.
- **AT05 — PASS:** run/library/day/strain variation is keyed. Remaining `rep()` calls are full `n_cells` constants or replication of a value selected from a unique keyed row; no dangerous short-vector recycling was found.
- **AT12 — PASS:** all 80,024 Hazzard2024 cells agree with the authoritative run table for library, animal, parasite lineage, day, and infection design.
- **AT14 — PENDING_DECISION:** the candidate follows the documented retain-and-flag policy, but the acceptance language requires the eventual DEC01/DEC02 decision, which does not exist yet.
- **AT15 — PASS:** exactly 1,438 cells use the validated complete author matrix; all 4,722 shared features are exactly equal, the D036 near-empty artifact is absent, and provenance is populated for all cells.
- **AT17 — FAIL:** all artifacts have identical unique cell-key sets. Raw, normalized, and scaled exports share one sorted order, while `cleaned_dataset.rds$mr_data` has a different order. Any positional app join is unsafe.
- **AT23 — FAIL:** no input checksums are consolidated in the manifest; build rows cite the build-start SHA; artifacts have repeated historical rows; the current atlas hash matches no row; raw/normalized/scaled exports have no rows. The supplemental manifest records current bytes only and does not reconstruct missing provenance.
- **AT28 — PASS:** for all 90,219 retained IDC labels, raw label, delta, and pruned flag are non-missing and no retained label is marked pruned.
- **AT29 — PASS with caveat:** the saved-graph clustering rerun produced identical numeric labels and ARI 1. It tests fixed-graph clustering, not end-to-end Harmony/UMAP, and relies on the clustering function's default seed rather than an explicit argument.
- **AT30 — PASS with caveat:** 22,905/23,689 (96.6904%) original gametocyte calls retain the same marker-threshold-defined sex call under resolution perturbation. This is biological marker-call stability, not numeric cluster-ID stability.

AT24, AT26, and AT27 have semantically relevant kNN/neighborhood evidence. AT25 fails because its output only tabulates four marker means by broad stage; it has no expected direction, quantitative threshold, or assertion, while the Phase-2 suite merely checked that the combined metrics file existed.

## Candidate binaries and contract

Eight required artifacts match an identifiable final manifest row: six study RDS files, the cleaned dataset, and corrected bibliography CSV. The atlas is a hash mismatch; three expression exports are absent from the manifest. See `candidate_hash_validation.tsv`.

Direct candidate facts:

- 105,963 cells across six studies.
- 88 cleaned metadata columns.
- Raw and normalized exports: 6,860 gene columns; scaled export: 2,000 gene columns.
- Gene columns are cleanly distinguished from the 88 metadata columns and all use `PVP01_` identifiers.
- Zero literal `"NA"` strings in desired app-visible fields.
- 90,219 IDC-labeled cells, 23,689 gametocyte calls, and zero liver-stage gametocyte calls.
- Ruberto feature categories are confirmed as 4,722 measured, 288 biotype-excluded, and 1,801 inferred release-gap genes. Integration uses the 5,203-gene six-study intersection before selecting 2,000 HVGs, so structural zero-padding does not enter PCA/Harmony.

The corrected bibliography has exactly one row for each of the six candidate study labels and matching cell counts. It has **no PMID column**, so candidate `study_pmid` agreement cannot be tested against the actual file. This is recorded as `LABEL_MATCH_PMID_NOT_IN_BIBLIOGRAPHY`, not PASS.

## Phase-2 report audit

Counts, Ruberto replacement, integration feature space, annotation totals, impossible-state claims, candidate locations, and the saved island-dispersion comparison are supported. The claims that all hashes are in the manifest, that every named acceptance criterion was exercised, and that full acceptance validation was complete are incorrect. Historical non-overwrite is only partially independently verifiable. See `report_claim_validation.tsv`.

## Application baseline and repair branch

Current master has obsolete config mappings, direct obsolete column references, old HPI vocabulary, stale `functions/preprocessing.R`, no bibliography foreign-key validation, hardcoded label/color maps, an unreachable `blood_stage` branch, and silent error fallbacks. It cannot consume the candidate unchanged.

At the final remote check, `github/plavisca-app-production-repair` did not exist. The separate app agent had a local branch with uncommitted work; QA deliberately did not review that mutable state. Therefore AT19 and AT20 remain `PENDING_APP`, no repaired-app SHA was reviewed, and no claim is made that the application is ready.

## Reproduction notes

- Run `Rscript audit/phase3_preflight/run_candidate_qa.R` from the QA worktree.
- Run `audit/phase3_preflight/verify_candidate_hashes.sh` for independent SHA256 output.
- Production R scripts parse successfully with base R. Full package-level tests require the Phase-2 R 4.3.3 environment (or rebuilt R 4.6-compatible packages); the current user library has an ABI mismatch.
