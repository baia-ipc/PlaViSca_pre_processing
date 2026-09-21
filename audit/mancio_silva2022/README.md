# Mancio-Silva2022 audit package

Start with [AUDIT.md](AUDIT.md). The repository includes the audit report,
duplicate and lineage tables, diagnostic R scripts, validation outputs, historical
script snapshots, accession records, and selected upstream count columns.

The approximately 2.5 GB of downloaded GEO archives and full extracted count
matrices remain local and are excluded from Git. Their exact filenames, retrieval
URLs, extraction members, sizes, and SHA-256 checksums are recorded in
`source_download_manifest.tsv`. No source files were deleted.

To restore these inputs, download the three rows whose kind is `download` from
their listed URLs into this directory. For each `archive_member` row, extract
only its named member from the listed archive into this directory. The smaller
gene and cell-ID files are already included in Git. Verify restored files against
`SHA256SUMS.txt` before running the diagnostic commands in AUDIT.md.

From this directory, `sha256sum --check --ignore-missing SHA256SUMS.txt` verifies
the available evidence; use `sha256sum --check SHA256SUMS.txt` once all source
archives and matrices have been restored. Production input checksums are separate
in `production_inputs_SHA256SUMS.txt`, with paths relative to the project root.

Files named `historical_*.R` are archived evidence, not runnable entrypoints:
they contain original production paths and save operations. Run only the
`inspect_*.R` diagnostics as documented, using the scripts Pixi environment.

This package does not contain repairs to production code or regenerated
production datasets.

## Shared-pipeline consolidation audit

Once all six study-specific audits (this package plus `hazzard2024/`,
`ruberto2022_1/`, and `ruberto2022_2/`, with the Sa2020 and Hazzard2022
sections recorded inside this package's `AUDIT.md`) were complete, a
separate cross-study consolidation audit was performed at
[`../shared_pipeline/AUDIT.md`](../shared_pipeline/AUDIT.md). It builds the
authoritative defect register, the shared-mechanism analysis for
`singleR.R`/`integration.R`/`flatten_data.R`, the target metadata schema, and
the staged repair plan. Start there before beginning any repair work.
