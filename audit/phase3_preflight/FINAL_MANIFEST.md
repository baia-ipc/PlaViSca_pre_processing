# Authoritative final-release manifest

> **DEC01/DEC02 update:** the authoritative TSV now selects the 105,949-cell
> rebuild from scientific commit `173431cfca8e19c716a96164ee2ea6053441855b`.
> Focused resolution classified the original AT25 result as a validation
> semantics defect; corrected AT25 passes with unchanged thresholds, yielding
> 30 PASS / 0 FAIL. `REPORT.md` contains the current verdict. Statements below
> describing the prior 105,963-cell
> closure are retained as historical context and are superseded.

`final_release_manifest.tsv` is the authoritative manifest for the PlaViSca v1 data release candidate. The historical root `build_manifest.tsv` is preserved unchanged as an append-only account of intermediate Phase-2 generations; its repeated rows are not authoritative selectors.

The manifest cites technical code commit `8ce8828985aa75ea625bb4e762465f6ab8689905`. The final documentation/acceptance commit is a descendant and is verified by the pushed branch SHA, avoiding an impossible self-referential commit hash inside its own contents.

Every artifact row marked `AVAILABLE_AND_VERIFIED` was hashed directly from the exact file named in `path`. Inputs earlier than the six final study objects are not asserted to be exhaustively recoverable here; the committed per-study forensic audits and historical manifest remain their provenance trail. No unavailable historical checksum is invented.

The scientific atlas was not recomputed in release closure. Only `cleaned_dataset.rds` was rewritten to canonicalize cell order, and `data_source.csv` was extended with verified PMIDs. Raw, normalized, scaled, study-object, and atlas hashes are unchanged.
