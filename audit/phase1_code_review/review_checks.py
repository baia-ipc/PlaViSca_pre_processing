#!/usr/bin/env python3
"""Small deterministic, read-only checks for the Phase-1 code review."""

import csv
import hashlib
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent


def read_tsv(path):
    with open(path, newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, rows, fields):
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


ref_rows = []


def ref(table, check, observed, expected, status, evidence):
    ref_rows.append({
        "table": table,
        "check": check,
        "observed": observed,
        "expected": expected,
        "status": status,
        "evidence": evidence,
    })


# Hazzard2024: compare fields having direct counterparts in the authoritative table.
h_auth = {r["run_id"]: r for r in read_tsv(ROOT / "audit/hazzard2024/authoritative_run_table.tsv")}
h_prod = {r["run_id"]: r for r in read_tsv(ROOT / "scripts/ref/hazzard2024_run_table.tsv")}
ref("hazzard2024_run_table.tsv", "run key set", len(h_prod), len(h_auth),
    "PASS" if h_prod.keys() == h_auth.keys() else "FAIL", "exact run_id set comparison")
mapping = {
    "library_id": "library_id",
    "animal_id": "animal_id",
    "day_post_infection": "day_post_infection_authoritative",
    "infection_design": "infection_design_authoritative",
    "host_species": "host_species_authoritative",
    "sequencer_registered": "instrument_model",
    "sequencer_protocol_publication": "sequencing_platform_authoritative",
    "host_sex": "host_sex_authoritative",
    "registry_location": "registry_geo_loc_name",
    "host_provenance": "host_provenance_note",
}
for prod_field, auth_field in mapping.items():
    bad = [rid for rid in h_prod if h_prod[rid][prod_field] != h_auth[rid][auth_field]]
    ref("hazzard2024_run_table.tsv", prod_field, len(bad), 0,
        "PASS" if not bad else "FAIL", ";".join(bad[:10]) or f"matches {auth_field} for all 24 runs")
for key in ("run_id", "library_id"):
    vals = [r[key] for r in h_prod.values()]
    ref("hazzard2024_run_table.tsv", f"unique {key}", len(vals) - len(set(vals)), 0,
        "PASS" if len(vals) == len(set(vals)) else "FAIL", "duplicate-key count")

# Mancio: compare the committed production crosswalk to the audited relation.
m_prod = read_tsv(ROOT / "scripts/ref/mancio_silva2022_run_table.tsv")
m_auth = {r["orig.ident"]: r for r in read_tsv(ROOT / "audit/mancio_silva2022/source_group_accession_map.tsv")}
mmap = {
    "source_group": "Sample",
    "geo_accession": "geo_accession",
    "geo_title": "geo_title",
    "sequencer_registered": "registered_instrument",
    "run_accession": "run_accessions",
    "experiment_accession": "experiment_accession",
    "source_cell_count": "cells",
}
ref("mancio_silva2022_run_table.tsv", "row/key count", len(m_prod), len(m_auth),
    "PASS" if len(m_prod) == len(m_auth) and len({r['prefix_key'] for r in m_prod}) == len(m_prod) else "FAIL",
    "15 unique prefix_key rows expected")
for prod_field, auth_field in mmap.items():
    bad = [r["orig_ident"] for r in m_prod if r[prod_field] != m_auth[r["orig_ident"]][auth_field]]
    ref("mancio_silva2022_run_table.tsv", prod_field, len(bad), 0,
        "PASS" if not bad else "FAIL", ";".join(bad) or f"matches {auth_field} for all 15 groups")
bad_day = [r["orig_ident"] for r in m_prod if r["day_post_infection"] != m_auth[r["orig_ident"]]["Day"].lstrip("d")]
bad_tx = [r["orig_ident"] for r in m_prod if r["source_treatment"] != m_auth[r["orig_ident"]]["Treatment"]]
ref("mancio_silva2022_run_table.tsv", "day_post_infection", len(bad_day), 0, "PASS" if not bad_day else "FAIL", "source Day comparison")
ref("mancio_silva2022_run_table.tsv", "source_treatment", len(bad_tx), 0, "PASS" if not bad_tx else "FAIL", "source Treatment comparison")
multi = sum(";" in r["run_accession"] for r in m_prod)
ref("mancio_silva2022_run_table.tsv", "multi-run relation representation", multi, 8, "PASS" if multi == 8 else "FAIL",
    "Eight Infection2 groups retain two associated runs; production misuse as scalar per-cell run_id is reviewed separately")

# Sa relation table.
sa = read_tsv(ROOT / "scripts/ref/sa2020_main_analysis_cell_lineage.tsv")
ref("sa2020_main_analysis_cell_lineage.tsv", "author rows", len(sa), 9215, "PASS" if len(sa) == 9215 else "FAIL", "audit main-analysis count")
present = sum(r["plavisca_present"] == "TRUE" for r in sa)
ref("sa2020_main_analysis_cell_lineage.tsv", "present in PlaViSca", present, 9018, "PASS" if present == 9018 else "FAIL", "audited overlap")
for key in ("author_cell_id", "plavisca_cell_id"):
    vals = [r[key] for r in sa]
    ref("sa2020_main_analysis_cell_lineage.tsv", f"unique {key}", len(vals) - len(set(vals)), 0,
        "PASS" if len(vals) == len(set(vals)) else "FAIL", "duplicate count")

for table in ("hazzard2024_run_table.tsv", "mancio_silva2022_run_table.tsv", "sa2020_main_analysis_cell_lineage.tsv"):
    ref(table, "sha256", sha256(ROOT / "scripts/ref" / table), "recorded", "INFO", "review snapshot checksum")

write_tsv(OUT / "ref_table_validation.tsv", ref_rows,
          ["table", "check", "observed", "expected", "status", "evidence"])

# Curated schema implementation assessment. "Partial" is used where a field is
# emitted only by some studies, under an obsolete alias, or at the wrong granularity.
schema = read_tsv(ROOT / "audit/shared_pipeline/metadata_schema.tsv")
assessment = {
    "study_id": ("NO", "study_label alias remains; canonical slug field study_id is not emitted"),
    "run_accession": ("PARTIAL", "run_id remains per-cell; Mancio multi-run sets are semicolon strings, not a relation"),
    "library_id": ("PARTIAL", "implemented for Hazzard2024 only"),
    "animal_id": ("PARTIAL", "implemented for Hazzard2024 and Sa2020; nullable elsewhere"),
    "donor_id": ("PARTIAL", "implemented for Ruberto2022_1; Mancio remains true NA"),
    "source_cell_id": ("PARTIAL", "implemented only for Mancio; other study source IDs remain barcode aliases"),
    "reconstructed_barcode": ("PARTIAL", "implemented only for Mancio; legacy barcode retained elsewhere"),
    "parasite_lineage_isolate": ("PARTIAL", "not emitted consistently by Mancio/Hazzard2022"),
    "per_cell_parasite_genotype": ("PARTIAL", "Hazzard2024 run table has values only for mono-lineage runs; processed per-cell source columns are dropped"),
    "host_species": ("YES", "present for all studies"),
    "host_taxid": ("PARTIAL", "proper integer where populated; Sa2020 remains NA despite known species"),
    "tissue_or_sample_type": ("YES", "canonical vocabulary emitted by study scripts"),
    "source_stage_annotation": ("NO", "no exact schema field; heterogeneous source_state/source_ap2g/liver_form fields remain"),
    "parasite_origin_location": ("PARTIAL", "field present, but Mancio hard-codes unsupported Ubon-Ratchathani province"),
    "host_provenance": ("PARTIAL", "not emitted consistently by all study scripts"),
    "experimental_site": ("PARTIAL", "not emitted consistently; Hazzard2024 table uses a curated institutional label"),
    "sequencing_site": ("PARTIAL", "not emitted consistently; true NA where unresolved"),
    "registry_location": ("PARTIAL", "not emitted consistently"),
    "day_post_infection": ("YES", "integer/true NA assignments; keyed where repaired"),
    "treatment": ("PARTIAL", "harmonized values mostly implemented, but no central vocabulary validator"),
    "source_treatment": ("PARTIAL", "preserved where available; absent in some scripts"),
    "infection_design": ("PARTIAL", "implemented for Hazzard2024 only"),
    "biological_replicate_id": ("YES", "implemented in all reviewed study scripts"),
    "biological_replicate_type": ("YES", "implemented in all reviewed study scripts"),
    "technology": ("PARTIAL", "Hazzard2024 uses technology; other studies use obsolete sc_technology alias"),
    "sequencer_registered": ("PARTIAL", "present in most scripts; Sa2020 retains sequencer alias"),
    "sequencer_protocol_publication": ("PARTIAL", "implemented for Hazzard studies; missing elsewhere"),
    "study_pmid": ("YES", "present and corrected"),
    "source_run_count_deposited": ("INCORRECT", "stored per-cell despite study-level schema"),
    "source_run_count_imported": ("INCORRECT", "stored per-cell despite study-level schema"),
    "source_run_count_retained": ("INCORRECT", "stored per-cell; Sa uses cell count 9766 rather than run count"),
    "source_stage_provenance": ("PARTIAL", "missing upstream for Hazzard2022/Hazzard2024/Sa2020 and synthesized late"),
    "idc_reference_similarity_label": ("BLOCKED", "implemented in singleR.R but eligibility vector is invalid for most studies"),
    "idc_similarity_raw_label": ("BLOCKED", "implemented but same execution blocker"),
    "idc_similarity_score_delta": ("BLOCKED", "implemented but same execution blocker"),
    "idc_similarity_pruned_flag": ("BLOCKED", "implemented but same execution blocker"),
    "idc_heuristic_stage_bin": ("BLOCKED", "implemented but same execution blocker"),
    "source_life_cycle_stage": ("PARTIAL", "only Ruberto studies populate it; Hazzard2022 source sporozoites do not"),
    "harmonized_life_cycle_stage": ("INCORRECT", "gametocyte inference only fills NA and therefore cannot supersede existing IDC bins"),
    "parasite_broad_stage": ("BLOCKED", "deterministic mapping is sound for characters, but upstream harmonized stages are incomplete"),
    "pred_gametocyte_sex": ("BLOCKED", "undefined writer and character-index eligibility bug prevent marker summary/calls"),
    "adjudicated_display_stage": ("YES", "present as nullable placeholder; no adjudications claimed"),
}
schema_rows = []
for row in schema:
    status, note = assessment.get(row["field_name"], ("UNREVIEWED", "No curated assessment entry"))
    schema_rows.append({
        "field_name": row["field_name"],
        "required_or_nullable": "nullable" if row["nullable"].lower().startswith("yes") else "required",
        "implementation_status": status,
        "type_granularity_provenance_review": note,
    })
write_tsv(OUT / "schema_implementation_check.tsv", schema_rows,
          ["field_name", "required_or_nullable", "implementation_status", "type_granularity_provenance_review"])

# Per-defect independent verification. Locations are intentionally concise;
# REVIEW.md contains the cross-cutting blockers and full consequences.
V = {}
def verdict(did, location, summary, matches, regression, classification):
    V[did] = (location, summary, matches, regression, classification)

verdict("D001", "scripts/silva.R:68", "Correct PMID constant.", "yes", "none found", "VERIFIED FIX")
verdict("D002", "scripts/silva.R:53-78", "Uses a keyed group table, but writes semicolon-delimited run sets into scalar per-cell run_id.", "partial", "violates run relation/scalar accession schema", "PARTIAL FIX")
verdict("D003", "scripts/silva.R:85-88", "Keyed source day converted to integer.", "yes", "none found", "VERIFIED FIX")
verdict("D004", "scripts/silva.R:90-98", "Preserves CTRL/PI4K and maps to study_control/PI4K_inhibitor.", "yes", "none found", "VERIFIED FIX")
verdict("D005", "scripts/silva.R:74-83", "Registered instrument is keyed; publication protocol/conflict is left NA without an explicit conflict flag.", "partial", "conflict provenance remains implicit", "PARTIAL FIX")
verdict("D006", "scripts/silva.R:100-105", "Splits host and parasite geography but assigns Ubon-Ratchathani to every cell despite country-only authority.", "no", "turns unresolved province into asserted fact", "INCORRECT FIX")
verdict("D007", "scripts/pipeline_lib.R:18-25; scripts/silva.R:71", "Shared preprocessing constant is used, but the app repository/config is outside this diff.", "partial", "end-to-end spelling consistency unverified", "PARTIAL FIX")
verdict("D008", "scripts/silva.R:29-38,124-138", "Only selected fields are copied; 14 .1 records fail exact Updated_names match and lose source metadata.", "partial", "incomplete source provenance", "PARTIAL FIX")
verdict("D009", "scripts/silva.R:31-37,133-155", "Exact match is attempted, but Replicative is blanket-converted to liver schizont and source topcell is not preserved.", "no", "changes established source/inference precedence", "INCORRECT FIX")
verdict("D010", "scripts/silva.R:130-131,157-159", "Source ID and reconstructed token are retained; final candidate export has not been rebuilt/checked.", "partial", "none beyond unverified propagation", "PARTIAL FIX")
verdict("D013", "scripts/silva.R:117-121", "num_srr is removed but three run counts remain replicated per cell; retained count is also 58 rather than the 23 associated runs.", "no", "schema granularity and value error", "INCORRECT FIX")
verdict("D014", "scripts/singleR.R:460-463", "Uses %in% for sex labels.", "yes", "execution blocked elsewhere", "VERIFIED FIX")
verdict("D015", "scripts/singleR.R:439-451", "Correct NIH/PB_MACS prefix map.", "yes", "execution blocked elsewhere", "VERIFIED FIX")
verdict("D016", "scripts/sa2020_pv_analysis_script.R:217-219", "Correct unqualified 10x_Chromium label.", "yes", "none found", "VERIFIED FIX")
verdict("D017", "scripts/sa2020_pv_analysis_script.R:213-219,263-275", "Separates registry/experimental site from strain origin.", "yes", "none found", "VERIFIED FIX")
verdict("D019", "scripts/sa2020_pv_analysis_script.R:233-248,306-309,405-408", "Animal pairing is explicit, but source run counts remain per-cell and retained is incorrectly set to 9766 cells.", "partial", "wrong study-provenance granularity/value", "PARTIAL FIX")
verdict("D020", "scripts/hazzard2022_pv_analysis_script.R:225-263; scripts/singleR.R:263-276", "singleR expects source sporozoite fields that Hazzard2022 never creates, so runs 498/499 are not overridden.", "no", "confirmed contradiction would recur/not be repaired", "INCORRECT FIX")
verdict("D021", "scripts/singleR.R:46-89", "Pre-gating intent is sound, but Hazzard2022 has no total_umi_count and becomes NA-eligible.", "no", "SingleR subset path fails before annotation", "INCORRECT FIX")
verdict("D023", "scripts/hazzard2022_pv_analysis_script.R:212-217", "Keeps registered and publication protocol instruments separately.", "yes", "none found", "VERIFIED FIX")
verdict("D024", "scripts/hazzard2022_pv_analysis_script.R:252-257", "Uses true NA; unresolved host/day remain unasserted rather than fabricated.", "partial", "authoritative population reconciliation still open", "PARTIAL FIX")
verdict("D025", "audit/shared_pipeline decision DEC08 vs Phase-1 implementation", "Required publication/accession/QC audit was not completed, yet Hazzard-specific metadata was finalized.", "no", "repairs precede required evidence closure", "INCORRECT FIX")
verdict("D026", "scripts/ref/hazzard2024_run_table.tsv", "All 24 library_id/animal_id mappings match the authoritative table.", "yes", "none found", "VERIFIED FIX")
verdict("D027", "scripts/ref/hazzard2024_run_table.tsv", "SRR27021993 is Chesson via keyed table.", "yes", "none found", "VERIFIED FIX")
verdict("D028", "scripts/ref/hazzard2024_run_table.tsv", "All 24 days match the authoritative table.", "yes", "none found", "VERIFIED FIX")
verdict("D029", "scripts/ref/hazzard2024_run_table.tsv", "Uses unqualified 10x_Chromium.", "yes", "none found", "VERIFIED FIX")
verdict("D030", "scripts/ref/hazzard2024_run_table.tsv", "Animal replicate identity and type are keyed for all retained runs.", "yes", "none found", "VERIFIED FIX")
verdict("D032", "scripts/ref/hazzard2024_run_table.tsv; chunk scripts:216-225", "Removes num_srr but repeats three study-level counts on every cell.", "no", "violates authoritative schema granularity", "INCORRECT FIX")
verdict("D033", "scripts/pipeline_lib.R:213-258; Hazzard2024 chunk scripts", "Replaces positional vectors with keyed run-table application.", "yes", "helper still accepts NA/duplicate data keys", "VERIFIED FIX")
verdict("D034", "scripts/singleR.R:46-89", "Hazzard2024 lacks total_umi_count and source PseudoGroup eligibility fields, so intended asexual gate is not operational.", "no", "atlas annotation build fails/misclassifies sexual cells", "INCORRECT FIX")
verdict("D037", "scripts/ruberto2022_1_pv_analysis_script.R:223,360-370; scripts/singleR.R:83-89", "Liver tissue and near-zero fields are present and excluded before inference.", "yes", "D036 itself remains excluded from scope", "VERIFIED FIX")
verdict("D038", "scripts/ruberto2022_1_pv_analysis_script.R:214-218", "Correct unqualified technology spelling.", "yes", "none found", "VERIFIED FIX")
verdict("D039", "scripts/ruberto2022_1_pv_analysis_script.R:226-229", "Preserves source None but harmonizes it to No_Treatment.", "yes", "none found", "VERIFIED FIX")
verdict("D040", "scripts/ruberto2022_2_pv_analysis_script.R:210-214", "Separates parasite origin, experimental site, and sequencing site.", "yes", "none found", "VERIFIED FIX")
verdict("D041", "scripts/ruberto2022_2_pv_analysis_script.R:215-218", "Correct verified V3 spelling.", "yes", "none found", "VERIFIED FIX")
verdict("D042", "scripts/ruberto2022_2_pv_analysis_script.R:226-227,250-251", "Separates field-isolate label from geography.", "yes", "none found", "VERIFIED FIX")
verdict("D043", "scripts/ruberto2022_2_pv_analysis_script.R:221-247", "Adds integer host_taxid; legacy host_id remains for compatibility.", "yes", "legacy alias remains but canonical field is correct", "VERIFIED FIX")
verdict("D044", "scripts/ruberto2022_2_pv_analysis_script.R:338-343; scripts/singleR.R:83-89", "Source-defined vector tissue is structurally ineligible for IDC.", "yes", "none found", "VERIFIED FIX")
verdict("D045", "scripts/singleR.R:399-419", "Source-selection-defined stages are excluded and asserted against gametocyte calls.", "yes", "marker path itself is broken", "VERIFIED FIX")
verdict("D046", "scripts/flatten_data.R:89-114", "New names are intended, but legacy aliases remain and no rebuilt export/app contract demonstrates reconciliation.", "partial", "schema drift remains likely", "PARTIAL FIX")
verdict("D047", "scripts/singleR.R:282-302", "Uses graph.name and asserts RNA_snn exists.", "yes", "biological effect remains Phase-2 validation", "VERIFIED FIX")
verdict("D048", "scripts/singleR.R:305-395", "Replaces numeric IDs, but the implementation has an unnamed-vector index bug, undefined writer, and unsupported fixed thresholds.", "no", "annotation script cannot complete", "INCORRECT FIX")
verdict("D049", "scripts/singleR.R:500-532", "Adds source gate and manifest intent, but writer is undefined and non-source Anopheles removal still depends on inferred class.", "partial", "latent unvalidated removal remains", "PARTIAL FIX")
verdict("D050", "scripts/integration.R:79-116", "Adds layer assertions, but indexes an unnamed metadata vector by cell IDs, producing NA study membership.", "no", "Harmony preparation fails before integration", "INCORRECT FIX")
verdict("D051", "scripts/integration.R:58-152", "Layer-derived grouping and JoinLayers timing are explicitly documented and asserted.", "yes", "mapping expression needs correction", "VERIFIED FIX")
verdict("D052", "scripts/flatten_data.R:232-283", "Constructs and saves accumulated top_genes_df with a six-study assertion.", "yes", "dense export performance remains", "VERIFIED FIX")
verdict("D053", "scripts/flatten_data.R:48-87", "Cleans metadata names before any_of and asserts survival.", "yes", "duplicate-clean-name collision is not separately rejected", "VERIFIED FIX")
verdict("D054", "scripts/flatten_data.R:31-38,89-114", "Candidate-only export is safe, but schema reconciliation is not verifiable until rebuild/app update.", "partial", "candidate schema may not satisfy current app", "PARTIAL FIX")
verdict("D055", "scripts/flatten_data.R:133-165", "Sorts and asserts exact row sets/order before cbind.", "yes", "none found", "VERIFIED FIX")
verdict("D056", "scripts/pipeline_lib.R:18-25; scripts/data_source_manipulation.R:40-50", "Shared preprocessing constants used; app consumers are not updated in this diff.", "partial", "end-to-end mismatch remains possible", "PARTIAL FIX")
verdict("D057", "scripts/data_source_manipulation.R:28-100", "Correct source code values, but no candidate-only path and no regenerated artifact is reviewed.", "partial", "running script overwrites data/data_source.csv", "PARTIAL FIX")
verdict("D058", "scripts/singleR.R:73-213", "Renames fields and pre-gates, but total_umi_count is absent for five studies and source sexual labels are imported after eligibility.", "no", "failure plus circular/late eligibility", "INCORRECT FIX")
verdict("D059", "scripts/singleR.R:96-149", "Logs dropped genes and records reference provenance/cross-platform semantics.", "yes", "none found", "VERIFIED FIX")
verdict("D060", "scripts/singleR.R:214-222", "Persists explicitly named heuristic bins and documents cutoffs as PlaViSca-defined.", "yes", "none found", "VERIFIED FIX")
verdict("D061", "scripts/singleR.R:182-208", "Persists raw/pruned-equivalent label, delta.next and pruning flag separately.", "yes", "execution blocked elsewhere", "VERIFIED FIX")
verdict("D062", "scripts/pipeline_lib.R:144-155; scripts/flatten_data.R:116-121", "Adds literal-NA assertions, but factor columns escape export checking and run-count/treatment vocabulary lacks central validation.", "partial", "literal NA factors could pass", "PARTIAL FIX")
verdict("D063", "PlaViSca app code (not changed in Phase-1 diff)", "Dead hard-coded columns were not repaired.", "no", "existing defect remains", "INCORRECT FIX")
verdict("D064", "PlaViSca app code (not changed in Phase-1 diff)", "Hard-coded label/color maps were not repaired.", "no", "existing defect remains", "INCORRECT FIX")
verdict("D065", "PlaViSca app code (not changed in Phase-1 diff)", "Unreachable blood_stage branch was not repaired.", "no", "existing defect remains", "INCORRECT FIX")
verdict("D066", "PlaViSca app code (not changed in Phase-1 diff)", "Male gametocyte/config color handling was not repaired.", "no", "existing defect remains", "INCORRECT FIX")
verdict("D067", "scripts/data_source_manipulation.R:103-117", "Validates bibliography labels against constants at build time, not the required app-startup join against mr_data.", "partial", "drift can still reach app", "PARTIAL FIX")
verdict("D068", "PlaViSca app code (not changed in Phase-1 diff)", "Biological lifecycle ordering was not implemented.", "no", "existing defect remains", "INCORRECT FIX")
verdict("D069", "all changed study scripts", "No production Chomium spelling remains.", "yes", "legacy sc_technology field name remains", "VERIFIED FIX")

claimed = []
excluded = {11, 12, 18, 22, 31, 35, 36}
for number in range(1, 70):
    if number in excluded:
        continue
    in_claim = (number <= 10 or 13 <= number <= 17 or 19 <= number <= 21 or
                23 <= number <= 30 or 32 <= number <= 34 or 37 <= number <= 69)
    if not in_claim:
        continue
    did = f"D{number:03d}"
    location, summary, matches, regression, classification = V[did]
    claimed.append({
        "defect_id": did,
        "claimed_fixed": "yes",
        "code_location": location,
        "implementation_summary": summary,
        "matches_audit_spec": matches,
        "potential_regression": regression,
        "classification": classification,
    })
write_tsv(OUT / "defect_fix_verification.tsv", claimed,
          ["defect_id", "claimed_fixed", "code_location", "implementation_summary",
           "matches_audit_spec", "potential_regression", "classification"])
