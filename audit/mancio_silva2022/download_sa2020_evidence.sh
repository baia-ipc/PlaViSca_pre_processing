#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"

# Audit-only authoritative snapshots for Sa et al. 2020 (PMCID PMC7224573).
curl -L --fail --silent --show-error \
  https://www.ebi.ac.uk/europepmc/webservices/rest/PMC7224573/fullTextXML \
  -o sa2020_article.xml
curl -L --fail --silent --show-error \
  'https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA603327&result=read_run&fields=run_accession,experiment_accession,sample_accession,study_accession,secondary_study_accession,sample_alias,experiment_alias,instrument_model,library_strategy,library_source,library_selection,scientific_name&format=tsv&download=true' \
  -o sa2020_ena_runs.tsv
curl -L --fail --silent --show-error \
  'https://journals.plos.org/plosbiology/article/file?type=supplementary&id=10.1371/journal.pbio.3000711.s030' \
  -o sa2020_s030.xlsx
curl -L --fail --silent --show-error \
  'https://journals.plos.org/plosbiology/article/file?type=supplementary&id=10.1371/journal.pbio.3000711.s034' \
  -o sa2020_s034.xlsx
curl -L --fail --silent --show-error \
  'https://journals.plos.org/plosbiology/article/file?type=supplementary&id=10.1371/journal.pbio.3000711.s036' \
  -o sa2020_s036.xlsx
sha256sum sa2020_article.xml sa2020_ena_runs.tsv sa2020_s030.xlsx sa2020_s034.xlsx sa2020_s036.xlsx \
  > sa2020_authoritative_sha256.txt
