#!/usr/bin/env bash
set -euo pipefail

candidate_root=${PLAVISCA_CANDIDATE_ROOT:-$(pwd)}

sha256sum \
  "$candidate_root/pv_all_studies.rds" \
  "$candidate_root/data/candidate_export/normalize_df.rds" \
  "$candidate_root/data/candidate_export/raw_df.rds" \
  "$candidate_root/data/candidate_export/scale_df.rds" \
  "$candidate_root/data/candidate_export/cleaned_dataset.rds" \
  "$candidate_root/hazzard2022.rds" \
  "$candidate_root/hazzard2024.rds" \
  "$candidate_root/ruberto2022_1.rds" \
  "$candidate_root/ruberto2022_2.rds" \
  "$candidate_root/sa2020.rds" \
  "$candidate_root/silva2022.rds" \
  "$candidate_root/data/data_source.csv"
