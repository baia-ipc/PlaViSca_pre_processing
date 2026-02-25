#!/usr/bin/bash

# List of R scripts
scripts=(
  "hazzard2024_pv_combined_55_56.R"
  "hazzard2024_pv_combined_59_65.R"
  "hazzard2024_pv_combined_66_72.R"
  "hazzard2024_pv_combined_75_76.R"
  "hazzard2024_pv_combined_80_87.R"
  "hazzard2024_pv_combined_91_94.R"
)

# Loop through and run each one
for script in "${scripts[@]}"; do
  echo "Running $script"
  Rscript "$script"
  if [ $? -ne 0 ]; then
    echo "Error in $script, continuing to next ${script}..."
  fi
done


if [ -f "hazzard2024_merge_all.R" ]; then
    echo -e "\nrunning hazzard2024_merge_all..."
    Rscript hazzard2024_merge_all.R
else
    echo "hazzard2024_merge_all.R not found!"
    exit 1
fi

if [ -f "hazzard2022_pv_analysis_script.R" ]; then
    echo -e "\nrunning hazzard2022..."
    Rscript hazzard2022_pv_analysis_script.R
else
    echo "hazzard2022_pv_analysis_script.R not found!"
    exit 1
fi

if [ -f "ruberto2022_1_pv_analysis_script.R" ]; then
    echo -e "\nrunning ruberto2022_1..."
    Rscript ruberto2022_1_pv_analysis_script.R
else
    echo "ruberto2022_1_pv_analysis_script.R not found!"
    exit 1
fi

if [ -f "ruberto2022_2_pv_analysis_script.R" ]; then
    echo -e "\nrunning ruberto2022_2..."
    Rscript ruberto2022_2_pv_analysis_script.R

else
    echo "ruberto2022_2_pv_analysis_script.R not found!"
    exit 1
fi

if [ -f "sa2020_pv_analysis_script.R" ]; then
    echo -e "\nrunning sa2020..."
    Rscript sa2020_pv_analysis_script.R
else
    echo "sa2020_pv_analysis_script.R not found!"
    exit 1
fi


if [ -f "silva.R" ]; then
    echo -e "\nruning silva..."
    Rscript silva.R
else
    echo "silva.R not found!"
    exit 1
fi

if [ -f "integration.R" ]; then
    echo -e "\nruning integration..."
    Rscript integration.R
else
    echo "integration.R not found!"
    exit 1
fi

if [ -f "singleR.R" ]; then
    echo -e "\nruning singleR..."
    Rscript singleR.R
else
    echo "SingleR.R not found!"
    exit 1
fi

echo "All scripts run completly."






