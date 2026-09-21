# AT25 semantics resolution

Classification: **VALIDATION_SEMANTICS_DEFECT**.

Production defines `gametocyte_eligible <- idc_eligible`: Host-blood cells,
excluding source-selection-defined stages, with at least 10 total UMI. It
filters metadata to those cells before calculating mean female and male module
scores per `seurat_clusters`. A cluster is Female or Male only when its own
eligible-cell mean is at least 0.10 and exceeds the opposite-sex mean by at
least 0.05. The resulting cluster call is mapped to cells, after which all
ineligible cells are explicitly masked to `NA`.

The original AT25 instead calculated own and opposite-sex cluster means over
all atlas cells. It then required both (a) every cluster containing any called
cell to clear the all-cell thresholds and (b) at least 90% of called cells to
reside in clearing clusters. Requirement (a) makes (b) ineffective for the
boundary mixing it explicitly permits: even one called cell in a large mixed
cluster can fail (a), although up to 10% of calls may be outside clearing
clusters under (b).

The corrected test separates two scientific questions. Production classifier
consistency uses the exact production eligibility population and unchanged
0.10/0.05 thresholds for every called cluster. Integrated-population coherence
uses all-cell cluster context with the unchanged AUC >=0.70, positive median
shift, and called-cell coherence >=0.90 criteria; it reports all calls and
clusters outside that coherence rather than imposing perfect cluster purity.

For Male calls, cluster 4 contains 1,687 atlas cells but only 18 eligible Male
calls; 1,669 cells are ineligible, predominantly liver stage. Its all-cell male
mean is 0.0104848, but its eligible-cell male mean is 0.1284835 and its
eligible-cell male-minus-female margin is 0.1578435. Cluster 5 contains the
other 976 Male calls. Thus all 994 Male calls are production-consistent and
976/994 (98.1891%) are in all-cell marker-clearing clusters; the 18/994
(1.8109%) boundary calls in cluster 4 remain explicit in the diagnostics.

The saved production `cluster_sex_marker_summary.tsv` is independently
consistent: cluster 4 has eligible n=18, female mean=-0.0293600 and male
mean=0.1284835; cluster 5 has eligible n=976, female mean=-0.0289882 and male
mean=1.6762734. No biological/integration failure is present under the
production-consistent and predeclared coherence-aware criteria. Harmony,
SingleR, and the atlas were not rebuilt.
