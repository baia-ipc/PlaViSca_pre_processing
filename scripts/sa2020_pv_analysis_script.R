##### --------------------------------------------------------------------------
##### -Ruberto, CTEGD, Institute of Bioinformatics, UGA
##### --------------------------------------------------------------------------

##### --------------------------------------------------------------------------
##### 1. Load libraries
##### --------------------------------------------------------------------------
library(tidyverse)
library(Seurat)
# DropletUtils/rtracklayer/scater/SingleCellExperiment are only needed by the
# slow STARsolo-reprocessing fallback path (Section 3b below) - loaded lazily
# there via requireNamespace() so a fast-path run (the common case on a
# machine that already has sa2020.rds) does not pay their load cost either.

# set current working directory
.plavisca_root <- local({
  candidates <- unique(c(
    Sys.getenv("PLAVISCA_PREPROCESS_ROOT", unset = NA_character_),
    getwd(),
    "/home/sopheap/pvsca_b/pre_process_data"
  ))
  candidates <- candidates[!is.na(candidates) & nzchar(candidates)]
  hit <- candidates[file.exists(file.path(candidates, "scripts", "pipeline_lib.R"))]
  if (length(hit) == 0) {
    stop(
      "Could not locate the pre_process_data project root (looked for scripts/",
      "pipeline_lib.R under $PLAVISCA_PREPROCESS_ROOT, the current working ",
      "directory, and the historical hard-coded path). Set the ",
      "PLAVISCA_PREPROCESS_ROOT environment variable to the pre_process_data ",
      "directory's absolute path, or run this script from that directory.",
      call. = FALSE
    )
  }
  normalizePath(hit[[1]])
})
setwd(.plavisca_root)
source("scripts/pipeline_lib.R")

##### --------------------------------------------------------------------------
##### 2-5. Build per-run expression objects (all_pv), one Seurat object per
##### run, colnames = raw 10x barcode (not yet prefixed).
#####
##### PHASE 2 efficiency fix (2026-09-20, at user request): expression
##### counts/cell-membership for this study are NOT changed by Phase 1/2
##### (only metadata-assignment logic was buggy - see D016,D017,D019,D069
##### below); Part II of the Phase 2 spec requires non-Ruberto2022_1 studies'
##### expression matrices to remain consistent with their audited production
##### source. So when a prior sa2020.rds with the expected 9766-cell,
##### `<3-digit-run-suffix>_<barcode>` colname structure already exists, its
##### RNA counts are reused directly (split back into a per-run all_pv list)
##### instead of re-reading STARsolo's raw/filtered barcode matrices - this
##### produces bit-identical counts to a full reprocessing run (the same
##### matrix, not a recomputation) while skipping the slow read10xCounts/
##### barcodeRanks/emptyDrops I/O entirely. A from-scratch checkout with no
##### prior artifact falls back to the original full STARsolo pipeline.
##### --------------------------------------------------------------------------
study_num <- "32365102"
samples <- paste0("SRR11008", 269:278)
run_suffixes <- str_extract(samples, "\\d{3}$")
existing_path <- "sa2020.rds"

all_pv <- NULL
if (file.exists(existing_path)) {
  reused <- tryCatch(readRDS(existing_path), error = function(e) NULL)
  if (!is.null(reused) && inherits(reused, "Seurat") && ncol(reused) == 9766L) {
    cell_suffix <- sub("^([0-9]{3})_.*$", "\\1", colnames(reused))
    if (!anyNA(cell_suffix) && setequal(unique(cell_suffix), run_suffixes)) {
      reused_counts <- GetAssayData(reused, assay = "RNA", layer = "counts")
      all_pv <- stats::setNames(
        lapply(seq_along(samples), function(i) {
          suf <- run_suffixes[[i]]
          cells_i <- colnames(reused)[cell_suffix == suf]
          mat_i <- reused_counts[, cells_i, drop = FALSE]
          colnames(mat_i) <- sub(paste0("^", suf, "_"), "", cells_i)
          CreateSeuratObject(counts = mat_i, min.cells = 0, min.features = 0)
        }),
        samples
      )
      message(sprintf(
        "sa2020: reused expression counts from existing %s (STARsolo re-read skipped); recomputing metadata from current script logic only.",
        existing_path
      ))
    }
  }
}

if (is.null(all_pv)) {
  message("sa2020: no reusable prior sa2020.rds found - running full STARsolo reprocessing (slow path).")
  library(DropletUtils)
  library(rtracklayer)
  library(scater)
  library(SingleCellExperiment)

  pv.gff <- import.gff3("./scripts/ref/PlasmoDB-68_PvivaxP01.gff")
  pv.gff <- as.data.frame(pv.gff)
  rRNA_pv <- pv.gff %>%
    as_tibble() %>%
    filter(type == "rRNA") %>%
    select(ID, description) %>%
    mutate(ID = gsub("\\.1", "", ID))

  all_pv <- list()
  for (sample_id in samples) {
    path <- file.path(
      paste0("./counts/", study_num, "/", sample_id, "_solo_out/Solo.out/GeneFull/filtered")
    )
    sce <- read10xCounts(path, col.names = TRUE)
    all_pv[[sample_id]] <- sce
  }

  # Since we are interested in assessing changes in protein-coding
  # transcripts, we can remove rRNA from the dataset.
  for (i in 1:length(all_pv)) {
    all_pv[[i]] <- all_pv[[i]][!rownames(all_pv[[i]]) %in% rRNA_pv$ID, ]
  }

  # Knee plots (QC visualization only, not used for droplet-calling - see
  # commented-out emptyDrops() below, unchanged from the original pipeline).
  pdf(paste0(study_num, "_knee_plots.pdf"), width = 12, height = 6)
  par(mfrow = c(1, 2))
  bcrank <- list()
  uniq <- list()
  for (i in seq_along(all_pv)) {
    bcrank[[i]] <- barcodeRanks(all_pv[[i]])
    uniq[[i]] <- !duplicated(bcrank[[i]]$rank)
    plot(
      bcrank[[i]]$rank[uniq[[i]]], bcrank[[i]]$total[uniq[[i]]],
      log = "xy", xlab = "Rank", ylab = "Total UMI Count", cex.lab = 1.2,
      main = names(all_pv)[i]
    )
    abline(h = metadata(bcrank[[i]])$inflection, col = "darkgreen", lty = 2)
    abline(h = metadata(bcrank[[i]])$knee, col = "dodgerblue", lty = 2)
    legend("bottomleft", legend = c("Inflection", "Knee"), col = c("darkgreen", "dodgerblue"), lty = 2, cex = 1.2)
    if (i %% 2 == 0 && i < length(all_pv)) par(mfrow = c(1, 2))
  }
  dev.off()
  rm(bcrank, uniq)

  # Another option is to use the emptyDrops.
  # emptyDrops performs Monte Carlo simulations to compute p-values, so we
  # need to set the seed to obtain reproducible results.
  # see PMID: 30902100 for rationale and statistical framework underlying
  # this method
  # set.seed(123456)
  # all_pv <- lapply(all_pv, function(x) {
  #   e.out <- emptyDrops(x, lower = 40) #change this value if needed
  #   x <- x[, which(e.out$FDR <= 0.001)]
  # })

  for (i in 1:length(all_pv)) {
    sce <- all_pv[[i]]
    mat <- counts(sce)
    all_pv[[i]] <- CreateSeuratObject(counts = mat, min.cells = 1, min.features = 1)
  }
}

# Add metadata via explicit per-run vectors (each already keyed 1:1 to
# `samples`, never recycled) plus keyed post-hoc fixes for D016,D017,D019.
run_id <- paste0("SRR11008", 269:278)
assert_cardinality(run_id, length(all_pv), label = "run_id")
study_label <- rep(STUDY_LABELS[["sa2020"]], length(all_pv))
pub_year <- rep(2020, length(all_pv))
# D017: sample/registry location (deposit site), kept separate from strain
# parasite-origin below - never combined into one "geographic_location".
registry_location <- rep("USA_Rockville", length(all_pv))
experimental_site <- rep("NIH (Rockville, MD)", length(all_pv))
# D016,D069: corrected spelling; chemistry version (V2) not independently
# confirmed, so left unqualified.
sc_technology <- rep("10x_Chromium", length(all_pv))
sequencer <- rep("Illumina HiSeq 4000", length(all_pv))
host_species <- c(
  "Saimiri boliviensis",
  "Aotus nancymaae",
  "Aotus nancymaae",
  "Saimiri boliviensis",
  "Saimiri boliviensis",
  "Aotus nancymaae",
  "Aotus nancymaae",
  "Saimiri boliviensis",
  "Aotus nancymaae",
  "Aotus nancymaae"
)
# D019: host_id already uniquely identifies the individual animal (three IDs
# repeat: 86574, 86436, 86416, each pairing one chloroquine-treated run with
# its untreated same-animal counterpart) - this is exactly the replicate
# structure the schema's biological_replicate_id/animal_id fields require.
animal_id <- c(
  "3879",
  "86574",
  "86574",
  "5541",
  "4215",
  "86436",
  "86436",
  "5107",
  "86416",
  "86416"
)
sample_type <- rep("Mammalian host: blood", length(all_pv))
tissue_or_sample_type <- rep("Host blood", length(all_pv))
strain <- c(
  "NIH-1993",
  "Indonesia-I",
  "Indonesia-I",
  "NIH-1993",
  "Chesson",
  "Chesson",
  "Chesson",
  "AMRU-I",
  "AMRU-I",
  "AMRU-I"
)
# D017: parasite strain-origin, kept as a separate field from the
# registry/experimental deposit location above.
parasite_origin_location <- c(
  "unresolved (NIH-1993; closely related to but distinct from Salvador-I)",
  "Indonesia",
  "Indonesia",
  "unresolved (NIH-1993; closely related to but distinct from Salvador-I)",
  "New Guinea",
  "New Guinea",
  "New Guinea",
  "Papua New Guinea",
  "Papua New Guinea",
  "Papua New Guinea"
)
day_post_infection <- c(15, 16, 12, 21, 24, 16, 12, 24, 16, 12)
source_treatment <- c(
  "None",
  "CQ_(16h_5mg/kg)",
  "None",
  "None",
  "None",
  "CQ_(16h_10mg/kg)",
  "None",
  "None",
  "CQ_(16h_10mg/kg)",
  "None"
)
treatment <- ifelse(source_treatment == "None", "No_Treatment", source_treatment)

for (i in seq_along(all_pv)) {
  n_cells <- ncol(all_pv[[i]])
  assert_cardinality(colnames(all_pv[[i]]), n_cells, label = paste0("colnames(", names(all_pv)[i], ")"))

  all_pv[[i]]$study_pmid <- rep(study_num, n_cells)
  all_pv[[i]]$run_id <- rep(run_id[[i]], n_cells)
  all_pv[[i]]$study_label <- rep(study_label[[i]], n_cells)
  all_pv[[i]]$pub_year <- rep(pub_year[[i]], n_cells)
  all_pv[[i]]$registry_location <- rep(registry_location[[i]], n_cells)
  all_pv[[i]]$experimental_site <- rep(experimental_site[[i]], n_cells)
  all_pv[[i]]$sequencing_site <- rep(NA_character_, n_cells)
  all_pv[[i]]$sc_technology <- rep(sc_technology[[i]], n_cells)
  all_pv[[i]]$sequencer <- rep(sequencer[[i]], n_cells)
  all_pv[[i]]$host_species <- rep(host_species[[i]], n_cells)
  all_pv[[i]]$host_taxid <- rep(NA_integer_, n_cells)
  all_pv[[i]]$animal_id <- rep(animal_id[[i]], n_cells)
  all_pv[[i]]$biological_replicate_id <- rep(animal_id[[i]], n_cells)
  all_pv[[i]]$biological_replicate_type <- rep("animal", n_cells)
  all_pv[[i]]$sample_type <- rep(sample_type[[i]], n_cells)
  all_pv[[i]]$tissue_or_sample_type <- rep(tissue_or_sample_type[[i]], n_cells)
  all_pv[[i]]$strain <- rep(strain[[i]], n_cells)
  all_pv[[i]]$parasite_lineage_isolate <- rep(strain[[i]], n_cells)
  all_pv[[i]]$parasite_origin_location <- rep(parasite_origin_location[[i]], n_cells)
  all_pv[[i]]$day_post_infection <- rep(day_post_infection[[i]], n_cells)
  all_pv[[i]]$source_treatment <- rep(source_treatment[[i]], n_cells)
  all_pv[[i]]$treatment <- rep(treatment[[i]], n_cells)

  all_pv[[i]]$barcode <- paste(
    str_extract(all_pv[[i]]$run_id, "\\d{3}$"),
    "_",
    colnames(all_pv[[i]]),
    sep = ""
  )
  colnames(all_pv[[i]]) <- all_pv[[i]]$barcode
}

# D013-style fix (num_srr): study-level provenance, never an ambiguous
# per-cell scalar (D019 explicitly forbids a "num_srr cell field" for Sa2020).
sa2020_source_run_count_deposited <- 10L
sa2020_source_run_count_imported <- 10L

##### --------------------------------------------------------------------------
##### 6. Processing step 3: data normalization, variable selection, scaling
##### and dimension reduction
##### --------------------------------------------------------------------------
# The next step is to normalize the data. At the moment, the only information
# in the Seurat object is 'count' data. In order to facilitate downstream
# analysis, we will want to normalize the data. In Seurat, the default is
# By default,a global-scaling normalization method “LogNormalize” is employed.
# This step normalizes the feature expression measurements for each cell by the
# total expression, multiplies this by a scale factor (10,000 by default), and
# log-transforms the result.
# In Seurat v5, normalized values are stored in <name_of_seurat_object>[["RNA"]]$data
# Next, we will calculate a subset of features that exhibit high cell-to-cell
# variation in the dataset (i.e, they are highly expressed in some cells, and
# lowly expressed in others). Focusing on these genes in downstream analysis
# helps to highlight biological signal in single-cell datasets (PMID: 24056876

#for (i in 1:length(all_pv)) {
# all_pv[[i]] <- NormalizeData(all_pv[[i]], verbose = FALSE)
#all_pv[[i]] <- FindVariableFeatures(all_pv[[i]], selection.method = "vst",
#                                      nfeatures = nrow(all_pv[[i]])*.3,verbose = FALSE)
#}

##### --------------------------------------------------------------------------
##### 7. Processing step 4: integration, re-normalization, scaling and data red
##### --------------------------------------------------------------------------

# Next, we will inegtate the data and apply a linear transformation (‘scaling’) that is a standard
# pre-processing step prior to dimensional reduction techniques like PCA.
# The ScaleData() function:
#  - Shifts the expression of each gene, so that the mean expression across cells is 0
#  - Scales the expression of each gene, so that the variance across cells is 1
# This step gives equal weight in downstream analyses, so that highly-expressed
# genes do not dominate
# The results of this are stored in <name_of_seurat_object>[["RNA"]]$scale.data
# By default, only variable features are scaled. You can specify the features
# argument to scale additional features. Since the # of transcripts detected is
# relatively low, we will scale the entire dataset.
# Next we perform PCA on the scaled data. By default, only the previously
# determined variable features are used as input, but can be defined using
# features argument if you wish to choose a different subset (if you do want to
# use a custom subset of features, make sure you pass these to ScaleData first).
# For the first principal components, Seurat outputs a list of genes with the
# most positive and negative loadings, representing modules of genes that
# exhibit either correlation (or anti-correlation) across single-cells in the
# dataset.

pv.combined.all <- Reduce(function(x, y) merge(x, y = y), all_pv)

pv.combined.all <- JoinLayers(pv.combined.all)

# pv.combined.all  <- NormalizeData(pv.combined.all)
# pv.combined.all  <- FindVariableFeatures(pv.combined.all, selection.method = "vst",
# 										nfeatures = nrow(pv.combined.all) * 0.3, verbose = FALSE)
# pv.combined.all  <- ScaleData(pv.combined.all)
# pv.combined.all  <- RunPCA(pv.combined.all)

# D018/DEC03: retain all 9766 cells (option b - conservative Phase 1 policy,
# no cell-membership change), and add main_analysis_member as a provenance
# covariate from the existing audited publication mapping (9018 TRUE / 748
# reprocessing-only / 197 publication-main-analysis cells absent from
# PlaViSca and therefore untaggable here).
main_analysis_crosswalk <- read.delim(
  "scripts/ref/sa2020_main_analysis_cell_lineage.tsv",
  stringsAsFactors = FALSE
)
main_analysis_ids <- main_analysis_crosswalk$plavisca_cell_id[
  main_analysis_crosswalk$plavisca_present
]
pv.combined.all$main_analysis_member <- colnames(pv.combined.all) %in% main_analysis_ids

# D013/D032-style fix: run counts are study-level provenance, never a
# per-cell scalar (D019 explicitly forbids a "num_srr" cell field).
pv.combined.all$source_run_count_deposited <- sa2020_source_run_count_deposited
pv.combined.all$source_run_count_imported <- sa2020_source_run_count_imported
pv.combined.all$source_run_count_retained <- ncol(pv.combined.all)

# Explicit count provenance (Phase 2 spec Part I.4).
# QC UMI fields required by singleR.R's IDC eligibility gate (Part IV) -
# must be populated for every study, not only Ruberto2022_1.
pv.combined.all <- set_qc_umi_fields(pv.combined.all)

pv.combined.all <- set_count_provenance(
  pv.combined.all,
  source_of_counts = "plavisca_starsolo_raw_umi",
  counting_pipeline = "STARsolo",
  count_reference_version = "PlasmoDB-68_PvivaxP01",
  count_provenance_status = "plavisca_production_pipeline"
)

# Conservative Phase 1 membership policy: preserve all 9766 cells.
assert_cell_count(ncol(pv.combined.all), 9766L, label = "sa2020.rds")
assert_unique_cell_ids(colnames(pv.combined.all), label = "sa2020.rds colnames")
assert_no_literal_na_string(pv.combined.all$treatment, label = "sa2020$treatment")
assert_no_literal_na_string(pv.combined.all$day_post_infection, label = "sa2020$day_post_infection")

# Save Seurat object
saveRDS(pv.combined.all, file = "sa2020.rds")
record_build_manifest(
  artifact_path = "sa2020.rds",
  script_path = "scripts/sa2020_pv_analysis_script.R",
  cell_count = ncol(pv.combined.all),
  notes = "Fixes D016,D017,D019,D069; D018 resolved as DEC03 option (b) - retain all 9766, add main_analysis_member covariate, no membership change; D014/D015 gametocyte-import fixes live in singleR.R"
)
# Remove all object
rm(list = ls())
gc()
