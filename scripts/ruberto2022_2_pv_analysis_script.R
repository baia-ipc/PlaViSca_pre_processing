##### --------------------------------------------------------------------------
##### -Ruberto, CTEGD, Institute of Bioinformatics, UGA
##### --------------------------------------------------------------------------

##### --------------------------------------------------------------------------
##### 1. Load libraries
##### --------------------------------------------------------------------------
library(tidyverse)
library(Seurat)
# DropletUtils/rtracklayer/scater/SingleCellExperiment are only needed by the
# slow STARsolo-reprocessing fallback path below - loaded lazily there.

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
##### 2-5. Build per-run expression objects (all_pv).
#####
##### PHASE 2 efficiency fix (2026-09-20, at user request): expression
##### counts/cell-membership are unchanged by Phase 1/2 for this study (only
##### metadata was buggy - D040-D043,D069 below); Part II of the Phase 2 spec
##### requires non-Ruberto2022_1 studies' expression matrices to remain
##### consistent with their audited production source. Membership itself is
##### independently determined below (Section 7) from
##### data/6_PvSPZ.BS.combined.rds's Stage=="Sporozoite" selection, not from
##### STARsolo. When a prior ruberto2022_2.rds with the expected 9947-cell,
##### full 6811-feature panel already exists, its RNA counts for the union of
##### the 3 runs are reused directly (split back into a per-run all_pv list,
##### then re-subset to the same sporozoite cell list below) instead of
##### re-reading STARsolo's raw barcode matrices - bit-identical counts,
##### skips the slow I/O. Falls back to full STARsolo reprocessing (which
##### then also needs the fuller, pre-subset raw barcode space) if no
##### reusable artifact exists.
##### --------------------------------------------------------------------------
study_num <- "35926062"
samples <- paste0("ERR5087", 438:440)
run_suffixes <- str_extract(samples, "\\d{3}$")
existing_path <- "ruberto2022_2.rds"

all_pv <- NULL
if (file.exists(existing_path)) {
  reused <- tryCatch(readRDS(existing_path), error = function(e) NULL)
  if (!is.null(reused) && inherits(reused, "Seurat") && ncol(reused) == 9947L && nrow(reused) == 6811L) {
    cell_suffix <- sub("^([0-9]{3})_.*$", "\\1", colnames(reused))
    if (!anyNA(cell_suffix) && all(unique(cell_suffix) %in% run_suffixes)) {
      reused_counts <- GetAssayData(reused, assay = "RNA", layer = "counts")
      all_pv <- stats::setNames(
        lapply(seq_along(samples), function(i) {
          suf <- run_suffixes[[i]]
          cells_i <- colnames(reused)[cell_suffix == suf]
          if (length(cells_i) == 0) {
            # This run's cells may already all be present (post-subset
            # objects only ever shrink, never gain runs), so an empty
            # per-run slice here is unexpected - fall back rather than
            # silently building a run with zero cells.
            return(NULL)
          }
          mat_i <- reused_counts[, cells_i, drop = FALSE]
          colnames(mat_i) <- sub(paste0("^", suf, "_"), "", cells_i)
          CreateSeuratObject(counts = mat_i, min.cells = 0, min.features = 0)
        }),
        samples
      )
      if (any(vapply(all_pv, is.null, logical(1)))) all_pv <- NULL
    }
  }
}
if (!is.null(all_pv)) {
  message(sprintf(
    "ruberto2022_2: reused expression counts from existing %s (STARsolo re-read skipped); recomputing metadata from current script logic only.",
    existing_path
  ))
}

if (is.null(all_pv)) {
  message("ruberto2022_2: no reusable prior ruberto2022_2.rds found - running full STARsolo reprocessing (slow path).")
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
    path <- file.path(paste0("./counts/", study_num, "/", sample_id, "_solo_out/Solo.out/GeneFull/raw"))
    sce <- read10xCounts(path, col.names = TRUE)
    all_pv[[sample_id]] <- sce
  }

  for (i in 1:length(all_pv)) {
    all_pv[[i]] <- all_pv[[i]][!rownames(all_pv[[i]]) %in% rRNA_pv$ID, ]
  }

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

  # set.seed(123456)
  # all_pv <- lapply(all_pv, function(x) {
  #   e.out <- emptyDrops(x, lower = 40) #change this value if needed
  #   x <- x[, which(e.out$FDR <= 0.001)]
  # })

  for (i in 1:length(all_pv)) {
    sce <- all_pv[[i]]
    mat <- counts(sce)
    all_pv[[i]] <- CreateSeuratObject(counts = mat, min.cells = 0, min.features = 0)
  }
}

# Add metadata via explicit per-run vectors, each keyed 1:1 to the three
# ERR runs/isolates - fixes D040-D043,D069.
run_id <- paste0("ERR5087", 438:440)
assert_cardinality(run_id, length(all_pv), label = "run_id")
study_label <- rep(STUDY_LABELS[["ruberto2022_2"]], length(all_pv))
pub_year <- rep(2022, length(all_pv))
# D040: parasite origin / experimental site / sequencing site kept separate,
# never a single combined "Cambodia" field.
parasite_origin_location <- rep("Cambodia_Mondulkiri Province", length(all_pv))
experimental_site <- rep("Institut Pasteur du Cambodge (Phnom Penh)", length(all_pv))
sequencing_site <- rep("Macrogen (Seoul)", length(all_pv))
# D041,D069: corrected spelling; V3 chemistry IS independently verified for
# this study, so it is retained (unlike the other studies' unqualified
# "10x_Chromium").
sc_technology <- rep("10x_Chromium_V3", length(all_pv))
sequencer_registered <- rep("HiSeq X Ten", length(all_pv))
host_species <- rep("Anopheles dirus", length(all_pv))
# D043: proper integer NCBI taxid field, not an ad hoc "tax_id-7168" string
# (the value itself, Anopheles dirus 7168, was already correct).
host_taxid <- rep(7168L, length(all_pv))
sample_type <- rep("Vector host: salivary gland", length(all_pv))
tissue_or_sample_type <- rep("Vector salivary gland", length(all_pv))
# D042: strain/field-isolate terminology separated from geography.
strain <- rep("Field isolate", length(all_pv))
day_post_infection <- c(16, 18, 17)
biological_replicate_id <- c("1", "2", "3")
biological_replicate_type <- rep("infection", length(all_pv))

for (i in 1:length(all_pv)) {
  n_cells <- ncol(all_pv[[i]])
  assert_cardinality(colnames(all_pv[[i]]), n_cells, label = paste0("colnames(", run_id[[i]], ")"))

  all_pv[[i]]$study_pmid <- rep(study_num, n_cells)
  all_pv[[i]]$run_id <- rep(run_id[[i]], n_cells)
  all_pv[[i]]$study_label <- rep(study_label[[i]], n_cells)
  all_pv[[i]]$pub_year <- rep(pub_year[[i]], n_cells)
  all_pv[[i]]$parasite_origin_location <- rep(parasite_origin_location[[i]], n_cells)
  all_pv[[i]]$experimental_site <- rep(experimental_site[[i]], n_cells)
  all_pv[[i]]$sequencing_site <- rep(sequencing_site[[i]], n_cells)
  all_pv[[i]]$sc_technology <- rep(sc_technology[[i]], n_cells)
  all_pv[[i]]$sequencer_registered <- rep(sequencer_registered[[i]], n_cells)
  all_pv[[i]]$host_species <- rep(host_species[[i]], n_cells)
  all_pv[[i]]$host_taxid <- rep(host_taxid[[i]], n_cells)
  all_pv[[i]]$host_id <- rep(paste0("tax_id-", host_taxid[[i]]), n_cells) # retained for backward compatibility with existing app schema
  all_pv[[i]]$sample_type <- rep(sample_type[[i]], n_cells)
  all_pv[[i]]$tissue_or_sample_type <- rep(tissue_or_sample_type[[i]], n_cells)
  all_pv[[i]]$strain <- rep(strain[[i]], n_cells)
  all_pv[[i]]$parasite_lineage_isolate <- rep(strain[[i]], n_cells)
  all_pv[[i]]$day_post_infection <- rep(day_post_infection[[i]], n_cells)
  all_pv[[i]]$treatment <- rep("No_Treatment", n_cells)
  all_pv[[i]]$biological_replicate_id <- rep(biological_replicate_id[[i]], n_cells)
  all_pv[[i]]$biological_replicate_type <- rep(biological_replicate_type[[i]], n_cells)

  all_pv[[i]]$barcode <- paste(
    str_extract(all_pv[[i]]$run_id, "\\d{3}$"),
    "_",
    colnames(all_pv[[i]]),
    sep = ""
  )
  colnames(all_pv[[i]]) <- all_pv[[i]]$barcode
}

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

# for (i in 1:length(all_pv)) {
#   all_pv[[i]] <- NormalizeData(all_pv[[i]], verbose = FALSE)
#   all_pv[[i]] <- FindVariableFeatures(all_pv[[i]], selection.method = "vst",
#                                         nfeatures = nrow(all_pv[[i]])*.3,
#                                       verbose = FALSE)
# }

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


# read data from ruberto2022_2 study
ruberto2022_2 <- readRDS("data/6_PvSPZ.BS.combined.rds")
spz_cell_ids <- WhichCells(ruberto2022_2, expression = Stage == "Sporozoite")

cells_to_subset <- data.frame(cells = spz_cell_ids) |>
  mutate(
    cells_name = sub("_.*", "", cells),
    cells_num = sub("^[^_]+_", "", cells),
    new_cells = case_when(
      cells_num == "Case_909" ~ "438",
      cells_num == "Case_922" ~ "439",
      cells_num == "Case_923" ~ "440"
    ),
    new_cells_name = paste0(new_cells, "_", cells_name)
  ) |>
  select(cells = new_cells_name)

pv.combined.all <- subset(pv.combined.all, cells = cells_to_subset$cells)

# Source-selection-defined Sporozoite provenance (D044,D045,D046): every
# cell in this study was selected because the source object's own Stage ==
# "Sporozoite" - this is authoritative and must never be silently overridden
# downstream by an incompatible blood-stage classifier (see singleR.R).
pv.combined.all$source_life_cycle_stage <- "Sporozoite"
pv.combined.all$source_stage_provenance <- "source_selection_defined"

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

# Conservative Phase 1 membership policy: preserve the exact 9,947-cell
# population; do not modify counts or membership.
assert_cell_count(ncol(pv.combined.all), 9947L, label = "ruberto2022_2.rds")
assert_unique_cell_ids(colnames(pv.combined.all), label = "ruberto2022_2.rds colnames")

# Save Seurat object
saveRDS(pv.combined.all, file = "ruberto2022_2.rds")
record_build_manifest(
  artifact_path = "ruberto2022_2.rds",
  script_path = "scripts/ruberto2022_2_pv_analysis_script.R",
  cell_count = ncol(pv.combined.all),
  notes = "Fixes D040-D043,D069; D044/D045/D046 HPI-gating and override-scoping fixes live in singleR.R"
)
# Remove all object
rm(list = ls())
gc()
