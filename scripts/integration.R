library(DropletUtils) # scRNA-seq empty droplet detection
library(Seurat) # contains scRNA-seq analysis functions

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

# Phase 1 (production repair) explicitly forbids the expensive atlas-wide
# Harmony/UMAP/clustering rebuild - this script's job in Phase 1 is to be
# correct and testable, not to execute. Phase 2 sets this env var once all
# per-study cell-membership/expression decisions (D011,D012,D018,D022,D036)
# are frozen.
if (!isTRUE(as.logical(Sys.getenv("PLAVISCA_ALLOW_FULL_INTEGRATION_REBUILD", "FALSE")))) {
  pipeline_fail(paste(
    "Full atlas-wide Harmony integration rebuild is a Phase 2 operation.",
    "Set PLAVISCA_ALLOW_FULL_INTEGRATION_REBUILD=TRUE explicitly to run this",
    "script. See audit/shared_pipeline/AUDIT.md Part 13 (staged repair",
    "strategy) - Ruberto2022_1's D036 expression defect must be resolved",
    "first."
  ))
}

# read all Seurat object
files <- c(
    "hazzard2022.rds",
    "hazzard2024.rds",
    "ruberto2022_1.rds",
    "ruberto2022_2.rds",
    "sa2020.rds",
    "silva2022.rds"
)

pv.combined.all <- lapply(files, readRDS)

# ============================================================================
# Phase 2 Part 6: integration feature universe. Capture each study's own
# feature panel BEFORE merge() zero-pads everything to the union - merging
# first and inspecting rownames(pv.combined.all) afterward would make every
# study look "measured" for every gene, hiding exactly the reference-
# incompatibility this step exists to avoid encoding into PCA/Harmony.
# ============================================================================
study_feature_sets <- lapply(pv.combined.all, rownames)
names(study_feature_sets) <- files

atlas_feature_universe <- Reduce(union, study_feature_sets)
integration_features <- Reduce(intersect, study_feature_sets)

cat(sprintf(
  "[integration.R] Feature universe: %d total (union across 6 studies), %d integration-eligible (intersection, measured in all 6 studies), %d lost to reference incompatibility (union - intersection)\n",
  length(atlas_feature_universe), length(integration_features),
  length(atlas_feature_universe) - length(integration_features)
))

dir.create("audit/phase2_rebuild", showWarnings = FALSE, recursive = TRUE)
feature_universe_tbl <- data.frame(
  feature = atlas_feature_universe,
  stringsAsFactors = FALSE
)
for (f in files) {
  feature_universe_tbl[[paste0("measured_in__", sub("\\.rds$", "", f))]] <-
    feature_universe_tbl$feature %in% study_feature_sets[[f]]
}
feature_universe_tbl$n_studies_measured <- rowSums(
  feature_universe_tbl[, paste0("measured_in__", sub("\\.rds$", "", files)), drop = FALSE]
)
feature_universe_tbl$integration_eligible <- feature_universe_tbl$feature %in% integration_features
write.table(
  feature_universe_tbl, "audit/phase2_rebuild/feature_universe.tsv",
  sep = "\t", quote = FALSE, row.names = FALSE
)

# Feature-availability table (Part 6/24): distinguishes "zero expression"
# from "gene unavailable/not quantified in that study" at study level - the
# per-cell/per-study-panel version of the same distinction
# phase2_feature_coverage.tsv already makes for Ruberto2022_1 specifically.
feature_availability_tbl <- do.call(rbind, lapply(files, function(f) {
  data.frame(
    study_file = f,
    feature = atlas_feature_universe,
    availability = ifelse(
      atlas_feature_universe %in% study_feature_sets[[f]],
      "measured", "not_measured_reference_incompatible"
    ),
    stringsAsFactors = FALSE
  )
}))
write.table(
  feature_availability_tbl, "audit/phase2_rebuild/feature_availability_by_study.tsv",
  sep = "\t", quote = FALSE, row.names = FALSE
)

# merge all layers
pv.combined.all <- Reduce(function(x, y) merge(x, y = y), pv.combined.all)

# ============================================================================
# D050/D051: explicit, testable integration-grouping diagnostics.
#
# Direct inspection of the installed Seurat 5.3.0 / SeuratObject 5.2.0
# source (audit/shared_pipeline/AUDIT.md Part 7) confirmed that neither
# IntegrateLayers() nor HarmonyIntegration() has a `group_by` formal
# parameter - a `group_by="study_label"` argument is silently absorbed into
# `...` and never read. The Harmony batch-grouping variable Seurat/Harmony
# actually use is derived by CreateIntegrationGroups() from the per-assay
# DATA LAYER structure at IntegrateLayers() time (one integration group per
# layer), not from any metadata column.
#
# The repaired mechanism therefore uses the layer mechanism deliberately and
# explicitly - never `group_by`, which is removed below since leaving it in
# would misleadingly imply it does something - and PROVES, before every
# Harmony run, that "one layer" genuinely means "one study": exact layer
# names, exact cell counts per layer, an explicit layer->study mapping, and
# zero ambiguous/mixed-study layers. This must be re-verified on every rerun
# (e.g. if `JoinLayers()` timing or merge order ever changes), never assumed
# from today's coincidental 6-layers-for-6-studies result.
# ============================================================================
rna_layers <- Layers(pv.combined.all, assay = "RNA", search = "counts")
cat(sprintf("[integration.R] RNA assay layers before JoinLayers: %s\n", paste(rna_layers, collapse = ", ")))

if (length(rna_layers) != length(files)) {
  pipeline_fail(sprintf(
    "Expected exactly %d RNA count layers (one per study), found %d (%s). The Harmony grouping mechanism (CreateIntegrationGroups, layer-derived) requires this invariant to hold.",
    length(files), length(rna_layers), paste(rna_layers, collapse = ", ")
  ))
}

layer_study_map <- lapply(rna_layers, function(layer_name) {
  cells_in_layer <- colnames(LayerData(pv.combined.all, layer = layer_name, assay = "RNA"))
  studies_in_layer <- unique(pv.combined.all$study_label[cells_in_layer])
  list(layer = layer_name, n_cells = length(cells_in_layer), studies = studies_in_layer)
})

for (entry in layer_study_map) {
  cat(sprintf(
    "[integration.R] layer '%s': %d cells, study/studies = %s\n",
    entry$layer, entry$n_cells, paste(entry$studies, collapse = ", ")
  ))
  if (length(entry$studies) != 1) {
    pipeline_fail(sprintf(
      "Ambiguous/mixed-study layer '%s' spans %d studies (%s) - the Harmony grouping mechanism requires exactly one study per layer. Do not proceed with IntegrateLayers().",
      entry$layer, length(entry$studies), paste(entry$studies, collapse = ", ")
    ))
  }
}

all_studies_represented <- sort(unique(vapply(layer_study_map, function(e) e$studies, character(1))))
expected_studies <- sort(unname(STUDY_LABELS))
if (!identical(all_studies_represented, expected_studies)) {
  pipeline_fail(sprintf(
    "Layer->study mapping does not cover exactly the 6 expected studies. Found: %s. Expected: %s.",
    paste(all_studies_represented, collapse = ", "), paste(expected_studies, collapse = ", ")
  ))
}
cat("[integration.R] Verified: one integration group per study, zero ambiguous/mixed-study layers, all 6 studies represented.\n")

pv.combined.all <- NormalizeData(pv.combined.all)

# Phase 2 Part 6: HVG selection is restricted to integration_features (the
# 5,203-gene intersection measured in all 6 studies), computed on a
# separate subsetted copy so the vst variance ranking is never distorted by
# genes that are structurally zero (reference-absent, not biologically
# zero) in one or more studies - exactly the "encode study/reference
# identity into the matrix" failure mode Part 6 warns against. The main
# object keeps its full union assay (needed for per-study raw/normalized
# exports in flatten_data.R); only VariableFeatures/scale.data/PCA are
# computed from the integration-eligible subset.
integration_subset <- subset(pv.combined.all, features = integration_features)
integration_subset <- FindVariableFeatures(
    integration_subset,
    selection.method = "vst",
    nfeatures = 2000
)
n_hvg <- length(VariableFeatures(integration_subset))
cat(sprintf(
  "[integration.R] Highly variable genes selected from the %d-gene integration-eligible universe: %d (contributed to PCA)\n",
  length(integration_features), n_hvg
))

VariableFeatures(pv.combined.all) <- VariableFeatures(integration_subset)
pv.combined.all <- ScaleData(pv.combined.all, features = VariableFeatures(pv.combined.all))
pv.combined.all <- RunPCA(pv.combined.all, features = VariableFeatures(pv.combined.all))
rm(integration_subset)
pv.combined.all <- RunUMAP(
    pv.combined.all,
    reduction = "pca",
    dims = 1:30,
    n.components = 3,
    n.threads = 10
)

pv.combined.all@reductions$pca_unintegrated <- pv.combined.all@reductions$pca
pv.combined.all@reductions$umap_unintegrated <- pv.combined.all@reductions$umap

# merge all seurat all study and run integration. NOTE: `group_by` is
# deliberately NOT passed here (D050) - it has no effect on installed
# Seurat/HarmonyIntegration and would misleadingly imply per-study grouping
# is happening via this argument. The actual grouping mechanism is the
# verified one-layer-per-study invariant asserted above; see the layer
# diagnostics printed before this call.
pv.combined.all <- IntegrateLayers(
    object = pv.combined.all,
    method = HarmonyIntegration,
    orig.reduction = "pca",
    new.reduction = "pca_integrated",
    verbose = FALSE
)

pv.combined.all <- JoinLayers(pv.combined.all)

# Run umap, findneighbors and find clusters
pv.combined.all <- RunUMAP(
    pv.combined.all,
    reduction = "pca_integrated",
    dims = 1:30,
    n.components = 3,
    n.threads = 10
)

pv.combined.all@reductions$umap_integrated <- pv.combined.all@reductions$umap

pv.combined.all <- FindNeighbors(
    pv.combined.all,
    reduction = "pca_integrated",
    dims = 1:30
)


for (res in seq(0.01, 0.1, by = 0.01)) {
    pv.combined.all <- FindClusters(
        pv.combined.all,
        resolution = res,
        graph.name = "RNA_snn"
    )
    colname <- paste0("RNA_snn_res.", res)
    pv.combined.all@meta.data[[colname]] <- Idents(pv.combined.all)
}


# Run t-SNE
pv.combined.all <- RunTSNE(
    object = pv.combined.all,
    reduction = "pca_integrated",
    dims = 1:30,
    seed.use = 123,
    dim.embed = 3,
    reduction.name = "tsne_integrated",
    reduction.key = "tSNE_",
    do.fast = TRUE
)

# Remove unneeded data
pv.combined.all@reductions$pca <- NULL
pv.combined.all@reductions$umap <- NULL
pv.combined.all@meta.data$seurat_clusters <- NULL
pv.combined.all@reductions$tsne <- NULL

saveRDS(pv.combined.all, file = "pv_all_studies.rds")
record_build_manifest(
    artifact_path = "pv_all_studies.rds",
    script_path = "scripts/integration.R",
    cell_count = ncol(pv.combined.all),
    seed = 123,
    notes = sprintf(
      "Fixes D050,D051: removed non-functional group_by argument, added explicit layer->study grouping assertions/diagnostics before Harmony. Phase 2 Part 6: HVG/PCA restricted to the %d-gene integration-eligible intersection (of %d-gene atlas union) - see audit/phase2_rebuild/feature_universe.tsv and feature_availability_by_study.tsv",
      length(integration_features), length(atlas_feature_universe)
    )
)
