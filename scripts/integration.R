library(DropletUtils) # scRNA-seq empty droplet detection
library(Seurat) # contains scRNA-seq analysis functions

# set current working directory
setwd("/home/sopheap/pvsca_b/pre_process_data")

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

# merge all layers
pv.combined.all <- Reduce(function(x, y) merge(x, y = y), pv.combined.all)

pv.combined.all <- NormalizeData(pv.combined.all)
pv.combined.all <- FindVariableFeatures(
    pv.combined.all,
    selection.method = "vst",
    nfeatures = 2000
    # nfeatures = nrow(pv.combined.all) * 0.3
)
pv.combined.all <- ScaleData(pv.combined.all)
pv.combined.all <- RunPCA(pv.combined.all)
pv.combined.all <- RunUMAP(
    pv.combined.all,
    reduction = "pca",
    dims = 1:30,
    n.components = 3,
    n.threads = 10
)

pv.combined.all@reductions$pca_unintegrated <- pv.combined.all@reductions$pca
pv.combined.all@reductions$umap_unintegrated <- pv.combined.all@reductions$umap

# merge all seurat all study and run integration
pv.combined.all <- IntegrateLayers(
    object = pv.combined.all,
    method = HarmonyIntegration,
    orig.reduction = "pca",
    new.reduction = "pca_integrated",
    group_by = "study_label",
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
