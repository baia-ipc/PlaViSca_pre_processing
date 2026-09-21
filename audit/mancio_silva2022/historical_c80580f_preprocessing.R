# This script is used to clean and manipulate data for pvsca project

# load library
library(Seurat)
library(tidyverse)
library(janitor)

setwd("/srv/baia/prj/pvsca")
# load gff data
gff_data <- readRDS("ref/PvivaxP01_gff_data.rds")

# load seurat object
so <- readRDS("./pre_process_data/pv_all_studies.rds")

metadata <- colnames(so@meta.data)

# subset reduction columns from seurat object
umap <- colnames(so@reductions$umap_unintegrated@cell.embeddings)
umap_int <- colnames(so@reductions$umap_integrated@cell.embeddings)
pca <- colnames(so@reductions$pca_unintegrated@cell.embeddings)[1:3]
pca_int <- colnames(so@reductions$pca_integrated@cell.embeddings)[1:3]
tsne <- colnames(so@reductions$tsne_integrated@cell.embeddings)

# Extract metadata and reduction data
mr_data <- FetchData(
  so,
  vars = c(
    metadata,
    umap,
    umap_int,
    pca,
    pca_int,
    tsne
  )
) %>%
  clean_names()

# modify life_cycle_stage to include parasite_stage for schizont stage
mr_data <- mr_data %>%
  mutate(
    life_cycle_stage = case_when(
      life_cycle_stage == "Schizont" ~
        paste0(life_cycle_stage, " (", parasite_stage, ")"),
      TRUE ~ life_cycle_stage
    )
  )


# clean mr_data
mr_data <- mr_data %>%
  dplyr::rename(
    umap_i_1 = cnqwt_1,
    umap_i_2 = cnqwt_2,
    umap_i_3 = cnqwt_3,
    umap_u_1 = umap_1,
    umap_u_2 = umap_2,
    umap_u_3 = umap_3,
    pca_i_1 = harmony_1,
    pca_i_2 = harmony_2,
    pca_i_3 = harmony_3,
    pca_u_1 = pc_1,
    pca_u_2 = pc_2,
    pca_u_3 = pc_3
  )


# reduction and metadata
red_df <- mr_data %>%
  select(
    any_of(metadata),
    any_of(contains("umap")),
    any_of(contains("pca")),
    any_of(contains("tsne"))
  )
# sort by row names
red_df <- red_df[order(rownames(red_df)), ]

# normalized data
norm_exp <- GetAssayData(so, layer = "data")
norm_exp <- Matrix::t(norm_exp)
norm_exp <- norm_exp[order(rownames(norm_exp)), ]

# combined normalized data
normalize_df <- cbind(red_df, as.matrix(norm_exp))


# raw data
raw_exp <- GetAssayData(so, layer = "counts")
raw_exp <- Matrix::t(raw_exp)
raw_exp <- raw_exp[order(rownames(raw_exp)), ]

# combined raw data
raw_df <- cbind(red_df, as.matrix(raw_exp))

# scale data
scale_exp <- GetAssayData(so, layer = "scale")
scale_exp <- Matrix::t(scale_exp)
scale_exp <- scale_exp[order(rownames(scale_exp)), ]

# combined raw data
scale_df <- cbind(red_df, as.matrix(scale_exp))

# save normalized, raw and scaled data
saveRDS(normalize_df, "data/normalize_df.rds")
saveRDS(raw_df, "data/raw_df.rds")
saveRDS(scale_df, "data/scale_df.rds")

# pca_std_df and hvf_data
# split data per study
so_list <- SplitObject(so, split.by = "study_label")

per_study_data <- list()
for (study_name in names(so_list)) {
  so_sub <- so_list[[study_name]]

  # Normalize
  so_sub <- NormalizeData(so_sub, assay = "RNA", verbose = FALSE)

  # PCA
  so_sub <- RunPCA(so_sub, assay = "RNA", verbose = FALSE)
  pca_std_df <- data.frame(
    pc = 1:length(so_sub@reductions$pca@stdev),
    std = so_sub@reductions$pca@stdev
  )

  # Highly variable genes (2000 by default)
  so_sub <- FindVariableFeatures(
    so_sub,
    assay = "RNA",
    selection.method = "vst",
    nfeatures = 2000
  )

  hiv_data <- HVFInfo(so_sub, assay = "RNA", method = "vst", status = TRUE) %>%
    rownames_to_column("id") %>%
    left_join(gff_data, by = "id") %>%
    select(id, name, mean, variable, rank, variance.standardized) %>%
    mutate(name = coalesce(name, id))

  # Top expressed genes
  top_genes <- rowMeans(GetAssayData(so_sub, assay = "RNA", layer = "data"))
  top_genes <- names(sort(top_genes, decreasing = TRUE)[1:10])

  top_genes_exp <- FetchData(so_sub, vars = c(top_genes)) %>%
    rownames_to_column("id") %>%
    pivot_longer(
      any_of(top_genes),
      names_to = "gene",
      values_to = "expression"
    ) %>%
    left_join(gff_data, by = c("gene" = "id")) %>%
    select(id, gene, name, expression) %>%
    mutate(name = coalesce(name, gene))

  # Store per-study results
  per_study_data[[study_name]] <- list(
    pca_std_df = pca_std_df,
    hiv_data = hiv_data,
    top_genes_exp = top_genes_exp
  )
}


pca_df <- bind_rows(
  lapply(names(per_study_data), function(study_name) {
    df <- per_study_data[[study_name]]$pca_std_df
    df$study <- study_name
    df
  })
)

hiv_df <- bind_rows(
  lapply(names(per_study_data), function(study_name) {
    df <- per_study_data[[study_name]]$hiv_data
    df$study <- study_name
    df
  })
)

top_genes_df <- bind_rows(
  lapply(names(per_study_data), function(study_name) {
    df <- per_study_data[[study_name]]$top_genes_exp
    df$study <- study_name
    df
  })
)

save_data <- list(
  mr_data = mr_data,
  pca_df = pca_df,
  hiv_data = hiv_df,
  top_genes_exp = top_genes_exp
)

# # save list of data
saveRDS(save_data, "data/cleaned_dataset.rds")

# clear objects
# rm(list = ls())
