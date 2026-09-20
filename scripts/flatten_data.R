# This script is used to clean and manipulate data for pvsca project

# load library
library(Seurat)
library(tidyverse)
library(janitor)

.plavisca_project_root <- local({
  candidates <- unique(c(
    Sys.getenv("PLAVISCA_PROJECT_ROOT", unset = NA_character_),
    getwd(),
    "/home/sopheap/pvsca_b"
  ))
  candidates <- candidates[!is.na(candidates) & nzchar(candidates)]
  hit <- candidates[file.exists(file.path(candidates, "pre_process_data", "scripts", "pipeline_lib.R"))]
  if (length(hit) == 0) {
    stop(
      "Could not locate the PlaViSca project root (looked for pre_process_data/",
      "scripts/pipeline_lib.R under $PLAVISCA_PROJECT_ROOT, the current working ",
      "directory, and the historical hard-coded path). Set the ",
      "PLAVISCA_PROJECT_ROOT environment variable to the directory containing ",
      "both pre_process_data/ and PlaViSca/.",
      call. = FALSE
    )
  }
  normalizePath(hit[[1]])
})
setwd(.plavisca_project_root)
source("pre_process_data/scripts/pipeline_lib.R")

# Section 13 (do not overwrite production artifacts): every export in this
# script writes to a clearly separated candidate/staging location, never
# directly to the deployed PlaViSca/data/*.rds files. Promoting a validated
# candidate export to the deployed location is a deliberate, separate,
# reviewed Phase 2 action - never an automatic side effect of running this
# script.
candidate_dir <- "pre_process_data/data/candidate_export"
dir.create(candidate_dir, showWarnings = FALSE, recursive = TRUE)

# load gff data (fixes the historical "./ref/..." path, which assumed a
# personal layout where ref/ lived directly under the project root; it
# actually ships inside the PlaViSca app repo)
gff_data <- readRDS("PlaViSca/ref/PvivaxP01_gff_data.rds")

# load seurat object
so <- readRDS("pre_process_data/pv_all_studies.rds")

# D053: capture metadata column names AFTER clean_names(), not before - the
# original bug captured `metadata` from the raw Seurat column names, then
# evaluated `any_of(metadata)` against the already-clean_names()-transformed
# mr_data, silently dropping any renamed column (orig.ident, nCount_RNA,
# nFeature_RNA, ...) since any_of() never errors on a missing name.
metadata_raw <- colnames(so@meta.data)

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
    metadata_raw,
    umap,
    umap_int,
    pca,
    pca_int,
    tsne
  )
) %>%
  clean_names()

# D053 fix continued: metadata names re-derived from the post-clean_names()
# columns that actually exist in mr_data, by cleaning the same raw name
# vector the same way clean_names() would - this is what any_of(metadata)
# below is actually evaluated against.
metadata <- janitor::make_clean_names(metadata_raw)
missing_metadata_cols <- setdiff(metadata, colnames(mr_data))
if (length(missing_metadata_cols) > 0) {
  pipeline_fail(sprintf(
    "D053 regression: expected metadata column(s) missing from mr_data after clean_names(): %s",
    paste(missing_metadata_cols, collapse = ", ")
  ))
}

# clean mr_data - drop internal/intermediate annotation-pipeline columns
# that are not part of the exported schema (D046,D054: reconciled against
# the current metadata_schema.tsv vocabulary, not the stale deployed
# artifact's 52-column schema).
internal_only_columns <- c(
  "pred_gametocyte", "liver_form", "refine_state",
  "femalegams_module_score1", "malegams_module_score1",
  grep("snn", colnames(mr_data), value = TRUE)
)

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
  ) %>%
  select(-any_of(internal_only_columns))

# AT17 release closure: one canonical cell order for every cell-indexed app
# export. The expression exports have always used lexicographically sorted
# cell keys via red_df; apply that same deterministic order to mr_data itself
# before any downstream extraction or serialization.
canonical_cell_ids <- sort(rownames(mr_data))
assert_unique_cell_ids(canonical_cell_ids, label = "canonical export cell IDs")
mr_data <- mr_data[canonical_cell_ids, , drop = FALSE]
assert_rowname_order_equal(
  rownames(mr_data), canonical_cell_ids,
  label_a = "mr_data", label_b = "canonical_cell_ids"
)

# D062: never export a literal "NA" string - this is a fail-loud build gate,
# not a silent app-side workaround.
character_cols <- names(mr_data)[vapply(mr_data, is.character, logical(1))]
for (col in character_cols) {
  assert_no_literal_na_string(mr_data[[col]], label = paste0("mr_data$", col))
}

# reduction and metadata
red_df <- mr_data %>%
  select(
    any_of(metadata),
    any_of(contains("umap")),
    any_of(contains("pca")),
    any_of(contains("t_sne"))
  )

# mr_data is already in the canonical sorted order; retain and assert it.
assert_rowname_order_equal(
  rownames(red_df), canonical_cell_ids,
  label_a = "red_df", label_b = "canonical_cell_ids"
)

# D055/AT16/AT17: every expression matrix is explicitly key-aligned to
# red_df before any cbind - exact rowname-set AND order equality asserted,
# never merely assumed from independent per-object sorting.
align_to_red_df <- function(expr_matrix, label) {
  expr_matrix <- expr_matrix[order(rownames(expr_matrix)), , drop = FALSE]
  assert_rowname_order_equal(
    rownames(red_df), rownames(expr_matrix),
    label_a = "red_df", label_b = label
  )
  expr_matrix
}

# normalized data
norm_exp <- Matrix::t(GetAssayData(so, layer = "data"))
norm_exp <- align_to_red_df(norm_exp, "norm_exp")
normalize_df <- cbind(red_df, as.matrix(norm_exp))

# raw data
raw_exp <- Matrix::t(GetAssayData(so, layer = "counts"))
raw_exp <- align_to_red_df(raw_exp, "raw_exp")
raw_df <- cbind(red_df, as.matrix(raw_exp))

# scale data
scale_exp <- Matrix::t(GetAssayData(so, layer = "scale"))
scale_exp <- align_to_red_df(scale_exp, "scale_exp")
scale_df <- cbind(red_df, as.matrix(scale_exp))

# AT01/AT17: no duplicated cell IDs, identical key sets across every export.
assert_unique_cell_ids(rownames(normalize_df), label = "normalize_df")
assert_rowname_set_equal(rownames(normalize_df), rownames(raw_df), "normalize_df", "raw_df")
assert_rowname_set_equal(rownames(normalize_df), rownames(scale_df), "normalize_df", "scale_df")
assert_rowname_order_equal(rownames(mr_data), rownames(normalize_df), "mr_data", "normalize_df")
assert_rowname_order_equal(rownames(mr_data), rownames(raw_df), "mr_data", "raw_df")
assert_rowname_order_equal(rownames(mr_data), rownames(scale_df), "mr_data", "scale_df")

# save normalized, raw and scaled data (candidate/staging paths - see
# section 13 note above; never the deployed PlaViSca/data/*.rds files)
saveRDS(normalize_df, file.path(candidate_dir, "normalize_df.rds"))
saveRDS(raw_df, file.path(candidate_dir, "raw_df.rds"))
saveRDS(scale_df, file.path(candidate_dir, "scale_df.rds"))

# pca_std_df and hvf_data
# split data per study
so_list <- SplitObject(so, split.by = "study_label")

# AT18: assert all six studies are represented before building the
# per-study export loop, rather than discovering a missing study only after
# top_genes_exp/hiv_data/pca_df silently omit it.
missing_studies <- setdiff(unname(STUDY_LABELS), names(so_list))
if (length(missing_studies) > 0) {
  pipeline_fail(sprintf(
    "pv_all_studies.rds is missing expected stud(y/ies) for per-study export: %s",
    paste(missing_studies, collapse = ", ")
  ))
}

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

# D052: save the correctly-accumulated top_genes_df (all six studies),
# never the stale last-loop-iteration top_genes_exp - the original bug
# reused the per-study loop variable name for the final saved object, so
# `save_data$top_genes_exp` silently contained only the LAST study
# processed (Mancio-Silva2022's 14,940 rows) instead of all six studies'
# accumulated data.
top_genes_df <- bind_rows(
  lapply(names(per_study_data), function(study_name) {
    df <- per_study_data[[study_name]]$top_genes_exp
    df$study <- study_name
    df
  })
)

# AT18: exactly six studies represented in the top-gene export.
n_top_gene_studies <- dplyr::n_distinct(top_genes_df$study)
if (n_top_gene_studies != length(STUDY_LABELS)) {
  pipeline_fail(sprintf(
    "AT18 violation: top_genes_df contains %d distinct studies, expected %d",
    n_top_gene_studies, length(STUDY_LABELS)
  ))
}

save_data <- list(
  mr_data = mr_data,
  pca_df = pca_df,
  hiv_data = hiv_df,
  top_genes_exp = top_genes_df
)

# save list of data (candidate/staging path - see section 13 note above)
saveRDS(save_data, file.path(candidate_dir, "cleaned_dataset.rds"))
record_build_manifest(
  artifact_path = file.path(candidate_dir, "cleaned_dataset.rds"),
  script_path = "pre_process_data/scripts/flatten_data.R",
  cell_count = nrow(mr_data),
  notes = "Fixes D052,D053,D054,D055,D062; candidate/staging export only, deployed PlaViSca/data/*.rds untouched (section 13)"
)

# clear objects
# rm(list = ls())
