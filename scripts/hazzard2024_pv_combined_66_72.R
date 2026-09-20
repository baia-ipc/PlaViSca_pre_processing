##### --------------------------------------------------------------------------
##### -Ruberto, CTEGD, Institute of Bioinformatics, UGA
##### --------------------------------------------------------------------------

##### --------------------------------------------------------------------------
##### 1. Load libraries
##### --------------------------------------------------------------------------
library(DropletUtils)
library(tidyverse)
library(rtracklayer)
library(Seurat)
library(scater)
library(SingleCellExperiment)

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
##### 2. Load annotation files
##### --------------------------------------------------------------------------
# Before beginning the data processing steps, let's upload the PvP01 gene
# annotation file. This will come in handy when we perform differential gene
# expression analyses and for the generation of the data tables.

pv.gff <- import.gff3("./scripts/ref/PlasmoDB-68_PvivaxP01.gff")
pv.gff <- as.data.frame(pv.gff)
gene.info <- pv.gff %>%
    mutate(seurat = gsub("_", "-", ID)) %>%
    filter(type == "protein_coding_gene")

gene.info <- as.data.frame(gene.info)
gene.info$GeneDescription <- paste(
    gene.info$ID,
    gene.info$description,
    sep = "::"
)

rRNA_pv <- pv.gff %>%
    as_tibble() %>%
    filter(type == "rRNA") %>%
    select(ID, description) %>%
    mutate(ID = gsub("\\.1", "", ID))

##### --------------------------------------------------------------------------
##### 3. Load scRNAseq data in to R
##### --------------------------------------------------------------------------
# Great, we will now  upload the aligned data to R. We will use the data from
# Sa et al. as an example, PMID:32365102. This dataset contains 10 single-cell
# RNAs-seq datasets from P. vivax parasites.
# Choose one or more of these single-cell RNA-seq datasets for processing.
# For the data that you choose to analyze, simply unhash the sample
# corresponding to the the one that you choose. Be sure to also change the
# path of the read_count_input() function so that it matches the location on your
# system

# For the sample(s) that you choose, please provide information linked to how
# the sample(s) were processed. Simply go to https://www.ncbi.nlm.nih.gov/sra,
# type in the SRR accession number and you will find this information.
study_num <- "39223117"
samples <- c(paste0("SRR27021", 966:972))

all_pv <- list()
for (sample_id in samples) {
    path <- file.path(
        paste0(
            "./counts/",
            study_num,
            "/",
            sample_id,
            "_solo_out/Solo.out/GeneFull/raw"
        )
    )

    sce <- read10xCounts(path, col.names = TRUE)
    assign(sample_id, sce)
    all_pv[[sample_id]] <- sce
}

rm(list = samples)

##### --------------------------------------------------------------------------
##### 4. Processing step 1: select droplets containing cells
##### --------------------------------------------------------------------------
# We will now process the *P. vivax* mixed blood stage parasite data.
# Note: When running Kallisto Bustools, we generated both an unfiltered
# countmatrix. This means that we need to run a droplet detection step to
# select droplets containing cells versus droplets did not capture any cells or
# droplets that captured ambient RNA

# Since we are interested in assessing changes in protein-coding transcripts,
# we can remove rRNA from the dataset.

for (i in 1:length(all_pv)) {
    all_pv[[i]] <- all_pv[[i]][!rownames(all_pv[[i]]) %in% rRNA_pv$ID, ]
}
# Let's now generate a 'knee plot' to visualize the distribution of RNA per cell.
# This provides a nice visual for determining droplet-containing versus
# empty versus ambient RNA containing cells.
# Visualizes the inflection point to filter empty droplets

# Create the PDF file of knee plot
pdf(paste0(study_num, "_knee_plots.pdf"), width = 12, height = 6)
par(mfrow = c(1, 2))

bcrank <- list()
uniq <- list()

for (i in seq_along(all_pv)) {
    bcrank[[i]] <- barcodeRanks(all_pv[[i]])
    uniq[[i]] <- !duplicated(bcrank[[i]]$rank)

    plot(
        bcrank[[i]]$rank[uniq[[i]]],
        bcrank[[i]]$total[uniq[[i]]],
        log = "xy",
        xlab = "Rank",
        ylab = "Total UMI Count",
        cex.lab = 1.2,
        main = names(all_pv)[i]
    )

    abline(h = metadata(bcrank[[i]])$inflection, col = "darkgreen", lty = 2)
    abline(h = metadata(bcrank[[i]])$knee, col = "dodgerblue", lty = 2)
    legend(
        "bottomleft",
        legend = c("Inflection", "Knee"),
        col = c("darkgreen", "dodgerblue"),
        lty = 2,
        cex = 1.2
    )

    # Start new page after every 2 plots
    if (i %% 2 == 0 && i < length(all_pv)) {
        par(mfrow = c(1, 2))
    }
}

dev.off()

rm(bcrank, uniq)

# Another option is to use the emptyDrops.
# emptyDrops performs Monte Carlo simulations to compute p-values, so we need to
# set the seed to obtain reproducible results.
# see PMID: 30902100 for  rationale and statistical framework underlying this
# method
# set.seed(123456)
# all_pv <- lapply(all_pv, function(x) {
#     e.out <- emptyDrops(x, lower = 150) #change this value if needed
#     x <- x[, which(e.out$FDR <= 0.001)]
# })

##### --------------------------------------------------------------------------
##### 5. Processing step 2: conversion of dgCMatrix to a Seurat object
##### --------------------------------------------------------------------------
# There are various tools
# available to process and analyze single-cell RNA-seq data. The Seurat suite
# https://satijalab.org/seurat/ is one of the most popular thanks to its regular
# updates, its ease of use, and streamlined data handling and processing
# commands. In what follows, we will follow the standard steps for handling,
# processing, visualizing, and comparing data in the Seurat suite
# (https://satijalab.org/seurat/articles/pbmc3k_tutorial).

# Let's now transform the data in to a Seurat object.

for (i in 1:length(all_pv)) {
    sce <- all_pv[[i]]
    mat <- counts(sce)
    all_pv[[i]] <- CreateSeuratObject(
        counts = mat,
        min.cells = 1,
        min.features = 1
    )
}

# Note the use of the min.cells and min.features aruguments in the
# CreateSeuratObject function. What do these arguments do?
# In the console run '?CreateSeuratObject' to access the description of the
# function and its arugments.
# After transforming the data, check out the dimensions of the new Seurat object.
# How has the noumber of cells and the number of features changed?

# Add metadata via a keyed authoritative run table - never positional
# vectors/rep() recycling (fixes D026-D030,D032,D033,D069; see
# scripts/pipeline_lib.R and scripts/ref/hazzard2024_run_table.tsv, derived
# from audit/hazzard2024/authoritative_run_table.tsv).
source("scripts/pipeline_lib.R")

run_table <- read.delim(
  "scripts/ref/hazzard2024_run_table.tsv",
  stringsAsFactors = FALSE
)

membership <- split_run_membership(all_pv, run_table)

retained_fields <- c(
  "library_id", "animal_id", "biological_replicate_id",
  "biological_replicate_type", "day_post_infection", "infection_design",
  "parasite_lineage_isolate", "per_cell_parasite_genotype",
  "parasite_origin_location", "host_provenance", "experimental_site",
  "sequencing_site", "registry_location", "host_species", "technology",
  "sequencer_registered", "sequencer_protocol_publication", "host_sex",
  "source_run_count_deposited", "source_run_count_imported",
  "source_run_count_retained"
)

if (length(membership$retained) > 0) {
  membership$retained <- apply_keyed_run_metadata(
    membership$retained,
    run_table,
    retained_fields
  )
  for (rid in names(membership$retained)) {
    n_cells <- ncol(membership$retained[[rid]])
    membership$retained[[rid]]$sample_type <- rep("Host blood", n_cells)
    membership$retained[[rid]]$tissue_or_sample_type <- rep("Host blood", n_cells) # Phase 2 fix: this field, not sample_type, is what singleR.R reads for IDC eligibility
    membership$retained[[rid]]$retained_by_plavisca_qc <- rep(TRUE, n_cells)
  }
}

# Runs imported for QC but not among the 24 authoritative retained blood
# libraries (see scripts/hazzard2024_merge_all.R's explicit subset() against
# data/Proccessed_Data.txt) never survive to the saved study object; assign
# explicit NA placeholders rather than any authoritative value, since none of
# the corrected fields above are established for a population never retained.
for (rid in names(membership$pending_filter)) {
  n_cells <- ncol(membership$pending_filter[[rid]])
  for (f in retained_fields) {
    membership$pending_filter[[rid]][[f]] <- rep(NA, n_cells)
  }
  membership$pending_filter[[rid]]$sample_type <- rep(NA_character_, n_cells)
  membership$pending_filter[[rid]]$tissue_or_sample_type <- rep(NA_character_, n_cells)
  membership$pending_filter[[rid]]$retained_by_plavisca_qc <- rep(FALSE, n_cells)
}

all_pv <- c(membership$retained, membership$pending_filter)[names(all_pv)]

for (i in seq_along(all_pv)) {
  n_cells <- ncol(all_pv[[i]])
  assert_cardinality(colnames(all_pv[[i]]), n_cells, label = paste0("colnames(", names(all_pv)[i], ")"))

  all_pv[[i]]$study_pmid <- rep(study_num, n_cells)
  all_pv[[i]]$run_id <- rep(names(all_pv)[i], n_cells)
  all_pv[[i]]$study_label <- rep(STUDY_LABELS[["hazzard2024"]], n_cells)
  all_pv[[i]]$pub_year <- rep(2024L, n_cells)
  all_pv[[i]]$sc_technology <- all_pv[[i]]$technology
  all_pv[[i]]$sequencer <- all_pv[[i]]$sequencer_registered
  # No drug-treatment arm exists in this study (infection design - mono/
  # consecutive/simultaneous/sporozoite - is a separate covariate carried in
  # infection_design, not a treatment arm); "No_Treatment" is the genuinely
  # correct value here, never a literal "NA" string.
  all_pv[[i]]$treatment <- rep("No_Treatment", n_cells)
  all_pv[[i]]$source_treatment <- rep(NA_character_, n_cells)
  all_pv[[i]]$host_id <- all_pv[[i]]$library_id

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
#     all_pv[[i]] <- NormalizeData(all_pv[[i]], verbose = FALSE)
#     all_pv[[i]] <- FindVariableFeatures(all_pv[[i]], selection.method = "vst",
#                                         nfeatures = nrow(all_pv[[i]])*.3,
#                                         verbose = FALSE)
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

# pv.combined.all <- NormalizeData(pv.combined.all)
# pv.combined.all <- FindVariableFeatures(pv.combined.all)
# pv.combined.all <- ScaleData(pv.combined.all)
# pv.combined.all <- RunPCA(pv.combined.all)

# Save Seurat object
saveRDS(pv.combined.all, file = "hazzard2024_pv_combined_66_72.rds")
# Remove all object
rm(list = ls())
gc()
