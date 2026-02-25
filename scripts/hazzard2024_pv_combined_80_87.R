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
setwd("/home/sopheap/pvsca_b/pre_process_data")

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
samples <- c(paste0("SRR27021", 980:987))

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

# Add metadata. Note: the same will need to be done for the other samples that we analyze.
study_pmid <- rep(study_num, length(all_pv))
run_id <- rep(paste0("SRR27021", 980:987), length(all_pv))
study_label <- rep("Hazzard2024", length(all_pv))
num_srr <- rep(24, length(all_pv))
pub_year <- rep(2024, length(all_pv))
geographic_location <- rep("Papua New Guinea/El Salvador", length(all_pv))
sc_technology <- rep("10x_Chomium_V3", length(all_pv))
sequencer <- rep("Illumina NovaSeq 6000", length(all_pv))
host_species <- rep("Saimiri boliviensis", length(all_pv))
host_id <- c(
    "5163_2",
    "5163_1",
    "5708_2",
    "5550_1",
    "5537_4",
    "5537_3",
    "5537_3",
    "5537_1"
)

sample_type <- rep("Mammalian host: blood", length(all_pv))

strain <- rep("NIH1993-F3&Chesson", length(all_pv))
day_post_infection <- c(16, 14, 21, 17, 26, 24, 21, 16)
treatment <- rep("No_Treatment", length(all_pv))
biological_replicate <- rep(NA, length(all_pv))


for (i in 1:length(all_pv)) {
    all_pv[[i]]$study_pmid <- paste0(study_pmid[[i]])
    all_pv[[i]]$run_id <- paste0(run_id[[i]])
    all_pv[[i]]$study_label <- paste0(study_label[[i]])
    all_pv[[i]]$num_srr <- paste0(num_srr[[i]])
    all_pv[[i]]$pub_year <- paste0(pub_year[[i]])
    all_pv[[i]]$geographic_location <- paste0(geographic_location[[i]])
    all_pv[[i]]$sc_technology <- paste0(sc_technology[[i]])
    all_pv[[i]]$sequencer <- paste0(sequencer[[i]])
    all_pv[[i]]$host_species <- paste0(host_species[[i]])
    all_pv[[i]]$host_id <- paste0(host_id[[i]])
    all_pv[[i]]$sample_type <- paste0(sample_type[[i]])
    all_pv[[i]]$strain <- paste0(strain[[i]])
    all_pv[[i]]$day_post_infection <- paste0(day_post_infection[[i]])
    all_pv[[i]]$treatment <- paste0(treatment[[i]])
    all_pv[[i]]$biological_replicate <- paste0(biological_replicate[[i]])

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
saveRDS(pv.combined.all, file = "hazzard2024_pv_combined_80_87.rds")
# Remove all object
rm(list = ls())
gc()
