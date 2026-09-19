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
study_num <- "36525464"
samples <- c("SRR20710498", "SRR20710499", "SRR20710500", "SRR20710501")

all_pv <- list()
for (sample_id in samples) {
  path <- file.path(
    paste0(
      "/srv/baia/prj/pvsca/pre_process_data/counts/",
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
#rm(list = samples)

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
set.seed(123456)
all_pv <- lapply(all_pv, function(x) {
  e.out <- emptyDrops(x, lower = 1) #change this value if needed
  x <- x[, which(e.out$FDR <= 0.001)]
})

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

# Add metadata via explicit per-run vectors, each keyed 1:1 to `samples`
# (never recycled via paste0(NA), which silently produced the literal string
# "NA" for day_post_infection/host_id - D024).
run_id <- samples
assert_cardinality(run_id, length(all_pv), label = "run_id")

study_label <- rep(STUDY_LABELS[["hazzard2022"]], length(all_pv))
pub_year <- rep(2022, length(all_pv))
# D005/D023-style fix: registered vs publication-protocol sequencer kept as
# two separate fields, never silently resolved. ENA registers HiSeq2500 for
# runs 500/501 and NovaSeq6000 for 498/499; the publication text states
# NovaSeq throughout - both are preserved, not merged.
sequencer_registered <- c("Illumina NovaSeq 6000", "Illumina NovaSeq 6000", "Illumina HiSeq 2500", "Illumina HiSeq 2500")
sequencer_protocol_publication <- rep("Illumina NovaSeq 6000", length(all_pv))
sc_technology <- rep("10x_Chromium", length(all_pv))
host_species <- c(
  "Anopheles freeborni",
  "Anopheles stephensi",
  "Saimiri boliviensis",
  "Saimiri boliviensis"
)
sample_type <- c(
  "Vector host: salivary gland",
  "Vector host: salivary gland",
  "Mammalian host: blood",
  "Mammalian host: blood"
)
tissue_or_sample_type <- c(
  "Vector salivary gland",
  "Vector salivary gland",
  "Host blood",
  "Host blood"
)
strain <- c("Sal1/Chesson", "Sal1/Chesson", "Sal1", "Sal1")
biological_replicate_id <- c("2", "1", "2", "1")

for (i in 1:length(all_pv)) {
  n_cells <- ncol(all_pv[[i]])
  assert_cardinality(colnames(all_pv[[i]]), n_cells, label = paste0("colnames(", i, ")"))

  all_pv[[i]]$study_pmid <- rep(study_num, n_cells)
  all_pv[[i]]$run_id <- rep(run_id[[i]], n_cells)
  all_pv[[i]]$study_label <- rep(study_label[[i]], n_cells)
  all_pv[[i]]$pub_year <- rep(pub_year[[i]], n_cells)
  all_pv[[i]]$sc_technology <- rep(sc_technology[[i]], n_cells)
  all_pv[[i]]$sequencer_registered <- rep(sequencer_registered[[i]], n_cells)
  all_pv[[i]]$sequencer_protocol_publication <- rep(sequencer_protocol_publication[[i]], n_cells)
  all_pv[[i]]$host_species <- rep(host_species[[i]], n_cells)
  # D024: true NA - host_id/day_post_infection are genuinely undetermined at
  # this preprocessing stage (D025: pending the remaining Hazzard2022
  # publication/accession/QC/provenance reconciliation), never the literal
  # string "NA".
  all_pv[[i]]$host_id <- rep(NA_character_, n_cells)
  all_pv[[i]]$day_post_infection <- rep(NA_integer_, n_cells)
  all_pv[[i]]$sample_type <- rep(sample_type[[i]], n_cells)
  all_pv[[i]]$tissue_or_sample_type <- rep(tissue_or_sample_type[[i]], n_cells)
  all_pv[[i]]$strain <- rep(strain[[i]], n_cells)
  all_pv[[i]]$treatment <- rep("No_Treatment", n_cells)
  all_pv[[i]]$biological_replicate_id <- rep(biological_replicate_id[[i]], n_cells)
  all_pv[[i]]$biological_replicate_type <- rep("infection", n_cells)

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
#for (i in 1:length(all_pv)) {
#  all_pv[[i]] <- NormalizeData(all_pv[[i]], verbose = FALSE)
#  all_pv[[i]] <- FindVariableFeatures(all_pv[[i]], selection.method = "vst",
#                                      nfeatures = nrow(all_pv[[i]])*.3,
#                                      verbose = FALSE)
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

# pv.combined.all <- NormalizeData(pv.combined.all)
# pv.combined.all <- FindVariableFeatures(pv.combined.all, selection.method = "vst",
# 										nfeatures = nrow(pv.combined.all) * 0.3)
# pv.combined.all <- ScaleData(pv.combined.all)
# pv.combined.all <- RunPCA(pv.combined.all)

# D006-style geography split (never a combined field); D022's QC/inclusion
# authority is explicitly deferred (DEC04/DEC08) - do not alter the QC
# population while the publication/QC reconciliation remains incomplete.
pv.combined.all$experimental_site <- "Walter Reed Army Institute of Research (Silver Spring, MD)"
pv.combined.all$registry_location <- "USA: District of Columbia"
pv.combined.all$sequencing_site <- NA_character_

# Conservative Phase 1 membership policy: preserve the exact 3,294-cell
# population; do not alter QC/inclusion while D022/DEC04/DEC08 remain open.
assert_cell_count(ncol(pv.combined.all), 3294L, label = "hazzard2022.rds")
assert_unique_cell_ids(colnames(pv.combined.all), label = "hazzard2022.rds colnames")
assert_no_literal_na_string(pv.combined.all$host_id, label = "hazzard2022$host_id")
assert_no_literal_na_string(pv.combined.all$day_post_infection, label = "hazzard2022$day_post_infection")

# Save Seurat object
saveRDS(pv.combined.all, file = "hazzard2022.rds")
record_build_manifest(
  artifact_path = "hazzard2022.rds",
  script_path = "scripts/hazzard2022_pv_analysis_script.R",
  cell_count = ncol(pv.combined.all),
  notes = "Fixes D023,D024,D069 (true NA, registered-vs-publication sequencer split, technology spelling); D022 QC-authority explicitly deferred per DEC04/DEC08; D020/D021 broad-stage/HPI fixes live in singleR.R"
)

# Remove all object
rm(list = ls())
gc()
