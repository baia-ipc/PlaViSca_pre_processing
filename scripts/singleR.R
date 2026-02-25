library(DropletUtils) # scRNA-seq empty droplet detection
library(readxl) # data import
library(rtracklayer) # contains functions to handle gff/gtf files
library(Seurat) # contains scRNA-seq analysis functions
library(SingleCellExperiment) # single-cell data storage
library(SingleR) # single-cell data identification tool
library(scCustomize) # single-cell data handling and plotting tools
library(tidyverse)

# set current working directory
setwd("/home/sopheap/pvsca_b/pre_process_data")

pv.combined.all <- readRDS("./pv_all_studies.rds")

# subset to see barcode of study silva2020
silva <- subset(pv.combined.all, subset = study_label == "Mancio Silva2022")


# First, upload (PMID:xxxx) data to R (data can be found in the OneDrive folder)
Zhu_HiDef_PVX_IDC_timecourse_1 <- read_excel(
  "./data/Zhu_SciReps_2016.xls",
  sheet = "SMRU1"
)

# Unfortunantely the genome used to analyze this data was the P. vivax Sal-I
# instead of the PvP01 (reference genome) that we used to align our data.
# We therefore need to convert the PVX_ gene ids to PVP01_ gene ids before
# we try to assign cellular identities to our data.

# Upload PVP01 orthologs
PVX_PVP01_orthologs <- read_excel("./data/GenesByOrthologs_Summary.xlsx")

# Process and tidy the timecourse data
Zhu_HiDef_PVX_IDC_timecourse_1 <- as.data.frame(Zhu_HiDef_PVX_IDC_timecourse_1)
Zhu_HiDef_PVX_IDC_timecourse_1 <- left_join(
  Zhu_HiDef_PVX_IDC_timecourse_1,
  PVX_PVP01_orthologs,
  by = c("Gid" = "InputOrtholog")
)
Zhu_HiDef_PVX_IDC_timecourse_1 <- drop_na(Zhu_HiDef_PVX_IDC_timecourse_1)
Zhu_HiDef_PVX_IDC_timecourse_1 <- Zhu_HiDef_PVX_IDC_timecourse_1[, c(2:9)]
rownames(
  Zhu_HiDef_PVX_IDC_timecourse_1
) <- Zhu_HiDef_PVX_IDC_timecourse_1$GeneID
Zhu_HiDef_PVX_IDC_timecourse_1 <- Zhu_HiDef_PVX_IDC_timecourse_1[, c(1:7)]
rownames(Zhu_HiDef_PVX_IDC_timecourse_1)
counts <- as.matrix(Zhu_HiDef_PVX_IDC_timecourse_1)
# Convert the dataframe to a summarized experiment
IDC_cycle <- SummarizedExperiment(
  assay = SimpleList(logcounts = counts),
  colData = colnames(Zhu_HiDef_PVX_IDC_timecourse_1)
)
IDC_cycle$X <- NULL

# rownames(IDC_cycle) <- gsub("-", "_", rownames(IDC_cycle))

# Convert Seurat object to single-cell experiment (sce) object for SingleR mapping
sce.all <- as.SingleCellExperiment(DietSeurat(
  pv.combined.all,
  layers = "data",
  assays = "RNA"
))
rownames(sce.all) <- gsub("_", "-", rownames(sce.all))


# Run SingleR
pred.IDC.all <- SingleR(
  test = sce.all,
  ref = IDC_cycle,
  assay.type.test = 1,
  labels = colnames(IDC_cycle)
)

# add pruned prediction time to seurat
pv.combined.all$hour_post_invasion <- pred.IDC.all$pruned.labels

# NA on hour post invasion for sample type rather than blood
pv.combined.all$hour_post_invasion[
  pv.combined.all$sample_type != "Mammalian host: blood"
] <- NA

# sopheap identify IDC time
pred.IDC.all$pruned.labels <- as.numeric(pred.IDC.all$pruned.labels)
pred.IDC.all$pruned.labels <- case_when(
  pred.IDC.all$pruned.labels >= 0 & pred.IDC.all$pruned.labels < 18 ~ "Ring",
  pred.IDC.all$pruned.labels >= 18 & pred.IDC.all$pruned.labels < 30 ~
    "Trophozoite",
  pred.IDC.all$pruned.labels >= 30 & pred.IDC.all$pruned.labels < 46 ~
    "Schizont",
  pred.IDC.all$pruned.labels >= 46 & pred.IDC.all$pruned.labels <= 48 ~
    "Merozoite",
  TRUE ~ NA
)

# add to seurat object
match_cell <- intersect(colnames(pv.combined.all), rownames(pred.IDC.all))

pv.combined.all$parasite_stages <- pred.IDC.all[match_cell, ]$pruned.labels

# NA for sample type rather than blood
pv.combined.all$parasite_stages[
  pv.combined.all$sample_type != "Mammalian host: blood"
] <- NA


plotDeltaDistribution(pred.IDC.all, ncol = 3)
# Add SingleR annotations to Seurat object
# pv.combined.all$parasite_stage <- pred.IDC.all$pruned.labels
# Replace NA vaues with "99"
# pv.combined.all$parasite_stage <- pv.combined.all$parasite_stage %>%
# replace_na("99")
# Convert annotations to factor and relevel
# pv.combined.all$IDC.pred <- factor(
#   pv.combined.all$IDC.pred,
#   levels = sort(as.numeric(unique(pv.combined.all$IDC.pred)))
# )

# overwrite life cycle stage with column liver form frm Ruberto2022_1
pv.combined.all$liver_form[
  pv.combined.all$liver_form == "Hypnozoites"
] <- "Hypnozoite"

pv.combined.all$liver_form[
  pv.combined.all$liver_form == "Schizonts"
] <- "Schizont"


liver_form <- pv.combined.all$liver_form %in% c("Hypnozoite", "Schizont")

pv.combined.all$parasite_stages[liver_form] <- pv.combined.all$liver_form[
  liver_form
]


# Update parasite_stages only where refine_state is not NA
pv.combined.all$parasite_stages[!is.na(pv.combined.all$refine_state)] <-
  pv.combined.all$refine_state[!is.na(pv.combined.all$refine_state)]


# overwrite for the sample type salivary gland to sporozoite. This is for Rubeto2022_2, hazzard2022 498 and 499

ruberto2022_2 <- pv.combined.all$study_label == "Ruberto2022_2"
pv.combined.all$parasite_stages[ruberto2022_2] <- "Sporozoite"

hazzard2022 <- pv.combined.all$run_id %in% c("SRR20710498", "SRR20710499")
pv.combined.all$parasite_stages[hazzard2022] <- "Sporozoite"

# visualize the number of cells in the dataset assigned to a parasite stage:
SCpubr::do_BarPlot(
  pv.combined.all,
  group.by = "parasite_stages",
  position = "stack",
  flip = T,
  order = T
)

# Output the exact # of cells predicted per IDC stage
# table(pv.combined.all$parasite_stages, useNA = "always")
# Plot the data in low-dimensional space
# before cell labelling
SCpubr::do_DimPlot(
  pv.combined.all,
  group.by = "parasite_stages",
  dims = c(1, 2),
  pt.size = 0.1,
  split.by = "parasite_stages"
)


# Let's now cluster the cell de novo
# Identify clusters of cells by a shared nearest neighbor (SNN) modularity
# optimization based clustering algorithm. First calculate k-nearest neighbors
# and construct the SNN graph. Then optimize the modularity function to determine
# clusters. For a full description of the algorithms, see Waltman and van Eck
# (2013) The European Physical Journal B.
# what is nearest-neighbour algorithm? https://www.ibm.com/topics/knn
# pv.combined.all <- FindNeighbors(
#   pv.combined.all,
#   dims = 1:30,
#   reduction = "pca_integrated"
# ) # change as necessary

pv.combined.all <- FindClusters(
  pv.combined.all,
  resolution = 0.05,
  graph_name = "new_clustering"
) # change as necessary

# ggarrange(
#   SCpubr::do_DimPlot(
#     pv.combined.all,
#     group.by = "seurat_clusters",
#     dims = c(1, 2),
#     pt.size = 0.1,
#     plot.axes = T,
#     reduction = "pca_integrated"
#   ),
#   SCpubr::do_DimPlot(
#     pv.combined.all,
#     group.by = "seurat_clusters",
#     dims = c(1, 3),
#     pt.size = 0.1,
#     plot.axes = T,
#     reduction = "pca_integrated"
#   ),
#   SCpubr::do_DimPlot(
#     pv.combined.all,
#     group.by = "seurat_clusters",
#     dims = c(2, 3),
#     pt.size = 0.1,
#     plot.axes = T,
#     reduction = "pca_integrated"
#   ),
#   align = "hv",
#   nrow = 1,
#   ncol = 3,
#   common.legend = T
# )

# pv.combined.all_stage_stats <- Cluster_Stats_All_Samples(
#   seurat_object = pv.combined.all,
#   group_by_var = "parasite_stages"
# )

# Do these clustering stats make sense? If you want to generate more or less
# clusters, rerun "FindClusters" with a different resolution.
# The closer to 0, the less clusters; closer to 1, more clusters.

# In Sa et al. (2020), they indicate that there are gametocytes in the samples.
# Can we detect these cells? Remember that the singleR cell identification
# strategy only contained asexual parasites. Let's try to identify gametocytes.

# We will assign a 'module score' to each cell based on their combined expression
# of known female and male gametocyte markers.

femaleGams <- c(
  "PVP01-1207200",
  "PVP01-0616100",
  "PVP01-1119300",
  "PVP01-1465500",
  "PVP01-1259400",
  "PVP01-1441000",
  "PVP01-0702600",
  "PVP01-1306800",
  "PVP01-0946800",
  "PVP01-0517400",
  "PVP01-1027600",
  "PVP01-1003000",
  "PVP01-1143200",
  "PVP01-1017500",
  "PVP01-1024300",
  "PVP01-0806000",
  "PVP01-1240300",
  "PVP01-0603600",
  "PVP01-0712800"
)

maleGams <- c(
  "PVP01-1262200",
  "PVP01-0530800",
  "PVP01-1412100",
  "PVP01-1025600",
  "PVP01-1266500",
  "PVP01-1229400"
)

pv.combined.all <- AddModuleScore(
  object = pv.combined.all,
  features = list(c(femaleGams)),
  name = "femaleGams_module_score"
)
pv.combined.all <- AddModuleScore(
  object = pv.combined.all,
  features = list(c(maleGams)),
  name = "maleGams_module_score"
)

VlnPlot_scCustom(
  seurat_object = pv.combined.all,
  features = "maleGams_module_score1"
)
VlnPlot_scCustom(
  seurat_object = pv.combined.all,
  features = "femaleGams_module_score1"
)

# table(pv.combined.all$RNA_snn_res.0.4, pv.combined.all$parasite_stages)

# Can you predict which clusters might be representative of gametocytes?
# Once you are happy with the clustering resolution, let's go ahead and rename
# the clusters identified.

# Below is an example of how to do so. Note that this need to be modified
# accordingly!

pv.combined.all <- RenameIdents(
  object = pv.combined.all,
  "0" = "Asexual 1",
  "1" = "Asexual 2",
  "2" = "Female gametocyte",
  "3" = "Asexual 3",
  "4" = "Male gametocyte",
  "5" = "Asexual 4",
  "6" = "Asexual 5",
  "7" = "Asexual 6",
  "8" = "Asexual 7",
  "9" = "Asexual 8",
  "10" = "Asexual 9",
  "11" = "Asexual 10",
  "12" = "Asexual 11",
  "13" = "Asexual 12",
  "14" = "Asexual 13",
  "15" = "Asexual 14",
  "16" = "Asexual 15",
  "17" = "Asexual 16",
  "18" = "Asexual 17",
  "19" = "Asexual 18",
  "20" = "Asexual 19",
  "21" = "Asexual 20"
)


pv.combined.all$pred_gametocyte <- as.character(Idents(pv.combined.all))

# add male and female gametocyte from prediction to life_cycle_stage
male_female <- pv.combined.all$pred_gametocyte %in%
  c("Male gametocyte", "Female gametocyte")

pv.combined.all$parasite_stages[
  male_female
] <- pv.combined.all$pred_gametocyte[male_female]


# DimPlot_scCustom(
#   pv.combined.all,
#   group.by = "RenamedClusters1",
#   reduction = "pca_integrated",
#   colors_use = ColorBlind_Pal()
# ) +
#   theme_grey() +
#   border()

# (SCpubr::do_DimPlot(
#   sample = pv.combined.all,
#   label = TRUE,
#   label.color = "black",
#   repel = T,
#   reduction = "pca_integrated",
#   dims = c(1, 2),
#   pt.size = 0.25,
#   colors.use = c(
#     "Asexual 1" = "orange",
#     "Asexual 2" = "#0072B2",
#     "Female Gametocyte 1" = "#009E73",
#     "Asexual 3" = "#CC79A7",
#     "Asexual 4" = "#F0E442",
#     "Female Gametocyte 2" = "firebrick2",
#     "Male Gametocyte" = "#56B4E9"
#   )
# ) +
#   NoLegend())

# read gametocyte metadata from Sa et al. and overwrite the life_cylce_stage
df <- read_excel("./data/pbio.3000711.s034.xlsx") |>
  select(sample = Sample, type = Type) |>
  mutate(
    prefix = sub("(_[^_]+)$", "", sample),
    cells = sub("^[^_]+_[^_]+_", "", sample)
  ) |>
  select(prefix, cells, type) |>
  mutate(
    prefix = case_when(
      prefix == "AMRU_Ao" ~ "278",
      prefix == "AMRU_CQ" ~ "277",
      prefix == "AMRU_Sa" ~ "276",
      prefix == "Ches_Ao" ~ "275",
      prefix == "Ches_CQ" ~ "274",
      prefix == "Ches_Sa" ~ "273",
      prefix == "NIH_Ao" ~ "269",
      prefix == "NIH_CQ" ~ "269",
      prefix == "NIH_Sa" ~ "272",
    ),
    new_cells = paste0(prefix, "_", cells),
    type = case_when(
      type == "Female" ~ "Female gametocyte",
      type == "Male" ~ "Male gametocyte"
    )
  ) |>
  filter(type == c("Male gametocyte", "Female gametocyte"))

# Find indices of Seurat cells that match new_cells
match_cell <- intersect(pv.combined.all$barcode, df$new_cells)


# Assign type using these indices
pv.combined.all$parasite_stages[match_cell] <- df$type[match(
  match_cell,
  df$new_cells
)]


# mutate parasite stages if sample type is liver call schizont (liver stage) and if blood call schizont (blood stage)
pv.combined.all$parasite_stages[
  pv.combined.all$sample_type == "Mammalian host: hepatocyte" &
    pv.combined.all$parasite_stages == "Schizont"
] <- "Schizont (liver stage)"

pv.combined.all$parasite_stages[
  pv.combined.all$sample_type == "Mammalian host: blood" &
    pv.combined.all$parasite_stages == "Schizont"
] <- "Schizont (blood stage)"


# ============================================================================
# DERIVE NEW 3-TIER TERMINOLOGY COLUMNS FROM LIFE_CYCLE_STAGE
# ============================================================================
# These columns provide a scientifically accurate 3-tier stage classification:
# 1. development_phase: 3 broad categories (Blood stages, Liver stages, Sporozoite)
# 2. blood_stage: 6 blood-specific forms + NA for non-blood
# 3. parasite_stages: 9 specific stage forms across all development phases
#
# Source of truth: life_cycle_stage (generated by SingleR predictions above)

pv.combined.all@meta.data <- pv.combined.all@meta.data %>%
  mutate(
    # 1. Development phase - 3 groups
    development_phase = case_when(
      parasite_stages == "Sporozoite" ~ "Sporozoite",
      parasite_stages == "Hypnozoite" ~ "Liver stages",
      parasite_stages == "Schizont (liver stage)" ~ "Liver stages",
      TRUE ~ "Blood stages"
    ),

    # 2. Blood stage - blood-specific forms only (NA for non-blood)
    blood_stage = case_when(
      parasite_stages == "Merozoite" ~ "Merozoite",
      parasite_stages == "Ring" ~ "Ring stage",
      parasite_stages == "Trophozoite" ~ "Trophozoite",
      parasite_stages == "Schizont (blood stage)" ~ "Schizont (blood)",
      parasite_stages == "Female gametocyte" ~ "Female gametocyte",
      parasite_stages == "Male gametocyte" ~ "Male gametocyte",
      TRUE ~ NA_character_
    ),

    # 3. Parasite stage - specific forms across all phases
    parasite_stages = case_when(
      parasite_stages == "Sporozoite" ~ "Sporozoite",
      parasite_stages == "Merozoite" ~ "Merozoite",
      parasite_stages == "Ring" ~ "Ring stage",
      parasite_stages == "Trophozoite" ~ "Trophozoite",
      parasite_stages == "Schizont (blood stage)" ~ "Schizont (blood)",
      parasite_stages == "Female gametocyte" ~ "Female gametocyte",
      parasite_stages == "Male gametocyte" ~ "Male gametocyte",
      parasite_stages == "Hypnozoite" ~ "Hypnozoite",
      parasite_stages == "Schizont (liver stage)" ~ "Schizont (liver)",
      TRUE ~ parasite_stages
    )
  )

# Remove Male and female gametocyte prediction with Anopheles
pv.combined.all <- subset(
  pv.combined.all,
  subset = !((parasite_stages %in% c("Male gametocyte", "Female gametocyte")) &
    grepl("Anopheles", host_species))
)


# Verify the new columns
cat("\n=== NEW TERMINOLOGY COLUMNS CREATED ===\n\n")
cat("Development phase (3 categories):\n")
print(table(pv.combined.all$development_phase, useNA = "ifany"))

cat("\n\nBlood stage (6 forms + NA for non-blood):\n")
print(table(pv.combined.all$blood_stage, useNA = "ifany"))

cat("\n\nParasite stage (9 specific forms):\n")
print(table(pv.combined.all$parasite_stages, useNA = "ifany"))

cat("\n=== VERIFICATION COMPLETE ===\n\n")


# clean up
# change features from - to _ to match with gff gene_id
features <- gsub("-", "_", rownames(pv.combined.all[["RNA"]]@features))
rownames(pv.combined.all[["RNA"]]@features) <- features


pv.combined.all$orig.ident <- NULL
pv.combined.all$femaleGams_module_score1 <- NULL
pv.combined.all$maleGams_module_score1 <- NULL
pv.combined.all@project.name <- "PlaViSca"

saveRDS(pv.combined.all, file = "./pv_all_studies.rds")
