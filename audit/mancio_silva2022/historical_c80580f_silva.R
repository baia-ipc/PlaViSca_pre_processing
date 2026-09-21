library(Seurat)
library(dplyr)
# set current working directory
setwd("/srv/baia/prj/pvsca/pre_process_data")

df <- readxl::read_xls("data/TableS2_Cell_Metadata.xls", skip = 1)
df <- df |>
  mutate(
    refine_state = case_when(
      State == "Hypnozoite" ~ "Hypnozoite",
      TRUE ~ topcell
    ),
    refine_state = case_when(
      refine_state == "Male" ~ "Male gametocyte",
      refine_state == "Female" ~ "Female gametocyte",
      TRUE ~ refine_state
    )
  )


so <- readRDS("data/PvData_CHM_Final.RDS")

raw <- GetAssayData(so, assay = "RNA", layer = "counts")
rownames(raw) <- sub("^([^ -]*-[^ -]*).*", "\\1", rownames(raw))


silva <- CreateSeuratObject(counts = raw)

silva$study_pmid <- "39223117"
silva$run_id <- paste0("SRR18134", 227:284)
silva$study_label <- "Mancio_silva2022"
silva$num_srr <- 58
silva$pub_year <- 2022
silva$goegraphic_location <- "USA_RockvilleThailand_Ubon-Ratchathani"
silva$sc_technology <- NA
silva$sequncer <- c(
  rep("NextSeq 500", 16),
  rep("Illumina NovaSeq 6000", 42)
)
silva$host_species <- "Homo sapiens"
silva$host_id <- NA
silva$sample_type <- "Infected host hepatocytes"
silva$parasite_stage <- "Liver stage"
silva$strain <- "patient isolate"

silva$day_post_infection <- c(
  rep(1, 2),
  rep(5, 2),
  rep(8, 2),
  rep(11, 2),
  rep(14, 2),
  rep(8, 2),
  rep(11, 2),
  rep(14, 2),
  rep(1, 6),
  rep(4, 6),
  rep(5, 6),
  rep(8, 6),
  rep(11, 6),
  rep(8, 6),
  rep(11, 6)
)

silva$treatment <- c(
  rep("No_Treatment", 10),
  rep("No_Treatment", 6),
  rep("No_Treatment", 30),
  rep("No_Treatment", 12)
)

silva$biological_replicate <- NA

silva$refine_state <- NA
match_cell <- intersect(colnames(silva), df$Updated_names)
silva$refine_state[match_cell] <- df$refine_state[match(
  match_cell,
  df$Updated_names
)]

saveRDS(silva, file = "silva2022.rds")
