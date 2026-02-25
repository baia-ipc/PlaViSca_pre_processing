library(readr)
library(dplyr)
library(stringr)

setwd("/home/sopheap/pvsca_b/pre_process_data")

# Publication dates
Publication_date <- c(
  "4-May-20",
  "4-Aug-22",
  "25-Aug-22",
  "16-Nov-22",
  "2-Sep-24",
  "19-Apr-22"
)

# Study labels
Study_label <- c(
  "Sar2020",
  "Ruberto2022_2",
  "Ruberto2022_1",
  "Hazzard2022",
  "Hazzard2024",
  "Mancio-Silva2022"
)

# Authors
Authors <- c(
  "Sà JM et al.",
  "Ruberto AA et al.",
  "Ruberto AA et al.",
  "Hazzard B et al.",
  "Hazzard B et al.",
  "Mancio-Silva L et al."
)

# Journals
Journal <- c(
  "PLoS Biol",
  "PLoS Negl Trop Dis",
  "Front. Cell. Infect. Microbiol.",
  "PLoS Negl Trop Dis",
  "Nat Commun",
  "Cell Host & Microbe"
)

# Titles
Title <- c(
  "Single-cell transcription analysis of <em>Plasmodium vivax</em> blood-stage parasites identifies stage- and species-specific profiles of expression",
  "Single-cell RNA sequencing of <em>Plasmodium vivax</em> sporozoites reveals stage- and species-specific transcriptomic signatures",
  "Single-cell RNA profiling of <em>Plasmodium vivax</em>-infected hepatocytes reveals parasite- and host- specific transcriptomic signatures and therapeutic targets",
  "Long read single cell RNA sequencing reveals the isoform diversity of <em>Plasmodium vivax</em> transcripts",
  "Single-cell analyses of polyclonal <em>Plasmodium vivax</em> infections and their consequences on parasite transmission",
  "A single-cell liver atlas of <em>Plasmodium vivax</em> infection"
)

# DOIs with links
DOI <- c(
  "<a href='https://doi.org/10.1371/journal.pbio.3000711' target='_blank'>https://doi.org/10.1371/journal.pbio.3000711</a>",
  "<a href='https://doi.org/10.1371/journal.pntd.0010633' target='_blank'>https://doi.org/10.1371/journal.pntd.0010633</a>",
  "<a href='https://doi.org/10.3389/fcimb.2022.986314' target='_blank'>https://doi.org/10.3389/fcimb.2022.986314</a>",
  "<a href='https://doi.org/10.1371/journal.pntd.0010991' target='_blank'>https://doi.org/10.1371/journal.pntd.0010991</a>",
  "<a href='https://doi.org/10.1038/s41467-024-51949-8' target='_blank'>https://doi.org/10.1038/s41467-024-51949-8</a>",
  "<a href='https://doi.org/10.1016/j.chom.2022.03.034' target='_blank'>https://doi.org/10.1016/j.chom.2022.03.034</a>"
)

# Number of cells
Number_of_cells <- c(
  9766,
  9947,
  1438,
  3294,
  80024,
  1494
)

df <- cbind(
  Publication_date,
  Study_label,
  Authors,
  Journal,
  Title,
  DOI,
  Number_of_cells
)
df <- as.data.frame(df)

# Ensure the folder exists, then write to CSV

write_csv(df, "data/data_source.csv")

# Print the data frame
print(df)
