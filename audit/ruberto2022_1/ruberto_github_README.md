# Workflow for Single-cell RNA profiling of *Plasmodium vivax*-infected hepatocytes reveals parasite- and host- specific transcriptomic signatures and therapeutic targets
#### Anthony A. Ruberto<sup>1</sup>, Steven P. Maher<sup>1</sup>, Amélie Vantaux<sup>2</sup>, Chester J Joyner<sup>1,3</sup>, Caitlin Bourke<sup>4</sup>, Balu Balan<sup>4,5</sup>, Aaron Jex<sup>4,5</sup>, Ivo Mueller<sup>4</sup>, Benoit Witkowski<sup>2</sup>, Dennis E Kyle<sup>1</sup>
#### 1 Center for Tropical and Emerging Tropical Diseases, University of Georgia, Athens, GA, USA
#### 2 Malaria Molecular Epidemiology Unit, Institut Pasteur du Cambodge, Phnom Penh, Cambodia
#### 3 Department of Infectious Diseases, College of Veterinary Medicine, University of Georgia, Athens, GA, USA
#### 4 Division of Population Health and Immunity, Walter and Eliza Hall Institute of Medical Research, Parkville, Victoria, Australia
#### 5 Faculty of Veterinary and Agricultural Sciences, University of Melbourne, Parkville, Victoria, Australia

![](https://www.ncbi.nlm.nih.gov/corehtml/pmc/pmcgifs/logo-fcimb.gif)

[![DOI](https://zenodo.org/badge/DOI/10.3389/fcimb.2022.986314.svg)](https://doi.org/10.3389/fcimb.2022.986314)


## Content

R Markdown files containing codes used to analyze *Plasmodium vivax*-infected hepatocytes.

To perform the analysis in its entirety, download and unzip "Pv_Hep_scRNAseq_Shell_scripts_outputs" from Zenodo.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6463338.svg)](https://doi.org/10.5281/zenodo.6463338)

Once complete you will find 6 folders that correspond to the each of the *Plasmodium vivax*-infected hepatocytes 10x Genomics scRNA-seq runs.

Modify the path of these folders in the Pv_Hep_scRNAseq_additional_file_1.Rmd so that they correspond to the directory that they are located on your system.

When saving the RDS files for each analysis step, be sure you have the correct output path specified.

## Markdown file descriptions

**Pv_Hep_scRNAseq_additional_file_1.Rmd**: In this document, we distinguish between droplets containing cells from those containing ambient RNA / dead/ dying cells.

**Pv_Hep_scRNAseq_additional_file_2.Rmd**: In this document, we first process the data to filter reads from human hepatocytes, leaving us with only *P. vivax* transcripts. Next, we perform cell and gene filtering for each of the samples. Last, we perform data reduction and clustering. 

**Pv_Hep_scRNAseq_additional_file_3.Rmd**: In this document, we merge the *P. vivax* data. We will first merge data derived from replicate 1 (day 5 and day 9) and replicate 2 (day 5 and day 9) separately. Second, we will merge all data (replicates 1 and 2). The rationale for assessing the data from replicate 1 and 2 separately gives us an opportunity to pick up any biases that may arise that are biological or technical in nature. Ways these biases would manifest include differences in clustering outputs or differential gene expression testing. 

**Pv_Hep_scRNAseq_additional_file_4.Rmd**: In this document, we perform the differential gene expression analysis used in the
manuscript for the parasite side.

**Pv_Hep_scRNAseq_additional_file_5.Rmd**: In this document, we first process the data to filter reads from *P.vivax*, leaving us with only human transcripts. Next, we perform cell and gene filtering for each of the samples. Last, we perform data reduction and clustering.

**Pv_Hep_scRNAseq_additional_file_6.Rmd**: In this document, we first merge the processed human hepatocyte scRNAseq data from each replicate. Next, we perform low dimensional reduction and clustering. Finally, the cells are parsed and differential gene expression analyses are performed.

## Outputs

In the event you do not want to perform the analysis, the .rds files (Pv_Hep_scRNAseq_R_scripts_outputs.zip) from each step can be downloaded from Zenodo.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6463338.svg)](https://doi.org/10.5281/zenodo.6463338)


## Support

If there are any issues or questions, do not hesitate to contact Anthony Ruberto [aruberto@uga.edu].






