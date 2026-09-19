#!/usr/bin/env Rscript

# Deterministic, read-only forensic inspection for Hazzard et al. 2024.
# Run from pre_process_data/scripts with:
# env -u R_LIBS_USER pixi run Rscript --vanilla ../audit/hazzard2024/inspect_hazzard2024.R

suppressPackageStartupMessages({
    library(digest)
    library(dplyr)
    library(readxl)
    library(Seurat)
    library(xml2)
})

audit_dir <- normalizePath("../audit/hazzard2024", mustWork = TRUE)
data_dir <- normalizePath("../data", mustWork = TRUE)
pre_dir <- normalizePath("..", mustWork = TRUE)
app_dir <- normalizePath("../../PlaViSca/data", mustWork = TRUE)

write_tsv <- function(x, filename) {
    write.table(
        x,
        file.path(audit_dir, filename),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE,
        na = "NA"
    )
}

stopifnot(
    tools::md5sum(file.path(data_dir, "Proccessed_Data.txt")) ==
        "6ede1f7da4e5d57841c9f90883b2b2c1",
    tools::md5sum(file.path(audit_dir, "supplemental_data_1.xlsx")) ==
        "702418299ebb4d4e3af441fec9966abc"
)

# The merger's explicit, auditable library-to-run crosswalk.
group_to_run <- c(
    `4881_1` = "SRR27021994", `5142_1` = "SRR27021992",
    `5142_2` = "SRR27021991", `5163_1` = "SRR27021981",
    `5163_2` = "SRR27021980", `5164_1` = "SRR27021976",
    `5164_2` = "SRR27021975", `5309_1` = "SRR27021956",
    `5309_2` = "SRR27021955", `5350_1` = "SRR27021972",
    `5350_2` = "SRR27021970", `5350_3` = "SRR27021969",
    `5370_1` = "SRR27021968", `5370_2` = "SRR27021967",
    `5370_3` = "SRR27021966", `5537_1` = "SRR27021987",
    `5537_2` = "SRR27021986", `5537_3` = "SRR27021985",
    `5537_4` = "SRR27021984", `5550_1` = "SRR27021983",
    `5708_1` = "SRR27021993", `5708_2` = "SRR27021982",
    `5708_3` = "SRR27021971", `5708_4` = "SRR27021960"
)

source <- read.delim(
    file.path(data_dir, "Proccessed_Data.txt"),
    check.names = FALSE,
    stringsAsFactors = FALSE
)
source$source_cell_id <- rownames(source)
source$source_group <- sub("(_[^_]+)$", "", source$source_cell_id)
source$barcode_token <- sub("^[^_]+_[^_]+_", "", source$source_cell_id)
source$run_id <- unname(group_to_run[source$source_group])
source$plavisca_cell_id <- paste0(sub("^SRR27021", "", source$run_id), "_", source$barcode_token)

stopifnot(
    nrow(source) == 80024L,
    !anyDuplicated(source$source_cell_id),
    !anyDuplicated(source$plavisca_cell_id),
    !anyNA(source$run_id),
    identical(as.character(source$Sample), source$source_group)
)

# Publication Supplemental Data 1 is the primary authority for retained blood
# library identities, days, infection design and published retained-cell counts.
supp_raw <- read_excel(
    file.path(audit_dir, "supplemental_data_1.xlsx"),
    sheet = "Sheet1",
    .name_repair = "minimal"
)
names(supp_raw)[1] <- "library_id"
names(supp_raw) <- make.unique(names(supp_raw), sep = "__duplicate_")
supp <- supp_raw |>
    filter(Library == "scRNA-seq", !grepl("^SPZ", library_id)) |>
    transmute(
        library_id,
        run_id = unname(group_to_run[library_id]),
        animal_id = sub("_.*$", "", library_id),
        day_post_infection_authoritative = as.integer(`Days Post Infection`),
        infection_design_authoritative = case_when(
            Condition == "Successive infection" ~ "Consecutive infection",
            TRUE ~ Condition
        ),
        inoculum_or_lineage_authoritative = Strain,
        publication_total_reads = as.numeric(`Total Reads (w/ barcode)`),
        publication_cells_gt1000 = as.integer(`Cells with >1000 reads within genes`),
        publication_chesson_cells = as.integer(`Chesson Cells`),
        publication_nih1993_cells = as.integer(`NIH 1993 Cells`)
    )

stopifnot(
    nrow(supp) == 24L,
    !anyDuplicated(supp$library_id),
    !anyDuplicated(supp$run_id),
    !anyNA(supp$run_id),
    sum(supp$publication_cells_gt1000) == 80024L
)

ena <- read.delim(
    file.path(audit_dir, "ena_PRJNA1047651_runs.tsv"),
    check.names = FALSE,
    stringsAsFactors = FALSE
)
stopifnot(nrow(ena) == 43L, !anyDuplicated(ena$run_accession))

# Flatten NCBI BioSample attributes without treating registry fields as
# experimental truth where they conflict with publication evidence.
biosample_doc <- read_xml(file.path(audit_dir, "ncbi_biosamples.xml"))
biosamples <- lapply(xml_find_all(biosample_doc, ".//BioSample"), function(node) {
    ids <- xml_find_all(node, ".//Id")
    attrs <- xml_find_all(node, ".//Attribute")
    get_attr <- function(harmonized_name) {
        hit <- attrs[xml_attr(attrs, "harmonized_name") == harmonized_name]
        if (length(hit) == 0L) NA_character_ else xml_text(hit[[1]])
    }
    data.frame(
        sample_accession = xml_attr(node, "accession"),
        biosample_title = xml_text(xml_find_first(node, "./Description/Title")),
        registry_organism = xml_text(xml_find_first(node, "./Description/Organism/OrganismName")),
        registry_geo_loc_name = get_attr("geo_loc_name"),
        registry_host = get_attr("host"),
        registry_isolate = get_attr("isolate"),
        registry_strain = get_attr("strain"),
        registry_tissue = get_attr("tissue"),
        stringsAsFactors = FALSE
    )
}) |>
    bind_rows()

study <- readRDS(file.path(pre_dir, "hazzard2024.rds"))
final_obj <- readRDS(file.path(app_dir, "cleaned_dataset.rds"))
final <- final_obj$mr_data |>
    as.data.frame() |>
    filter(study_label == "Hazzard2024")

stopifnot(
    ncol(study) == 80024L,
    nrow(final) == 80024L,
    !anyDuplicated(colnames(study)),
    !anyDuplicated(rownames(final)),
    setequal(colnames(study), rownames(final)),
    setequal(source$plavisca_cell_id, colnames(study))
)

study_meta <- study@meta.data
study_meta$plavisca_cell_id <- rownames(study_meta)
final$plavisca_cell_id <- rownames(final)

shared_fields <- intersect(names(study_meta), names(final))
shared_fields <- setdiff(shared_fields, "plavisca_cell_id")
final_aligned <- final[match(study_meta$plavisca_cell_id, final$plavisca_cell_id), , drop = FALSE]
metadata_concordance <- lapply(shared_fields, function(field) {
    a <- as.character(study_meta[[field]])
    b <- as.character(final_aligned[[field]])
    same <- (is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & a == b)
    data.frame(
        field,
        compared_cells = length(a),
        equal_cells = sum(same),
        unequal_cells = sum(!same),
        stringsAsFactors = FALSE
    )
}) |>
    bind_rows()
write_tsv(metadata_concordance, "study_to_deployed_metadata_concordance.tsv")

cell_lineage <- source |>
    select(
        source_cell_id, source_group, barcode_token, run_id,
        plavisca_cell_id, source_genotype = Genotype,
        source_condition = Condition, source_pseudotime = pseudotime,
        source_chesson_alleles = Chesson, source_nih_alleles = NIH
    ) |>
    left_join(supp, by = c("source_group" = "library_id", "run_id")) |>
    mutate(
        study_position = match(plavisca_cell_id, colnames(study)),
        deployed_position = match(plavisca_cell_id, rownames(final)),
        study_present = !is.na(study_position),
        deployed_present = !is.na(deployed_position)
    ) |>
    left_join(
        study_meta |>
            transmute(
                plavisca_cell_id,
                deployed_run_id = as.character(run_id),
                deployed_host_id = as.character(host_id),
                deployed_strain = as.character(strain),
                deployed_day_post_infection = as.integer(day_post_infection),
                deployed_geographic_location = as.character(geographic_location),
                deployed_host_species = as.character(host_species),
                deployed_technology = as.character(sc_technology),
                deployed_sequencer = as.character(sequencer),
                deployed_biological_replicate = as.character(biological_replicate)
            ),
        by = "plavisca_cell_id"
    ) |>
    mutate(
        expected_library_id = source_group,
        expected_animal_id = animal_id,
        host_library_id_error = deployed_host_id != expected_library_id,
        day_error = deployed_day_post_infection != day_post_infection_authoritative,
        source_genotype_disagrees_deployed_strain = case_when(
            is.na(source_genotype) ~ NA,
            source_genotype == "Chesson" ~ !grepl("Chesson", deployed_strain),
            source_genotype == "NIH" ~ !grepl("NIH", deployed_strain),
            source_genotype == "Outcross" ~ !grepl("Chesson", deployed_strain) |
                !grepl("NIH", deployed_strain),
            TRUE ~ NA
        )
    )

stopifnot(
    all(cell_lineage$study_present),
    all(cell_lineage$deployed_present),
    sum(cell_lineage$host_library_id_error) == 19810L,
    sum(cell_lineage$source_genotype_disagrees_deployed_strain, na.rm = TRUE) == 16321L
)
write_tsv(cell_lineage, "cell_lineage.tsv")

source_counts <- source |>
    count(source_group, Genotype, Condition, name = "source_cells")

run_table <- supp |>
    left_join(
        ena |>
            select(
                run_id = run_accession, experiment_accession, sample_accession,
                study_accession, secondary_study_accession, sample_alias,
                experiment_alias, instrument_model, library_strategy,
                library_source, library_selection, scientific_name
            ),
        by = "run_id"
    ) |>
    left_join(biosamples, by = "sample_accession") |>
    left_join(
        source |>
            count(source_group, run_id, name = "deposited_processed_cells"),
        by = c("library_id" = "source_group", "run_id")
    ) |>
    left_join(
        study_meta |>
            count(
                run_id,
                deployed_host_id = host_id,
                deployed_strain = strain,
                deployed_day = day_post_infection,
                deployed_geography = geographic_location,
                deployed_host_species = host_species,
                deployed_technology = sc_technology,
                deployed_sequencer = sequencer,
                deployed_biological_replicate = biological_replicate,
                name = "deployed_cells"
            ),
        by = "run_id"
    ) |>
    mutate(
        host_library_id_correct = deployed_host_id == library_id,
        day_correct = as.integer(deployed_day) == day_post_infection_authoritative,
        cell_count_correct = publication_cells_gt1000 == deposited_processed_cells &
            deposited_processed_cells == deployed_cells,
        retained_by_plavisca = TRUE,
        sample_material = "P. vivax-infected Saimiri boliviensis blood",
        capture_technology_authoritative = "10X Genomics 3-prime scRNA-seq; kit 1000268",
        sequencing_platform_authoritative = "Illumina NovaSeq",
        host_species_authoritative = "Saimiri boliviensis",
        host_sex_authoritative = "female",
        treatment_authority_note = "No drug-treatment arm described; infection design and splenectomy are separate concepts",
        parasite_provenance_note = case_when(
            inoculum_or_lineage_authoritative == "coinfection" ~
                "Descended from Chesson plus NIH-1993-F3 mosquito-stage coinfection; includes parental and outcross genotypes",
            grepl("Chesson", inoculum_or_lineage_authoritative) &
                grepl("NIH", inoculum_or_lineage_authoritative) ~
                "Chesson historical lineage (Papua New Guinea) plus NIH-1993-F3; exact NIH-1993 geographic origin not established here",
            grepl("Chesson", inoculum_or_lineage_authoritative) ~
                "Chesson historical lineage; Papua New Guinea origin supported",
            TRUE ~
                "NIH-1993-F3 lineage; closely related to but distinct from Salvador-I; exact geographic origin not established here"
        ),
        host_provenance_note = "Female Saimiri boliviensis from NIH-approved sources; donor/birth geography not reported",
        experimental_site_note = "NIAID/NIH animal protocol; do not infer animal origin from institutional location",
        sequencing_site_note = "Not stated in publication or Supplemental Data 1"
    ) |>
    arrange(run_id)
write_tsv(run_table, "authoritative_run_table.tsv")
write_tsv(source_counts, "source_genotype_condition_counts.tsv")

imported_runs <- c(
    paste0("SRR27021", 955:956), paste0("SRR27021", 959:972),
    paste0("SRR27021", 975:976), paste0("SRR27021", 980:987),
    paste0("SRR27021", 991:994)
)
run_membership <- ena |>
    transmute(
        run_id = run_accession, experiment_accession, sample_accession,
        sample_alias, scientific_name, instrument_model,
        deposited = TRUE,
        imported_by_plavisca_chunk = run_id %in% imported_runs,
        retained_by_plavisca = run_id %in% supp$run_id,
        retained_cells = as.integer(table(factor(source$run_id, levels = run_id)))
    ) |>
    arrange(run_id)
stopifnot(
    sum(run_membership$deposited) == 43L,
    sum(run_membership$imported_by_plavisca_chunk) == 30L,
    sum(run_membership$retained_by_plavisca) == 24L,
    sum(run_membership$retained_cells) == 80024L
)
write_tsv(run_membership, "all_deposited_run_membership.tsv")

discrepancies <- bind_rows(
    cell_lineage |>
        filter(host_library_id_error) |>
        count(
            issue = "wrong host/library ID", run_id,
            observed = deployed_host_id, expected = expected_library_id,
            name = "affected_cells"
        ),
    cell_lineage |>
        filter(source_genotype_disagrees_deployed_strain %in% TRUE) |>
        count(
            issue = "deployed strain excludes deposited cell genotype",
            run_id, observed = deployed_strain, expected = source_genotype,
            name = "affected_cells"
        ),
    cell_lineage |>
        filter(day_error) |>
        count(
            issue = "wrong day post infection", run_id,
            observed = as.character(deployed_day_post_infection),
            expected = as.character(day_post_infection_authoritative),
            name = "affected_cells"
        ),
    study_meta |>
        count(
            issue = "misspelled technology label", run_id,
            observed = sc_technology,
            expected = "10X Genomics 3-prime scRNA-seq; kit 1000268",
            name = "affected_cells"
        ),
    cell_lineage |>
        filter(is.na(deployed_biological_replicate) |
            deployed_biological_replicate == "NA") |>
        count(
            issue = "missing biological replicate despite known animal",
            run_id, observed = deployed_biological_replicate,
            expected = animal_id,
            name = "affected_cells"
        )
) |>
    arrange(issue, run_id, observed, expected)
write_tsv(discrepancies, "metadata_discrepancies.tsv")

hpi_stage <- final |>
    count(
        run_id, host_id, strain, day_post_infection,
        hour_post_invasion, life_cycle_stage, pred_gametocyte,
        name = "cells"
    ) |>
    arrange(run_id, hour_post_invasion, life_cycle_stage, pred_gametocyte)
write_tsv(hpi_stage, "hpi_lifecycle_impact.tsv")

source_vs_lifecycle <- source |>
    transmute(
        plavisca_cell_id,
        author_pseudotime_group = PseudoGroup,
        author_pseudotime = pseudotime
    ) |>
    left_join(
        final |>
            transmute(
                plavisca_cell_id,
                deployed_hpi = hour_post_invasion,
                deployed_lifecycle = life_cycle_stage,
                deployed_cluster_interpretation = pred_gametocyte
            ),
        by = "plavisca_cell_id"
    ) |>
    count(
        author_pseudotime_group, deployed_hpi, deployed_lifecycle,
        deployed_cluster_interpretation, name = "cells"
    ) |>
    arrange(
        author_pseudotime_group, deployed_hpi, deployed_lifecycle,
        deployed_cluster_interpretation
    )
write_tsv(source_vs_lifecycle, "source_vs_deployed_lifecycle.tsv")

stage_summary <- final |>
    summarise(
        deployed_cells = n(),
        hpi_nonmissing = sum(!is.na(hour_post_invasion)),
        hpi_distinct_labels = n_distinct(hour_post_invasion, na.rm = TRUE),
        final_female_gametocyte = sum(life_cycle_stage == "Female gametocyte", na.rm = TRUE),
        final_male_gametocyte = sum(life_cycle_stage == "Male gametocyte", na.rm = TRUE),
        cluster_female_gametocyte = sum(pred_gametocyte == "Female gametocyte", na.rm = TRUE),
        cluster_male_gametocyte = sum(pred_gametocyte == "Male gametocyte", na.rm = TRUE),
        author_sexual_group = sum(source$PseudoGroup == "Sexual", na.rm = TRUE),
        author_sexual_not_final_female = sum(
            source$PseudoGroup == "Sexual" &
                final$life_cycle_stage[match(source$plavisca_cell_id, rownames(final))] !=
                    "Female gametocyte",
            na.rm = TRUE
        ),
        final_female_not_author_sexual = sum(
            source$PseudoGroup != "Sexual" &
                final$life_cycle_stage[match(source$plavisca_cell_id, rownames(final))] ==
                    "Female gametocyte",
            na.rm = TRUE
        )
    )
write_tsv(stage_summary, "annotation_summary.tsv")

# Exact sparse count-profile fingerprints within the study object. Hash matches
# are verified elementwise before they are reported as exact duplicates.
counts <- GetAssayData(study, assay = "RNA", layer = "counts")
fingerprints <- vapply(seq_len(ncol(counts)), function(i) {
    lo <- counts@p[[i]] + 1L
    hi <- counts@p[[i + 1L]]
    if (lo > hi) {
        digest(list(integer(), numeric()), algo = "sha256")
    } else {
        digest(list(counts@i[lo:hi], counts@x[lo:hi]), algo = "sha256")
    }
}, character(1))
candidate_groups <- split(seq_along(fingerprints), fingerprints)
candidate_groups <- candidate_groups[lengths(candidate_groups) > 1L]
duplicate_rows <- list()
if (length(candidate_groups) > 0L) {
    for (hash in names(candidate_groups)) {
        idx <- candidate_groups[[hash]]
        anchor <- idx[[1]]
        exact <- idx[vapply(idx, function(j) {
            isTRUE(all(counts[, anchor, drop = FALSE] == counts[, j, drop = FALSE]))
        }, logical(1))]
        if (length(exact) > 1L) {
            duplicate_rows[[length(duplicate_rows) + 1L]] <- data.frame(
                hash_sha256 = hash,
                cell_id = colnames(counts)[exact],
                run_id = study_meta$run_id[match(colnames(counts)[exact], rownames(study_meta))],
                group_size = length(exact),
                stringsAsFactors = FALSE
            )
        }
    }
}
exact_duplicates <- bind_rows(duplicate_rows)
if (nrow(exact_duplicates) == 0L) {
    exact_duplicates <- data.frame(
        hash_sha256 = character(), cell_id = character(),
        run_id = character(), group_size = integer()
    )
}
write_tsv(exact_duplicates, "exact_count_duplicates.tsv")

inventory <- data.frame(
    metric = c(
        "deposited_runs", "chunk_imported_runs", "retained_runs",
        "publication_cells_gt1000", "deposited_processed_rows",
        "study_rds_cells", "deployed_cells", "unique_source_ids",
        "unique_plavisca_ids", "study_features", "study_count_nnz",
        "host_id_error_cells", "strain_error_cells", "day_error_cells",
        "technology_typo_cells", "missing_biological_replicate_cells",
        "hpi_nonmissing_cells", "female_gametocyte_cells",
        "exact_duplicate_profile_groups", "exact_duplicate_profile_cells"
    ),
    value = c(
        nrow(ena), length(unique(imported_runs)), nrow(supp),
        sum(supp$publication_cells_gt1000), nrow(source),
        ncol(study), nrow(final), n_distinct(source$source_cell_id),
        n_distinct(source$plavisca_cell_id), nrow(study), length(counts@x),
        sum(cell_lineage$host_library_id_error),
        sum(cell_lineage$source_genotype_disagrees_deployed_strain, na.rm = TRUE),
        sum(cell_lineage$day_error), ncol(study),
        sum(is.na(study_meta$biological_replicate) |
            as.character(study_meta$biological_replicate) == "NA"),
        sum(!is.na(final$hour_post_invasion)),
        sum(final$life_cycle_stage == "Female gametocyte", na.rm = TRUE),
        n_distinct(exact_duplicates$hash_sha256), nrow(exact_duplicates)
    )
)
write_tsv(inventory, "inventory_summary.tsv")

validation <- c(
    "PASS source Proccessed_Data.txt MD5 matches Zenodo-advertised checksum",
    "PASS publication Supplemental Data 1 checksum fixed and readable",
    "PASS 43 unique deposited runs in ENA snapshot",
    "PASS 30 unique runs imported by chunk scripts",
    "PASS 24 unique blood libraries retained by deposited processed-cell list",
    "PASS publication retained-cell counts sum to 80024",
    "PASS deposited processed table has 80024 unique source IDs",
    "PASS transformed source IDs have 80024 unique values",
    "PASS transformed source IDs equal the study-RDS cell set",
    "PASS study-RDS IDs equal deployed Hazzard2024 row IDs",
    paste0("CHECK deployed day mismatch cells=", sum(cell_lineage$day_error),
        "; see metadata_discrepancies.tsv and authoritative_run_table.tsv"),
    "PASS wrong host/library IDs affect exactly 19810 deployed cells",
    "PASS deposited-genotype/deployed-strain conflict affects exactly 16321 cells",
    paste0("PASS exact count duplicate inspection completed; groups=", n_distinct(exact_duplicates$hash_sha256),
        "; cells=", nrow(exact_duplicates)),
    "LIMIT source count-table ZIP could not be retrieved from Zenodo during this audit; raw/source-count equivalence is not established",
    "LIMIT central integrated Seurat object is absent; Harmony and SingleR were not rerun",
    "LIMIT exact executed historical build remains unproven without a build manifest"
)
writeLines(validation, file.path(audit_dir, "validation.txt"))

rm(counts, study, final_obj)
gc()

manifest_files <- list.files(audit_dir, full.names = TRUE)
manifest_files <- manifest_files[
    file.info(manifest_files)$isdir %in% FALSE &
        basename(manifest_files) != "evidence_manifest.tsv"
]
evidence_manifest <- data.frame(
    filename = basename(manifest_files),
    bytes = unname(file.info(manifest_files)$size),
    sha256 = vapply(manifest_files, function(path) {
        digest(path, algo = "sha256", file = TRUE, serialize = FALSE)
    }, character(1)),
    stringsAsFactors = FALSE
) |>
    arrange(filename)
write_tsv(evidence_manifest, "evidence_manifest.tsv")
