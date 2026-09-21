#!/usr/bin/env Rscript

# Base-R-only helper tests. The dplyr import is removed solely so these pure
# helper functions can be exercised despite the baseline Pixi post-link defect.
source_lines <- readLines("scripts/pipeline_lib.R")
source_lines <- source_lines[!grepl(
  "suppressPackageStartupMessages\\(library\\(dplyr\\)\\)", source_lines
)]
eval(parse(text = source_lines), envir = .GlobalEnv)

result <- data.frame(test = character(), observed = character(), expected = character(), status = character())
add <- function(test, observed, expected, status) {
  result <<- rbind(result, data.frame(test, observed, expected, status))
}
did_error <- function(expr) inherits(tryCatch(expr, error = identity), "error")

add("assert_keyed_join duplicate data keys", did_error(assert_keyed_join(c("a", "a"), "a")), "TRUE", "FAIL")
add("assert_keyed_join NA key", did_error(assert_keyed_join(NA_character_, NA_character_)), "TRUE", "FAIL")
factor_observed <- paste(map_broad_stage(factor(c("Sporozoite", "Ring"))), collapse = ";")
add("map_broad_stage factor input", factor_observed, "Sporozoite stage;Blood stage",
    ifelse(factor_observed == "Sporozoite stage;Blood stage", "PASS", "FAIL"))
character_observed <- paste(map_broad_stage(c("Sporozoite", "Ring")), collapse = ";")
add("map_broad_stage character input", character_observed, "Sporozoite stage;Blood stage",
    ifelse(character_observed == "Sporozoite stage;Blood stage", "PASS", "FAIL"))
eligibility <- c(TRUE, FALSE, TRUE)
cell_id <- c("c1", "c2", "c3")
indexed <- paste(eligibility[cell_id], collapse = ";")
add("unnamed eligibility indexed by cell ID", indexed, "TRUE;FALSE;TRUE", "FAIL")

write.table(result, "audit/phase1_code_review/helper_edge_case_results.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
