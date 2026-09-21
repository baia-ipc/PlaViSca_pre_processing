# Environment smoke test - confirms the reproducible Pixi environment can
# load every package the preprocessing/annotation scripts depend on,
# including the Bioconductor "data" packages whose conda post-link.sh hook
# pixi does not execute (see scripts/setup/ensure_bioc_data_packages.sh).
#
# Run with:
#   env -u R_LIBS_USER pixi run --manifest-path scripts/pixi.toml Rscript --vanilla scripts/tests/test_environment.R
# or:
#   env -u R_LIBS_USER pixi run --manifest-path scripts/pixi.toml test-env

required_packages <- c(
  "tidyverse",
  "Seurat",
  "SeuratObject",
  "Matrix",
  "DropletUtils",
  "scater",
  "SingleCellExperiment",
  "rtracklayer",
  "GenomeInfoDb",
  "GenomeInfoDbData",
  "SingleR",
  "scCustomize"
)

failures <- character(0)
for (pkg in required_packages) {
  ok <- requireNamespace(pkg, quietly = TRUE)
  cat(sprintf("[%s] %s\n", if (ok) "OK" else "MISSING", pkg))
  if (!ok) failures <- c(failures, pkg)
}

if (length(failures) > 0) {
  stop(sprintf(
    "Environment smoke test FAILED: could not load %s. See scripts/setup/ensure_bioc_data_packages.sh for the known GenomeInfoDbData/post-link.sh issue.",
    paste(failures, collapse = ", ")
  ))
}

cat("Environment smoke test PASSED: all required packages load.\n")
