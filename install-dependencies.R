## Install everything viroscape needs, then the package itself.
## Run once from the project directory:  Rscript install-dependencies.R

cran <- c(
  # core
  "ape", "phangorn", "tibble", "dplyr", "tidyr", "purrr", "rlang", "cli",
  "jsonlite", "ggplot2", "scales",
  # modelling extras
  "lme4", "car",
  # retrieval
  "rentrez",
  # app
  "shiny", "bslib", "DT",
  # pipeline, docs, tests
  "targets", "testthat", "knitr", "rmarkdown", "roxygen2", "devtools", "patchwork"
)

missing <- cran[!cran %in% rownames(installed.packages())]
if (length(missing)) {
  message("Installing from CRAN: ", paste(missing, collapse = ", "))
  install.packages(missing, repos = "https://cloud.r-project.org")
} else {
  message("All CRAN dependencies already present.")
}

## Optional: MUSCLE alignment via Bioconductor, to reproduce the published
## alignment exactly. Everything else works without it - viroscape falls back to
## DECIPHER, an external mafft/muscle binary, or pre-aligned input.
if (interactive()) {
  ans <- readline("Install the Bioconductor aligners (msa, Biostrings, DECIPHER)? [y/N] ")
  if (tolower(trimws(ans)) %in% c("y", "yes")) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager", repos = "https://cloud.r-project.org")
    }
    BiocManager::install(c("msa", "Biostrings", "DECIPHER"), ask = FALSE)
  }
} else {
  message("\nOptional Bioconductor aligners (skipped in a non-interactive run):")
  message('  BiocManager::install(c("msa", "Biostrings", "DECIPHER"))')
}

## Install viroscape itself from this directory
if (requireNamespace("devtools", quietly = TRUE)) {
  devtools::install(".", dependencies = FALSE, upgrade = "never")
} else {
  install.packages(".", repos = NULL, type = "source")
}

message("\nDone. Try:  library(viroscape); example_analysis()")
