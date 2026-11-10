library(testthat)

script_dir <- if (exists("snakemake")) snakemake@scriptdir else getwd()

source(file.path(script_dir, "../scripts/gwas_functions.R"))

test_dir(script_dir, reporter = "summary")
