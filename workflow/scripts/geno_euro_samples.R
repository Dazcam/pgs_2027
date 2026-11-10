#--------------------------------------------------------------------------------------
#
#    Extract European-ancestry sample IDs
#
#--------------------------------------------------------------------------------------
#
# Pipeline:  01GENOTYPES | Rule: euro_samples
#            Upstream:   sample sheet (copied in by hand)
#            Downstream: subset_vcf
#
# Purpose:   - Read the sample sheet (QC'd samples only)
#            - Keep samples of the configured ancestry
#            - Write one sample ID per line for bcftools --samples-file
#
# Inputs:    sample_sheet  Excel sample sheet with sample ID and ancestry columns
#
# Outputs:   euro_samples  Sample IDs, one per line, no header
#
#--------------------------------------------------------------------------------------

## Set up logging for Snakemake
if (exists("snakemake")) {
  log_smk <- function() {
    if (exists("snakemake") && length(snakemake@log) != 0) {
      log <- file(snakemake@log[[1]], open = "wt")
      sink(log, append = TRUE)
      sink(log, append = TRUE, type = "message")
    }
  }
  log_smk()
}

message("\n\nExtracting European-ancestry sample IDs ...")

##  Load packages and variables  ------------------------------------------------------
library(tidyverse)
library(readxl)

# Input and output paths
sheet_in    <- snakemake@input[[1]]
samples_out <- snakemake@output[[1]]

# Sheet settings
sheet         <- snakemake@params[["sheet"]]
id_col        <- snakemake@params[["id_col"]]
ancestry_col  <- snakemake@params[["ancestry_col"]]
keep_ancestry <- snakemake@params[["ancestry"]]

# Print what each variable is set to
message("\nVariables")
message("============================")
tibble(
  variable = c("sheet_in", "samples_out", "sheet", "id_col", "ancestry_col", "keep_ancestry"),
  value    = c(sheet_in, samples_out, sheet, id_col, ancestry_col, keep_ancestry)) |>
  knitr::kable(format = "simple") |>
  cat(sep = "\n")
message("============================\n")

##  Read sheet and filter  ------------------------------------------------------------
# Read every column as text so numeric IDs are not written as e.g. 1102.0
# .env$ below: the sheet's own ancestry column would otherwise mask the variable
anns <- read_excel(sheet_in, sheet = sheet, col_types = "text")

missing_cols <- setdiff(c(id_col, ancestry_col), colnames(anns))
if (length(missing_cols) > 0) {
  stop("Columns not found in sheet: ", paste(missing_cols, collapse = ", "))
}

anns <- anns |>
  mutate(id = str_trim(.data[[id_col]]),
         anc = str_trim(.data[[ancestry_col]]))

if (anyNA(anns$id) || any(duplicated(anns$id))) {
  stop("Sample IDs in the sheet are missing or duplicated")
}

message("Samples in sheet: ", nrow(anns))
message("Samples by ancestry:")
count(anns, anc) |>
  knitr::kable(format = "simple") |>
  cat(sep = "\n")

samples <- anns |>
  dplyr::filter(anc == .env$keep_ancestry) |>
  pull(id)

if (length(samples) == 0) {
  stop("No samples with ancestry '", keep_ancestry, "' found in the sheet")
}

message("\nSamples kept: ", length(samples))

##  Write output  ---------------------------------------------------------------------
dir.create(dirname(samples_out), recursive = TRUE, showWarnings = FALSE)
write_lines(samples, samples_out)

message("Done.")

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
