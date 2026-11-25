#--------------------------------------------------------------------------------------
#
#    Pseudobulk counts for one cell type, for sharing
#
#--------------------------------------------------------------------------------------
#
# Pipeline:  04DATA-SHARING | Rule: share_pseudobulk
#            Upstream:   eQTL_study_2025 02SCANPY pseudobulk CSVs; share_scores
#            Downstream: share_readme, share_tar
#
# Purpose:   - Read the Scanpy pseudobulk CSV (donors x genes; first column unnamed)
#            - Keep the donors in the shared PGS file, in the same order, matched
#              exactly; stop if any is missing, unless it is listed for this cell type
#              in data_sharing: allowed_missing
#            - Name the first column iid; gene IDs unchanged
#            - Write tab-separated, gzipped, under the shared cell type name
#
# Inputs:    csv      <cell_type>_pseudobulk.csv from the eQTL pipeline
#            scores   shared pgs_scores.tsv (donor order)
#
# Outputs:   counts   <label>.tsv.gz
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

message("\n\nPreparing pseudobulk counts for sharing ...")

##  Load packages, functions and variables  -------------------------------------------
library(tidyverse)

source(file.path(snakemake@scriptdir, "share_functions.R"))

# Input and output paths
csv_in     <- snakemake@input[["csv"]]
scores_in  <- snakemake@input[["scores"]]
counts_out <- snakemake@output[["counts"]]

# Settings
cell_type <- snakemake@params[["cell_type"]]
allowed   <- as.character(unlist(snakemake@params[["allowed_missing"]]))
label     <- snakemake@wildcards[["label"]]

# Print what each variable is set to
message("\nVariables")
message("============================")
tibble(
  variable = c("csv_in", "scores_in", "counts_out", "cell_type", "label", "allowed"),
  value    = c(csv_in, scores_in, counts_out, cell_type, label,
               ifelse(length(allowed) == 0, "none", paste(allowed, collapse = ", ")))) |>
  knitr::kable(format = "simple") |>
  cat(sep = "\n")
message("============================\n")

##  Read  -----------------------------------------------------------------------------
# First column (donor ID) as text, counts as numbers; name_repair = "minimal" keeps
# gene IDs exactly as written, even if repeated
n_cols <- length(str_split_1(read_lines(csv_in, n_max = 1), ","))
counts <- read_csv(csv_in, name_repair = "minimal", progress = FALSE,
                   col_types = paste0("c", strrep("d", n_cols - 1)))
colnames(counts)[1] <- "iid"
donors <- read_tsv(scores_in, col_types = cols(iid = col_character()), col_select = iid)$iid

message("Donors in CSV: ", nrow(counts), "; genes: ", ncol(counts) - 1)
dup_genes <- sum(duplicated(colnames(counts)[-1]))
if (dup_genes > 0) message("Note: ", dup_genes, " repeated gene IDs (kept as in the CSV)")

##  Match donors  ---------------------------------------------------------------------
expected <- share_expected_donors(donors, counts$iid, allowed)
rows     <- share_match_donors(counts$iid, expected, paste0(cell_type, " pseudobulk"))
dropped  <- setdiff(counts$iid, donors)
absent   <- setdiff(donors, expected)
message("Donors kept: ", length(rows), "; dropped (not in the PGS file): ", length(dropped))
if (length(absent) > 0) {
  message("PGS donors absent from this cell type (allowed in the config): ",
          paste(absent, collapse = ", "))
}

out <- counts[rows, ]
if (anyNA(select(out, -iid))) stop("Missing counts for kept donors")

##  Write  ----------------------------------------------------------------------------
dir.create(dirname(counts_out), recursive = TRUE, showWarnings = FALSE)
write_tsv(out, counts_out)
message("\nWritten: ", counts_out, " (", nrow(out), " donors x ", ncol(out) - 1, " genes)")
message("Done.")

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
