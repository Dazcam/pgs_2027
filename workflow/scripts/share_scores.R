#--------------------------------------------------------------------------------------
#
#    Check the PGS file against its sources and copy it for sharing
#
#--------------------------------------------------------------------------------------
#
# Pipeline:  04DATA-SHARING | Rule: share_scores
#            Upstream:   03PRS combine_scores (primary N set), 01GENOTYPES genotype_pcs,
#                        euro_samples
#            Downstream: share_pseudobulk, share_readme, share_tar
#
# Purpose:   - Donors: exactly the European sample list, no duplicates
#            - PC1-PCn: equal to the .eigenvec for every donor
#            - sex and PCW: equal to the sample sheet for every donor; sex is M or F;
#              PCW is never missing
#            - Write the file unchanged for the tarball
#
# Inputs:    scores   pgs_scores_<nset>_<tag>.tsv (iid, pgs_*, pgs_raw_*, PCs, sex, PCW)
#            pcs      PLINK 1.9 .eigenvec (no header: FID IID PC1 ...)
#            samples  European-ancestry sample IDs
#            sheet    sample sheet (sex, PCW)
#
# Outputs:   scores   pgs_scores.tsv for the tarball
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

message("\n\nChecking the PGS file for sharing ...")

##  Load packages, functions and variables  -------------------------------------------
library(tidyverse)
library(readxl)

source(file.path(snakemake@scriptdir, "share_functions.R"))

# Input and output paths
scores_in  <- snakemake@input[["scores"]]
pcs_in     <- snakemake@input[["pcs"]]
samples_in <- snakemake@input[["samples"]]
sheet_in   <- snakemake@input[["sheet"]]
scores_out <- snakemake@output[["scores"]]

# Settings
sheet   <- snakemake@params[["sheet"]]
id_col  <- snakemake@params[["id_col"]]
sex_col <- snakemake@params[["sex_col"]]
pcw_col <- snakemake@params[["pcw_col"]]
sex_ok  <- unlist(snakemake@params[["sex_values"]])

# Print what each variable is set to
message("\nVariables")
message("============================")
tibble(
  variable = c("scores_in", "pcs_in", "samples_in", "sheet_in", "scores_out", "sheet",
               "id_col", "sex_col", "pcw_col", "sex_ok"),
  value    = c(scores_in, pcs_in, samples_in, sheet_in, scores_out, sheet, id_col, sex_col,
               pcw_col, paste(sex_ok, collapse = ", "))) |>
  knitr::kable(format = "simple") |>
  cat(sep = "\n")
message("============================\n")

##  Read  -----------------------------------------------------------------------------
scores  <- read_tsv(scores_in, col_types = cols(iid = col_character(), sex = col_character(),
                                                .default = col_double()))
samples <- read_lines(samples_in)

pcs <- read_table(pcs_in, col_names = FALSE, col_types = cols(.default = col_character()),
                  progress = FALSE)
colnames(pcs) <- c("fid", "iid", paste0("PC", seq_len(ncol(pcs) - 2)))
pcs <- pcs |>
  select(-fid) |>
  mutate(across(starts_with("PC"), as.numeric))

anns <- read_excel(sheet_in, sheet = sheet, col_types = "text") |>
  transmute(iid = str_trim(.data[[id_col]]),
            sex = str_trim(.data[[sex_col]]),
            PCW = as.numeric(.data[[pcw_col]]))

message("Donors in PGS file: ", nrow(scores), "; in European sample list: ", length(samples))

##  Checks  ---------------------------------------------------------------------------
# Donors: exactly the sample list
share_match_donors(scores$iid, samples, "PGS file")
if (nrow(scores) != length(samples)) {
  stop("PGS file has ", nrow(scores), " donors; the European sample list has ", length(samples))
}

# Genotype PCs
pc_cols <- str_subset(colnames(scores), "^PC[0-9]+$")
pcs_src <- pcs[share_match_donors(pcs$iid, scores$iid, ".eigenvec"), ]
walk(pc_cols, \(pc) share_check_equal(scores[[pc]], pcs_src[[pc]], scores$iid, pc))
message("PCs checked against the .eigenvec: ", pc_cols[1], "-", pc_cols[length(pc_cols)])

# Sex and PCW
anns_src <- anns[share_match_donors(anns$iid, scores$iid, "sample sheet"), ]
share_check_equal(scores$sex, anns_src$sex, scores$iid, "sex")
share_check_equal(scores$PCW, anns_src$PCW, scores$iid, "PCW")
message("Sex and PCW checked against the sample sheet")

bad_sex <- scores$iid[is.na(scores$sex) | !scores$sex %in% sex_ok]
if (length(bad_sex) > 0) {
  stop("Sex not one of ", paste(sex_ok, collapse = "/"), " for: ", paste(bad_sex, collapse = ", "))
}
if (anyNA(scores$PCW)) {
  stop("PCW missing for: ", paste(scores$iid[is.na(scores$PCW)], collapse = ", "))
}

message("\nSex:")
count(scores, sex) |>
  knitr::kable(format = "simple") |>
  cat(sep = "\n")
message("PCW range: ", min(scores$PCW), "-", max(scores$PCW))

##  Write  ----------------------------------------------------------------------------
dir.create(dirname(scores_out), recursive = TRUE, showWarnings = FALSE)
file.copy(scores_in, scores_out, overwrite = TRUE)
message("\nWritten: ", scores_out)
message("Done.")

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
