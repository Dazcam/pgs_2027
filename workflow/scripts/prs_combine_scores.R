#--------------------------------------------------------------------------------------
#
#    Combine PGS across GWAS with genotype PCs, sex and PCW
#
#--------------------------------------------------------------------------------------
#
# Pipeline:  03PRS | Rule: combine_scores
#            Upstream:   score (one .sscore per GWAS), 01GENOTYPES genotype_pcs, euro_samples
#            Downstream: deliverable for the collaborators; report
#
# Purpose:   - Read the PLINK 2 score for each GWAS (sum of dosage x posterior effect)
#            - Standardise each score across donors (mean 0, SD 1)
#            - Add genotype PCs and sex and PCW from the sample sheet
#            - Summarise SNPs used and raw score distribution per GWAS
#
# Inputs:    sscores  PLINK 2 .sscore files, one per GWAS (same order as params gwas)
#            pcs      PLINK 1.9 .eigenvec (no header: FID IID PC1 ...)
#            samples  European-ancestry sample IDs used for the genotypes
#            sheet    sample sheet (sex, PCW)
#
# Outputs:   scores   one row per donor: iid, pgs_<gwas> (standardised), pgs_raw_<gwas>,
#                     PC1-PCn, sex, PCW
#            summary  per GWAS: SNPs used, raw score mean and SD
#
# Notes:     Scores are not adjusted for covariates; that is left to the downstream model
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

message("\n\nCombining PGS across GWAS ...")

##  Load packages, functions and variables  -------------------------------------------
library(tidyverse)
library(readxl)

# Load functions
source(file.path(snakemake@scriptdir, "prs_functions.R"))

# Input and output paths
sscores_in  <- unlist(snakemake@input[["sscores"]])
pcs_in      <- snakemake@input[["pcs"]]
samples_in  <- snakemake@input[["samples"]]
sheet_in    <- snakemake@input[["sheet"]]
scores_out  <- snakemake@output[["scores"]]
summary_out <- snakemake@output[["summary"]]

# Settings
gwas_ids <- unlist(snakemake@params[["gwas"]])
sheet    <- snakemake@params[["sheet"]]
id_col   <- snakemake@params[["id_col"]]
sex_col  <- snakemake@params[["sex_col"]]
pcw_col  <- snakemake@params[["pcw_col"]]

# Print what each variable is set to
message("\nVariables")
message("============================")
tibble(
  variable = c("sscores_in", "pcs_in", "samples_in", "sheet_in", "scores_out", "summary_out",
               "gwas_ids", "sheet", "id_col", "sex_col", "pcw_col"),
  value    = c(paste(basename(sscores_in), collapse = ", "), pcs_in, samples_in, sheet_in,
               scores_out, summary_out, paste(gwas_ids, collapse = ", "), sheet, id_col,
               sex_col, pcw_col)) |>
  knitr::kable(format = "simple") |>
  cat(sep = "\n")
message("============================\n")

if (length(sscores_in) != length(gwas_ids)) {
  stop("Number of score files (", length(sscores_in), ") does not match number of GWAS (",
       length(gwas_ids), ")")
}

##  Read scores  ----------------------------------------------------------------------
samples <- read_lines(samples_in)
message("Donors expected: ", length(samples))

scores_long <- map2(sscores_in, gwas_ids, \(path, gwas) {
  sc <- prs_read_sscore(path)
  if (!setequal(sc$iid, samples) || nrow(sc) != length(samples)) {
    stop(gwas, ": donors in the score file do not match the European sample list")
  }
  sc |> mutate(gwas = gwas)
}) |>
  list_rbind()

summary_tbl <- scores_long |>
  group_by(gwas) |>
  summarise(snps_used = first(allele_ct) / 2,
            allele_ct_varies = n_distinct(allele_ct) > 1,
            raw_mean = mean(score_sum),
            raw_sd   = sd(score_sum),
            .groups = "drop") |>
  mutate(gwas = factor(gwas, levels = gwas_ids)) |>
  arrange(gwas) |>
  mutate(gwas = as.character(gwas))

message("Per-GWAS summary:")
summary_tbl |>
  knitr::kable(format = "simple", digits = 6) |>
  cat(sep = "\n")

if (any(summary_tbl$allele_ct_varies)) {
  message("\nNote: ALLELE_CT differs between donors for: ",
          paste(summary_tbl$gwas[summary_tbl$allele_ct_varies], collapse = ", "))
}

scores_wide <- scores_long |>
  group_by(gwas) |>
  mutate(pgs = prs_standardise(score_sum)) |>
  ungroup() |>
  select(iid, gwas, pgs, pgs_raw = score_sum) |>
  pivot_wider(names_from = gwas, values_from = c(pgs, pgs_raw), names_vary = "slowest")

##  Covariates  -----------------------------------------------------------------------
pcs <- read_table(pcs_in, col_names = FALSE, col_types = cols(.default = col_character()),
                  progress = FALSE)
colnames(pcs) <- c("fid", "iid", paste0("PC", seq_len(ncol(pcs) - 2)))
pcs <- pcs |>
  select(-fid) |>
  mutate(across(starts_with("PC"), as.numeric))

anns <- read_excel(sheet_in, sheet = sheet, col_types = "text") |>
  transmute(iid = str_trim(.data[[id_col]]),
            sex = .data[[sex_col]],
            PCW = as.numeric(.data[[pcw_col]]))

out <- scores_wide |>
  left_join(pcs, by = "iid") |>
  left_join(anns, by = "iid") |>
  arrange(match(iid, samples))

if (anyNA(out$PC1)) stop("Genotype PCs missing for some donors")
if (anyNA(out$sex) || anyNA(out$PCW)) message("Note: sex or PCW missing for some donors")

message("\nDonors written: ", nrow(out), "; columns: ", ncol(out))

##  Write outputs  --------------------------------------------------------------------
dir.create(dirname(scores_out), recursive = TRUE, showWarnings = FALSE)
write_tsv(out, scores_out)
write_tsv(summary_tbl |> select(-allele_ct_varies), summary_out)

message("Done.")

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
