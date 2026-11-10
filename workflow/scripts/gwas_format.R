#--------------------------------------------------------------------------------------
#
#    Format GWAS summary statistics for PRS-CS
#
#--------------------------------------------------------------------------------------
#
# Pipeline:  02PREP-GWAS | Rule: format_gwas (one job per GWAS)
#            Upstream:   00SETUP get_gwas_sumstats
#            Downstream: 03PRS (PRS-CS)
#
# Purpose:   - Read the columns set in the GWAS block of the config
#            - Standardise to SNP, A1 (effect allele), A2, BETA (log(OR) for OR), SE
#            - Per-SNP N by the configured method
#            - QC: valid values, GenomicSEM error/warning, INFO, strand-ambiguous SNPs,
#              duplicate IDs, optional low-N filter
#            - Single N for PRS-CS (median or mean of per-SNP N after QC)
#
# Inputs:    raw_gwas  GWAS file from ../resources/gwas/<gwas>/
#
# Outputs:   sumstats  PRS-CS input: SNP A1 A2 BETA SE (tab-separated)
#            n_gwas    single value for PRS-CS --n_gwas
#            qc        SNPs remaining after each QC step, and the per-SNP N summary
#
# Notes:     PRS-CS reads the sumstats by column position, so the column order is fixed
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

message("\n\nFormatting GWAS summary statistics for PRS-CS ...")

##  Load packages, functions and variables  -------------------------------------------
library(tidyverse)

# Load functions
source(file.path(snakemake@scriptdir, "gwas_functions.R"))

# Input and output paths
gwas_in      <- snakemake@input[["raw_gwas"]]
sumstats_out <- snakemake@output[["sumstats"]]
n_gwas_out   <- snakemake@output[["n_gwas"]]
qc_out       <- snakemake@output[["qc"]]

# GWAS settings
gwas        <- snakemake@wildcards[["gwas"]]
gwas_cfg    <- snakemake@params[["gwas_cfg"]]
cols        <- unlist(gwas_cfg[["cols"]])
delim       <- gwas_cfg[["delim"]]
effect_type <- gwas_cfg[["effect_type"]]
trait_type  <- gwas_cfg[["trait_type"]]
qc_cols     <- unlist(gwas_cfg[["qc_cols"]])
n_cfg       <- gwas_cfg[["n"]]

# QC settings
info_min        <- snakemake@params[["info_min"]]
drop_ambiguous  <- isTRUE(snakemake@params[["drop_ambiguous"]])
drop_duplicates <- isTRUE(snakemake@params[["drop_duplicates"]])
n_min_frac      <- snakemake@params[["n_min_frac"]]
n_summary       <- snakemake@params[["n_summary"]]
auto_maf        <- unlist(snakemake@params[["auto_maf"]])

if (is.null(qc_cols)) qc_cols <- character(0)
if (!is.null(n_min_frac) && n_cfg[["method"]] %in% c("auto", "fixed")) {
  message("n_min_frac is ignored for n method '", n_cfg[["method"]], "'")
  n_min_frac <- NULL
}

n_desc <- switch(n_cfg[["method"]],
                 column       = paste0(n_cfg[["col"]], " x ", n_cfg[["scale"]]),
                 case_control = paste0("4 / (1/", n_cfg[["cases"]], " + 1/", n_cfg[["controls"]], ")"),
                 fixed        = as.character(n_cfg[["value"]]),
                 auto         = paste0("auto (", trait_type, ", MAF ", paste(auto_maf, collapse = "-"), ")"))

# Print what each variable is set to
message("\nVariables")
message("============================")
tibble(
  variable = c("gwas", "gwas_in", "sumstats_out", "n_gwas_out", "qc_out", "delim", "cols",
               "effect_type", "trait_type", "qc_cols", "n", "info_min", "drop_ambiguous",
               "drop_duplicates", "n_min_frac", "n_summary"),
  value    = c(gwas, gwas_in, sumstats_out, n_gwas_out, qc_out, delim,
               paste(names(cols), cols, sep = "=", collapse = ", "),
               effect_type, trait_type, paste(qc_cols, collapse = ", "), n_desc,
               ifelse(is.null(info_min), "off", info_min), drop_ambiguous, drop_duplicates,
               ifelse(is.null(n_min_frac), "off", n_min_frac), n_summary)) |>
  knitr::kable(format = "simple") |>
  cat(sep = "\n")
message("============================\n")

##  Read and standardise  -------------------------------------------------------------
n_raw_cols <- unlist(n_cfg[c("col", "cases", "controls")])
read_cols  <- unique(c(cols, qc_cols, n_raw_cols))

raw <- gwas_read(gwas_in, delim, read_cols)
message("SNPs read: ", nrow(raw))

gwas_df <- tibble(
  snp  = raw[[cols[["snp"]]]],
  a1   = toupper(raw[[cols[["a1"]]]]),
  a2   = toupper(raw[[cols[["a2"]]]]),
  beta = gwas_to_beta(as.numeric(raw[[cols[["effect"]]]]), effect_type),
  se   = as.numeric(raw[[cols[["se"]]]]),
  p    = as.numeric(raw[[cols[["p"]]]]))

if (!is.null(gwas_cfg[["cols"]][["freq"]])) gwas_df$freq <- as.numeric(raw[[cols[["freq"]]]])
if (!is.null(gwas_cfg[["cols"]][["info"]])) gwas_df$info <- as.numeric(raw[[cols[["info"]]]])
for (col in qc_cols) gwas_df[[col]] <- raw[[col]]

gwas_df$n <- gwas_snp_n(bind_cols(gwas_df, raw[setdiff(n_raw_cols, colnames(gwas_df))]),
                        n_cfg, trait_type, auto_maf)
rm(raw)

##  QC  -------------------------------------------------------------------------------
qc_res <- gwas_qc(gwas_df, info_min = info_min, drop_ambiguous = drop_ambiguous,
                  drop_duplicates = drop_duplicates, qc_cols = qc_cols, n_min_frac = n_min_frac)
gwas_df <- qc_res$df

qc_tbl <- qc_res$qc |>
  mutate(removed = lag(snps) - snps)

message("QC steps:")
qc_tbl |>
  knitr::kable(format = "simple") |>
  cat(sep = "\n")

##  N for PRS-CS  ---------------------------------------------------------------------
n_gwas <- round(gwas_n_summary(gwas_df$n, n_summary))
n_stats <- tibble(
  stat  = c("n_min", "n_median", "n_mean", "n_max", "n_snps_used", "n_gwas"),
  value = c(min(gwas_df$n, na.rm = TRUE), median(gwas_df$n, na.rm = TRUE),
            mean(gwas_df$n, na.rm = TRUE), max(gwas_df$n, na.rm = TRUE),
            sum(!is.na(gwas_df$n)), n_gwas))

message("\nPer-SNP N after QC (", n_desc, "):")
n_stats |>
  mutate(value = format(round(value), big.mark = ",")) |>
  knitr::kable(format = "simple") |>
  cat(sep = "\n")

##  Write outputs  --------------------------------------------------------------------
dir.create(dirname(sumstats_out), recursive = TRUE, showWarnings = FALSE)

gwas_df |>
  transmute(SNP = snp, A1 = a1, A2 = a2, BETA = beta, SE = se) |>
  write_tsv(sumstats_out)

write_lines(n_gwas, n_gwas_out)

bind_rows(
  qc_tbl |> transmute(gwas = gwas, type = "qc_step", name = step, value = snps),
  n_stats |> transmute(gwas = gwas, type = "n", name = stat, value = value)) |>
  write_tsv(qc_out)

message("\nSNPs written: ", nrow(gwas_df))
message("n_gwas: ", n_gwas)
message("Done.")

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
