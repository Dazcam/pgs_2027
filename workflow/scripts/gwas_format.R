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
#            - Flag SNPs PRS-CS would use: rsID with compatible alleles in the LD
#              reference and the genotype .bim (PRS-CS's own matching)
#            - QC: rows with an invalid SNP, BETA, SE or P are always removed. Every
#              other filter (biallelic SNVs, GenomicSEM error/warning, INFO,
#              strand-ambiguous SNPs, duplicate IDs, low N) is counted, in the whole
#              file and among the SNPs PRS-CS would use, and removes SNPs only if
#              switched on in the config
#            - Single N for PRS-CS (median or mean of per-SNP N after QC)
#
# Inputs:    raw_gwas  GWAS file from ../resources/gwas/<gwas>/
#            ref       PRS-CS LD reference (snpinfo_ukbb_hm3 is read)
#            bed       genotypes (the .bim at the same prefix is read)
#
# Outputs:   sumstats  PRS-CS input: SNP A1 A2 BETA SE (tab-separated)
#            n_gwas    single value for PRS-CS --n_gwas
#            qc        per-SNP N summary
#            qc_steps  per QC step: applied, SNPs failing and remaining, in the whole
#                      file and among the SNPs PRS-CS would use
#
# Notes:     - PRS-CS reads the sumstats by column position, so the column order is fixed
#            - The SNPs remaining among those PRS-CS would use, after the applied
#              steps, should equal the SNPs PRS-CS reports in common (snps_used)
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
qc_steps_out <- snakemake@output[["qc_steps"]]
ref_snpinfo  <- snakemake@params[["ref_snpinfo"]]
bim_file     <- snakemake@params[["bim"]]

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
filters   <- snakemake@params[["filters"]]
n_summary <- snakemake@params[["n_summary"]]
auto_maf  <- unlist(snakemake@params[["auto_maf"]])

if (is.null(qc_cols)) qc_cols <- character(0)

# The low-N filter needs per-SNP N that varies between SNPs
low_n_ok <- !n_cfg[["method"]] %in% c("auto", "fixed")
if (!low_n_ok) message("Low-N filter not evaluated for n method '", n_cfg[["method"]], "'")

filter_desc <- names(filters) |>
  map_chr(\(f) {
    x     <- filters[[f]]
    extra <- setdiff(names(x), "apply")
    paste0(f, ": ", ifelse(isTRUE(x[["apply"]]), "applied", "counted only"),
           if (length(extra) > 0) paste0(" (", paste(extra, unlist(x[extra]), sep = " ", collapse = ", "), ")") else "")
  }) |>
  paste(collapse = "; ")

n_desc <- switch(n_cfg[["method"]],
                 column       = paste0(n_cfg[["col"]], " x ", n_cfg[["scale"]]),
                 case_control = paste0("4 / (1/", n_cfg[["cases"]], " + 1/", n_cfg[["controls"]], ")"),
                 fixed        = as.character(n_cfg[["value"]]),
                 auto         = paste0("auto (", trait_type, ", MAF ", paste(auto_maf, collapse = "-"), ")"))

# Print what each variable is set to
message("\nVariables")
message("============================")
tibble(
  variable = c("gwas", "gwas_in", "ref_snpinfo", "bim_file", "sumstats_out", "n_gwas_out",
               "qc_out", "qc_steps_out", "delim", "cols", "effect_type", "trait_type",
               "qc_cols", "n", "filters", "n_summary"),
  value    = c(gwas, gwas_in, ref_snpinfo, bim_file, sumstats_out, n_gwas_out, qc_out,
               qc_steps_out, delim, paste(names(cols), cols, sep = "=", collapse = ", "),
               effect_type, trait_type, paste(qc_cols, collapse = ", "), n_desc,
               filter_desc, n_summary)) |>
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

##  SNPs PRS-CS would use  ------------------------------------------------------------
ref <- read_table(ref_snpinfo, col_types = cols(.default = col_character()), progress = FALSE) |>
  transmute(chr = CHR, snp = SNP, a1 = A1, a2 = A2)

bim <- read_tsv(bim_file, col_names = FALSE, col_select = c(1, 2, 5, 6),
                col_types = cols(.default = col_character()), progress = FALSE) |>
  set_names(c("chr", "snp", "a1", "a2")) |>
  filter(snp %in% ref$snp)

gwas_df$overlap <- gwas_in_overlap(gwas_df$snp, gwas_df$a1, gwas_df$a2, ref, bim)
message("\nSNPs in the LD reference: ", nrow(ref),
        "; of these in the genotype .bim: ", n_distinct(bim$snp),
        "; GWAS rows matching both: ", sum(gwas_df$overlap))
rm(ref, bim)

##  QC  -------------------------------------------------------------------------------
qc_res  <- gwas_qc(gwas_df, filters = filters, qc_cols = qc_cols, low_n_ok = low_n_ok)
gwas_df <- qc_res$df
qc_tbl  <- qc_res$qc

message("\nQC steps (failing: SNPs failing the step on its own; remaining: after the applied steps;")
message("overlap: among the SNPs PRS-CS would use):")
qc_tbl |>
  knitr::kable(format = "simple", format.args = list(big.mark = ",")) |>
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

n_stats |>
  transmute(gwas = gwas, type = "n", name = stat, value = value) |>
  write_tsv(qc_out)

qc_tbl |>
  mutate(gwas = gwas, .before = 1) |>
  write_tsv(qc_steps_out)

message("\nSNPs written: ", nrow(gwas_df), " (of these PRS-CS would use: ", sum(gwas_df$overlap), ")")
message("n_gwas: ", n_gwas)
message("Done.")

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
