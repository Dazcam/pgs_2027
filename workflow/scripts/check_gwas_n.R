#--------------------------------------------------------------------------------------
#
#    Check the sample size (N) of each GWAS
#
#--------------------------------------------------------------------------------------
#
# Standalone check, not part of the Snakemake pipeline. Run from workflow/:
#   singularity exec ../resources/containers/prs.sif Rscript scripts/check_gwas_n.R [gwas ...]
# With no arguments, every GWAS in the config is checked.
#
# Purpose:   For each raw GWAS file, report
#            - ## header lines that mention sample sizes, and the column names
#            - per-SNP N as set in the config (method and scale factor):
#              min, median, max, and SNPs below 0.3 and 0.5 of the 90th percentile,
#              with the chromosomes those SNPs are on
#            - where the file has case and control columns, 4/(1/cases + 1/controls)
#              and its ratio to the configured N, to check the scale factor
#            - an independent estimate from SE and allele frequency (Mallard et al.,
#              as n: auto), and its ratio to the configured N
#
# Notes:     - All rows are used, except that the estimate from SE uses SNPs with
#              INFO >= info_min (where there is an INFO column) and MAF in auto_maf
#            - A ratio near 1 supports the configured N; near 0.5 or 2 suggests a
#              wrong scale factor. The SE-based estimate is approximate (SE
#              corrections, imputation and frequency source all shift it), so treat
#              differences of 10-20% as normal
#            - Output: ../results/02PREP-GWAS/check_gwas_n.tsv, and the log on screen
#
#--------------------------------------------------------------------------------------

message("\n\nChecking GWAS sample sizes ...")

##  Load packages, functions and variables  -------------------------------------------
library(tidyverse)

source("scripts/gwas_functions.R")

# R's YAML reader treats the key n: as FALSE (YAML 1.1 boolean); keep y/n as text
yes_no   <- list("bool#no"  = \(x) if (x %in% c("n", "N")) x else FALSE,
                 "bool#yes" = \(x) if (x %in% c("y", "Y")) x else TRUE)
config   <- yaml::read_yaml("../config/config.yaml", handlers = yes_no)
gwas_dir <- config$setup$get_gwas_sumstats$dir
prep_cfg <- config$prep_gwas$format_gwas
info_min <- prep_cfg$filters$info$min
auto_maf <- unlist(prep_cfg$auto_maf)
out_file <- "../results/02PREP-GWAS/check_gwas_n.tsv"

args     <- commandArgs(trailingOnly = TRUE)
gwas_ids <- if (length(args) > 0) args else names(config$gwas)

chr_cols  <- c("CHR", "CHROM", "#CHROM", "chr", "Chr")
case_cols <- c("NCAS", "Nca", "Ncases", "N_CAS", "Nca_eff")
ctrl_cols <- c("NCON", "Nco", "Ncontrols", "N_CON", "Nco_eff")

check_gwas <- function(gwas) {

  cfg  <- config$gwas[[gwas]]
  path <- file.path(str_replace(gwas_dir, fixed("{gwas}"), gwas), cfg$file)
  message("\n============================================================")
  message(gwas, ": ", path)
  message("============================================================")
  if (!file.exists(path)) {
    message("File not found, skipped")
    return(NULL)
  }

  # Header lines about sample size, and column names
  skip     <- gwas_meta_lines(path)
  meta     <- read_lines(path, n_max = skip)
  n_meta   <- meta[str_detect(meta, regex("ncas|ncon|nca|nco|neff|ncase|ncontrol|neffective|sample", ignore_case = TRUE))]
  n_meta   <- str_trunc(n_meta, 200)
  header   <- read_lines(path, skip = skip, n_max = 1) |> str_split(if (cfg$delim == "tab") "\t" else "\\s+") |> unlist()
  header   <- header[header != ""]

  message("\nHeader lines mentioning sample size:")
  if (length(n_meta) == 0) message("  none") else walk(n_meta, \(x) message("  ", x))
  message("\nColumns: ", paste(header, collapse = " "))

  # Columns to read
  cols      <- unlist(cfg$cols)
  n_cfg     <- cfg$n_calc
  n_raw     <- unlist(n_cfg[c("col", "cases", "controls")])
  chr_col   <- intersect(chr_cols, header)[1]
  case_col  <- intersect(case_cols, header)[1]
  ctrl_col  <- intersect(ctrl_cols, header)[1]
  read_cols <- unique(na.omit(c(cols, n_raw, chr_col, case_col, ctrl_col)))

  raw <- gwas_read(path, cfg$delim, read_cols)
  message("\nSNPs in file: ", format(nrow(raw), big.mark = ","))

  df <- tibble(
    se   = as.numeric(raw[[cols[["se"]]]]),
    freq = if (!is.null(cfg$cols$freq)) as.numeric(raw[[cols[["freq"]]]]) else NA_real_,
    info = if (!is.null(cfg$cols$info)) as.numeric(raw[[cols[["info"]]]]) else NA_real_,
    chr  = if (!is.na(chr_col)) raw[[chr_col]] else NA_character_)
  df$n <- gwas_snp_n(bind_cols(df, raw[n_raw]), n_cfg, cfg$trait_type, auto_maf)

  # Configured per-SNP N
  n_desc <- switch(n_cfg$method,
                   column       = paste0(n_cfg$col, " x ", n_cfg$scale),
                   case_control = paste0("4/(1/", n_cfg$cases, " + 1/", n_cfg$controls, ")"),
                   fixed        = paste0("fixed ", n_cfg$value),
                   auto         = paste0("auto (", cfg$trait_type, ")"))
  n_med  <- median(df$n, na.rm = TRUE)
  n_ref  <- quantile(df$n, 0.9, na.rm = TRUE, names = FALSE)
  below  <- map_int(c(0.3, 0.5), \(f) sum(df$n < f * n_ref, na.rm = TRUE))

  message("\nConfigured N (", n_desc, "):")
  message("  min ", format(round(min(df$n, na.rm = TRUE)), big.mark = ","),
          " | median ", format(round(n_med), big.mark = ","),
          " | max ", format(round(max(df$n, na.rm = TRUE)), big.mark = ","))
  message("  SNPs below 0.3 / 0.5 x 90th percentile: ",
          format(below[1], big.mark = ","), " / ", format(below[2], big.mark = ","))

  if (below[2] > 0 && !is.na(chr_col)) {
    low_chr <- df |> filter(n < 0.5 * n_ref) |> count(chr) |>
      mutate(chr = fct_inseq(factor(chr), ordered = FALSE)) |> arrange(chr)
    message("  Chromosomes of SNPs below 0.5: ",
            paste0(low_chr$chr, ":", low_chr$n, collapse = " "))
  }

  # Case and control columns, if present
  cc_med <- NA_real_
  if (!is.na(case_col) && !is.na(ctrl_col)) {
    cases    <- as.numeric(raw[[case_col]])
    controls <- as.numeric(raw[[ctrl_col]])
    cc_med    <- median(4 / (1 / cases + 1 / controls), na.rm = TRUE)
    message("\nCase/control columns (", case_col, ", ", ctrl_col, "):")
    message("  median cases ", format(round(median(cases, na.rm = TRUE)), big.mark = ","),
            " | median controls ", format(round(median(controls, na.rm = TRUE)), big.mark = ","))
    message("  median 4/(1/cases + 1/controls): ", format(round(cc_med), big.mark = ","),
            " (ratio to configured N ", round(cc_med / n_med, 2), ")")
  }

  # Independent estimate from SE and allele frequency
  se_med <- NA_real_
  if (!all(is.na(df$freq))) {
    k   <- if (cfg$trait_type == "binary") 4 else 1
    maf <- pmin(df$freq, 1 - df$freq)
    keep <- maf >= auto_maf[1] & maf <= auto_maf[2] & !is.na(df$se) & df$se > 0
    if (!is.null(info_min) && !all(is.na(df$info))) keep <- keep & !is.na(df$info) & df$info >= info_min
    se_med <- median(k / (2 * df$freq[keep] * (1 - df$freq[keep]) * df$se[keep]^2), na.rm = TRUE)
    message("\nEstimate from SE and frequency (k = ", k, ", MAF ", paste(auto_maf, collapse = "-"),
            ", ", format(sum(keep, na.rm = TRUE), big.mark = ","), " SNPs): ",
            format(round(se_med), big.mark = ","), " (ratio to configured N ", round(se_med / n_med, 2), ")")
  } else {
    message("\nEstimate from SE and frequency: not possible (no frequency column)")
  }

  tibble(gwas = gwas, n_method = n_desc, snps = nrow(raw),
         n_min = min(df$n, na.rm = TRUE), n_median = n_med, n_max = max(df$n, na.rm = TRUE),
         below_0.3 = below[1], below_0.5 = below[2],
         cc_median = cc_med, cc_ratio = cc_med / n_med,
         se_median = se_med, se_ratio = se_med / n_med)
}

##  Run  ------------------------------------------------------------------------------
res <- map(gwas_ids, check_gwas) |> list_rbind()

message("\n\nSummary")
message("============================")
res |>
  mutate(across(c(n_min, n_median, n_max, cc_median, se_median), \(x) round(x)),
         across(ends_with("ratio"), \(x) round(x, 2))) |>
  knitr::kable(format = "simple") |>
  cat(sep = "\n")

dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
write_tsv(res, out_file)
message("\nWritten: ", out_file)
message("Done.")

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
