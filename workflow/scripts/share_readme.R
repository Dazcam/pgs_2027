#--------------------------------------------------------------------------------------
#
#    README for the shared files
#
#--------------------------------------------------------------------------------------
#
# Pipeline:  04DATA-SHARING | Rule: share_readme
#            Upstream:   share_scores, share_pseudobulk, 03PRS combine_scores summary and
#                        set_n_gwas (primary N set)
#            Downstream: share_tar
#
# Purpose:   - Check every shared file has the same donors in the same order, apart from
#              donors allowed to be absent from a cell type (data_sharing: allowed_missing)
#            - Write README.txt: files, cell types, covariates (sex counts, PCW range,
#              PC settings) and PGS settings with a per-GWAS table (file, n_gwas and its
#              source, SNPs in the score), from the config and this run's outputs
#
# Inputs:    scores   shared pgs_scores.tsv
#            counts   shared pseudobulk files, in config cell type order
#            summary  03PRS pgs_summary for the primary N set (snps_used per GWAS)
#            n_gwas   n_gwas passed to PRS-CS per GWAS (primary N set, GWAS order)
#            qc       format_gwas per-SNP N summary per GWAS (GWAS order)
#
# Outputs:   readme   README.txt
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

message("\n\nWriting the README for the shared files ...")

##  Load packages, functions and variables  -------------------------------------------
library(tidyverse)

source(file.path(snakemake@scriptdir, "share_functions.R"))

# Input and output paths
scores_in  <- snakemake@input[["scores"]]
counts_in  <- unlist(snakemake@input[["counts"]])
summary_in <- snakemake@input[["summary"]]
n_gwas_in  <- unlist(snakemake@input[["n_gwas"]])
qc_in      <- unlist(snakemake@input[["qc"]])
readme_out <- snakemake@output[["readme"]]

# Settings
labels      <- unlist(snakemake@params[["labels"]])
allowed     <- snakemake@params[["allowed_missing"]]     # per label: donors allowed to be absent
gwas_ids    <- unlist(snakemake@params[["gwas"]])
gwas_files  <- unlist(snakemake@params[["gwas_files"]])
n_sources   <- unlist(snakemake@params[["n_sources"]])
prscs       <- snakemake@params[["prscs"]]
ld_ref      <- snakemake@params[["ld_reference"]]
pca         <- snakemake@params[["pca"]]
sex_values  <- unlist(snakemake@params[["sex_values"]])
n_summary   <- snakemake@params[["n_summary"]]
title       <- snakemake@params[["title"]]
code        <- snakemake@params[["code"]]

# Print what each variable is set to
message("\nVariables")
message("============================")
tibble(
  variable = c("scores_in", "counts_in", "summary_in", "n_gwas_in", "readme_out", "labels",
               "gwas_ids", "ld_ref", "title", "code"),
  value    = c(scores_in, paste(basename(counts_in), collapse = ", "), summary_in,
               paste(basename(n_gwas_in), collapse = ", "), readme_out,
               paste(labels, collapse = ", "), paste(gwas_ids, collapse = ", "), ld_ref,
               title, code)) |>
  knitr::kable(format = "simple") |>
  cat(sep = "\n")
message("============================\n")

##  Check donors across files  --------------------------------------------------------
scores <- read_tsv(scores_in, col_types = cols(iid = col_character(), sex = col_character(),
                                                .default = col_double()))
n_donors <- nrow(scores)

absent <- map2(counts_in, labels, \(f, lab) {
  # Base R: the gene header line is too long for readr's connection buffer
  ids      <- sub("\t.*$", "", readLines(f))[-1]
  ok       <- as.character(unlist(allowed[[lab]]))
  expected <- share_expected_donors(scores$iid, ids, ok)
  if (!identical(ids, expected)) {
    stop(lab, ": donors differ from pgs_scores.tsv in identity or order")
  }
  tibble(label = lab, iid = setdiff(scores$iid, ids))
}) |>
  list_rbind()
message("All ", length(counts_in), " pseudobulk files have the PGS donors in the same order",
        if (nrow(absent) > 0) paste0("; absent as allowed: ", nrow(absent), " donor-cell type pairs") else "")

# README sentence on donors: absences listed per donor
donor_lines <- if (nrow(absent) == 0) {
  paste0("All files contain the same ", n_donors, " donors in the same order.")
} else {
  c(paste0("All files contain the same ", n_donors, " donors in the same order, except that these donors"),
    "failed QC for some cell types and are absent from those files:",
    absent |>
      mutate(label = factor(label, levels = labels)) |>
      arrange(label) |>
      summarise(types = paste(label, collapse = ", "), .by = iid) |>
      transmute(line = paste0("  ", iid, ": ", types)) |>
      pull(line))
}

##  Values for the README  ------------------------------------------------------------
pc_cols  <- str_subset(colnames(scores), "^PC[0-9]+$")
pc_range <- paste0(pc_cols[1], "-", pc_cols[length(pc_cols)])
sex_n    <- map_chr(sex_values, \(s) paste0(s, " (n = ", sum(scores$sex == s), ")"))
pcw      <- format(range(scores$PCW), trim = TRUE)
prune    <- str_split_1(str_squish(pca$indep_pairwise), " ")

summary_tbl <- read_tsv(summary_in, show_col_types = FALSE)

# Per-SNP N range after QC: an N column that varies between SNPs was summarised
# (prep_gwas: n_summary) to the single value PRS-CS takes
n_varies <- map_lgl(qc_in, \(f) {
  n <- read_tsv(f, show_col_types = FALSE) |> filter(type == "n")
  n$value[n$name == "n_min"] != n$value[n$name == "n_max"]
})

gwas_tbl <- tibble(gwas = gwas_ids, file = gwas_files[gwas_ids], source = n_sources[gwas_ids],
                   n_gwas = map_dbl(n_gwas_in, \(f) as.numeric(read_lines(f)[1])),
                   varies = n_varies) |>
  mutate(source = if_else(source == "N column in file" & varies,
                          paste(n_summary, "of per-SNP N in file"), source)) |>
  select(-varies) |>
  left_join(summary_tbl |> select(gwas, snps_used), by = "gwas")
if (anyNA(gwas_tbl)) stop("Missing values in the GWAS table")

fmt_n <- \(x) format(x, big.mark = ",", trim = TRUE, scientific = FALSE)
cols  <- list(c("GWAS", gwas_tbl$gwas), c("File", gwas_tbl$file),
              c("n_gwas", fmt_n(gwas_tbl$n_gwas)), c("n_gwas source", gwas_tbl$source),
              c("SNPs in score", fmt_n(gwas_tbl$snps_used)))
widths <- map_int(cols, \(x) max(nchar(x)))
table_lines <- pmap_chr(cols, \(...) {
  x <- c(...)
  paste0("  ", paste(str_pad(x, widths + 3, side = "right"), collapse = "") |> str_trim("right"))
})

##  README  ---------------------------------------------------------------------------
readme <- c(
  str_replace(title, fixed("{n_donors}"), as.character(n_donors)),
  "",
  "Files",
  "  pgs_scores.tsv             One row per donor: iid; pgs_<gwas> and pgs_raw_<gwas> for each GWAS;",
  paste0("                             genotype PCs ", pc_range, "; sex; PCW"),
  "  pseudobulk/<cell_type>.tsv.gz",
  paste0("                             Raw pseudobulk counts, one file per cell type (",
         length(labels), "); one row per donor"),
  "                             (iid), one column per gene (mixed encoding, gene symbols or Ensembl IDs, as output by Scanpy)",
  "  README.txt                 This file",
  "  md5sums.txt                MD5 checksums of the files above; check with: md5sum -c md5sums.txt",
  "",
  donor_lines,
  "",
  "Cell types",
  share_cell_type_lines(labels),
  "",
  "Covariates",
  paste0("  sex           ", paste(sex_n, collapse = " or ")),
  paste0("  PCW           post-conception weeks, range ", pcw[1], "-", pcw[2]),
  paste0("  ", str_pad(pc_range, 14, side = "right"), "genotype principal components: ",
         pca$software, " PCA in the ", n_donors, " donors;"),
  paste0("                SNPs with MAF >= ", pca$maf, ", LD-pruned (window ", prune[1],
         ", step ", prune[2], ", r2 ", prune[3], ")"),
  "",
  "Polygenic scores",
  paste0("  Method        PRS-CS: a = ", prscs$a, ", b = ", prscs$b, ", phi = ",
         prscs$phi %||% "learnt from the data (auto)", "; ", prscs$n_iter,
         " MCMC iterations, ", prscs$n_burnin, " burn-in,"),
  paste0("                thinning ", prscs$thin, ", seed ", prscs$seed),
  paste0("  LD reference  ", ld_ref),
  "  pgs_raw_      sum of imputed dosage x posterior effect size over SNPs (PLINK 2, autosomes)",
  paste0("  pgs_          pgs_raw_ standardised across the ", n_donors,
         " donors: (score - mean) / SD"),
  "",
  table_lines,
  "",
  "  n_gwas: the single GWAS sample size passed to PRS-CS (--n_gwas); Supplementary Table 1 N is the",
  "  total sample size as published (cases + controls for case/control GWAS)",
  "",
  paste0("Code available at ", code))

##  Write  ----------------------------------------------------------------------------
dir.create(dirname(readme_out), recursive = TRUE, showWarnings = FALSE)
write_lines(readme, readme_out)
message("\nREADME:\n")
cat(readme, sep = "\n")
message("\nWritten: ", readme_out)
message("Done.")

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
