#--------------------------------------------------------------------------------------
#
#    PGS functions
#
#--------------------------------------------------------------------------------------
#
# Functions: prs_read_sscore()   read a PLINK 2 .sscore file (with or without FID)
#            prs_standardise()   z-score across donors
#
# Packages:  caller must have readr and dplyr available (pkg:: calls used throughout)
#
# Used in:   prs_combine_scores.R
#
#--------------------------------------------------------------------------------------

prs_read_sscore <- function(path) {

  # PLINK 2 writes #IID (no FID) or #FID IID; returns iid, allele_ct, score_sum
  df <- readr::read_tsv(path, col_types = readr::cols(.default = readr::col_character()),
                        progress = FALSE)
  colnames(df) <- sub("^#", "", colnames(df))

  needed <- c("IID", "ALLELE_CT", "SCORE1_SUM")
  missing_cols <- setdiff(needed, colnames(df))
  if (length(missing_cols) > 0) {
    stop("Columns not found in ", path, ": ", paste(missing_cols, collapse = ", "),
         " (was --score run with cols=+scoresums?)")
  }

  tibble::tibble(iid       = df$IID,
                 allele_ct = as.numeric(df$ALLELE_CT),
                 score_sum = as.numeric(df$SCORE1_SUM))
}

prs_standardise <- function(x) {

  s <- stats::sd(x)
  if (is.na(s) || s == 0) stop("Cannot standardise: score has no variation across donors")
  (x - mean(x)) / s
}

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
