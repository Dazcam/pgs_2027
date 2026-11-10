#--------------------------------------------------------------------------------------
#
#    GWAS formatting functions
#
#--------------------------------------------------------------------------------------
#
# Functions: gwas_read()          read only the columns needed from a raw GWAS file
#            gwas_to_beta()       effect column to BETA (log for OR)
#            gwas_is_ambiguous()  flag strand-ambiguous SNPs (A/T, C/G)
#            gwas_snp_n()         per-SNP N by method (column, case_control, fixed, auto)
#            gwas_qc()            sequential QC with SNP counts per step
#            gwas_n_summary()     single N for PRS-CS from per-SNP N
#
# Packages:  caller must have readr and dplyr available (pkg:: calls used throughout)
#
# Used in:   gwas_format.R
#
#--------------------------------------------------------------------------------------

gwas_read <- function(path, delim, cols) {

  # cols: character vector of raw column names to keep
  # All columns are read as text; the caller converts numeric fields
  col_spec <- readr::cols(.default = readr::col_character())

  if (delim == "tab") {
    df <- readr::read_tsv(path, col_types = col_spec, col_select = dplyr::all_of(cols),
                          progress = FALSE)
  } else if (delim == "space") {
    # read_table has no col_select in older readr, so read all columns then select
    df <- readr::read_table(path, col_types = col_spec, progress = FALSE) |>
      dplyr::select(dplyr::all_of(cols))
  } else {
    stop("delim must be 'tab' or 'space', not '", delim, "'")
  }

  df
}

gwas_to_beta <- function(effect, effect_type) {

  if (effect_type == "beta") return(effect)
  if (effect_type == "or")   return(log(effect))
  stop("effect_type must be 'beta' or 'or', not '", effect_type, "'")
}

gwas_is_ambiguous <- function(a1, a2) {

  pair <- paste0(toupper(a1), toupper(a2))
  pair %in% c("AT", "TA", "CG", "GC")
}

gwas_snp_n <- function(df, n_cfg, trait_type, auto_maf) {

  # df has standard columns (se, and freq for auto) plus any raw N columns named in n_cfg
  # Returns per-SNP N; for auto, NA outside the MAF window
  method <- n_cfg[["method"]]

  if (method == "column") {
    return(as.numeric(df[[n_cfg[["col"]]]]) * n_cfg[["scale"]])
  }

  if (method == "case_control") {
    cases    <- as.numeric(df[[n_cfg[["cases"]]]])
    controls <- as.numeric(df[[n_cfg[["controls"]]]])
    return(4 / (1 / cases + 1 / controls))
  }

  if (method == "fixed") {
    return(rep(as.numeric(n_cfg[["value"]]), nrow(df)))
  }

  if (method == "auto") {
    if (!"freq" %in% colnames(df)) stop("n: auto needs an allele frequency column (cols: freq)")
    k   <- switch(trait_type, binary = 4, continuous = 1,
                  stop("trait_type must be 'binary' or 'continuous', not '", trait_type, "'"))
    maf <- pmin(df$freq, 1 - df$freq)
    n   <- k / (2 * df$freq * (1 - df$freq) * df$se^2)
    n[!(maf >= auto_maf[1] & maf <= auto_maf[2])] <- NA
    return(n)
  }

  stop("n method must be column, case_control, fixed or auto, not '", method, "'")
}

gwas_qc <- function(df, info_min = NULL, drop_ambiguous = TRUE, drop_duplicates = TRUE,
                    qc_cols = character(0), n_min_frac = NULL) {

  # df has standard columns snp, a1, a2, beta, se, p, n (and info if present)
  # Returns list(df = filtered data, qc = SNPs remaining after each step)
  steps <- list()
  log_step <- function(step, d) steps[[length(steps) + 1]] <<- tibble::tibble(step = step, snps = nrow(d))

  log_step("input", df)

  df <- df |>
    dplyr::filter(!is.na(snp), !is.na(beta), is.finite(beta), !is.na(se), se > 0,
                  !is.na(p), p > 0, p <= 1)
  log_step("valid effect, SE and P", df)

  for (col in qc_cols) {
    df <- df[!is.na(df[[col]]) & as.numeric(df[[col]]) == 0, ]
  }
  if (length(qc_cols) > 0) log_step(paste0(paste(qc_cols, collapse = " & "), " == 0"), df)

  if (!is.null(info_min) && "info" %in% colnames(df)) {
    df <- df |> dplyr::filter(!is.na(info), info >= info_min)
    log_step(paste0("INFO >= ", info_min), df)
  }

  if (isTRUE(drop_ambiguous)) {
    df <- df[!gwas_is_ambiguous(df$a1, df$a2), ]
    log_step("not strand-ambiguous", df)
  }

  if (isTRUE(drop_duplicates)) {
    dup <- df$snp[duplicated(df$snp)]
    df  <- df[!df$snp %in% dup, ]
    log_step("unique SNP ID", df)
  }

  if (!is.null(n_min_frac)) {
    n_cut <- n_min_frac * stats::quantile(df$n, 0.9, na.rm = TRUE, names = FALSE)
    df    <- df |> dplyr::filter(!is.na(n), n >= n_cut)
    log_step(paste0("N >= ", n_min_frac, " x 90th percentile"), df)
  }

  if (nrow(df) == 0) stop("No SNPs left after QC")

  list(df = df, qc = dplyr::bind_rows(steps))
}

gwas_n_summary <- function(n, n_summary = "median") {

  n <- n[!is.na(n)]
  if (length(n) == 0) stop("No per-SNP N values to summarise")
  switch(n_summary,
         median = stats::median(n),
         mean   = mean(n),
         stop("n_summary must be 'median' or 'mean', not '", n_summary, "'"))
}

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
