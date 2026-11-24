#--------------------------------------------------------------------------------------
#
#    GWAS formatting functions
#
#--------------------------------------------------------------------------------------
#
# Functions: gwas_meta_lines()    count leading ## metadata lines (e.g. PGC VCF-style files)
#            gwas_read()          read only the columns needed from a raw GWAS file
#            gwas_to_beta()       effect column to BETA (log for OR)
#            gwas_is_snv()        flag biallelic single-base SNVs (A/C/G/T, alleles differ)
#            gwas_is_ambiguous()  flag strand-ambiguous SNPs (A/T, C/G)
#            gwas_snp_n()         per-SNP N by method (column, case_control, fixed, auto)
#            gwas_allele_key()    allele-pair key invariant to allele order and strand
#            gwas_in_overlap()    flag SNPs PRS-CS would use (LD reference and .bim)
#            gwas_qc()            QC: every filter counted, switched-on ones applied
#            gwas_n_summary()     single N for PRS-CS from per-SNP N
#
# Packages:  caller must have readr and dplyr available (pkg:: calls used throughout)
#
# Used in:   gwas_format.R
#
#--------------------------------------------------------------------------------------

gwas_meta_lines <- function(path, n_check = 1000) {

  # Number of lines at the top of the file that start with ##
  head_lines <- readr::read_lines(path, n_max = n_check, progress = FALSE)
  not_meta   <- which(!startsWith(head_lines, "##"))
  if (length(not_meta) == 0) stop("No header line found in the first ", n_check, " lines of ", path)
  not_meta[1] - 1
}

gwas_read <- function(path, delim, cols) {

  # cols: character vector of raw column names to keep
  # Leading ## metadata lines are skipped; a single # (e.g. #CHROM) stays in the header
  # All columns are read as text; the caller converts numeric fields
  if (!delim %in% c("tab", "space")) stop("delim must be 'tab' or 'space', not '", delim, "'")

  col_spec <- readr::cols(.default = readr::col_character())
  skip     <- gwas_meta_lines(path)

  if (delim == "tab") {
    df <- readr::read_tsv(path, col_types = col_spec, col_select = dplyr::all_of(cols),
                          skip = skip, progress = FALSE)
  } else {
    # read_table has no col_select in older readr, so read all columns then select
    df <- readr::read_table(path, col_types = col_spec, skip = skip, progress = FALSE) |>
      dplyr::select(dplyr::all_of(cols))
  }

  df
}

gwas_to_beta <- function(effect, effect_type) {

  if (effect_type == "beta") return(effect)
  if (effect_type == "or")   return(log(effect))
  stop("effect_type must be 'beta' or 'or', not '", effect_type, "'")
}

gwas_is_snv <- function(a1, a2) {

  # PRS-CS reads only SNPs whose alleles are single A, C, G or T; indels and
  # other allele codes are dropped there, so drop them here to match
  a1 <- toupper(a1)
  a2 <- toupper(a2)
  bases <- c("A", "C", "G", "T")
  a1 %in% bases & a2 %in% bases & a1 != a2
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

gwas_allele_key <- function(a1, a2) {

  # Key for an allele pair that is the same whichever allele is listed first and
  # on either strand (A/G, G/A, T/C and C/T all give "AG"), as PRS-CS matches
  # alleles. NA unless both alleles are single A/C/G/T bases and differ
  comp <- c(A = "T", C = "G", G = "C", T = "A")
  a1   <- toupper(a1)
  a2   <- toupper(a2)
  ok   <- a1 %in% names(comp) & a2 %in% names(comp) & a1 != a2

  sorted_pair <- function(x, y) ifelse(x < y, paste0(x, y), paste0(y, x))
  key <- pmin(sorted_pair(a1, a2), sorted_pair(unname(comp[a1]), unname(comp[a2])))
  key[!ok] <- NA_character_
  key
}

gwas_in_overlap <- function(snp, a1, a2, ref, bim) {

  # TRUE for GWAS rows PRS-CS would use: rsID with compatible alleles in both the
  # LD reference and the genotype .bim, on the same chromosome in those two
  # ref, bim: tibbles with chr, snp, a1, a2
  site_key <- function(d) paste(sub("^chr", "", d$chr), d$snp, gwas_allele_key(d$a1, d$a2))
  ref_ok   <- !is.na(gwas_allele_key(ref$a1, ref$a2))
  shared   <- intersect(site_key(ref)[ref_ok], site_key(bim))

  # Drop the chromosome: the GWAS is matched on rsID and alleles only
  shared_snp <- sub("^\\S+ ", "", shared)
  key <- gwas_allele_key(a1, a2)
  !is.na(key) & paste(snp, key) %in% shared_snp
}

gwas_qc <- function(df, filters = list(), qc_cols = character(0), low_n_ok = TRUE) {

  # df has standard columns snp, a1, a2, beta, se, p, n, plus info and qc_cols where
  # present, and overlap (TRUE if PRS-CS would use the SNP) if known
  # filters: named list, one entry per filter, each list(apply = TRUE/FALSE, ...)
  #   snv, qc_cols, info (min), ambiguous, duplicates, low_n (min_frac)
  # Rows with a missing or invalid SNP, BETA, SE or P are always removed. Every other
  # filter is evaluated on the valid rows, on its own (independent of the others),
  # and only the applied ones remove SNPs
  # Returns list(df = filtered data, qc = one row per step) with columns
  #   failing / failing_overlap:     SNPs failing the step on its own
  #   remaining / remaining_overlap: SNPs left after the applied steps up to this one
  applied <- function(f) isTRUE(filters[[f]][["apply"]])
  setting <- function(f, key, default) {
    x <- filters[[f]][[key]]
    if (is.null(x)) default else x
  }

  # overlap NA (not known) gives NA overlap counts
  if (!"overlap" %in% colnames(df)) df$overlap <- NA
  n_ovl <- function(x) sum(x)

  valid <- !is.na(df$snp) & !is.na(df$beta) & is.finite(df$beta) &
           !is.na(df$se) & df$se > 0 & !is.na(df$p) & df$p > 0 & df$p <= 1

  steps <- list(
    tibble::tibble(step = "input", label = "input", applied = NA,
                   failing = NA_integer_, failing_overlap = NA_integer_,
                   remaining = nrow(df), remaining_overlap = n_ovl(df$overlap)),
    tibble::tibble(step = "valid", label = "valid SNP, effect, SE and P", applied = TRUE,
                   failing = sum(!valid), failing_overlap = n_ovl(!valid & df$overlap),
                   remaining = sum(valid), remaining_overlap = n_ovl(valid & df$overlap)))
  df <- df[valid, ]

  # Pass (TRUE) / fail for each filter that applies to this GWAS
  info_min <- setting("info", "min", 0.8)
  min_frac <- setting("low_n", "min_frac", 0.5)
  pass   <- list()
  labels <- list()

  pass$snv   <- gwas_is_snv(df$a1, df$a2)
  labels$snv <- "biallelic SNV (A/C/G/T)"

  if (length(qc_cols) > 0) {
    pass$qc_cols   <- Reduce(`&`, lapply(qc_cols, \(col) !is.na(df[[col]]) & as.numeric(df[[col]]) == 0))
    labels$qc_cols <- paste0(paste(qc_cols, collapse = " & "), " == 0")
  }

  if ("info" %in% colnames(df)) {
    pass$info   <- !is.na(df$info) & df$info >= info_min
    labels$info <- paste0("INFO >= ", info_min)
  }

  pass$ambiguous   <- !gwas_is_ambiguous(df$a1, df$a2)
  labels$ambiguous <- "not strand-ambiguous (A/T, C/G)"

  dup <- df$snp[duplicated(df$snp)]
  pass$duplicates   <- !df$snp %in% dup
  labels$duplicates <- "unique SNP ID"

  if (low_n_ok) {
    n_cut <- min_frac * stats::quantile(df$n, 0.9, na.rm = TRUE, names = FALSE)
    pass$low_n   <- !is.na(df$n) & df$n >= n_cut
    labels$low_n <- paste0("N >= ", min_frac, " x 90th percentile N")
  }

  # Count each filter on its own, then apply the switched-on ones in order
  keep <- rep(TRUE, nrow(df))
  for (f in names(pass)) {
    on <- applied(f)
    if (on) keep <- keep & pass[[f]]
    steps[[length(steps) + 1]] <- tibble::tibble(
      step = f, label = labels[[f]], applied = on,
      failing = sum(!pass[[f]]), failing_overlap = n_ovl(!pass[[f]] & df$overlap),
      remaining         = if (on) sum(keep) else NA_integer_,
      remaining_overlap = if (on) n_ovl(keep & df$overlap) else NA_integer_)
  }

  df <- df[keep, ]
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
