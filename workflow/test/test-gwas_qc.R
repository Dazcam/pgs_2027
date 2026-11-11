# test-gwas_qc.R
# - Used in gwas_format.R

toy <- tibble::tibble(
  snp  = c("rs1", "rs2", "rs3", "rs4", "rs5", "rs6", "rs6", "rs7", "rs8", "rs9", "rs10"),
  a1   = c("A",   "A",   "A",   "A",   "A",   "C",   "C",   "C",   "G",   "AT",  "G"),
  a2   = c("G",   "G",   "G",   "T",   "C",   "T",   "T",   "T",   "A",   "A",   "G"),
  beta = c(0.1,   NA,    0.1,   0.1,   0.1,   0.1,   0.1,   0.1,   0.1,   0.1,   0.1),
  se   = c(0.01,  0.01,  0.01,  0.01,  0.01,  0.01,  0.01,  0.01,  0.01,  0.01,  0.01),
  p    = c(0.5,   0.5,   0.5,   0.5,   0.5,   0.5,   0.5,   0.5,   0.5,   0.5,   0.5),
  info = c(0.9,   0.9,   0.7,   0.9,   0.9,   0.9,   0.9,   0.9,   0.95,  0.9,   0.9),
  n    = c(100,   100,   100,   100,   100,   100,   100,   100,   10,    100,   100),
  error = c("0",  "0",   "0",   "0",   "0",   "0",   "0",   "1",   "0",   "0",   "0"))

run_qc <- function(...) {
  suppressMessages(gwas_qc(toy, ...))
}

test_that("default QC removes invalid, non-SNV, low INFO, ambiguous and duplicated SNPs", {
  res <- run_qc(info_min = 0.8)
  # rs2 NA beta; rs9 indel, rs10 G/G (not SNVs); rs3 INFO 0.7; rs4 A/T ambiguous;
  # both rs6 duplicated
  expect_equal(res$df$snp, c("rs1", "rs5", "rs7", "rs8"))
  expect_equal(res$qc$snps, c(11, 10, 8, 7, 6, 4))
  expect_equal(res$qc$step[3], "biallelic SNV (A/C/G/T)")
})

test_that("each filter can be switched off", {
  res <- run_qc(info_min = NULL, drop_ambiguous = FALSE, drop_duplicates = FALSE)
  # the SNV filter always runs (PRS-CS drops non-SNVs anyway)
  expect_equal(nrow(res$df), 8)
  expect_equal(res$qc$step, c("input", "valid effect, SE and P", "biallelic SNV (A/C/G/T)"))
})

test_that("qc_cols keep only rows equal to 0", {
  res <- run_qc(info_min = NULL, drop_ambiguous = FALSE, drop_duplicates = FALSE, qc_cols = "error")
  expect_false("rs7" %in% res$df$snp)
})

test_that("n_min_frac drops low-N SNPs", {
  res <- run_qc(info_min = NULL, n_min_frac = 0.5)
  expect_false("rs8" %in% res$df$snp)
})

test_that("invalid P and SE values are removed", {
  bad <- toy[c(1, 1, 1, 1), ]
  bad$snp <- c("rs1", "rs2", "rs3", "rs4")
  bad$p   <- c(0.5, 0,    1.2,  0.5)
  bad$se  <- c(0.01, 0.01, 0.01, -1)
  res <- gwas_qc(bad, info_min = NULL, drop_ambiguous = FALSE, drop_duplicates = FALSE)
  expect_equal(res$df$snp, "rs1")
})

test_that("QC that leaves no SNPs stops", {
  bad <- toy[1:2, ]
  bad$p <- c(0, 2)
  expect_error(gwas_qc(bad), "No SNPs left")
})
