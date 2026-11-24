# test-gwas_qc.R
# - Used in gwas_format.R

toy <- tibble::tibble(
  snp     = c("rs1", "rs2", "rs3", "rs4", "rs5", "rs6", "rs6", "rs7", "rs8", "rs9", "rs10"),
  a1      = c("A",   "A",   "A",   "A",   "A",   "C",   "C",   "C",   "G",   "AT",  "G"),
  a2      = c("G",   "G",   "G",   "T",   "C",   "T",   "T",   "T",   "A",   "A",   "G"),
  beta    = c(0.1,   NA,    0.1,   0.1,   0.1,   0.1,   0.1,   0.1,   0.1,   0.1,   0.1),
  se      = c(0.01,  0.01,  0.01,  0.01,  0.01,  0.01,  0.01,  0.01,  0.01,  0.01,  0.01),
  p       = c(0.5,   0.5,   0.5,   0.5,   0.5,   0.5,   0.5,   0.5,   0.5,   0.5,   0.5),
  info    = c(0.9,   0.9,   0.7,   0.9,   0.9,   0.9,   0.9,   0.9,   0.95,  0.9,   0.9),
  n       = c(100,   100,   100,   100,   100,   100,   100,   100,   10,    100,   100),
  error   = c("0",   "0",   "0",   "0",   "0",   "0",   "0",   "1",   "0",   "0",   "0"),
  overlap = c(TRUE,  TRUE,  TRUE,  FALSE, TRUE,  TRUE,  TRUE,  FALSE, TRUE,  FALSE, FALSE))

all_on <- list(snv = list(apply = TRUE), qc_cols = list(apply = TRUE),
               info = list(apply = TRUE, min = 0.8), ambiguous = list(apply = TRUE),
               duplicates = list(apply = TRUE), low_n = list(apply = TRUE, min_frac = 0.5))
all_off <- lapply(all_on, \(f) modifyList(f, list(apply = FALSE)))

run_qc <- function(filters, ...) {
  suppressMessages(gwas_qc(toy, filters = filters, qc_cols = "error", ...))
}
step_row <- function(res, s) res$qc[res$qc$step == s, ]

test_that("all filters on removes invalid, non-SNV, flagged, low INFO, ambiguous, duplicated and low-N SNPs", {
  res <- run_qc(all_on)
  # rs2 NA beta; rs9 indel, rs10 G/G; rs7 error flag; rs3 INFO 0.7; rs4 A/T;
  # both rs6 duplicated; rs8 low N
  expect_equal(res$df$snp, c("rs1", "rs5"))
  expect_equal(res$qc$step, c("input", "valid", "snv", "qc_cols", "info", "ambiguous", "duplicates", "low_n"))
  expect_equal(res$qc$remaining, c(11, 10, 8, 7, 6, 5, 3, 2))
})

test_that("filters switched off are counted but remove nothing", {
  res <- run_qc(all_off)
  # only the invalid row is removed
  expect_equal(nrow(res$df), 10)
  expect_equal(res$qc$applied[-1], c(TRUE, rep(FALSE, 6)))
  expect_true(all(is.na(res$qc$remaining[-(1:2)])))
  expect_equal(step_row(res, "snv")$failing, 2)
  expect_equal(step_row(res, "info")$failing, 1)
  expect_equal(step_row(res, "ambiguous")$failing, 1)
  expect_equal(step_row(res, "duplicates")$failing, 2)
  expect_equal(step_row(res, "low_n")$failing, 1)
})

test_that("failing counts each filter on its own, whatever else is applied", {
  on  <- run_qc(all_on)
  off <- run_qc(all_off)
  expect_equal(on$qc$failing, off$qc$failing)
})

test_that("overlap counts are restricted to SNPs PRS-CS would use", {
  res <- run_qc(all_off)
  expect_equal(step_row(res, "input")$remaining_overlap, 7)
  expect_equal(step_row(res, "valid")$failing_overlap, 1)    # rs2
  expect_equal(step_row(res, "snv")$failing_overlap, 0)      # rs9, rs10 not in overlap
  expect_equal(step_row(res, "info")$failing_overlap, 1)     # rs3
  expect_equal(step_row(res, "duplicates")$failing_overlap, 2)
  expect_equal(step_row(run_qc(all_on), "low_n")$remaining_overlap, 2)  # rs1, rs5
})

test_that("only the switched-on filters remove SNPs", {
  dup_only <- all_off
  dup_only$duplicates$apply <- TRUE
  res <- run_qc(dup_only)
  expect_false("rs6" %in% res$df$snp)
  expect_true(all(c("rs3", "rs4", "rs7", "rs9") %in% res$df$snp))
})

test_that("info and qc_cols rows only appear where they apply", {
  no_info <- dplyr::select(toy, -info)
  res <- suppressMessages(gwas_qc(no_info, filters = all_off))
  expect_false(any(c("info", "qc_cols") %in% res$qc$step))
})

test_that("low-N filter is not evaluated when low_n_ok is FALSE", {
  res <- run_qc(all_on, low_n_ok = FALSE)
  expect_false("low_n" %in% res$qc$step)
  expect_true("rs8" %in% res$df$snp)
})

test_that("overlap counts are NA when overlap is not known", {
  res <- suppressMessages(gwas_qc(dplyr::select(toy, -overlap), filters = all_off))
  expect_true(all(is.na(res$qc$remaining_overlap[1:2])))
})

test_that("invalid P and SE values are always removed", {
  bad <- toy[c(1, 1, 1, 1), ]
  bad$snp <- c("rs1", "rs2", "rs3", "rs4")
  bad$p   <- c(0.5, 0,    1.2,  0.5)
  bad$se  <- c(0.01, 0.01, 0.01, -1)
  res <- gwas_qc(bad, filters = all_off)
  expect_equal(res$df$snp, "rs1")
})

test_that("QC that leaves no SNPs stops", {
  bad <- toy[1:2, ]
  bad$p <- c(0, 2)
  expect_error(gwas_qc(bad, filters = all_off), "No SNPs left")
})
