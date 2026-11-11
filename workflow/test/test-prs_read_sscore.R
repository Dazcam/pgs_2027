# test-prs_read_sscore.R
# - Used in prs_combine_scores.R

test_that("reads PLINK 2 sscore with #IID or #FID IID headers", {
  f1 <- tempfile(fileext = ".sscore")
  f2 <- tempfile(fileext = ".sscore")
  writeLines(c("#IID\tALLELE_CT\tNAMED_ALLELE_DOSAGE_SUM\tSCORE1_AVG\tSCORE1_SUM",
               "1102\t200\t95.5\t0.001\t0.2",
               "1103\t200\t101.2\t-0.002\t-0.4"), f1)
  writeLines(c("#FID\tIID\tALLELE_CT\tNAMED_ALLELE_DOSAGE_SUM\tSCORE1_AVG\tSCORE1_SUM",
               "1102\t1102\t200\t95.5\t0.001\t0.2"), f2)

  s1 <- prs_read_sscore(f1)
  expect_equal(s1$iid, c("1102", "1103"))
  expect_equal(s1$allele_ct, c(200, 200))
  expect_equal(s1$score_sum, c(0.2, -0.4))
  expect_equal(prs_read_sscore(f2)$iid, "1102")
})

test_that("missing SCORE1_SUM gives a clear error", {
  f <- tempfile(fileext = ".sscore")
  writeLines(c("#IID\tALLELE_CT\tSCORE1_AVG", "1102\t200\t0.001"), f)
  expect_error(prs_read_sscore(f), "cols=\\+scoresums")
})
