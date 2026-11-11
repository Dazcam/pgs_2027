# test-prs_standardise.R
# - Used in prs_combine_scores.R

test_that("standardised scores have mean 0 and SD 1", {
  z <- prs_standardise(c(2, 4, 6, 8))
  expect_equal(mean(z), 0)
  expect_equal(sd(z), 1)
  expect_equal(z, (c(2, 4, 6, 8) - 5) / sd(c(2, 4, 6, 8)))
})

test_that("constant scores stop", {
  expect_error(prs_standardise(c(1, 1, 1)), "no variation")
})
