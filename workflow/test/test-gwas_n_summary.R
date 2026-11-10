# test-gwas_n_summary.R
# - Used in gwas_format.R

test_that("median and mean ignore NA", {
  expect_equal(gwas_n_summary(c(1, 2, 10, NA), "median"), 2)
  expect_equal(gwas_n_summary(c(1, 2, 9, NA), "mean"), 4)
})

test_that("all-NA or unknown summary stops", {
  expect_error(gwas_n_summary(c(NA, NA)), "No per-SNP N")
  expect_error(gwas_n_summary(1, "max"), "n_summary")
})
