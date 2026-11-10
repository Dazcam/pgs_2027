# test-gwas_to_beta.R
# - Used in gwas_format.R

test_that("BETA passes through and OR becomes log(OR)", {
  expect_equal(gwas_to_beta(c(0.1, -0.2), "beta"), c(0.1, -0.2))
  expect_equal(gwas_to_beta(c(1, exp(0.5), exp(-0.3)), "or"), c(0, 0.5, -0.3))
})

test_that("unknown effect type stops", {
  expect_error(gwas_to_beta(1, "z"), "effect_type")
})
