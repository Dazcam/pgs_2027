# test-gwas_is_ambiguous.R
# - Used in gwas_format.R

test_that("A/T and C/G pairs are ambiguous in either order and case", {
  expect_equal(gwas_is_ambiguous(c("A", "T", "c", "G"), c("T", "A", "g", "C")), rep(TRUE, 4))
})

test_that("other pairs and indels are not ambiguous", {
  expect_equal(gwas_is_ambiguous(c("A", "A", "C", "AT"), c("G", "C", "T", "A")), rep(FALSE, 4))
})
