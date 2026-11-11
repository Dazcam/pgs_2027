# test-gwas_is_snv.R
# - Used in gwas_format.R

test_that("single-base A/C/G/T pairs with different alleles are SNVs, in any case", {
  expect_equal(gwas_is_snv(c("A", "c", "G", "T"), c("G", "t", "A", "A")), rep(TRUE, 4))
})

test_that("indels, other codes, same alleles and missing values are not SNVs", {
  expect_equal(gwas_is_snv(c("AT", "I", "A", "A", "N", NA), c("A", "D", "A", "<DEL>", "G", "C")),
               rep(FALSE, 6))
})
