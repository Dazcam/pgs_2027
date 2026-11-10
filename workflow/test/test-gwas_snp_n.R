# test-gwas_snp_n.R
# - Used in gwas_format.R

df <- tibble::tibble(se = c(0.02, 0.02, 0.01), freq = c(0.25, 0.05, 0.6),
                     NEFF = c("1000", "2000", "3000"), Nca = c("100", "200", "300"),
                     Nco = c("300", "200", "900"))

test_that("column method multiplies by scale", {
  expect_equal(gwas_snp_n(df, list(method = "column", col = "NEFF", scale = 2), "binary", c(0.1, 0.4)),
               c(2000, 4000, 6000))
})

test_that("case_control gives effective N 4/(1/cases + 1/controls)", {
  # 100 cases, 300 controls -> 300; balanced 200/200 -> 400
  expect_equal(gwas_snp_n(df, list(method = "case_control", cases = "Nca", controls = "Nco"),
                          "binary", c(0.1, 0.4)),
               c(300, 400, 900))
})

test_that("fixed repeats the value", {
  expect_equal(gwas_snp_n(df, list(method = "fixed", value = 44367), "binary", c(0.1, 0.4)),
               rep(44367, 3))
})

test_that("auto uses k/(2f(1-f)SE^2) within the MAF window, NA outside", {
  # f = 0.25: 2f(1-f) = 0.375; SE^2 = 4e-4 -> continuous 6666.67, binary 26666.67
  # f = 0.05 (MAF 0.05) outside window; f = 0.6 (MAF 0.4) inside: 2f(1-f) = 0.48, SE^2 = 1e-4
  cont <- gwas_snp_n(df, list(method = "auto"), "continuous", c(0.1, 0.4))
  bin  <- gwas_snp_n(df, list(method = "auto"), "binary", c(0.1, 0.4))
  expect_equal(cont, c(1 / (0.375 * 4e-4), NA, 1 / (0.48 * 1e-4)))
  expect_equal(bin, 4 * cont)
})

test_that("auto without a frequency column stops", {
  expect_error(gwas_snp_n(df[, c("se")], list(method = "auto"), "binary", c(0.1, 0.4)), "frequency")
})
