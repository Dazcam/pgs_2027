# test-gwas_read.R
# - Used in gwas_format.R

test_that("tab and space files read only the requested columns, as text", {
  tab <- tempfile(fileext = ".tsv.gz")
  spc <- tempfile(fileext = ".txt")
  con <- gzfile(tab, "w")
  writeLines(c("#CHROM\tID\tEA\tNEA\tBETA", "1\trs1\tA\tG\t0.10"), con); close(con)
  writeLines(c("CHR SNP A1 A2 OR", "8 rs2 C T 0.99631"), spc)

  t1 <- gwas_read(tab, "tab", c("ID", "EA", "BETA"))
  t2 <- gwas_read(spc, "space", c("SNP", "OR"))
  expect_equal(colnames(t1), c("ID", "EA", "BETA"))
  expect_equal(t1$BETA, "0.10")
  expect_equal(t2$OR, "0.99631")
})

test_that("unknown delimiter stops", {
  expect_error(gwas_read("x", "comma", "a"), "delim")
})
