# test-gwas_in_overlap.R
# - Used in gwas_format.R

test_that("allele key ignores allele order and strand", {
  expect_equal(gwas_allele_key(c("A", "G", "T", "C"), c("G", "A", "C", "T")), rep("AG", 4))
  expect_equal(gwas_allele_key("A", "T"), gwas_allele_key("T", "A"))
  expect_equal(gwas_allele_key("a", "g"), "AG")
})

test_that("allele key is NA for non-SNV alleles", {
  expect_true(all(is.na(gwas_allele_key(c("AT", "G", "N", "-"), c("A", "G", "A", "C")))))
})

ref <- tibble::tibble(chr = c("1",   "1",   "2",   "2"),
                      snp = c("rs1", "rs2", "rs3", "rs4"),
                      a1  = c("A",   "A",   "G",   "A"),
                      a2  = c("G",   "C",   "T",   "C"))
bim <- tibble::tibble(chr = c("1",   "1",   "2",   "1"),
                      snp = c("rs1", "rs2", "rs3", "rs4"),
                      a1  = c("G",   "A",   "A",   "A"),
                      a2  = c("A",   "G",   "C",   "C"))

test_that("SNPs need the rsID with compatible alleles in the reference and .bim", {
  # rs1: .bim swapped, GWAS on the other strand -> used
  # rs2: .bim alleles (A/G) incompatible with reference (A/C) -> not used
  # rs3: .bim on the other strand to the reference -> used
  # rs4: on different chromosomes in reference and .bim -> not used
  # rs5: not in reference -> not used
  res <- gwas_in_overlap(snp = c("rs1", "rs2", "rs3", "rs4", "rs5"),
                         a1  = c("T",   "A",   "G",   "A",   "A"),
                         a2  = c("C",   "C",   "T",   "C",   "G"),
                         ref = ref, bim = bim)
  expect_equal(res, c(TRUE, FALSE, TRUE, FALSE, FALSE))
})

test_that("GWAS alleles must match the SNP, not just the rsID", {
  res <- gwas_in_overlap(snp = c("rs1", "rs1"), a1 = c("A", "A"), a2 = c("G", "C"),
                         ref = ref, bim = bim)
  expect_equal(res, c(TRUE, FALSE))
})

test_that("chr prefixes in the .bim are ignored", {
  bim_chr <- dplyr::mutate(bim, chr = paste0("chr", chr))
  expect_true(gwas_in_overlap("rs1", "A", "G", ref, bim_chr))
})
