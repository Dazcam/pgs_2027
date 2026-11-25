# test-share_functions.R
# - Used in share_scores.R, share_pseudobulk.R, share_readme.R

labels <- list("Glu-UL" = "Glu-A", "Glu-DL" = "Glu-B")

test_that("relabel renames L1 types and their subclusters", {
  expect_equal(share_relabel(c("Glu-UL", "Glu-DL-2", "GABA-0", "NPC"), labels),
               c("Glu-A", "Glu-B-2", "GABA-0", "NPC"))
})

test_that("relabel needs named labels", {
  expect_error(share_relabel("Glu-UL", c("Glu-A")), "named")
})

test_that("expected donors keep order and drop only allowed absences", {
  donors <- c("a", "b", "c", "d")
  expect_equal(share_expected_donors(donors, c("a", "c", "d"), "b"), c("a", "c", "d"))
  expect_equal(share_expected_donors(donors, donors, "b"), donors)          # allowed but present
  expect_equal(share_expected_donors(donors, c("a", "d"), "b"), c("a", "c", "d"))  # c not allowed
  expect_equal(share_expected_donors(donors, donors), donors)
})

test_that("match_donors returns positions in expected order, allowing extra donors", {
  expect_equal(share_match_donors(c("13008_A", "1102", "12971_FC", "99_FC"),
                                  c("1102", "12971_FC", "13008_A"), "file"),
               c(2, 3, 1))
})

test_that("match_donors needs exact IDs", {
  expect_error(share_match_donors(c("1102_FC", "12971_FC"), c("1102", "12971_FC"), "file"),
               "1 expected donors missing: 1102")
})

test_that("match_donors stops on duplicated IDs", {
  expect_error(share_match_donors(c("1102", "1102"), "1102", "file"), "duplicated")
})

test_that("check_equal passes equal values and tiny numeric differences", {
  expect_true(share_check_equal(c(0.1, NA, -0.003), c(0.1 + 1e-15, NA, -0.003), 1:3, "PC1"))
  expect_true(share_check_equal(c("M", "F"), c("M", "F"), 1:2, "sex"))
})

test_that("check_equal stops and names donors that differ", {
  expect_error(share_check_equal(c(1, 2, 3), c(1, 2.1, 3), c("a", "b", "c"), "PCW"),
               "1 donors: b")
  expect_error(share_check_equal(c("M", NA), c("M", "F"), c("a", "b"), "sex"), "b")
})

test_that("cell type lines put L1 first, then two parent types per line", {
  ct <- c("Glu-A", "Glu-B", "NPC", "Glu-A-0", "Glu-A-1", "Glu-B-0", "NPC-0")
  expect_equal(share_cell_type_lines(ct),
               c("  Glu-A, Glu-B, NPC,", "  Glu-A-0, Glu-A-1, Glu-B-0,", "  NPC-0"))
})
