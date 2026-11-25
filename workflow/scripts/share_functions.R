#--------------------------------------------------------------------------------------
#
#    Data sharing functions
#
#--------------------------------------------------------------------------------------
#
# Functions: share_relabel()       pipeline cell type names -> manuscript names
#            share_expected_donors()  donors expected in a file, allowing listed absences
#            share_match_donors()  position of each expected donor in a file's IDs
#            share_check_equal()   stop if two value vectors differ, listing donors
#            share_cell_type_lines()  cell types laid out for the README
#
# Packages:  base R only
#
# Used in:   share_scores.R, share_pseudobulk.R, share_readme.R
#
#--------------------------------------------------------------------------------------

share_relabel <- function(x, labels) {

  # Apply each substitution in labels (pipeline text -> shared text), in order, to every
  # occurrence in each string: "Glu-UL" -> "Glu-A" and "Glu-UL-0" -> "Glu-A-0"
  labels <- unlist(labels)
  if (is.null(names(labels)) || any(names(labels) == "")) {
    stop("share_relabel: labels must be a named list or vector")
  }
  x <- as.character(x)
  for (old in names(labels)) {
    x <- gsub(old, labels[[old]], x, fixed = TRUE)
  }
  x
}

share_expected_donors <- function(donors, present, allowed_missing = character(0)) {

  # Donors expected in a file, in the order of donors: every donor, except those in
  # allowed_missing that are absent from present
  donors[donors %in% present | !donors %in% allowed_missing]
}

share_match_donors <- function(ids, expected, what) {

  # Position in ids of each expected donor, matched exactly. Stops if any expected donor
  # is missing or if ids has duplicates; donors in ids but not expected are allowed
  ids <- as.character(ids)
  dup <- unique(ids[duplicated(ids)])
  if (length(dup) > 0) {
    stop(what, ": duplicated donor IDs: ", paste(dup, collapse = ", "))
  }
  missing <- setdiff(expected, ids)
  if (length(missing) > 0) {
    stop(what, ": ", length(missing), " expected donors missing: ",
         paste(utils::head(missing, 20), collapse = ", "))
  }
  match(expected, ids)
}

share_check_equal <- function(x, y, ids, what, tolerance = 1e-8) {

  # Stop if x and y differ for any donor (numeric: relative tolerance; otherwise exact,
  # with NA only equal to NA)
  if (is.numeric(x) && is.numeric(y)) {
    same <- (is.na(x) & is.na(y)) |
            (!is.na(x) & !is.na(y) & abs(x - y) <= tolerance * pmax(1, abs(y)))
  } else {
    x <- as.character(x)
    y <- as.character(y)
    same <- (is.na(x) & is.na(y)) | (!is.na(x) & !is.na(y) & x == y)
  }
  if (!all(same)) {
    stop(what, ": values differ from the source for ", sum(!same), " donors: ",
         paste(utils::head(ids[!same], 20), collapse = ", "))
  }
  invisible(TRUE)
}

share_cell_type_lines <- function(cell_types, per_line = 2, indent = "  ") {

  # L1 types (no "-<n>" suffix) on the first line, then the L2 subclusters with
  # per_line parent types on each line, in config order
  is_l2  <- grepl("-[0-9]+$", cell_types)
  l1     <- cell_types[!is_l2]
  l2     <- cell_types[is_l2]
  parent <- sub("-[0-9]+$", "", l2)
  groups <- split(l2, factor(parent, levels = unique(parent)))
  chunks <- split(groups, ceiling(seq_along(groups) / per_line))
  lines  <- c(paste(l1, collapse = ", "),
              vapply(chunks, \(g) paste(unlist(g), collapse = ", "), character(1)))
  lines  <- unname(lines[lines != ""])
  paste0(indent, lines, c(rep(",", length(lines) - 1), ""))
}

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
