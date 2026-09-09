# Roxygen blocks are written in markdown; without `Roxygen: list(markdown =
# TRUE)` in DESCRIPTION the markers pass into the Rd files verbatim, which is
# how the 0.5.2 CRAN manual came to show literal "**bold**" and "|---|" tables.
# Reads the source man/ directory when running from the package tree, the
# installed Rd database under R CMD check.

test_that("no raw markdown survives in the Rd files", {
  root <- test_path("..", "..")
  db <- if (dir.exists(file.path(root, "man"))) {
    tools::Rd_db(dir = root)
  } else {
    tools::Rd_db("bifactory")
  }
  expect_gt(length(db), 0L)
  raw <- vapply(db, function(rd) paste(as.character(rd), collapse = ""), "")
  bad <- names(raw)[grepl("\\*\\*|\n\\| |\\|-{2,}\\|", raw)]
  expect_identical(bad, character(0))
})
