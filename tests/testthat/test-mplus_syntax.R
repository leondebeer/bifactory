# generate_mplus_syntax() and generate_mplus_besem_syntax() must write input
# that Mplus accepts: every VARIABLE statement ends in a semicolon, no line
# exceeds Mplus's 90-character limit, and ordinal items get a CATEGORICAL
# statement with WLSMV.

items19 <- paste0("item_name", sprintf("%02d", 1:19))
sf19    <- list(SpecA = items19[1:10], SpecB = items19[11:19])

statement <- function(syntax, keyword) {
  # the full statement from `keyword` to its terminating semicolon
  regmatches(syntax, regexpr(paste0(keyword, "[^;]*;"), syntax))
}

test_that("B-ESEM syntax: terminated statements and lines within 90 characters", {
  s <- generate_mplus_besem_syntax(sf19, items19)
  lines <- strsplit(s, "\n")[[1]]
  expect_true(all(nchar(lines) <= 90))
  # NAMES must stop before USEVARIABLES, USEVARIABLES before MISSING
  expect_false(grepl("USEVARIABLES", statement(s, "NAMES =")))
  expect_false(grepl("MISSING", statement(s, "USEVARIABLES =")))
  expect_match(s, "ESTIMATOR = MLR;", fixed = TRUE)
  expect_false(grepl("CATEGORICAL", s))
})

test_that("B-ESEM syntax for ordinal items declares them CATEGORICAL under WLSMV", {
  s <- generate_mplus_besem_syntax(sf19, items19, ordered = TRUE)
  expect_true(all(nchar(strsplit(s, "\n")[[1]]) <= 90))
  cat_stmt <- statement(s, "CATEGORICAL =")
  expect_true(all(vapply(items19, grepl, logical(1), x = cat_stmt, fixed = TRUE)))
  expect_match(s, "ESTIMATOR = WLSMV;", fixed = TRUE)
  expect_match(s, "ROTATION = TARGET (orthogonal);", fixed = TRUE)
  # besem_ordered() fits correlations with residual variances 1 - communality,
  # which is Mplus's default delta parameterization: under theta the
  # residual variances are fixed and Mplus cannot report a negative one.
  expect_false(grepl("PARAMETERIZATION", s))
})

test_that("ESEM syntax: terminated statements, wrapped lines, ordinal option", {
  s <- generate_mplus_syntax(sf19)
  expect_true(all(nchar(strsplit(s, "\n")[[1]]) <= 90))
  expect_false(grepl("USEVARIABLES", statement(s, "NAMES =")))
  expect_false(grepl("MISSING", statement(s, "USEVARIABLES =")))

  # esem_ordered()'s default lavaan path uses the theta parameterization
  s_ord <- generate_mplus_syntax(sf19, ordered = TRUE)
  expect_match(s_ord, "CATEGORICAL =", fixed = TRUE)
  expect_match(s_ord, "ESTIMATOR = WLSMV;", fixed = TRUE)
  expect_match(s_ord, "PARAMETERIZATION = THETA;", fixed = TRUE)
})

test_that("generated input tightens Mplus's convergence criterion", {
  # Mplus's default criterion stops early enough to move chi-square and
  # loadings away from the R solution; same settings as the invariance syntax
  for (s in list(generate_mplus_besem_syntax(sf19, items19, ordered = TRUE),
                 generate_mplus_syntax(sf19))) {
    expect_match(s, "CONVERGENCE = 0.000001;", fixed = TRUE)
    expect_match(s, "ITERATIONS = 10000;", fixed = TRUE)
  }
})

test_that("missing_code = '.' writes MISSING = .;", {
  expect_match(generate_mplus_besem_syntax(sf19, items19, missing_code = "."),
               "MISSING = .;", fixed = TRUE)
  expect_match(generate_mplus_syntax(sf19, missing_code = "."),
               "MISSING = .;", fixed = TRUE)
  expect_match(generate_mplus_syntax(sf19), "MISSING ARE ALL (999);", fixed = TRUE)
})
