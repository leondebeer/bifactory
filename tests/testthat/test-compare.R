# esem_compare(): the comparison CFA must be fitted with the ESEM's own
# estimator settings, and the table must show the same (robust) statistics the
# rest of the package reports.

skip_on_cran()

hs    <- lavaan::HolzingerSwineford1939[, paste0("x", 1:9)]
quiet <- function(expr) suppressWarnings(suppressMessages(expr))
fit   <- quiet(esem(hs, nfactors = 3))
cfa_model <- "
  Visual  =~ x1 + x2 + x3
  Textual =~ x4 + x5 + x6
  Speed   =~ x7 + x8 + x9
"

test_that("esem_compare() tabulates the robust MLR chi-square for both models", {
  cmp <- quiet(esem_compare(fit, cfa_model, data = hs))
  tb  <- cmp$fit_table
  x2_esem <- unname(lavaan::fitMeasures(fit$lavaan_fit, "chisq.scaled"))
  expect_equal(tb$ESEM[tb$index == "chisq"], x2_esem, tolerance = 1e-4)
  # the CFA inherits the ESEM's robust test, so it has a scaled chi-square too
  expect_true(any(grepl("yuan.bentler", lavaan::lavInspect(cmp$cfa_fit, "options")$test)))
  x2_cfa <- unname(lavaan::fitMeasures(cmp$cfa_fit, "chisq.scaled"))
  expect_equal(tb$CFA[tb$index == "chisq"], x2_cfa, tolerance = 1e-4)
  expect_equal(tb$CFA[tb$index == "cfi"],
               unname(lavaan::fitMeasures(cmp$cfa_fit, "cfi.scaled")), tolerance = 1e-4)
})
