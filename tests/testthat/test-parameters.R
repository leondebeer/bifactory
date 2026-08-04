# parameters() classed return + print.bifactory_parameters ------------------
# One small ESEM fit (HolzingerSwineford1939, 9 items, 3 factors) shared
# across tests; fitting is the slow part, so do it once.

hs <- lavaan::HolzingerSwineford1939[, paste0("x", 1:9)]
fit <- suppressWarnings(suppressMessages(esem(hs, nfactors = 3)))

test_that("parameters() returns a visible classed data frame", {
  expect_visible(parameters(fit))
  p <- parameters(fit)
  expect_s3_class(p, "bifactory_parameters")
  expect_s3_class(p, "data.frame")
  expect_true(all(c("factor", "item", "std", "se", "z", "p", "loading_type")
                  %in% names(p)))
  expect_gt(nrow(p), 0L)
})

test_that("parameters() carries the attributes print() renders from", {
  p <- parameters(fit, digits = 2, highlight_primary = FALSE)
  expect_type(attr(p, "blocks"), "list")
  expect_identical(attr(p, "digits"), 2)
  expect_false(attr(p, "highlight_primary"))
})

test_that("print.bifactory_parameters renders the table and returns invisibly", {
  p <- parameters(fit)
  expect_output(print(p), "Standardized Parameters (STDYX)", fixed = TRUE)
  expect_invisible(print(p))
  expect_identical(withVisible(print(p))$value, p)
})

test_that("row subsetting keeps the class and render attributes", {
  p  <- parameters(fit)
  ps <- p[seq_len(min(3L, nrow(p))), ]
  expect_s3_class(ps, "bifactory_parameters")
  expect_type(attr(ps, "blocks"), "list")
})

test_that("print falls back to print.data.frame when attributes are lost", {
  ps <- parameters(fit)[seq_len(3L), ]
  attr(ps, "blocks") <- NULL   # simulate loss (e.g. rbind/merge strip attrs)
  expect_output(print(ps), "std")   # fallback renders without error
})
