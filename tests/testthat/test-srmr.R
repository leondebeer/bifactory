# Single-group SRMR by the Mplus definition (Asparouhov & Muthen 2018), checked
# against Mplus 9 on psych::bfi (validation/_validate_ewc_mplus.R): ESEM 0.026
# and B-ESEM 0.018 under MLR (listwise, N = 2436); ESEM 0.026 under WLSMV.
# Mplus always models the means, so the continuous denominator is
# p(p+1)/2 + p even when lavaan fits no mean structure.

skip_on_cran()

quiet <- function(expr) suppressWarnings(suppressMessages(expr))
fl <- list(Agree = paste0("A", 1:5), Consc = paste0("C", 1:5), Extra = paste0("E", 1:5),
           Neuro = paste0("N", 1:5), Open = paste0("O", 1:5))
items <- unname(unlist(fl))
bfi   <- psych::bfi[, items]
bfi_c <- bfi[stats::complete.cases(bfi), ]

test_that("continuous ESEM and B-ESEM SRMR match Mplus at 3 dp", {
  spec <- do.call(specify_model, c(fl, list(data = bfi_c, ordered = FALSE, label = "bfi")))
  res  <- quiet(run_comparison(spec, run_alignment = FALSE))
  expect_identical(fit_indices(res$fit_esem)[["SRMR"]],  "0.026")
  expect_identical(fit_indices(res$fit_besem)[["SRMR"]], "0.018")
})

test_that("ordered ESEM (lavaan WLSMV path) SRMR matches Mplus at 3 dp", {
  spec <- do.call(specify_model, c(fl, list(data = bfi, ordered = TRUE, label = "bfi ord")))
  fit  <- quiet(esem_ordered(data = bfi, nfactors = 5L, indicators = items, rotation = "target",
                             target = spec$target, factor_names = spec$factor_names,
                             missing = "pairwise"))
  expect_identical(fit_indices(fit)[["SRMR"]], "0.026")
})
