# ESEM-within-CFA / B-ESEM-within-CFA: exactness of the re-expression ---------
# HolzingerSwineford1939 (9 items, 3 factors) plus a quartile-cut copy for the
# ordered path. The CFA/ESEM/B-ESEM trio is slow, so each estimator's trio is
# fitted once and shared across tests.

skip_on_cran()

hs <- lavaan::HolzingerSwineford1939
fl <- list(Vis = c("x1", "x2", "x3"), Tex = c("x4", "x5", "x6"), Spd = c("x7", "x8", "x9"))
hs_ord <- hs
for (it in paste0("x", 1:9))
  hs_ord[[it]] <- cut(hs[[it]], quantile(hs[[it]], c(0, .25, .5, .75, 1)),
                      include.lowest = TRUE, labels = FALSE)

quiet  <- function(expr) suppressWarnings(suppressMessages(expr))
spec_c <- do.call(specify_model, c(fl, list(data = hs,     ordered = FALSE, label = "HS")))
spec_o <- do.call(specify_model, c(fl, list(data = hs_ord, ordered = TRUE,  label = "HS ord")))
res_c  <- quiet(run_comparison(spec_c, mplus_folder = NULL, run_alignment = FALSE))
res_o  <- quiet(run_comparison(spec_o, mplus_folder = NULL, run_alignment = FALSE))
k      <- length(fl) + 1L            # G + specific factors in the B-ESEM

fm <- function(fit, ordered) {
  keys <- if (ordered) c("chisq.scaled", "df.scaled") else c("chisq", "df")
  unname(lavaan::fitMeasures(fit, keys))
}
fvar <- function(fit, factors) {
  pe <- lavaan::parameterEstimates(fit)
  pe[pe$op == "~~" & pe$lhs == pe$rhs & pe$lhs %in% factors, c("lhs", "est", "se")]
}
fcov <- function(fit, factors) {
  pe <- lavaan::parameterEstimates(fit)
  pe[pe$op == "~~" & pe$lhs != pe$rhs & pe$lhs %in% factors & pe$rhs %in% factors, ]
}
ref_loading_free <- function(fit, factor, item) {
  pt <- lavaan::parameterTable(fit)
  pt$free[pt$op == "=~" & pt$lhs == factor & pt$rhs == item] != 0L
}

# -- Plain EWC: both identification modes must return the ESEM's chi2 and df --

test_that("fixed-variance EWC re-expresses the ESEM exactly (guard)", {
  for (est in list(list(res_c, spec_c, FALSE), list(res_o, spec_o, TRUE))) {
    ewc <- quiet(fit_ewc(est[[1]]$fit_esem, est[[2]], var_fixed = TRUE))
    expect_equal(fm(ewc$lavaan_fit, est[[3]]), fm(est[[1]]$fit_esem$lavaan_fit, est[[3]]),
                 tolerance = 1e-4)
  }
})

test_that("free-variance EWC re-expresses a continuous ESEM exactly", {
  ewc <- quiet(fit_ewc(res_c$fit_esem, spec_c, var_fixed = FALSE))
  expect_equal(fm(ewc$lavaan_fit, FALSE), fm(res_c$fit_esem$lavaan_fit, FALSE), tolerance = 1e-4)
  v <- fvar(ewc$lavaan_fit, names(fl))
  expect_true(all(v$se > 0))                          # variances estimated, not fixed
  expect_equal(v$est, rep(1, 3), tolerance = 0.01)    # same metric as the std.lv ESEM
  for (f in names(fl)) expect_false(ref_loading_free(ewc$lavaan_fit, f, ewc$referents[[f]]))
})

test_that("free-variance EWC re-expresses an ordered ESEM exactly, with unit variances", {
  ewc <- quiet(fit_ewc(res_o$fit_esem, spec_o, var_fixed = FALSE))
  expect_equal(fm(ewc$lavaan_fit, TRUE), fm(res_o$fit_esem$lavaan_fit, TRUE), tolerance = 1e-4)
  v <- fvar(ewc$lavaan_fit, names(fl))
  expect_true(all(v$se > 0))
  expect_equal(v$est, rep(1, 3), tolerance = 0.01)
})

# -- B-EWC: df is the B-ESEM df + k(k-1)/2 by construction (guard) -------------

test_that("B-EWC has k(k-1)/2 more df than the B-ESEM and the same chi-square under ML", {
  fb <- res_c$fit_besem
  b  <- quiet(fit_ewc(fb, spec_c))
  expect_equal(fm(b$lavaan_fit, FALSE)[2], fm(fb$lavaan_fit, FALSE)[2] + k * (k - 1) / 2)
  expect_equal(fm(b$lavaan_fit, FALSE)[1], fm(fb$lavaan_fit, FALSE)[1], tolerance = 1e-3)
  cv <- fcov(b$lavaan_fit, c(fb$g_name, names(fl)))
  expect_equal(nrow(cv), k * (k - 1) / 2)
  expect_true(all(cv$est == 0 & cv$se == 0))
})

# -- B-EWC: var_fixed = FALSE must be honoured, not silently ignored ----------

test_that("B-EWC honours var_fixed = FALSE", {
  fb   <- res_c$fit_besem
  facs <- c(fb$g_name, names(fl))
  b1   <- quiet(fit_ewc(fb, spec_c))
  b2   <- quiet(fit_ewc(fb, spec_c, var_fixed = FALSE))
  expect_false(b2$var_fixed)
  expect_false(identical(b1$syntax, b2$syntax))
  v <- fvar(b2$lavaan_fit, facs)
  expect_true(all(v$se > 0))
  expect_equal(v$est, rep(1, k), tolerance = 0.01)
  cv <- fcov(b2$lavaan_fit, facs)
  expect_true(all(cv$est == 0 & cv$se == 0))          # orthogonality kept
  expect_equal(fm(b2$lavaan_fit, FALSE), fm(b1$lavaan_fit, FALSE), tolerance = 1e-4)  # same model
  for (f in facs) expect_false(ref_loading_free(b2$lavaan_fit, f, b2$referents[[f]]))
})

# -- Ordered B-EWC: report the B-EWC's own fit, not the source B-ESEM's -------

test_that("ordered B-EWC reports its own df rather than the source B-ESEM's", {
  fb     <- res_o$fit_besem
  b      <- quiet(fit_ewc(fb, spec_o))
  own_df <- fb$wlsmv_stats$df + k * (k - 1) / 2
  expect_null(b$wlsmv_stats)
  expect_equal(fm(b$lavaan_fit, TRUE)[2], own_df)
  expect_output(quiet(print(b)), sprintf("X2\\(%d\\)", own_df))
})

# -- WLSMV SRMR: EWC must report the Mplus-denominator SRMR, like the pipeline

test_that("ordered EWC reports the same SRMR as the ESEM it re-expresses", {
  ewc <- quiet(fit_ewc(res_o$fit_esem, spec_o))
  ct  <- compare_ewc(res_o, ewc)
  srmr_esem <- as.numeric(ct$ESEM_R[ct$Index == "SRMR"])
  expect_equal(ct$EWC_varfix[ct$Index == "SRMR"], srmr_esem, tolerance = 1e-3)
  expect_output(quiet(print(ewc)), sprintf("SRMR = %.3f", srmr_esem))
})

# -- MLR: EWC must report the robust (scaled) fit that the pipeline and Mplus report

test_that("continuous EWC reports the robust MLR fit, like the pipeline", {
  ewc <- quiet(fit_ewc(res_c$fit_esem, spec_c))
  ct  <- compare_ewc(res_c, ewc)
  for (i in c("CFI", "TLI", "RMSEA", "SRMR", "X2"))
    expect_equal(ct$EWC_varfix[ct$Index == i], as.numeric(ct$ESEM_R[ct$Index == i]),
                 tolerance = 1e-3, label = i)
  x2 <- as.numeric(ct$ESEM_R[ct$Index == "X2"]); df <- as.numeric(ct$ESEM_R[ct$Index == "df"])
  expect_output(quiet(print(ewc)), sprintf("X2\\(%d\\) = %.3f", df, x2))
})

test_that("continuous ESEM and B-ESEM print the robust MLR chi-square the pipeline reports", {
  ct <- res_c$comparison_table
  x2 <- function(col) as.numeric(ct[[col]][ct$Index == "X2"])
  expect_output(quiet(print(res_c$fit_esem)),  sprintf("= %.3f", x2("ESEM_R")),  fixed = TRUE)
  expect_output(quiet(print(res_c$fit_besem)), sprintf("= %.3f", x2("BESEM_R")), fixed = TRUE)
})

# -- Custom ordered B-ESEM: SRMR by the Mplus definition (Asparouhov & Muthen
#    2018): polychoric residuals over n_pairs + total number of categories
#    (thresholds are saturated here, so their residuals are 0).

test_that("ordered B-ESEM (custom DWLS path) reports the Mplus-definition SRMR", {
  fb  <- res_o$fit_besem
  L   <- as.matrix(fb$rotated_loadings)
  Phi <- if (is.null(fb$factor_correlations)) diag(ncol(L)) else as.matrix(fb$factor_correlations)
  R   <- as.matrix(fb$polychoric)[rownames(L), rownames(L)]
  imp <- L %*% Phi %*% t(L)
  lt  <- lower.tri(R)
  n_cat <- sum(vapply(hs_ord[, rownames(L)], function(v) length(unique(v)), 1L))
  srmr_mplus <- sqrt(sum((R[lt] - imp[lt])^2) / (sum(lt) + n_cat))
  expect_equal(fb$wlsmv_stats$srmr, srmr_mplus, tolerance = 1e-6)
})
