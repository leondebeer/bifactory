# The custom DWLS paths (besem_ordered(), esem_ordered(method = "rotation"))
# fit the polychoric correlations that lavaan computed W and Gamma for.  Up to
# 0.6.0 they fitted psych::polychoric() instead, which differs from lavaan and
# Mplus when items have unequal numbers of categories (0.013 to 0.021 on the
# data below), so the chi-square did not reproduce in either program.  The
# B-ESEM and ESEM are rotations of the unrestricted k-factor model, so their
# fit must equal lavaan's own WLSMV fit of that model.

skip_on_cran()

quiet <- function(expr) suppressWarnings(suppressMessages(expr))

# G .6 and specific .5 loadings, 15 items, three specific factors; items
# alternate between 4 and 5 floor-skewed categories, so some bivariate
# tables have empty cells.
sim_mixed_cats <- function(seed = 1L, n = 800L) {
  set.seed(seed)
  g <- stats::rnorm(n); s <- matrix(stats::rnorm(3 * n), n, 3)
  c4 <- list(c(0.3, 1.2, 2.0), c(0.4, 1.3, 2.1), c(0.2, 1.1, 1.9))
  c5 <- list(c(0.3, 1.0, 1.7, 2.4), c(0.4, 1.1, 1.8, 2.5), c(0.2, 0.9, 1.6, 2.3))
  X <- list()
  for (j in 1:3) for (i in 1:5) {
    y    <- 0.6 * g + 0.5 * s[, j] + stats::rnorm(n, sd = sqrt(1 - 0.36 - 0.25))
    cuts <- if (i %% 2) c4[[(i %% 3) + 1]] else c5[[(i %% 3) + 1]]
    X[[paste0(letters[j], i)]] <- findInterval(y, cuts)
  }
  as.data.frame(X)
}
sf <- list(S1 = paste0("a", 1:5), S2 = paste0("b", 1:5), S3 = paste0("c", 1:5))
X  <- sim_mixed_cats()
R_lav <- lavaan::lavCor(X, ordered = names(X))

lav_efa_fit <- function(X, k) {
  lhs <- paste(sprintf('efa("e")*F%d', seq_len(k)), collapse = " + ")
  f   <- quiet(lavaan::cfa(paste(lhs, "=~", paste(names(X), collapse = " + ")),
                           data = X, ordered = names(X), estimator = "WLSMV",
                           parameterization = "theta", rotation = "geomin"))
  lavaan::fitMeasures(f, c("chisq.scaled", "df.scaled", "cfi.scaled", "rmsea.scaled"))
}

test_that("besem_ordered() fits lavaan's polychorics and reproduces lavaan's fit", {
  fit <- quiet(besem_ordered(X, sf, n_starts = 10L))
  expect_identical(fit$polychoric_source, "lavaan")
  expect_equal(unname(fit$polychoric), unname(unclass(R_lav)), tolerance = 1e-6)

  ref <- lav_efa_fit(X, 4L)
  expect_equal(fit$wlsmv_stats$df, unname(ref[["df.scaled"]]))
  expect_lt(abs(fit$wlsmv_stats$chisq - ref[["chisq.scaled"]]), 0.05)
  expect_lt(abs(fit$wlsmv_stats$cfi   - ref[["cfi.scaled"]]),   1e-3)
  expect_lt(abs(fit$wlsmv_stats$rmsea - ref[["rmsea.scaled"]]), 1e-3)
})

test_that("esem_ordered(method = 'rotation') fits lavaan's polychorics", {
  tgt <- make_target(sf, item_names = names(X))
  fit <- quiet(esem_ordered(X, nfactors = 3L, target = tgt, method = "rotation",
                            n_starts = 10L))
  expect_identical(fit$polychoric_source, "lavaan")
  expect_equal(unname(fit$polychoric), unname(unclass(R_lav)), tolerance = 1e-6)

  ref <- lav_efa_fit(X, 3L)
  expect_equal(fit$wlsmv_stats$df, unname(ref[["df.scaled"]]))
  expect_lt(abs(fit$wlsmv_stats$chisq - ref[["chisq.scaled"]]), 0.05)
})

test_that("r_obs_override is used as given and warned about when it departs from lavaan", {
  R_psych <- suppressWarnings(psych::polychoric(X)$rho)
  expect_warning(
    fit <- suppressMessages(besem_ordered(X, sf, n_starts = 2L, r_obs_override = R_psych)),
    "differs from lavaan")
  expect_identical(fit$polychoric_source, "r_obs_override")
  expect_equal(unname(fit$polychoric), unname(R_psych))

  # lavaan's own matrix passes without a warning about the source
  expect_no_warning(
    suppressMessages(besem_ordered(X, sf, n_starts = 2L, r_obs_override = unclass(R_lav))),
    message = "differs from lavaan")
})

test_that("multi-group besem_ordered() pools the polychorics and keeps the groups", {
  Xg <- cbind(X, grp = rep(c("m", "f"), length.out = nrow(X)))
  fit <- quiet(besem_ordered(Xg, sf, n_starts = 2L, group = "grp"))
  expect_equal(unname(fit$polychoric), unname(unclass(R_lav)), tolerance = 1e-6)
  expect_equal(lavaan::lavInspect(fit$lavaan_fit, "ngroups"), 2L)
  fs <- quiet(factor_scores(fit))
  expect_true("group" %in% names(fs))
})
