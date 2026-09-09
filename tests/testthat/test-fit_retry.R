# .fit_with_retry(): lavaan is started twice (default and "simple" start values)
# and the converged fit with the lower fit function is kept.  On the BAT 2023
# data the default start left the two-country configural B-ESEM in a well 9 %
# shallower than the one Mplus reaches (NL loadings 0.21 off); the simple start
# finds the deeper well.
# .mplus_conv_lines(): the generator's first Mplus attempt runs at a criterion
# tighter than Mplus's default (.00005 stopped the eight-country means model 15
# chi-square units early) and the formatter must write such criteria in full.

hs         <- lavaan::HolzingerSwineford1939
saturated  <- "f =~ x1 + x2 + x3"         # three indicators: fit function 0
misfit     <- "f =~ x1 + x2 + x3 + x7"    # x7 belongs to another factor: fit function > 0
fx         <- function(fit) fit$lavaan_fit@optim$fx

# A fit function that returns the saturated model for one start and the
# misfitting model for the other (attempt 0 passes no `start`, i.e. default).
fit_by_start <- function(saturated_start) {
  function(ctrl = NULL, opts = NULL) {
    st <- if (is.null(opts$start)) "default" else opts$start
    m  <- if (identical(st, saturated_start)) saturated else misfit
    list(lavaan_fit = lavaan::cfa(m, data = hs))
  }
}

test_that(".fit_with_retry keeps the start with the lower fit function", {
  expect_equal(fx(.fit_with_retry(fit_by_start("simple"),  verbose = FALSE)), 0, tolerance = 1e-8)
  expect_equal(fx(.fit_with_retry(fit_by_start("default"), verbose = FALSE)), 0, tolerance = 1e-8)
})

# Equally deep wells: the eight-country BAT 2023 means model has several optima whose fit
# functions differ by 0.003 to 0.02 percent (section 5c of the BAT 2023 note) and whose
# loadings differ by 0.11.  Switching to the alternative start for such a gain only moves
# the reported solution away from the default start's well (which is where Mplus lands);
# the alternative must be materially deeper to replace the default fit.
test_that(".fit_with_retry ignores an alternative start that is not materially deeper", {
  f <- lavaan::cfa(misfit, data = hs)
  deeper_by <- function(rel) { g <- f; g@optim$fx <- f@optim$fx * (1 - rel); g }
  fit_fn <- function(rel) function(ctrl = NULL, opts = NULL) {
    if (is.null(opts$start)) list(lavaan_fit = f, tag = "default")
    else list(lavaan_fit = deeper_by(rel), tag = "alt")
  }
  expect_identical(.fit_with_retry(fit_fn(1e-4), verbose = FALSE)$tag, "default")   # 0.01 %: a tie
  expect_identical(.fit_with_retry(fit_fn(1e-2), verbose = FALSE)$tag, "alt")       # 1 %: a deeper well
})

test_that(".mplus_conv_lines writes a tight first criterion and formats small criteria in full", {
  expect_match(.mplus_conv_lines(),        "CONVERGENCE = 0.000001;", fixed = TRUE)
  expect_match(.mplus_conv_lines(1e-7),    "CONVERGENCE = 0.0000001;", fixed = TRUE)
  expect_match(.mplus_conv_lines(2.5e-2),  "CONVERGENCE = 0.025;", fixed = TRUE)
})

test_that(".fit_with_retry(alt_start = NULL) fits once", {
  n <- 0L
  fit_fn <- function(ctrl = NULL, opts = NULL) { n <<- n + 1L; list(lavaan_fit = lavaan::cfa(saturated, data = hs)) }
  .fit_with_retry(fit_fn, verbose = FALSE, alt_start = NULL)
  expect_identical(n, 1L)
})

# The configural model has no cross-group constraint, so its optimum is the sum of the
# per-group optima; the joint optimizer's path depends on the number of groups (four BAT
# 2023 countries: NL stuck in a well 9% shallower from both starts, alone it finds the
# deeper one in seconds).  The ordered configural fit therefore starts from the per-group
# fits (each from both starts, deeper kept).
hs_ord <- hs
for (it in paste0("x", 1:9))
  hs_ord[[it]] <- cut(hs[[it]], quantile(hs[[it]], c(0, .25, .5, .75, 1)), include.lowest = TRUE, labels = FALSE)
spec_o <- specify_model(Vis = c("x1", "x2", "x3"), Tex = c("x4", "x5", "x6"), Spd = c("x7", "x8", "x9"),
                        data = hs_ord, group = "school", ordered = TRUE, label = "HS ord")

test_that("ordered configural fit starts from the per-group optima and is never shallower", {
  x0 <- .configural_starts(spec_o, model = "besem", missing = "pairwise")
  fit <- suppressWarnings(.fit_invariance_model(spec_o, group_equal = NULL, is_ordered = TRUE,
                                                model = "besem", missing = "pairwise", verbose = FALSE))
  lf <- fit$lavaan_fit
  # A numeric `start` is capped at 1000 entries by lavaan (four bfi education tiers need
  # 1040, eight BAT countries 1576), so the starts travel as a parameter table (lhs/op/rhs/
  # group/est, matched row by row) holding the unrotated per-group optima.
  expect_s3_class(x0, "data.frame")
  expect_setequal(x0$group, seq_along(spec_o$group_levels))
  expect_identical(sum(x0$free > 0L), length(lf@optim$x))     # unrotated free parameters, all groups
  expect_true(all(is.finite(x0$est[x0$free > 0L])))
  expect_true(is.data.frame(lf@Options$start))                 # the joint fit really started there
  # HS quartile-cut B-ESEM is degenerate in Pasteur (a flat valley), so the joint fit may go
  # deeper than the per-group optimum there; it must never end shallower in any group.
  expect_true(all(lf@optim$fx.group <= attr(x0, "fx") + 1e-8))
})

# When 1e-6 fails (eight-country weak/strong/strict), the first retry must not be looser
# than the Mplus default at which those models used to converge.
test_that("the Mplus retry sequence starts at Mplus's default criterion", {
  expect_equal(.MPLUS_CONV_RETRY[1], 5e-5)
  expect_match(.mplus_conv_lines(.MPLUS_CONV_RETRY[1]), "CONVERGENCE = 0.00005;", fixed = TRUE)
})

# A level whose Mplus run printed no fit (eight-country var/cov: singular information
# matrix) saves no DIFFTEST derivatives file; the next level's input must then omit
# DIFFTEST instead of failing with "The file specified for the DIFFTEST option cannot be found".
test_that(".mplus_difftest_line omits DIFFTEST when the previous level saved no file", {
  dat <- file.path(tempdir(), "besem_inv_varcov.dat")
  unlink(dat)
  expect_message(line <- .mplus_difftest_line(dat), "DIFFTEST")
  expect_identical(line, "")
  writeLines("x", dat)
  expect_silent(line <- .mplus_difftest_line(dat))
  expect_identical(line, "\n  DIFFTEST = besem_inv_varcov.dat;")
  unlink(dat)
})
