# esem_invariance(): the default sequence stops at strict; `through` adds the
# latent variance/covariance and latent mean levels (Morin's sequence).
# HolzingerSwineford1939 by school (2 groups), continuous and quartile-cut.
# Each invariance run is slow, so the six-level runs are fitted once and shared.

skip_on_cran()

hs <- lavaan::HolzingerSwineford1939
fl <- list(Vis = c("x1", "x2", "x3"), Tex = c("x4", "x5", "x6"), Spd = c("x7", "x8", "x9"))
hs_ord <- hs
for (it in paste0("x", 1:9))
  hs_ord[[it]] <- cut(hs[[it]], quantile(hs[[it]], c(0, .25, .5, .75, 1)),
                      include.lowest = TRUE, labels = FALSE)

quiet  <- function(expr) suppressWarnings(suppressMessages(expr))
spec_c <- do.call(specify_model, c(fl, list(data = hs,     group = "school", ordered = FALSE, label = "HS")))
spec_o <- do.call(specify_model, c(fl, list(data = hs_ord, group = "school", ordered = TRUE,  label = "HS ord")))
k  <- length(fl)        # ESEM factors
kb <- k + 1L            # B-ESEM: G + specific factors

inv_c  <- quiet(esem_invariance(spec_c, through = "means", verbose = FALSE))
inv_o  <- quiet(esem_invariance(spec_o, through = "means", verbose = FALSE))
binv_c <- quiet(esem_invariance(spec_c, model = "besem", through = "means", verbose = FALSE))
binv_o <- quiet(esem_invariance(spec_o, model = "besem", through = "means", verbose = FALSE))

df_of <- function(inv, lv) unname(lavaan::fitMeasures(inv$models[[lv]]$lavaan_fit, "df"))
ptab  <- function(inv, lv) lavaan::parameterTable(inv$models[[lv]]$lavaan_fit)
facs  <- function(pt) unique(pt$lhs[pt$op == "=~"])
row_of <- function(inv, label) inv$table[inv$table$Model == label, ]

# -- Default unchanged ----------------------------------------------------------

test_that("the default sequence still stops at strict (guard)", {
  inv <- quiet(esem_invariance(spec_c, verbose = FALSE))
  expect_identical(names(inv$models), c("configural", "weak", "strong", "strict"))
  expect_equal(nrow(inv$table), 4L)
})

test_that("an unknown `through` value is rejected", {
  expect_error(esem_invariance(spec_c, through = "latent", verbose = FALSE), "should be one of")
})

# -- Level 5: latent variances and covariances ---------------------------------

test_that("through = 'varcov' adds one level with k(k+1)/2 more df than strict", {
  inv <- quiet(esem_invariance(spec_c, through = "varcov", verbose = FALSE))
  expect_identical(names(inv$models), c("configural", "weak", "strong", "strict", "varcov"))
  expect_equal(nrow(inv$table), 5L)
  expect_equal(df_of(inv, "varcov") - df_of(inv, "strict"), k * (k + 1) / 2)
  expect_equal(row_of(inv, "5. Latent var/cov")$delta_df, k * (k + 1) / 2)   # LRT vs strict
})

test_that("ESEM varcov level: non-reference variances 1 and covariances equal to group 1", {
  for (inv in list(inv_c, inv_o)) {
    expect_equal(df_of(inv, "varcov") - df_of(inv, "strict"), k * (k + 1) / 2)
    pt <- ptab(inv, "varcov"); f <- facs(pt)
    cv <- pt[pt$op == "~~" & pt$lhs %in% f & pt$rhs %in% f, ]
    v2 <- cv[cv$lhs == cv$rhs & cv$group == 2L, ]
    expect_equal(v2$est, rep(1, k))
    c1 <- cv[cv$lhs != cv$rhs & cv$group == 1L, ]
    c2 <- cv[cv$lhs != cv$rhs & cv$group == 2L, ]
    expect_equal(nrow(c2), k * (k - 1) / 2)
    # equal up to optimizer precision (lavaan imposes the efa-block equality numerically)
    expect_equal(c2$est[order(c2$lhs, c2$rhs)], c1$est[order(c1$lhs, c1$rhs)], tolerance = 1e-4)
  }
})

test_that("B-ESEM varcov level: non-reference variances 1 and all covariances 0", {
  for (inv in list(binv_c, binv_o)) {
    expect_equal(df_of(inv, "varcov") - df_of(inv, "strict"), kb * (kb + 1) / 2)
    pt <- ptab(inv, "varcov"); f <- facs(pt)
    cv <- pt[pt$op == "~~" & pt$lhs %in% f & pt$rhs %in% f & pt$group == 2L, ]
    expect_equal(cv$est[cv$lhs == cv$rhs], rep(1, kb))
    expect_equal(cv$est[cv$lhs != cv$rhs], rep(0, kb * (kb - 1) / 2))
    expect_true(all(cv$free == 0L))
  }
})

# -- Level 6: latent means ---------------------------------------------------------

test_that("through = 'means' adds a sixth level with k more df than varcov", {
  for (inv in list(inv_c, inv_o)) {
    expect_identical(names(inv$models),
                     c("configural", "weak", "strong", "strict", "varcov", "means"))
    expect_equal(df_of(inv, "means") - df_of(inv, "varcov"), k)
    expect_equal(row_of(inv, "6. Latent means")$delta_df, k)
    pt <- ptab(inv, "means"); f <- facs(pt)
    m2 <- pt[pt$op == "~1" & pt$lhs %in% f & pt$group == 2L, ]
    expect_equal(m2$est, rep(0, k))
  }
  for (inv in list(binv_c, binv_o)) {
    expect_equal(df_of(inv, "means") - df_of(inv, "varcov"), kb)
    pt <- ptab(inv, "means"); f <- facs(pt)
    m2 <- pt[pt$op == "~1" & pt$lhs %in% f & pt$group == 2L, ]
    expect_equal(m2$est, rep(0, kb))
  }
})

test_that("print lists the extra levels", {
  expect_output(print(inv_c), "5. Latent var/cov", fixed = TRUE)
  expect_output(print(inv_c), "6. Latent means",   fixed = TRUE)
  expect_output(print(inv_c), "inv$models$means",  fixed = TRUE)
})

# -- Downstream consumers --------------------------------------------------------

test_that("factor_scores() auto-selects among all fitted levels and accepts the new levels", {
  sel <- quiet(bifactory:::.fs_select_invariance_level(inv_c, "auto"))
  expect_identical(sel$level, bifactory:::.inv_supported_level(inv_c$table, inv_c$models, notes = inv_c$notes))
  expect_true(sel$level %in% names(inv_c$models))
  fs <- quiet(factor_scores(inv_c, level = "means"))
  expect_s3_class(fs, "data.frame")
  expect_true("group" %in% names(fs))
})

# -- SRMR: Mplus definition (Asparouhov & Muthen 2018) ---------------------------
# Mplus 9 values for the same six models (validation/_validate_through_mplus.R).
# Continuous: correlation + variance + standardized-mean residuals; ordered:
# polychoric + category-probability residuals; groups pooled by sample size.

test_that("multi-group SRMR matches Mplus at every level (continuous, MLR)", {
  mplus <- c(0.021, 0.044, 0.057, 0.079, 0.094, 0.130)
  expect_equal(round(inv_c$table$SRMR, 3), mplus)
})

test_that("multi-group SRMR matches Mplus at every level (ordered, WLSMV)", {
  mplus <- c(0.026, 0.035, 0.042, 0.044, 0.055, 0.056)
  expect_equal(round(inv_o$table$SRMR, 3), mplus)
})

# -- Reference group -------------------------------------------------------------
# lavaan orders groups by first appearance in the data; specify_model() sorts
# the levels and the Mplus generator takes the first sorted level as the
# reference group. The lavaan fits must use the same order, otherwise R and
# Mplus fix orthogonality in different groups and report differently oriented
# loadings (HS rows start with Pasteur; sorted levels start with Grant-White).
test_that("lavaan group order follows spec$group_levels, not data order", {
  for (inv in list(inv_c, inv_o, binv_c, binv_o))
    for (lv in names(inv$models))
      expect_identical(lavaan::lavInspect(inv$models[[lv]]$lavaan_fit, "group.label"),
                       as.character(spec_c$group_levels))
})

# -- Ordered configural identification ------------------------------------------
# The configural model uses the standard ordinal identification in every group
# (residual variances 1, factor means 0, all thresholds free). Freeing the
# non-reference residuals and means against one equated threshold per item
# (the weak-level scheme) leaves each item's scale nearly unidentified when that
# threshold sits near zero in a non-reference group: on bfi (G = 3) the O5
# residual variance ran to 333 and the chi-square drifted 14 units from Mplus.
test_that("ordered configural model fixes residuals and means in every group, thresholds free", {
  for (inv in list(inv_o, binv_o)) {
    pt <- ptab(inv, "configural")
    rv <- pt[pt$op == "~~" & pt$lhs == pt$rhs & pt$lhs %in% spec_o$all_items, ]
    expect_true(all(rv$free == 0L & rv$ustart == 1))
    fm <- pt[pt$op == "~1" & !pt$lhs %in% spec_o$all_items, ]
    expect_true(all(fm$free == 0L & fm$ustart == 0))
    expect_true(all(pt$free[pt$op == "|"] > 0L))
    expect_identical(sum(pt$op == "=="), 0L)
  }
})
