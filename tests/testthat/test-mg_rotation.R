# From the weak level on the loadings are equal across groups and the model
# is rotated once.  lavaan's default scales the rows of the common loading
# matrix by diag(AA') plus the residual variances averaged over the groups
# before rotating; Mplus rotates them unscaled.  Same fit, a different
# rotation of it.  .mg_rotation_args() passes std.ov = FALSE in that case.

test_that(".mg_rotation_args() drops row scaling only for equal loadings across groups", {
  base <- list(rotation.args = list(target = matrix(0, 2, 2)))
  mg   <- .mg_rotation_args(c(base, list(group = "g", group.equal = "loadings")))
  expect_false(mg$rotation.args$std.ov)
  expect_true(is.matrix(mg$rotation.args$target))
  expect_null(.mg_rotation_args(c(base, list(group = "g")))$rotation.args$std.ov)        # configural
  expect_null(.mg_rotation_args(c(base, list(group.equal = "loadings")))$rotation.args$std.ov)
  expect_null(.mg_rotation_args(c(base, list(group = "g", group.equal = "intercepts")))$rotation.args$std.ov)
  user <- list(rotation.args = list(std.ov = TRUE), group = "g", group.equal = "loadings")
  expect_true(.mg_rotation_args(user)$rotation.args$std.ov)                               # caller's choice wins
  snake <- list(rotation.args = list(std_ov = TRUE), group = "g", group.equal = "loadings")
  expect_identical(.mg_rotation_args(snake), snake)                                      # no duplicate spelling
  one <- list(data = data.frame(g = c("a", "a", NA)), group = "g", group.equal = "loadings")
  expect_null(.mg_rotation_args(one)$rotation.args$std.ov)                                # one group in lavaan
  two <- list(data = data.frame(g = c("a", "b")), group = "g", group.equal = "loadings")
  expect_false(.mg_rotation_args(two)$rotation.args$std.ov)
  mat <- list(data = cbind(x = 1:2, g = 1:2), group = "g", group.equal = "loadings")          # lavaan takes a matrix
  expect_false(.mg_rotation_args(mat)$rotation.args$std.ov)
})

# Oblique target rotation: delta-method rotation SEs with a tighter GPA stop;
# orthogonal target (B-ESEM) and other rotations keep lavaan's bordered default.
test_that(".mg_rotation_args() picks delta rotation SEs for oblique target rotation only", {
  tgt <- list(rotation = "target", rotation.args = list(target = matrix(0, 2, 2)))
  ob  <- .mg_rotation_args(tgt)
  expect_identical(ob$rotation.se, "delta")
  expect_identical(ob$rotation.args$gpa.tol, 1e-6)
  expect_identical(.mg_rotation_args(utils::modifyList(tgt, list(rotation = "Target")))$rotation.se, "delta")
  expect_null(.mg_rotation_args(utils::modifyList(tgt, list(rotation.args = list(orthogonal = TRUE))))$rotation.se)
  expect_null(.mg_rotation_args(list(rotation = "geomin"))$rotation.se)
  expect_identical(.mg_rotation_args(c(tgt, list(rotation.se = "bordered")))$rotation.se, "bordered")
  snake <- .mg_rotation_args(utils::modifyList(tgt, list(rotation.args = list(gpa_tol = 1e-5))))
  expect_null(snake$rotation.args$gpa.tol)                                                # no duplicate spelling
})

# lavaan accepts group.partial on an efa loading but reports group 1's rotated
# loading in every group, so a freed loading would be shown as equal.
test_that("group.partial on a factor loading is refused", {
  expect_error(.mg_rotation_args(list(group = "g", group.equal = "loadings", group.partial = "F1=~x1")),
               "group.partial")
  expect_silent(.check_group_partial(c("x1|t1", "x2~~x2", "x3~1")))
  hs <- lavaan::HolzingerSwineford1939
  spec <- suppressMessages(specify_model(visual = c("x1", "x2", "x3"), textual = c("x4", "x5", "x6"),
                                         speed = c("x7", "x8", "x9"), data = hs, group = "school"))
  expect_error(esem_invariance(spec, group.partial = "visual=~x1"), "group.partial")
  expect_false(.mg_rotation_args(list(group = "g", group.equal = c("loadings", "thresholds")))$rotation.args$std.ov)
})

# Mplus 9 STDYX loadings, weak level, HS by school (Grant-White = group 1),
# ESEM with target rotation (validation/through_mplus/esem_{mlr,wlsmv}_weak.out).
# With lavaan's default rotation R was 0.0107 (MLR) and 0.0261 (WLSMV) away.
hs_mplus_weak <- list(
  mlr = list(
    matrix(c( 0.537,  0.182,  0.090,   0.540, -0.014, -0.035,   0.732, -0.021,  0.031,
              0.035,  0.848, -0.026,  -0.052,  0.880,  0.022,   0.101,  0.735,  0.000,
             -0.243,  0.097,  0.785,   0.041, -0.109,  0.844,   0.304,  0.004,  0.548), 9, 3, byrow = TRUE),
    matrix(c( 0.579,  0.192,  0.075,   0.508, -0.013, -0.026,   0.693, -0.019,  0.023,
              0.036,  0.831, -0.020,  -0.051,  0.849,  0.017,   0.115,  0.818,  0.000,
             -0.242,  0.095,  0.604,   0.046, -0.120,  0.729,   0.325,  0.005,  0.452), 9, 3, byrow = TRUE)),
  wlsmv = list(
    matrix(c( 0.538,  0.189,  0.109,   0.516,  0.105, -0.050,   0.679, -0.052,  0.052,
              0.052,  0.760,  0.013,   0.007,  0.899,  0.008,  -0.031,  0.846, -0.025,
             -0.247,  0.052,  0.772,  -0.002, -0.055,  0.834,   0.346,  0.030,  0.535), 9, 3, byrow = TRUE),
    matrix(c( 0.443,  0.270,  0.095,   0.356,  0.126, -0.037,   0.856, -0.115,  0.070,
              0.032,  0.811,  0.009,   0.004,  0.858,  0.005,  -0.019,  0.916, -0.017,
             -0.216,  0.079,  0.718,  -0.001, -0.063,  0.584,   0.295,  0.044,  0.486), 9, 3, byrow = TRUE)))

# Mplus 9 STDYX loading SEs, configural level, same models
# (validation/through_mplus/esem_{mlr,wlsmv}_configural.out).  With lavaan's
# bordered rotation SEs R was up to 0.014 away; with delta 0.0024.
hs_mplus_cfg_se <- list(
  mlr = list(
    matrix(c(0.133, 0.069, 0.070,  0.129, 0.109, 0.073,  0.132, 0.077, 0.068,
             0.058, 0.050, 0.045,  0.072, 0.058, 0.056,  0.057, 0.050, 0.048,
             0.069, 0.049, 0.107,  0.082, 0.050, 0.088,  0.107, 0.074, 0.073), 9, 3, byrow = TRUE),
    matrix(c(0.099, 0.078, 0.062,  0.109, 0.075, 0.069,  0.076, 0.062, 0.057,
             0.058, 0.057, 0.080,  0.080, 0.063, 0.049,  0.055, 0.040, 0.043,
             0.106, 0.091, 0.393,  0.181, 0.095, 0.299,  0.119, 0.099, 0.251), 9, 3, byrow = TRUE)),
  wlsmv = list(
    matrix(c(0.127, 0.090, 0.076,  0.124, 0.110, 0.092,  0.116, 0.089, 0.085,
             0.081, 0.061, 0.061,  0.077, 0.069, 0.058,  0.072, 0.064, 0.053,
             0.078, 0.057, 0.111,  0.103, 0.060, 0.109,  0.089, 0.058, 0.085), 9, 3, byrow = TRUE),
    matrix(c(0.113, 0.085, 0.078,  0.112, 0.090, 0.092,  0.112, 0.068, 0.070,
             0.061, 0.054, 0.061,  0.054, 0.056, 0.055,  0.051, 0.048, 0.050,
             0.085, 0.070, 0.139,  0.105, 0.067, 0.124,  0.088, 0.073, 0.114), 9, 3, byrow = TRUE)))

test_that("weak-level ESEM loadings match Mplus without alignment (MLR and WLSMV)", {
  skip_on_cran()
  hs <- lavaan::HolzingerSwineford1939
  hs$school <- as.integer(hs$school)
  fl <- list(Vis = c("x1", "x2", "x3"), Tex = c("x4", "x5", "x6"), Spd = c("x7", "x8", "x9"))
  items <- unname(unlist(fl))
  hs_ord <- hs
  for (it in items)
    hs_ord[[it]] <- cut(hs[[it]], stats::quantile(hs[[it]], c(0, .25, .5, .75, 1)),
                        include.lowest = TRUE, labels = FALSE)
  for (est in c("mlr", "wlsmv")) {
    spec <- do.call(specify_model, c(fl, list(data = if (est == "mlr") hs else hs_ord,
                                              group = "school", ordered = est == "wlsmv")))
    inv <- suppressWarnings(suppressMessages(esem_invariance(spec, verbose = FALSE, cores = 1L)))
    ss  <- lavaan::standardizedSolution(inv$models$weak$lavaan_fit, se = FALSE)
    ss  <- ss[ss$op == "=~", ]
    for (g in 1:2) {
      s <- ss[ss$group == g, ]
      R <- matrix(NA_real_, 9, 3, dimnames = list(items, names(fl)))
      R[cbind(s$rhs, s$lhs)] <- s$est.std
      expect_lt(max(abs(R - hs_mplus_weak[[est]][[g]])), 0.003)
    }
    cs <- lavaan::standardizedSolution(inv$models$configural$lavaan_fit)
    cs <- cs[cs$op == "=~", ]
    for (g in 1:2) {
      s <- cs[cs$group == g, ]
      R <- matrix(NA_real_, 9, 3, dimnames = list(items, names(fl)))
      R[cbind(s$rhs, s$lhs)] <- s$se
      expect_lt(max(abs(R - hs_mplus_cfg_se[[est]][[g]])), 0.003)
    }
  }
})
