# partial_invariance(): only partial strong and partial strict exist for ESEM
# (the rotated loading block is invariant as a unit); the result table must
# assemble without error; and a "partial" refit must actually release the
# reported constraints (fewer df, lower chi-square than the full model) on
# every path: continuous and ordered, ESEM and B-ESEM.

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
inv    <- quiet(esem_invariance(spec_c, verbose = FALSE))

df_of  <- function(fit) unname(lavaan::fitMeasures(fit$lavaan_fit, "df"))
x2_of  <- function(fit) unname(lavaan::fitMeasures(fit$lavaan_fit, "chisq"))
expect_released <- function(p, full) {
  n_freed <- nrow(p$freed_params)
  expect_gt(n_freed, 0L)
  expect_equal(df_of(full) - df_of(p$partial_fit), n_freed)
  expect_lt(x2_of(p$partial_fit), x2_of(full))
}

test_that("partial weak is refused: ESEM loadings are invariant as a block", {
  expect_error(partial_invariance(inv, level = "weak", verbose = FALSE), "strong")
})

test_that("partial strong returns a table with the invariance table's columns", {
  p <- quiet(partial_invariance(inv, level = "strong", max_free = 2L, verbose = FALSE))
  expect_s3_class(p, "esem_partial_invariance")
  expect_identical(names(p$table), c(names(inv$table), "pass"))
  expect_true(any(grepl("partial", p$table$Model)))
  expect_output(quiet(print(p)), "partial", fixed = TRUE)
})

test_that("continuous ESEM: the partial refit really releases the freed intercepts", {
  p <- quiet(partial_invariance(inv, level = "strong", max_free = 2L, verbose = FALSE))
  expect_released(p, inv$models$strong)
  expect_true(all(grepl("intercept", p$freed_params$label)))
})

inv_o  <- quiet(esem_invariance(spec_o, verbose = FALSE))

test_that("ordered ESEM: freed thresholds and residuals are really released", {
  p_strong <- quiet(partial_invariance(inv_o, level = "strong", max_free = 2L, verbose = FALSE))
  expect_released(p_strong, inv_o$models$strong)
  expect_true(all(grepl("threshold", p_strong$freed_params$label)))
  # the downstream strict refit keeps the released thresholds free
  expect_lt(df_of(p_strong$downstream$strict), df_of(inv_o$models$strict))
  p_strict <- quiet(partial_invariance(inv_o, level = "strict", max_free = 1L, verbose = FALSE))
  expect_released(p_strict, inv_o$models$strict)
  expect_true(all(grepl("residual", p_strict$freed_params$label)))
})

test_that("B-ESEM is supported on both paths", {
  for (spec in list(spec_c, spec_o)) {
    binv <- quiet(esem_invariance(spec, model = "besem", verbose = FALSE))
    p <- quiet(partial_invariance(binv, level = "strong", max_free = 1L, verbose = FALSE))
    expect_s3_class(p, "esem_partial_invariance")
    expect_released(p, binv$models$strong)
  }
})

test_that("two invariant anchors per factor are kept and heavy release warns", {
  expect_warning(
    p <- suppressMessages(partial_invariance(inv_o, level = "strong", max_free = 9L,
                                             delta_cfi_cutoff = 1, verbose = FALSE)),
    "20%")
  freed_items <- sub(" .*$", "", p$freed_params$label)
  for (f in fl) expect_gte(sum(!(f %in% freed_items)), 2L)
  expect_released(p, inv_o$models$strong)
})

test_that("a hand-specified group.partial releases an item's thresholds and re-fixes its residual", {
  gp  <- paste0("x9|t", 1:3)
  inv <- quiet(esem_invariance(spec_o, verbose = FALSE, group.partial = gp))
  pt  <- lavaan::parameterTable(inv$models$strong$lavaan_fit)
  thr <- pt[pt$op == "|" & pt$lhs == "x9", ]
  expect_true(all(thr$free > 0))                      # free in both groups
  expect_equal(sum(pt$op == "==" & pt$lhs %in% thr$plabel), 0L)
  res <- pt[pt$op == "~~" & pt$lhs == "x9" & pt$rhs == "x9", ]
  expect_equal(res$free, c(0L, 0L)); expect_equal(res$ustart, c(1, 1))
  expect_equal(df_of(inv_o$models$strong) - df_of(inv$models$strong), 3L - 1L)
})

# -- Recovery of planted non-invariance ------------------------------------------
# Two groups, 3 factors x 4 items; group 2 has the intercepts of y2 (+0.6) and y7
# (-0.6) shifted. The release rule must free exactly those items (as intercepts on
# the continuous path, as thresholds on the ordered path) and nothing else.
test_that("planted non-invariance is recovered on both paths", {
  set.seed(2026)
  n <- 600
  L <- matrix(0.15, 12, 3); for (f in 1:3) L[(f - 1) * 4 + 1:4, f] <- 0.7
  gen <- function(n, shift) {
    eta <- MASS::mvrnorm(n, rep(0, 3), matrix(c(1, .3, .3, .3, 1, .3, .3, .3, 1), 3))
    y <- eta %*% t(L) + matrix(rnorm(n * 12, sd = sqrt(1 - rowSums(L^2))), n, 12, byrow = TRUE)
    y <- sweep(y, 2, shift, "+"); colnames(y) <- paste0("y", 1:12); as.data.frame(y)
  }
  shift2 <- rep(0, 12); shift2[2] <- 0.6; shift2[7] <- -0.6
  d  <- rbind(cbind(gen(n, rep(0, 12)), g = 1L), cbind(gen(n, shift2), g = 2L))
  fs <- list(F1 = paste0("y", 1:4), F2 = paste0("y", 5:8), F3 = paste0("y", 9:12))
  item_of <- function(p) sub(" .*$", "", p$freed_params$label)

  inv_s <- quiet(esem_invariance(do.call(specify_model, c(fs, list(data = d, group = "g", label = "sim"))), verbose = FALSE))
  p_s <- quiet(partial_invariance(inv_s, level = "strong", verbose = FALSE))
  expect_true(p_s$converged)
  expect_setequal(item_of(p_s), c("y2", "y7"))
  expect_true(all(grepl("intercept", p_s$freed_params$label)))

  do <- d; for (v in paste0("y", 1:12)) do[[v]] <- as.integer(cut(d[[v]], c(-Inf, -1.2, -0.4, 0.4, 1.2, Inf)))
  inv_so <- quiet(esem_invariance(do.call(specify_model, c(fs, list(data = do, group = "g", ordered = TRUE, label = "sim ord"))), verbose = FALSE))
  p_so <- quiet(partial_invariance(inv_so, level = "strong", verbose = FALSE))
  expect_true(p_so$converged)
  expect_setequal(unique(item_of(p_so)), c("y2", "y7"))
  expect_true(all(grepl("threshold", p_so$freed_params$label)))
})
