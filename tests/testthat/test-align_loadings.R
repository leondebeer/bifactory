# extract_mplus_loadings(): when Mplus cannot compute standard errors it prints the
# standardized results as a single "StdYX Estimate" column, which MplusAutomation
# returns as character; the loading matrices must still be numeric.

test_that("extract_mplus_loadings returns numeric matrices when Mplus prints estimates as text", {
  ps <- data.frame(paramHeader = c("G.BY", "G.BY", "EX.BY", "Residual.Variances"),
                   param       = c("EX1", "EX2", "EX1", "CFA MODEL"),
                   est         = c("0.870", "0.822", "0.103", "g"),
                   Group       = "BE", stringsAsFactors = FALSE)
  M <- extract_mplus_loadings(list(parameters = list(stdyx.standardized = ps)),
                              items = c("EX1", "EX2"), factors = c("G", "EX"))
  expect_true(is.numeric(M[["1"]]))
  expect_equal(M[["1"]]["EX1", "G"], 0.870)
  expect_equal(M[["1"]]["EX2", "EX"], 0)
})

# Orientation transfer: from the weak level on the non-reference groups have free factor
# covariances, and the standardization divides by sqrt(diag(Q' Psi_g Q)), which no rotation
# of the standardized matrix can undo. Comparing against another solution therefore has to
# rotate Psi_g along with the shared loadings before standardizing.
.hs_besem_inv <- function() {
  data(HolzingerSwineford1939, package = "lavaan"); hs <- HolzingerSwineford1939
  for (it in paste0("x", 1:9))
    hs[[it]] <- cut(hs[[it]], quantile(hs[[it]], c(0, .25, .5, .75, 1)), include.lowest = TRUE, labels = FALSE)
  spec <- specify_model(Vis = c("x1", "x2", "x3"), Tex = c("x4", "x5", "x6"), Spd = c("x7", "x8", "x9"),
                        data = hs, group = "school", ordered = TRUE, label = "HS ord")
  suppressWarnings(esem_invariance(spec, model = "besem", through = "strict", verbose = FALSE, cores = 1L))
}

test_that("a solution as target is compared by orientation transfer, exact in every group", {
  skip_on_cran()
  inv <- .hs_besem_inv()
  sol <- .extract_solution(inv, level = "weak")
  k <- ncol(sol[[1]]$lambda)
  set.seed(1); Q <- qr.Q(qr(matrix(rnorm(k * k), k)))
  # the same solution in another orientation: shared loadings rotated, Psi_g transformed
  rotated <- lapply(sol, function(s) list(lambda = s$lambda %*% Q, psi = t(Q) %*% s$psi %*% Q, theta = s$theta))
  # old yardstick: group-wise Procrustes on the standardized loadings of the rotated solution
  old <- align_loadings(inv, target = lapply(rotated, .stdyx_from_solution), level = "weak", se_method = "none")
  new <- align_loadings(inv, target = rotated, level = "weak", se_method = "none")
  expect_identical(attr(old, "method"), "procrustes")
  expect_identical(attr(new, "method"), "transfer")
  expect_lt(old[[1]]$residual_max, 1e-6)       # reference group (Psi = I): exact either way
  expect_gt(old[[2]]$residual_max, 1e-3)       # non-reference group: Procrustes cannot undo it
  expect_lt(new[[1]]$residual_max, 1e-8)
  expect_lt(new[[2]]$residual_max, 1e-8)
  expect_equal(new[[2]]$Q, new[[1]]$Q)         # one rotation for the shared loadings
  # a fit compared with itself: identity rotation, zero residual
  self <- align_loadings(inv, target = inv, level = "weak", se_method = "none")
  expect_identical(attr(self, "method"), "transfer")
  expect_equal(unname(self[[2]]$Q), diag(k), tolerance = 1e-8)
  expect_lt(self[[2]]$residual_max, 1e-10)
})
