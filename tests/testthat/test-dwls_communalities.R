# The custom DWLS paths (besem_ordered(), esem_ordered(method = "rotation"))
# fit only the polychoric correlations between items, so nothing keeps an
# item's communality below 1.  A communality above 1 (a Heywood case, negative
# residual variance) gives a warning; the communalities are returned.

test_that(".dwls_communalities() warns only when a communality exceeds 1", {
  L <- rbind(c(0.6, 0.5), c(0.7, 0.3), c(0.9, 0.6))
  expect_warning(h <- .dwls_communalities(L, c("x1", "x2", "x3")), "x3")
  expect_equal(unname(h), rowSums(L^2))
  expect_silent(.dwls_communalities(L[1:2, ], c("x1", "x2")))
})

# G .6 and specific .5 loadings, two specific factors of five items, five
# categories.  Seed 4 gives a communality of 1.026 for a5; seed 1 is proper.
sim_g2 <- function(seed, n = 800L) {
  set.seed(seed)
  g <- stats::rnorm(n); s <- matrix(stats::rnorm(2 * n), n, 2)
  X <- list()
  for (j in 1:2) for (i in 1:5) {
    y <- 0.6 * g + 0.5 * s[, j] + stats::rnorm(n, sd = sqrt(1 - 0.36 - 0.25))
    X[[paste0(letters[j], i)]] <- findInterval(y, c(-0.8, -0.2, 0.4, 1.0) + 0.1 * (i %% 3))
  }
  as.data.frame(X)
}

test_that("besem_ordered() flags a Heywood case and stays silent on a proper solution", {
  skip_on_cran()
  sf <- list(S1 = paste0("a", 1:5), S2 = paste0("b", 1:5))
  expect_warning(f <- suppressMessages(besem_ordered(sim_g2(4L), specific_factors = sf)), "a5")
  expect_gt(f$communalities[["a5"]], 1)
  w <- character(0)
  f1 <- withCallingHandlers(suppressMessages(besem_ordered(sim_g2(1L), specific_factors = sf)),
                            warning = function(e) { w <<- c(w, conditionMessage(e)); invokeRestart("muffleWarning") })
  expect_false(any(grepl("Heywood", w)))
  expect_lt(max(f1$communalities), 1)
})

# Seed 3: from the rotated EFA solution with its echelon cells set to 0 the
# optimizer ran off to a communality in the millions, in a region where the
# DWLS discrepancy is lower than at the proper local minimum that Mplus and
# lavaan report (chi-square 18.79 / 18.77).  The second start, the EFA solution
# rotated into echelon form, reaches that minimum, and it is kept.
test_that("besem_ordered() keeps the proper local solution over a lower improper one", {
  skip_on_cran()
  sf <- list(S1 = paste0("a", 1:5), S2 = paste0("b", 1:5))
  w <- m <- character(0)
  f <- withCallingHandlers(besem_ordered(sim_g2(3L), specific_factors = sf),
                           warning = function(e) { w <<- c(w, conditionMessage(e)); invokeRestart("muffleWarning") },
                           message = function(e) { m <<- c(m, conditionMessage(e)); invokeRestart("muffleMessage") })
  expect_length(w, 0L)
  expect_true(any(grepl("lower DWLS discrepancy", m)))
  expect_lt(max(f$communalities), 1)
  expect_lt(abs(f$wlsmv_stats$chisq - 18.769), 0.01)
})
