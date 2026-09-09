# compute_indices(): composite reliabilities use absolute loadings (Morin, Arens & Marsh
# 2016), so reversing an item or reflecting a factor cannot change them, and residual
# variances of correlated-factor models come from the model-implied variance
# (1 - diag(Lambda Phi Lambda')), not from 1 - rowSums(Lambda^2).
# HolzingerSwineford1939 serves the CFA / ESEM checks; its nine-item bifactor is
# degenerate (Heywood cases), so the B-ESEM checks use psych::bfi (ordered path).
data("HolzingerSwineford1939", package = "lavaan")
hs <- HolzingerSwineford1939
.indices <- function(spec) {
  res <- suppressWarnings(suppressMessages(run_comparison(spec, run_alignment = FALSE, n_starts = 5L)))
  list(res = res, ind = suppressWarnings(suppressMessages(compute_indices(res))))
}
.hs_spec <- function(d) specify_model(Visual = c("x1", "x2", "x3"), Textual = c("x4", "x5", "x6"),
                                      Speed = c("x7", "x8", "x9"), data = d, label = "HS")
.bfi_spec <- function(d) {
  fl <- list(A = paste0("A", 1:5), C = paste0("C", 1:5), E = paste0("E", 1:5), N = paste0("N", 1:5), O = paste0("O", 1:5))
  do.call(specify_model, c(fl, list(data = d, ordered = TRUE, label = "BFI")))
}
std_L <- function(lav) {
  ss <- lavaan::standardizedSolution(lav, se = FALSE); ss <- ss[ss$op == "=~", ]
  L <- matrix(0, length(unique(ss$rhs)), length(unique(ss$lhs)), dimnames = list(unique(ss$rhs), unique(ss$lhs)))
  for (i in seq_len(nrow(ss))) L[ss$rhs[i], ss$lhs[i]] <- ss$est.std[i]
  L
}

test_that("reversing items leaves the CFA and ESEM composite reliabilities unchanged", {
  skip_on_cran()
  a <- .indices(.hs_spec(hs))
  hs_rev <- hs; hs_rev$x1 <- -hs_rev$x1; hs_rev$x7 <- -hs_rev$x7
  b <- .indices(.hs_spec(hs_rev))
  for (m in c("cfa", "esem")) {
    expect_equal(b$ind[[m]]$omega_total, a$ind[[m]]$omega_total, tolerance = 2e-3, label = paste(m, "omega_total"))
    for (s in c("Visual", "Textual", "Speed"))
      expect_equal(b$ind[[m]]$subscales[[s]]$omega_sub, a$ind[[m]]$subscales[[s]]$omega_sub, tolerance = 2e-3,
                   label = paste(m, s, "omega_sub"))
  }
  expect_equal(b$ind$alpha$G, a$ind$alpha$G, tolerance = 2e-3)
})

test_that("ESEM omega total uses residual variances implied by correlated factors", {
  skip_on_cran()
  a <- .indices(.hs_spec(hs))
  lav <- a$res$fit_esem$lavaan_fit
  L <- std_L(lav); Phi <- lavaan::lavInspect(lav, "cor.lv")[colnames(L), colnames(L)]
  theta <- 1 - diag(L %*% Phi %*% t(L))
  # factors reflected so that every column sum is positive, then absolute loadings
  D <- diag(sign(colSums(L))); Phi_r <- D %*% Phi %*% D
  cs <- colSums(abs(L)); common <- as.numeric(t(cs) %*% Phi_r %*% cs)
  expect_equal(a$ind$esem$omega_total, round(common / (common + sum(theta)), 3))
})

test_that("B-ESEM omega total and omega-H follow the absolute-loading formula and survive item reversal", {
  skip_on_cran()
  data(bfi, package = "psych")
  a <- .indices(.bfi_spec(bfi))
  fb <- a$res$fit_besem
  L <- fb$std_rotated_loadings; g <- fb$g_name
  theta <- pmax(1 - rowSums(L^2), 1e-6)
  cs <- colSums(abs(L)); total <- sum(cs^2) + sum(theta)
  expect_equal(a$ind$besem$omega_total, round(sum(cs^2) / total, 3))
  expect_equal(a$ind$besem$omega_h_g, round(cs[[g]]^2 / total, 3))
  # reverse-score the reverse-worded items: every composite index must stay the same
  bfi_rev <- bfi
  for (it in c("A1", "C4", "C5", "E1", "E2", "O2", "O5")) bfi_rev[[it]] <- 7 - bfi_rev[[it]]
  b <- .indices(.bfi_spec(bfi_rev))
  expect_equal(b$ind$besem$omega_total, a$ind$besem$omega_total, tolerance = 2e-3)
  expect_equal(b$ind$besem$omega_h_g, a$ind$besem$omega_h_g, tolerance = 2e-3)
  expect_equal(b$ind$besem$ecv, a$ind$besem$ecv, tolerance = 2e-3)
  for (s in c("A", "C", "E", "N", "O"))
    expect_equal(b$ind$besem$subscales[[s]]$omega_sub, a$ind$besem$subscales[[s]]$omega_sub, tolerance = 2e-3, label = s)
})

test_that("subscale omega footnotes say the other factor is treated as error, not partialled out", {
  skip_on_cran()
  data(bfi, package = "psych")
  a <- .indices(.bfi_spec(bfi))
  txt <- paste(capture.output(print(a$ind)), collapse = " ")
  expect_false(grepl("partialled", txt, fixed = TRUE))
  expect_match(txt, "additional source of error", fixed = TRUE)
})

test_that("reliability and alignment printouts follow the console width", {
  skip_on_cran()
  data(bfi, package = "psych")
  a <- .indices(.bfi_spec(bfi))
  al <- alignment_check(na.omit(bfi), list(A = paste0("A", 1:5), C = paste0("C", 1:5), E = paste0("E", 1:5),
                                          N = paste0("N", 1:5), O = paste0("O", 1:5)))
  out <- withr::with_options(list(width = 80), c(capture.output(print(a$ind)), capture.output(print(al))))
  expect_lte(max(nchar(out)), 84)   # fixed-format tables are up to 82 wide
})
