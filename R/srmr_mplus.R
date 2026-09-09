# Internal: SRMR as Mplus computes it (Asparouhov & Muthen, 2018, "SRMR in
# Mplus", statmodel.com/download/SRMR2.pdf).
#
# Continuous items: squared residuals of the correlations (r_jk - rho_jk), of
# the variances ((s_jj - sigma_jj) / s_jj) and of the standardized means
# (m_j / sqrt(s_jj) - mu_j / sqrt(sigma_jj)), over p(p+1)/2 + p.  Mplus always
# models the means; when lavaan fits no mean structure the mean residuals are
# 0 and only the denominator counts them.
# Ordered items: squared residuals of the polychoric correlations plus the
# category-probability residuals implied by the model vs the sample thresholds
# (theta parameterization: thresholds standardized by the implied SD), over
# p(p-1)/2 + the total number of categories.
# Several groups: the sums of squares and the denominators are pooled with the
# group sample-size weights before the square root.  This reproduces Mplus 9
# to three decimals at all six invariance levels of the validation runs
# (validation/_validate_through_*.R); the n-weighted mean of the group SRMRs
# (eq. 11 of the note) is off by up to 0.002 there.
.srmr_mplus <- function(lav) {
  G   <- lavaan::lavInspect(lav, "ngroups")
  nob <- lavaan::lavInspect(lav, "nobs")
  ss  <- lavaan::lavInspect(lav, "sampstat")
  im  <- lavaan::lavInspect(lav, "implied")
  if (G == 1L) { ss <- list(ss); im <- list(im) }
  ordered <- length(lavaan::lavNames(lav, "ov.ord")) > 0L

  parts <- vapply(seq_len(G), function(g) {
    S <- ss[[g]]$cov; Sig <- im[[g]]$cov; p <- ncol(S); lt <- lower.tri(S)
    ssq <- sum((stats::cov2cor(S)[lt] - stats::cov2cor(Sig)[lt])^2)
    if (ordered) {
      th_s <- ss[[g]]$th
      th_i <- im[[g]]$th
      sd_i <- sqrt(diag(Sig)); names(sd_i) <- colnames(Sig)
      item <- sub("\\|.*$", "", names(th_i))
      th_i <- th_i / sd_i[item]
      n_cat <- 0L
      for (v in unique(item)) {
        pi <- diff(c(0, stats::pnorm(th_i[item == v]), 1))
        qi <- diff(c(0, stats::pnorm(th_s[item == v]), 1))
        ssq   <- ssq + sum((pi - qi)^2)
        n_cat <- n_cat + length(pi)
      }
      return(c(ssq, p * (p - 1) / 2 + n_cat))
    }
    s <- diag(S); v <- diag(Sig)
    ssq <- ssq + sum(((s - v) / s)^2)
    if (!is.null(ss[[g]]$mean))
      ssq <- ssq + sum((ss[[g]]$mean / sqrt(s) - im[[g]]$mean / sqrt(v))^2)
    c(ssq, p * (p + 1) / 2 + p)
  }, numeric(2L))

  w <- nob / sum(nob)
  sqrt(sum(w * parts[1L, ] / parts[2L, ]))
}
