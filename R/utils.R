# -- Internal utilities --------------------------------------------------------

#' Null-coalescing operator
#' @noRd
`%||%` <- function(a, b) if (!is.null(a)) a else b


# Pre-flight check used by ordered = TRUE paths (esem_ordered, besem_ordered).
# psych::polychoric only supports up to 8 distinct response categories per item;
# raise an informative error pointing the user to ordered = FALSE for items
# above that bound, before any heavy computation kicks off.
.assert_polychoric_compatible <- function(data, indicators) {
  ncats <- vapply(indicators, function(it)
    length(unique(stats::na.omit(data[[it]]))), integer(1L))
  bad <- which(ncats > 8L)
  if (length(bad)) {
    stop(
      "ordered = TRUE supports up to 8 response categories. ",
      "Item(s) with too many categories: ",
      paste(sprintf("%s (%d)", indicators[bad], ncats[bad]), collapse = ", "),
      ".\n  Use ordered = FALSE (MLR + Pearson correlations) for items with ",
      "more than 8 categories.",
      call. = FALSE
    )
  }
  invisible(NULL)
}


# FIML factor analysis: handles missing data case-by-case.
# Retained for future use if lavaan's FIML EFA convergence improves.
# Currently unused: continuous B-ESEM uses listwise + lavaan EFA instead.
# @keywords internal
.fiml_efa <- function(X, k, n_starts = 50L, max_iter = 5000L) {
  # X: N x p data matrix (may contain NAs)
  # k: number of factors
  # Returns: loadings, uniquenesses, mu, chi-squared, df, etc.

  N <- nrow(X); p <- ncol(X)

  # -- Pre-compute per-pattern sufficient statistics --------------------------
  # Group cases by missing-data pattern for efficiency
  obs_pattern <- !is.na(X)
  pat_key     <- apply(obs_pattern, 1, paste, collapse = "")
  pat_groups  <- split(seq_len(N), pat_key)

  patterns <- lapply(pat_groups, function(rows) {
    obs   <- obs_pattern[rows[1], ]
    idx   <- which(obs)
    p_i   <- length(idx)
    x_mat <- X[rows, idx, drop = FALSE]
    list(idx = idx, p_i = p_i, n_i = length(rows), x_mat = x_mat)
  })

  # -- Echelon loading index map ----------------------------------------------
  loading_idx <- matrix(0L, nrow = p, ncol = k)
  counter <- 0L
  for (i in seq_len(p))
    for (j in seq_len(k))
      if (i > k || j <= i) { counter <- counter + 1L; loading_idx[i, j] <- counter }
  n_load <- counter
  # Parameters: loading params + log-uniquenesses + means
  n_theta <- n_load + p + p

  unpack <- function(theta) {
    L <- matrix(0, p, k)
    lp <- theta[seq_len(n_load)]
    for (i in seq_len(p)) for (j in seq_len(k)) {
      idx <- loading_idx[i, j]; if (idx > 0L) L[i, j] <- lp[idx]
    }
    list(L   = L,
         psi = exp(theta[(n_load + 1L):(n_load + p)]),
         mu  = theta[(n_load + p + 1L):n_theta])
  }
  pack <- function(L, log_psi, mu) {
    lp <- numeric(n_load)
    for (i in seq_len(p)) for (j in seq_len(k)) {
      idx <- loading_idx[i, j]; if (idx > 0L) lp[idx] <- L[i, j]
    }
    c(lp, log_psi, mu)
  }

  # -- FIML negative log-likelihood (sum over patterns) -----------------------
  fn <- function(theta) {
    u     <- unpack(theta)
    Sigma <- tcrossprod(u$L) + diag(u$psi, nrow = p)
    total <- 0
    for (pat in patterns) {
      idx  <- pat$idx; p_i <- pat$p_i; n_i <- pat$n_i
      S_i  <- Sigma[idx, idx, drop = FALSE]
      R_i  <- tryCatch(chol(S_i), error = function(e) NULL)
      if (is.null(R_i)) return(1e10)
      Si_i <- chol2inv(R_i)
      ld_i <- 2 * sum(log(diag(R_i)))
      mu_i <- u$mu[idx]
      # Sum over cases in this pattern
      for (r in seq_len(n_i)) {
        d <- pat$x_mat[r, ] - mu_i
        total <- total + ld_i + sum(d * (Si_i %*% d))
      }
    }
    total / (2 * N)  # average per case (Mplus convention for optimisation)
  }

  # -- Analytical gradient ----------------------------------------------------
  gr <- function(theta) {
    u     <- unpack(theta)
    Sigma <- tcrossprod(u$L) + diag(u$psi, nrow = p)

    # Accumulate gradient contributions across patterns
    dL_full   <- matrix(0, p, k)
    dpsi_full <- numeric(p)
    dmu_full  <- numeric(p)

    for (pat in patterns) {
      idx <- pat$idx; p_i <- pat$p_i; n_i <- pat$n_i
      S_i <- Sigma[idx, idx, drop = FALSE]
      R_i <- tryCatch(chol(S_i), error = function(e) NULL)
      if (is.null(R_i)) return(rep(0, n_theta))
      Si_i <- chol2inv(R_i)
      mu_i <- u$mu[idx]

      # Accumulate outer products and mean residuals for this pattern
      sum_dd <- matrix(0, p_i, p_i)
      sum_d  <- numeric(p_i)
      for (r in seq_len(n_i)) {
        d <- pat$x_mat[r, ] - mu_i
        sum_dd <- sum_dd + tcrossprod(d)
        sum_d  <- sum_d + d
      }
      # D_i = Si_i - Si_i %*% (sum_dd / n_i) %*% Si_i  (per-case average)
      # But for total gradient, use full sums
      D_i <- n_i * Si_i - Si_i %*% sum_dd %*% Si_i

      # Gradient w.r.t. Lambda (submatrix)
      dL_full[idx, ] <- dL_full[idx, ] + D_i %*% u$L[idx, , drop = FALSE]

      # Gradient w.r.t. psi (diagonal of D_i)
      dpsi_full[idx] <- dpsi_full[idx] + diag(D_i)

      # Gradient w.r.t. mu
      dmu_full[idx] <- dmu_full[idx] - Si_i %*% sum_d
    }

    # Pack gradient: dF/dL, dF/d(log_psi), dF/dmu
    gl <- numeric(n_load)
    for (i in seq_len(p)) for (j in seq_len(k)) {
      idx <- loading_idx[i, j]
      if (idx > 0L) gl[idx] <- dL_full[i, j]
    }
    c(gl / N, u$psi * dpsi_full / N, dmu_full / N)
  }

  # -- EM algorithm for starting values ---------------------------------------
  # Compute EM-estimated covariance and means (handles missing data)
  mu_em  <- colMeans(X, na.rm = TRUE)
  S_em   <- cov(X, use = "pairwise.complete.obs")
  # Replace any NAs in S_em with 0 (rare edge case)
  S_em[is.na(S_em)] <- 0

  # Use factanal on EM covariance for initial loading estimates
  fa_init <- tryCatch(
    factanal(covmat = S_em, factors = k, n.obs = N, rotation = "none",
             control = list(nstart = 10L)),
    error = function(e) NULL
  )

  # -- Multi-start optimisation -----------------------------------------------
  best_val <- Inf; best_par <- NULL
  # Seed scoped to this block via withr (RNG state restored afterwards) so
  # multi-start is reproducible without disturbing the caller's RNG stream.
  withr::with_seed(42L, for (s in seq_len(n_starts)) {
    if (s == 1L && !is.null(fa_init)) {
      # First start: use factanal/EM solution (warm start)
      L0   <- unclass(fa_init$loadings)
      psi0 <- fa_init$uniquenesses
      mu0  <- mu_em
    } else {
      # Random starts
      L0 <- matrix(rnorm(p * k, 0, 0.3), p, k)
      psi0 <- runif(p, 0.3, 0.9)
      mu0  <- mu_em + rnorm(p, 0, 0.1)
    }
    # Apply echelon constraints
    for (i in seq_len(min(k, p))) if (i < k) L0[i, (i + 1L):k] <- 0

    theta0 <- pack(L0, log(pmax(psi0, 0.005)), mu0)
    res <- tryCatch(nlminb(theta0, fn, gr,
      control = list(iter.max = max_iter, eval.max = max_iter * 5L,
                     rel.tol = 1e-10, x.tol = 1e-10)),
      error = function(e) NULL)
    if (!is.null(res) && is.finite(res$objective) && res$objective < best_val) {
      best_val <- res$objective; best_par <- res$par
    }
  })
  if (is.null(best_par)) stop(".fiml_efa: all starts failed", call. = FALSE)

  # -- Extract solution -------------------------------------------------------
  final <- unpack(best_par)
  L   <- final$L; psi <- final$psi; mu <- final$mu

  # FIML chi-squared: 2*(LL_H1 - LL_H0) where both use the FIML formula.
  # H1 (saturated): unrestricted mean + covariance (EM-estimated).
  # H0 (model): Lambda*Lambda' + Psi with estimated mu.
  # Estimate H1 via EM algorithm on the raw data with missing values.
  mu_h1 <- colMeans(X, na.rm = TRUE)
  S_h1  <- cov(X, use = "pairwise.complete.obs")
  S_h1[is.na(S_h1)] <- 0
  # Simple EM to refine H1 estimates
  for (em_iter in seq_len(200L)) {
    T_sum  <- matrix(0, p, p); mu_sum <- numeric(p); n_eff <- 0
    for (pat in patterns) {
      idx <- pat$idx; p_i <- pat$p_i; n_i <- pat$n_i
      mis <- setdiff(seq_len(p), idx)
      S_oo <- S_h1[idx, idx, drop = FALSE]
      R_oo <- tryCatch(chol(S_oo), error = function(e) NULL)
      if (is.null(R_oo)) next
      Si_oo <- chol2inv(R_oo)
      if (length(mis) > 0) {
        S_mo <- S_h1[mis, idx, drop = FALSE]
        B    <- S_mo %*% Si_oo
        S_mm_cond <- S_h1[mis, mis, drop = FALSE] - B %*% S_h1[idx, mis, drop = FALSE]
      }
      for (r in seq_len(n_i)) {
        x_o <- pat$x_mat[r, ]
        x_imp <- mu_h1
        x_imp[idx] <- x_o
        if (length(mis) > 0) {
          x_imp[mis] <- mu_h1[mis] + B %*% (x_o - mu_h1[idx])
          T_sum[mis, mis] <- T_sum[mis, mis] + S_mm_cond
        }
        mu_sum <- mu_sum + x_imp
        T_sum  <- T_sum + tcrossprod(x_imp)
        n_eff  <- n_eff + 1L
      }
    }
    mu_new <- mu_sum / n_eff
    S_new  <- T_sum / n_eff - tcrossprod(mu_new)
    if (max(abs(mu_new - mu_h1)) < 1e-8 && max(abs(S_new - S_h1)) < 1e-8) break
    mu_h1 <- mu_new; S_h1 <- S_new
  }

  # Compute FIML LL for H1 (saturated) using EM estimates
  ll_h1 <- 0
  for (pat in patterns) {
    idx <- pat$idx; p_i <- pat$p_i; n_i <- pat$n_i
    S_i <- S_h1[idx, idx, drop = FALSE]
    R_i <- tryCatch(chol(S_i), error = function(e) NULL)
    if (is.null(R_i)) next
    Si_i <- chol2inv(R_i)
    ld_i <- 2 * sum(log(diag(R_i)))
    mu_i <- mu_h1[idx]
    for (r in seq_len(n_i)) {
      d <- pat$x_mat[r, ] - mu_i
      ll_h1 <- ll_h1 + ld_i + sum(d * (Si_i %*% d))
    }
  }
  ll_h1 <- ll_h1 / (2 * N)   # same scale as best_val

  chi  <- 2 * N * (best_val - ll_h1)  # should be positive (model fits worse than saturated)
  df   <- ((p - k)^2 - p - k) %/% 2

  # Model-implied covariance (for fit indices)
  Sigma <- tcrossprod(L) + diag(psi, nrow = p)

  list(loadings     = L,
       uniquenesses = psi,
       mu           = mu,
       Sigma        = Sigma,
       STATISTIC    = chi,
       dof          = df,
       PVAL         = if (df > 0) pchisq(chi, df, lower.tail = FALSE) else NA_real_,
       f_min        = best_val,
       f_h1         = ll_h1,
       N            = N)
}


# Pre-compute standardized loadings (STDYX) and delta-method SEs from a
# single-group lavaan rotation fit and return as named matrices.
# Called at the end of esem() and besem() to avoid expensive on-demand
# recomputation in parameters() / std_loadings().
# Returns NULL silently on any failure (caller falls back to standardizedsolution).
.cache_std_loadings <- function(fit) {
  tryCatch({
    ss    <- lavaan::standardizedsolution(fit)
    ss_lv <- ss[ss$op == "=~", , drop = FALSE]
    if (nrow(ss_lv) == 0L) return(NULL)
    items   <- unique(ss_lv$rhs)
    factors <- unique(ss_lv$lhs)
    L  <- matrix(NA_real_, nrow = length(items), ncol = length(factors),
                 dimnames = list(items, factors))
    SE <- matrix(NA_real_, nrow = length(items), ncol = length(factors),
                 dimnames = list(items, factors))
    for (i in seq_len(nrow(ss_lv))) {
      ri <- ss_lv$rhs[i]; fi <- ss_lv$lhs[i]
      if (ri %in% items && fi %in% factors) {
        L[ri, fi]  <- ss_lv$est.std[i]
        SE[ri, fi] <- ss_lv$se[i]
      }
    }
    list(L = L, SE = SE)
  }, error = function(e) NULL)
}


# Extract factor correlation matrix from lavaan standardized solution.
# Returns a symmetric nf x nf matrix with 1s on the diagonal.
.extract_phi_from_ss <- function(ss, factor_names) {
  nf <- length(factor_names)
  Phi <- diag(nf)
  dimnames(Phi) <- list(factor_names, factor_names)
  ss_phi <- ss[ss$op == "~~" & ss$lhs %in% factor_names &
                 ss$rhs %in% factor_names & ss$lhs != ss$rhs, ]
  for (r in seq_len(nrow(ss_phi))) {
    i <- match(ss_phi$lhs[r], factor_names)
    j <- match(ss_phi$rhs[r], factor_names)
    if (!is.na(i) && !is.na(j))
      Phi[i, j] <- Phi[j, i] <- ss_phi$est.std[r]
  }
  Phi
}


#' Parse Mplus Polychoric Correlation Matrix from Output File
#'
#' Reads a Mplus \code{.out} file and assembles the polychoric (tetrachoric)
#' correlation matrix that Mplus prints under \code{SAMPLE STATISTICS} when
#' \code{SAMPSTAT} is requested in the \code{OUTPUT} section.
#'
#' The matrix can be passed directly to \code{\link{besem_ordered}} via the
#' \code{r_obs_override} argument so that R's DWLS optimisation uses the
#' same polychoric estimates as Mplus, eliminating rotational discrepancies
#' caused by differences in the polychoric estimators.
#'
#' @param out_file Character. Full path to a Mplus \code{.out} file that
#'   contains a \code{CORRELATION MATRIX (WITH VARIANCES ON THE DIAGONAL)}
#'   section. This section is written when \code{SAMPSTAT} appears in the
#'   \code{OUTPUT} block of the Mplus input file.
#'
#' @return A named, symmetric numeric matrix of polychoric correlations with
#'   1s on the diagonal.  Row and column names are the item names as reported
#'   by Mplus (typically upper-case).
#'
#' @examples
#' \dontrun{
#' # Reads a polychoric matrix from a Mplus .out file (SAMPSTAT output),
#' # so it needs a .out produced by a licensed Mplus run.
#' R_mplus <- parse_mplus_polychoric(
#'   "path/to/besem_measurement.out"
#' )
#' dim(R_mplus)   # 18 x 18
#'
#' # Use as input to besem_ordered so R rotates from the same matrix as Mplus
#' fit_b <- besem_ordered(
#'   data             = mydata,
#'   specific_factors = factor_items,
#'   r_obs_override   = R_mplus
#' )
#' }
#'
#' @export
parse_mplus_polychoric <- function(out_file) {

  if (!file.exists(out_file))
    stop("File not found: ", out_file, call. = FALSE)

  lines  <- readLines(out_file, warn = FALSE)
  marker <- "CORRELATION MATRIX (WITH VARIANCES ON THE DIAGONAL)"
  block_starts <- grep(marker, lines, fixed = TRUE)

  if (length(block_starts) == 0)
    stop(
      "No '", marker, "' section found in:\n  ", out_file, "\n",
      "  Add SAMPSTAT to the OUTPUT section of your Mplus input file.",
      call. = FALSE
    )

  # Collect entries from all blocks: list(row, cols, vals)
  entries        <- list()
  all_item_names <- character(0)   # preserves encounter order

  for (bs in block_starts) {

    i <- bs + 1L

    # Skip blank lines to the column-header line
    while (i <= length(lines) && trimws(lines[i]) == "") i <- i + 1L
    if (i > length(lines)) next

    col_names <- strsplit(trimws(lines[i]), "\\s+")[[1]]
    i <- i + 1L

    # Skip underscore separator line(s)
    while (i <= length(lines) && grepl("^\\s*_+", lines[i])) i <- i + 1L

    # Read data rows until a blank line or a non-item line
    while (i <= length(lines)) {
      ln    <- lines[i]
      ln_tr <- trimws(ln)
      if (ln_tr == "") break                     # blank line ends block

      parts <- strsplit(ln_tr, "\\s+")[[1]]
      if (length(parts) == 0) break

      # First token must look like an item name (starts with a letter)
      if (!grepl("^[A-Za-z]", parts[1L])) break

      row_name <- parts[1L]
      vals     <- suppressWarnings(as.numeric(parts[-1L]))
      vals     <- vals[!is.na(vals)]             # drop any stray non-numerics

      # Accumulate ordered item names
      all_item_names <- union(all_item_names,
                              c(col_names, row_name))

      if (length(vals) > 0L)
        entries[[length(entries) + 1L]] <- list(
          row  = row_name,
          cols = col_names[seq_along(vals)],
          vals = vals
        )

      i <- i + 1L
    }
  }

  if (length(all_item_names) == 0L)
    stop("Could not parse any item names from the correlation matrix blocks.",
         call. = FALSE)

  # Build symmetric matrix with 1s on diagonal
  n   <- length(all_item_names)
  mat <- diag(n)
  dimnames(mat) <- list(all_item_names, all_item_names)

  for (e in entries) {
    ri <- match(e$row, all_item_names)
    for (ci_seq in seq_along(e$vals)) {
      cj <- match(e$cols[[ci_seq]], all_item_names)
      v  <- e$vals[[ci_seq]]
      if (!is.na(ri) && !is.na(cj) && is.finite(v)) {
        mat[ri, cj] <- v
        mat[cj, ri] <- v
      }
    }
  }

  mat
}


# -- Heywood case detection and iterative rotation correction ------------------
#
# Called by esem() and esem_ordered() when heywood_fix = TRUE.
# Detects standardised loadings |STDYX| >= 1.0 and resolves them by
# reconstructing unrotated loadings via Cholesky decomposition of Phi,
# then re-rotating with GPArotation::targetQ() using a Mplus-style target
# matrix (NA for free primary cells, 0 for cross-loadings, step-down values
# for Heywood items). Each round is a rotation call (~1ms), not a full
# lavaan::cfa() re-estimation.
#
# Arguments:
#   L_rot             -- items x factors matrix of STDYX loadings from lavaan
#   Phi_rot           -- factors x factors correlation matrix from lavaan
#   initial_target    -- target matrix used in the initial fit, or NULL
#   indicators        -- character vector of item names
#   factor_names      -- character vector of factor names
#   original_rotation -- rotation string used for the initial fit
#   max_rounds        -- integer, default 15
#
# Returns list(L_stdyx = <corrected loading matrix>,
#              Phi_cor = <corrected factor correlations>,
#              log     = <heywood_log list or NULL>)
.heywood_retry_loop <- function(L_rot, Phi_rot, initial_target,
                                 indicators, factor_names,
                                 original_rotation = "geomin",
                                 max_rounds = 15L) {

  max_rounds <- as.integer(max_rounds)
  n_items    <- length(indicators)
  n_factors  <- length(factor_names)

  # Single-factor models have no rotation freedom -- nothing to fix
  if (n_factors < 2L) return(list(L_stdyx = L_rot, Phi_cor = Phi_rot, log = NULL))

  # -- 1. L_rot from lavaan's standardizedSolution() IS already STDYX ---------
  # Phi_rot is a correlation matrix (diagonal = 1), so this sweep is a no-op
  # in practice. Kept for correctness if Phi ever has non-unit diagonal.
  # Ensure dimnames are set (lavInspect may strip them).
  if (is.null(rownames(Phi_rot))) dimnames(Phi_rot) <- list(factor_names, factor_names)
  fsd     <- sqrt(diag(Phi_rot))
  L_stdyx <- sweep(L_rot, 2, fsd, "*")

  # -- 2. Detect Heywood loadings ---------------------------------------------
  hw_mask <- abs(L_stdyx) >= 1.0
  if (!any(hw_mask)) return(list(L_stdyx = L_stdyx, Phi_cor = Phi_rot, log = NULL))

  heywood_pairs <- which(hw_mask, arr.ind = TRUE)
  heywood_items   <- unique(indicators[heywood_pairs[, 1L]])
  heywood_factors <- factor_names[heywood_pairs[, 2L]]
  heywood_lambdas <- L_stdyx[hw_mask]

  # -- 3. Report detection ----------------------------------------------------
  message("\nHeywood case(s) detected -- standardised loading(s) above 1.0:")
  for (k in seq_along(heywood_lambdas))
    message(sprintf("  %-6s-> %-12s lam = %.3f",
                heywood_factors[k],
                indicators[heywood_pairs[k, 1L]],
                heywood_lambdas[k]))

  detected_df <- data.frame(
    factor = heywood_factors,
    item   = indicators[heywood_pairs[, 1L]],
    lambda = heywood_lambdas,
    stringsAsFactors = FALSE
  )

  # -- 4. Reconstruct unrotated loadings via Cholesky -------------------------
  # L_rot = L_unrot  x  T^{-1}, Phi_rot = T  x  T'
  # Therefore: L_unrot = L_rot  x  t(chol(Phi_rot))
  message("\nRe-rotating via GPArotation (bypassing lavaan rotation)...")
  Phi_chol <- t(chol(Phi_rot))
  L_unrot  <- L_rot %*% Phi_chol

  # -- 5. Determine primary factor per item -----------------------------------
  if (!is.null(initial_target)) {
    missing_rows <- setdiff(indicators, rownames(initial_target))
    if (length(missing_rows))
      stop(".heywood_retry_loop: initial_target is missing rows for: ",
           paste(missing_rows, collapse = ", "), call. = FALSE)
    pf_idx <- apply(initial_target[indicators, factor_names, drop = FALSE],
                    1, function(row) {
                      w <- which(row == 1 | is.na(row))
                      if (length(w)) w[1L] else NA_integer_
                    })
  } else {
    pf_idx <- apply(abs(L_stdyx), 1, function(row) {
      if (all(is.na(row))) NA_integer_ else which.max(row)
    })
  }
  names(pf_idx) <- indicators

  # -- 6. Build Mplus-style target matrix -------------------------------------
  # NA = free (not targeted), 0 = cross-loading target
  tmat <- matrix(0, n_items, n_factors, dimnames = list(indicators, factor_names))
  for (item in indicators) {
    pf <- pf_idx[item]
    if (!is.na(pf)) tmat[item, pf] <- NA   # primary cells are FREE
  }

  # -- 7. Initialise per-item targets and tracking ----------------------------
  # Preserve sign: a negative Heywood loading should target -0.95, not +0.95
  cur_target <- setNames(vapply(heywood_items, function(item) {
    pf <- pf_idx[item]
    sgn <- if (!is.na(pf)) sign(L_stdyx[item, pf]) else 1
    sgn * 0.95
  }, numeric(1L)), heywood_items)
  still_above       <- heywood_items
  resolved_at_round <- setNames(rep(NA_integer_, length(heywood_items)), heywood_items)
  rounds_done       <- 0L

  best_L_stdyx <- L_stdyx
  best_Phi_cor <- Phi_rot
  best_hw_max  <- max(abs(heywood_lambdas))

  # -- 8. Retry loop (rotation only -- no re-estimation) ----------------------
  for (round in seq_len(max_rounds)) {
    rounds_done <- round

    for (item in still_above) {
      pf <- pf_idx[item]
      if (!is.na(pf)) tmat[item, pf] <- cur_target[item]
    }

    tgt_str <- paste(
      vapply(still_above, function(it) sprintf("%s = %.2f", it, cur_target[it]),
             character(1L)),
      collapse = ", ")
    message(sprintf("  Targets -> %s", tgt_str))

    # Use targetT for orthogonal, targetQ for oblique
    # Note: "T$" regex must not match "target" (oblique) -- only explicit "targetT"
    is_orth <- tolower(original_rotation) %in% c("targett", "varimax", "geomint")
    rot_fn  <- if (is_orth) GPArotation::targetT else GPArotation::targetQ
    rot <- tryCatch(
      rot_fn(L_unrot, Target = tmat),
      error = function(e) {
        message(paste("  GPArotation failed:", conditionMessage(e)))
        NULL
      })
    if (is.null(rot)) break

    L_new   <- rot$loadings
    Phi_new <- rot$Phi
    if (is.null(Phi_new)) Phi_new <- diag(n_factors)

    fsd_new     <- sqrt(diag(Phi_new))
    L_stdyx_new <- sweep(L_new, 2, fsd_new, "*")
    Phi_cor_new <- Phi_new / (fsd_new %o% fsd_new)
    dimnames(L_stdyx_new) <- list(indicators, factor_names)
    dimnames(Phi_cor_new) <- list(factor_names, factor_names)

    new_hw_max <- max(vapply(still_above, function(item) {
      pf <- pf_idx[item]
      if (!is.na(pf)) abs(L_stdyx_new[item, pf]) else 0
    }, numeric(1L)))

    if (new_hw_max < best_hw_max) {
      best_L_stdyx <- L_stdyx_new
      best_Phi_cor <- Phi_cor_new
      best_hw_max  <- new_hw_max
    }

    floor_hit        <- FALSE
    still_above_next <- character(0L)

    for (item in still_above) {
      pf      <- pf_idx[item]
      pf_name <- factor_names[pf]
      lam     <- if (!is.na(pf)) abs(L_stdyx_new[item, pf]) else NA_real_

      # Acceptance threshold 0.994: ensures lam rounds to <1.00 at 2-dp display
      # (APA) while leaving a buffer below the 0.9995 3-dp boundary.
      if (is.na(lam) || lam >= 0.994) {
        message(sprintf("  %-6s-> %-12s lam = %.3f  still >= 0.994",
                    pf_name, item, lam))
        still_above_next <- c(still_above_next, item)
        # Step toward zero, preserving sign
        sgn   <- sign(cur_target[item])
        new_t <- cur_target[item] - sgn * 0.05
        if (abs(new_t) < 0.10) {
          cur_target[item] <- sgn * 0.10
          floor_hit <- TRUE
        } else {
          cur_target[item] <- new_t
        }
      } else {
        message(sprintf("  %-6s-> %-12s lam = %.3f  resolved",
                    pf_name, item, lam))
        resolved_at_round[item] <- round
      }
    }

    still_above <- still_above_next

    if (length(still_above) == 0L) {
      message("\nAll Heywood case(s) resolved. Corrected rotation accepted.")
      break
    }
    if (floor_hit) {
      message("Stopping: target floor (0.10) reached. Returning best attempt.")
      break
    }
  }

  # -- 9. Build log -----------------------------------------------------------
  fin_rows <- do.call(rbind, lapply(unique(detected_df$item), function(item) {
    pf <- pf_idx[item]
    pf_name <- if (!is.na(pf)) factor_names[pf] else NA_character_
    data.frame(
      factor         = pf_name,
      item           = item,
      lambda         = if (!is.na(pf)) best_L_stdyx[item, pf] else NA_real_,
      resolved_round = resolved_at_round[item],
      stringsAsFactors = FALSE)
  }))

  log <- list(
    detected          = detected_df,
    rounds            = rounds_done,
    original_rotation = original_rotation,
    final_rotation    = "target",
    final_targets     = cur_target[unique(detected_df$item)],
    resolved          = length(still_above) == 0L,
    unresolved_items  = if (length(still_above) > 0L) still_above else NULL,
    delta_cfi         = 0,
    final_loadings    = fin_rows,
    method            = "rotation"
  )

  list(L_stdyx = best_L_stdyx, Phi_cor = best_Phi_cor, log = log)
}
