#' ESEM for Ordered-Categorical (Likert) Data
#'
#' Fits ESEM on ordered-categorical indicators (WLSMV / Theta parameterization).
#' Called automatically from \code{\link{esem}} when \code{ordered} is set.
#'
#' **Default** (\code{method = "lavaan"}): lavaan's \code{efa()} block with WLSMV
#' (unrestricted model, post-hoc rotation, delta-method SEs)---aligned with Mplus
#' \code{ESTIMATOR = WLSMV; ROTATION = TARGET (oblique)}.
#'
#' **Alternate** (\code{method = "rotation"}): custom DWLS pipeline (polychoric
#' correlations via \pkg{psych}, \pkg{GPArotation}, sandwich + rotation Jacobian SEs).
#'
#' @inheritSection doc_estimator_paths Estimator paths
#' @inheritSection doc_estimator_paths Missing-data defaults
#' @inheritSection doc_estimator_paths The lavaan_fit slot on custom WLSMV fits
#' @inheritSection doc_estimator_paths Heywood fix and loadings
#'
#' @param data A \code{data.frame} of observed ordered indicators.
#' @param nfactors Integer. Number of latent factors.
#' @param indicators Character vector of item names. If \code{NULL}, all
#'   columns of \code{data} except \code{group} are used.
#' @param rotation Character. Rotation method. Default \code{"target"}.
#'   Supports \code{"target"} (oblique), \code{"targetT"} (orthogonal),
#'   \code{"geomin"}, \code{"geominT"} (orthogonal geomin),
#'   \code{"oblimin"}, \code{"varimax"}.
#' @param target A target matrix from \code{\link{make_target}}. Required
#'   when \code{rotation} contains \code{"target"}.
#'   Convention: \code{NA} = free (no penalty), \code{0} = targeted toward zero.
#' @param factor_names Optional character vector of factor names.
#' @param method Character. \code{"lavaan"} (default): lavaan \code{efa()} + WLSMV.
#'   \code{"rotation"}: custom DWLS + post-hoc target rotation (see
#'   \link[=doc_estimator_paths]{The lavaan_fit slot on custom WLSMV fits}).
#' @param n_starts Integer. Number of random rotation starts when
#'   \code{method = "rotation"}. Default \code{100L}.
#' @param r_obs_override Optional observed correlation matrix to use instead
#'   of polychoric estimation (for reproducibility / testing).
#' @param std.lv Logical. Fix factor variances to 1. Default \code{TRUE}.
#' @param missing Character. Missing data handling passed to lavaan.
#'   Default \code{"pairwise"} (Mplus-style polychoric pairs). Use
#'   \code{"listwise"} for complete cases. When comparing with
#'   \code{\link{besem_ordered}} on the same data, set \code{missing} explicitly
#'   on both calls (\code{besem_ordered} defaults to \code{"listwise"}).
#' @param group Character. Grouping variable name for multi-group models.
#' @param group_equal Character vector of lavaan equality constraints.
#' @param heywood_fix Logical. Retry rotation with Cholesky unrotation when a
#'   standardised loading exceeds 1. Default \code{TRUE}.
#' @param n_obs Ignored (kept for backward compatibility).
#' @param ... Additional arguments passed to \code{lavaan::cfa()}.
#'
#' @return An object of class \code{c("esem_fit_ordered", "esem_fit")}. For
#' \code{method = "lavaan"}, \code{lavaan_fit} is the fitted ESEM model. For
#' \code{method = "rotation"}, \code{lavaan_fit} is an auxiliary one-factor CFA;
#' use \code{std_rotated_loadings}, \code{wlsmv_stats}, and
#' \code{\link{std_loadings}}. Slot \code{estimator} is \code{"WLSMV"} or
#' \code{"DWLS"} depending on path.
#'
#' @details
#' ## Mplus equivalence
#'
#' \preformatted{
#' Mplus:
#'   ANALYSIS:
#'     ESTIMATOR = WLSMV;
#'     ROTATION  = TARGET;
#'     PARAMETERIZATION = THETA;
#'   MODEL:
#'     EX MD CI BY item1-item18 (*1);
#'
#' R (this function):
#'   esem_ordered(data, nfactors = 3, rotation = "target", target = tgt,
#'                factor_names = c("EX","MD","CI"))
#' }
#'
#' ## Algorithm (\code{method = "lavaan"})
#' lavaan's \code{efa()} block fits an unrestricted k-factor model under WLSMV,
#' applies post-hoc rotation, and propagates standard errors through the rotation
#' transformation (Asparouhov & Muthen, 2009).
#'
#' @references
#' Asparouhov, T., & Muthen, B. (2009). Exploratory structural equation
#' modeling. \emph{Structural Equation Modeling, 16}(3), 397--438.
#'
#' @seealso \code{\link{esem}} for continuous data,
#'   \code{\link{besem_ordered}} for bifactor ordered ESEM.
#'
#' @examples
#' \dontrun{
#' tgt <- make_target(
#'   list(EX = items_EX, MD = items_MD, CI = items_CI),
#'   item_names = all_items
#' )
#'
#' fit_ord <- esem_ordered(
#'   data         = Rdata,
#'   nfactors     = 3,
#'   indicators   = all_items,
#'   rotation     = "target",
#'   target       = tgt,
#'   factor_names = c("EX","MD","CI")
#' )
#'
#' summary(fit_ord, fit.measures = TRUE, standardized = TRUE)
#' std_loadings(fit_ord)
#' factor_correlations(fit_ord)
#' }
#'
#' @importFrom lavaan cfa
#' @export
esem_ordered <- function(data,
                          nfactors,
                          indicators    = NULL,
                          rotation      = "target",
                          target        = NULL,
                          factor_names  = NULL,
                          method        = c("lavaan", "rotation"),
                          n_starts      = 100L,
                          r_obs_override = NULL,
                          std.lv        = TRUE,
                          missing       = "pairwise",
                          group         = NULL,
                          group_equal   = NULL,
                          heywood_fix   = TRUE,
                          n_obs         = NULL,
                          ...) {

  mc     <- match.call()
  method <- match.arg(method)

  # -- Validate -----------------------------------------------------------------
  nfactors <- as.integer(nfactors)

  all_cols   <- colnames(data)
  indicators <- if (is.null(indicators)) setdiff(all_cols, group) else indicators
  bad        <- setdiff(indicators, all_cols)
  if (length(bad))
    stop("Indicators not found in data: ", paste(bad, collapse = ", "), call. = FALSE)

  if (is.null(r_obs_override))
    .assert_polychoric_compatible(data, indicators)

  if (is.null(factor_names))
    factor_names <- paste0("F", seq_len(nfactors))
  if (length(factor_names) != nfactors)
    stop("`factor_names` length must equal `nfactors`.", call. = FALSE)
  factor_names <- make.names(factor_names)

  if (grepl("target", rotation, ignore.case = TRUE) && is.null(target))
    stop("Provide a `target` matrix when rotation = \"target\".", call. = FALSE)

  n_i <- length(indicators)
  k   <- nfactors
  n   <- nrow(data)

  # =======================================--====================================
  # METHOD: "rotation"
  # Custom pipeline: psych polychorics + DWLS extraction + oblique targetQ
  # rotation.  Bypasses lavaan's polychoric floor bug (sqrt(.Machine$double.eps)
  # clamping in lav_bvord_noexo_pi_cache) which biases rho upward by up to
  # 0.024 on highly correlated item pairs with near-zero probability cells.
  # Matches lavaan's rotation settings: normalize=FALSE (row.weights="none"),
  # 100 random starts, GPA algorithm.
  # Validated: loadings max|d|=0.001 vs Mplus, SEs max|d|=0.002.
  # ============--============================---==================================
  if (method == "rotation") {

    if (!requireNamespace("psych",       quietly = TRUE))
      stop("Package 'psych' is required for method='rotation'. ",
           "Install: install.packages('psych')", call. = FALSE)
    if (!requireNamespace("GPArotation", quietly = TRUE))
      stop("Package 'GPArotation' is required for method='rotation'. ",
           "Install: install.packages('GPArotation')", call. = FALSE)

    # Convert target from make_target format (1=primary, 0=cross) to
    # GPArotation format (NA=free, 0=target-zero) if needed.
    # make_target() produces 1/0; manual targets may already use NA/0.
    if (!anyNA(target) && any(target == 1)) {
      target[target == 1] <- NA_real_
    }
    cross_mask <- !is.na(target) & target == 0

    # -- Stage 1: Polychoric correlations -------------------------------------
    if (!is.null(r_obs_override)) {
      message("Stage 1: Using supplied polychoric matrix (r_obs_override)...")
      rn  <- rownames(r_obs_override)
      idx <- match(tolower(indicators), tolower(rn))
      if (anyNA(idx))
        stop("r_obs_override item names do not match indicators.\n",
             "  Missing: ", paste(indicators[is.na(idx)], collapse = ", "),
             call. = FALSE)
      R_poly <- r_obs_override[idx, idx]
      rownames(R_poly) <- colnames(R_poly) <- indicators
    } else {
      message("Stage 1: Polychoric correlations via psych::polychoric()...")
      poly_result <- tryCatch(
        psych::polychoric(data[, indicators, drop = FALSE]),
        error = function(e) stop("polychoric() failed: ", conditionMessage(e), call. = FALSE)
      )
      R_poly <- poly_result$rho
    }

    # -- Stage 2: Unrotated EFA starting values --------------------------------
    message("Stage 2: Unrotated EFA (factanal) + GPArotation::targetQ for starting values...")
    fa_unrotated <- tryCatch(
      factanal(covmat = R_poly, factors = k, n.obs = n, rotation = "none"),
      error = function(e) stop("Unrotated EFA failed: ", conditionMessage(e), call. = FALSE)
    )
    A <- unclass(fa_unrotated$loadings)
    # Use normalize=TRUE for starting values (better global optimum search).
    # Only Stage 5's final rotation uses normalize=FALSE to match lavaan.
    efa_rot <- tryCatch({
      suppressWarnings(
        GPArotation::targetQ(A, Target = target, maxit = 10000L, normalize = TRUE)
      )$loadings
    }, error = function(e) stop("Initial EFA rotation failed: ", conditionMessage(e), call. = FALSE))

    # -- Stage 3: DWLS weight matrices via auxiliary CFA ----------------------
    message("Stage 3: Extracting DWLS weight matrices (W, Gamma) from auxiliary lavaan CFA...")
    simple_syntax <- paste0(factor_names[1L], " =~ ", paste(indicators, collapse = " + "))
    fit_aux <- tryCatch(
      lavaan::cfa(simple_syntax, data = data, ordered = indicators,
                  std.lv = TRUE, parameterization = "theta", missing = missing),
      error = function(e) stop("Auxiliary WLSMV CFA failed: ", conditionMessage(e), call. = FALSE)
    )

    wls_obs_raw <- lavaan::lavInspect(fit_aux, "wls.obs")
    W_raw       <- lavaan::lavInspect(fit_aux, "wls.v")
    Gamma_raw   <- lavaan::lavInspect(fit_aux, "gamma")

    # Multi-group: use first group's matrices
    if (is.list(wls_obs_raw) && !is.null(names(wls_obs_raw)) &&
        !any(grepl("~~", names(wls_obs_raw)))) {
      wls_obs_raw <- wls_obs_raw[[1L]]
      W_raw       <- W_raw[[1L]]
      Gamma_raw   <- Gamma_raw[[1L]]
    }

    wls_obs_full <- wls_obs_raw
    W_full       <- W_raw
    Gamma_full   <- as.matrix(Gamma_raw)

    obs_names <- names(wls_obs_full)
    n_wls     <- length(wls_obs_full)
    cor_idx   <- grep("~~", obs_names)

    if (length(cor_idx) == 0)
      stop("Cannot identify polychoric correlation elements in wls.obs.", call. = FALSE)

    W_cor     <- diag(as.matrix(W_full))[cor_idx]
    Gamma_cor <- Gamma_full[cor_idx, cor_idx, drop = FALSE]
    n_pairs   <- length(cor_idx)

    pair_names <- obs_names[cor_idx]
    pairs <- do.call(rbind, lapply(strsplit(pair_names, "~~"), function(nm)
      c(match(nm[1], indicators), match(nm[2], indicators))))
    if (anyNA(pairs))
      stop("Polychoric pair names do not match indicators.", call. = FALSE)
    swap          <- pairs[, 1L] < pairs[, 2L]
    pairs[swap, ] <- pairs[swap, c(2L, 1L), drop = FALSE]
    r_obs <- R_poly[cbind(pairs[, 1L], pairs[, 2L])]
    message("  ", n_pairs, " polychoric pairs, ", n_wls, " WLS elements.")

    # -- Echelon parameterization ------------------------------------------------
    fixed_ij <- list()
    for (jj in seq_len(k))
      for (ii in seq_len(jj - 1L))
        fixed_ij[[length(fixed_ij) + 1L]] <- c(ii, jj)
    fixed_vec_idx <- vapply(fixed_ij, function(ij)
      as.integer((ij[2L] - 1L) * n_i + ij[1L]), integer(1L))
    free_vec_idx <- setdiff(seq_len(n_i * k), fixed_vec_idx)
    n_free       <- length(free_vec_idx)

    # -- Stage 4: DWLS optimization --------------------------------------------
    message("Stage 4: DWLS optimization via nlminb (", n_free, " free echelon loadings)...")

    sigma_fn <- function(L) tcrossprod(L)[cbind(pairs[, 1L], pairs[, 2L])]

    f_dwls_ech <- function(free_params) {
      lam               <- rep(0, n_i * k)
      lam[free_vec_idx] <- free_params
      res               <- r_obs - sigma_fn(matrix(lam, n_i, k))
      sum(W_cor * res * res)
    }

    grad_dwls_ech <- function(free_params) {
      lam               <- rep(0, n_i * k)
      lam[free_vec_idx] <- free_params
      L                 <- matrix(lam, n_i, k)
      wr                <- W_cor * (r_obs - sigma_fn(L))
      WR                <- matrix(0, n_i, n_i)
      WR[cbind(pairs[, 1L], pairs[, 2L])] <- wr
      WR[cbind(pairs[, 2L], pairs[, 1L])] <- wr
      as.vector(-2 * WR %*% L)[free_vec_idx]
    }

    start_vec                <- as.vector(efa_rot)
    start_vec[fixed_vec_idx] <- 0

    opt <- tryCatch(
      nlminb(start     = start_vec[free_vec_idx],
             objective = f_dwls_ech,
             gradient  = grad_dwls_ech,
             control   = list(iter.max = 3000L, eval.max = 8000L,
                              rel.tol  = 1e-10, x.tol   = 1e-10)),
      error = function(e) stop("DWLS nlminb failed: ", conditionMessage(e), call. = FALSE)
    )
    if (opt$convergence != 0)
      warning("DWLS optimization may not have fully converged (code=", opt$convergence,
              "): ", opt$message, call. = FALSE)

    lam_unrot               <- rep(0, n_i * k)
    lam_unrot[free_vec_idx] <- opt$par
    L_unrot <- matrix(lam_unrot, n_i, k)
    message("  Converged (code=", opt$convergence, ")  F = ", round(opt$objective, 4))

    # -- Stage 5: Oblique targetQ rotation --------------------------------------
    # Key: normalize=FALSE matches lavaan's row.weights="none" for target rotation.
    # Using Kaiser normalization (normalize=TRUE) causes a systematic loading gap.
    message("Stage 5: Oblique targetQ rotation via GPArotation (", n_starts, " random starts)...")
    # Seed local to this call (restored on exit) so random-start rotation is
    # reproducible without disturbing the caller's RNG stream.
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv))
                  get(".Random.seed", envir = .GlobalEnv) else NULL
    on.exit({
      if (!is.null(old_seed))
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      else
        suppressWarnings(rm(".Random.seed", envir = .GlobalEnv))
    }, add = TRUE)
    set.seed(42L)
    best_fn  <- Inf
    best_rot <- NULL

    for (s in seq_len(n_starts)) {
      if (s == 1L) {
        L_start <- L_unrot
      } else {
        T0      <- qr.Q(qr(matrix(stats::rnorm(k * k), k, k)))
        L_start <- L_unrot %*% T0
      }
      rot_s <- tryCatch(
        suppressWarnings(
          GPArotation::targetQ(L_start, Target = target, maxit = 10000L, normalize = FALSE)
        ),
        error = function(e) NULL
      )
      fn_s <- if (!is.null(rot_s) && all(is.finite(rot_s$loadings))) {
        0.5 * sum(rot_s$loadings[cross_mask]^2)
      } else NA_real_
      if (isTRUE(is.finite(fn_s) && fn_s < best_fn)) {
        best_fn  <- fn_s
        best_rot <- rot_s
      }
    }

    if (is.null(best_rot))
      stop("targetQ rotation failed on all ", n_starts, " starts.", call. = FALSE)

    message("  Best criterion = ", round(best_fn, 6))
    L_rot   <- best_rot$loadings
    Phi_rot <- best_rot$Phi
    colnames(L_rot) <- factor_names
    rownames(L_rot) <- indicators

    # Sign correction: flip each column so primary loadings sum positive
    for (j in seq_len(k)) {
      prim_idx <- which(is.na(target[, j]))
      if (sum(L_rot[prim_idx, j]) < 0) {
        L_rot[, j]   <- -L_rot[, j]
        Phi_rot[j, ] <- -Phi_rot[j, ]
        Phi_rot[, j] <- -Phi_rot[, j]
      }
    }
    colnames(Phi_rot) <- rownames(Phi_rot) <- factor_names

    # --- Stage 6: SEs (sandwich + rotation delta method) ----------------------
    message("Stage 6: Standard errors (sandwich ACM + numDeriv rotation-delta Jacobian)...")

    # 6a: Jacobian at L_unrot (free columns) for sandwich ACM
    G_unrot <- matrix(0, n_pairs, n_i * k)
    for (p in seq_len(n_pairs)) {
      mi <- pairs[p, 1L]; ni <- pairs[p, 2L]
      for (j in seq_len(k)) {
        G_unrot[p, (j - 1L) * n_i + mi] <- L_unrot[ni, j]
        G_unrot[p, (j - 1L) * n_i + ni] <- L_unrot[mi, j]
      }
    }
    G_free    <- G_unrot[, free_vec_idx, drop = FALSE]
    W_diag    <- diag(W_cor)
    GfWGf     <- crossprod(G_free, W_diag %*% G_free)
    GfWGf_inv <- MASS::ginv(GfWGf)
    ACM_free  <- GfWGf_inv %*%
      crossprod(G_free, W_diag %*% Gamma_cor %*% W_diag %*% G_free) %*%
      GfWGf_inv / (n - 1L)

    # 6b: Rotation Jacobian (numerical)
    SE_r <- matrix(NA_real_, n_i, k, dimnames = list(indicators, factor_names))
    J_rot <- NULL
    if (requireNamespace("numDeriv", quietly = TRUE)) {
      rot_fn_oblique <- function(free_params) {
        lam               <- rep(0, n_i * k)
        lam[free_vec_idx] <- free_params
        L <- matrix(lam, n_i, k)
        rot <- suppressWarnings(
          GPArotation::targetQ(L, Target = target, maxit = 10000L, normalize = FALSE)
        )
        L_r <- rot$loadings
        for (jj in seq_len(k)) {
          prim <- which(is.na(target[, jj]))
          if (sum(L_r[prim, jj]) < 0) L_r[, jj] <- -L_r[, jj]
        }
        as.vector(L_r)
      }
      J_rot <- tryCatch(
        numDeriv::jacobian(rot_fn_oblique, opt$par, method = "Richardson"),
        error = function(e) { warning("numDeriv failed: ", conditionMessage(e)); NULL }
      )
    } else {
      message("  Install 'numDeriv' for rotation-corrected standard errors.")
    }

    if (!is.null(J_rot)) {
      ACM_rot <- J_rot %*% ACM_free %*% t(J_rot)
      SE_r[]  <- sqrt(pmax(diag(ACM_rot), 0))
    }

    # -- Stage 7: WLSMV chi-square (full WLS, correlations + thresholds) ------
    message("Stage 7: WLSMV chi-square (scaled-and-shifted, Asparouhov & Muthen 2010)...")

    resid_v   <- r_obs - sigma_fn(L_unrot)
    T_naive   <- (n - 1L) * sum(W_cor * resid_v * resid_v)
    df_wlsmv  <- n_pairs - n_free

    # Full WLS approach: embed model Jacobian in n_wls-dimensional space.
    # Mplus treats thresholds as free model parameters, so the Jacobian Delta
    # has columns for both loading params AND threshold params.  Without the
    # threshold columns, the hat matrix H projects onto too few dimensions,
    # U = I - H is too large, and the scaling factors a/b are wrong.
    W_full_diag <- diag(as.matrix(W_full))
    W_full_sqrt <- diag(sqrt(W_full_diag))

    n_thr   <- n_wls - n_pairs
    thr_idx <- setdiff(seq_len(n_wls), cor_idx)

    G_full_wls <- matrix(0, n_wls, n_free + n_thr)
    G_full_wls[cor_idx, seq_len(n_free)]        <- G_free          # loading params
    G_full_wls[thr_idx, n_free + seq_len(n_thr)] <- diag(n_thr)    # threshold params

    GtWG_fwls     <- crossprod(G_full_wls, diag(W_full_diag) %*% G_full_wls)
    GtWG_fwls_inv <- solve(GtWG_fwls)

    H_fwls          <- W_full_sqrt %*% G_full_wls %*% GtWG_fwls_inv %*%
                       t(G_full_wls) %*% W_full_sqrt
    Gamma_full_tilde <- W_full_sqrt %*% Gamma_full %*% W_full_sqrt
    U_fwls          <- diag(n_wls) - H_fwls
    UG_fwls         <- U_fwls %*% Gamma_full_tilde

    a_scale <- sum(diag(UG_fwls))
    b_scale <- sum(diag(UG_fwls %*% UG_fwls))

    if (!is.finite(b_scale) || b_scale <= 0) {
      warning("WLSMV b_scale <= 0; falling back to mean-adjusted.", call. = FALSE)
      b_scale <- a_scale^2 / df_wlsmv
    }
    c_scale <- sqrt(df_wlsmv / b_scale)
    T_wlsmv <- (T_naive - a_scale) * c_scale + df_wlsmv
    pval    <- pchisq(T_wlsmv, df = df_wlsmv, lower.tail = FALSE)

    # Null model (independence: all correlations = 0, thresholds free)
    # Null Jacobian: only threshold params -> H_null zeros out correlation block,
    # passes threshold block -> U_null = I for correlations, 0 for thresholds.
    resid_null           <- wls_obs_full
    resid_null[cor_idx]  <- r_obs  # correlations are the residuals (null = 0)
    resid_null[-cor_idx] <- 0      # thresholds: model = observed
    T_null_naive <- (n - 1L) * sum(W_full_diag * resid_null * resid_null)
    df_null      <- n_pairs

    U_null <- diag(n_wls)
    U_null[thr_idx, thr_idx] <- 0            # threshold dimensions projected out
    UG_null <- U_null %*% Gamma_full_tilde
    a_null  <- sum(diag(UG_null))
    b_null  <- sum(diag(UG_null %*% UG_null))
    c_null  <- sqrt(df_null / b_null)
    T_null_wlsmv <- (T_null_naive - a_null) * c_null + df_null

    # -- Stage 8: Fit indices -----------------------------------------------------
    denom_cfi <- max(T_null_wlsmv - df_null, 0)
    CFI   <- if (denom_cfi == 0) 1 else 1 - max(T_wlsmv - df_wlsmv, 0) / denom_cfi
    TLI   <- if (df_null == 0 || T_null_wlsmv == 0) NA_real_ else
      (T_null_wlsmv / df_null - T_wlsmv / df_wlsmv) / (T_null_wlsmv / df_null - 1)
    RMSEA <- sqrt(max((T_wlsmv - df_wlsmv) / (df_wlsmv * (n - 1L)), 0))
    SRMR  <- sqrt(sum(resid_v * resid_v) / n_wls)

    message("  chi2(", df_wlsmv, ") = ", round(T_wlsmv, 3),
            "  p = ", round(pval, 4),
            "  CFI = ", round(CFI, 3),
            "  TLI = ", round(TLI, 3),
            "  RMSEA = ", round(RMSEA, 3),
            "  SRMR = ", round(SRMR, 3))

    wlsmv_stats <- list(
      chisq      = T_wlsmv,  df       = df_wlsmv, pvalue  = pval,
      cfi        = CFI,      tli      = TLI,      rmsea   = RMSEA,
      srmr       = SRMR,
      chisq_null = T_null_wlsmv, df_null = df_null,
      T_naive    = T_naive,  a_scale  = a_scale,
      b_scale    = b_scale,  c_scale  = c_scale
    )

    message("Done.")
    return(structure(
      list(
        lavaan_fit           = fit_aux,
        syntax               = simple_syntax,
        nfactors             = k,
        rotation             = "Target oblique (DWLS + post-hoc targetQ, A&M 2009)",
        factor_names         = factor_names,
        indicators           = indicators,
        polychoric           = R_poly,
        efa_loadings         = efa_rot,
        unrotated_loadings   = L_unrot,
        rotation_criterion   = best_fn,
        rotated_loadings     = L_rot,
        std_rotated_loadings = L_rot,  # polychoric -> already STDYX
        se_loadings          = SE_r,
        factor_correlations  = Phi_rot,
        wlsmv_stats          = wlsmv_stats,
        estimator            = "WLSMV",
        call                 = mc
      ),
      class = c("esem_fit_ordered", "esem_fit")
    ))
  }

  # =================--==========================---===============================
  # METHOD: "lavaan" (default)
  # ======================--=====================================================

  # -- Build efa() block model syntax --------------------------------------------
  # All k factors in one EFA block: lavaan handles identification internally
  # (echelon constraints) and computes rotation-corrected SEs via delta method.
  # This is the Asparouhov & Muthen (2009) algorithm, same as Mplus.
  lhs          <- paste0("efa('esem')*", factor_names, collapse = " + ")
  model_syntax <- paste0(lhs, " =~ ", paste(indicators, collapse = " + "))

  # -- Map rotation string to lavaan arguments -----------------------------------
  is_orth    <- grepl("T$|varimax", rotation)
  lav_rot    <- sub("T$|Q$", "", tolower(rotation))   # "target","geomin","oblimin","varimax"
  rot_args   <- list(orthogonal = is_orth)
  if (!is.null(target)) rot_args$target <- target     # NA=free, 0=target-zero

  # -- Fit: efa() + WLSMV + rotation --------------------------------------------
  message("Stage 1: Fitting ESEM via lavaan efa() + WLSMV + ", rotation, " rotation...")
  message("  (lavaan computes polychoric correlations and rotation-corrected SEs internally)")

  cfa_args <- list(
    model            = model_syntax,
    data             = data,
    ordered          = indicators,
    std.lv           = std.lv,
    parameterization = "theta",
    missing          = missing,
    rotation         = lav_rot,
    rotation.args    = rot_args
  )
  if (!is.null(group))       cfa_args$group       <- group
  if (!is.null(group_equal)) cfa_args$group.equal <- group_equal
  cfa_args <- c(cfa_args, list(...))

  fit <- tryCatch(
    do.call(lavaan::cfa, cfa_args),
    error = function(e) stop("ESEM efa() WLSMV failed: ", conditionMessage(e), call. = FALSE)
  )

  conv <- lavaan::lavInspect(fit, "converged")
  message("Done. Estimator: ", fit@Options$estimator, " | Converged: ", conv)
  if (!conv)
    warning("ESEM efa() WLSMV did not converge. Check model specification.", call. = FALSE)

  # -- Heywood case correction (rotation-only) ---------------------------------
  heywood_log <- NULL
  corrected_L <- NULL
  if (isTRUE(heywood_fix) && is.null(group)) {
    cached_pre <- .cache_std_loadings(fit)
    if (!is.null(cached_pre)) {
      Phi_rot <- tryCatch(lavaan::lavInspect(fit, "cor.lv"), error = function(e) {
        ss_tmp <- lavaan::standardizedsolution(fit)
        .extract_phi_from_ss(ss_tmp, factor_names)
      })
      if (is.list(Phi_rot)) Phi_rot <- Phi_rot[[1L]]

      result <- .heywood_retry_loop(
        L_rot             = cached_pre$L,
        Phi_rot           = Phi_rot,
        initial_target    = target,
        indicators        = indicators,
        factor_names      = factor_names,
        original_rotation = rotation,
        max_rounds        = 15L
      )
      heywood_log <- result$log
      corrected_L <- result$L_stdyx
    }
  }

  # -- Extract standardised loading matrix --------------------------------------
  if (!is.null(corrected_L)) {
    loading_matrix <- corrected_L
  } else {
    loading_matrix <- matrix(NA_real_, n_i, nfactors,
                             dimnames = list(indicators, factor_names))
    ss <- tryCatch(lavaan::standardizedsolution(fit), error = function(e) NULL)
    if (!is.null(ss)) {
      ss_lv <- ss[ss$op == "=~", , drop = FALSE]
      for (rr in seq_len(nrow(ss_lv))) {
        fi <- match(ss_lv$lhs[rr], factor_names)
        ii <- match(ss_lv$rhs[rr], indicators)
        if (!is.na(fi) && !is.na(ii)) loading_matrix[ii, fi] <- ss_lv$est.std[rr]
      }
    }
  }

  # -- Return --------------------------------------------------------------------
  structure(
    list(
      lavaan_fit   = fit,
      syntax       = model_syntax,
      nfactors     = nfactors,
      rotation     = if (!is.null(heywood_log)) heywood_log$final_rotation else rotation,
      factor_names = factor_names,
      indicators   = indicators,
      polychoric   = NULL,        # handled internally by lavaan
      efa_loadings = loading_matrix,
      estimator    = fit@Options$estimator,
      heywood_log  = heywood_log,
      call         = mc
    ),
    class = c("esem_fit_ordered", "esem_fit")
  )
}


# -- Analytical rotation Jacobian (implicit function theorem) -----------------
# Computes d(vec(L_rotated)) / d(free_echelon_params) exactly for orthogonal
# targetT rotation, avoiding numerical differentiation through the iterative
# GPArotation optimizer.  Eliminates local-minimum sensitivity and numerical
# noise in the delta-method SEs.
#
# At the rotation solution, L'G is symmetric (optimality condition) where
# L = A %*% T (rotated loadings) and G = 2*W*(L - Target).  A perturbation
# dA in the unrotated loadings induces dT = T %*% Phi (Phi skew-symmetric)
# to maintain optimality.  Phi is determined by a k(k-1)/2 linear system
# from the implicit function theorem.
#
# Reference: Asparouhov & Muthen (2009), SEM, 16(3), 397-438.
.analytical_rot_jacobian <- function(L_unrot, T_mat, Target, free_vec_idx) {
  p <- nrow(L_unrot)
  k <- ncol(L_unrot)
  n_free <- length(free_vec_idx)
  n_rot  <- (k * (k - 1L)) %/% 2L

  L   <- L_unrot %*% T_mat
  W   <- ifelse(is.na(Target), 0, 1)
  Tgt <- Target; Tgt[is.na(Tgt)] <- 0
  G   <- 2 * W * (L - Tgt)

  # Verify rotation optimality: L'G should be symmetric
  LtG <- crossprod(L, G)
  skew_residual <- max(abs(LtG - t(LtG))) / 2
  if (skew_residual > 0.05) {
    warning("Rotation optimality poorly satisfied (max|skew| = ",
            round(skew_residual, 4), "); falling back to numerical Jacobian.",
            call. = FALSE)
    return(NULL)
  }

  # Upper-triangle indices for skew-symmetric parameterization
  ut <- which(upper.tri(matrix(0, k, k)), arr.ind = TRUE)

  # Rotation Hessian F on the orthogonal manifold:
  # F(Phi) = skew(2 L'(W * (L Phi)) - Phi (L'G))
  S <- LtG
  F_mat <- matrix(0, n_rot, n_rot)
  for (r in seq_len(n_rot)) {
    Phi_r <- matrix(0, k, k)
    Phi_r[ut[r, 1L], ut[r, 2L]] <-  1
    Phi_r[ut[r, 2L], ut[r, 1L]] <- -1
    LPhi   <- L %*% Phi_r
    result <- 2 * crossprod(L, W * LPhi) - Phi_r %*% S
    skew_r <- (result - t(result)) / 2
    for (s in seq_len(n_rot))
      F_mat[s, r] <- skew_r[ut[s, 1L], ut[s, 2L]]
  }

  cond <- tryCatch(rcond(F_mat), error = function(e) 0)
  if (cond < 1e-12) {
    warning("Rotation Hessian near-singular (rcond = ", signif(cond, 3),
            "); falling back to numerical Jacobian.", call. = FALSE)
    return(NULL)
  }
  F_inv <- solve(F_mat)

  # Build Jacobian column by column: one column per free echelon parameter
  J <- matrix(0, p * k, n_free)
  for (j in seq_len(n_free)) {
    idx   <- free_vec_idx[j]
    col_j <- (idx - 1L) %/% p + 1L
    row_j <- (idx - 1L) %% p + 1L

    # b = skew(T' dA' G + 2 L' (W * (dA T)))
    # dA has a single 1 at (row_j, col_j); exploit sparsity.
    TtdAtG   <- outer(T_mat[col_j, ], G[row_j, ])
    WdAT_row <- W[row_j, ] * T_mat[col_j, ]
    LtWdAT   <- 2 * outer(L[row_j, ], WdAT_row)
    b_full   <- TtdAtG + LtWdAT
    b_skew   <- (b_full - t(b_full)) / 2
    b_vec    <- vapply(seq_len(n_rot), function(s)
      b_skew[ut[s, 1L], ut[s, 2L]], numeric(1))

    phi_vec <- -F_inv %*% b_vec
    Phi <- matrix(0, k, k)
    for (s in seq_len(n_rot)) {
      Phi[ut[s, 1L], ut[s, 2L]] <-  phi_vec[s]
      Phi[ut[s, 2L], ut[s, 1L]] <- -phi_vec[s]
    }

    # dL = dA T + L Phi  (dA T only touches row row_j)
    dL <- L %*% Phi
    dL[row_j, ] <- dL[row_j, ] + T_mat[col_j, ]
    J[, j] <- as.vector(dL)
  }
  J
}

#' Bifactor ESEM for Ordered-Categorical Data
#'
#' Fits B-ESEM on ordered-categorical indicators. Called from \code{\link{besem}}
#' when \code{ordered} is set. For Mplus-aligned loadings and fit indices use
#' \code{method = "rotation"} (default; also used by \code{\link{run_comparison}}).
#'
#' @inheritSection doc_estimator_paths Estimator paths
#' @inheritSection doc_estimator_paths Missing-data defaults
#' @inheritSection doc_estimator_paths The lavaan_fit slot on custom WLSMV fits
#'
#' @param data A \code{data.frame} of observed ordered indicators.
#' @param specific_factors Named list of specific factor -> item assignments.
#' @param indicators Character vector of all indicator names. If \code{NULL},
#'   derived from \code{specific_factors}.
#' @param g_name Character. General factor name. Default \code{"G"}.
#' @param method Character. \code{"rotation"} (default): custom DWLS/WLSMV +
#'   orthogonal \code{targetT} rotation (Mplus-aligned B-ESEM). \code{"set-esem"}:
#'   lavaan WLSMV **bifactor CFA** with non-primary specific loadings fixed at
#'   \code{0*} and starts from polychoric EFA---does **not** match Mplus B-ESEM
#'   loadings (restricted model; use only for debugging or when you want zero
#'   specific-factor cross-loadings).
#' @param n_starts Integer. Random rotation starts. Default \code{30L}
#'   (matches Mplus).
#' @param r_obs_override Optional observed correlation matrix to use instead
#'   of polychoric estimation.
#' @param group Character. Grouping variable for multi-group models.
#' @param group_equal Character vector of lavaan equality constraints.
#' @param missing Character. Missing data handling. Default \code{"listwise"}.
#'   \code{\link{run_comparison}} passes \code{"pairwise"} from
#'   \code{\link{specify_model}} when fitting via the pipeline. Set explicitly
#'   when comparing with \code{\link{esem_ordered}} (default \code{"pairwise"}).
#' @param std.lv Logical. Default \code{TRUE}.
#' @param ... Additional arguments passed to \code{lavaan::cfa()}.
#'
#' @return An object of class \code{c("besem_fit_ordered","besem_fit","esem_fit")}.
#' For \code{method = "rotation"}, \code{lavaan_fit} is auxiliary (one-factor CFA);
#' use \code{std_rotated_loadings}, \code{wlsmv_stats}, \code{\link{std_loadings}},
#' and \code{\link[lavaan:fitMeasures]{fitMeasures}(x)}. For \code{method = "set-esem"},
#' \code{lavaan_fit} is the fitted bifactor CFA.
#'
#' @section Method \code{"set-esem"} vs \code{"rotation"}:
#' \code{rotation} estimates cross-loadings (targeted toward zero), then rotates---
#' same estimand as Mplus B-ESEM WLSMV. \code{set-esem} **fixes** non-primary
#' specific loadings at zero in lavaan; G and primary loadings are re-optimized
#' under that harder constraint, so loadings and fit differ from Mplus.
#'
#' @seealso \code{\link{esem_ordered}}, \code{\link{besem}},
#'   \code{\link{run_comparison}}
#'
#' @examples
#' \dontrun{
#' fit_b_ord <- besem_ordered(
#'   data = Rdata,
#'   specific_factors = list(
#'     EX = items_EX, MD = items_MD, CI = items_CI
#'   )
#' )
#' summary(fit_b_ord, fit.measures = TRUE, standardized = TRUE)
#' }
#'
#' @importFrom psych polychoric fa
#' @importFrom lavaan cfa
#' @export
besem_ordered <- function(data,
                           specific_factors,
                           indicators     = NULL,
                           g_name         = "G",
                           method         = c("rotation", "set-esem"),
                           n_starts       = 30L,
                           r_obs_override = NULL,
                           group          = NULL,
                           group_equal    = NULL,
                           missing        = "listwise",
                           std.lv         = TRUE,
                           ...) {

  mc     <- match.call()
  method <- match.arg(method)

  if (!requireNamespace("psych",       quietly = TRUE))
    stop("Package 'psych' is required. Install: install.packages('psych')", call. = FALSE)
  if (!requireNamespace("GPArotation", quietly = TRUE))
    stop("Package 'GPArotation' is required. Install: install.packages('GPArotation')", call. = FALSE)

  if (is.null(indicators))
    indicators <- unlist(specific_factors, use.names = FALSE)

  specific_names   <- names(specific_factors)
  all_factor_names <- c(g_name, specific_names)
  k                <- length(all_factor_names)   # total factors (G + specific)
  n_i              <- length(indicators)
  n                <- nrow(data)

  if (is.null(r_obs_override))
    .assert_polychoric_compatible(data, indicators)

  # Build bifactor target: G column = NA (free target), specific = 1/NA
  btgt     <- make_bifactor_target(specific_factors = specific_factors,
                                   indicators       = indicators,
                                   g_name           = g_name)
  # make_bifactor_target produces no NAs (all 0s/1s), so btgt_rot == btgt
  btgt_rot <- btgt

  # -- Stage 1: Polychoric correlations (shared by both methods) ----------------
  if (!is.null(r_obs_override)) {
    message("Stage 1: Using supplied polychoric matrix (r_obs_override, e.g. parsed from Mplus .out)...")
    # Match Mplus row/col names (often UPPER case) to R's mixed-case indicators
    rn  <- rownames(r_obs_override)
    idx <- match(tolower(indicators), tolower(rn))
    if (anyNA(idx))
      stop("r_obs_override item names do not match indicators.\n",
           "  Missing: ", paste(indicators[is.na(idx)], collapse = ", "),
           call. = FALSE)
    R_poly <- r_obs_override[idx, idx]
    rownames(R_poly) <- colnames(R_poly) <- indicators
    message("  Matched ", length(indicators), " items from supplied matrix.")
  } else {
    message("Stage 1: Polychoric correlations via psych::polychoric()...")
    poly_result <- tryCatch(
      psych::polychoric(data[, indicators, drop = FALSE]),
      error = function(e) stop("polychoric() failed: ", conditionMessage(e), call. = FALSE)
    )
    R_poly <- poly_result$rho
  }

  # -- Stage 2: Unrotated EFA + bifactor targetT rotation (shared) --------------
  message("Stage 2: Bifactor EFA (factanal + GPArotation::targetT, ",
          k, " factors: 1 general + ", k - 1L, " specific)...")
  fa_unrotated <- tryCatch(
    factanal(covmat = R_poly, factors = k, n.obs = n, rotation = "none"),
    error = function(e) stop("Unrotated EFA failed: ", conditionMessage(e), call. = FALSE)
  )
  A <- unclass(fa_unrotated$loadings)

  efa_rot <- tryCatch({
    suppressWarnings(
      GPArotation::targetT(A, Target = btgt_rot, maxit = 10000L, eps = 1e-4)
    )$loadings
  }, error = function(e) stop("Bifactor EFA rotation failed: ", conditionMessage(e), call. = FALSE))
  colnames(efa_rot) <- all_factor_names
  rownames(efa_rot) <- indicators

  # -- Orthogonality constraints (shared by both methods) -----------------------
  orth_lines <- unlist(lapply(seq_len(k - 1), function(i)
    lapply(seq(i + 1L, k), function(j)
      paste0(all_factor_names[i], " ~~ 0*", all_factor_names[j]))))

  # -- Shared lavaan CFA fitting helper -----------------------------------------
  .fit_wlsmv <- function(syntax) {
    args <- list(model            = syntax,
                 data             = data,
                 ordered          = indicators,
                 std.lv           = TRUE,
                 parameterization = "theta",
                 missing          = missing)
    if (!is.null(group))       args$group       <- group
    if (!is.null(group_equal)) args$group.equal <- group_equal
    args <- c(args, list(...))
    do.call(lavaan::cfa, args)
  }

  # ============================================================================
  # METHOD: "rotation"
  # Fit unrestricted WLSMV model (echelon identification), then apply
  # post-hoc bifactor targetT rotation with numDeriv SE propagation.
  # Fit statistics are rotation-invariant -> match Mplus BESEM df/CFI/RMSEA.
  # Reference: Asparouhov & Muthen (2009); Ogasawara (2000).
  # ============================================================================
  if (method == "rotation") {

    # ==========================================================================
    # Asparouhov & Muthen (2009) DWLS/WLSMV EFA algorithm -- implemented from
    # scratch to match Mplus BESEM WLSMV output.
    #
    # The 4-factor bifactor model cannot be estimated directly by lavaan (the
    # G factor's large loadings on all items create near-singular Hessian with
    # 60+ free parameters). Instead we:
    #   Stage 3: Extract DWLS weight matrices (W, Gamma) from an auxiliary CFA.
    #   Stage 4: Minimize the DWLS objective over all n_ixk loadings via nlminb.
    #   Stage 5: Apply bifactor targetT rotation to the unrotated DWLS solution.
    #   Stage 6: Sandwich ACM + numDeriv rotation delta method -> SEs.
    #   Stage 7: WLSMV chi-square with mean-variance adjustment.
    #   Stage 8: CFI, TLI, RMSEA, SRMR from WLSMV statistics.
    #   Stage 9: Return (loadings from polychoric fit are already STDYX).
    # ==========================================================================

    # -- Stage 3: Extract DWLS weight matrices via auxiliary 1-factor CFA ------
    # W (wls.v) and Gamma depend on data + thresholds only, not factor structure.
    # Any converging WLSMV fit on the same data yields the same matrices.
    message("Stage 3: Extracting DWLS weight matrices (W, Gamma) from auxiliary lavaan CFA...")
    simple_syntax <- paste0(g_name, " =~ ", paste(indicators, collapse = " + "))
    fit_aux <- tryCatch(
      .fit_wlsmv(simple_syntax),
      error = function(e) stop("Auxiliary WLSMV CFA failed: ", conditionMessage(e), call. = FALSE)
    )

    wls_obs_raw  <- lavaan::lavInspect(fit_aux, "wls.obs")
    W_raw        <- lavaan::lavInspect(fit_aux, "wls.v")
    Gamma_raw    <- lavaan::lavInspect(fit_aux, "gamma")

    # Multi-group auxiliary CFA: lavInspect returns a list (one per group).
    # W and Gamma depend only on marginal thresholds, so any group's matrices
    # are equivalent for the DWLS objective.  Use the first group.
    if (is.list(wls_obs_raw) && !is.null(names(wls_obs_raw)) &&
        !any(grepl("~~", names(wls_obs_raw)))) {
      wls_obs_raw <- wls_obs_raw[[1L]]
      W_raw       <- W_raw[[1L]]
      Gamma_raw   <- Gamma_raw[[1L]]
    }

    wls_obs_full <- wls_obs_raw
    W_full       <- W_raw
    Gamma_full   <- as.matrix(Gamma_raw)

    obs_names    <- names(wls_obs_full)
    n_wls        <- length(wls_obs_full)   # n_pairs + n_thresholds (Mplus SRMR denominator)
    cor_idx      <- grep("~~", obs_names)

    if (length(cor_idx) == 0)
      stop("Cannot identify polychoric correlation elements in wls.obs. ",
           "This is a lavaan version issue -- please report.", call. = FALSE)

    # wls.v is a square diagonal matrix (not a vector) -- extract the diagonal.
    # as.vector() would serialize the FULL matrix and give wrong indices.
    # W and Gamma depend on thresholds (marginal distributions) not on the
    # bivariate polychoric estimates, so they are taken from lavaan throughout.
    W_cor     <- diag(as.matrix(W_full))[cor_idx]
    Gamma_cor <- Gamma_full[cor_idx, cor_idx, drop = FALSE]
    n_pairs   <- length(cor_idx)

    # Parse (row_item, col_item) index pairs for each polychoric correlation
    pair_names <- obs_names[cor_idx]
    pairs <- do.call(rbind, lapply(strsplit(pair_names, "~~"), function(nm)
      c(match(nm[1], indicators), match(nm[2], indicators))))
    if (anyNA(pairs))
      stop("Polychoric pair names in wls.obs do not match `indicators`. ",
           "Ensure all indicators are listed correctly.", call. = FALSE)

    # Ensure row index > col index (lower-triangle convention)
    swap          <- pairs[, 1L] < pairs[, 2L]
    pairs[swap, ] <- pairs[swap, c(2L, 1L), drop = FALSE]

    # r_obs: polychoric correlations for the DWLS objective.
    # We use R_poly (psych::polychoric or r_obs_override) rather than lavaan's
    # wls.obs polychorics.  Empirically, psych and Mplus agree to < 0.001 on
    # all pairs, whereas lavaan can differ by up to ~0.024 on highly-correlated
    # items (e.g., batCI1~~batCI2 = 0.825 in lavaan vs 0.801 in Mplus/psych).
    # That discrepancy shifts the DWLS objective surface and can push the
    # rotation to a different local minimum.  Using psych/Mplus polychorics
    # throughout gives rotation solutions that match Mplus to within rounding.
    # W and Gamma still come from lavaan -- they depend on thresholds, not on
    # the bivariate correlation estimates.
    # R_poly is either psych::polychoric (default) or r_obs_override (user-supplied).
    # Both paths use the same polychoric matrix set in Stage 1.
    r_obs <- R_poly[cbind(pairs[, 1L], pairs[, 2L])]
    if (!is.null(r_obs_override))
      message("  DWLS objective uses r_obs_override polychoric correlations.")

    message("  ", n_pairs, " polychoric correlation pairs extracted.")

    # -- Echelon parameterization -----------------------------------------------
    # ML EFA (factanal) satisfies LL' + Psi = R_poly exactly, so the DWLS
    # objective is trivially 0 at the ML solution -- giving T_naive = 0 and no
    # meaningful chi-square.  Mplus avoids this by using an echelon
    # parameterization: fix the k*(k-1)/2 upper-triangular elements of the
    # first k rows of L to 0.  This breaks the exact-reproduction property and
    # gives T_naive > 0, from which the WLSMV adjustment can be computed.
    # The minimum value of T_naive is invariant to which valid echelon is used.
    message("Stage 3b: Setting up echelon parameterization (", k*(k-1L)/2L,
            " loading constraints)...")

    # Fixed indices: upper triangle of first k rows, column-major order
    # (i < j, i <= k): lambda[1,2]=0, lambda[1,3]=0, lambda[2,3]=0, lambda[1,4]=0, lambda[2,4]=0, lambda[3,4]=0
    fixed_ij <- list()
    for (jj in seq_len(k)) {
      for (ii in seq_len(jj - 1L)) fixed_ij[[length(fixed_ij) + 1L]] <- c(ii, jj)
    }
    fixed_vec_idx <- vapply(fixed_ij, function(ij) (ij[2L] - 1L) * n_i + ij[1L], integer(1L))
    free_vec_idx  <- setdiff(seq_len(n_i * k), fixed_vec_idx)
    n_free        <- length(free_vec_idx)   # n_i*k - k*(k-1)/2  (= 66 for 18 items, 4 factors)

    # -- Stage 4: DWLS optimization (echelon-constrained, 66 free params) -------
    # Start from bifactor EFA rotation with constrained elements zeroed out.
    # With 6 zeros enforced the ML exact-reproduction no longer holds and
    # T_naive > 0 at the constrained minimum.
    message("Stage 4: DWLS optimization via nlminb (", n_free, " free echelon loadings)...")

    sigma_fn <- function(L) {
      LLt <- tcrossprod(L)
      LLt[cbind(pairs[, 1L], pairs[, 2L])]
    }

    f_dwls_ech <- function(free_params) {
      lam               <- rep(0, n_i * k)
      lam[free_vec_idx] <- free_params
      L                 <- matrix(lam, n_i, k)
      res               <- r_obs - sigma_fn(L)
      sum(W_cor * res * res)
    }

    grad_dwls_ech <- function(free_params) {
      lam               <- rep(0, n_i * k)
      lam[free_vec_idx] <- free_params
      L                 <- matrix(lam, n_i, k)
      wr                <- W_cor * (r_obs - sigma_fn(L))
      WR                <- matrix(0, n_i, n_i)
      WR[cbind(pairs[, 1L], pairs[, 2L])] <- wr
      WR[cbind(pairs[, 2L], pairs[, 1L])] <- wr
      grad_full <- as.vector(-2 * WR %*% L)
      grad_full[free_vec_idx]
    }

    # Start: bifactor EFA solution with 6 constrained elements zeroed
    start_vec                    <- as.vector(efa_rot)
    start_vec[fixed_vec_idx]     <- 0
    start_free                   <- start_vec[free_vec_idx]

    opt <- tryCatch(
      nlminb(start     = start_free,
             objective = f_dwls_ech,
             gradient  = grad_dwls_ech,
             control   = list(iter.max = 3000L, eval.max = 8000L,
                              rel.tol  = 1e-10, x.tol   = 1e-10)),
      error = function(e) stop("DWLS nlminb failed: ", conditionMessage(e), call. = FALSE)
    )
    if (opt$convergence != 0)
      warning("DWLS optimization may not have fully converged (code=", opt$convergence,
              "): ", opt$message, call. = FALSE)

    lam_unrot                  <- rep(0, n_i * k)
    lam_unrot[free_vec_idx]    <- opt$par
    L_unrot                    <- matrix(lam_unrot, n_i, k)
    message("  Converged (code=", opt$convergence, ")  F = ", round(opt$objective, 4))

    # -- Stage 5: Bifactor targetT rotation with random restarts ----------------
    # Mplus TARGET convention: NA = free (G column + primary loadings),
    #                          0  = penalise toward zero (cross-loadings only).
    # btgt_rot has 1s for "primary/G" and 0s for cross-loadings.
    # We convert 1s -> NA so GPForth only penalises cross-loadings, exactly
    # matching Mplus.  btgt_rot (0s/1s) is kept for the cross-loading criterion.
    btgt_na      <- btgt_rot
    btgt_na[btgt_rot == 1L] <- NA_real_
    cross_mask   <- btgt_rot == 0L      # elements targeted toward 0 (cross-loadings)

    # Mplus uses 30 random orthogonal starting matrices to escape local minima
    # of the rotation criterion.  We match that default (controllable via n_starts).
    # Approach: pre-rotate L_unrot by a random orthogonal T0, then call
    # targetT from the identity -- equivalent to a random Tmat start but avoids
    # numerical instability from passing large-norm Tmat directly to GPForth.
    # The composite rotation of the best start is stored as best_Th for Stage 6.
    message("Stage 5: Bifactor targetT rotation via GPArotation (",
            n_starts, " random starts)...")
    # Seed local to this call (restored on exit) so random-start rotation is
    # reproducible without disturbing the caller's RNG stream.
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv))
                  get(".Random.seed", envir = .GlobalEnv) else NULL
    on.exit({
      if (!is.null(old_seed))
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      else
        suppressWarnings(rm(".Random.seed", envir = .GlobalEnv))
    }, add = TRUE)
    set.seed(42L)   # reproducible random starts
    best_fn  <- Inf
    best_rot <- NULL
    best_Th  <- diag(k)   # composite rotation of winning start (for Stage 6)

    for (s in seq_len(n_starts)) {
      if (s == 1L) {
        L_start <- L_unrot          # deterministic: identity start
        T0      <- diag(k)
      } else {
        T0      <- qr.Q(qr(matrix(stats::rnorm(k * k), k, k)))
        L_start <- L_unrot %*% T0  # pre-rotate; targetT refines from identity
      }

      rot_s <- tryCatch(
        suppressWarnings(
          GPArotation::targetT(L_start, Target = btgt_na, maxit = 10000L, eps = 1e-8)
        ),
        error = function(e) NULL
      )
      # Criterion: sum of squared cross-loadings (elements targeted toward 0).
      # Using NA for primary/G elements means GPForth only penalises cross-loadings,
      # matching Mplus's TARGET criterion.  No fn.value dependency needed.
      fn_s <- if (!is.null(rot_s) && all(is.finite(rot_s$loadings))) {
        0.5 * sum(rot_s$loadings[cross_mask]^2)
      } else NA_real_
      if (isTRUE(is.finite(fn_s) && fn_s < best_fn)) {
        best_fn  <- fn_s
        best_rot <- rot_s
        best_Th  <- T0 %*% rot_s$Th   # composite: random pre-rotation + refinement
      }
    }

    if (is.null(best_rot))
      stop("targetT rotation failed on all ", n_starts, " starts.", call. = FALSE)

    message("  Best criterion = ", round(best_fn, 6),
            " (", n_starts, " starts, seed 42)")
    L_rot  <- best_rot$loadings
    colnames(L_rot) <- all_factor_names
    rownames(L_rot) <- indicators

    # Sign correction: with NA targets the sign of each factor is arbitrary.
    # Convention: flip each column so its primary loadings sum to positive.
    # G factor (col 1): sum over all items; specific factors: sum over primary items.
    for (j in seq_len(k)) {
      prim_idx <- if (j == 1L) seq_len(n_i) else which(btgt_rot[, j] == 1L)
      if (sum(L_rot[prim_idx, j]) < 0) L_rot[, j] <- -L_rot[, j]
    }

    # -- Stage 6: SEs via sandwich estimator + rotation delta method -----------
    # Two-step procedure following Asparouhov & Muthen (2009):
    #   (a) Sandwich ACM for the n_free echelon parameters at L_unrot
    #       (the DWLS optimum), scaled by 1/(n-1).
    #   (b) Delta method via J_rot (Jacobian of rotation map) to propagate
    #       the ACM from echelon space to rotated loading space.
    #
    # G_mat (at L_rot, rank-66 pseudo-inverse) is kept for Stage 7's hat matrix.
    # For SEs, a separate G_free is built at L_unrot using the free columns only.
    message("Stage 6: Standard errors (sandwich ACM + numDeriv rotation-delta Jacobian)...")

    # -- 6a: Jacobian G at L_rot (full 72 columns) -- needed for Stage 7 H_mat -
    #   d[LL']_{mn}/dlambda_{mj} = lambda_{nj},  d[LL']_{mn}/dlambda_{nj} = lambda_{mj}
    G_mat <- matrix(0, n_pairs, n_i * k)
    for (p in seq_len(n_pairs)) {
      mi <- pairs[p, 1L]; ni <- pairs[p, 2L]
      for (j in seq_len(k)) {
        G_mat[p, (j - 1L) * n_i + mi] <- L_rot[ni, j]
        G_mat[p, (j - 1L) * n_i + ni] <- L_rot[mi, j]
      }
    }
    W_diag   <- diag(W_cor)
    GtWG     <- crossprod(G_mat, W_diag %*% G_mat)   # 72x72, rank 66
    GtWG_inv <- MASS::ginv(GtWG)                      # pseudo-inverse for Stage 7

    # -- 6b: Sandwich ACM for free echelon parameters at L_unrot --------------
    # Evaluate Jacobian at the DWLS optimum (L_unrot), not at the rotated
    # solution.  Only free columns (free_vec_idx) are included; the 6 fixed
    # columns contribute zero gradient and are excluded.
    G_unrot <- matrix(0, n_pairs, n_i * k)
    for (p in seq_len(n_pairs)) {
      mi <- pairs[p, 1L]; ni <- pairs[p, 2L]
      for (j in seq_len(k)) {
        G_unrot[p, (j - 1L) * n_i + mi] <- L_unrot[ni, j]
        G_unrot[p, (j - 1L) * n_i + ni] <- L_unrot[mi, j]
      }
    }
    G_free    <- G_unrot[, free_vec_idx, drop = FALSE]   # n_pairs x n_free
    GfWGf     <- crossprod(G_free, W_diag %*% G_free)    # n_free x n_free
    GfWGf_inv <- MASS::ginv(GfWGf)
    # Lavaan's Gamma is the large-sample ACM: sqrt(N)(rho_hat-rho) -> N(0,Gamma).
    # Divide by (n-1) to recover Var(rho_hat); lavaan uses (N-1) convention.
    ACM_free  <- GfWGf_inv %*%
      crossprod(G_free, W_diag %*% Gamma_cor %*% W_diag %*% G_free) %*%
      GfWGf_inv / (n - 1L)                             # n_free x n_free

    # -- 6c: Delta method -- rotate ACM from echelon to rotated loading space --
    # J_rot: Jacobian of (free echelon params -> vectorised rotated loadings),
    # shape (n_i*k) x n_free = 72 x 66.
    # Analytical Jacobian via implicit function theorem (preferred); falls back
    # to numerical Jacobian (numDeriv) if the analytical path fails.
    # ACM_rot = J_rot %*% ACM_free %*% t(J_rot)  ->  (n_i*k)x(n_i*k).
    SE_r <- matrix(NA_real_, n_i, k, dimnames = list(indicators, all_factor_names))

    J_rot <- .analytical_rot_jacobian(L_unrot, best_Th, btgt_na, free_vec_idx)

    if (is.null(J_rot)) {
      # Fallback: numerical Jacobian through iterative rotation
      message("  Falling back to numerical Jacobian (numDeriv)...")
      if (requireNamespace("numDeriv", quietly = TRUE)) {
        rot_fn_ech <- function(free_params) {
          lam               <- rep(0, n_i * k)
          lam[free_vec_idx] <- free_params
          as.vector(suppressWarnings(
            GPArotation::targetT(matrix(lam, n_i, k),
                                 Target = btgt_na,
                                 maxit  = 10000L)
          )$loadings)
        }
        J_rot <- tryCatch(
          numDeriv::jacobian(rot_fn_ech, opt$par, method = "Richardson"),
          error = function(e) { warning("numDeriv failed: ", conditionMessage(e)); NULL }
        )
      } else {
        message("  Install 'numDeriv' for rotation-corrected standard errors.")
      }
    }

    if (!is.null(J_rot)) {
      ACM_rot <- J_rot %*% ACM_free %*% t(J_rot)
      SE_r[]  <- sqrt(pmax(diag(ACM_rot), 0))
    } else {
      se_full <- rep(NA_real_, n_i * k)
      se_full[free_vec_idx] <- sqrt(pmax(diag(ACM_free), 0))
      SE_r[] <- matrix(se_full, n_i, k)
    }

    # -- Stage 7: WLSMV chi-square (mean-and-variance adjusted, "scaled.shifted") -
    # Matches Mplus WLSMV formula: T* = (T_naive - a) x c + df
    # where a = tr(U Gamma),  b = tr[(U Gamma)^2],  c = sqrt(df/b)
    # Reference: Satorra & Bentler (2010); Muthen (1993).
    # df = n_pairs - n_free  (= 153 - 66 = 87 for 18 items, 4 factors)
    #
    # CRITICAL: The WLSMV scaling must be computed in the FULL WLS space
    # (correlations + thresholds), not just the correlation space.  Mplus treats
    # thresholds as free model parameters.  Without threshold columns in the
    # Jacobian, the hat matrix H projects onto too few dimensions, U = I - H is
    # too large, and the scaling factors a/b are biased -> wrong chi-square.
    message("Stage 7: WLSMV chi-square (scaled-and-shifted, Asparouhov & Muthen 2010)...")

    sigma_rot <- sigma_fn(L_rot)
    resid_rot <- r_obs - sigma_rot
    T_naive   <- (n - 1L) * sum(W_cor * resid_rot * resid_rot)

    df_wlsmv  <- n_pairs - n_free

    # Full WLS space: Jacobian with loading params + threshold params
    W_full_diag <- diag(as.matrix(W_full))
    n_thr       <- n_wls - n_pairs
    thr_idx     <- setdiff(seq_len(n_wls), cor_idx)
    n_load_cols <- ncol(G_mat)   # n_i*k (may be > n_free due to echelon)

    G_full_wls <- matrix(0, n_wls, n_load_cols + n_thr)
    G_full_wls[cor_idx, seq_len(n_load_cols)]          <- G_mat       # loading params
    G_full_wls[thr_idx, n_load_cols + seq_len(n_thr)]  <- diag(n_thr) # threshold params

    GtWG_fwls     <- crossprod(G_full_wls, diag(W_full_diag) %*% G_full_wls)
    GtWG_fwls_inv <- MASS::ginv(GtWG_fwls)  # pseudo-inverse (G_mat is rank-deficient)

    W_full_sqrt <- diag(sqrt(W_full_diag))

    H_fwls          <- W_full_sqrt %*% G_full_wls %*% GtWG_fwls_inv %*%
                       t(G_full_wls) %*% W_full_sqrt
    Gamma_full_tilde <- W_full_sqrt %*% Gamma_full %*% W_full_sqrt
    U_fwls          <- diag(n_wls) - H_fwls
    UG              <- U_fwls %*% Gamma_full_tilde
    a_scale     <- sum(diag(UG))               # tr(U Gamma)  -- trace.UGamma
    b_scale     <- sum(diag(UG %*% UG))        # tr[(U Gamma)^2] -- trace.UGamma2

    if (!is.finite(b_scale) || b_scale <= 0) {
      warning("WLSMV b_scale <= 0; falling back to mean-adjusted statistic.", call. = FALSE)
      b_scale <- a_scale^2 / df_wlsmv
    }
    c_scale  <- sqrt(df_wlsmv / b_scale)       # = 1 / scaling.factor
    T_wlsmv  <- (T_naive - a_scale) * c_scale + df_wlsmv
    pval     <- pchisq(T_wlsmv, df = df_wlsmv, lower.tail = FALSE)

    # Null model (independence: all correlations = 0, thresholds free)
    # Null Jacobian has only threshold params -> H_null zeros out correlation
    # block, passes threshold block -> U_null = I for correlations, 0 for thresholds.
    T_null_naive <- (n - 1L) * sum(W_cor * r_obs * r_obs)
    df_null      <- n_pairs
    U_null <- diag(n_wls)
    U_null[thr_idx, thr_idx] <- 0            # threshold dimensions projected out
    UG_null <- U_null %*% Gamma_full_tilde
    a_null  <- sum(diag(UG_null))
    b_null  <- sum(diag(UG_null %*% UG_null))
    c_null  <- sqrt(df_null / b_null)
    T_null_wlsmv <- (T_null_naive - a_null) * c_null + df_null

    # -- Stage 8: Fit indices ---------------------------------------------------
    denom_cfi <- max(T_null_wlsmv - df_null, 0)
    CFI  <- if (denom_cfi == 0) 1 else
      1 - max(T_wlsmv - df_wlsmv, 0) / denom_cfi
    TLI  <- if (df_null == 0 || T_null_wlsmv == 0) NA_real_ else
      (T_null_wlsmv / df_null - T_wlsmv / df_wlsmv) / (T_null_wlsmv / df_null - 1)
    RMSEA <- sqrt(max((T_wlsmv - df_wlsmv) / (df_wlsmv * (n - 1L)), 0))
    # Mplus SRMR denominator: total WLS observed statistics = n_pairs + n_thresholds.
    # Mplus includes the full wls.obs vector in the denominator; threshold residuals
    # are identically 0 (thresholds are free params) so only the denominator changes.
    # This matches Mplus SRMR convention; remaining gap vs Mplus is from rotation.
    SRMR  <- sqrt(sum(resid_rot * resid_rot) / n_wls)

    message("  chi2(", df_wlsmv, ") = ", round(T_wlsmv, 3),
            "  p = ", round(pval, 4),
            "  CFI = ", round(CFI, 3),
            "  TLI = ", round(TLI, 3),
            "  RMSEA = ", round(RMSEA, 3),
            "  SRMR = ", round(SRMR, 3))

    wlsmv_stats <- list(
      chisq      = T_wlsmv,  df       = df_wlsmv, pvalue  = pval,
      cfi        = CFI,      tli      = TLI,      rmsea   = RMSEA,
      srmr       = SRMR,
      chisq_null = T_null_wlsmv, df_null = df_null,
      T_naive    = T_naive,  a_scale  = a_scale,
      b_scale    = b_scale,  c_scale  = c_scale
    )

    # -- Stage 9: STDYX loading matrix -----------------------------------------
    # Fitting to polychoric CORRELATIONS (unit item variance) means L_rot
    # elements are already fully standardised (STDYX equivalent).
    L_r_std <- L_rot
    message("Done.")

    return(structure(
      list(
        lavaan_fit           = fit_aux,   # auxiliary 1-factor CFA (weight extraction)
        syntax               = simple_syntax,
        nfactors             = k,
        rotation             = "TargetT bifactor (DWLS + post-hoc rotation, A&M 2009)",
        factor_names         = all_factor_names,
        g_name               = g_name,
        specific_factors     = specific_factors,
        indicators           = indicators,
        polychoric           = R_poly,
        efa_loadings         = efa_rot,
        unrotated_loadings   = L_unrot,    # echelon-parameterised DWLS solution pre-rotation
        rotation_criterion   = best_fn,    # best sum-of-squared cross-loadings achieved
        rotated_loadings     = L_rot,
        std_rotated_loadings = L_r_std,
        se_loadings          = SE_r,
        wlsmv_stats          = wlsmv_stats,
        estimator            = "WLSMV",
        call                 = mc
      ),
      class = c("besem_fit_ordered", "besem_fit", "esem_fit")
    ))
  }

  # ============================================================================
  # METHOD: "set-esem"  (fallback / comparison)
  # G loads on all items; specific factors load only on primary items (=0 cross).
  # Simpler but does not match Mplus B-ESEM loadings.
  # ============================================================================

  message("Stage 3: Building bifactor set-ESEM lavaan syntax...")
  get_bf_primary <- function(fn) which(btgt[, fn] == 1)

  syntax_lines_se <- vapply(seq_len(k), function(j) {
    fn        <- all_factor_names[j]
    primary_j <- get_bf_primary(fn)
    terms <- vapply(seq_len(n_i), function(i) {
      sv <- round(efa_rot[i, j], 4)
      if (fn == g_name || i %in% primary_j) paste0("start(", sv, ")*", indicators[i])
      else                                   paste0("0*", indicators[i])
    }, character(1))
    paste0(fn, " =~ ", paste(terms, collapse = " + "))
  }, character(1))

  syntax_se <- paste(
    c("# bifactory: Bifactor set-ESEM (polychoric + WLSMV)",
      "# G: all items free; Specific: primary free, non-primary = 0",
      "", syntax_lines_se, "", "# Orthogonality", orth_lines),
    collapse = "\n"
  )

  message("Stage 4: Fitting bifactor set-ESEM via lavaan::cfa() (WLSMV + theta)...")
  fit <- tryCatch(
    .fit_wlsmv(syntax_se),
    error = function(e) stop("Bifactor set-ESEM failed: ", conditionMessage(e), call. = FALSE)
  )
  message("Done. Estimator: ", fit@Options$estimator,
          " | Converged: ", lavaan::lavInspect(fit, "converged"))

  structure(
    list(
      lavaan_fit       = fit,
      syntax           = syntax_se,
      nfactors         = k,
      rotation         = "TargetT (orthogonal, set-ESEM)",
      factor_names     = all_factor_names,
      g_name           = g_name,
      specific_factors = specific_factors,
      indicators       = indicators,
      polychoric       = R_poly,
      efa_loadings     = efa_rot,
      estimator        = fit@Options$estimator,
      call             = mc
    ),
    class = c("besem_fit_ordered", "besem_fit", "esem_fit")
  )
}
