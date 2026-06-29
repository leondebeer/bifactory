#' Exploratory Structural Equation Modeling
#'
#' Fits ESEM using \pkg{lavaan}'s native \code{efa()} block. For **continuous**
#' data (default), estimation and rotation are integrated in lavaan's ML
#' optimiser---the same single-step approach as Mplus \code{(*1)} syntax. For
#' **ordered** indicators, supply \code{ordered} and the call is routed to
#' \code{\link{esem_ordered}} (set-ESEM / WLSMV; not identical to continuous
#' \code{esem()}).
#'
#' @inheritSection doc_estimator_paths Estimator paths
#' @inheritSection doc_estimator_paths Missing-data defaults
#' @inheritSection doc_estimator_paths Heywood fix and loadings
#'
#' @section How this differs from a two-stage EFA + CFA workaround:
#' Many R implementations run EFA first, extract the loading matrix, then paste
#' it into a CFA as starting values. That is an approximation. \pkg{bifactory}
#' instead uses \pkg{lavaan}'s \code{efa()} block syntax (available since
#' lavaan 0.6-12), which estimates the rotation and the SEM parameters
#' simultaneously -- exactly as Mplus does with the \code{(*1)} syntax.
#'
#' @param data A \code{data.frame} containing the observed indicators (and
#'   optionally a grouping variable).
#' @param nfactors Integer. Number of latent factors.
#' @param indicators Character vector of item names to include. If \code{NULL}
#'   (default), all columns of \code{data} except \code{group} are used.
#' @param rotation Character. Rotation criterion. Passed to
#'   \code{lavaan::cfa(rotation = ...)}. Common choices:
#'   \itemize{
#'     \item \code{"geomin"} (default) -- oblique geomin, the Mplus default
#'     \item \code{"target"} -- requires \code{target} matrix
#'     \item \code{"varimax"} -- orthogonal varimax
#'     \item \code{"oblimin"} -- oblique oblimin
#'     \item \code{"none"} -- no rotation (confirmatory EFA)
#'   }
#'   The full list of lavaan-supported rotations is in
#'   \code{?lavaan::efaRotate}.
#' @param target A numeric matrix (items x factors) for target rotation.
#'   Create with \code{\link{make_target}}. Required when
#'   \code{rotation = "target"}.
#' @param estimator Character. lavaan estimator. Default \code{"MLR"}
#'   (ML with Huber-White robust SEs and Satorra-Bentler scaled chi-square).
#'   Use \code{"WLSMV"} for ordered indicators.
#' @param std.lv Logical. Fix factor variances to 1 for identification?
#'   Default \code{TRUE} (recommended for ESEM).
#' @param ordered Character vector of ordered-categorical item names.
#'   When non-\code{NULL}, routes to \code{\link{esem_ordered}} (WLSMV). Default
#'   \code{method} there is \code{"lavaan"}; see \code{?esem_ordered} for the
#'   custom \code{"rotation"} path.
#' @param group Character. Name of a grouping variable in \code{data} for
#'   multi-group ESEM (configural by default).
#' @param group_equal Character vector of lavaan parameter labels to constrain
#'   equal across groups, e.g. \code{"loadings"} (metric) or
#'   \code{c("loadings", "intercepts")} (scalar).
#' @param missing Character. Missing data handling for **continuous** fits.
#'   Default \code{"listwise"}; use \code{"fiml"} for full-information ML.
#'   When \code{ordered} is set, passed to \code{\link{esem_ordered}} (default
#'   \code{"pairwise"} there unless you override).
#' @param factor_names Optional character vector of length \code{nfactors}.
#'   Defaults to \code{F1, F2, ...}.
#' @param rotation_args Named list of extra arguments passed to the rotation
#'   function (e.g. \code{list(geomin.epsilon = 0.001)}).  When
#'   \code{rotation = "geomin"} (or any geomin variant) and
#'   \code{geomin.epsilon} is not supplied, it defaults to
#'   \code{0.0001} / \code{0.001} / \code{0.01} for 2 / 3 / 4+ factors
#'   respectively, matching Mplus's defaults.
#' @param heywood_fix Logical. Retry rotation with Cholesky unrotation if a
#'   standardised loading exceeds 1 (single-group continuous fits only).
#'   Default \code{TRUE}. If correction runs, \code{\link{std_loadings}} may
#'   differ from \code{lavaan::standardizedsolution(lavaan_fit(x))}; see
#'   \link[=doc_estimator_paths]{Heywood fix and loadings}.
#' @param ... Additional arguments passed to \code{lavaan::cfa()}.
#'
#' @return An object of class \code{"esem_fit"} containing:
#' \describe{
#'   \item{\code{lavaan_fit}}{The \pkg{lavaan} fit object for continuous ESEM.
#'     For ordered / custom WLSMV paths see \code{\link{esem_ordered}} and
#'     \link[=doc_estimator_paths]{The lavaan_fit slot on custom WLSMV fits}.}
#'   \item{\code{syntax}}{The lavaan model string that was estimated.}
#'   \item{\code{nfactors}}{Number of factors.}
#'   \item{\code{rotation}}{Rotation method used.}
#'   \item{\code{factor_names}}{Factor names.}
#'   \item{\code{indicators}}{Item names included in the model.}
#'   \item{\code{call}}{The matched call.}
#' }
#'
#' @details
#' ## Mplus equivalence
#'
#' | Mplus syntax | bifactory equivalent |
#' |---|---|
#' | `F1-F3 BY y1-y15 (*1);` | `esem(data, nfactors = 3)` |
#' | `(*1)` with geomin (default) | `rotation = "geomin"` |
#' | `(*1)` with target | `rotation = "target", target = tgt` |
#' | `GROUPING = g;` | `group = "g"` |
#' | Metric invariance | `group_equal = "loadings"` |
#' | Scalar invariance | `group_equal = c("loadings", "intercepts")` |
#'
#' ## Why not automate cross-loadings from a CFA?
#' Adding cross-loadings stepwise from modification indices is a different
#' (exploratory CFA) approach and is not recommended: it capitalises on chance,
#' inflates Type I error, and produces a different model on every dataset.
#' ESEM instead estimates \emph{all} cross-loadings simultaneously, with
#' rotation acting as a mathematical penalty for complexity -- giving a
#' reproducible, theory-neutral solution.
#'
#' @seealso
#' \code{\link{make_target}} for target matrices,
#' \code{\link{esem_ordered}} for ordered-categorical data,
#' \code{\link{esem_compare}} for ESEM vs CFA comparison,
#' \code{\link{std_loadings}} for the standardised loading matrix.
#'
#' @examples
#' data("HolzingerSwineford1939", package = "lavaan")
#' d <- HolzingerSwineford1939[, paste0("x", 1:9)]
#'
#' # Basic ESEM: 3 factors, geomin rotation (Mplus default)
#' fit <- esem(d, nfactors = 3)
#' round(std_loadings(fit), 2)
#' lavaan::fitMeasures(lavaan_fit(fit), c("cfi", "tli", "rmsea", "srmr"))
#'
#' # Named factors
#' fit2 <- esem(d, nfactors = 3,
#'              factor_names = c("Visual", "Textual", "Speed"))
#'
#' \donttest{
#' # Target rotation
#' tgt <- make_target(list(Vis = 1:3, Txt = 4:6, Spd = 7:9), nitems = 9)
#' fit3 <- esem(d, nfactors = 3, rotation = "target", target = tgt)
#'
#' # Multi-group configural ESEM
#' fit_mg <- esem(HolzingerSwineford1939, nfactors = 3,
#'                indicators = paste0("x", 1:9), group = "sex")
#' }
#'
#' @importFrom lavaan cfa
#' @export
esem <- function(data,
                 nfactors,
                 indicators    = NULL,
                 rotation      = "geomin",
                 target        = NULL,
                 estimator     = "MLR",
                 std.lv        = TRUE,
                 ordered       = NULL,
                 group         = NULL,
                 group_equal   = NULL,
                 missing       = "listwise",
                 factor_names  = NULL,
                 rotation_args = list(),
                 heywood_fix   = TRUE,
                 ...) {

  mc <- match.call()

  # -- 0. Input validation ---------------------------------------------------
  if (!is.data.frame(data) && !is.matrix(data))
    stop("`data` must be a data.frame or matrix.", call. = FALSE)

  nfactors <- as.integer(nfactors)
  if (nfactors < 1L)
    stop("`nfactors` must be a positive integer.", call. = FALSE)

  if (grepl("target", rotation, ignore.case = TRUE) && is.null(target))
    stop("Provide a `target` matrix when using target rotation. ",
         "See ?make_target.", call. = FALSE)

  # Resolve indicator names
  all_cols   <- colnames(data)
  indicators <- if (is.null(indicators)) {
    setdiff(all_cols, group)
  } else {
    bad <- setdiff(indicators, all_cols)
    if (length(bad))
      stop("Indicator(s) not found in `data`: ", paste(bad, collapse = ", "),
           call. = FALSE)
    indicators
  }

  nitems <- length(indicators)
  if (nitems < nfactors * 2L)
    stop("Too few indicators (", nitems, ") for ", nfactors, " factors. ",
         "Need at least ", nfactors * 2L, ".", call. = FALSE)

  # Factor names
  if (is.null(factor_names)) {
    factor_names <- paste0("F", seq_len(nfactors))
  } else {
    if (length(factor_names) != nfactors)
      stop("`factor_names` length must equal `nfactors` (", nfactors, ").",
           call. = FALSE)
    factor_names <- make.names(factor_names)
  }

  # -- Ordered indicators: lavaan's efa() block does not support WLSMV ---------
  # Route automatically to esem_ordered() which implements the set-ESEM pipeline
  # (polychoric correlations -> EFA rotation -> WLSMV CFA), matching Mplus behaviour.
  if (!is.null(ordered)) {
    message("Ordered indicators detected. lavaan's efa() block does not support ",
            "WLSMV rotation.\nRouting to esem_ordered() (set-ESEM + WLSMV pipeline).")
    return(esem_ordered(
      data         = data,
      nfactors     = nfactors,
      indicators   = indicators,
      rotation     = rotation,
      target       = target,
      factor_names = factor_names,
      group        = group,
      group_equal  = group_equal,
      missing      = missing,
      std.lv       = std.lv,
      heywood_fix  = heywood_fix,
      ...
    ))
  }

  # -- 1. Build native lavaan ESEM syntax ------------------------------------
  # Uses lavaan's efa() block -- a single integrated estimation with rotation
  # built into the ML optimiser. This is NOT a two-stage workaround.
  #
  # Equivalent Mplus syntax:
  #   F1-Fn BY y1-ym (*1);
  syntax <- .build_esem_syntax(
    factor_names = factor_names,
    indicators   = indicators,
    nfactors     = nfactors
  )

  # -- 2. Fit in lavaan ------------------------------------------------------
  cfa_args <- list(
    model     = syntax,
    data      = data,
    std.lv    = std.lv,
    estimator = estimator,
    missing   = missing,
    rotation  = rotation
  )

  if (!is.null(ordered))    cfa_args$ordered      <- ordered
  if (!is.null(group))      cfa_args$group        <- group
  if (!is.null(group_equal)) cfa_args$group.equal <- group_equal

  # Auto-set Geomin epsilon from nfactors to match Mplus convention
  # (Asparouhov & Muthen 2009, footnote 6): 0.0001 for 2 factors, 0.001 for
  # 3 factors, 0.01 for 4+ factors.  lavaan's default is 0.001 regardless of
  # nfactors; only inject when the user did not supply one explicitly.
  if (grepl("geomin", rotation, ignore.case = TRUE) &&
      is.null(rotation_args$geomin.epsilon)) {
    rotation_args$geomin.epsilon <- if (nfactors == 2L) 0.0001
                                    else if (nfactors == 3L) 0.001
                                    else 0.01
  }

  # Merge rotation_args, injecting target matrix if needed
  if (!is.null(target)) rotation_args$target <- target
  if (length(rotation_args)) cfa_args$rotation.args <- rotation_args

  # User overrides (e.g. se = "robust.huber.white")
  cfa_args <- c(cfa_args, list(...))

  fit <- tryCatch(
    do.call(lavaan::cfa, cfa_args),
    error = function(e) {
      stop(
        "lavaan::cfa() failed:\n  ", conditionMessage(e),
        "\n\nRequires lavaan >= 0.6-12. Check: packageVersion('lavaan')",
        call. = FALSE
      )
    }
  )

  # -- 3. Heywood case correction (rotation-only) ----------------------------
  heywood_log    <- NULL
  corrected_L    <- NULL
  if (isTRUE(heywood_fix) && is.null(group)) {
    cached_pre <- .cache_std_loadings(fit)
    if (!is.null(cached_pre)) {
      # Use lavInspect for Phi -- parameterization-invariant
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

  # -- 4. Cache standardized solution ------------------------------------------
  cached <- if (!is.null(corrected_L)) {
    se_mat <- .cache_std_loadings(fit)$SE
    list(L = corrected_L, SE = se_mat)
  } else if (is.null(group)) {
    .cache_std_loadings(fit)
  } else NULL

  # -- 5. Assemble and return --------------------------------------------------
  structure(
    list(
      lavaan_fit           = fit,
      syntax               = syntax,
      nfactors             = nfactors,
      rotation             = if (!is.null(heywood_log)) heywood_log$final_rotation else rotation,
      factor_names         = factor_names,
      indicators           = indicators,
      std_rotated_loadings = if (!is.null(cached)) cached$L  else NULL,
      se_loadings          = if (!is.null(cached)) cached$SE else NULL,
      heywood_log          = heywood_log,
      call                 = mc
    ),
    class = "esem_fit"
  )
}


# -- Internal helpers ----------------------------------------------------------

.build_esem_syntax <- function(factor_names, indicators, nfactors) {
  # The efa("label") prefix groups factors into one EFA block.
  # lavaan estimates ALL cross-loadings and applies rotation internally.
  #
  # Generated form:
  #   efa("esem")*F1 +
  #   efa("esem")*F2 +
  #   efa("esem")*F3 =~
  #       y1 + y2 + y3 + y4 + y5 + y6 +
  #       y7 + y8 + ...

  lhs <- paste(
    paste0('efa("esem")*', factor_names),
    collapse = " +\n    "
  )

  # Wrap RHS at 6 items per line for readability
  chunks   <- split(indicators, ceiling(seq_along(indicators) / 6))
  rhs_rows <- vapply(chunks, paste, character(1L), collapse = " + ")
  rhs      <- paste(rhs_rows, collapse = " +\n    ")

  paste0(
    "# bifactory: ESEM via lavaan native efa() block\n",
    "# Mplus equivalent: ",
    paste(factor_names, collapse = "-"),
    " BY ",
    indicators[1], "-", indicators[length(indicators)], " (*1);\n\n",
    lhs, " =~\n    ", rhs, "\n"
  )
}
