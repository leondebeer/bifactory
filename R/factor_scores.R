# -- Factor Score Extraction ---------------------------------------------------

#' Extract Latent Factor Scores
#'
#' Returns a data frame of estimated latent factor scores, one column per
#' factor and one row per observation.  Works with any \code{esem_fit},
#' plain \code{lavaan} object, or \code{esem_invariance} result.
#'
#' @param x An \code{esem_fit}, \code{lavaan} S4, or \code{esem_invariance}
#'   object.
#' @param method Character. Estimation method: \code{"regression"} (default,
#'   BLUP -- minimises MSE, scores are correlated when factors are correlated)
#'   or \code{"bartlett"} (unbiased, correct factor variance).  For WLSMV
#'   B-ESEM models with orthogonal rotation both methods return identical
#'   scores (\eqn{\Phi = I}).
#' @param level Character. Only used when \code{x} is an \code{esem_invariance}
#'   object.  \code{"auto"} (default) selects the most constrained fitted
#'   invariance level whose own transition \eqn{\Delta}CFI \eqn{\ge -0.010},
#'   including \code{"varcov"} and \code{"means"} when
#'   \code{esem_invariance()} was run with \code{through}.  Override
#'   with \code{"configural"}, \code{"weak"}, \code{"strong"},
#'   \code{"strict"}, \code{"varcov"} or \code{"means"}.  \code{"auto"}
#'   skips a level whose reported solution
#'   is inadmissible (a negative variance or a non-positive-definite matrix,
#'   recorded in \code{inv$notes}) even when its \eqn{\Delta}CFI passes.
#'   Silently ignored for plain fit objects.
#' @param dCFI_cutoff Numeric. The \eqn{\Delta}CFI threshold used by
#'   \code{level = "auto"} to accept an invariance level (a level qualifies
#'   when its transition \eqn{\Delta}CFI \eqn{\ge} \code{dCFI_cutoff}).
#'   Defaults to \code{-0.010} (Cheung & Rensvold, 2002).  Only used when
#'   \code{x} is an \code{esem_invariance} object and \code{level = "auto"}.
#' @param align Optional alignment of the score columns. \code{NULL} (default)
#'   returns scores in the model's native rotation orientation -- which can
#'   sign-flip or column-permute across reruns or solvers (lavaan vs Mplus)
#'   for B-ESEM target rotation. Supply one of:
#'   \itemize{
#'     \item \code{"canonical"} -- deterministic rule (columns sorted by
#'       \eqn{\Sigma\lambda^2} descending, sign-flipped so the largest
#'       \eqn{|\lambda|} per column is positive). No external reference;
#'       reproducible across reruns.
#'     \item \code{"group1"} -- multi-group only; rotate each group's scores
#'       to match group 1's loading orientation.
#'     \item a numeric matrix (single-group) or list of matrices (multi-group)
#'       of reference standardized loadings (rownames = items, colnames =
#'       factors). Useful for aligning to an external reference such as
#'       Mplus output via \code{extract_mplus_loadings()}.
#'     \item another fit object (\code{esem_fit}, \code{esem_invariance})
#'       -- aligns to that fit's loadings.
#'   }
#'   Internally calls \code{align_loadings()} to compute the per-group
#'   rotation \eqn{Q}, then rotates the score matrix by \eqn{Q}. Note: under
#'   \code{"canonical"} columns may be permuted, in which case the column
#'   names retain their original ordering -- treat them as positional after
#'   alignment. Not supported for \code{besem_fit_ordered} fits produced
#'   by \code{besem_ordered()} with custom WLSMV polychoric scoring.
#' @param ... Reserved for future use.
#'
#' @return A \code{data.frame} with one column per latent factor (named by
#'   the model's factor labels) and one row per observation.  A \code{group}
#'   column is prepended automatically for multi-group models.  Rows for
#'   missing observations contain \code{NA} in all factor columns.
#'
#' @details
#' For multi-group ordered B-ESEM scoring, a single pooled polychoric
#' correlation matrix is applied to every group. A warning is issued because
#' scores for non-reference groups may be biased when their correlation
#' structures differ. If the polychoric correlation matrix or the
#' score-information matrix is singular, \code{factor_scores()} raises an
#' informative error identifying the singular matrix rather than exposing a
#' bare \code{solve()} failure.
#'
#' @examples
#' data("HolzingerSwineford1939", package = "lavaan")
#' d <- HolzingerSwineford1939[, paste0("x", 1:9)]
#'
#' \donttest{
#' # Single fit
#' fit    <- esem(d, nfactors = 3)
#' scores <- factor_scores(fit)
#' head(scores)
#'
#' # Merge back to original data
#' d_with_scores <- cbind(d, factor_scores(fit))
#'
#' # Bartlett method
#' scores_b <- factor_scores(fit, method = "bartlett")
#' }
#'
#' @export
factor_scores <- function(x,
                          method      = c("regression", "bartlett"),
                          level       = "auto",
                          align       = NULL,
                          dCFI_cutoff = -0.010,
                          ...) {
  method <- match.arg(method)

  valid_levels <- c("auto", "configural", "weak", "strong", "strict",
                    "varcov", "means")
  if (!level %in% valid_levels)
    stop("level must be one of: auto, configural, weak, strong, strict, ",
         "varcov, means.", call. = FALSE)

  if (!is.numeric(dCFI_cutoff) || length(dCFI_cutoff) != 1L || is.na(dCFI_cutoff))
    stop("dCFI_cutoff must be a single numeric value.", call. = FALSE)

  # besem_fit_ordered comes in two flavours:
  #   (a) Custom from-scratch WLSMV (besem_ordered, method = "rotation"):
  #       lavaan_fit is echelon-parameterized, so lavPredict() would return
  #       wrong-factor scores. std_rotated_loadings + polychoric are stored
  #       so we do the math manually (Path 2).
  #   (b) Set-esem / invariance pipeline (lavaan::cfa with rotation="target"):
  #       lavaan_fit is a real rotated fit, lavPredict() returns correct
  #       scores. No custom fields stored -- route to Path 1.
  # besem_fit_continuous (bifactor.R:220) inherits esem_fit -> falls to Path 1 below.
  scores <- if (inherits(x, "esem_invariance")) {
    .fs_invariance(x, method = method, level = level, dCFI_cutoff = dCFI_cutoff)
  } else if (inherits(x, "besem_fit_ordered")) {
    if (!is.null(x$std_rotated_loadings))
      .fs_manual_wlsmv(x, method = method)
    else
      .fs_from_lavaan(x, method = method)
  } else if (inherits(x, "esem_fit") || inherits(x, "lavaan")) {
    .fs_from_lavaan(x, method = method)
  } else {
    stop("x must be an esem_fit, lavaan, or esem_invariance object.",
         call. = FALSE)
  }

  if (!is.null(align))
    scores <- .fs_apply_alignment(scores, x, target = align, level = level,
                                  dCFI_cutoff = dCFI_cutoff)

  scores
}
# -- Path 1: Standard lavaan (CFA, ESEM, continuous BESEM) --------------------

.fs_from_lavaan <- function(x, method) {
  # Unwrap esem_fit wrapper; plain lavaan S4 passes through as-is
  lav <- if (inherits(x, "esem_fit")) x$lavaan_fit else x

  # Factor names: prefer stored names, fall back to lavaan's own list
  fnames <- if (inherits(x, "esem_fit") && !is.null(x$factor_names))
    x$factor_names
  else
    lavaan::lavNames(lav, "lv")

  if (length(fnames) == 0L)
    stop("factor_scores: could not determine factor names from the model.",
         call. = FALSE)

  # lavPredict returns a matrix (single-group) or named list of matrices (multi-group)
  # lavaan 0.6.21 accepts "regression", "bartlett", "Bartlett", and "EBM" -- pass through directly
  pred <- lavaan::lavPredict(lav, method = method)

  if (is.list(pred)) {
    # Multi-group: bind rows, prepend group indicator
    grp_labels <- lavaan::lavInspect(lav, "group.label")
    if (is.null(grp_labels) || length(grp_labels) != length(pred))
      grp_labels <- seq_along(pred)

    dfs <- mapply(function(mat, grp) {
      df        <- as.data.frame(mat)
      if (ncol(df) != length(fnames))
        warning("factor_scores: lavPredict returned ", ncol(df),
                " factor column(s) but ", length(fnames),
                " factor name(s) found; names may be misaligned.", call. = FALSE)
      names(df) <- fnames[seq_len(ncol(df))]
      df$group  <- grp
      df[, c("group", setdiff(names(df), "group")), drop = FALSE]
    }, pred, grp_labels, SIMPLIFY = FALSE)

    result <- do.call(rbind, dfs)

    # Sanity check: compare against expected row count from original data
    expected_rows <- tryCatch(
      sum(sapply(lavaan::lavInspect(lav, "data"), nrow)),
      error = function(e) NULL
    )
    if (!is.null(expected_rows) && nrow(result) != expected_rows)
      warning("factor_scores: row count (", nrow(result), ") does not match ",
              "original data (", expected_rows, " rows) -- check for missing ",
              "observations.", call. = FALSE)

    return(result)
  }

  # Single-group
  df        <- as.data.frame(pred)
  if (ncol(df) != length(fnames))
    warning("factor_scores: lavPredict returned ", ncol(df),
            " factor column(s) but ", length(fnames),
            " factor name(s) found; names may be misaligned.", call. = FALSE)
  names(df) <- fnames[seq_len(ncol(df))]

  expected_rows <- tryCatch(nrow(lavaan::lavInspect(lav, "data")),
                            error = function(e) NULL)
  if (!is.null(expected_rows) && nrow(df) != expected_rows)
    warning("factor_scores: row count (", nrow(df), ") does not match ",
            "original data (", expected_rows, " rows) -- check for missing ",
            "observations.", call. = FALSE)
  df
}

# -- Path 2: WLSMV BESEM with custom rotation (manual scoring) ----------------

.fs_manual_wlsmv <- function(x, method) {
  # Guard: set-esem fallback path does not store rotated loadings
  Lambda <- x$std_rotated_loadings
  if (is.null(Lambda))
    stop("std_rotated_loadings not available -- this besem_fit_ordered was ",
         "produced by the set-esem fallback path, which does not support ",
         "manual factor score extraction.", call. = FALSE)

  R_poly <- x$polychoric
  if (is.null(R_poly))
    stop("Polychoric correlation matrix (x$polychoric) is NULL.", call. = FALSE)

  # Check positive definiteness of R before inversion
  eigs <- tryCatch(eigen(R_poly, symmetric = TRUE, only.values = TRUE)$values,
                   error = function(e) NULL)
  if (!is.null(eigs) && any(eigs <= 0))
    warning("factor_scores: polychoric correlation matrix is not positive ",
            "definite -- matrix inversion may be unstable.", call. = FALSE)

  # Raw data: lavInspect returns matrix (single-group) or list (multi-group)
  raw <- lavaan::lavInspect(x$lavaan_fit, "data")
  if (is.null(raw))
    stop("Raw data not available in the lavaan fit object (lavaan internal ",
         "state issue).", call. = FALSE)

  # Normalise to list for uniform handling
  if (!is.list(raw)) raw <- list(raw)

  fnames    <- x$factor_names
  if (is.null(fnames) || length(fnames) == 0L)
    stop("factor_scores: factor_names is NULL or empty on this besem_fit_ordered object.",
         call. = FALSE)

  # Note: x$polychoric is the pooled polychoric matrix from Stage 1 of besem_ordered().
  # For multi-group models, this pooled matrix is applied to all groups -- a standard
  # approximation when per-group polychorics are not separately estimated.
  if (length(raw) > 1L)
    warning("factor_scores: one pooled polychoric matrix is applied to all groups; ",
            "scores for non-reference groups may be biased if correlation structures differ.",
            call. = FALSE)
  R_inv <- tryCatch(
    solve(R_poly),
    error = function(e) stop(
      "factor_scores: the polychoric correlation matrix is singular; ",
      "this usually indicates redundant indicators or an ill-conditioned ",
      "correlation structure.", call. = FALSE))
  A         <- t(Lambda) %*% R_inv %*% Lambda   # k  x  k shared intermediate

  # Bartlett: W_bart = solve(A) %*% t(Lambda) %*% R_inv  [k  x  p]
  #           scores = X_std %*% t(W_bart)
  # Regression: W_reg = R_inv %*% Lambda %*% solve(A)    [p  x  k]
  #             scores = X_std %*% W_reg
  # For orthogonal B-ESEM (Phi = I): t(W_bart) == W_reg exactly,
  # so both methods produce identical scores. This is verified in the
  # validation script (test_factor_scores.R, Test 7).
  W_bart <- tryCatch(
    solve(A) %*% t(Lambda) %*% R_inv,
    error = function(e) stop(
      "factor_scores: the score-information matrix A is singular for ",
      "Bartlett scoring; this usually indicates redundant or unidentified ",
      "factor loadings.", call. = FALSE))                  # k  x  p
  W_reg <- tryCatch(
    R_inv %*% Lambda %*% solve(A),
    error = function(e) stop(
      "factor_scores: the score-information matrix A is singular for ",
      "regression scoring; this usually indicates redundant or unidentified ",
      "factor loadings.", call. = FALSE))                  # p  x  k
  W      <- if (method == "bartlett") t(W_bart) else W_reg  # always p  x  k

  score_group <- function(X_raw) {
    # Standardize each column (polychoric scale)
    X_std <- scale(X_raw, center = TRUE, scale = TRUE)
    # Note: scale() uses sample SD (n-1 denominator). For large n this is negligible
    # relative to the polychoric correlation metric used in Lambda.

    # Align columns: Lambda rows correspond to x$indicators order
    indicators <- x$indicators
    if (!is.null(indicators)) {
      missing_ind <- setdiff(indicators, colnames(X_std))
      if (length(missing_ind) > 0L)
        stop("factor_scores: indicator(s) not found in data: ",
             paste(missing_ind, collapse = ", "), call. = FALSE)
      X_std <- X_std[, indicators, drop = FALSE]
    }

    # Rows with any NA get all-NA scores
    complete   <- complete.cases(X_std)
    scores_mat <- matrix(NA_real_, nrow = nrow(X_std), ncol = ncol(Lambda))
    if (any(complete))
      scores_mat[complete, ] <- X_std[complete, , drop = FALSE] %*% W

    df        <- as.data.frame(scores_mat)
    if (ncol(df) != length(fnames))
      warning("factor_scores: Lambda has ", ncol(Lambda), " column(s) but ",
              length(fnames), " factor name(s) found; names may be misaligned.",
              call. = FALSE)
    names(df) <- fnames[seq_len(ncol(Lambda))]
    df
  }

  if (length(raw) == 1L) {
    # Single-group
    return(score_group(raw[[1]]))
  }

  # Multi-group
  grp_labels <- lavaan::lavInspect(x$lavaan_fit, "group.label")
  if (is.null(grp_labels) || length(grp_labels) != length(raw))
    grp_labels <- seq_along(raw)

  dfs <- mapply(function(X_raw, grp) {
    df       <- score_group(X_raw)
    df$group <- grp
    df[, c("group", setdiff(names(df), "group")), drop = FALSE]
  }, raw, grp_labels, SIMPLIFY = FALSE)

  do.call(rbind, dfs)
}

# -- Path 3 helpers: invariance level selection --------------------------------

.fs_select_invariance_level <- function(inv, level, dCFI_cutoff = -0.010) {
  models <- inv$models
  table  <- inv$table

  # Map from bare level keys to the decorated Model column labels used in inv$table
  level_labels <- c(
    configural = "1. Configural",
    weak       = "2. Weak (metric)",
    strong     = "3. Strong (scalar)",
    strict     = "4. Strict",
    varcov     = "5. Latent var/cov",
    means      = "6. Latent means"
  )

  if (is.null(table))
    stop("inv$table is NULL -- the invariance object is incomplete.", call. = FALSE)

  if (all(vapply(models, is.null, logical(1))))
    stop("No invariance models available -- all levels failed to fit.",
         call. = FALSE)

  if (level != "auto") {
    # User-specified level
    fit <- models[[level]]
    if (is.null(fit))
      stop(sprintf("Invariance model '%s' is NULL (failed to fit).", level),
           call. = FALSE)
    # Look up dCFI for the message (configural has no dCFI -- show NA)
    dcfi_val <- if (level == "configural") NA_real_ else {
      lbl <- level_labels[level]
      row <- table[table$Model == lbl, , drop = FALSE]
      if (nrow(row) > 0 && "dCFI" %in% names(row)) row$dCFI[1] else NA_real_
    }
    message(sprintf(
      "factor_scores: using %s invariance model (dCFI = %s)",
      level,
      if (is.na(dcfi_val)) "NA" else sprintf("%.4f", dcfi_val)
    ))
    return(list(level = level, fit = fit))
  }

  # Auto-selection: every fitted level, most constrained first (means, varcov,
  # strict, strong, weak), each judged on its own transition; configural is the
  # unconditional terminal fallback (dCFI is NA -- never tested)
  cutoff      <- dCFI_cutoff
  level_order <- c("means", "varcov", "strict", "strong", "weak")

  # Extract dCFI by matching Model label (case-insensitive)
  .get_dcfi <- function(lv) {
    lbl <- level_labels[lv]
    row <- table[table$Model == lbl, , drop = FALSE]
    if (nrow(row) == 0 || !("dCFI" %in% names(row))) return(NA_real_)
    row$dCFI[1]
  }

  # a level with an inadmissible solution (inv$notes: negative variance, non-PD
  # matrix) is skipped even when its dCFI passes
  selected <- .inv_supported_level(table, models, cutoff, inv$notes)

  dcfi_selected <- if (selected == "configural") NA_real_ else .get_dcfi(selected)

  if (selected == "configural" && any(!vapply(models[intersect(level_order, names(models))], is.null, logical(1)))) {
    warning(sprintf(
      "factor_scores: no invariance level passed the dCFI >= %.4f cutoff. %s",
      cutoff, "Falling back to configural."), call. = FALSE
    )
  }

  message(sprintf(
    "factor_scores: using %s invariance model (dCFI = %s)",
    selected,
    if (is.na(dcfi_selected)) "NA" else sprintf("%.4f", dcfi_selected)
  ))

  list(level = selected, fit = models[[selected]])
}


.fs_invariance <- function(x, method, level, dCFI_cutoff = -0.010) {
  sel <- .fs_select_invariance_level(x, level, dCFI_cutoff = dCFI_cutoff)
  fit <- sel$fit

  # Invariance B-ESEM uses the lavaan-rotated path (no std_rotated_loadings
  # stored), so lavPredict() is the right call. Only the single-group custom
  # WLSMV pipeline needs the manual math.
  if (inherits(fit, "besem_fit_ordered") && !is.null(fit$std_rotated_loadings)) {
    .fs_manual_wlsmv(fit, method = method)
  } else {
    .fs_from_lavaan(fit, method = method)
  }
}

# -- Path 4: Score column alignment via Procrustes / canonical Q ---------------

# Detect column permutation in Q (used by canonical's signed-permutation Q).
# Returns NULL if Q is a non-permutation rotation (e.g. genuine Procrustes
# with off-diagonals large enough that no column is dominated by a single row).
# When non-NULL, the caller relabels score columns to reflect the new factor
# identity so scores$A really contains factor A's scores under the alignment.
.fs_detect_perm <- function(Q, names_old) {
  sources <- vapply(seq_len(ncol(Q)), function(j) which.max(abs(Q[, j])),
                    integer(1L))
  if (length(unique(sources)) != length(sources)) return(NULL)
  is_perm <- all(vapply(seq_len(ncol(Q)), function(j) {
    max(abs(Q[, j])) > 0.95
  }, logical(1L)))
  if (!is_perm) return(NULL)
  names_old[sources]
}

# Apply Q to score matrix M. When Q is a (signed) permutation, returned column
# names reflect the permutation. When Q is a genuine rotation, returned column
# names keep the original ordering — the rotated values land in the original
# factor-named columns. Keeping colnames non-NULL is what lets downstream
# `for (cn in colnames(Mq))` loops actually write the rotated values back.
.fs_apply_to_group <- function(M, Q) {
  Mq <- M %*% Q
  new_names <- .fs_detect_perm(Q, colnames(M))
  colnames(Mq) <- if (!is.null(new_names)) new_names else colnames(M)
  Mq
}

.fs_apply_alignment <- function(scores, x, target, level, dCFI_cutoff = -0.010) {
  if (inherits(x, "besem_fit_ordered") && !is.null(x$std_rotated_loadings))
    stop("factor_scores: align is not yet supported for besem_fit_ordered ",
         "produced by besem_ordered() with custom WLSMV polychoric scoring. ",
         "Compute Q from align_loadings() and rotate the score columns ",
         "manually.", call. = FALSE)

  effective_level <- level
  if (inherits(x, "esem_invariance") && level == "auto") {
    sel <- suppressMessages(.fs_select_invariance_level(x, level,
                                                        dCFI_cutoff = dCFI_cutoff))
    effective_level <- sel$level
  }

  al <- align_loadings(x, target = target, level = effective_level,
                       se_method = "none")
  Qs <- lapply(al, `[[`, "Q")

  has_grp <- "group" %in% names(scores)
  fac_cols <- setdiff(names(scores), "group")
  if (length(fac_cols) == 0L) return(scores)

  if (!has_grp) {
    if (length(Qs) != 1L)
      stop("factor_scores: alignment returned ", length(Qs),
           " rotation matrices but scores are single-group.", call. = FALSE)
    M  <- as.matrix(scores[, fac_cols, drop = FALSE])
    Mq <- .fs_apply_to_group(M, Qs[[1L]])
    # Assign by name -- under permutation, Mq column "X" goes into scores$X
    for (cn in colnames(Mq)) scores[[cn]] <- Mq[, cn]
    return(scores)
  }

  # multi-group: match by group label (both align_loadings and factor_scores
  # use lavInspect(lav, "group.label") so names align)
  group_perms <- vector("list", length(Qs))
  names(group_perms) <- names(Qs)

  for (g_label in names(Qs)) {
    sel <- as.character(scores$group) == g_label
    if (!any(sel)) {
      warning("factor_scores: alignment has no rows for group '",
              g_label, "'; skipping.", call. = FALSE)
      next
    }
    M  <- as.matrix(scores[sel, fac_cols, drop = FALSE])
    Mq <- .fs_apply_to_group(M, Qs[[g_label]])
    # Assign by column name so a permuted Mq writes its "X" column into scores$X
    for (cn in colnames(Mq)) scores[sel, cn] <- Mq[, cn]
    group_perms[[g_label]] <- colnames(Mq)
  }

  perms <- Filter(Negate(is.null), group_perms)
  if (length(perms) > 1L &&
      !all(vapply(perms, identical, logical(1L), perms[[1L]]))) {
    warning("factor_scores: alignment produced different column permutations ",
            "across groups (likely a configural / weak fit with per-group ",
            "loadings of differing Sigma lambda^2 ordering). Each group's ",
            "scores are correctly named, but the per-row meaning of column ",
            "j differs across groups -- interpret with care.", call. = FALSE)
  }
  scores
}
