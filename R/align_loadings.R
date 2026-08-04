# -- Loading Orientation Alignment ---------------------------------------------
#
# Orthogonal bifactor target rotation under NA cross-loading targets is
# non-unique: two implementations (e.g., lavaan's GPA vs Mplus's pairwise
# target) can converge to global optima of identical fit and identical
# column-space, but with loadings related by a near-identity orthogonal
# rotation Q. This file provides:
#
#   - align_loadings():  rotate a fit's standardized loadings to match a
#                        reference (Procrustes), or apply a canonical
#                        deterministic orientation rule.
#   - extract_mplus_loadings(): convenience helper to pull per-group
#                               STDYX matrices from MplusAutomation output.
#
# Substantive quantities (Sigma_lambda^2, omega/ECV, fit indices, df) are
# invariant under orthogonal rotation, so alignment changes only the
# presentation, not the inference.

#' Align ESEM/B-ESEM standardized loadings to a reference
#'
#' Rotates a fitted solution's standardized loadings by an orthogonal Q so
#' they match a reference matrix (Procrustes) or a deterministic canonical
#' orientation. Useful when comparing bifactor solutions across software,
#' or when you want a reproducible orientation across reruns.
#'
#' Multi-group fits and \code{esem_invariance} objects are aligned per group.
#' Item communalities (\eqn{\Sigma\lambda^2}), reliability indices (omega,
#' ECV), and model fit are invariant under orthogonal rotation, so alignment
#' affects only the per-loading partition, not substantive conclusions.
#'
#' @param x An \code{esem_fit}, \code{besem_fit}, \code{esem_invariance},
#'   or raw lavaan S4 fit.
#' @param target One of:
#'   \itemize{
#'     \item a numeric matrix (single-group) or list of matrices (multi-group)
#'       of reference standardized loadings, with rownames = items, colnames
#'       = factors;
#'     \item another fit object of compatible structure;
#'     \item \code{"group1"} -- align all groups to the first group's loadings
#'       within \code{x} (within-fit consistency, no external reference);
#'     \item \code{"canonical"} (default) -- apply deterministic sign + column
#'       order rule (sort columns by \eqn{\Sigma\lambda^2} desc, sign-flip so
#'       the largest \eqn{|\lambda|} per column is positive). No external
#'       reference required.
#'   }
#' @param level Only used for \code{esem_invariance} objects. One of
#'   \code{"configural"} (default), \code{"weak"}, \code{"strong"},
#'   \code{"strict"}.
#' @param se_method \code{"approx"} (default) returns row-wise diagonal
#'   delta-method SEs (ignores within-row covariance between loadings on
#'   different factors). \code{"none"} returns aligned point estimates only.
#'
#' @return An object of class \code{aligned_loadings}: a list with one
#'   element per group, each containing
#'   \describe{
#'     \item{loadings}{aligned standardized loading matrix (items x factors)}
#'     \item{se}{aligned SEs (matching shape) or \code{NULL}}
#'     \item{Q}{the orthogonal rotation matrix applied}
#'     \item{residual_max, residual_mean}{max / mean abs residual vs target
#'       (only when an external reference was supplied)}
#'   }
#'
#' @examples
#' data("HolzingerSwineford1939", package = "lavaan")
#'
#' spec <- specify_model(
#'   Visual  = c("x1", "x2", "x3"),
#'   Textual = c("x4", "x5", "x6"),
#'   Speed   = c("x7", "x8", "x9"),
#'   data  = HolzingerSwineford1939,
#'   group = "school"
#' )
#'
#' \donttest{
#'   inv <- esem_invariance(spec)
#'
#'   # Canonical orientation (no external reference, deterministic)
#'   aligned <- align_loadings(inv, target = "canonical", level = "configural")
#'   print(aligned)
#'
#'   # Within-fit: align all groups to group 1's orientation
#'   aligned <- align_loadings(inv, target = "group1", level = "configural")
#' }
#'
#' \dontrun{
#'   # Align to Mplus output (requires a Mplus .out file + MplusAutomation)
#'   mp  <- MplusAutomation::readModels("besem_inv_configural.out")
#'   tgt <- extract_mplus_loadings(mp)
#'   aligned <- align_loadings(inv, target = tgt, level = "configural")
#' }
#' @export
align_loadings <- function(x,
                           target    = "canonical",
                           level     = "configural",
                           se_method = c("approx", "none")) {
  se_method <- match.arg(se_method)

  std    <- .extract_loadings_with_se(x, level = level)
  groups <- names(std)
  tgt    <- .resolve_target(target, std)
  method <- attr(tgt, "method")

  out <- vector("list", length(groups))
  names(out) <- groups

  for (g in groups) {
    L_R  <- std[[g]]$loadings
    SE_R <- std[[g]]$se

    Q <- if (method == "canonical") {
      .canonical_rotation(L_R)
    } else {
      L_T <- .check_target_compat(tgt[[g]], L_R, group = g)
      .procrustes_rotation(L_R, L_T)
    }

    L_aligned  <- L_R %*% Q
    dimnames(L_aligned) <- dimnames(L_R)
    SE_aligned <- if (se_method == "approx") .se_rotate_approx(SE_R, Q) else NULL
    if (!is.null(SE_aligned)) dimnames(SE_aligned) <- dimnames(SE_R)

    res <- list(loadings = L_aligned, se = SE_aligned, Q = Q)
    if (method == "procrustes" && !is.null(tgt[[g]])) {
      d <- abs(L_aligned - tgt[[g]])
      res$residual_max  <- max(d, na.rm = TRUE)
      res$residual_mean <- mean(d, na.rm = TRUE)
    }
    out[[g]] <- res
  }

  attr(out, "method")    <- method
  attr(out, "se_method") <- se_method
  attr(out, "level")     <- if (inherits(x, "esem_invariance")) level else NA_character_
  class(out) <- "aligned_loadings"
  out
}

# Orthogonal Procrustes: argmin_Q ||L_R Q - L_T||_F  s.t. Q'Q = I.
# Closed form via SVD: Q = U V', where SVD(L_R' L_T) = U D V'.
.procrustes_rotation <- function(L_R, L_T) {
  M  <- crossprod(L_R, L_T)
  sv <- svd(M)
  sv$u %*% t(sv$v)
}

# Deterministic canonical orientation:
#   1. permute columns by Sigma_lambda^2 descending
#   2. sign-flip so the largest |lambda| in each column is positive
# Produces a unique orientation for any rotated solution (up to ties).
.canonical_rotation <- function(L_R) {
  M    <- ncol(L_R)
  ord  <- order(colSums(L_R^2, na.rm = TRUE), decreasing = TRUE)
  P    <- diag(M)[, ord, drop = FALSE]
  Lp   <- L_R %*% P
  sgn  <- vapply(seq_len(M), function(j) {
    col <- Lp[, j]
    idx <- which.max(abs(col))
    s   <- sign(col[idx])
    if (s == 0) 1 else s
  }, numeric(1L))
  P %*% diag(sgn, nrow = M)
}

# Approximate delta-method SE propagation under orthogonal rotation:
# (Lambda Q)_ij = sum_k Lambda_ik Q_kj
# Var((Lambda Q)_ij) ~= sum_k Q_kj^2 * Var(Lambda_ik)
#                       + 2 sum_{k<l} Q_kj Q_lj Cov(Lambda_ik, Lambda_il)
# We retain only the diagonal (first) term. The cross-covariance term is
# bounded by |Q_kj| |Q_lj| sqrt(Var_ik Var_il); for Q near-identity (off-
# diagonals < 0.15 in practice) the omitted contribution is small.
.se_rotate_approx <- function(SE_R, Q) {
  sqrt(SE_R^2 %*% Q^2)
}

# Resolve user-facing target argument to a list of per-group reference
# matrices, plus a "method" attribute = "procrustes" or "canonical".
.resolve_target <- function(target, std) {
  groups <- names(std)

  if (is.character(target) && length(target) == 1L) {
    if (target == "canonical") {
      out <- vector("list", length(groups))
      names(out) <- groups
      attr(out, "method") <- "canonical"
      return(out)
    }
    if (target == "group1") {
      L1  <- std[[groups[1L]]]$loadings
      out <- replicate(length(groups), L1, simplify = FALSE)
      names(out) <- groups
      attr(out, "method") <- "procrustes"
      return(out)
    }
    stop("Unknown character target: '", target,
         "'. Use 'canonical', 'group1', or supply a matrix / list / fit.",
         call. = FALSE)
  }

  if (is.matrix(target)) {
    if (length(groups) != 1L)
      stop("For multi-group fit, target must be a list of matrices ",
           "(one per group), not a single matrix.", call. = FALSE)
    out <- list(target)
    names(out) <- groups
    attr(out, "method") <- "procrustes"
    return(out)
  }

  if (is.list(target) && !inherits(target, c("esem_fit", "esem_invariance"))) {
    if (length(target) != length(groups))
      stop("target list length (", length(target),
           ") must match number of groups (", length(groups), ")",
           call. = FALSE)
    if (is.null(names(target))) names(target) <- groups
    attr(target, "method") <- "procrustes"
    return(target)
  }

  # Treat as another fit -- extract its loadings.
  ext <- .extract_loadings_with_se(target)
  if (length(ext) != length(groups))
    stop("Reference fit has ", length(ext), " group(s); current fit has ",
         length(groups), ".", call. = FALSE)
  out <- lapply(ext, `[[`, "loadings")
  names(out) <- groups
  attr(out, "method") <- "procrustes"
  out
}

# Verify a candidate target matrix matches Lambda_R's items / factors.
# Reorders rows / columns by name when possible. Errors with a useful
# message if the structure is incompatible.
.check_target_compat <- function(L_T, L_R, group) {
  if (is.null(L_T))
    stop("No target supplied for group '", group, "'.", call. = FALSE)
  if (!identical(dim(L_T), dim(L_R))) {
    stop("Group '", group, "': target dim is ",
         paste(dim(L_T), collapse = " x "),
         " but fit Lambda is ", paste(dim(L_R), collapse = " x "), ".",
         call. = FALSE)
  }
  if (!is.null(rownames(L_T)) && !is.null(rownames(L_R))) {
    if (!setequal(rownames(L_T), rownames(L_R)))
      stop("Group '", group, "': target rownames don't match fit items.",
           call. = FALSE)
    L_T <- L_T[rownames(L_R), , drop = FALSE]
  }
  if (!is.null(colnames(L_T)) && !is.null(colnames(L_R))) {
    if (!setequal(colnames(L_T), colnames(L_R)))
      stop("Group '", group, "': target colnames don't match fit factors.",
           call. = FALSE)
    L_T <- L_T[, colnames(L_R), drop = FALSE]
  }
  L_T
}

# Extract per-group standardized loadings + SEs from any supported fit type.
# Returns a named list (one element per group), each a list with:
#   loadings: items x factors matrix
#   se:       items x factors matrix (zeros where standardizedSolution
#             reports no SE -- e.g., fixed parameters)
.extract_loadings_with_se <- function(x, level = "configural") {
  if (inherits(x, "esem_invariance")) {
    fit <- x$models[[level]]
    if (is.null(fit))
      stop("Invariance object has no '", level, "' level fit.", call. = FALSE)
    lav <- if (isS4(fit)) fit else fit$lavaan_fit
  } else if (inherits(x, c("esem_fit", "besem_fit"))) {
    lav <- x$lavaan_fit
  } else if (isS4(x) && methods::is(x, "lavaan")) {
    lav <- x
  } else {
    stop("Unsupported class for align_loadings: ",
         paste(class(x), collapse = "/"), call. = FALSE)
  }

  ss <- lavaan::standardizedSolution(lav, type = "std.all", se = TRUE)
  ss <- ss[ss$op == "=~", , drop = FALSE]

  has_grp <- "group" %in% names(ss) && length(unique(ss$group)) > 1L
  groups  <- if (has_grp) sort(unique(ss$group)) else 1L

  grp_labels <- if (has_grp) {
    lab <- tryCatch(lavaan::lavInspect(lav, "group.label"),
                    error = function(e) NULL)
    if (is.null(lab) || length(lab) != length(groups)) as.character(groups)
    else lab
  } else "1"

  factors <- unique(ss$lhs)
  items   <- unique(ss$rhs)

  out <- vector("list", length(groups))
  names(out) <- grp_labels

  for (k in seq_along(groups)) {
    g  <- groups[k]
    sg <- if (has_grp) ss[ss$group == g, , drop = FALSE] else ss
    L  <- matrix(0, nrow = length(items), ncol = length(factors),
                 dimnames = list(items, factors))
    SE <- matrix(0, nrow = length(items), ncol = length(factors),
                 dimnames = list(items, factors))
    for (i in seq_len(nrow(sg))) {
      L[sg$rhs[i],  sg$lhs[i]] <- sg$est.std[i]
      se_val <- sg$se[i]
      if (!is.na(se_val)) SE[sg$rhs[i], sg$lhs[i]] <- se_val
    }
    out[[k]] <- list(loadings = L, se = SE)
  }
  out
}

#' Extract standardized loadings from MplusAutomation output
#'
#' Convenience helper for use as the \code{target} of
#' \code{align_loadings()}: pulls the STDYX-standardized loadings from a
#' \code{MplusAutomation::readModels()} result and returns one matrix per
#' group, in items-by-factors form with names matching common conventions.
#'
#' @param mp_out An object returned by \code{MplusAutomation::readModels()}.
#' @param standardized Which Mplus section to read. Default
#'   \code{"stdyx.standardized"}; set to \code{"unstandardized"} to read raw
#'   point estimates (note: comparison with lavaan's
#'   \code{parameterEstimates()} is not meaningful for ESEM \code{efa()}
#'   blocks because lavaan reports the unrotated optimization basis).
#' @param items Optional character vector. If supplied, the returned
#'   matrices use this row order; otherwise the order from Mplus is used.
#' @param factors Optional character vector. If supplied, the returned
#'   matrices use this column order.
#'
#' @return A named list of matrices, one per group, each with rownames
#'   (items) and colnames (factors).
#'
#' @export
extract_mplus_loadings <- function(mp_out,
                                   standardized = "stdyx.standardized",
                                   items   = NULL,
                                   factors = NULL) {
  ps <- mp_out$parameters[[standardized]]
  if (is.null(ps))
    stop("Mplus output does not contain section '", standardized, "'.",
         call. = FALSE)
  ld <- ps[grepl("BY$", ps$paramHeader), , drop = FALSE]
  if (nrow(ld) == 0L)
    stop("No BY rows found in Mplus output section '", standardized, "'.",
         call. = FALSE)
  ld$factor <- sub(".BY$", "", ld$paramHeader)
  ld$item   <- toupper(ld$param)

  has_grp <- "Group" %in% names(ld) && length(unique(ld$Group)) > 1L
  groups  <- if (has_grp) unique(ld$Group) else "1"

  if (is.null(items))   items   <- unique(ld$item)
  if (is.null(factors)) factors <- unique(ld$factor)
  items_u <- toupper(items)

  out <- vector("list", length(groups))
  names(out) <- groups
  for (g in groups) {
    sg <- if (has_grp) ld[ld$Group == g, , drop = FALSE] else ld
    M  <- matrix(0, nrow = length(items_u), ncol = length(factors),
                 dimnames = list(items_u, factors))
    for (i in seq_len(nrow(sg))) {
      it <- sg$item[i]
      fc <- sg$factor[i]
      if (it %in% items_u && fc %in% factors)
        M[it, fc] <- sg$est[i]
    }
    out[[g]] <- M
  }
  out
}

#' Print method for aligned_loadings
#' @param x An \code{aligned_loadings} object.
#' @param digits Integer; number of decimal places to print. Default 3.
#' @param ... Ignored.
#' @return Invisibly returns \code{x}; called for the side effect of printing
#'   the aligned loading matrices.
#' @export
#' @keywords internal
print.aligned_loadings <- function(x, digits = 3, ...) {
  meth <- attr(x, "method")
  lvl  <- attr(x, "level")
  cat("aligned_loadings (method = '", meth, "'", sep = "")
  if (!is.na(lvl)) cat(", level = '", lvl, "'", sep = "")
  cat(", groups = ", length(x), ")\n", sep = "")
  for (g in names(x)) {
    cat("\n-- group: ", g, " --\n", sep = "")
    print(round(x[[g]]$loadings, digits))
    if (!is.null(x[[g]]$residual_max))
      cat(sprintf("residual: max = %.4f, mean = %.4f\n",
                  x[[g]]$residual_max, x[[g]]$residual_mean))
  }
  invisible(x)
}
