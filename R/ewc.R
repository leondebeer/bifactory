# R/ewc.R ---------------------------------------------------------------------
# ESEM-within-CFA (EWC): converts a fitted ESEM solution into explicit lavaan
# CFA syntax following Marsh et al. (2014) / esemComp approach.
#
# Algorithm (from the reference .inp files in within/):
#   1. Extract UNSTANDARDIZED loadings from the fitted ESEM (parameterEstimates)
#   2. Select one referent per factor -- item with highest |unstd loading| on
#      its primary factor
#   3. Generate lavaan syntax:
#        - Referent item's cross-loadings on OTHER factors  -> FIXED: value*item
#        - var_fixed=FALSE also fixes referent's own primary loading -> value*item
#        - Everything else                                  -> FREE:  start(v)*item
#   4. Fit via lavaan::cfa() -- no rotation needed, identification is algebraic
#
# Two identification modes (matching the two Mplus .inp reference files):
#   var_fixed = TRUE  (default) : factor variances fixed to 1 (std.lv=TRUE)
#                                  referent's own primary loading stays free
#   var_fixed = FALSE           : factor variances freely estimated;
#                                  referent's own primary loading is FIXED
# Both modes apply to the B-ESEM-within-CFA path too (G and every specific
# factor). Identification comes entirely from the generated syntax: lavaan's
# auto.fix.first stays OFF, because with it on lavaan also pins the first
# indicator of any factor whose referent is not that first item (at its
# start() value), and the EWC stops being a re-expression of the ESEM.
#
# These functions are self-contained and do not modify any pipeline file.
# They use public lavaan functions and .rmsea_ci() from pipeline.R (same namespace).
# -----------------------------------------------------------------------------


# == 1. Referent selection ======================================================

#' Auto-Select EWC Referent Items
#'
#' For each factor, selects the primary indicator with the largest absolute
#' unstandardised loading from the ESEM solution.  This is the item that will
#' anchor that factor's identification in the EWC model.
#'
#' @param esem_fit An \code{esem_fit} object (e.g. \code{results$fit_esem}
#'   from \code{\link{run_comparison}}).
#' @param spec An \code{esem_spec} from \code{\link{specify_model}}.
#'
#' @return A named character vector: names are factor names, values are the
#'   selected referent item names (in the original case from \code{spec}).
#'
#' @details
#' Referents are chosen solely by maximum |loading| on the primary factor.
#' You can override automatic selection by supplying a custom named vector
#' to the \code{referents} argument of \code{\link{ewc_syntax}} or
#' \code{\link{fit_ewc}}.
#'
#' @seealso \code{\link{ewc_syntax}}, \code{\link{fit_ewc}}
#' @export
find_ewc_referents <- function(esem_fit, spec) {
  if (!inherits(esem_fit, "esem_fit"))
    stop("`esem_fit` must be an esem_fit object.", call. = FALSE)
  if (!inherits(spec, "esem_spec"))
    stop("`spec` must be an esem_spec from specify_model().", call. = FALSE)

  # B-ESEM: return specific-factor referents + one additional G referent.
  # G referent is chosen among items NOT used as a specific-factor referent
  # (follows Mplus convention for BIFACTOR-ESEM-WITHIN-CFA).
  if (inherits(esem_fit, "besem_fit"))
    return(.besem_ewc_referents(esem_fit, spec))

  ld <- .ewc_unstd_loadings(esem_fit)

  factors <- names(spec$factors)
  refs    <- setNames(character(length(factors)), factors)

  for (f in factors) {
    prim    <- tolower(spec$factors[[f]])
    f_rows  <- ld[tolower(ld$lhs) == tolower(f) &
                  tolower(ld$rhs) %in% prim, , drop = FALSE]
    if (nrow(f_rows) == 0L)
      stop("No loadings found for factor '", f,
           "'. Check that esem_fit and spec are compatible.", call. = FALSE)
    best_idx <- which.max(abs(f_rows$est))
    # Return item name in original case (from spec, not lavaan's stored case)
    matched  <- spec$factors[[f]][tolower(spec$factors[[f]]) == tolower(f_rows$rhs[best_idx])]
    refs[f]  <- if (length(matched) > 0L) matched[1L] else f_rows$rhs[best_idx]
  }
  refs
}


# == 2. Syntax generator =======================================================

#' Generate ESEM-within-CFA Lavaan Syntax
#'
#' Converts a fitted ESEM solution into explicit lavaan CFA syntax following
#' the ESEM-within-CFA approach (Marsh et al. 2014).  The syntax uses
#' unstandardised ESEM loadings as starting values (\code{start(v)*item}) or
#' fixed values (\code{v*item}) depending on the referent scheme.
#'
#' @param esem_fit An \code{esem_fit} object from \code{\link{run_comparison}}
#'   (\code{results$fit_esem}) or \code{\link{esem}} / \code{\link{esem_ordered}}.
#' @param spec An \code{esem_spec} from \code{\link{specify_model}}.
#' @param referents Named character vector of referent items, one per factor.
#'   Names = factor names; values = item names.  If \code{NULL} (default),
#'   referents are selected automatically via \code{\link{find_ewc_referents}}.
#' @param var_fixed Logical.  Identification mode:
#'   \describe{
#'     \item{\code{TRUE} (default)}{Factor variances fixed to 1 (\code{std.lv = TRUE}).
#'       Only referent cross-loadings on other factors are fixed.
#'       Mirrors Mplus \code{EX\@1; MD\@1; CI\@1;}.}
#'     \item{\code{FALSE}}{Factor variances freely estimated.
#'       Referent's own primary loading is \emph{also} fixed for identification.
#'       Mirrors Mplus \code{EX*; MD*; CI*;}.}
#'   }
#'   For a B-ESEM source the same two schemes apply to G and to every
#'   specific factor; all factor covariances stay fixed at zero in both.
#'
#' @return A character string of lavaan model syntax.  Use \code{cat()} to
#'   inspect, or pass directly to \code{\link{fit_ewc}} via
#'   \code{custom_syntax}.
#'
#' @references
#' Marsh, H. W., Morin, A. J. S., Parker, P. D., & Kaur, G. (2014).
#' Exploratory structural equation modeling. \emph{Annual Review of Clinical
#' Psychology, 10}, 85-110.
#'
#' @seealso \code{\link{find_ewc_referents}}, \code{\link{fit_ewc}}
#' @export
ewc_syntax <- function(esem_fit, spec,
                       referents = NULL,
                       var_fixed = TRUE) {
  if (!inherits(esem_fit, "esem_fit"))
    stop("`esem_fit` must be an esem_fit object.", call. = FALSE)
  if (!inherits(spec, "esem_spec"))
    stop("`spec` must be an esem_spec from specify_model().", call. = FALSE)

  # B-ESEM path: route to the bifactor EWC generator (different convention --
  # G referent + orthogonality constraints).
  if (inherits(esem_fit, "besem_fit"))
    return(.besem_ewc_syntax(esem_fit, spec, referents = referents,
                             var_fixed = var_fixed))

  if (is.null(referents))
    referents <- find_ewc_referents(esem_fit, spec)

  factors <- names(spec$factors)

  # Validate referents
  for (f in factors) {
    if (is.null(referents[f]) || is.na(referents[f]) || referents[f] == "")
      stop("Missing referent for factor '", f, "'.", call. = FALSE)
    if (!tolower(referents[f]) %in% tolower(spec$factors[[f]]))
      stop("Referent '", referents[f], "' is not a primary item of factor '", f, "'.",
           call. = FALSE)
  }

  # Build loading lookup: factor (lower)  x  item (lower) -> unstd est
  ld  <- .ewc_unstd_loadings(esem_fit)
  lkp <- function(f, it) {
    r <- ld[tolower(ld$lhs) == tolower(f) & tolower(ld$rhs) == tolower(it), , drop = FALSE]
    if (nrow(r) == 0L) return(NA_real_)
    r$est[1L]
  }

  # referent_owner: lower(item) -> factor that owns it as referent
  ref_owner <- setNames(factors, tolower(unname(referents)))

  # Generate per-factor loading lines
  load_lines <- vapply(factors, function(f) {
    ref_lower <- tolower(referents[f])

    # Primary items first, cross-loading items after (matches Mplus convention)
    prim_items  <- spec$factors[[f]]
    cross_items <- setdiff(spec$all_items, prim_items)
    ordered_items <- c(prim_items, cross_items)

    terms <- vapply(ordered_items, function(it) {
      val      <- lkp(f, it)
      it_lower <- tolower(it)

      # Decide: FIXED or FREE?
      is_fixed <- FALSE
      if (it_lower %in% names(ref_owner)) {
        owner <- ref_owner[it_lower]
        if (owner != f) {
          is_fixed <- TRUE          # referent of ANOTHER factor -> FIXED (cross-loading)
        } else if (!var_fixed) {
          is_fixed <- TRUE          # own referent, var_free mode -> FIXED (primary loading)
        }
      }

      if (!is.na(val)) {
        if (is_fixed)
          sprintf("%.3f*%s", val, it)           # FIXED  (no start())
        else
          sprintf("start(%.3f)*%s", val, it)    # FREE with starting value
      } else {
        if (is_fixed) it else sprintf("start(0)*%s", it)
      }
    }, character(1L))

    # Wrap at 4 items per line for readability
    chunks   <- split(terms, ceiling(seq_along(terms) / 4L))
    rhs_rows <- vapply(chunks, paste, character(1L), collapse = " + ")
    rhs      <- paste(rhs_rows, collapse = " +\n       ")
    paste0(f, " =~ ", rhs)
  }, character(1L))

  header <- paste0(
    "# ESEM-within-CFA (EWC) lavaan syntax\n",
    "# Marsh et al. (2014) / esemComp approach\n",
    "# var_fixed = ", var_fixed, "\n",
    "# Referents: ",
    paste(paste0(names(referents), " = '", referents, "'"), collapse = ", "), "\n",
    "#   value*item      = FIXED to ESEM unstd loading\n",
    "#   start(v)*item   = FREE, ESEM unstd loading as starting value\n",
    "# Generated by bifactory::ewc_syntax() -- edit freely\n"
  )

  if (var_fixed) {
    id_note <- "# Identification: std.lv = TRUE (factor variances = 1) in fit_ewc()"
  } else {
    id_note <- "# Identification: factor variances free, referent primary loadings fixed above"
  }

  paste(c(header, id_note, "", load_lines), collapse = "\n")
}


# == 3. Model fitter ===========================================================

#' Fit an ESEM-within-CFA Model
#'
#' Generates EWC lavaan syntax from a fitted ESEM solution and estimates it as
#' a standard \code{lavaan::cfa()} model -- no rotation required.
#'
#' @param esem_fit An \code{esem_fit} object (\code{results$fit_esem}).
#' @param spec An \code{esem_spec} from \code{\link{specify_model}}.
#' @param referents Named character vector of referent items or \code{NULL}
#'   (auto-selected).  See \code{\link{find_ewc_referents}}.
#' @param var_fixed Logical.  Identification mode.  Default \code{TRUE}
#'   (factor variances = 1).  See \code{\link{ewc_syntax}} for details.
#' @param missing Character.  Missing data handling.  Defaults to
#'   \code{"pairwise"} for ordered, \code{"listwise"} for continuous.
#' @param custom_syntax Character or \code{NULL}.  Supply a hand-edited
#'   syntax string (from \code{\link{ewc_syntax}}) instead of auto-generating.
#'   When non-\code{NULL}, \code{referents} and \code{var_fixed} still control
#'   the \code{lavaan::cfa()} options (\code{std.lv}).
#' @param ... Additional arguments forwarded to \code{lavaan::cfa()}.
#'
#' @return An object of class \code{"ewc_fit"}:
#' \describe{
#'   \item{\code{lavaan_fit}}{lavaan fit object; all lavaan generics work on it.}
#'   \item{\code{syntax}}{lavaan model string used.}
#'   \item{\code{referents}}{Named vector of referent items.}
#'   \item{\code{var_fixed}}{Identification mode used.}
#'   \item{\code{estimator}}{Estimator (\code{"MLR"} for continuous, \code{"DWLS"} for ordered).}
#'   \item{\code{spec}}{The spec object.}
#' }
#'
#' @seealso \code{\link{ewc_syntax}}, \code{\link{find_ewc_referents}},
#'   \code{\link{compare_ewc}}
#' @export
fit_ewc <- function(esem_fit,
                    spec,
                    referents     = NULL,
                    var_fixed     = TRUE,
                    missing       = NULL,
                    custom_syntax = NULL,
                    ...) {
  if (!inherits(esem_fit, "esem_fit"))
    stop("`esem_fit` must be an esem_fit object.", call. = FALSE)
  if (!inherits(spec, "esem_spec"))
    stop("`spec` must be an esem_spec from specify_model().", call. = FALSE)

  is_besem <- inherits(esem_fit, "besem_fit")

  if (is.null(referents))
    referents <- find_ewc_referents(esem_fit, spec)

  is_ordered <- isTRUE(spec$ordered) ||
    (is.character(spec$ordered) && length(spec$ordered) > 0L)
  estimator  <- if (is_ordered) "WLSMV" else "MLR"
  if (is.null(missing)) missing <- if (is_ordered) "pairwise" else "listwise"

  syntax <- if (!is.null(custom_syntax)) {
    custom_syntax
  } else {
    ewc_syntax(esem_fit, spec, referents = referents, var_fixed = var_fixed)
  }

  # Identification comes entirely from the generated syntax (referent loadings
  # pinned; variances fixed to 1 through std.lv when var_fixed). auto.fix.first
  # must stay off: with it on, lavaan also pins the first indicator of every
  # factor whose referent is not that first item, at its start() value, and
  # the EWC gains df it should not have.
  cfa_args <- list(
    model          = syntax,
    data           = spec$data,
    std.lv         = isTRUE(var_fixed),
    auto.fix.first = FALSE,
    estimator      = estimator,
    missing        = missing
  )

  # B-EWC writes the orthogonality constraints (F ~~ 0*F2) directly into the
  # model syntax; we do NOT also pass orthogonal = TRUE to lavaan because it
  # behaves unexpectedly in multi-group contexts (see CLAUDE.md). The explicit
  # constraints are enough for identification in single-group B-EWC.

  if (is_ordered) {
    cfa_args$ordered          <- spec$all_items
    # Pinned values must live on the scale the EWC is fitted on. A custom-DWLS
    # B-ESEM stores rotated loadings on the standardised (delta) scale; a
    # lavaan-path ESEM stores estimates on its own parameterization (theta by
    # default in esem_ordered()), so follow the source.
    cfa_args$parameterization <- if (is_besem) "delta" else
      .ewc_source_parameterization(esem_fit)
  }

  dots <- list(...)
  if (length(dots)) cfa_args <- c(cfa_args, dots)

  lav_fit <- tryCatch(
    do.call(lavaan::cfa, cfa_args),
    error = function(e)
      stop("lavaan::cfa() failed in fit_ewc():\n  ", conditionMessage(e),
           call. = FALSE)
  )

  # The B-EWC is a CFA with k(k-1)/2 more constraints than the orthogonal
  # B-ESEM (all referent cross-loadings pinned, not just the k(k-1)/2 the
  # rotation needs), so its df is larger and its fit is its own -- reported
  # from this lavaan fit, never copied from the source.
  structure(
    list(
      lavaan_fit  = lav_fit,
      syntax      = syntax,
      referents   = referents,
      var_fixed   = isTRUE(var_fixed),
      estimator   = estimator,
      spec        = spec,
      besem       = is_besem
    ),
    class = "ewc_fit"
  )
}


# == 4. Comparison =============================================================

#' Compare Pipeline Results with an EWC Model
#'
#' Appends fit indices from a \code{\link{fit_ewc}} result to the pipeline
#' comparison table, producing a unified data frame with CFA, ESEM, B-ESEM,
#' and EWC columns side by side.
#'
#' @param results An \code{esem_comparison_pipeline} from
#'   \code{\link{run_comparison}}.
#' @param ewc An \code{ewc_fit} from \code{\link{fit_ewc}}.
#'
#' @return An object of class \code{c("ewc_comparison","data.frame")}.
#'   Print with \code{print()} for a formatted table.
#'
#' @seealso \code{\link{fit_ewc}}
#' @export
compare_ewc <- function(results, ewc) {
  if (!inherits(results, "esem_comparison_pipeline"))
    stop("`results` must be an esem_comparison_pipeline from run_comparison().",
         call. = FALSE)
  if (!inherits(ewc, "ewc_fit"))
    stop("`ewc` must be an ewc_fit from fit_ewc().", call. = FALSE)

  base  <- results$comparison_table
  ewc_v <- .ewc_fi(ewc$lavaan_fit, ewc$estimator)

  idx_labels <- c("CFI", "TLI", "RMSEA",
                  "RMSEA [L90%CI]", "RMSEA [U90%CI]",
                  "SRMR", "X2", "df", "p",
                  "AIC", "BIC", "SABIC")

  ewc_ordered <- unname(ewc_v[match(base[[1L]], idx_labels)])

  id_tag    <- if (ewc$var_fixed) "varfix" else "varfree"
  col_label <- paste0("EWC_", id_tag)

  out <- cbind(base,
               setNames(data.frame(ewc_ordered, stringsAsFactors = FALSE),
                        col_label))
  attr(out, "ewc_var_fixed")  <- ewc$var_fixed
  attr(out, "ewc_referents")  <- ewc$referents
  attr(out, "ewc_estimator")  <- ewc$estimator
  class(out) <- c("ewc_comparison", "data.frame")
  out
}


# == 5. S3 print methods =======================================================

#' Print Method for ewc_fit
#'
#' Compact fit summary (CFI/TLI/RMSEA/SRMR and \eqn{\chi^2}) for an
#' ESEM-within-CFA fit.
#'
#' @param x An \code{ewc_fit} object from \code{\link{fit_ewc}}.
#' @param ... Ignored.
#'
#' @return Invisibly returns \code{x}.
#' @export
print.ewc_fit <- function(x, ...) {
  obj_name <- deparse(substitute(x))
  if (length(obj_name) != 1L || !nzchar(obj_name) || grepl("[^A-Za-z0-9._]", obj_name))
    obj_name <- "x"
  fi <- .ewc_fi(x$lavaan_fit, x$estimator)

  id_label <- if (x$var_fixed) "factor variances fixed to 1"
              else              "factor variances free (referent primary loadings fixed)"

  variant_label <- if (isTRUE(x$besem)) "B-ESEM-within-CFA"
                   else                 "ESEM-within-CFA"

  cat(sprintf("\n== EWC fit: %s =============================\n", variant_label))
  cat(" Identification  :", id_label, "\n")
  cat(" Estimator       :", x$estimator, "\n")
  cat(" Referents       :",
      paste(paste0(names(x$referents), "=", x$referents), collapse = ", "), "\n\n")

  if (!is.na(fi[["X2"]])) {
    cat(sprintf("  CFI = %.3f   TLI = %.3f   RMSEA = %.3f   SRMR = %.3f\n",
                fi[["CFI"]], fi[["TLI"]], fi[["RMSEA"]], fi[["SRMR"]]))
    cat(sprintf("  X2(%g) = %.3f\n", fi[["df"]], fi[["X2"]]))
  }

  cat("\n Use summary(", obj_name, ") for full output (adds colour-coded loadings table)\n", sep = "")
  cat(" Use parameters(", obj_name, ") for just the loadings table\n", sep = "")
  cat(" Use cat(", obj_name, "$syntax) to view/copy the model string\n", sep = "")
  invisible(x)
}


#' Summary for ESEM-within-CFA Fits
#'
#' Forwards to \code{lavaan::lavaan-class} \code{summary} on the underlying
#' lavaan fit, so you do not need to \code{library(lavaan)} separately.
#'
#' @param object An \code{ewc_fit} object from \code{\link{fit_ewc}}.
#' @param fit.measures Logical. Include fit indices. Default \code{TRUE}.
#' @param standardized Logical. Include standardised estimates. Default \code{TRUE}.
#' @param show_loadings Logical. After the lavaan summary, also print the
#'   \code{\link{parameters}} table with primary loadings colour-highlighted.
#'   Default \code{TRUE}.
#' @param ... Additional arguments passed to lavaan's summary method.
#'
#' @return Invisibly returns the lavaan summary object.
#' @export
summary.ewc_fit <- function(object,
                             fit.measures   = TRUE,
                             standardized   = TRUE,
                             show_loadings  = TRUE,
                             ...) {
  lav_summary <- methods::getMethod("summary", "lavaan")
  out <- lav_summary(object$lavaan_fit,
                     fit.measures = fit.measures,
                     standardized = standardized,
                     ...)
  # lavaan's S4 summary returns a summary object that prints on auto-display.
  # When called from inside a function, auto-print is suppressed -- force it.
  print(out)
  if (isTRUE(show_loadings))
    print(parameters(object))
  invisible(out)
}


#' Print Method for ewc_comparison
#'
#' Displays the fit-index comparison table returned by \code{\link{compare_ewc}}.
#'
#' @param x An \code{ewc_comparison} object.
#' @param ... Ignored.
#'
#' @return Invisibly returns \code{x}.
#' @export
print.ewc_comparison <- function(x, ...) {
  cat("\n== EWC vs Pipeline -- Fit Comparison =======================\n\n")

  val_cols <- setdiff(colnames(x), "Index")
  idx_col  <- x[["Index"]]

  # Numeric display
  disp <- as.data.frame(lapply(x, as.character), stringsAsFactors = FALSE)
  for (col in val_cols) {
    vals <- suppressWarnings(as.numeric(x[[col]]))
    disp[[col]] <- ifelse(is.na(vals), "-", sprintf("%.3f", vals))
  }

  # Merge RMSEA + CI into one row
  rmsea_r <- which(idx_col == "RMSEA")
  li       <- which(idx_col == "RMSEA [L90%CI]")
  hi       <- which(idx_col == "RMSEA [U90%CI]")
  if (length(rmsea_r) == 1L && length(li) == 1L && length(hi) == 1L) {
    disp$Index[rmsea_r] <- "RMSEA [90% CI]"
    for (col in val_cols) {
      rv <- suppressWarnings(as.numeric(x[[col]][rmsea_r]))
      lv <- suppressWarnings(as.numeric(x[[col]][li]))
      hv <- suppressWarnings(as.numeric(x[[col]][hi]))
      if (!is.na(rv) && !is.na(lv) && !is.na(hv)) {
        disp[[col]][rmsea_r] <- sprintf("%.3f [%s, %s]", rv,
          sub("^0", "", sprintf("%.3f", lv)),
          sub("^0", "", sprintf("%.3f", hv)))
      }
    }
    disp <- disp[-c(li, hi), , drop = FALSE]
  }

  print(disp, row.names = FALSE)

  refs <- attr(x, "ewc_referents")
  cat("\n Referents: ",
      paste(paste0(names(refs), "=", refs), collapse = ", "), "\n",
      " Identification: ",
      if (attr(x, "ewc_var_fixed")) "factor variances fixed to 1"
      else "factor variances free", "\n",
      sep = "")
  invisible(x)
}


# == 6. Internal helpers =======================================================

# Extract unstandardised loading estimates from any esem_fit.
# Returns a data frame with columns lhs (factor), rhs (item), est.
.ewc_unstd_loadings <- function(esem_fit) {
  lav <- esem_fit$lavaan_fit
  pe  <- tryCatch(
    lavaan::parameterEstimates(lav),
    error = function(e)
      stop("Could not extract parameter estimates from esem_fit. ",
           "Is lavaan_fit a valid lavaan object?", call. = FALSE)
  )
  ld <- pe[pe$op == "=~", c("lhs", "rhs", "est"), drop = FALSE]
  if (nrow(ld) == 0L)
    stop("No loading parameters (op == '=~') found in esem_fit.", call. = FALSE)
  ld
}

# Extract 12 fit values (matching the comparison_table row order) through the
# shared helper, so the EWC column agrees with the pipeline columns on scaled
# (MLR / WLSMV) statistics and on the Mplus SRMR denominator.
.ewc_fi <- function(lav, estimator) {
  labels <- c("CFI","TLI","RMSEA","RMSEA [L90%CI]","RMSEA [U90%CI]",
              "SRMR","X2","df","p","AIC","BIC","SABIC")
  fi <- .fit_indices_lav(lav)
  if (is.null(fi)) return(setNames(rep(NA_real_, 12L), labels))

  lo <- fi[["rmsea.ci.lower"]]; hi <- fi[["rmsea.ci.upper"]]
  if (is.na(lo) || is.na(hi)) {
    n  <- tryCatch(sum(lavaan::lavInspect(lav, "nobs")),
                   error = function(e) NA_integer_)
    ci <- .rmsea_ci(fi[["chisq"]], fi[["df"]], n)
    lo <- ci[1L]; hi <- ci[2L]
  }
  setNames(round(c(fi[["cfi"]], fi[["tli"]], fi[["rmsea"]], lo, hi, fi[["srmr"]],
                   fi[["chisq"]], fi[["df"]], fi[["pvalue"]],
                   fi[["aic"]], fi[["bic"]], fi[["bic2"]]), 3L), labels)
}


# == 7. B-ESEM-within-CFA helpers =============================================
# Follows Mplus convention (see BIFACTOR-ESEM-WITHIN-CFA reference):
#   - One referent per specific factor (item with highest |specific loading|)
#   - One additional referent for G (highest |G loading| among items not already
#     used as a specific-factor referent)
#   - Every referent has ALL its other-factor loadings FIXED (@value) --
#     including the G loading for specific referents and the specific
#     loadings for the G referent.
#   - var_fixed = TRUE : referent's own-factor loading FREE, variances @1.
#     var_fixed = FALSE: referent's own-factor loading FIXED, variances free.
#   - All factor covariances fixed to 0 (orthogonal bifactor structure).

.besem_ewc_referents <- function(besem_fit, spec) {
  L <- besem_fit$rotated_loadings
  if (is.null(L))
    stop("besem_fit$rotated_loadings is missing -- cannot build B-EWC referents.",
         call. = FALSE)
  g_name  <- besem_fit$g_name %||% "G"
  factors <- names(spec$factors)

  # Row name lookup is case-insensitive to match how find_ewc_referents handles it
  row_lc <- tolower(rownames(L))
  refs   <- setNames(character(length(factors)), factors)

  for (f in factors) {
    prim_lc <- tolower(spec$factors[[f]])
    idx     <- which(row_lc %in% prim_lc)
    if (!length(idx))
      stop("No rows in rotated_loadings matched the primary items of factor '",
           f, "'.", call. = FALSE)
    col <- match(f, colnames(L))
    if (is.na(col))
      stop("rotated_loadings has no column named '", f, "'.", call. = FALSE)
    best       <- idx[which.max(abs(L[idx, col]))]
    # Return in spec's original case
    matched    <- spec$factors[[f]][tolower(spec$factors[[f]]) == row_lc[best]]
    refs[f]    <- if (length(matched)) matched[1L] else rownames(L)[best]
  }

  # G referent: highest |G loading| among items NOT already a specific referent
  g_col       <- match(g_name, colnames(L))
  if (is.na(g_col))
    stop("rotated_loadings has no G column ('", g_name, "').", call. = FALSE)
  used_lc     <- tolower(unname(refs))
  candidates  <- which(!row_lc %in% used_lc)
  if (!length(candidates))
    stop("Cannot pick a G referent: all items are already specific referents.",
         call. = FALSE)
  g_best      <- candidates[which.max(abs(L[candidates, g_col]))]
  g_ref_item  <- rownames(L)[g_best]
  # Try to recover spec-case for the G referent too
  all_spec_lc <- tolower(spec$all_items)
  matched_g   <- spec$all_items[all_spec_lc == tolower(g_ref_item)]
  g_ref_item  <- if (length(matched_g)) matched_g[1L] else g_ref_item

  out <- c(setNames(g_ref_item, g_name), refs)
  out
}

# Parameterization the source lavaan fit was estimated under ("theta" /
# "delta"); "delta" when the slot is absent (custom-DWLS fits, continuous).
.ewc_source_parameterization <- function(esem_fit) {
  p <- tryCatch(lavaan::lavInspect(esem_fit$lavaan_fit, "options")$parameterization,
                error = function(e) NULL)
  if (is.character(p) && length(p) == 1L && p %in% c("theta", "delta")) p else "delta"
}

.besem_ewc_syntax <- function(besem_fit, spec, referents = NULL, var_fixed = TRUE) {
  L <- besem_fit$rotated_loadings
  if (is.null(L))
    stop("besem_fit$rotated_loadings is missing -- cannot build B-EWC syntax.",
         call. = FALSE)

  g_name    <- besem_fit$g_name %||% "G"
  specifics <- names(spec$factors)
  all_fact  <- c(g_name, specifics)

  if (is.null(referents))
    referents <- .besem_ewc_referents(besem_fit, spec)

  # Validate
  if (!all(all_fact %in% names(referents)))
    stop("B-EWC referents must include G ('", g_name, "') and all specific ",
         "factors (", paste(specifics, collapse = ", "), ").", call. = FALSE)

  # owner_of[tolower(item)] -> factor name whose referent is that item
  owner_of <- setNames(names(referents), tolower(unname(referents)))

  # Loading lookup: factor, item (both in original case)
  row_lc <- tolower(rownames(L))
  col_of <- function(f) match(f, colnames(L))

  lkp <- function(f, it) {
    ci <- col_of(f); ri <- match(tolower(it), row_lc)
    if (is.na(ci) || is.na(ri)) return(NA_real_)
    L[ri, ci]
  }

  # Item order: primaries of each specific factor first (A1..A5, C1..C5, ...)
  # then any unclaimed items. Match Mplus reference which groups by factor.
  ordered_items <- unique(c(unlist(spec$factors, use.names = FALSE),
                            spec$all_items))

  build_factor_line <- function(f) {
    terms <- vapply(ordered_items, function(it) {
      val       <- lkp(f, it)
      it_lower  <- tolower(it)
      own       <- owner_of[it_lower]
      is_owner  <- !is.na(own) && own == f
      is_refent <- !is.na(own)   # item is *some* factor's referent

      # Rule:
      #   * non-referent item                              -> FREE  start(v)*item
      #   * referent item, own factor                      -> FREE  (var_fixed)
      #                                                       FIXED (!var_fixed)
      #   * referent item, another factor (cross-loading)  -> FIXED val*item
      is_fixed <- is_refent && (!is_owner || !var_fixed)
      if (!is.na(val)) {
        if (is_fixed)
          sprintf("%.5f*%s", val, it)
        else
          sprintf("start(%.5f)*%s", val, it)
      } else {
        if (is_fixed) it else sprintf("start(0)*%s", it)
      }
    }, character(1L))

    chunks <- split(terms, ceiling(seq_along(terms) / 4L))
    rhs    <- paste(vapply(chunks, paste, character(1L), collapse = " + "),
                    collapse = " +\n       ")
    paste0(f, " =~ ", rhs)
  }

  load_lines <- vapply(all_fact, build_factor_line, character(1L))

  # Orthogonality: every pair of factors fixed to 0
  pairs      <- utils::combn(all_fact, 2L)
  ortho_lines <- apply(pairs, 2L, function(p) sprintf("%s ~~ 0*%s", p[1L], p[2L]))

  header <- paste0(
    "# Bifactor ESEM-within-CFA (B-EWC) lavaan syntax\n",
    "# Mplus BIFACTOR-ESEM-WITHIN-CFA convention\n",
    "# Referents: ",
    paste(paste0(names(referents), " = '", referents, "'"), collapse = ", "),
    "\n",
    "#   start(v)*item   = FREE, rotated loading as starting value\n",
    "#   value*item      = FIXED to rotated loading\n",
    "# var_fixed = ", var_fixed, "\n",
    if (var_fixed)
      "# Factor variances = 1 (std.lv=TRUE); all factor covariances fixed to 0.\n"
    else
      "# Factor variances free (referent own-factor loadings fixed above); all factor covariances fixed to 0.\n"
  )

  paste(c(header, "", load_lines, "", ortho_lines), collapse = "\n")
}
