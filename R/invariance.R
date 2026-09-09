# -- Measurement Invariance Testing --------------------------------------------

#' Measurement Invariance Testing for ESEM Models
#'
#' Tests configural, weak (metric), strong (scalar), and strict invariance for
#' an ESEM model across groups, optionally followed by latent
#' variance/covariance and latent mean invariance (\code{through}). Returns a
#' formatted table of fit indices and chi-square difference tests, analogous
#' to Mplus's multi-group output.
#'
#' @inheritSection doc_estimator_paths Estimator paths
#'
#' @param spec An \code{esem_spec} object from \code{\link{specify_model}} that
#'   includes a \code{group} variable. If \code{spec$group} is \code{NULL}, prints
#'   a message and returns \code{invisible(NULL)} (not an error).
#' @param model Character. Which model type to test: \code{"esem"} (default) for
#'   standard ESEM, or \code{"besem"} for Bifactor ESEM. Ordered B-ESEM uses lavaan
#'   multi-group \code{efa()} with explicit syntax patches (orthogonal rotation,
#'   Theta identification, threshold labels); see \strong{Scope} below.
#' @param missing Character. Missing data handling passed to \code{lavaan::cfa()}.
#'   Default \code{NULL} uses \code{spec$missing} from \code{\link{specify_model}}
#'   (\code{"pairwise"} for ordered data, \code{"listwise"} for continuous).
#' @param verbose Logical. Print progress messages. Default \code{TRUE}. The
#'   `Fitting <level> ...` line is shown in bold blue where the console
#'   supports ANSI colour (RStudio, or a terminal whose \code{TERM} is not
#'   "dumb"); \code{options(bifactory.ansi_colour = FALSE)} turns the colour
#'   off, \code{TRUE} forces it.
#' @param through Character. Last level of the sequence to fit.
#'   \code{"strict"} (default) stops at the four measurement-invariance levels.
#'   \code{"varcov"} adds latent variance/covariance invariance (level 5:
#'   factor variances fixed to 1 in every group and factor covariances equal
#'   across groups, i.e. 0 in every group for B-ESEM). \code{"means"} adds
#'   latent mean invariance (level 6: factor means fixed to 0 in every group).
#'   Levels are cumulative, so \code{"means"} fits all six models.
#' @param cores Integer. Worker processes for the two starting-value fits of
#'   each level (lavaan's default start and \code{start = "simple"}), which are
#'   independent and run at the same time when \code{cores >= 2} (default 2, the
#'   most that helps). The results are identical; a level then takes as long as
#'   its slower start instead of both together (\code{psych::bfi}, three groups,
#'   six levels: 13 instead of 18 minutes). \code{1} fits sequentially. The
#'   workers load this package and the calling session's library path.
#' @param ... Additional arguments passed to the underlying fit function. Do not
#'   pass \code{group}, \code{group_equal}, \code{ordered}, or
#'   \code{parameterization} here -- these are managed internally.
#'
#' @return An object of class \code{"esem_invariance"} containing:
#' \describe{
#'   \item{\code{table}}{Data frame with fit indices and D-statistics.
#'     \code{print()} renders it formatted.}
#'   \item{\code{models}}{Named list of fit objects, one per fitted level:
#'     \code{configural}, \code{weak}, \code{strong}, \code{strict}, and, with
#'     \code{through}, \code{varcov} and \code{means}.
#'     \code{NULL} entries indicate a model that failed to fit.}
#'   \item{\code{lrt}}{Named list of \code{lavTestLRT()} outputs, one per
#'     fitted level after configural (each compared with the level before it).}
#'   \item{\code{spec}}{The original model specification.}
#'   \item{\code{model}}{Character: \code{"esem"} or \code{"besem"}.}
#'   \item{\code{notes}}{Named list, one character vector per level that raised a
#'     lavaan warning while being fitted or compared (singular information matrix
#'     for the robust test statistic, non-positive-definite latent or residual
#'     covariance matrix, negative variance). The warnings are caught, reworded
#'     and shown by \code{print()} instead of being raised, so a run that
#'     concludes does not end with a lavaan warning. Empty when nothing was raised.}
#'   \item{\code{fallback_from}, \code{fallback_note}}{If B-ESEM configural fails to
#'     converge, the result may be from an **ESEM** fallback (no general factor);
#'     not comparable df-for-df to Mplus B-ESEM.}
#' }
#'
#' @section Scope and caveats:
#' \itemize{
#'   \item **Ordered B-ESEM:** Validated against Mplus on BFI for
#'     \eqn{G \in \{2, 3, 4\}} (max \eqn{|\Delta\mathrm{CFI}| \le 0.001},
#'     \eqn{|\Delta\mathrm{SRMR}| \le 0.002}, df match). At \eqn{G \ge 5} with
#'     \eqn{\ge 18} items a warning is issued; verify externally before reporting.
#'   \item **WLSMV SRMR:** The invariance table applies an Mplus-style SRMR
#'     denominator correction for ordered models where implemented.
#'   \item **B-ESEM configural failure:** After retries, may fall back to ESEM
#'     (\code{fallback_from = "besem"}). The printed table is structurally different
#'     from B-ESEM; do not compare to Mplus B-ESEM output.
#'   \item **Empty categories:** Ordered multi-group fits may recode categories
#'     absent in one group (see verbose output).
#' }
#' See \code{system.file("VALIDATION.md", package = "bifactory")} for a summary.
#'
#' @section Constraint mapping:
#'
#' | Level | Continuous (MLR) | Ordered (WLSMV/Theta) |
#' |---|---|---|
#' | Configural | free | free |
#' | Weak | loadings | loadings |
#' | Strong | loadings + intercepts | loadings + thresholds |
#' | Strict | + residuals | + residuals |
#' | Latent var/cov (\code{through = "varcov"}) | + factor variances (1) and covariances | + factor variances (1) and covariances |
#' | Latent means (\code{through = "means"}) | + factor means (0) | + factor means (0) |
#'
#' The last two levels follow Morin, Arens and Marsh (2016): the degrees of
#' freedom grow by \eqn{k(k+1)/2} (variances and covariances of \eqn{k}
#' factors) and then by \eqn{k} (means) per non-reference group.
#'
#' @section Scaled chi-square difference tests:
#'
#' Simple subtraction of scaled chi-square values is not valid for MLR or
#' WLSMV. \code{lavaan::lavTestLRT()} is used, which applies the
#' Satorra-Bentler (2001) correction for MLR and a mean-variance-adjusted
#' difference test for WLSMV. WLSMV difference test results may differ
#' numerically from Mplus's \code{DIFFTEST} procedure.
#'
#' @section Identification across groups:
#'
#' At the configural level \code{std.lv = TRUE} fixes factor variances to 1
#' in all groups. Under weak invariance lavaan automatically frees factor
#' variances in Group 2+ and keeps them at 1 in Group 1. Under strong
#' invariance factor means are freed in Group 2+ and fixed to 0 in Group 1.
#' The optional latent var/cov level fixes the Group 2+ variances back to 1
#' and equates the covariances; the latent mean level fixes the Group 2+
#' means back to 0.
#'
#' With ordered indicators (theta parameterization) the configural model
#' uses the standard identification in every group: residual variances 1,
#' factor means 0, all thresholds free. From the weak level on, the
#' non-reference residual variances and factor means are free and identified
#' through the equal loadings plus partial threshold equalities (two
#' thresholds of one referent item per factor, one threshold of every other
#' item), which is the Mplus convention. The same scheme is not used at the
#' configural level: with free loadings an item's single equated threshold
#' is its only link across groups, and when that threshold lies near zero in
#' a non-reference group the item's residual variance is empirically
#' unidentified (it ran to several hundred on \code{psych::bfi}, with the
#' chi-square drifting by 14 to 36 units between programs). The standard
#' identification is the equivalent, well-conditioned model.
#'
#' Group 1 (the reference group) is the first sorted level of the grouping
#' variable, \code{spec$group_levels}, which is also the reference group in
#' the Mplus syntax written by \code{\link{run_mplus_besem_invariance}}. The
#' fits pass this order to lavaan as \code{group.label}; lavaan's own default
#' would be the order in which the groups appear in the data, and a different
#' reference group gives an equivalent model with the same fit but differently
#' oriented loadings and factor covariances.
#'
#' @seealso \code{\link{specify_model}}, \code{\link{esem}},
#'   \code{\link{esem_ordered}}
#'
#' @examples
#' data("HolzingerSwineford1939", package = "lavaan")
#'
#' spec <- specify_model(
#'   Visual  = c("x1", "x2", "x3"),
#'   Textual = c("x4", "x5", "x6"),
#'   Speed   = c("x7", "x8", "x9"),
#'   data  = HolzingerSwineford1939,
#'   group = "school",
#'   label = "Holzinger-Swineford"
#' )
#'
#' \donttest{
#' # Continuous ESEM measurement invariance across schools
#' inv <- esem_invariance(spec)
#' print(inv)
#'
#' # Access individual model fits
#' summary(inv$models$strong, fit.measures = TRUE, standardized = TRUE)
#' lavaan::lavTestScore(inv$models$strong$lavaan_fit)
#'
#' # Ordered data (WLSMV) and B-ESEM invariance follow the same pattern:
#' # add `ordered = TRUE` to specify_model(), then optionally `model = "besem"`.
#' items <- paste0("x", 1:9)
#' ord <- as.data.frame(lapply(HolzingerSwineford1939[, items], function(v) {
#'   as.integer(cut(v, breaks = quantile(v, probs = seq(0, 1, 0.2)),
#'                  include.lowest = TRUE))
#' }))
#' names(ord) <- items
#' ord$school <- HolzingerSwineford1939$school
#'
#' spec_ord <- specify_model(
#'   Visual = c("x1", "x2", "x3"), Textual = c("x4", "x5", "x6"),
#'   Speed = c("x7", "x8", "x9"),
#'   data = ord, group = "school", ordered = TRUE
#' )
#' inv_besem <- esem_invariance(spec_ord, model = "besem")
#' print(inv_besem)
#' }
#'
#' @export
esem_invariance <- function(spec,
                             model   = c("esem", "besem"),
                             missing = NULL,
                             verbose = TRUE,
                             through = c("strict", "varcov", "means"),
                             cores   = 2L,
                             ...) {

  if (!inherits(spec, "esem_spec"))
    stop("`spec` must be an esem_spec from specify_model().", call. = FALSE)
  if (is.null(spec$group)) {
    message("Invariance testing skipped: no group variable in spec.\n",
            "  Set group in specify_model(), e.g.: specify_model(..., group = 'sex')")
    return(invisible(NULL))
  }

  model   <- match.arg(model)
  through <- match.arg(through)

  cl <- .make_cluster(cores)
  if (!is.null(cl)) on.exit(parallel::stopCluster(cl), add = TRUE)

  if (model == "besem" && is.null(spec$bifactor_target))
    stop(
      "B-ESEM invariance requires `spec$bifactor_target`. ",
      "Ensure specify_model() was called correctly.",
      call. = FALSE
    )

  is_ordered <- !is.null(spec$ordered) && length(spec$ordered) > 0

  if (is.null(missing)) missing <- spec$missing

  # -- Collapse empty response categories (WLSMV multi-group) -----------------
  # lavaan's polychoric estimator fails when a response category is observed in
  # one group but not another (e.g. top category used by Group 1 but never by
  # Group 2).  Remap each such category to the nearest category that is present
  # in all groups.  The modified data is local to this call -- the caller's spec
  # is not affected.
  if (is_ordered && !is.null(spec$group)) {
    recode <- .collapse_empty_cats(spec$data, spec$all_items, spec$group)
    if (length(recode$recoded) > 0L) {
      if (verbose) {
        message("  [NOTE] Empty response categories in >=1 group -- collapsing to nearest non-empty:")
        for (it in names(recode$recoded)) {
          m <- recode$recoded[[it]]
          message(sprintf("    %-10s: %s", it,
                      paste0(names(m), "->", as.character(m), collapse = ", ")))
        }
      }
      spec$data <- recode$data
    }
  }

  # Constraint vectors differ between continuous and ordered estimator paths.
  # For ordered (WLSMV + theta parameterization) Mplus equates thresholds at
  # scalar level, not intercepts.
  # Levels 5 and 6 (Morin's sequence) fix the non-reference groups' latent
  # variances back to 1 and their covariances equal to the reference group's
  # (0 for B-ESEM), then their latent means to 0.
  constraints <- if (is_ordered) {
    list(
      configural = NULL,
      weak       = "loadings",
      strong     = c("loadings", "thresholds"),
      strict     = c("loadings", "thresholds", "residuals"),
      varcov     = c("loadings", "thresholds", "residuals", "lv.variances", "lv.covariances"),
      means      = c("loadings", "thresholds", "residuals", "lv.variances", "lv.covariances", "means")
    )
  } else {
    list(
      configural = NULL,
      weak       = "loadings",
      strong     = c("loadings", "intercepts"),
      strict     = c("loadings", "intercepts", "residuals"),
      varcov     = c("loadings", "intercepts", "residuals", "lv.variances", "lv.covariances"),
      means      = c("loadings", "intercepts", "residuals", "lv.variances", "lv.covariances", "means")
    )
  }

  model_label <- if (model == "besem") "B-ESEM" else "ESEM"

  if (verbose) {
    message("======================================================")
    message(sprintf(" %s Measurement Invariance: %s", model_label, spec$label))
    message(paste(
      " Group:", spec$group,
      " (", length(spec$group_levels), "groups:",
      paste(spec$group_levels, collapse = ", "), ")"
    ))
    message(paste(
      " Estimator:", if (is_ordered) "WLSMV (ordinal)" else paste0(spec$estimator_esem, " (continuous)")
    ))
    message("======================================================\n")
  }

  # Scope notice: the lavaan B-ESEM invariance path has been validated against
  # Mplus on BFI for G in {2, 3, 4} (max |dCFI| = 0.001, max |dSRMR| = 0.002,
  # df match exactly). At G >= 5 with >= 18 items the lavaan optimizer can
  # stall on B-ESEM target rotation; the fit may still complete but is outside
  # the package's validated range.
  if (verbose && model == "besem" && is_ordered) {
    G <- length(spec$group_levels)
    n_items <- length(spec$all_items)
    if (G >= 5L && n_items >= 18L) {
      warning(sprintf(
        "esem_invariance(model = 'besem'): %d groups x %d items is outside ",
        G, n_items),
        "the validated range (G in {2, 3, 4}). Convergence or SE instability ",
        "may occur; verify against an external benchmark before reporting.",
        call. = FALSE)
    }
  }

  all_levels <- c("configural", "weak", "strong", "strict", "varcov", "means")
  levels     <- all_levels[seq_len(match(through, all_levels))]
  labels <- c(
    configural = "1. Configural",
    weak       = "2. Weak (metric)",
    strong     = "3. Strong (scalar)",
    strict     = "4. Strict",
    varcov     = "5. Latent var/cov",
    means      = "6. Latent means"
  )

  # -- Fit the models up to `through` ----------------------------------------------
  fits  <- list()
  notes <- list()
  for (lv in levels) {
    if (verbose) message(.fit_header(labels[lv]))

    r <- .with_notes(tryCatch(
      .fit_invariance_model(
        spec        = spec,
        group_equal = constraints[[lv]],
        is_ordered  = is_ordered,
        model       = model,
        missing     = missing,
        verbose     = verbose,
        cl          = cl,
        ...
      ),
      error = function(e) {
        # Suppress warning for B-ESEM configural: failure is handled
        # gracefully below with an informative fallback to ESEM.
        if (!(model == "besem" && lv == "configural")) {
          warning(
            sprintf("esem_invariance: %s model failed: %s", lv, conditionMessage(e)),
            call. = FALSE
          )
        }
        NULL
      }
    ))
    fit <- r$value
    lv_notes <- unique(c(attr(fit, "notes"), r$notes))
    # replace lavaan's generic admissibility warnings by group-specific findings
    # read from the reported solution (which group, which parameter)
    if (!is.null(fit)) {
      adm <- .inv_admissibility(fit$lavaan_fit)
      lv_notes <- c(lv_notes[!grepl("negative variance|not positive definite", lv_notes)], adm)
    }
    if (length(lv_notes)) notes[[lv]] <- lv_notes

    conv <- if (!is.null(fit)) {
      tryCatch(lavaan::lavInspect(fit$lavaan_fit, "converged"), error = function(e) NA)
    } else {
      FALSE
    }

    if (verbose) {
      if (is.null(fit)) {
        message(" FAILED")
      } else {
        message(if (isTRUE(conv)) " OK" else " WARNING (did not converge)")
      }
      for (n in lv_notes) message("      note: ", n, " (noted, continuing)")
    }

    fits[[lv]] <- fit

    # -- B-ESEM configural fallback -----------------------------------------
    # If the B-ESEM configural model fails or does not converge, the remaining
    # invariance levels cannot be meaningfully compared against the baseline.
    # Fall back to ESEM, which is the most suitable specification that
    # converges across groups -- even though B-ESEM provided the best
    # single-group fit.
    if (model == "besem" && lv == "configural" && !isTRUE(conv)) {
      if (verbose) {
        message("")
        message("  -------------------------------------------------------")
        message("  B-ESEM configural did not converge after 7 retries.")
        message("  Falling back to ESEM (n_factors - 1, no general factor).")
        message("")
        message("  WARNING: ESEM is structurally NOT equivalent to B-ESEM:")
        message("    - ESEM has k specific factors only.")
        message("    - B-ESEM has k specific factors + 1 general factor.")
        message("    - df / chi2 will differ by k * (n_items) per group.")
        message("    - Direct numerical comparison vs Mplus B-ESEM is invalid.")
        message("")
        message("  Suggested next steps:")
        message("    1. Reduce n_groups (the convergence problem grows with groups).")
        message("    2. Inspect single-group fits per country to find the failing one.")
        message("    3. Run Mplus directly via run_comparison(..., mplus_folder=...) to")
        message("       cross-check the lavaan retry sequence.")
        message("  -------------------------------------------------------\n")
      }
      warning(
        "B-ESEM configural did not converge in lavaan after 7 retries. Falling back ",
        "to ESEM (k specific factors, no general factor). The returned invariance ",
        "table is structurally different from a B-ESEM solution and cannot be ",
        "compared df-for-df against Mplus B-ESEM output.",
        call. = FALSE
      )
      inv_esem <- esem_invariance(spec, model = "esem", missing = missing,
                                  verbose = verbose, through = through, ...)
      inv_esem$fallback_from <- "besem"
      inv_esem$fallback_note <- paste0(
        "B-ESEM configural failed to converge after 7 retries. Returned model ",
        "is ESEM (k specific factors, no general factor) -- structurally NOT ",
        "equivalent to B-ESEM. df, chi-square, and fit indices CANNOT be compared ",
        "directly against Mplus B-ESEM output (df gap = n_groups * n_items)."
      )
      return(inv_esem)
    }
  }

  # -- Chi-square difference tests -----------------------------------------------
  # Each comparison: less-constrained vs more-constrained.
  # We extract raw lavaan objects for lavTestLRT() -- it requires lavaan S4
  # objects, not the esem_fit wrappers.
  comparison_pairs <- list(
    weak   = c("configural", "weak"),
    strong = c("weak",       "strong"),
    strict = c("strong",     "strict"),
    varcov = c("strict",     "varcov"),
    means  = c("varcov",     "means")
  )[levels[-1L]]

  lrt_results <- list()
  for (comp in names(comparison_pairs)) {
    lo <- comparison_pairs[[comp]][1]
    hi <- comparison_pairs[[comp]][2]

    if (is.null(fits[[lo]]) || is.null(fits[[hi]])) next

    r <- .with_notes(tryCatch(
      lavaan::lavTestLRT(fits[[lo]]$lavaan_fit, fits[[hi]]$lavaan_fit),
      error = function(e) {
        warning(
          sprintf("lavTestLRT failed (%s vs %s): %s", lo, hi, conditionMessage(e)),
          call. = FALSE
        )
        NULL
      }
    ))
    lrt_results[[comp]] <- r$value
    if (length(r$notes)) notes[[hi]] <- c(notes[[hi]], paste("difference test:", r$notes))
  }

  # -- Build table ---------------------------------------------------------------
  table_out <- .build_invariance_table(fits, lrt_results, levels, labels, is_ordered)

  out <- structure(
    list(
      table  = table_out,
      models = fits,
      lrt    = lrt_results,
      spec   = spec,
      model  = model,
      notes  = notes
    ),
    class = "esem_invariance"
  )

  if (verbose) {
    message("")
    for (line in .inv_conclusion(out)) message("  ", line)
    if (length(notes)) message("  Estimation notes are listed by print(); see inv$notes.")
    message("  Remember: print(inv)")
  }
  out
}


# -- Print method --------------------------------------------------------------

#' Print Method for esem_invariance
#'
#' @param x An \code{esem_invariance} object.
#' @param ... Ignored.
#' @return Invisibly returns \code{x}; called for the side effect of printing
#'   the invariance comparison table.
#' @export
print.esem_invariance <- function(x, ...) {
  is_ordered  <- !is.null(x$spec$ordered) && length(x$spec$ordered) > 0
  model_label <- if (!is.null(x$model) && x$model == "besem") "B-ESEM" else "ESEM"

  cat("======================================================\n")
  cat(sprintf(" %s Measurement Invariance: %s\n", model_label, x$spec$label))
  cat(
    " Group:", x$spec$group,
    " (", length(x$spec$group_levels), "groups:",
    paste(x$spec$group_levels, collapse = ", "), ")\n"
  )
  cat("======================================================\n\n")

  if (!is.null(x$fallback_from) && x$fallback_from == "besem") {
    cat("  WARNING: B-ESEM configural failed to converge after 7 retries.\n")
    cat("  Results below are from ESEM (k specific factors, NO general factor).\n")
    cat("  This is structurally NOT equivalent to B-ESEM. df / chi-square /\n")
    cat("  CFI / RMSEA cannot be compared df-for-df against Mplus B-ESEM output.\n\n")
  }

  .print_invariance_table(x$table, is_ordered)

  failed <- names(Filter(is.null, x$models))
  if (length(failed))
    cat("  NOTE: Models that failed to fit:", paste(failed, collapse = ", "), "\n\n")

  # lavaan's estimation warnings, reworded (see .translate_lavaan_warning()).
  if (length(x$notes)) {
    cat("  Estimation notes (lavaan messages for the reported fit, reworded):\n")
    for (lv in names(x$notes))
      for (n in x$notes[[lv]]) cat(sprintf("    %s: %s\n", lv, n))
    cat("  These describe the fit at that level; the fit indices above are reported\n")
    cat("  as computed. A non-positive-definite matrix at the most constrained level\n")
    cat("  means that level's solution is inadmissible and should not be interpreted.\n\n")
  }
  for (line in .inv_conclusion(x)) cat("  ", line, "\n", sep = "")
  cat("\n")

  cat("Access results:\n")
  for (lv in names(x$models))
    cat(sprintf("  inv$models$%-11s -- esem_fit (%s)\n", lv, lv))
  cat("  inv$lrt$weak           -- lavTestLRT() output (weak vs configural)\n")
  cat("  inv$table              -- data frame of all fit indices\n\n")

  invisible(x)
}


# -- Internal helpers ----------------------------------------------------------

# Convergence retry strategy.
#
# History: Earlier versions tried `c(1e-4, 5e-4, ...)` as progressively LOOSER
# rel.tol values, mirroring Mplus's CONVERGENCE setting.  This is wrong for
# nlminb: lavaan's `rel.tol` is a relative *function-change* tolerance, not a
# parameter-change tolerance like Mplus's CONVERGENCE.  Loosening rel.tol makes
# nlminb stop EARLY at flat saddles (fx values orders of magnitude above the
# true minimum), and lavaan's post-fit gradient check then rejects the result
# with "the optimizer (NLMINB) claimed the model converged, but not all elements
# of the gradient are (near) zero".
#
# Each entry: list(label, ctrl = list passed as `control`, opts = top-level
# lavaan options merged into the cfa() call).
#
# Order: tighter (Mplus-equivalent) first, then progressively allow gradient
# escape valves.  The escape valves (check.gradient = FALSE,
# optim.force.converged = TRUE) accept fits where nlminb itself converged but
# the gradient norm is above lavaan's 0.001 cutoff -- common for large
# multi-group B-ESEM where the function landscape is very flat near the
# minimum.
.CONV_RETRY_SEQ <- list(
  list(label = "rel.tol=1e-4",
       ctrl  = list(iter.max = 10000L, rel.tol = 1e-4),
       opts  = NULL),
  list(label = "rel.tol=1e-3",
       ctrl  = list(iter.max = 10000L, rel.tol = 1e-3),
       opts  = NULL),
  list(label = "rel.tol=1e-3 + check.gradient=FALSE",
       ctrl  = list(iter.max = 10000L, rel.tol = 1e-3),
       opts  = list(check.gradient = FALSE)),
  list(label = "rel.tol=1e-2 + check.gradient=FALSE",
       ctrl  = list(iter.max = 10000L, rel.tol = 1e-2),
       opts  = list(check.gradient = FALSE)),
  list(label = "rel.tol=5e-2 + check.gradient=FALSE",
       ctrl  = list(iter.max = 10000L, rel.tol = 5e-2),
       opts  = list(check.gradient = FALSE)),
  list(label = "rel.tol=5e-2 + force-converged",
       ctrl  = list(iter.max = 10000L, rel.tol = 5e-2),
       opts  = list(check.gradient = FALSE, optim.force.converged = TRUE))
)

# Convergence criteria for Mplus: first attempt 1e-6 (Mplus's default .00005
# stopped the eight-country BAT 2023 means model 15 chi-square units and the
# two-country strict model 71 units above its own optimum), then the retry
# sequence starting at Mplus's default (the eight-country weak/strong/strict
# models converge there but not at 1e-6).  Mplus's CONVERGENCE is a
# parameter-change criterion (not function-change like lavaan's rel.tol), and
# Mplus's optimizer uses 20 steepest-descent iterations before switching to
# quasi-Newton -- so progressive loosening genuinely helps Mplus escape saddles
# (unlike for lavaan's nlminb; see .CONV_RETRY_SEQ above).
.MPLUS_CONV_RETRY <- c(5e-5, 1e-4, 5e-4, 1e-3, 5e-3, 1e-2, 5e-2, 2.5e-2)

.mplus_conv_lines <- function(conv_crit = NULL) {
  if (is.null(conv_crit)) conv_crit <- 1e-6
  crit <- format(conv_crit, scientific = FALSE, trim = TRUE)
  paste0(
    "\n  ITERATIONS = 10000;",
    "\n  H1ITERATIONS = 10000;",
    "\n  CONVERGENCE = ", crit, ";",
    "\n  H1CONVERGENCE = ", crit, ";"
  )
}

# DIFFTEST line for the next level, or "" when the previous level saved no
# derivatives file (Mplus prints no fit and saves nothing when the information
# matrix is singular, e.g. the eight-country BAT 2023 var/cov model); an input
# naming a missing file fails with "*** ERROR in ANALYSIS command".
.mplus_difftest_line <- function(dat_path) {
  if (file.exists(dat_path))
    return(paste0("\n  DIFFTEST = ", basename(dat_path), ";"))
  message(sprintf("  (no %s from the previous level: DIFFTEST omitted)", basename(dat_path)))
  ""
}

.mplus_converged <- function(out_path) {
  if (!file.exists(out_path)) return(FALSE)
  txt <- readLines(out_path, warn = FALSE)
  any(grepl("THE MODEL ESTIMATION TERMINATED NORMALLY", txt))
}


.inv_converged <- function(fit) {
  if (is.null(fit)) return(FALSE)
  isTRUE(tryCatch(
    lavaan::lavInspect(fit$lavaan_fit, "converged"),
    error = function(e) FALSE
  ))
}

.inv_fx <- function(fit) fit$lavaan_fit@optim$fx

# lavaan's own warnings ("lavaan->lav_test_sb(): could not invert ...") read
# like failures to users; esem_invariance() catches them per level, on the
# master and on the workers, and reports them reworded (inv$notes, print()).
# Package warnings are not lavaan-prefixed and pass through untouched.
.translate_lavaan_warning <- function(msg) {
  m <- gsub("\\s+", " ", msg)
  if (grepl("could not invert information matrix", m))
    return("robust test statistic could not be computed (singular information matrix)")
  if (grepl("latent variables is not positive definite", m))
    return("latent covariance matrix (psi) is not positive definite in at least one group")
  if (grepl("residual", m) && grepl("not positive definite", m))
    return("residual covariance matrix (theta) is not positive definite in at least one group")
  if (grepl("negative", m) && grepl("variance", m))
    return("a negative variance estimate (Heywood case)")
  sub("^lavaan(->[^:]*\\(\\))?[ :]*(WARNING:)?\\s*", "lavaan: ", m)
}

.with_notes <- function(expr) {
  notes <- character()
  value <- withCallingHandlers(expr, warning = function(w) {
    msg <- conditionMessage(w)
    if (grepl("^lavaan", msg)) {
      notes <<- c(notes, .translate_lavaan_warning(msg))
      invokeRestart("muffleWarning")
    }
  })
  list(value = value, notes = unique(notes))
}

# Most constrained level (strict, strong, weak) whose own transition dCFI passes
# the cutoff; configural otherwise. Shared by factor_scores(level = "auto")
# and the closing line of esem_invariance() / print().
# Admissibility of the reported fit, read from the solution itself so the
# note can say which group and which parameter: every variance estimate below
# zero (latent or residual) and every group whose latent covariance matrix is
# not positive definite. Returns character(0) when the solution is admissible.
.inv_admissibility_pt <- function(pt, group_labels, psi) {
  out <- character()
  v <- pt[pt$op == "~~" & pt$lhs == pt$rhs & !is.na(pt$est) & pt$est < 0, , drop = FALSE]
  if (nrow(v)) {
    lat <- unique(pt$lhs[pt$op == "=~"])
    by_g <- vapply(split(v, v$group), function(d) {
      kind <- ifelse(d$lhs %in% lat, "latent variance", "residual variance")
      paste0(group_labels[d$group[1]], " (", paste(sprintf("%s %s %.3f", d$lhs, kind, d$est), collapse = "; "), ")")
    }, character(1))
    out <- c(out, paste0("a negative variance estimate (Heywood case): ", paste(by_g, collapse = ", ")))
  }
  npd <- vapply(psi, function(m) {
    m <- as.matrix(m)
    nrow(m) > 0 && any(eigen(m, symmetric = TRUE, only.values = TRUE)$values < -1e-8)
  }, logical(1))
  if (any(npd))
    out <- c(out, paste0("latent covariance matrix (psi) is not positive definite: ",
                         paste(group_labels[which(npd)], collapse = ", ")))
  out
}

.inv_admissibility <- function(lav) {
  tryCatch({
    pt  <- lavaan::parTable(lav)
    lab <- lavaan::lavInspect(lav, "group.label")
    if (!length(lab)) lab <- "1"
    psi <- lavaan::lavInspect(lav, "cov.lv")
    if (!is.list(psi)) psi <- list(psi)
    .inv_admissibility_pt(pt, lab, psi)
  }, error = function(e) character())
}

# A level whose reported solution is inadmissible (negative variance,
# non-positive-definite psi or theta) is skipped even if its dCFI passes; a
# note about the test statistic alone does not disqualify it.
.inv_inadmissible <- function(notes) {
  notes[grepl("negative variance|not positive definite", notes)]
}

.inv_supported_level <- function(table, models, cutoff = -0.010, notes = NULL) {
  # every fitted level counts, most constrained first (levels 5 and 6 when
  # esem_invariance() was run with through = "varcov" / "means")
  labels <- c(means = "6. Latent means", varcov = "5. Latent var/cov", strict = "4. Strict",
              strong = "3. Strong (scalar)", weak = "2. Weak (metric)")
  for (lv in names(labels)) {
    if (is.null(models[[lv]]) || length(.inv_inadmissible(notes[[lv]]))) next
    row <- table[table$Model == labels[[lv]], , drop = FALSE]
    dcfi <- if (nrow(row) && "dCFI" %in% names(row)) row$dCFI[1] else NA_real_
    if (!is.na(dcfi) && dcfi >= cutoff) return(lv)
  }
  "configural"
}

.inv_conclusion <- function(x, cutoff = -0.010) {
  fitted <- sum(!vapply(x$models, is.null, logical(1)))
  by_fit <- .inv_supported_level(x$table, x$models, cutoff)          # dCFI only
  lv     <- .inv_supported_level(x$table, x$models, cutoff, x$notes)  # admissible too
  out <- sprintf("Invariance testing concluded: %d of %d levels fitted.", fitted, length(x$models))
  out <- c(out, paste("Levels considered:", paste(names(x$models), collapse = ", ")),
           "The invariance fit table can be interpreted.")
  if (by_fit != lv) {
    out <- c(out, sprintf("%s is supported by dCFI but its solution is inadmissible (%s); its parameters and factor scores should not be used.",
                          by_fit, paste(.inv_inadmissible(x$notes[[by_fit]]), collapse = "; ")))
    # At strict the residual variances are already equal, so the offending
    # variance is a latent one in some group; equating the latent variances
    # (level 5) often removes it, and partial_invariance() is the other route.
    if (by_fit == "strict" && is.null(x$models$varcov))
      out <- c(out, paste0('Tip: fit through = "varcov" or through = "means" (levels 5 and 6 equate the latent ',
                           'variances/covariances and means across groups, which often removes a negative latent ',
                           'variance in one group), or release the offending item with partial_invariance(); ',
                           'then choose that level explicitly, e.g. factor_scores(inv, level = "varcov").'))
  }
  out <- c(out,
    if (lv == "configural")
      sprintf("No %slevel beyond configural is supported by dCFI >= %.3f; parameters differ by group (inv$models$configural), or try partial_invariance().",
              if (by_fit != lv) "admissible " else "", cutoff)
    else
      sprintf("The most constrained %slevel supported by dCFI >= %.3f is: %s. Use its parameters (loadings, thresholds, ...): parameters(inv$models$%s).",
              if (by_fit != lv) "admissible " else "", cutoff, lv, lv))
  if (length(x$notes[[lv]]))
    out <- c(out, sprintf("Note at that level: %s. Inspect the solution before reporting it.",
                          paste(x$notes[[lv]], collapse = "; ")))
  out
}

# "Fitting <level> ..." in bold blue where the console supports it, so each
# level's block stands out in a long log.
.fit_header <- function(label) {
  txt <- sprintf("  Fitting %-22s", paste0(label, " ..."))
  if (.ansi_colour()) paste0("\033[1;34m", txt, "\033[0m") else txt
}

.set_notes <- function(fit, notes) {
  if (!is.null(fit) && length(notes)) attr(fit, "notes") <- unique(c(attr(fit, "notes"), notes))
  fit
}


# Starting values for the ordered multi-group configural fit: each group fitted
# alone (default and simple start, deeper kept) and the unrotated free-parameter
# values collected in group order.  The configural model has no cross-group
# constraint, so its optimum is the sum of the per-group optima, but the joint
# optimizer's path depends on the number of groups: on four BAT 2023 countries
# both joint starts left NL in a well 9 % shallower than the one NL alone
# reaches in five seconds.  Returned as a parameter table (lhs/op/rhs/group/est,
# which lav_start() matches row by row) rather than a numeric vector: lavaan caps
# a numeric `start` at 1000 entries (lav_options_check), and four bfi education
# tiers already need 1040, eight BAT countries 1576.  Returns NULL when any group
# fails; attribute "fx" holds the per-group fit function values.
.configural_starts <- function(spec, model, missing, cl = NULL, ...) {
  if (is.null(spec$group) || length(spec$group_levels) < 2L) return(NULL)
  builder <- if (model == "besem") .fit_besem_inv_ordered else .fit_esem_inv_ordered
  dots <- list(...); dots$group.label <- NULL
  fits <- lapply(spec$group_levels, function(g) {
    sg <- spec
    sg$data <- spec$data[spec$data[[spec$group]] %in% g, , drop = FALSE]
    sg$group <- NULL; sg$group_levels <- NULL
    tryCatch(.fit_with_retry(function(ctrl = NULL, opts = NULL) {
      args <- c(list(spec = sg, group_equal = NULL, missing = missing), dots)
      if (!is.null(ctrl)) args$control <- ctrl
      args[names(opts)] <- opts
      suppressWarnings(do.call(builder, args))
    }, verbose = FALSE, cl = cl), error = function(e) NULL)
  })
  if (any(vapply(fits, is.null, logical(1)))) return(NULL)
  tabs <- lapply(seq_along(fits), function(i) {
    lf <- fits[[i]]$lavaan_fit
    pt <- lavaan::parTable(lf)
    # parTable() holds the rotated solution; the optimizer works on the unrotated
    # one (the efa-fixed loadings dropped), which lavaan keeps in est.unrotated.
    if (!is.null(pt$est.unrotated)) { pt$est <- pt$est.unrotated; pt$free <- pt$free.unrotated }
    if (sum(pt$free > 0L) != length(lf@optim$x)) return(NULL)
    data.frame(lhs = pt$lhs, op = pt$op, rhs = pt$rhs, group = i, free = pt$free, est = pt$est)
  })
  if (any(vapply(tabs, is.null, logical(1)))) return(NULL)
  x0 <- do.call(rbind, tabs)
  attr(x0, "fx") <- vapply(fits, .inv_fx, numeric(1))
  x0
}


# Worker processes for the two starts of each level (esem_invariance(cores = )).
# NULL (sequential) for cores < 2 or when no cluster can be started.  The workers
# get the master's library path and this package: the installed namespace, or the
# load_all() source tree while developing (pkgload marks such a namespace).
.make_cluster <- function(cores) {
  cores <- suppressWarnings(as.integer(cores)[1])
  if (is.na(cores) || cores < 2L) return(NULL)
  n <- parallel::detectCores()
  if (!is.na(n)) cores <- min(cores, n)
  cl <- tryCatch(parallel::makeCluster(cores), error = function(e) NULL)
  if (is.null(cl)) {
    message("  (no worker processes available: starts run sequentially)")
    return(NULL)
  }
  # Nothing written in this file goes over the wire here: such a closure carries
  # this namespace as its environment, which a worker cannot resolve before it
  # has loaded the package.  .libPaths goes by name, because the function object
  # would bring its own copy of base's .lib.loc environment along and set that
  # copy instead of the worker's library path.
  ns <- asNamespace("bifactory")
  ok <- tryCatch({
    parallel::clusterCall(cl, ".libPaths", .libPaths())
    if (exists(".__DEVTOOLS__", envir = ns, inherits = FALSE) &&
        requireNamespace("pkgload", quietly = TRUE)) {
      parallel::clusterCall(cl, pkgload::load_all, getNamespaceInfo(ns, "path"),
                            quiet = TRUE)
    } else {
      parallel::clusterCall(cl, "loadNamespace", "bifactory")
    }
    TRUE
  }, error = function(e) FALSE)
  if (!ok) {
    parallel::stopCluster(cl)
    message("  (workers could not load bifactory: starts run sequentially)")
    return(NULL)
  }
  cl
}


# Relative fit-function gain an alternative start needs before it replaces the
# default start's converged fit (0.1 %; see the comment inside .fit_with_retry).
.ALT_START_GAIN <- 1e-3

.fit_with_retry <- function(fit_fn, verbose = TRUE, alt_start = "simple", cl = NULL) {

  # Attempt 0: default lavaan settings (tight rel.tol = 1e-10, check.gradient = TRUE).
  # Best fit quality if it converges; matches non-invariance esem/besem behaviour.
  parallel_starts <- !is.null(cl) && !is.null(alt_start)
  if (verbose) message("\n    [attempt 0/", length(.CONV_RETRY_SEQ), "] default",
                       if (parallel_starts) " and simple start on 2 workers", " ...")
  alt <- NULL
  if (parallel_starts) {
    # The two starts are independent fits: both at once on the worker processes.
    both <- parallel::parLapply(
      cl, list(NULL, list(start = alt_start)),
      function(o, fit_fn) tryCatch(.with_notes(suppressMessages(fit_fn(ctrl = NULL, opts = o))),
                                   error = function(e) NULL),
      fit_fn = fit_fn)
    fit <- .set_notes(both[[1]]$value, both[[1]]$notes)
    alt <- .set_notes(both[[2]]$value, both[[2]]$notes)
  } else {
    r <- tryCatch(.with_notes(fit_fn(ctrl = NULL, opts = NULL)), error = function(e) {
      if (verbose) message(sprintf(" [lavaan ERROR] %s", conditionMessage(e)))
      NULL
    })
    fit <- .set_notes(r$value, r$notes)
  }
  if (.inv_converged(fit)) {
    # Second start: the multi-group fit function can have several wells (BAT
    # 2023, two-country configural B-ESEM: the default start stopped 9 % above
    # the well Mplus reaches).  Keep the deeper of the two converged fits.
    # alt_start = NULL skips this (configural fits started from the per-group
    # optima, .configural_starts()).
    if (!parallel_starts && !is.null(alt_start)) {
      r <- tryCatch(.with_notes(suppressMessages(fit_fn(ctrl = NULL, opts = list(start = alt_start)))),
                    error = function(e) NULL)
      alt <- .set_notes(r$value, r$notes)
    }
    # Only a materially deeper alternative replaces the default fit.  The
    # eight-country BAT 2023 means and strong models have several optima whose
    # fit functions differ by 0.003 % (section 5c of the BAT 2023 note) while
    # their loadings differ by 0.11; switching for such a gain only moves the
    # reported solution away from the default start's well, which is where
    # Mplus's default start lands.  Real wells differ by 1 % or more (two- and
    # four-country configural: 9 %, eight-country configural: 1.2 %).
    if (.inv_converged(alt) &&
        isTRUE(.inv_fx(alt) < .inv_fx(fit) * (1 - .ALT_START_GAIN))) fit <- alt
    if (verbose) message(" OK")
    return(fit)
  }

  # Attempts 1..N from .CONV_RETRY_SEQ: progressively looser rel.tol with
  # gradient-check escape valves (see comments on .CONV_RETRY_SEQ).
  n_retry <- length(.CONV_RETRY_SEQ)
  for (i in seq_along(.CONV_RETRY_SEQ)) {
    a <- .CONV_RETRY_SEQ[[i]]
    if (verbose)
      message(sprintf("\n    [attempt %d/%d] %s ...", i, n_retry, a$label))

    r <- tryCatch(.with_notes(suppressMessages(fit_fn(ctrl = a$ctrl, opts = a$opts))),
                  error = function(e) NULL)
    fit <- .set_notes(r$value, r$notes)
    if (.inv_converged(fit)) {
      if (verbose) message(" OK")
      return(fit)
    }
    if (verbose) message(" failed")
  }

  stop(
    n_retry + 1L, " convergence strategies were tried (default tight rel.tol, ",
    "progressive loosening, check.gradient = FALSE escape valves) and the ",
    "model likely will not converge with this package.  For 5+ group B-ESEM, ",
    "consider running run_comparison(..., mplus_folder = ...) and using Mplus ",
    "output for invariance comparisons.",
    call. = FALSE
  )
}


.fit_invariance_model <- function(spec, group_equal, is_ordered, model,
                                   missing, verbose = TRUE, cl = NULL, ...) {

  # Closure that performs one fit attempt.
  # `ctrl` is either NULL or a named list passed as lavaan's `control` argument
  # (e.g. list(iter.max = 10000L, rel.tol = 1e-3)).
  # `opts` is either NULL or a named list of top-level lavaan options merged
  # into the cfa() call (e.g. list(check.gradient = FALSE)).
  make_fit <- function(ctrl = NULL, opts = NULL) {
    extra <- list()
    if (!is.null(ctrl)) extra$control <- ctrl
    # Reference group = first sorted level (specify_model()), as in the Mplus
    # generator. lavaan's default is the order of appearance in the data.
    if (!is.null(spec$group) && !is.null(spec$group_levels) &&
        is.null(list(...)$group.label))
      extra$group.label <- as.character(spec$group_levels)
    # `opts` (retry escape valves, start values) override anything in `...`.
    with_opts <- function(args) { args <- c(args, extra); args[names(opts)] <- opts; args }

    if (model == "besem") {
      if (is_ordered) {
        do.call(
          .fit_besem_inv_ordered,
          with_opts(c(list(spec = spec, group_equal = group_equal, missing = missing),
                      list(...)))
        )
      } else {
        do.call(
          .fit_besem_inv_continuous,
          with_opts(c(list(spec = spec, group_equal = group_equal, missing = missing),
                      list(...)))
        )
      }
    } else if (is_ordered) {
      # Multi-group ordered ESEM goes through .fit_esem_inv_ordered(), which
      # adds the same Theta-parameterization identification fixes that the
      # B-ESEM helper applies (free residuals/means in non-reference groups +
      # partial threshold equality at configural/weak; explicit threshold
      # labels + residual control at strong/strict).  Without these, lavaan
      # silently keeps residuals = 1 and factor means = 0 in every group,
      # producing df mismatches vs Mplus on multi-group fits.  The
      # single-group case (no group, or n_groups == 1) bypasses identification
      # fixes entirely (needs_theta_id = FALSE), which preserves the existing
      # behaviour for non-invariance callers.
      do.call(
        .fit_esem_inv_ordered,
        with_opts(c(list(spec = spec, group_equal = group_equal, missing = missing),
                    list(...)))
      )
    } else {
      suppressMessages(do.call(
        esem,
        with_opts(c(list(
            data         = spec$data,
            nfactors     = spec$nfactors,
            indicators   = spec$all_items,
            rotation     = "target",
            target       = spec$target,
            factor_names = spec$factor_names,
            estimator    = spec$estimator_esem,
            group        = spec$group,
            group_equal  = group_equal,
            missing      = missing,
            heywood_fix  = FALSE          # invariance models must stay nested
          ),
          list(...)))
      ))
    }
  }

  # Ordered multi-group configural: start the joint fit from the per-group
  # optima (.configural_starts(), a parameter table); no second joint start is
  # needed then.  Any failure falls back to the plain two-start fit.
  if (is_ordered && is.null(group_equal) && !is.null(spec$group) &&
      length(spec$group_levels) > 1L && model %in% c("besem", "esem")) {
    x0 <- .configural_starts(spec, model = model, missing = missing, cl = cl, ...)
    if (!is.null(x0)) {
      fit <- tryCatch(
        .fit_with_retry(function(ctrl = NULL, opts = NULL)
                          make_fit(ctrl, c(list(start = x0), opts)),
                        verbose = verbose, alt_start = NULL),
        error = function(e) NULL)
      if (!is.null(fit)) return(fit)
    }
  }

  .fit_with_retry(make_fit, verbose = verbose, cl = cl)
}


.collapse_empty_cats <- function(data, items, group) {
  # Force a plain data.frame copy so data.table / tibble reference semantics
  # cannot propagate changes back to the caller's spec$data.
  data    <- as.data.frame(data)
  groups  <- sort(unique(na.omit(data[[group]])))
  recoded <- list()

  # Helper: extract numeric category values from a column that may be a factor
  # whose levels are numeric labels (e.g. "1","2",..."6").  as.integer() on a
  # factor gives level *indices*, not label values -- use the labels directly.
  .cat_vals <- function(col) {
    if (is.factor(col)) {
      lv <- suppressWarnings(as.integer(levels(col)))
      if (!anyNA(lv)) return(lv[as.integer(col)])   # numeric labels: use them
    }
    as.integer(col)   # numeric/integer column or non-numeric factor -> indices
  }

  for (item in items) {
    cats_per_group <- lapply(groups, function(g) {
      vals <- data[[item]][!is.na(data[[group]]) & data[[group]] == g]
      sort(unique(.cat_vals(na.omit(vals))))
    })

    common_cats <- Reduce(intersect, cats_per_group)
    all_cats    <- sort(unique(unlist(cats_per_group)))
    problem     <- setdiff(all_cats, common_cats)
    if (length(problem) == 0L || length(common_cats) == 0L) next

    old_col <- .cat_vals(data[[item]])
    new_col <- old_col
    changed <- integer(0)
    for (v in problem) {
      target <- common_cats[which.min(abs(common_cats - v))]
      new_col[!is.na(old_col) & old_col == v] <- target
      changed <- c(changed, setNames(target, as.character(v)))
    }
    data[[item]] <- new_col
    recoded[[item]] <- changed
  }

  list(data = data, recoded = recoded)
}


# Mute lavaan's "vcov not positive definite" warning around target-rotated
# ESEM/B-ESEM cfa() calls. Rotation indeterminacy makes the parameter vcov
# structurally near-singular; the warning is benign here and confuses users.
.muffle_rotated_vcov <- function(expr) {
  withCallingHandlers(
    expr,
    warning = function(w) {
      if (grepl("lav_model_vcov|positive definite",
                conditionMessage(w), ignore.case = TRUE))
        invokeRestart("muffleWarning")
    }
  )
}


.fit_besem_inv_ordered <- function(spec, group_equal, missing, ...) {
  # lavaan-style group.partial labels ("x9|t3", "x9~~x9") are honoured by the
  # threshold/residual writer below: the explicit equality labels written there
  # would otherwise defeat lavaan's own group.partial.
  dots <- list(...)
  group_partial <- dots$group.partial
  dots$group.partial <- NULL

  # Build orthogonal target: 1 -> NA (free), 0 -> 0 (target toward zero).
  # make_bifactor_target() produces a 0/1 numeric matrix; lavaan's rotation uses
  # NA to mean "no penalty" and 0 to mean "rotate toward zero". This matches
  # Mplus: primary/G loadings are free, cross-loadings are targeted to 0.
  btgt    <- unclass(spec$bifactor_target)   # strip esem_target S3 class
  btgt_na <- btgt
  btgt_na[btgt != 0] <- NA_real_

  all_factor_names <- colnames(btgt)   # c("G", "EX", "MD", "CI")
  lhs          <- paste(paste0('efa("besem")*', all_factor_names), collapse = " + ")
  model_syntax <- paste0(lhs, " =~ ", paste(spec$all_items, collapse = " + "))

  # -- Rotation strategy --------------------------------------------------------
  # Configural: orthogonal=TRUE enforces zero factor covariances in ALL groups
  #   via the rotation algorithm.  Both groups have covariances = 0. Matches Mplus.
  #
  # Weak / Strong / Strict: Mplus applies the orthogonal rotation only to the
  #   SHARED loading matrix (reference group identification).  The non-reference
  #   group's k*(k-1)/2 factor covariances are ADDITIONAL FREE parameters (not
  #   constrained by the rotation), adding 6 free params and reducing df by 6.
  #
  #   lavaan's orthogonal=TRUE alone fixes covariances to 0 in ALL groups (over-
  #   constraining by 6, hence df_R = df_M + 6 at the weak level).
  #
  #   Fix: use orthogonal=TRUE (so the rotation enforces zero covariances in the
  #   reference group) and add explicit c(0, NA)* covariance syntax.  The explicit
  #   NA for the non-reference group overrides the rotation's global orthogonality
  #   constraint, allowing those covariances to be freely estimated.  This produces
  #   the same parameterisation as Mplus: reference group truly orthogonal
  #   (rotation-enforced), non-reference group has free covariances.
  #
  #   NOTE: using orthogonal=FALSE (oblique) with c(0, NA)* does NOT work — the
  #   oblique rotation determines factor covariances implicitly, and the explicit
  #   constraint is ignored.  Group 1 specific-specific covariances end up ~0.3–0.4
  #   instead of 0.  This produces a correct df but a mis-parameterised model.

  is_configural <- is.null(group_equal)
  # Levels 5/6 of esem_invariance(): latent var/cov and latent mean invariance.
  fix_varcov <- !is.null(group_equal) && "lv.variances" %in% group_equal
  fix_means  <- !is.null(group_equal) && "means"        %in% group_equal
  # Use spec$group_levels (na.omit'd in specify_model) so NA values in the
  # grouping column are not counted as an extra group.  Fall back to an
  # na.omit computation if group_levels is absent for any reason.
  n_groups <- if (!is.null(spec$group_levels)) {
    length(spec$group_levels)
  } else if (!is.null(spec$group)) {
    length(unique(na.omit(spec$data[[spec$group]])))
  } else {
    1L
  }

  if (!is_configural && !is.null(spec$group)) {
    # Latent var/cov invariance: the non-reference covariances go back to 0.
    grp_spec  <- paste(c("0", rep(if (fix_varcov) "0" else "NA", n_groups - 1L)),
                       collapse = ", ")
    fac_pairs <- combn(all_factor_names, 2L, simplify = FALSE)
    cov_lines <- vapply(fac_pairs,
      function(p) paste0(p[1L], " ~~ c(", grp_spec, ")*", p[2L]),
      character(1L))
    model_syntax <- paste(c(model_syntax, cov_lines), collapse = "\n")
    rot_args <- list(orthogonal = TRUE, target = btgt_na, rstarts = 30L)

  } else {
    rot_args <- list(orthogonal = TRUE, target = btgt_na, rstarts = 30L)
  }

  # -- Theta identification: match Mplus at weak -------------------------------
  #
  # At weak, Mplus frees the non-reference residuals and factor means and
  # identifies them through partial threshold equalities; lavaan's default
  # keeps residuals = 1 and means = 0 in ALL groups, a different model with
  # the same df (chi-square gap of ~18 when group thresholds truly differ).
  #
  # Not at configural: with free loadings the equated threshold is an item's
  # only cross-group link, and when it sits near zero in a non-reference group
  # the item's residual variance is empirically unidentified (bfi O5: 333 in
  # R, 145 in Mplus at a tight tolerance; chi-square 14 units apart at the
  # default tolerances). The standard identification (residuals 1, means 0,
  # thresholds free in every group) is the equivalent, well-conditioned model;
  # both programs then agree to about one chi-square unit.
  #
  # Applied at weak (strong/strict handled separately below):
  #   1. Free residuals in non-reference group: item ~~ c(1, NA)*item
  #   2. Free factor means in non-reference group: F ~ c(0, NA)*1
  #   3. Partial threshold constraints for identification:
  #      - Referent items (1 per factor): 2 thresholds constrained equal
  #      - Non-referent items: 1 threshold constrained equal
  #
  # Net df change = 0: freed residuals + freed means = threshold constraints.
  needs_theta_id <- !is_configural && !is.null(group_equal) &&
    "loadings" %in% group_equal && !("thresholds" %in% group_equal)

  if (needs_theta_id && n_groups >= 2L) {

    # Referent indicators: first item of each specific factor + 2nd item of
    # first factor (for G).  Matches the Mplus syntax generator (.besem_thresh_lines).
    factor_names_spec_w <- names(spec$factors)
    refs2 <- unique(c(
      vapply(spec$factors, `[`, character(1), 1),    # 1st item of each specific factor
      spec$factors[[1]][2]                            # G referent = 2nd item of factor 1
    ))
    refs2 <- refs2[!is.na(refs2)]

    # 1. Free factor means in non-reference group(s)
    grp_spec_mean <- paste(c("0", rep("NA", n_groups - 1L)), collapse = ", ")
    mean_lines <- vapply(all_factor_names, function(f)
      paste0(f, " ~ c(", grp_spec_mean, ")*1"), character(1L))
    model_syntax <- paste(c(model_syntax, mean_lines), collapse = "\n")

    # 2. Free residual variances in non-reference group(s)
    res_vals <- paste(c("1", rep("NA", n_groups - 1L)), collapse = ", ")
    res_lines <- vapply(spec$all_items, function(it)
      paste0(it, " ~~ c(", res_vals, ")*", it), character(1L))
    model_syntax <- paste(c(model_syntax, res_lines), collapse = "\n")

    # 3. Partial threshold constraints for identification
    lbl_n <- 1L
    thresh_lines <- character(0)
    for (it in spec$all_items) {
      n_thr <- length(unique(na.omit(spec$data[[it]]))) - 1L
      if (n_thr == 0L) next
      n_constrained <- if (it %in% refs2) min(2L, n_thr) else 1L
      parts <- character(n_thr)
      for (k in seq_len(n_thr)) {
        if (k <= n_constrained) {
          lbl     <- paste0("wthr", lbl_n)
          grp_lbl <- paste(rep(lbl, n_groups), collapse = ", ")
          parts[k] <- paste0("c(", grp_lbl, ")*t", k)
          lbl_n   <- lbl_n + 1L
        } else {
          parts[k] <- paste0("t", k)
        }
      }
      thresh_lines <- c(thresh_lines,
        paste0(it, " | ", paste(parts, collapse = " + ")))
    }
    model_syntax <- paste(c(model_syntax, thresh_lines), collapse = "\n")
  }

  # -- Strong/Strict: replace group.equal="thresholds" with full explicit labeling --
  #
  # Using group.equal="thresholds" alongside the c(0,NA)* factor-covariance syntax
  # causes lavaan to add only 36 of the expected 50 threshold equality constraints
  # (df_R=266 vs df_M=280 at strong).
  #
  # The correct approach (per besem-proof/ Mplus files and multigroup.pdf):
  # specify ALL threshold equality via explicit c(lbl, lbl)*tN labels in the model
  # syntax -- the same way Mplus uses identical threshold labels in MODEL and MODEL G2.
  # Removing "thresholds" from group.equal avoids the conflict.
  #
  # df accounting (18 items  x  4 thresholds = 72 total):
  #   - lavaan WEAK identification fixes some thresholds implicitly.
  #   - Explicit equality for all 72 thresholds adds constraints for the ones not
  #     yet fixed, targeting df_strong = 280 (matching Mplus).
  #   - NOTE: lavaan fixes factor means at 0 in all groups (std.lv=TRUE), while
  #     Mplus frees factor means in group 2+ at strong.  This difference in mean
  #     identification may produce a small df offset (+/-4) that should be verified
  #     empirically after running the model.
  #   - Strict: group.equal="residuals" adds 18 more -> target df_strict = 298.
  if (!is_configural && !is.null(group_equal) && "thresholds" %in% group_equal) {

    # Capture whether strict residual equality is requested BEFORE removing it.
    # We handle residuals explicitly in the model string (see below) so that we
    # can precisely control which groups have fixed vs free residuals -- lavaan's
    # Theta parameterization silently fixes ALL groups' residuals to 1 when all
    # thresholds are explicitly constrained equal, which adds 18 hidden constraints
    # at strong and makes strict indistinguishable from strong.
    # Covariances and means are written explicitly below, so their group.equal
    # tokens are dropped too; "lv.variances" stays for lavaan to pin the
    # non-reference variances back at 1.
    has_residuals <- "residuals" %in% group_equal
    group_equal   <- setdiff(group_equal,
                             c("thresholds", "residuals", "lv.covariances", "means"))

    # All thresholds equal: item | c(lbl,lbl)*t1 + c(lbl,lbl)*t2 + ...
    # The part AFTER '*' must be a POSITIONAL indicator (t1, t2, ...) -- not an
    # arbitrary name.  Using any other token causes lavaan to fail with a C-level
    # "subscript out of bounds" error when it tries to look up the threshold index.
    pl <- .ordered_partial_lines(spec, n_groups, has_residuals, group_partial)
    model_syntax <- paste(c(model_syntax, pl$thresholds), collapse = "\n")

    # Free factor means in group 2+ (Mplus strong: MODEL FEMALE has [G*] [EX*] etc.)
    # With explicit threshold equality the scale is fully identified, so factor means
    # in the non-reference group(s) can be freely estimated for mean comparison.
    # lavaan's std.lv=TRUE + group.equal="thresholds" would do this automatically,
    # but since we bypassed "thresholds", we add it explicitly.
    # Syntax: "F ~ c(0, NA, NA, ...)*1" -- 0 in reference group, NA (free) elsewhere.
    # Latent mean invariance: the non-reference means go back to 0.
    grp_spec_mean <- paste(c("0", rep(if (fix_means) "0" else "NA", n_groups - 1L)),
                           collapse = ", ")
    mean_lines <- vapply(all_factor_names, function(f) {
      paste0(f, " ~ c(", grp_spec_mean, ")*1")
    }, character(1L))
    model_syntax <- paste(c(model_syntax, mean_lines), collapse = "\n")

    # Explicit residual variance constraints -- prevents lavaan's Theta from silently
    # fixing ALL groups' residuals to 1 when all thresholds are explicitly equal.
    #
    # Strong (Mplus): MALE residuals = 1 (ref group, Theta), FEMALE residuals FREE.
    #   -> c(1, NA) keeps non-reference group residuals estimable.
    #   -> Net: -18 free params released vs lavaan's default "all Theta = 1".
    #
    # Strict (Mplus): residuals equal across groups (both = 1, from Theta).
    #   -> c(1, 1) explicitly fixes both = 1, equivalent to residual equality.
    #   -> +18 constraints vs strong -> df_strict = df_strong + 18.
    model_syntax <- paste(c(model_syntax, pl$residuals), collapse = "\n")
  }

  cfa_args <- list(
    model            = model_syntax,
    data             = spec$data,
    ordered          = spec$all_items,
    std.lv           = TRUE,
    parameterization = "theta",
    missing          = missing,
    rotation         = "target",
    rotation.args    = rot_args
  )
  if (!is.null(spec$group))  cfa_args$group       <- spec$group
  if (!is.null(group_equal)) cfa_args$group.equal <- group_equal
  cfa_args <- c(cfa_args, dots)

  fit <- .muffle_rotated_vcov(tryCatch(
    do.call(lavaan::cfa, cfa_args),
    error = function(e) {
      # Print the generated syntax so it can be inspected when debugging crashes.
      message("\n--- [DEBUG] B-ESEM model syntax that caused the error ---")
      message(model_syntax)
      message("--- [END DEBUG] ---")
      stop("B-ESEM efa() WLSMV failed: ", conditionMessage(e), call. = FALSE)
    }
  ))

  structure(
    list(
      lavaan_fit   = fit,
      syntax       = model_syntax,
      factor_names = all_factor_names,
      nfactors     = length(all_factor_names)
    ),
    class = c("besem_fit_ordered", "besem_fit", "esem_fit")
  )
}


.fit_esem_inv_ordered <- function(spec, group_equal, missing, ...) {
  # lavaan-style group.partial labels ("x9|t3", "x9~~x9") are honoured by the
  # threshold/residual writer below: the explicit equality labels written there
  # would otherwise defeat lavaan's own group.partial.
  dots <- list(...)
  group_partial <- dots$group.partial
  dots$group.partial <- NULL

  tgt    <- unclass(spec$target)        # strip esem_target S3 class
  tgt_na <- tgt
  tgt_na[tgt != 0] <- NA_real_

  all_factor_names <- colnames(tgt)
  lhs          <- paste(paste0("efa('esem')*", all_factor_names), collapse = " + ")
  model_syntax <- paste0(lhs, " =~ ", paste(spec$all_items, collapse = " + "))

  is_configural <- is.null(group_equal)
  # Levels 5/6 of esem_invariance(): latent var/cov and latent mean invariance.
  fix_varcov <- !is.null(group_equal) && "lv.variances" %in% group_equal
  fix_means  <- !is.null(group_equal) && "means"        %in% group_equal
  n_groups <- if (!is.null(spec$group_levels)) {
    length(spec$group_levels)
  } else if (!is.null(spec$group)) {
    length(unique(na.omit(spec$data[[spec$group]])))
  } else {
    1L
  }

  # Plain ESEM: oblique target rotation in every group, no orthogonality
  # override needed.  rstarts = 30L matches the B-ESEM helper for stable
  # reproducibility on complex multi-group fits.
  rot_args <- list(orthogonal = FALSE, target = tgt_na, rstarts = 30L)

  # -- Theta identification fix at weak (why not configural: see the same
  # block in .fit_besem_inv_ordered) --------------------------------------------
  needs_theta_id <- !is_configural && !is.null(group_equal) &&
    "loadings" %in% group_equal && !("thresholds" %in% group_equal)

  if (needs_theta_id && n_groups >= 2L) {

    # Referent items: first item of each specific factor (no G referent).
    refs2 <- unique(vapply(spec$factors, `[`, character(1), 1))
    refs2 <- refs2[!is.na(refs2)]

    grp_spec_mean <- paste(c("0", rep("NA", n_groups - 1L)), collapse = ", ")
    mean_lines <- vapply(all_factor_names, function(f)
      paste0(f, " ~ c(", grp_spec_mean, ")*1"), character(1L))
    model_syntax <- paste(c(model_syntax, mean_lines), collapse = "\n")

    res_vals <- paste(c("1", rep("NA", n_groups - 1L)), collapse = ", ")
    res_lines <- vapply(spec$all_items, function(it)
      paste0(it, " ~~ c(", res_vals, ")*", it), character(1L))
    model_syntax <- paste(c(model_syntax, res_lines), collapse = "\n")

    lbl_n <- 1L
    thresh_lines <- character(0)
    for (it in spec$all_items) {
      n_thr <- length(unique(na.omit(spec$data[[it]]))) - 1L
      if (n_thr == 0L) next
      n_constrained <- if (it %in% refs2) min(2L, n_thr) else 1L
      parts <- character(n_thr)
      for (k in seq_len(n_thr)) {
        if (k <= n_constrained) {
          lbl     <- paste0("wthr", lbl_n)
          grp_lbl <- paste(rep(lbl, n_groups), collapse = ", ")
          parts[k] <- paste0("c(", grp_lbl, ")*t", k)
          lbl_n   <- lbl_n + 1L
        } else {
          parts[k] <- paste0("t", k)
        }
      }
      thresh_lines <- c(thresh_lines,
        paste0(it, " | ", paste(parts, collapse = " + ")))
    }
    model_syntax <- paste(c(model_syntax, thresh_lines), collapse = "\n")
  }

  # -- Strong/Strict: full explicit threshold labelling + residual control -----
  if (!is_configural && !is.null(group_equal) && "thresholds" %in% group_equal) {

    # Means are written explicitly below; lv.variances / lv.covariances stay
    # for lavaan (they pin the non-reference unrotated Phi back to identity).
    has_residuals <- "residuals" %in% group_equal
    group_equal   <- setdiff(group_equal, c("thresholds", "residuals", "means"))

    pl <- .ordered_partial_lines(spec, n_groups, has_residuals, group_partial)
    model_syntax <- paste(c(model_syntax, pl$thresholds), collapse = "\n")

    # Latent mean invariance: the non-reference means go back to 0.
    grp_spec_mean <- paste(c("0", rep(if (fix_means) "0" else "NA", n_groups - 1L)),
                           collapse = ", ")
    mean_lines <- vapply(all_factor_names, function(f) {
      paste0(f, " ~ c(", grp_spec_mean, ")*1")
    }, character(1L))
    model_syntax <- paste(c(model_syntax, mean_lines), collapse = "\n")

    model_syntax <- paste(c(model_syntax, pl$residuals), collapse = "\n")
  }

  cfa_args <- list(
    model            = model_syntax,
    data             = spec$data,
    ordered          = spec$all_items,
    std.lv           = TRUE,
    parameterization = "theta",
    missing          = missing,
    rotation         = "target",
    rotation.args    = rot_args
  )
  if (!is.null(spec$group))  cfa_args$group       <- spec$group
  if (!is.null(group_equal)) cfa_args$group.equal <- group_equal
  cfa_args <- c(cfa_args, dots)

  fit <- .muffle_rotated_vcov(tryCatch(
    do.call(lavaan::cfa, cfa_args),
    error = function(e) {
      message("\n--- [DEBUG] ESEM model syntax that caused the error ---")
      message(model_syntax)
      message("--- [END DEBUG] ---")
      stop("ESEM efa() WLSMV failed: ", conditionMessage(e), call. = FALSE)
    }
  ))

  structure(
    list(
      lavaan_fit   = fit,
      syntax       = model_syntax,
      factor_names = all_factor_names,
      nfactors     = length(all_factor_names),
      indicators   = spec$all_items,
      estimator    = fit@Options$estimator
    ),
    class = c("esem_fit_ordered", "esem_fit")
  )
}


.fit_besem_inv_continuous <- function(spec, group_equal, missing, ...) {

  btgt    <- unclass(spec$bifactor_target)
  btgt_na <- btgt
  btgt_na[btgt != 0] <- NA_real_

  all_factor_names <- colnames(btgt)   # c("G", "EX", ...)
  lhs          <- paste(paste0('efa("besem")*', all_factor_names), collapse = " + ")
  model_syntax <- paste0(lhs, " =~ ", paste(spec$all_items, collapse = " + "))

  is_configural <- is.null(group_equal)
  # Levels 5/6 of esem_invariance(): latent var/cov and latent mean invariance.
  fix_varcov <- !is.null(group_equal) && "lv.variances" %in% group_equal
  fix_means  <- !is.null(group_equal) && "means"        %in% group_equal
  n_groups <- if (!is.null(spec$group_levels)) {
    length(spec$group_levels)
  } else if (!is.null(spec$group)) {
    length(unique(na.omit(spec$data[[spec$group]])))
  } else {
    1L
  }

  if (!is_configural && !is.null(spec$group)) {
    # Fix: Mplus's orthogonal rotation only enforces zero covariances in the
    # reference group. Non-reference group's k*(k-1)/2 covariances are free.
    # lavaan's orthogonal=TRUE alone over-constrains by fixing all groups (+6 df).
    # Fix: keep orthogonal=TRUE and add explicit c(0, NA)* for every factor pair.
    # The explicit NA for group 2+ overrides the rotation's global constraint.
    # NOTE: oblique rotation does NOT work — the rotation silently overrides c(0,NA)*.
    # Latent var/cov invariance: the non-reference covariances go back to 0
    # (explicit, so the group.equal token is dropped); latent means are left
    # to lavaan's "means" handling.
    grp_spec  <- paste(c("0", rep(if (fix_varcov) "0" else "NA", n_groups - 1L)),
                       collapse = ", ")
    fac_pairs <- combn(all_factor_names, 2L, simplify = FALSE)
    cov_lines <- vapply(fac_pairs,
      function(p) paste0(p[1L], " ~~ c(", grp_spec, ")*", p[2L]),
      character(1L))
    model_syntax <- paste(c(model_syntax, cov_lines), collapse = "\n")
    group_equal  <- setdiff(group_equal, "lv.covariances")
    rot_args <- list(orthogonal = TRUE, target = btgt_na)
  } else {
    rot_args <- list(orthogonal = TRUE, target = btgt_na)
  }

  cfa_args <- list(
    model         = model_syntax,
    data          = spec$data,
    std.lv        = TRUE,
    estimator     = if (!is.null(spec$estimator_esem)) spec$estimator_esem else "MLR",
    missing       = missing,
    rotation      = "target",
    rotation.args = rot_args
  )
  if (!is.null(spec$group))   cfa_args$group       <- spec$group
  if (!is.null(group_equal))  cfa_args$group.equal <- group_equal
  cfa_args <- c(cfa_args, list(...))

  fit <- .muffle_rotated_vcov(tryCatch(
    do.call(lavaan::cfa, cfa_args),
    error = function(e) {
      message("\n--- [DEBUG] B-ESEM MLR syntax that caused the error ---")
      message(model_syntax)
      message("--- [END DEBUG] ---")
      stop("B-ESEM MLR CFA failed: ", conditionMessage(e), call. = FALSE)
    }
  ))

  structure(
    list(
      lavaan_fit   = fit,
      syntax       = model_syntax,
      factor_names = all_factor_names,
      nfactors     = length(all_factor_names)
    ),
    class = c("besem_fit_continuous", "besem_fit", "esem_fit")
  )
}


.extract_fit_inv <- function(fit_obj) {
  na_row <- c(chisq = NA_real_, df = NA_real_, cfi = NA_real_,
              tli   = NA_real_, rmsea = NA_real_, rmsea_lo = NA_real_,
              rmsea_hi = NA_real_, srmr = NA_real_)

  if (is.null(fit_obj)) return(na_row)

  # Always use the raw lavaan S4 object -- fitMeasures() on the wrapper may not
  # return the full named vector needed for scaled-index detection.
  lav <- fit_obj$lavaan_fit
  fm  <- tryCatch(lavaan::fitMeasures(lav), error = function(e) NULL)
  if (is.null(fm)) return(na_row)

  # SRMR by the Mplus definition (.srmr_mplus()); lavaan's value as fallback.
  srmr_corrected <- tryCatch(.srmr_mplus(lav), error = function(e) unname(fm["srmr"]))

  # RMSEA 90% CI -- read from lavaan fitMeasures, fallback to .rmsea_ci().
  # `lo_keys`/`hi_keys` are character vectors so we try both lavaan spellings:
  # newer versions emit `rmsea.scaled.ci.lower`, older ones `rmsea.ci.lower.scaled`.
  # Without this, the lookup misses and the displayed CI comes from the naive
  # NCP inversion of chisq.scaled, which is a different formula than lavaan's
  # robust RMSEA -- the CI then no longer brackets the point estimate.
  .get_ci <- function(lo_keys, hi_keys, chisq_val, df_val) {
    .first <- function(keys) {
      hit <- intersect(keys, names(fm))
      if (length(hit)) as.numeric(fm[hit[1L]]) else NA_real_
    }
    lo <- .first(lo_keys); hi <- .first(hi_keys)
    if (is.na(lo) || is.na(hi)) {
      n_obs <- tryCatch(sum(lavaan::lavInspect(lav, "nobs")), error = function(e) NA_integer_)
      ci    <- .rmsea_ci(chisq_val, df_val, n_obs)
      lo <- ci[1L]; hi <- ci[2L]
    }
    c(lo, hi)
  }

  # Prefer scaled indices (MLR Satorra-Bentler; WLSMV mean-variance adjusted)
  if ("cfi.scaled" %in% names(fm) && !is.na(fm["cfi.scaled"])) {
    ci <- .get_ci(c("rmsea.ci.lower.scaled", "rmsea.scaled.ci.lower"),
                  c("rmsea.ci.upper.scaled", "rmsea.scaled.ci.upper"),
                  unname(fm["chisq.scaled"]), unname(fm["df.scaled"]))
    c(
      chisq    = unname(fm["chisq.scaled"]),
      df       = unname(fm["df.scaled"]),
      cfi      = unname(fm["cfi.scaled"]),
      tli      = unname(fm["tli.scaled"]),
      rmsea    = unname(fm["rmsea.scaled"]),
      rmsea_lo = ci[1L],
      rmsea_hi = ci[2L],
      srmr     = srmr_corrected
    )
  } else {
    ci <- .get_ci("rmsea.ci.lower", "rmsea.ci.upper",
                  unname(fm["chisq"]), unname(fm["df"]))
    c(
      chisq    = unname(fm["chisq"]),
      df       = unname(fm["df"]),
      cfi      = unname(fm["cfi"]),
      tli      = unname(fm["tli"]),
      rmsea    = unname(fm["rmsea"]),
      rmsea_lo = ci[1L],
      rmsea_hi = ci[2L],
      srmr     = srmr_corrected
    )
  }
}


.extract_lrt_row <- function(lrt) {
  if (is.null(lrt))
    return(c(delta_chisq = NA_real_, delta_df = NA_real_, pvalue = NA_real_))

  # lavTestLRT() returns a data.frame with class "anova".
  # Row 1 = less constrained; Row 2 = more constrained (the diff row).
  # Column names vary slightly across lavaan versions -- detect robustly.
  cn <- colnames(lrt)

  chisq_col <- .match_col(cn, c("Chisq diff", "Chisq.diff", "Scaled Chisq diff"))
  df_col    <- .match_col(cn, c("Df diff",    "Df.diff"))
  pval_col  <- .match_col(cn, c("Pr(>Chisq)", "Pr..Chisq."))

  c(
    delta_chisq = if (!is.na(chisq_col)) as.numeric(lrt[2L, chisq_col]) else NA_real_,
    delta_df    = if (!is.na(df_col))    as.numeric(lrt[2L, df_col])    else NA_real_,
    pvalue      = if (!is.na(pval_col))  as.numeric(lrt[2L, pval_col])  else NA_real_
  )
}


.match_col <- function(col_names, candidates) {
  hit <- intersect(candidates, col_names)
  if (length(hit)) hit[1L] else NA_character_
}


.build_invariance_table <- function(fits, lrt_results, levels, labels, is_ordered) {

  rows <- list()
  prev <- c(chisq = NA_real_, df = NA_real_, cfi = NA_real_, tli = NA_real_,
            rmsea = NA_real_, rmsea_lo = NA_real_, rmsea_hi = NA_real_, srmr = NA_real_)

  # Map from level name to the lrt_results key that compares it vs previous level
  lrt_key <- c(weak = "weak", strong = "strong", strict = "strict",
               varcov = "varcov", means = "means")

  for (lv in levels) {
    fi <- .extract_fit_inv(fits[[lv]])

    # Deltas computed from rounded (3 d.p.) values so that displayed D statistics
    # are consistent with the rounded fit indices shown in the table.
    # e.g. CFI 0.98524 -> 0.986, CFI 0.98521 -> 0.985, D = -0.001 not -0.00003.
    fi_r   <- round(fi, 3)
    prev_r <- round(prev, 3)
    .d <- function(key) {
      if (!is.na(prev_r[key]) && !is.na(fi_r[key])) unname(fi_r[key] - prev_r[key]) else NA_real_
    }

    delta_cfi   <- .d("cfi")
    delta_tli   <- .d("tli")
    delta_rmsea <- .d("rmsea")
    delta_srmr  <- .d("srmr")

    # DX2, Ddf and p from lavTestLRT() -- correct for both ML and WLSMV.
    # For WLSMV, lavTestLRT() applies the mean-variance adjusted difference test;
    # naive pchisq(DX2_scaled, Ddf) is not valid for scaled statistics.
    lrt_row     <- .extract_lrt_row(lrt_results[[lrt_key[lv]]])
    delta_chisq <- unname(lrt_row["delta_chisq"])
    delta_df    <- unname(lrt_row["delta_df"])
    pvalue      <- unname(lrt_row["pvalue"])

    rows[[lv]] <- data.frame(
      Model       = labels[lv],
      chisq       = unname(fi["chisq"]),
      df          = unname(fi["df"]),
      delta_chisq = delta_chisq,
      delta_df    = delta_df,
      pvalue      = pvalue,
      CFI         = unname(fi["cfi"]),
      dCFI        = delta_cfi,
      TLI         = unname(fi["tli"]),
      dTLI        = delta_tli,
      RMSEA       = unname(fi["rmsea"]),
      RMSEA_lo    = unname(fi["rmsea_lo"]),
      RMSEA_hi    = unname(fi["rmsea_hi"]),
      dRMSEA      = delta_rmsea,
      SRMR        = unname(fi["srmr"]),
      dSMR        = delta_srmr,
      row.names   = NULL,
      stringsAsFactors = FALSE
    )

    prev <- fi_r[c("chisq", "df", "cfi", "tli", "rmsea", "rmsea_lo", "rmsea_hi", "srmr")]
  }

  do.call(rbind, rows)
}


.print_invariance_table <- function(tbl, is_ordered) {

  .fmt <- function(x, digits = 3) {
    s <- ifelse(is.na(x), "", formatC(round(x, digits), format = "f", digits = digits))
    # formatC("-0.000") is confusing; collapse to "0.000"
    zero_pat <- paste0("^-0\\.", paste(rep("0", digits), collapse = ""), "$")
    gsub(zero_pat, paste0("0.", paste(rep("0", digits), collapse = "")), s)
  }
  .fmt_int <- function(x) ifelse(is.na(x), "", as.character(round(x, 0)))
  .fmt_p   <- function(x) {
    ifelse(is.na(x), "",
      ifelse(x < .001, "<.001", formatC(x, format = "f", digits = 3)))
  }

  display <- data.frame(
    "Model"      = tbl$Model,
    "chi2"       = .fmt(tbl$chisq),
    "df"         = .fmt_int(tbl$df),
    "dchi2"      = .fmt(tbl$delta_chisq),
    "ddf"        = .fmt_int(tbl$delta_df),
    "p"          = .fmt_p(tbl$pvalue),
    "CFI"        = .fmt(tbl$CFI),
    "dCFI"       = .fmt(tbl$dCFI),
    "TLI"        = .fmt(tbl$TLI),
    "dTLI"       = .fmt(tbl$dTLI),
    "RMSEA"      = mapply(function(r, lo, hi) {
      if (!is.na(r) && !is.na(lo) && !is.na(hi))
        sprintf("%.3f [%s, %s]", r,
                sub("^0", "", sprintf("%.3f", lo)),
                sub("^0", "", sprintf("%.3f", hi)))
      else .fmt(r)
    }, tbl$RMSEA, tbl$RMSEA_lo, tbl$RMSEA_hi),
    "dRMSEA"     = .fmt(tbl$dRMSEA),
    "SRMR"       = .fmt(tbl$SRMR),
    "dSRMR"      = .fmt(tbl$dSMR),
    stringsAsFactors = FALSE,
    check.names  = FALSE
  )
  names(display) <- c(
    "Model", "chi2(s)", "df",
    "dchi2", "ddf", "p",
    "CFI", "dCFI", "TLI", "dTLI",
    "RMSEA [90% CI]    ", "dRMSEA", "SRMR", "dSRMR"
  )

  print(display, row.names = FALSE)

  chisq_note <- if (is_ordered)
    "dchi2 = scaled difference test (lavTestLRT) -- NOT chi2_B minus chi2_A; WLSMV requires separate adjustment so dchi2 will not equal the naive difference. p is derived from dchi2 and ddf."
  else
    "dchi2 = scaled difference test (lavTestLRT) -- NOT chi2_B minus chi2_A; scaling corrections mean dchi2 will not equal the naive difference. p is derived from dchi2 and ddf."

  cat(
    "\nNotes:\n",
    " d = delta (more-constrained minus less-constrained).\n",
    " dCFI >= -0.010 suggests non-invariance (Cheung & Rensvold, 2002; Putnick & Bornstein, 2016).\n",
    " dRMSEA >= +0.015 suggests non-invariance (Chen, 2007; Putnick & Bornstein, 2016).\n",
    " p = significance of the scaled chi-square difference test (lavTestLRT).\n",
    "   A significant p (< .05) indicates the more-constrained model fits significantly\n",
    "   worse, but practical criteria (dCFI, dRMSEA) should take precedence over p alone.\n",
    " chi2(s) = ", if (is_ordered) "WLSMV mean-variance adjusted statistic. " else "Satorra-Bentler scaled chi-square (MLR). ",
    chisq_note, "\n\n",
    sep = ""
  )
}


# -- Mplus comparison ----------------------------------------------------------

#' Generate, Run, and Compare B-ESEM Invariance Models Against Mplus
#'
#' Creates complete Mplus \code{.inp} files for every level fitted in
#' \code{inv} (configural, weak, strong, strict and, when
#' \code{esem_invariance()} was run with \code{through}, latent var/cov and
#' latent means), runs them via \pkg{MplusAutomation}, then prints a
#' side-by-side comparison of fit statistics against the R results from
#' \code{\link{esem_invariance}}.
#'
#' @param inv An \code{esem_invariance} object (from \code{model = "besem"}).
#' @param output_folder Character. Folder where \code{data.dat}, \code{.inp},
#'   and \code{.out} files will be written. Created if it does not exist.
#' @param mplus_command Character. Full path to the Mplus executable.
#'   Default \code{"mplus"} (works if Mplus is on the system PATH).
#' @param group_labels Named character vector mapping group values to Mplus
#'   labels, e.g. \code{c("1" = "MALE", "2" = "FEMALE")}. If \code{NULL}
#'   (default) labels are auto-derived as \code{G1}, \code{G2}, ...
#' @param missing_code Numeric. Missing value sentinel written to the data
#'   file. Default \code{-999}.
#' @param difftest Logical. Generate \code{DIFFTEST} constraint files for
#'   chained nested-model comparisons (strong vs weak, strict vs strong).
#'   Default \code{TRUE}.
#'
#' @return A data frame (invisibly) with \code{_R}, \code{_Mplus}, and
#'   \code{delta_} columns for each fit statistic. Printed as a table.
#'
#' @seealso \code{\link{esem_invariance}}
#'
#' @examples
#' \dontrun{
#' # Requires a licensed Mplus installation reachable via `mplus_command`.
#' inv <- esem_invariance(spec, model = "besem")
#' cmp <- run_mplus_besem_invariance(
#'   inv,
#'   output_folder = tempfile("mplus_inv_"),
#'   mplus_command = "C:/Program Files/Mplus/mplus.exe",
#'   group_labels  = c("1" = "MALE", "2" = "FEMALE")
#' )
#' }
#'
#' @export
run_mplus_besem_invariance <- function(inv,
                                        output_folder,
                                        mplus_command = "mplus",
                                        group_labels  = NULL,
                                        missing_code  = -999,
                                        difftest      = TRUE) {

  if (!requireNamespace("MplusAutomation", quietly = TRUE))
    stop("MplusAutomation is required. Install: install.packages('MplusAutomation')",
         call. = FALSE)
  if (!inherits(inv, "esem_invariance"))
    stop("`inv` must be an esem_invariance object.", call. = FALSE)
  if (is.null(inv$model) || inv$model != "besem")
    stop("This function is for B-ESEM invariance (model = 'besem').", call. = FALSE)

  spec   <- inv$spec
  folder <- normalizePath(output_folder, mustWork = FALSE)
  dir.create(folder, showWarnings = FALSE, recursive = TRUE)

  # -- Sanitize variable names for Mplus (dots -> underscores) -----------------
  .mp  <- function(x) gsub(".", "_", x, fixed = TRUE)
  items   <- spec$all_items
  mp_items <- .mp(items)
  mp_group <- .mp(spec$group)

  # -- Group labels ------------------------------------------------------------
  g_levels <- spec$group_levels
  if (is.null(group_labels)) {
    group_labels <- setNames(paste0("G", seq_along(g_levels)), as.character(g_levels))
  }
  grouping_str <- paste(
    vapply(as.character(g_levels), function(v)
      paste0(v, " = ", group_labels[v]), character(1)),
    collapse = " "
  )
  # Non-reference group labels -- one MODEL <label>: block is generated per
  # non-reference group.  For K groups, K-1 blocks are produced.  Must match
  # exactly what is declared in the GROUPING block above.
  nonref_labels <- unname(group_labels[as.character(g_levels[-1])])

  # -- Write data file (items + group, missing -> missing_code) -----------------
  dat_path <- file.path(folder, "data.dat")
  df_out   <- spec$data[, c(items, spec$group), drop = FALSE]
  df_out[is.na(df_out)] <- missing_code
  # prepareMplusData writes a headerless space-delimited file
  MplusAutomation::prepareMplusData(df_out, filename = dat_path, overwrite = TRUE,
                                     inpfile = FALSE)
  message(paste("  Data written:", dat_path))

  is_ordered <- !is.null(spec$ordered) && length(spec$ordered) > 0

  # -- Number of thresholds per item (max unique levels - 1) -------------------
  n_thresh <- max(vapply(items, function(v) {
    length(unique(na.omit(spec$data[[v]]))) - 1L
  }, integer(1)))

  # -- Determine referent indicators -------------------------------------------
  # EX referent: items[1];  G referent: items[2]
  # Each subsequent specific factor: first item of that factor
  factor_names_spec <- names(spec$factors)   # e.g. c("EX","MD","CI")
  refs2 <- character(0)  # items that get BOTH $1 and $2 labeled (partial model)
  for (fi in seq_along(factor_names_spec)) {
    refs2 <- c(refs2, spec$factors[[fi]][1])
  }
  refs2 <- c(refs2, spec$factors[[1]][2])  # G referent = 2nd item of 1st factor
  refs2 <- unique(refs2[!is.na(refs2)])

  # -- BY block (shared across all levels) -------------------------------------
  by_block <- .besem_by_lines(spec, mp_items, factor_names_spec)

  # -- Generate and write one .inp per level -----------------------------------
  # Mirror the levels fitted in R (through = "varcov" / "means" adds two).
  levels    <- names(inv$models)
  inp_paths <- setNames(file.path(folder, paste0("besem_inv_", levels, ".inp")), levels)
  diff_files <- setNames(c(NA, paste0("besem_inv_", levels[-length(levels)], ".dat")),
                         levels)
  save_files <- paste0("besem_inv_", levels, ".dat")

  # Build wrapped variable list lines (Mplus 90-char limit; use 85 to be safe)
  .mplus_varline <- function(keyword, vars) {
    prefix  <- paste0("  ", keyword, " ARE ")
    indent  <- paste(rep(" ", nchar(prefix)), collapse = "")
    lines   <- character(0)
    current <- prefix
    for (v in vars) {
      candidate <- paste0(current, v, " ")
      if (nchar(candidate) > 85 && current != prefix) {
        lines   <- c(lines, current)
        current <- paste0(indent, v, " ")
      } else {
        current <- candidate
      }
    }
    lines <- c(lines, paste0(trimws(current, "right"), ";"))
    paste(lines, collapse = "\n")
  }

  names_line <- .mplus_varline("NAMES",        c(mp_items, mp_group))
  use_line   <- .mplus_varline("USEVARIABLES", c(mp_items, mp_group))
  miss_line  <- paste0("  MISSING ARE ALL (", missing_code, ");")
  group_line <- paste0("  GROUPING IS ", mp_group, " (", grouping_str, ");")
  output_sec <- "OUTPUT:\n  SAMPSTAT STANDARDIZED STDYX TECH1 TECH4;"

  # Estimator-specific VARIABLE and ANALYSIS settings
  if (is_ordered) {
    cat_line <- .mplus_varline("CATEGORICAL", mp_items)
    analysis <- paste(
      "ANALYSIS:",
      "  ESTIMATOR = WLSMV;",
      "  ROTATION = TARGET(ORTHOGONAL);",
      "  PARAMETERIZATION = THETA;",
      sep = "\n"
    )
  } else {
    cat_line <- ""   # MLR: continuous items, no CATEGORICAL statement
    ml_est <- if (!is.null(spec$estimator_esem) && nchar(spec$estimator_esem) > 0)
                spec$estimator_esem else "MLR"
    analysis <- paste(
      "ANALYSIS:",
      paste0("  ESTIMATOR = ", ml_est, ";"),
      "  ROTATION = TARGET(ORTHOGONAL);",
      sep = "\n"
    )
  }

  # -- Generate .inp, run, and retry on non-convergence (per level) -------------
  # WLSMV: run in order so DIFFTEST chain works.
  # ML: DIFFTEST is WLSMV-only; chi-square differences are computed by subtraction.
  message(sprintf("\n  Running Mplus (%d models)...", length(levels)))
  for (lv in levels) {
    message(sprintf("    %-12s ...", lv))

    full_thresh  <- lv %in% c("strong", "strict", "varcov", "means")
    strict_model <- lv %in% c("strict", "varcov", "means")
    fix_varcov   <- lv %in% c("varcov", "means")   # non-ref variances @1, covariances @0
    fix_means    <- lv == "means"                  # non-ref factor means @0
    configural   <- lv == "configural"

    fac_all      <- c("G", factor_names_spec)
    fac_var_ref  <- paste(paste0("  ", fac_all, "@1;"), collapse = "\n")
    fac_mean_ref <- paste(paste0("  [", fac_all, "@0];"), collapse = "\n")
    fac_var_fem  <- paste(paste0("  ", fac_all, "*;"),   collapse = "\n")
    fac_mean_fem <- paste(paste0("  [", fac_all, "*];"), collapse = "\n")
    if (fix_varcov) {
      # Mplus rejects "@0" (and bare "@1") on EFA-set parameters in a
      # non-reference group; the technical-supplement form is equality labels
      # on every covariance in both groups plus variances fixed at 1, which
      # equates the non-reference covariances to the reference group's
      # rotation-fixed zeros. (fac_var_ref is never used for a non-reference
      # block at this level: fix_varcov and configural are exclusive.)
      prs     <- combn(fac_all, 2L)
      cov_lab <- sprintf("  %s WITH %s (c%d);", prs[1L, ], prs[2L, ], seq_len(ncol(prs)))
      fac_var_ref <- paste(c(fac_var_ref, cov_lab), collapse = "\n")
      fac_var_fem <- paste(c(paste0("  ", fac_all, "@1;"), cov_lab), collapse = "\n")
    }
    if (fix_means) fac_mean_fem <- fac_mean_ref

    if (is_ordered) {
      # -- WLSMV path (original logic) -----------------------------------------
      thresh_ref <- .besem_thresh_lines(mp_items, n_thresh, refs2, full_thresh, grp = "ref",
                                        configural = configural)
      thresh_fem <- .besem_thresh_lines(mp_items, n_thresh, refs2, full_thresh, grp = "fem",
                                        configural = configural)

      uniq_ref <- paste0("  ", mp_items[1], "-", mp_items[length(mp_items)], "@1;")
      uniq_fem <- if (strict_model)
                    paste0("  ", mp_items[1], "-", mp_items[length(mp_items)], "@1;")
                  else
                    paste0("  ", mp_items[1], "-", mp_items[length(mp_items)], "*;")

      # One MODEL <label>: block per non-reference group; same body in each
      # (constraints are identical across non-reference groups for measurement
      # invariance — they differ only in which group's parameters they bind).
      make_nonref_block <- function(g_label) {
        if (configural) {
          # Standard configural identification in every group: factor
          # variances 1, factor means 0, residuals 1, all thresholds free.
          paste(c(paste0("MODEL ", g_label, ":"), by_block, fac_var_ref, fac_mean_ref,
                  thresh_fem, uniq_ref), collapse = "\n")
        } else {
          paste(c(paste0("MODEL ", g_label, ":"), fac_var_fem, fac_mean_fem,
                  thresh_fem, uniq_fem), collapse = "\n")
        }
      }
      model_female <- paste(vapply(nonref_labels, make_nonref_block, character(1)),
                            collapse = "\n\n")

      make_inp <- function(conv_crit = NULL) {
        analysis_str <- paste0(analysis, .mplus_conv_lines(conv_crit))
        analysis_with_diff <- if (difftest && !is.na(diff_files[lv]))
          paste0(analysis_str,
                 .mplus_difftest_line(file.path(dirname(inp_paths[lv]), diff_files[lv])))
        else
          analysis_str
        var_lines <- Filter(nchar, c(names_line, use_line, miss_line, group_line, cat_line))
        inp_parts <- c(
          paste0("TITLE: B-ESEM Invariance - ", lv, " (", spec$label, ");"),
          "", "DATA:", "  FILE = \"data.dat\";", "",
          "VARIABLE:", var_lines, "",
          analysis_with_diff, "",
          "MODEL:", by_block, fac_var_ref, fac_mean_ref, thresh_ref, uniq_ref, "",
          model_female, "",
          output_sec
        )
        if (difftest)
          inp_parts <- c(inp_parts, "", paste0("SAVEDATA:\n  DIFFTEST = besem_inv_", lv, ".dat;"))
        paste(inp_parts, collapse = "\n")
      }

    } else {
      # -- MLR path -------------------------------------------------------------
      # Pattern from besem-proof/ML/:
      #   Configural  -- separate BY blocks in MODEL FEMALE; fac var @1 both groups;
      #                 fac means @0 both groups; intercepts free (unlabeled);
      #                 uniquenesses free (unlabeled).
      #   Weak        -- no BY in MODEL FEMALE (loadings equal); fac var * in G2;
      #                 fac means @0 in G2 (still fixed!); intercepts/uniq free.
      #   Strong      -- fac var * + fac means * in G2; intercepts equal (labeled);
      #                 uniquenesses free (unlabeled).
      #   Strict      -- same as strong + uniquenesses equal (labeled).

      # Per-item line builders
      make_int_lines <- function(labeled) {
        if (labeled)
          paste(sprintf("  [%s] (i%d);", mp_items, seq_along(mp_items)), collapse = "\n")
        else
          paste(sprintf("  [%s];", mp_items), collapse = "\n")
      }
      make_uniq_lines <- function(labeled) {
        if (labeled)
          paste(sprintf("  %s (u%d);", mp_items, seq_along(mp_items)), collapse = "\n")
        else
          paste(sprintf("  %s;", mp_items), collapse = "\n")
      }

      int_lines_ref  <- make_int_lines(full_thresh)   # labeled at strong/strict
      int_lines_fem  <- make_int_lines(full_thresh)   # same labels -> equal
      uniq_lines_ref <- make_uniq_lines(strict_model) # labeled at strict
      uniq_lines_fem <- make_uniq_lines(strict_model) # same labels -> equal

      # Factor variances in MODEL FEMALE
      fac_var_fem_lv <- if (configural) fac_var_ref else fac_var_fem
      # Factor means in MODEL FEMALE: @0 at configural/weak; free at strong/strict
      fac_mean_fem_lv <- if (full_thresh) fac_mean_fem else fac_mean_ref

      # One MODEL <label>: block per non-reference group; same body in each.
      make_nonref_block <- function(g_label) {
        fem_parts <- c(paste0("MODEL ", g_label, ":"), "")
        if (configural) fem_parts <- c(fem_parts, by_block, "")
        fem_parts <- c(fem_parts,
                       fac_var_fem_lv, "",
                       fac_mean_fem_lv, "",
                       int_lines_fem, "",
                       uniq_lines_fem)
        paste(fem_parts, collapse = "\n")
      }
      model_female <- paste(vapply(nonref_labels, make_nonref_block, character(1)),
                            collapse = "\n\n")

      make_inp <- function(conv_crit = NULL) {
        analysis_str <- paste0(analysis, .mplus_conv_lines(conv_crit))
        var_lines <- Filter(nchar, c(names_line, use_line, miss_line, group_line))
        paste(c(
          paste0("TITLE: B-ESEM Invariance - ", lv, " (", spec$label, ");"),
          "", "DATA:", "  FILE = \"data.dat\";", "",
          "VARIABLE:", var_lines, "",
          analysis_str, "",
          "MODEL:", by_block, "",
          fac_var_ref, "",
          fac_mean_ref, "",
          int_lines_ref, "",
          uniq_lines_ref, "",
          model_female, "",
          output_sec
        ), collapse = "\n")
      }
    }

    out_path <- sub("\\.inp$", ".out", inp_paths[lv])

    # Attempt 0: CONVERGENCE = 1e-6 (.mplus_conv_lines default)
    writeLines(make_inp(), inp_paths[lv])
    MplusAutomation::runModels(target = inp_paths[lv], Mplus_command = mplus_command,
                                showOutput = FALSE, replaceOutfile = "always")

    if (.mplus_converged(out_path)) { message("  OK"); next }

    # If .out has a hard ERROR (syntax/spec), retrying won't help
    if (file.exists(out_path) &&
        any(grepl("\\*\\*\\* ERROR", readLines(out_path, warn = FALSE)))) {
      message("  ERROR"); next
    }

    # Retries: Mplus default, then progressively relaxed (mirrors R .fit_with_retry)
    converged_mplus <- FALSE
    for (i in seq_along(.MPLUS_CONV_RETRY)) {
      conv <- .MPLUS_CONV_RETRY[i]
      message(sprintf("\n      [retry %d/%d conv=%s] ...", i, length(.MPLUS_CONV_RETRY),
                      format(conv, scientific = FALSE)))
      writeLines(make_inp(conv), inp_paths[lv])
      MplusAutomation::runModels(target = inp_paths[lv], Mplus_command = mplus_command,
                                  showOutput = FALSE, replaceOutfile = "always")
      if (.mplus_converged(out_path)) {
        message("  OK")
        converged_mplus <- TRUE
        break
      }
    }
    if (!converged_mplus) message("  FAILED (all retries exhausted)")
  }

  # -- Read .out files ---------------------------------------------------------
  out_paths <- sub("\\.inp$", ".out", inp_paths)
  mplus_fi  <- lapply(setNames(levels, levels), function(lv) {
    p <- out_paths[lv]
    if (!file.exists(p)) return(NULL)
    tryCatch(MplusAutomation::readModels(p, quiet = TRUE),
             error = function(e) NULL)
  })

  # -- Build and print comparison table ----------------------------------------
  r_tbl <- inv$table
  rows  <- lapply(seq_along(levels), function(i) {
    mp  <- mplus_fi[[levels[i]]]
    r_chi <- r_tbl$chisq[i]; r_df  <- r_tbl$df[i]
    r_cfi <- r_tbl$CFI[i];   r_tli <- r_tbl$TLI[i]
    r_rms <- r_tbl$RMSEA[i]; r_srm <- r_tbl$SRMR[i]
    if (is.null(mp)) {
      m_chi <- m_df <- m_cfi <- m_tli <- m_rms <- m_srm <- NA_real_
    } else {
      s     <- mp$summaries
      m_chi <- .mp_val(s, c("ChiSqM_Value", "Chi2Value"))
      m_df  <- .mp_val(s, c("ChiSqM_DF",    "Chi2DF"))
      m_cfi <- .mp_val(s, c("CFI"))
      m_tli <- .mp_val(s, c("TLI"))
      m_rms <- .mp_val(s, c("RMSEA_Estimate"))
      m_srm <- .mp_val(s, c("SRMR", "SRMSR"))
    }
    data.frame(Level = r_tbl$Model[i],
      chi2_R = r_chi, chi2_M = m_chi, dchi2 = r_chi - m_chi,
      df_R   = r_df,  df_M   = m_df,  ddf   = r_df  - m_df,
      CFI_R  = r_cfi, CFI_M  = m_cfi, dCFI  = r_cfi - m_cfi,
      TLI_R  = r_tli, TLI_M  = m_tli, dTLI  = r_tli - m_tli,
      RMS_R  = r_rms, RMS_M  = m_rms, dRMS  = r_rms - m_rms,
      SRM_R  = r_srm, SRM_M  = m_srm, dSRM  = r_srm - m_srm,
      row.names = NULL, stringsAsFactors = FALSE)
  })
  cmp <- do.call(rbind, rows)

  message("\n==============================================================")
  message(" B-ESEM Invariance: R vs Mplus (d = R - Mplus)")
  message("==============================================================\n")

  .f3 <- function(x) ifelse(is.na(x), "   NA  ", formatC(round(x,3), format="f", digits=3, width=7))
  .fi <- function(x) ifelse(is.na(x), "  NA ", formatC(round(x,0),  format="d", width=5))

  message(sprintf("%-20s %7s %7s %7s  %5s %5s %5s  %6s %6s %7s  %6s %6s %7s  %6s %6s %7s  %6s %6s %7s",
    "Level","chi2_R","chi2_M","dchi2","df_R","df_M","ddf",
    "CFI_R","CFI_M","dCFI","TLI_R","TLI_M","dTLI",
    "RMS_R","RMS_M","dRMS","SRM_R","SRM_M","dSRM"))
  message(strrep("-", 142))
  for (i in seq_len(nrow(cmp))) {
    r <- cmp[i, ]
    message(sprintf("%-20s %7s %7s %7s  %5s %5s %5s  %6s %6s %7s  %6s %6s %7s  %6s %6s %7s  %6s %6s %7s",
      r$Level,
      .f3(r$chi2_R), .f3(r$chi2_M), .f3(r$dchi2),
      .fi(r$df_R),   .fi(r$df_M),   .fi(r$ddf),
      .f3(r$CFI_R),  .f3(r$CFI_M),  .f3(r$dCFI),
      .f3(r$TLI_R),  .f3(r$TLI_M),  .f3(r$dTLI),
      .f3(r$RMS_R),  .f3(r$RMS_M),  .f3(r$dRMS),
      .f3(r$SRM_R),  .f3(r$SRM_M),  .f3(r$dSRM)))
  }
  message("\n  d = R - Mplus. Target: |d| <= 0.001\n")

  invisible(cmp)
}


# -- Syntax helpers for B-ESEM invariance .inp generation ---------------------

.besem_by_lines <- function(spec, mp_items, factor_names_spec) {
  lines <- character(0)

  # Mplus requires (*N) labels to appear on the same line as some loading;
  # standalone "(*1);" is rejected with "should be declared on the same line
  # with some of the factor loadings".  We therefore append " (*1);" to the
  # final item token before wrapping, guaranteeing it stays attached.
  attach_label <- function(tokens, label = "(*1);") {
    n <- length(tokens)
    if (n == 0) return(label)
    tokens[n] <- paste0(tokens[n], " ", label)
    tokens
  }

  # G BY: all items load freely -- no cross-loading targets (~0)
  lines <- c(lines, "  G BY", .wrap85(attach_label(mp_items)))

  # Specific factors: primary items first, then all others targeted toward 0
  for (fac in factor_names_spec) {
    primary   <- .mp_names(spec$factors[[fac]])
    cross     <- setdiff(mp_items, primary)
    cross_tgt <- paste0(cross, "~0")   # e.g. "batMD1~0" -- no space, single token
    lines <- c(lines, paste0("  ", fac, " BY"),
               .wrap85(attach_label(c(primary, cross_tgt))))
  }

  lines
}

.wrap85 <- function(tokens, indent = "    ") {
  lines   <- character(0)
  current <- indent
  for (tok in tokens) {
    candidate <- paste0(current, tok, " ")
    if (nchar(candidate) > 85 && current != indent) {
      lines   <- c(lines, trimws(current, "right"))
      current <- paste0(indent, tok, " ")
    } else {
      current <- candidate
    }
  }
  if (nchar(trimws(current)) > 0)
    lines <- c(lines, trimws(current, "right"))
  lines
}

.mp_names <- function(x) gsub(".", "_", x, fixed = TRUE)

.besem_thresh_lines <- function(mp_items, n_thresh, refs2, full_equal, grp,
                                configural = FALSE) {
  lines  <- character(0)
  t_num  <- 1L
  for (it in mp_items) {
    for (k in seq_len(n_thresh)) {
      lbl  <- paste0("t", t_num)
      tag  <- paste0("[", it, "$", k, "]")
      # Decide whether this threshold gets a label (none at configural: the
      # standard identification keeps every threshold free in every group)
      labeled <- !configural &&
                 (full_equal ||
                  k == 1L ||
                  (k == 2L && it %in% refs2))
      line <- if (labeled) paste0("  ", tag, " (", lbl, ");")
              else          paste0("  ", tag, ";")
      lines  <- c(lines, line)
      t_num  <- t_num + 1L
    }
  }
  lines
}


.mp_val <- function(s, candidates) {
  for (nm in candidates) {
    v <- s[[nm]]
    if (!is.null(v) && length(v) == 1L && is.numeric(v) && !is.na(v))
      return(as.numeric(v))
  }
  NA_real_
}


# Internal: threshold and residual lines of the strong/strict ordered models.
# All thresholds are equal across groups through explicit labels (the part
# after "*" must stay the positional t1, t2, ...); lavaan-style group.partial
# labels release them: "x9|t3" leaves that threshold free in every group,
# "x9~~x9" (strict) frees the item's residual in the non-reference groups.
# An item whose thresholds are all released gets its residual fixed to 1 in
# every group (theta identification).
.ordered_partial_lines <- function(spec, n_groups, has_residuals, group_partial = NULL) {
  thresh_lines <- character(0)
  res_lines    <- character(0)
  lbl_n <- 1L
  for (it in spec$all_items) {
    n_thr <- length(unique(na.omit(spec$data[[it]]))) - 1L
    if (n_thr == 0L) next
    parts <- character(n_thr)
    freed <- paste0(it, "|t", seq_len(n_thr)) %in% group_partial
    for (k in seq_len(n_thr)) {
      lbl     <- paste0("ethr", lbl_n)
      grp_lbl <- paste(rep(lbl, n_groups), collapse = ", ")
      parts[k] <- if (freed[k]) paste0("t", k) else paste0("c(", grp_lbl, ")*t", k)
      lbl_n <- lbl_n + 1L
    }
    thresh_lines <- c(thresh_lines, paste0(it, " | ", paste(parts, collapse = " + ")))
    res_free <- paste0(it, "~~", it) %in% group_partial
    non_ref  <- if (all(freed) || (has_residuals && !res_free)) "1" else "NA"
    res_lines <- c(res_lines,
      paste0(it, " ~~ c(", paste(c("1", rep(non_ref, n_groups - 1L)), collapse = ", "), ")*", it))
  }
  list(thresholds = thresh_lines, residuals = res_lines)
}
