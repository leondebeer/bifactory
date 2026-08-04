#' Specify an ESEM Model Structure
#'
#' Creates a model specification object that flows automatically through
#' all subsequent steps -- CFA, ESEM, B-ESEM, alignment check, Mplus syntax,
#' and comparison tables. Define your factor structure once at the top;
#' everything else is derived automatically.
#'
#' @param ... Named character vectors, one per factor. The name becomes the
#'   factor name and the vector contains the indicator names.
#'   Example: \code{EX = c("y1","y2","y3"), MD = c("y4","y5","y6")}.
#' @param data A \code{data.frame} containing the indicators.
#' @param label Optional character string labelling the model (used in
#'   output headers and Mplus titles). Default \code{"ESEM Model"}.
#'
#' @return An object of class \code{"esem_spec"} containing:
#' \describe{
#'   \item{\code{factors}}{Named list of factor -> indicator assignments.}
#'   \item{\code{factor_names}}{Character vector of factor names.}
#'   \item{\code{all_items}}{Character vector of all indicators in order.}
#'   \item{\code{nfactors}}{Number of specific factors.}
#'   \item{\code{data}}{The supplied data frame.}
#'   \item{\code{label}}{Model label.}
#'   \item{\code{cfa_syntax}}{Ready-to-use lavaan CFA model string.}
#'   \item{\code{target}}{Target matrix for ESEM target rotation.}
#'   \item{\code{bifactor_target}}{Target matrix for B-ESEM.}
#' }
#'
#' @param ordered Logical or character vector. If \code{TRUE}, all indicators
#'   are treated as ordered-categorical and the estimator is automatically
#'   switched to \code{"WLSMV"} (both R and Mplus). If a character vector of
#'   item names is supplied, only those items are treated as ordered.
#'   Default \code{FALSE} (continuous).
#' @param group Character. Name of a grouping variable in \code{data} for
#'   multi-group models. When supplied, a configural model is fitted by
#'   default. Use \code{group_equal} in \code{run_comparison()} to test
#'   metric or scalar invariance.
#' @param estimator Character. Override the auto-selected lavaan estimator.
#'   If \code{NULL} (default), uses \code{"WLSMV"} when \code{ordered} is
#'   supplied and \code{"MLR"} otherwise.
#' @param missing Character. Missing-data handling. If \code{NULL} (default),
#'   \code{"pairwise"} for ordered data and \code{"listwise"} for continuous.
#'   Ordered accepts \code{"pairwise"} or \code{"listwise"}. Continuous accepts
#'   \code{"listwise"} or any FIML alias (\code{"fiml"}, \code{"ml"},
#'   \code{"direct"}).  The chosen method is applied uniformly to CFA, ESEM,
#'   and B-ESEM so that fit statistics are computed on the same sample (fair
#'   comparison).  Note that the FIML + GPArotation::targetT path used for
#'   B-ESEM with 4+ specific factors can settle into a local optimum (G
#'   absorbs specific-factor variance); a warning is emitted when that
#'   combination is detected.
#'
#' @examples
#' data("HolzingerSwineford1939", package = "lavaan")
#'
#' # Continuous indicators (MLR estimator by default)
#' spec <- specify_model(
#'   Visual  = c("x1", "x2", "x3"),
#'   Textual = c("x4", "x5", "x6"),
#'   Speed   = c("x7", "x8", "x9"),
#'   data    = HolzingerSwineford1939,
#'   label   = "Holzinger-Swineford 3-factor"
#' )
#' print(spec)
#'
#' # Multi-group spec (for invariance testing)
#' spec_mg <- specify_model(
#'   Visual  = c("x1", "x2", "x3"),
#'   Textual = c("x4", "x5", "x6"),
#'   Speed   = c("x7", "x8", "x9"),
#'   data    = HolzingerSwineford1939,
#'   group   = "sex",
#'   label   = "HS 3-factor (multi-group)"
#' )
#'
#' \donttest{
#' # Run the CFA/ESEM/B-ESEM comparison from a spec
#' results <- run_comparison(spec)
#' print(results)
#' }
#'
#' @export
specify_model <- function(..., data, label = "ESEM Model",
                           ordered = FALSE, group = NULL, estimator = NULL,
                           missing = NULL) {

  factors <- list(...)

  if (length(factors) == 0)
    stop("Supply at least one named factor, e.g. EX = c('y1','y2','y3').",
         call. = FALSE)

  if (is.null(names(factors)) || any(names(factors) == ""))
    stop("All factors must be named.", call. = FALSE)

  if (missing(data))
    stop("Supply a data frame via the `data` argument.", call. = FALSE)

  factor_names <- names(factors)
  all_items    <- unlist(factors, use.names = FALSE)

  # Validate items exist in data
  missing_items <- setdiff(all_items, colnames(data))
  if (length(missing_items))
    stop("Items not found in data: ",
         paste(missing_items, collapse = ", "), call. = FALSE)

  # Mplus truncates variable names to 8 chars (case-insensitive; dots -> _).
  # Warn if two factor names would alias after truncation -- silent collisions
  # in Mplus output are hard to diagnose later. Harmless for R-only workflows.
  mplus_keys <- substr(tolower(gsub(".", "_", factor_names, fixed = TRUE)), 1L, 8L)
  if (anyDuplicated(mplus_keys)) {
    coll     <- unique(mplus_keys[duplicated(mplus_keys)])
    culprits <- vapply(coll, function(k)
      paste(factor_names[mplus_keys == k], collapse = " / "),
      character(1))
    warning(
      "Factor names share an 8-character prefix and will collide in Mplus ",
      "output:\n  ", paste(culprits, collapse = "\n  "),
      "\nRename to unique 8-char prefixes if you plan to use mplus_folder.",
      call. = FALSE
    )
  }

  # Auto-build CFA syntax
  cfa_lines <- vapply(factor_names, function(f) {
    paste0("  ", f, " =~ ", paste(factors[[f]], collapse = " + "))
  }, FUN.VALUE = character(1))
  cfa_syntax <- paste(cfa_lines, collapse = "\n")

  # Auto-build ESEM target matrix
  tgt <- make_target(
    keys       = factors,
    item_names = all_items
  )

  # Auto-build B-ESEM target matrix
  btgt <- make_bifactor_target(
    specific_factors = factors,
    indicators       = all_items
  )

  # Resolve ordered argument
  ordered_items <- if (isTRUE(ordered)) {
    all_items   # all items ordered
  } else if (is.character(ordered) && length(ordered) > 0) {
    bad <- setdiff(ordered, all_items)
    if (length(bad))
      warning("ordered items not in indicators, ignoring: ",
              paste(bad, collapse = ", "), call. = FALSE)
    intersect(ordered, all_items)
  } else {
    NULL
  }

  # CFA: WLSMV + ordered items natively in lavaan.
  # ESEM/BESEM ordered: routed to esem_ordered() / besem_ordered() (set-ESEM /
  # custom DWLS paths), not the continuous efa()+MLR path in esem()/besem().
  # Continuous data: MLR default.
  estimator_cfa  <- if (!is.null(ordered_items)) "WLSMV" else "MLR"
  estimator_esem <- if (!is.null(ordered_items)) "WLSMV" else "MLR"
  if (!is.null(estimator)) {
    estimator_cfa  <- estimator
    estimator_esem <- estimator
  }

  # Resolve `missing`:
  #   ordered     -> "pairwise" (all available pairs for polychoric) [default]
  #                  "listwise" is also valid if complete-case is required.
  #   continuous  -> "listwise" (complete case) [default]
  #                  "fiml" / "ml" / "direct" are accepted (CFA + ESEM only;
  #                  B-ESEM internally forces listwise to avoid a lavaan EFA
  #                  local-optimum trap documented in bifactor.R:183-185).
  fiml_syn <- c("fiml", "ml", "ml.x", "direct", "direct.ml")
  if (is.null(missing)) {
    missing <- if (!is.null(ordered_items)) "pairwise" else "listwise"
  } else {
    missing <- tolower(missing)
    if (!is.null(ordered_items)) {
      if (!missing %in% c("pairwise", "listwise"))
        stop("For ordered data, `missing` must be \"pairwise\" or \"listwise\" ",
             "(you gave \"", missing, "\").", call. = FALSE)
    } else {
      if (missing %in% fiml_syn) {
        missing <- "ml"  # lavaan canonical name for FIML
      } else if (!missing %in% c("listwise")) {
        stop("For continuous data, `missing` must be \"listwise\" or one of ",
             paste0('"', fiml_syn, '"', collapse = ", "),
             " (you gave \"", missing, "\").", call. = FALSE)
      }
    }
  }

  # Summarise missing-data exposure on the analysis items.
  n_total    <- nrow(data)
  n_complete <- sum(stats::complete.cases(data[, all_items, drop = FALSE]))
  n_missing  <- n_total - n_complete

  if (!is.null(ordered_items)) {
    miss_line <- if (missing == "pairwise") {
      sprintf("  Missing   -> \"pairwise\" (all available pairs per polychoric; %d/%d complete, %d with >=1 NA)",
              n_complete, n_total, n_missing)
    } else {
      sprintf("  Missing   -> \"listwise\" (%d/%d complete cases used; %d dropped)",
              n_complete, n_total, n_missing)
    }
    message(
      "Note: ordered = TRUE detected.\n",
      "  CFA   (R)  -> WLSMV with polychoric correlations\n",
      "  ESEM  (R)  -> efa() block + WLSMV + target rotation\n",
      "  BESEM (R)  -> Bifactor WLSMV (bifactory custom polychoric + orthogonal target rotation)\n",
      "  All Mplus  -> WLSMV with CATEGORICAL block\n",
      "  All models use WLSMV consistently.\n",
      miss_line)
  } else {
    if (missing == "ml") {
      risk_line <- if (length(factors) >= 4L) {
        sprintf(
          "\n  Warning: %d specific factors + FIML puts B-ESEM in the local-optimum risk zone for the FIML + GPArotation::targetT path used here (G absorbs specific-factor variance). Inspect B-ESEM loadings and verify against other software such as Mplus. Use missing = \"listwise\" to avoid.",
          length(factors))
      } else ""
      message(sprintf(
        "Continuous MLR: missing = \"fiml\" applied to CFA, ESEM, and B-ESEM (all %d rows; %d have >=1 NA on analysis items).%s",
        n_total, n_missing, risk_line))
    } else {
      message(sprintf(
        "Continuous MLR: missing = \"listwise\" (%d/%d complete cases used; %d dropped).",
        n_complete, n_total, n_missing))
    }
  }

  # Validate group variable
  if (!is.null(group)) {
    if (!group %in% colnames(data))
      stop("Group variable '", group, "' not found in data.", call. = FALSE)
    group_levels <- sort(unique(na.omit(data[[group]])))
    message("Multi-group model: group = '", group, "' (",
            length(group_levels), " groups: ",
            paste(group_levels, collapse = ", "), ")")
  } else {
    group_levels <- NULL
  }

  structure(
    list(
      factors           = factors,
      factor_names      = factor_names,
      all_items         = all_items,
      nfactors          = length(factors),
      data              = data,
      label             = label,
      cfa_syntax        = cfa_syntax,
      target            = tgt,
      bifactor_target   = btgt,
      ordered           = ordered_items,
      group             = group,
      group_levels      = group_levels,
      estimator_cfa     = estimator_cfa,
      estimator_esem    = estimator_esem,
      missing           = missing
    ),
    class = "esem_spec"
  )
}


#' Print an esem_spec Object
#' @param x An \code{esem_spec} object.
#' @param ... Ignored.
#' @return Invisibly returns \code{x}; called for the side effect of printing
#'   the model specification.
#' @export
print.esem_spec <- function(x, ...) {
  cat("\n==============================================\n")
  cat(" Model Specification:", x$label, "\n")
  cat("==============================================\n\n")
  cat("Factors      :", x$nfactors, "\n")
  cat("Items total  :", length(x$all_items), "\n")
  cat("Estimator (CFA)   :", x$estimator_cfa, "\n")
  cat("Estimator (ESEM)  :", x$estimator_esem, "\n")
  cat("Missing data      :", x$missing, "\n")
  if (!is.null(x$ordered)) {
    cat("Note: ordered data -> WLSMV for all models.\n")
  }
  cat("Ordered items     :", if (is.null(x$ordered)) "none" else
        paste0(length(x$ordered), " items"), "\n")
  cat("Group        :", if (is.null(x$group)) "none (single group)" else
        paste0(x$group, " (", length(x$group_levels), " groups: ",
               paste(x$group_levels, collapse = ", "), ")"), "\n\n")
  for (f in x$factor_names) {
    cat(sprintf("  %-10s (%d items): %s\n",
                f,
                length(x$factors[[f]]),
                paste(x$factors[[f]], collapse = ", ")))
  }
  cat("\nReady to run\n\n")
  invisible(x)
}


#' Run the Full CFA / ESEM / B-ESEM Comparison Pipeline
#'
#' Takes a model specification from \code{\link{specify_model}} and
#' automatically fits CFA, ESEM, and B-ESEM in R, optionally runs all
#' three in Mplus, and returns a comparison table plus all fitted objects.
#'
#' @param spec An \code{esem_spec} object from \code{\link{specify_model}}.
#' @param mplus_folder Character. Path to a folder for Mplus files. If
#'   \code{NULL} (default), Mplus models are skipped.
#' @param mplus_command Character. Path or command used to invoke Mplus.
#'   Default \code{"Mplus"} (assumes it is on PATH). Requires
#'   \pkg{MplusAutomation}.
#' @param run_alignment Logical. Run alignment check? Default \code{TRUE}.
#' @param group_equal Character vector of lavaan equality constraints
#'   (e.g. \code{"loadings"} for metric, \code{c("loadings","intercepts")}
#'   for scalar invariance).
#' @param n_starts Integer. Random rotation starts for the B-ESEM rotation
#'   search. Default \code{30L} (matches Mplus).
#'
#' @return An object of class \code{"esem_comparison_pipeline"} containing:
#' \describe{
#'   \item{\code{spec}}{The original model specification.}
#'   \item{\code{fit_cfa}}{lavaan CFA fit object.}
#'   \item{\code{fit_esem}}{esem_fit object.}
#'   \item{\code{fit_besem}}{besem_fit object.}
#'   \item{\code{alignment}}{alignment_check result (if run).}
#'   \item{\code{comparison_table}}{Data frame of fit indices.}
#'   \item{\code{mplus_results}}{List of Mplus readModels results (if run).}
#' }
#'
#' @examples
#' data("HolzingerSwineford1939", package = "lavaan")
#'
#' spec <- specify_model(
#'   Visual  = c("x1", "x2", "x3"),
#'   Textual = c("x4", "x5", "x6"),
#'   Speed   = c("x7", "x8", "x9"),
#'   data    = HolzingerSwineford1939,
#'   label   = "Holzinger-Swineford"
#' )
#'
#' \donttest{
#' # Fit CFA, ESEM, and B-ESEM in R and build the comparison table
#' results <- run_comparison(spec, n_starts = 5L)
#' print(results)
#'
#' # Access individual fits
#' summary(results$fit_esem, fit.measures = TRUE, standardized = TRUE)
#' std_loadings(results$fit_besem, suppress = 0.10)
#' factor_correlations(results$fit_esem)
#' }
#'
#' \dontrun{
#' # Also run all three models in Mplus and compare side by side.
#' # Requires a licensed Mplus installation reachable via `mplus_command`.
#' results <- run_comparison(spec, mplus_folder = tempfile("mplus_"))
#' print(results)
#' }
#'
#' @export
run_comparison <- function(spec,
                           mplus_folder   = NULL,
                           mplus_command  = "Mplus",
                           run_alignment  = TRUE,
                           group_equal    = NULL,
                           n_starts       = 30L) {

  if (!inherits(spec, "esem_spec"))
    stop("`spec` must be an esem_spec object from specify_model().",
         call. = FALSE)

  message("======================================================")
  message(paste(" Running comparison pipeline:", spec$label))
  message("======================================================\n")

  data           <- spec$data
  estimator_cfa  <- spec$estimator_cfa
  estimator_esem <- spec$estimator_esem
  ordered        <- spec$ordered
  missing <- spec$missing
  # run_comparison() always fits single-group models (CFA vs ESEM vs B-ESEM).
  # Multi-group analysis is handled by esem_invariance().  Ignore spec$group here.
  group          <- NULL
  group_equal    <- NULL

  # -- 1. Alignment check ------------------------------------------------------
  alignment <- NULL
  if (run_alignment) {
    # Fit minimal CFA first for alignment check
    align_args <- list(model = spec$cfa_syntax, data = data,
                       std.lv = TRUE, estimator = estimator_cfa)
    if (!is.null(ordered) && length(ordered) > 0) align_args$ordered <- ordered
    cfa_quick <- tryCatch(
      do.call(lavaan::cfa, align_args),
      error = function(e) NULL
    )
    if (!is.null(cfa_quick)) {
      clusters <- setNames(
        rep(spec$factor_names, lengths(spec$factors)),
        spec$all_items
      )
      alignment <- tryCatch(
        # For alignment check always use Pearson correlations (not polychoric)
        # to avoid dependency on ordered estimator
        alignment_check(data, clusters, cfa_fit = cfa_quick,
                        is_cor = FALSE),
        error = function(e) {
          message("Alignment check failed: ", conditionMessage(e))
          NULL
        }
      )
      # alignment stored in results$alignment; printed separately by the caller
    }
  }


  # -- 2. Fit R models ---------------------------------------------------------
  is_ordered <- !is.null(ordered) && length(ordered) > 0

  # CFA
  message("Fitting 1/3: CFA via lavaan::cfa() ",
      if (is_ordered) "(WLSMV + theta)" else "(MLR)",
      "... ")
  cfa_args <- list(
    model     = spec$cfa_syntax,
    data      = data,
    std.lv    = TRUE,
    estimator = estimator_cfa,
    missing   = missing
  )
  if (is_ordered)             cfa_args$ordered          <- ordered
  if (is_ordered)             cfa_args$parameterization <- "theta"
  if (!is.null(group))        cfa_args$group            <- group
  if (!is.null(group_equal))  cfa_args$group.equal      <- group_equal
  fit_cfa <- do.call(lavaan::cfa, cfa_args)
  message("done")

  # ESEM: lavaan efa() path matches Mplus exactly (same rotation local minimum).
  # The custom DWLS path (method="dwls") is available as a fallback for Heywood

  # cases but finds different rotation minima -- not used by default.
  if (is_ordered) {
    message("Fitting 2/3: ESEM via lavaan efa() block (WLSMV + targetQ rotation)... ")
    fit_esem <- esem_ordered(
      data         = data,
      nfactors     = spec$nfactors,
      indicators   = spec$all_items,
      rotation     = "target",
      target       = spec$target,
      factor_names = spec$factor_names,
      group        = group,
      group_equal  = group_equal,
      missing      = missing
    )
    message("done")
  } else {
    message("Fitting 2/3: ESEM via lavaan efa() block (MLR + targetQ rotation)... ")
    fit_esem <- esem(
      data         = data,
      nfactors     = spec$nfactors,
      indicators   = spec$all_items,
      rotation     = "target",
      target       = spec$target,
      factor_names = spec$factor_names,
      estimator    = estimator_esem,
      missing      = missing,
      group        = group,
      group_equal  = group_equal
    )
    message("done")
  }

  # B-ESEM: post-hoc rotation method for ordered (Asparouhov & Muthen 2009);
  # unrestricted WLSMV fit + bifactor targetT rotation + numDeriv SE correction.
  if (is_ordered) {
    message("Fitting 3/3: B-ESEM via custom DWLS + targetT rotation ",
        "(Asparouhov & Muthen 2009, ", n_starts, " random starts)...")
    fit_besem <- besem_ordered(
      data             = data,
      specific_factors = spec$factors,
      method           = "rotation",
      n_starts         = n_starts,
      group            = group,
      group_equal      = group_equal,
      missing          = missing
    )
  } else {
    message("Fitting 3/3: B-ESEM (MLR + targetT rotation, ",
        n_starts, " random starts)... ")
    fit_besem <- besem(
      data             = data,
      specific_factors = spec$factors,
      estimator        = estimator_esem,
      missing          = missing,
      group            = group,
      group_equal      = group_equal,
      n_starts         = n_starts
    )
  }
  message("done")

  # For ML/MLR: pre-cache standardized solution so parameters() is instant.
  # standardizedsolution() on EFA rotation models with MLR/FIML can be slow;
  # paying the cost once here avoids repeated computation on every parameters() call.
  if (!is_ordered) {
    for (.nm in c("esem", "besem")) {
      .fo <- if (.nm == "esem") fit_esem else fit_besem
      if (is.null(.fo$std_rotated_loadings)) {
        message(sprintf("  Caching standardized solution (%s)... ", toupper(.nm)))
        .lav <- .fo$lavaan_fit
        .cc  <- .cache_std_loadings(.lav)
        if (!is.null(.cc)) {
          if (.nm == "esem")   { fit_esem$std_rotated_loadings  <- .cc$L; fit_esem$se_loadings  <- .cc$SE }
          if (.nm == "besem")  { fit_besem$std_rotated_loadings <- .cc$L; fit_besem$se_loadings <- .cc$SE }
          message("done")
        } else {
          message("skipped (will compute on demand)")
        }
      }
    }
    suppressWarnings(rm(.nm, .fo, .lav, .cc))
  }

  # -- 3. Mplus models (optional) ----------------------------------------------
  mplus_results <- NULL

  if (!is.null(mplus_folder)) {
    message("\nRunning Mplus models...")
    mplus_results <- tryCatch(
      .run_mplus_pipeline(.spec_nogroup(spec), mplus_folder,
                          mplus_command = mplus_command),
      error = function(e) {
        message("Mplus pipeline failed: ", conditionMessage(e))
        message("Continuing with R results only.")
        NULL
      }
    )
  }

  # -- 4. Build comparison table -----------------------------------------------
  message("\nBuilding comparison table...")
  comparison_table <- .build_comparison_table(
    fit_cfa, fit_esem, fit_besem, mplus_results
  )

  # -- 4b. R vs Mplus loading delta (only when Mplus results present) -----------
  mplus_delta <- if (!is.null(mplus_results))
    .mplus_delta(fit_cfa, fit_esem, fit_besem, mplus_results, spec$all_items)
  else
    NULL

  # -- 5. Print summary --------------------------------------------------------
  if (is_ordered) {
    ct        <- comparison_table
    cfi_cfa   <- ct[ct[, 1] == "CFI", "CFA_R",  drop = TRUE]
    cfi_besem <- ct[ct[, 1] == "CFI", "BESEM_R", drop = TRUE]
    if (length(cfi_cfa) && length(cfi_besem) &&
        !is.na(cfi_cfa) && !is.na(cfi_besem) && cfi_besem < cfi_cfa) {
      message(
        "\nNOTE: BESEM_R fit (CFI=", cfi_besem, ") is lower than CFA_R (CFI=", cfi_cfa, ").\n",
        "  This is unexpected -- check for convergence warnings in besem_ordered()."
      )
    }
  }
  message("")
  .print_pipeline_summary(spec, fit_cfa, fit_esem, fit_besem, comparison_table)

  # Attach spec to each S3 fit so standalone parameters() calls can build the
  # primary/cross colour map without the pipeline wrapper.
  if (!is.null(fit_esem))  fit_esem$spec  <- spec
  if (!is.null(fit_besem)) fit_besem$spec <- spec

  structure(
    list(
      spec             = spec,
      fit_cfa          = fit_cfa,
      fit_esem         = fit_esem,
      fit_besem        = fit_besem,
      alignment        = alignment,
      comparison_table = comparison_table,
      mplus_results    = mplus_results,
      mplus_delta      = mplus_delta
    ),
    class = "esem_comparison_pipeline"
  )
}


#' Print Method for esem_comparison_pipeline
#' @param x An \code{esem_comparison_pipeline} object.
#' @param hints Logical. Print a short help block listing the accessors
#'   available on the object (e.g. \code{x$fit_esem}). Default \code{TRUE}.
#' @param ... Ignored.
#' @return Invisibly returns \code{x}; called for the side effect of printing
#'   the comparison-pipeline overview.
#' @export
print.esem_comparison_pipeline <- function(x, hints = TRUE, ...) {
  cat("\n======================================================\n")
  cat(" Pipeline Results:", x$spec$label, "\n")
  cat("======================================================\n\n")
  cat("Fit Index Comparison:\n\n")
  print(x$comparison_table, row.names = FALSE)
  if (isTRUE(attr(x$comparison_table, "wlsmv_note"))) {
    cat(paste0(
      "\n* WLSMV dX2 computed via lavaan scaled difference test;",
      " may not match Mplus DIFFTEST.\n",
      "  Mplus dp omitted for WLSMV (statistic is not X2-distributed).\n",
      "  B-ESEM vs ESEM dX2 suppressed on WLSMV path",
      " (auxiliary lavaan fit not comparable).\n"
    ))
  }

  cat("\nFactor Correlations:\n")
  cat("  CFA:\n")
  cor_cfa_print <- lavaan::lavInspect(x$fit_cfa, "cor.lv")
  if (is.list(cor_cfa_print)) cor_cfa_print <- cor_cfa_print[[1L]]
  cor_cfa_mat <- round(cor_cfa_print, 3)
  print(cor_cfa_mat)
  cat("\n  ESEM:\n")
  cor_esem_mat <- factor_correlations(x$fit_esem)
  print(cor_esem_mat)
  # Delta: change in |correlation| (lower-triangle pairs only)
  fnames <- rownames(cor_cfa_mat)
  if (!is.null(fnames) && identical(fnames, rownames(cor_esem_mat))) {
    pairs <- which(lower.tri(cor_cfa_mat), arr.ind = TRUE)
    if (nrow(pairs) > 0L) {
      delta_parts <- character(nrow(pairs))
      for (k in seq_len(nrow(pairs))) {
        i <- pairs[k, 1L]; j <- pairs[k, 2L]
        d <- abs(cor_esem_mat[i, j]) - abs(cor_cfa_mat[i, j])
        delta_parts[k] <- sprintf("%s-%s = %+.3f", fnames[i], fnames[j], d)
      }
      cat(sprintf("\n  d|r| CFA->ESEM: %s\n",
                  paste(delta_parts, collapse = ", ")))
      cat("  (negative = reduction in |correlation|; expected when cross-loadings absorb shared variance)\n")
    }
  }
  cat("\n  B-ESEM: all fixed to 0 (orthogonal)\n")

  # Heywood notification
  hl <- if (!is.null(x$fit_esem)) x$fit_esem$heywood_log else NULL
  if (!is.null(hl) && length(hl$detected$item) > 0L) {
    hw_df <- hl$detected
    if (isTRUE(hl$resolved)) {
      cat("\n(i) ESEM Heywood case(s) detected and resolved via re-rotation:\n")
      for (i in seq_len(nrow(hw_df))) {
        fin_lam <- NA_real_
        if (!is.null(hl$final_loadings)) {
          fl_row <- hl$final_loadings[hl$final_loadings$item == hw_df$item[i], ]
          if (nrow(fl_row)) fin_lam <- fl_row$lambda[1L]
        }
        cat(sprintf("  %s -> %s  lam = %.3f -> %.3f\n",
                    hw_df$factor[i], hw_df$item[i], hw_df$lambda[i], fin_lam))
      }
      cat(sprintf("  Method: %s rotation (%d round%s). Fit statistics unchanged.\n",
                  hl$final_rotation, hl$rounds,
                  if (hl$rounds != 1L) "s" else ""))
    } else {
      cat("\n(!) ESEM Heywood case(s) -- standardised loading > 1.0 (inadmissible):\n")
      for (i in seq_len(nrow(hw_df)))
        cat(sprintf("  %s -> %s  lam = %.3f\n",
                    hw_df$factor[i], hw_df$item[i], hw_df$lambda[i]))
      cat("  Not resolved. Interpret with caution.\n")
      cat("  Consider: fewer factors, removing item, or B-ESEM specification.\n")
    }
  }

  if (!is.null(x$mplus_delta)) {
    cat("\nR vs Mplus equivalence (mean |D| across all loadings):\n\n")
    print(x$mplus_delta, row.names = FALSE)
    cat("\n")
  }

  if (hints) {
    cat("\nAccess results:\n")
    cat("  results$fit_cfa          -- lavaan CFA object\n")
    cat("  results$fit_esem         -- esem_fit object\n")
    cat("  results$fit_besem        -- besem_fit object\n")
    cat("  results$alignment        -- alignment check\n")
    cat("  results$comparison_table -- fit index table\n")
    cat("  results$mplus_delta      -- R vs Mplus loading delta table\n")
    cat("  plot(results$alignment)  -- alignment ratio plots\n\n")
  }
  invisible(x)
}


# -- Internal pipeline helpers -------------------------------------------------

#' @keywords internal
# Return a copy of spec with group/group_levels set to NULL.
# Used so run_comparison() always generates single-group Mplus .inp files --
# multi-group analysis is handled by esem_invariance() / run_mplus_besem_invariance().
.spec_nogroup <- function(spec) {
  spec$group        <- NULL
  spec$group_levels <- NULL
  spec
}

.mplus_delta <- function(fit_cfa, fit_esem, fit_besem, mplus_results, all_items) {

  # Mplus truncates variable names to 8 chars and uppercases them in output.
  # Build a lookup: truncated-lowercase key -> original R item name.
  mp_key <- substr(tolower(gsub(".", "_", all_items, fixed = TRUE)), 1, 8)
  lookup <- setNames(all_items, mp_key)

  # Extract Mplus BY loadings + SEs from a readModels object
  mp_loads <- function(mf) {
    p        <- mf$parameters$stdyx.standardized
    p        <- p[grepl("\\.BY$", p$paramHeader), ]
    p$factor <- sub("\\.BY$", "", p$paramHeader)
    p$item   <- lookup[substr(tolower(p$param), 1, 8)]
    p[!is.na(p$item), c("factor", "item", "est", "se")]
  }

  rows <- list()

  # -- CFA --------------------------------------------------------------------
  if (!is.null(mplus_results$cfa)) {
    tryCatch({
      # standardizedSolution() gives SEs for STDYX estimates (comparable to Mplus)
      ss    <- lavaan::standardizedSolution(fit_cfa, type = "std.all")
      ss    <- ss[ss$op == "=~", ]
      r     <- data.frame(factor = ss$lhs, item = ss$rhs,
                          r_est  = ss$est.std, r_se = ss$se,
                          stringsAsFactors = FALSE)
      mp    <- mp_loads(mplus_results$cfa)
      m     <- merge(r, mp, by = c("factor","item"))
      dl    <- abs(m$r_est - m$est)
      ds    <- if (all(is.na(m$r_se))) rep(NA_real_, nrow(m)) else abs(m$r_se - m$se)
      rows[["CFA"]] <- data.frame(
        Model       = "CFA",
        N           = nrow(m),
        mean_d_load = round(mean(dl, na.rm = TRUE), 4),
        max_d_load  = round(max(dl,  na.rm = TRUE), 4),
        mean_d_SE   = round(mean(ds, na.rm = TRUE), 4),
        max_d_SE    = round(max(ds,  na.rm = TRUE), 4)
      )
    }, error = function(e) NULL)
  }

  # -- ESEM -------------------------------------------------------------------
  if (!is.null(mplus_results$esem)) {
    tryCatch({
      if (!is.null(fit_esem$std_rotated_loadings)) {
        # DWLS path: loadings and SEs stored directly (auxiliary lavaan_fit is
        # a 1-factor CFA for weight extraction only -- cannot use parameterEstimates)
        L  <- fit_esem$std_rotated_loadings
        SE <- fit_esem$se_loadings
        r  <- data.frame(
          factor = colnames(L)[col(L)],
          item   = rownames(L)[row(L)],
          r_est  = as.vector(L),
          r_se   = if (!is.null(SE)) as.vector(SE) else NA_real_,
          stringsAsFactors = FALSE
        )
      } else {
        # efa() path: standardizedSolution gives STDYX SEs comparable to Mplus
        ss <- lavaan::standardizedSolution(fit_esem$lavaan_fit, type = "std.all")
        ss <- ss[ss$op == "=~", ]
        r  <- data.frame(factor = ss$lhs, item = ss$rhs,
                         r_est  = ss$est.std, r_se = ss$se,
                         stringsAsFactors = FALSE)
      }
      mp  <- mp_loads(mplus_results$esem)
      m   <- merge(r, mp, by = c("factor", "item"))
      dl  <- abs(m$r_est - m$est)
      ds  <- if (all(is.na(m$r_se))) rep(NA_real_, nrow(m)) else abs(m$r_se - m$se)
      rows[["ESEM"]] <- data.frame(
        Model       = "ESEM",
        N           = nrow(m),
        mean_d_load = round(mean(dl, na.rm = TRUE), 4),
        max_d_load  = round(max(dl,  na.rm = TRUE), 4),
        mean_d_SE   = round(mean(ds, na.rm = TRUE), 4),
        max_d_SE    = round(max(ds,  na.rm = TRUE), 4)
      )
    }, error = function(e) NULL)
  }

  # -- BESEM ------------------------------------------------------------------
  if (!is.null(mplus_results$besem) &&
      !is.null(fit_besem$std_rotated_loadings)) {
    tryCatch({
      L   <- fit_besem$std_rotated_loadings
      SE  <- fit_besem$se_loadings
      r   <- data.frame(
        factor = colnames(L)[col(L)],
        item   = rownames(L)[row(L)],
        r_est  = as.vector(L),
        r_se   = as.vector(SE),
        stringsAsFactors = FALSE
      )
      mp  <- mp_loads(mplus_results$besem)
      m   <- merge(r, mp, by = c("factor", "item"))
      dl  <- abs(m$r_est - m$est)
      ds  <- abs(m$r_se  - m$se)
      rows[["BESEM"]] <- data.frame(
        Model       = "BESEM",
        N           = nrow(m),
        mean_d_load = round(mean(dl, na.rm = TRUE), 4),
        max_d_load  = round(max(dl,  na.rm = TRUE), 4),
        mean_d_SE   = round(mean(ds, na.rm = TRUE), 4),
        max_d_SE    = round(max(ds,  na.rm = TRUE), 4)
      )
    }, error = function(e) NULL)
  }

  if (length(rows) == 0) return(NULL)
  do.call(rbind, rows)
}


.run_mplus_pipeline <- function(spec, folder, mplus_command = "Mplus") {

  if (!requireNamespace("MplusAutomation", quietly = TRUE))
    stop("MplusAutomation is required when 'mplus_folder' is supplied. ",
         "Install with: install.packages('MplusAutomation')", call. = FALSE)

  dir.create(folder, showWarnings = FALSE, recursive = TRUE)

  # Prepare data
  all_vars_for_data <- spec$all_items
  if (!is.null(spec$group)) all_vars_for_data <- c(all_vars_for_data, spec$group)
  mplus_df <- spec$data[, all_vars_for_data, drop = FALSE]
  mplus_df[is.na(mplus_df)] <- 999
  data_file <- file.path(folder, "mydata.dat")
  MplusAutomation::prepareMplusData(mplus_df, filename = data_file)

  items   <- spec$all_items
  fnames  <- spec$factor_names
  factors <- spec$factors
  ordered <- spec$ordered
  group   <- spec$group
  estimator_cfa  <- spec$estimator_cfa
  estimator_esem <- spec$estimator_esem

  # Sanitize variable names for Mplus: replace dots with underscores.
  # Mplus rejects dots in variable names. The data file is headerless so
  # the NAMES line is just positional labels -- they don't need to match R names.
  .mp <- function(x) gsub(".", "_", x, fixed = TRUE)
  mp_items   <- .mp(items)
  mp_factors <- lapply(factors, .mp)
  mp_ordered <- .mp(ordered)
  mp_group   <- if (!is.null(group)) .mp(group) else character(0)
  mp_vars    <- c(mp_items, mp_group)

  # Mplus CATEGORICAL block (if ordered items)
  # Mplus GROUPING block (if multi-group)
  group_block <- if (!is.null(group)) {
    levels <- sort(unique(na.omit(spec$data[[group]])))
    level_str <- paste(vapply(levels, function(l)
      paste0(l, "=", l), character(1)), collapse = " ")
    paste0("  GROUPING = ", mp_group, " (", level_str, ");\n")
  } else ""

  # Mplus estimator line (unused - kept for reference)
  estimator_line <- paste0("  ESTIMATOR = ", estimator_cfa, ";")

  # Mplus group_equal equivalent (for metric/scalar invariance)
  model_constraint_block <- ""

  # Build NAMES and USEVARIABLES lines -- wrapped at 70 chars, semicolon on same last line
  .inp_varlist <- function(prefix, vars, width = 70) {
    indent  <- paste(rep(" ", nchar(prefix)), collapse = "")
    lines   <- character(0)
    current <- prefix
    for (v in vars) {
      candidate <- paste0(current, v, " ")
      if (nchar(candidate) > width && current != prefix) {
        lines   <- c(lines, current)
        current <- paste0(indent, v, " ")
      } else {
        current <- candidate
      }
    }
    # Last line gets the semicolon -- trim trailing space first
    current <- trimws(current, "right")
    lines   <- c(lines, paste0(current, ";"))
    paste(lines, collapse = "\n")
  }

  cat_block <- if (!is.null(ordered) && length(ordered) > 0) {
    # Wrap CATEGORICAL line at 80 chars (Mplus limit is 90)
    cat_line  <- .inp_varlist("  CATEGORICAL = ", mp_ordered, width = 80)
    paste0(cat_line, "\n")  # .inp_varlist already appends semicolon
  } else ""



  names_line   <- .inp_varlist("  NAMES = ",        mp_vars)
  usevars_line <- .inp_varlist("  USEVARIABLES = ", mp_vars)

  # -- CFA input
  # Add asterisk (*) to first item per factor to free it from Mplus's implicit
  # referent-loading=1 default. Combined with F@1, this matches lavaan std.lv=TRUE:
  # factor variance = 1, all loadings freely estimated.
  cfa_by <- paste(vapply(fnames, function(f) {
    items_f      <- mp_factors[[f]]
    items_starred <- c(paste0(items_f[1], "*"), items_f[-1])
    .inp_varlist(paste0("  ", f, " BY "), items_starred, width = 85)
  }, character(1)), collapse = "\n")

  # Factor variance = 1 constraints (matching R's std.lv=TRUE)
  cfa_var <- paste(vapply(fnames, function(f) paste0("  ", f, "@1;"),
                          character(1)), collapse = "\n")

  # Listwise for continuous (all models use same complete-case data);
  # pairwise for ordered (WLSMV convention).  If the spec requested FIML,
  # drop LISTWISE = ON so Mplus's default ML-family FIML handling applies
  # (parallels missing = "fiml"/"ml"/"direct" on the R side).
  is_fiml <- !is.null(spec$missing) &&
              tolower(spec$missing) %in% c("ml", "ml.x", "fiml",
                                            "direct", "direct.ml")
  mplus_data_block <- if (is.null(ordered) || length(ordered) == 0) {
    if (is_fiml)
      "DATA:\n  FILE = mydata.dat;\n\n"
    else
      "DATA:\n  FILE = mydata.dat;\n  LISTWISE = ON;\n\n"
  } else {
    "DATA:\n  FILE = mydata.dat;\n\n"
  }

  cfa_inp <- paste0(
    "TITLE: CFA - ", spec$label, ";\n\n",
    mplus_data_block,
    "VARIABLE:\n", names_line, "\n",
    usevars_line, "\n",
    cat_block,
    group_block,
    "  MISSING ARE ALL (999);\n\n",
    if (!is.null(ordered) && length(ordered) > 0)
      paste0("ANALYSIS:\n  ESTIMATOR = ", estimator_cfa, ";\n  PARAMETERIZATION = THETA;\n\n")
    else
      paste0("ANALYSIS:\n  ESTIMATOR = ", estimator_cfa, ";\n\n"),
    "MODEL:\n", cfa_by, "\n", cfa_var, "\n\n",
    "OUTPUT:\n  STDYX;\n  SAMPSTAT;\n"
  )
  writeLines(cfa_inp, file.path(folder, "cfa_model.inp"))

  # -- ESEM input
  # For both ordered and continuous: true ESEM with rotation (*1).
  # Mplus WLSMV supports (*1) rotation syntax; ROTATION = TARGET added to ANALYSIS.
  # Lines are wrapped at 85 chars to stay under Mplus 90-char limit.
  esem_by <- paste(vapply(fnames, function(f) {
    primary <- mp_factors[[f]]
    cross   <- setdiff(mp_items, primary)
    primary_str <- .inp_varlist(paste0("  ", f, " BY "), primary, width = 85)
    primary_str <- sub(";$", "", primary_str)
    cross_str   <- .inp_varlist("    ", paste0(cross, "~0"), width = 85)
    cross_str   <- sub(";$", " (*1);", cross_str)
    paste0(primary_str, "\n", cross_str)
  }, character(1)), collapse = "\n\n")
  esem_var_block <- ""  # rotation handles identification for both ordered and continuous

  esem_inp <- paste0(
    "TITLE: ESEM - ", spec$label, ";\n\n",
    mplus_data_block,
    "VARIABLE:\n", names_line, "\n",
    usevars_line, "\n",
    cat_block,
    group_block,
    "  MISSING ARE ALL (999);\n\n",
    if (!is.null(ordered) && length(ordered) > 0)
      paste0("ANALYSIS:\n  ESTIMATOR = ", estimator_esem, ";\n  ROTATION = TARGET;\n  PARAMETERIZATION = THETA;\n\n")
    else
      paste0("ANALYSIS:\n  ESTIMATOR = ", estimator_esem, ";\n  ROTATION = TARGET;\n\n"),
    "MODEL:\n",
    esem_by,
    if (nchar(esem_var_block) > 0) paste0("\n", esem_var_block) else "",
    "\n\n",
    "OUTPUT:\n  STDYX;\n  SAMPSTAT;\n"
  )
  writeLines(esem_inp, file.path(folder, "esem_measurement.inp"))

  # -- B-ESEM input
  # G factor BY line -- all items; (*1) + rotation handles identification
  g_items_line <- .inp_varlist("  G BY ", mp_items, width = 80)
  g_items_line <- sub(";$", "", g_items_line)
  g_by <- paste0(g_items_line, " (*1);")  # rotation for both ordered and continuous

  # B-ESEM: true bifactor ESEM with rotation for both ordered and continuous.
  # (*1) with ROTATION = TARGET (orthogonal) in ANALYSIS handles identification
  # and orthogonality between G and specific factors.
  # Lines wrapped at 85 chars (Mplus 90-char limit).
  besem_by <- paste(vapply(fnames, function(f) {
    primary <- mp_factors[[f]]
    cross   <- setdiff(mp_items, primary)
    primary_str <- .inp_varlist(paste0("  ", f, " BY "), primary, width = 85)
    primary_str <- sub(";$", "", primary_str)
    cross_str   <- .inp_varlist("    ", paste0(cross, "~0"), width = 85)
    cross_str   <- sub(";$", " (*1);", cross_str)
    paste0(primary_str, "\n", cross_str)
  }, character(1)), collapse = "\n\n")

  besem_inp <- paste0(
    "TITLE: B-ESEM - ", spec$label, ";\n\n",
    mplus_data_block,
    "VARIABLE:\n", names_line, "\n",
    usevars_line, "\n",
    cat_block,
    group_block,
    "  MISSING ARE ALL (999);\n\n",
    if (!is.null(ordered) && length(ordered) > 0)
      "ANALYSIS:\n  ESTIMATOR = WLSMV;\n  ROTATION = TARGET (orthogonal);\n  PARAMETERIZATION = THETA;\n\n"
    else
      paste0("ANALYSIS:\n  ESTIMATOR = ", estimator_esem, ";\n  ROTATION = TARGET (orthogonal);\n\n"),
    "MODEL:\n", g_by, "\n\n", besem_by, "\n\n",
    "OUTPUT:\n  STDYX;\n  SAMPSTAT;\n"
  )
  writeLines(besem_inp, file.path(folder, "besem_measurement.inp"))

  # Delete stale .out files so Mplus is forced to re-run
  for (out_f in c("cfa_model.out","esem_measurement.out","besem_measurement.out")) {
    f_path <- file.path(folder, out_f)
    if (file.exists(f_path)) file.remove(f_path)
  }

  # Run all three models in Mplus
  for (inp in c("cfa_model.inp", "esem_measurement.inp", "besem_measurement.inp")) {
    message(paste("  Running Mplus:", sub("\\.inp$","", inp), "... "))
    MplusAutomation::runModels(file.path(folder, inp),
                               Mplus_command = mplus_command)
    out_file <- file.path(folder, sub("\\.inp$", ".out", inp))
    if (file.exists(out_file)) {
      out_lines <- readLines(out_file)
      status    <- if (any(grepl("\\*\\*\\* ERROR", out_lines))) "ERROR" else "OK"
    } else {
      status <- "NO OUTPUT"
    }
    message(status)
  }

  # Read results
  list(
    cfa   = MplusAutomation::readModels(file.path(folder, "cfa_model.out")),
    esem  = MplusAutomation::readModels(file.path(folder, "esem_measurement.out")),
    besem = MplusAutomation::readModels(file.path(folder, "besem_measurement.out"))
  )
}


.mplus_only <- function(mplus_folder) {
  if (!requireNamespace("MplusAutomation", quietly = TRUE))
    stop("MplusAutomation is required for run_mode 2.", call. = FALSE)
  read_out <- function(name) {
    tryCatch(MplusAutomation::readModels(file.path(mplus_folder, name)),
             error = function(e) NULL)
  }
  mp <- list(
    cfa   = read_out("cfa_model.out"),
    esem  = read_out("esem_measurement.out"),
    besem = read_out("besem_measurement.out")
  )
  fi <- function(m, label) {
    if (is.null(m)) return(data.frame(Model=label, CFI=NA, TLI=NA, RMSEA=NA, SRMR=NA))
    s <- m$summaries
    data.frame(Model = label,
               CFI   = round(s$CFI, 3),
               TLI   = round(s$TLI, 3),
               RMSEA = round(s$RMSEA_Estimate, 3),
               SRMR  = round(s$SRMR, 3))
  }
  tbl <- rbind(fi(mp$cfa,"CFA"), fi(mp$esem,"ESEM"), fi(mp$besem,"BESEM"))
  message("\nMplus fit indices:")
  message(paste(utils::capture.output(print(tbl, row.names = FALSE)), collapse = "\n"))
  structure(list(mplus_results = mp, comparison_table = tbl,
                 fit_cfa = NULL, fit_esem = NULL, fit_besem = NULL,
                 spec = NULL, alignment = NULL),
            class = "esem_comparison_pipeline")
}


.rmsea_ci <- function(chisq, df, n) {
  # RMSEA 90% CI via non-central chi-squared inversion (MacCallum et al. 1996).
  # As NCP increases, P(X >= chisq | ncp) increases (distribution shifts right).
  # Lower bound: smallest NCP where P = 0.05  -> small NCP -> small RMSEA.
  # Upper bound: largest  NCP where P = 0.95  -> large NCP -> large RMSEA.
  if (anyNA(c(chisq, df, n)) || df <= 0 || n <= 1) return(c(NA_real_, NA_real_))
  lo_ncp <- tryCatch(suppressWarnings({
    if (pchisq(chisq, df, ncp = 0, lower.tail = FALSE) >= 0.05) 0
    else uniroot(function(l) pchisq(chisq, df, ncp = l, lower.tail = FALSE) - 0.05,
                 lower = 0, upper = max(chisq * 4, 100))$root
  }), error = function(e) NA_real_)
  hi_ncp <- tryCatch(suppressWarnings(
    uniroot(function(l) pchisq(chisq, df, ncp = l, lower.tail = FALSE) - 0.95,
            lower = 0, upper = max(chisq * 10, 500))$root
  ), error = function(e) NA_real_)
  c(sqrt(lo_ncp / (df * (n - 1L))),   # L90%CI
    sqrt(hi_ncp / (df * (n - 1L))))   # U90%CI
}


#' Mplus-Matched Fit Indices for a Single Model
#'
#' Returns CFI, TLI, RMSEA (with 90\% CI) and SRMR formatted as a named
#' character vector. Uses the \code{.scaled} variants when the fit was
#' estimated with MLR/WLSMV, matching Mplus's default output. For B-ESEM /
#' ESEM WLSMV fits that carry pre-computed Mplus-matched values in
#' \code{$wlsmv_stats}, those values are used directly and the RMSEA CI is
#' computed from the scaled chi-square via non-central chi-squared inversion
#' (MacCallum, Browne & Sugawara, 1996).
#'
#' @param fit An \code{esem_fit} / \code{besem_fit} object, or a raw
#'   \code{lavaan} S4 fit.
#'
#' @return A named character vector with elements \code{"CFI"}, \code{"TLI"},
#'   \code{"RMSEA [90\% CIs]"}, and \code{"SRMR"}. Stack rows with
#'   \code{rbind()} (wrapped in \code{noquote()}) to build a comparison
#'   table.
#'
#' @examples
#' data("HolzingerSwineford1939", package = "lavaan")
#' d <- HolzingerSwineford1939[, paste0("x", 1:9)]
#'
#' \donttest{
#' fit_e <- esem(d, nfactors = 3)
#' fit_b <- besem(d, specific_factors = list(
#'   Visual  = c("x1", "x2", "x3"),
#'   Textual = c("x4", "x5", "x6"),
#'   Speed   = c("x7", "x8", "x9")
#' ), n_starts = 5L)
#'
#' noquote(rbind(
#'   ESEM  = fit_indices(fit_e),
#'   BESEM = fit_indices(fit_b)
#' ))
#' }
#' @export
fit_indices <- function(fit) {
  lf <- if (inherits(fit, "esem_fit"))      fit$lavaan_fit
        else if (isS4(fit))                  fit
        else stop("`fit` must be an esem_fit / besem_fit or a lavaan S4 object.",
                  call. = FALSE)
  ws <- if (inherits(fit, "esem_fit")) fit$wlsmv_stats else NULL
  if (!is.null(ws)) {
    n  <- lavaan::lavInspect(lf, "ntotal")
    ci <- .rmsea_ci(ws$chisq, ws$df, n)
    out <- c(sprintf("%.3f", ws$cfi),
             sprintf("%.3f", ws$tli),
             sprintf("%.3f [%.3f, %.3f]", ws$rmsea, ci[1], ci[2]),
             sprintf("%.3f", ws$srmr))
  } else {
    # Route through the shared helper so print() and fit_indices() agree
    # on scaled indices and the Mplus SRMR correction.
    fi <- .fit_indices_lav(lf)
    if (is.null(fi)) {
      out <- c("NA", "NA", "NA [NA, NA]", "NA")
    } else {
      out <- c(sprintf("%.3f", fi[["cfi"]]),
               sprintf("%.3f", fi[["tli"]]),
               sprintf("%.3f [%.3f, %.3f]",
                       fi[["rmsea"]], fi[["rmsea.ci.lower"]],
                       fi[["rmsea.ci.upper"]]),
               sprintf("%.3f", fi[["srmr"]]))
    }
  }
  setNames(out, c("CFI", "TLI", "RMSEA [90% CIs]", "SRMR"))
}


.build_comparison_table <- function(fit_cfa, fit_esem, fit_besem,
                                     mplus_results = NULL) {
  # For WLSMV models (ordered data), lavaan reports two sets of fit indices:
  #   "cfi"        = naive DWLS chi-square (NOT comparable to Mplus)
  #   "cfi.scaled" = mean-variance adjusted WLSMV chi-square (matches Mplus)
  # We detect which to use by checking whether the scaled variant is available.
  #
  # SRMR: corrected to match Mplus convention (see comment inside r_fi()).
  #
  # CFI/TLI residual gap vs Mplus: historically ~0.001-0.004 for large item
  # sets; root cause was an upward bias in lavaan's polychoric cell-probability
  # floor (`lav_bvord_noexo_pi()`), fixed upstream in lavaan 0.6-22.2568
  # (April 2026).  On the fixed version, CFA / ESEM / B-ESEM fit indices
  # agree with Mplus to 3 dp across the validation test set; any remaining
  # delta is at numerical-precision level.
  r_fi <- function(fit) {
    # B-ESEM ordered (DWLS from-scratch): use custom WLSMV stats directly
    if (inherits(fit, "esem_fit") && !is.null(fit$wlsmv_stats)) {
      ws <- fit$wlsmv_stats
      n  <- tryCatch(nrow(fit$spec$data), error = function(e) NULL)
      if (!isTRUE(n > 0L))
        n <- tryCatch(as.integer(lavaan::lavInspect(fit$lavaan_fit, "nobs")),
                     error = function(e) NA_integer_)
      ci <- .rmsea_ci(ws$chisq, ws$df, n)
      return(round(c(ws$cfi, ws$tli, ws$rmsea, ci[1L], ci[2L],
                     ws$srmr, ws$chisq, ws$df, ws$pvalue,
                     NA_real_, NA_real_, NA_real_), 3))
    }
    obj    <- if (inherits(fit, "esem_fit")) lavaan_fit(fit) else fit
    all_fm <- tryCatch(lavaan::fitMeasures(obj), error = function(e) NULL)
    if (is.null(all_fm)) return(rep(NA_real_, 12L))
    if ("cfi.scaled" %in% names(all_fm) && !is.na(all_fm["cfi.scaled"])) {
      # WLSMV: use scaled indices. Read from all_fm -- no second fitMeasures() call.
      vals <- as.numeric(all_fm[c("cfi.scaled", "tli.scaled", "rmsea.scaled", "srmr")])
      # Mplus WLSMV SRMR divides by n_pairs + n_thresholds (threshold residuals
      # are 0 but inflate the denominator).  lavaan uses only n_pairs.
      # Correct: SRMR_mplus = SRMR_lavaan * sqrt(n_pairs / n_wls).
      wls_obs <- tryCatch(lavaan::lavInspect(obj, "wls.obs"), error = function(e) NULL)
      if (!is.null(wls_obs)) {
        if (is.list(wls_obs)) wls_obs <- wls_obs[[1L]]
        n_wls  <- length(wls_obs)
        cor_ov <- lavaan::lavInspect(obj, "cor.ov")
        if (is.list(cor_ov)) cor_ov <- cor_ov[[1L]]
        n_items <- nrow(cor_ov)
        n_pairs <- n_items * (n_items - 1L) / 2L
        if (n_wls > n_pairs) vals[4L] <- vals[4L] * sqrt(n_pairs / n_wls)
      }
      ci_lo    <- as.numeric(all_fm["rmsea.scaled.ci.lower"])
      ci_hi    <- as.numeric(all_fm["rmsea.scaled.ci.upper"])
      chi_vals <- as.numeric(all_fm[c("chisq.scaled", "df.scaled", "pvalue.scaled")])
      # AIC/BIC: available for MLR but not WLSMV
      ic_vals <- if (!is.na(all_fm["aic"]))
        as.numeric(all_fm[c("aic", "bic", "bic2")])
      else
        c(NA_real_, NA_real_, NA_real_)
    } else {
      vals     <- as.numeric(all_fm[c("cfi", "tli", "rmsea", "srmr")])
      ci_lo    <- as.numeric(all_fm["rmsea.ci.lower"])
      ci_hi    <- as.numeric(all_fm["rmsea.ci.upper"])
      chi_vals <- as.numeric(all_fm[c("chisq", "df", "pvalue")])
      ic_vals  <- as.numeric(all_fm[c("aic", "bic", "bic2")])
    }
    # Fallback: if lavaan didn't populate CI keys, compute from scratch
    if (is.na(ci_lo) || is.na(ci_hi)) {
      n_obs <- tryCatch(sum(lavaan::lavInspect(obj, "nobs")), error = function(e) NA_integer_)
      ci    <- .rmsea_ci(chi_vals[1L], chi_vals[2L], n_obs)
      ci_lo <- ci[1L]; ci_hi <- ci[2L]
    }
    round(c(vals[1L], vals[2L], vals[3L], ci_lo, ci_hi, vals[4L], chi_vals, ic_vals), 3)
  }

  tbl <- data.frame(
    Index   = c("CFI", "TLI", "RMSEA", "RMSEA [L90%CI]", "RMSEA [U90%CI]",
                "SRMR", "X2", "df", "p", "AIC", "BIC", "SABIC"),
    CFA_R   = r_fi(fit_cfa),
    ESEM_R  = r_fi(fit_esem),
    BESEM_R = r_fi(fit_besem),
    stringsAsFactors = FALSE
  )

  mp_fi <- function(res) {
    if (is.null(res)) return(rep(NA_real_, 12L))
    tryCatch({
      s <- res$summaries
      if (is.null(s) || nrow(s) == 0L) return(rep(NA_real_, 12L))
      round(c(
        .mp_val(s, "CFI"),
        .mp_val(s, "TLI"),
        .mp_val(s, "RMSEA_Estimate"),
        .mp_val(s, "RMSEA_90CI_LB",    "RMSEA_90CI_Lower"),
        .mp_val(s, "RMSEA_90CI_UB",    "RMSEA_90CI_Upper"),
        .mp_val(s, "SRMR"),
        .mp_val(s, "ChiSqM_Value",   "Chi_Square_Value"),
        .mp_val(s, "ChiSqM_DF",      "Chi_Square_DF_Value"),
        .mp_val(s, "ChiSqM_PValue",  "Chi_Square_P_Value"),
        .mp_val(s, "AIC"),
        .mp_val(s, "BIC"),
        .mp_val(s, "aBIC")
      ), 3L)
    }, error = function(e) rep(NA_real_, 12L))
  }

  # -- Mplus columns ---------------------------------------------------------
  if (!is.null(mplus_results)) {
    mp_cfa   <- mp_fi(mplus_results$cfa)
    mp_esem  <- mp_fi(mplus_results$esem)
    mp_besem <- mp_fi(mplus_results$besem)
    tbl$CFA_Mplus   <- mp_cfa
    tbl$ESEM_Mplus  <- mp_esem
    tbl$BESEM_Mplus <- mp_besem
    if (!all(is.na(tbl$CFA_Mplus))) {
      tbl <- tbl[, c("Index", "CFA_R", "CFA_Mplus",
                     "ESEM_R", "ESEM_Mplus",
                     "BESEM_R", "BESEM_Mplus")]
    } else {
      message("Mplus results unavailable -- showing R models only.")
      tbl$CFA_Mplus <- tbl$ESEM_Mplus <- tbl$BESEM_Mplus <- NULL
    }
  }

  tbl
}


.print_comparison_table <- function(tbl) {
  # Collapses RMSEA + L90%CI + U90%CI rows into one formatted string for display.
  # The underlying comparison_table keeps all rows numeric for xlsx/CSV export.
  val_cols <- setdiff(names(tbl), "Index")

  # Build character display frame
  disp <- as.data.frame(
    lapply(tbl, as.character), stringsAsFactors = FALSE)
  for (col in val_cols) {
    vals <- suppressWarnings(as.numeric(tbl[[col]]))
    disp[[col]] <- ifelse(is.na(vals), "-", sprintf("%.3f", vals))
  }

  ri <- which(tbl$Index == "RMSEA")
  li <- which(tbl$Index == "RMSEA [L90%CI]")
  hi <- which(tbl$Index == "RMSEA [U90%CI]")

  if (length(ri) == 1L && length(li) == 1L && length(hi) == 1L) {
    disp$Index[ri] <- "RMSEA [90% CI]"
    for (col in val_cols) {
      rv <- as.numeric(tbl[[col]][ri])
      lv <- as.numeric(tbl[[col]][li])
      hv <- as.numeric(tbl[[col]][hi])
      if (!is.na(rv) && !is.na(lv) && !is.na(hv)) {
        disp[[col]][ri] <- sprintf("%.3f [%s, %s]", rv,
          sub("^0", "", sprintf("%.3f", lv)),
          sub("^0", "", sprintf("%.3f", hv)))
      } else if (!is.na(rv)) {
        disp[[col]][ri] <- sprintf("%.3f", rv)
      } else {
        disp[[col]][ri] <- "-"
      }
    }
    disp <- disp[-c(li, hi), , drop = FALSE]
  }

  for (col in val_cols) disp[[col]] <- format(disp[[col]], justify = "right")
  message(paste(utils::capture.output(print(disp, row.names = FALSE)), collapse = "\n"))
}


.mp_val <- function(s, ...) {
  # Extract a single numeric value from an MplusAutomation summaries data frame.
  # Tries each supplied column name in order; returns NA_real_ if none found.
  # Needed because res$summaries[missing_col] returns NULL silently (not NA).
  for (nm in c(...)) {
    v <- s[[nm]]
    if (!is.null(v) && length(v) == 1L && !is.na(v)) return(as.numeric(v))
  }
  NA_real_
}


.print_pipeline_summary <- function(spec, fit_cfa, fit_esem,
                                     fit_besem, comparison_table) {
  message("======================================================")
  message(paste(" Results:", spec$label))
  message("======================================================\n")
  message("Fit Index Comparison:\n")
  message(paste(utils::capture.output(print(comparison_table, row.names = FALSE)), collapse = "\n"))

  message("\nFactor Correlations:")
  cor_cfa <- lavaan::lavInspect(fit_cfa, "cor.lv")
  if (is.list(cor_cfa)) cor_cfa <- cor_cfa[[1L]]   # multi-group: use group 1
  cor_cfa <- round(cor_cfa, 3)
  pairs   <- utils::combn(rownames(cor_cfa), 2, simplify = FALSE)
  message("  CFA: ", paste(vapply(pairs, function(p)
    paste0(p[1],"-",p[2]," = ", cor_cfa[p[1],p[2]]),
    character(1)), collapse = ", "))

  cor_esem <- factor_correlations(fit_esem)
  pairs_e  <- utils::combn(rownames(cor_esem), 2, simplify = FALSE)
  message("  ESEM: ", paste(vapply(pairs_e, function(p)
    paste0(p[1],"-",p[2]," = ", cor_esem[p[1],p[2]]),
    character(1)), collapse = ", "))

  message("  B-ESEM: all 0 (orthogonal)\n")
}
