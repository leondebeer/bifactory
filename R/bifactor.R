#' Bifactor Exploratory Structural Equation Modeling (B-ESEM)
#'
#' Fits a Bifactor ESEM model -- a general factor (G) loading on all
#' indicators plus domain-specific factors each targeting a subset of
#' indicators. All factors are orthogonal (uncorrelated), matching the
#' Mplus B-ESEM specification on **continuous** data via lavaan's \code{efa()}
#' block. With \code{ordered} indicators, the call routes to
#' \code{\link{besem_ordered}} (default custom DWLS/WLSMV path; see
#' \code{?besem_ordered}).
#'
#' @inheritSection doc_estimator_paths Estimator paths
#' @inheritSection doc_estimator_paths Missing-data defaults
#' @inheritSection doc_estimator_paths The lavaan_fit slot on custom WLSMV fits
#'
#' @section Mplus target syntax:
#' \preformatted{
#' Mplus continuous / ordered target syntax:
#'
#' ROTATION = TARGET (orthogonal);
#' MODEL:
#'   G BY batEX1-batCI5 (*1);
#'   EX BY batEX1~1 ... batMD1~0 ... (*1);
#'   MD BY batEX1~0 ... batMD1~1 ... (*1);
#'   CI BY batEX1~0 ... batCI1~1 ... (*1);
#' }
#'
#' @param data A \code{data.frame} of observed indicators.
#' @param specific_factors A named list mapping specific factor names to
#'   their primary indicator names. The general factor G is added
#'   automatically.
#'   Example: \code{list(EX = c("y1","y2"), MD = c("y3","y4"), CI = c("y5","y6"))}.
#' @param indicators Optional character vector of all indicator names. If
#'   \code{NULL} (default), all items from \code{specific_factors} are used.
#' @param g_name Character. Name for the general factor. Default \code{"G"}.
#' @param estimator Character. Default \code{"MLR"}
#'   (ML with Huber-White robust SEs and Satorra-Bentler scaled chi-square).
#' @param std.lv Logical. Fix factor variances to 1. Default \code{TRUE}.
#' @param ordered Character vector of ordinal item names.
#'   When non-\code{NULL}, routes to \code{\link{besem_ordered}}. For Mplus-aligned
#'   ordered B-ESEM use \code{method = "rotation"} there (default in
#'   \code{\link{run_comparison}}).
#' @param group Character. Grouping variable for multi-group B-ESEM.
#' @param group_equal Character vector of lavaan equality constraints.
#' @param missing Character. Missing data handling for **continuous** B-ESEM.
#'   Default \code{"listwise"}. When \code{ordered} is set, passed to
#'   \code{\link{besem_ordered}} (default \code{"listwise"} there; the pipeline
#'   passes \code{"pairwise"} from \code{\link{specify_model}} when used via
#'   \code{\link{run_comparison}}).
#' @param n_starts Integer. Number of random orthogonal starting matrices for
#'   the target rotation (forwarded as \code{rstarts} to lavaan's rotation
#'   engine). Multiple random starts help escape local optima of the rotation
#'   criterion, particularly for bifactor models with many specific factors.
#'   Default \code{30L} (matches Mplus).
#' @param ... Additional arguments passed to \code{lavaan::cfa()}.
#'
#' @return An object of class \code{c("besem_fit", "esem_fit")} with the
#'   same structure as \code{\link{esem}}, plus:
#' \describe{
#'   \item{\code{g_name}}{Name of the general factor.}
#'   \item{\code{specific_factors}}{Named list of specific factor assignments.}
#' }
#'
#' @details
#' ## What makes B-ESEM different from ESEM
#'
#' | | ESEM | B-ESEM |
#' |---|---|---|
#' | Factor structure | k oblique specific factors | 1 general + k orthogonal specific |
#' | Factor correlations | Freely estimated | All fixed to zero (orthogonal) |
#' | Cross-loadings | Estimated via rotation | Estimated via orthogonal target rotation |
#' | G factor | None | Loads freely on all items |
#' | Rotation | Oblique target/geomin | Orthogonal target |
#'
#' ## Orthogonality
#' B-ESEM uses orthogonal target rotation (\code{"targetT"} in lavaan),
#' which constrains all factors to be uncorrelated. This means:
#' \itemize{
#'   \item G is uncorrelated with EX, MD, CI
#'   \item EX, MD, CI are uncorrelated with each other
#' }
#' This matches Mplus \code{ROTATION = TARGET (orthogonal)}.
#'
#' ## Interpreting results
#' \itemize{
#'   \item **G loadings**: variance shared across all items regardless of domain
#'   \item **Specific loadings**: domain-specific variance after accounting for G
#'   \item **omega_h (omega hierarchical)**: reliability of G (use \code{psych::omega()})
#'   \item **omega_s (omega specific)**: reliability of each specific factor
#' }
#'
#' ## Target matrix structure
#' \preformatted{
#'        G   EX   MD   CI
#' y1     1    1    0    0   <- G free on all; EX primary; MD/CI targeted to 0
#' y2     1    1    0    0
#' y6     1    0    1    0   <- MD item
#' y9     1    0    0    1   <- CI item
#' }
#'
#' @seealso \code{\link{esem}} for standard oblique ESEM,
#'   \code{\link{besem_ordered}} for ordinal B-ESEM,
#'   \code{\link{make_bifactor_target}} for the target matrix,
#'   \code{\link{generate_mplus_besem_syntax}} for Mplus comparison.
#'
#' @examples
#' data("HolzingerSwineford1939", package = "lavaan")
#' d <- HolzingerSwineford1939[, paste0("x", 1:9)]
#'
#' \donttest{
#' fit_b <- besem(
#'   data = d,
#'   specific_factors = list(
#'     Visual  = c("x1", "x2", "x3"),
#'     Textual = c("x4", "x5", "x6"),
#'     Speed   = c("x7", "x8", "x9")
#'   ),
#'   n_starts = 5L
#' )
#'
#' summary(fit_b, fit.measures = TRUE, standardized = TRUE)
#' std_loadings(fit_b)        # rows = items, cols = G + specific factors
#' factor_correlations(fit_b) # should all be ~0 (orthogonal)
#'
#' # Compare B-ESEM vs standard ESEM
#' fit_esem <- esem(d, nfactors = 3)
#' lavaan::fitMeasures(lavaan_fit(fit_b),    c("cfi", "rmsea", "aic"))
#' lavaan::fitMeasures(lavaan_fit(fit_esem), c("cfi", "rmsea", "aic"))
#' }
#'
#' @importFrom lavaan cfa
#' @export
besem <- function(data,
                  specific_factors,
                  indicators  = NULL,
                  g_name      = "G",
                  estimator   = "MLR",
                  std.lv      = TRUE,
                  ordered     = NULL,
                  group       = NULL,
                  group_equal = NULL,
                  missing     = "listwise",
                  n_starts    = 30L,
                  ...) {

  mc <- match.call()

  # -- 0. Validate -----------------------------------------------------------
  if (!is.list(specific_factors) || is.null(names(specific_factors)))
    stop("`specific_factors` must be a named list.", call. = FALSE)

  if (g_name %in% names(specific_factors))
    stop("`g_name` ('", g_name, "') cannot be the same as a specific factor name.",
         call. = FALSE)

  if (is.null(indicators))
    indicators <- unlist(specific_factors, use.names = FALSE)

  missing_items <- setdiff(indicators, colnames(data))
  if (length(missing_items))
    stop("Items not found in data: ", paste(missing_items, collapse = ", "),
         call. = FALSE)

  # -- Ordered indicators: lavaan's efa() block does not support WLSMV ---------
  # Route automatically to besem_ordered() which implements the bifactor
  # set-ESEM pipeline (polychoric -> EFA -> WLSMV CFA), matching Mplus behaviour.
  if (!is.null(ordered)) {
    message("Ordinal indicators detected. lavaan's efa() block does not support ",
            "WLSMV rotation.\nRouting to besem_ordered() (bifactor set-ESEM + WLSMV pipeline).")
    return(besem_ordered(
      data             = data,
      specific_factors = specific_factors,
      indicators       = indicators,
      g_name           = g_name,
      group            = group,
      group_equal      = group_equal,
      missing          = missing,
      std.lv           = std.lv,
      ...
    ))
  }

  specific_names <- names(specific_factors)
  all_factor_names <- c(g_name, specific_names)
  nfactors_total   <- length(all_factor_names)

  # -- 1. Build bifactor target matrix ---------------------------------------
  # G: all items are primary (1) -- G loads freely on everything
  # Specific factors: primary items = 1, all others = 0
  tgt <- make_bifactor_target(
    specific_factors = specific_factors,
    indicators       = indicators,
    g_name           = g_name
  )

  # -- 2. Build lavaan efa() block syntax ------------------------------------
  syntax <- .build_besem_syntax(
    factor_names = all_factor_names,
    indicators   = indicators,
    g_name       = g_name
  )

  # -- 3. Fit via lavaan with orthogonal target rotation ---------------------
  # Convert target 1s -> NA for GPArotation: NA = free, 0 = target toward zero.
  tgt_na <- tgt; tgt_na[tgt == 1] <- NA_real_

  # Under FIML, the FIML ML step for B-ESEM with 4+ specific factors can
  # converge to a degenerate unrotated loading matrix (G absorbs specific-
  # factor variance).  Multi-start target rotation (n_starts) addresses
  # rotation-stage local optima, but cannot recover variance already lost
  # at the ML stage.  Honour the caller's choice and warn -- silently
  # switching to listwise here would desync B-ESEM from CFA / ESEM
  # (different Ns), invalidating the fit comparison.
  besem_missing <- missing
  if (is.null(ordered) && tolower(besem_missing) %in%
        c("ml", "ml.x", "fiml", "direct", "direct.ml") &&
      length(specific_factors) >= 4L) {
    warning(
      "B-ESEM under FIML with ", length(specific_factors), " specific factors: ",
      "multi-start rotation (n_starts = ", n_starts, ") addresses rotation ",
      "local optima, but the preceding FIML ML step can still converge to ",
      "a degenerate unrotated solution (G absorbs specific-factor variance).  ",
      "Inspect the rotated loadings and, if possible, verify against other ",
      "software such as Mplus.  Use missing = \"listwise\" to avoid the risk ",
      "(at the cost of complete-case deletion).",
      call. = FALSE
    )
  }

  cfa_args <- list(
    model         = syntax,
    data          = data,
    std.lv        = std.lv,
    estimator     = estimator,
    missing       = besem_missing,
    rotation      = "target",
    rotation.args = list(target    = tgt_na,
                         orthogonal = TRUE,
                         rstarts    = as.integer(n_starts))
  )
  if (!is.null(ordered))     cfa_args$ordered      <- ordered
  if (!is.null(group))       cfa_args$group        <- group
  if (!is.null(group_equal)) cfa_args$group.equal  <- group_equal
  cfa_args <- c(cfa_args, list(...))

  # Seed scoped to this call via withr (RNG state restored afterwards) so
  # random-start rotation is reproducible without disturbing the caller's
  # RNG stream.
  fit <- tryCatch(
    withr::with_seed(42L, do.call(lavaan::cfa, cfa_args)),
    error = function(e) {
      stop("lavaan::cfa() failed for B-ESEM:\n  ", conditionMessage(e),
           "\n\nCheck lavaan >= 0.6-12 and that indicators are correctly specified.",
           call. = FALSE)
    }
  )

  # -- 4. Cache standardized solution ----------------------------------------
  cached <- if (is.null(group)) .cache_std_loadings(fit) else NULL

  # Pull the unstandardized rotated loadings out of lavaan for B-EWC
  # (and any caller that needs raw lambda). These are on the y*/cov scale
  # the CFA target needs when pinning referent cross-loadings.
  L_rot_unstd <- if (is.null(group)) {
    lam <- tryCatch(lavaan::lavInspect(fit, "est")$lambda,
                    error = function(e) NULL)
    if (!is.null(lam) && !is.null(cached) &&
        all(rownames(cached$L) %in% rownames(lam)) &&
        all(colnames(cached$L) %in% colnames(lam))) {
      lam[rownames(cached$L), colnames(cached$L), drop = FALSE]
    } else lam
  } else NULL

  # Sign correction: flip columns so primary loadings sum positive.
  # Apply the same flips to the unstandardized matrix so the two stay in sync.
  if (!is.null(cached)) {
    if (!is.null(colnames(cached$L)) && !is.null(colnames(tgt)) &&
        all(colnames(tgt) %in% colnames(cached$L))) {
      cached$L  <- cached$L[, colnames(tgt), drop = FALSE]
      cached$SE <- cached$SE[, colnames(tgt), drop = FALSE]
      if (!is.null(L_rot_unstd) && !is.null(colnames(L_rot_unstd)) &&
          all(colnames(tgt) %in% colnames(L_rot_unstd)))
        L_rot_unstd <- L_rot_unstd[, colnames(tgt), drop = FALSE]
    }
    for (j in seq_len(ncol(cached$L))) {
      prim_idx <- if (j == 1L) seq_len(nrow(cached$L)) else which(tgt[, j] == 1)
      if (sum(cached$L[prim_idx, j], na.rm = TRUE) < 0) {
        cached$L[, j] <- -cached$L[, j]
        if (!is.null(L_rot_unstd) && j <= ncol(L_rot_unstd))
          L_rot_unstd[, j] <- -L_rot_unstd[, j]
      }
    }
  }

  # -- 5. Return -------------------------------------------------------------
  structure(
    list(
      lavaan_fit           = fit,
      syntax               = syntax,
      target               = tgt,
      nfactors             = nfactors_total,
      rotation             = "target (orthogonal)",
      factor_names         = all_factor_names,
      g_name               = g_name,
      specific_factors     = specific_factors,
      indicators           = indicators,
      rotated_loadings     = L_rot_unstd,
      std_rotated_loadings = if (!is.null(cached)) cached$L  else NULL,
      se_loadings          = if (!is.null(cached)) cached$SE else NULL,
      call                 = mc
    ),
    class = c("besem_fit", "esem_fit")
  )
}


#' Create a Bifactor Target Matrix
#'
#' Builds the target loading matrix for B-ESEM. The general factor G has
#' target 1 for all items (free to load everywhere). Each specific factor
#' has target 1 for its primary items and 0 for all others.
#'
#' @param specific_factors Named list of specific factor -> item assignments.
#' @param indicators Character vector of all item names.
#' @param g_name Character. Name for the general factor. Default \code{"G"}.
#'
#' @return A numeric matrix (items  x  factors) of class \code{"esem_target"}.
#'   Column order: G first, then specific factors in list order.
#'
#' @examples
#' tgt <- make_bifactor_target(
#'   specific_factors = list(EX = c("y1","y2","y3"),
#'                           MD = c("y4","y5","y6")),
#'   indicators       = paste0("y", 1:6)
#' )
#' print(tgt)
#'
#' @export
make_bifactor_target <- function(specific_factors,
                                  indicators,
                                  g_name = "G") {

  nfactors     <- length(specific_factors) + 1   # +1 for G
  factor_names <- c(g_name, names(specific_factors))
  nitems       <- length(indicators)

  tmat <- matrix(0, nrow = nitems, ncol = nfactors,
                 dimnames = list(indicators, factor_names))

  # G loads on everything -- target 1 for all items
  tmat[, g_name] <- 1

  # Specific factors -- 1 for primary items, 0 for cross-loadings
  for (f in names(specific_factors)) {
    primary_items <- specific_factors[[f]]
    bad <- setdiff(primary_items, indicators)
    if (length(bad))
      stop("Items in specific_factors[['", f, "']] not in indicators: ",
           paste(bad, collapse = ", "), call. = FALSE)
    tmat[primary_items, f] <- 1
  }

  class(tmat) <- c("esem_target", "matrix", "array")
  tmat
}


#' Generate Mplus B-ESEM Syntax
#'
#' Generates the correct Mplus syntax for Bifactor ESEM with orthogonal
#' target rotation, using separate BY statements per factor with (*1).
#'
#' @param specific_factors Named list of specific factor -> item assignments.
#' @param all_indicators Character vector of all ESEM indicator names.
#' @param g_name Character. General factor name. Default \code{"G"}.
#' @param cfa_factors Optional named list of additional CFA factors.
#' @param regressions Optional character vector of regression statements.
#' @param covariances Optional character vector of covariance statements.
#' @param data_file Character. Data file name. Default \code{"mydata.dat"}.
#' @param missing_code Numeric. Missing value code. Default \code{999}.
#' @param output_path Character. Path to write .inp file. Returns syntax
#'   invisibly if \code{NULL}.
#'
#' @return Mplus syntax string (invisibly). Writes file if \code{output_path}
#'   supplied.
#'
#' @export
generate_mplus_besem_syntax <- function(specific_factors,
                                         all_indicators,
                                         g_name       = "G",
                                         cfa_factors  = NULL,
                                         regressions  = NULL,
                                         covariances  = NULL,
                                         data_file    = "mydata.dat",
                                         missing_code = 999,
                                         output_path  = NULL) {

  specific_names <- names(specific_factors)

  # CFA items and covariates
  cfa_items  <- if (!is.null(cfa_factors)) unique(unlist(cfa_factors)) else character(0)
  all_factor_names <- c(g_name, specific_names, names(cfa_factors))

  covariates <- character(0)
  if (!is.null(regressions)) {
    for (reg in regressions) {
      parts <- strsplit(reg, " ON ")[[1]]
      if (length(parts) == 2) {
        rhs_vars <- trimws(strsplit(parts[2], " ")[[1]])
        rhs_vars <- rhs_vars[rhs_vars != ""]
        covariates <- unique(c(covariates, setdiff(rhs_vars, all_factor_names)))
      }
    }
  }

  all_vars <- c(all_indicators, cfa_items, covariates)

  # VARIABLE section
  names_line   <- .mplus_wrap("  NAMES = ", all_vars)
  usevars_line <- .mplus_wrap("  USEVARIABLES = ", all_vars)

  # G factor -- loads on all items, no target restrictions
  g_items_str <- paste(all_indicators, collapse = " ")
  g_block <- paste0("  ! General factor -- loads freely on all items\n",
                    "  ", g_name, " BY\n    ", g_items_str, " (*1);")

  # Specific factors -- same as regular ESEM target rotation
  specific_blocks <- vapply(specific_names, function(f) {
    primary <- specific_factors[[f]]
    cross   <- setdiff(all_indicators, primary)
    # Primary items first, then cross-loadings targeted to zero
    primary_str <- paste(primary, collapse = " ")
    cross_str   <- paste(paste0(cross, "~0"), collapse = " ")
    block <- paste0("  ", f, " BY\n    ", primary_str)
    if (nchar(cross_str) > 0)
      block <- paste0(block, "\n    ", cross_str, " (*1);")
    else
      block <- paste0(block, " (*1);")
    block
  }, FUN.VALUE = character(1))

  by_section <- paste(
    c(g_block, specific_blocks),
    collapse = "\n\n"
  )

  # CFA factors
  cfa_section <- ""
  if (!is.null(cfa_factors)) {
    cfa_blocks <- vapply(names(cfa_factors), function(f) {
      paste0("  ", f, " BY ", paste(cfa_factors[[f]], collapse = " "), ";")
    }, FUN.VALUE = character(1))
    cfa_section <- paste0("\n  ! CFA factors\n", paste(cfa_blocks, collapse = "\n"))
  }

  # Regressions
  reg_section <- ""
  if (!is.null(regressions))
    reg_section <- paste0("\n  ! Structural paths\n",
                          paste(paste0("  ", regressions, ";"), collapse = "\n"))

  # Covariances
  cov_section <- ""
  if (!is.null(covariances))
    cov_section <- paste0("\n  ! Residual covariances\n",
                          paste(paste0("  ", covariances, ";"), collapse = "\n"))

  syntax <- paste0(
    "TITLE: Bifactor ESEM (B-ESEM) with Orthogonal Target Rotation\n\n",
    "DATA:\n  FILE = ", data_file, ";\n\n",
    "VARIABLE:\n",
    names_line, "\n",
    usevars_line, "\n",
    "  MISSING ARE ALL (", missing_code, ");\n\n",
    "ANALYSIS:\n",
    "  ESTIMATOR = MLR;\n",
    "  ROTATION = TARGET (orthogonal); !Use for B-ESEM\n\n",
    "MODEL:\n",
    "  ! G loads on all items (no target restrictions)\n",
    "  ! Specific factors orthogonal to G and each other\n\n",
    by_section,
    cfa_section,
    reg_section,
    cov_section,
    "\n\nOUTPUT:\n  STDYX;\n  MODINDICES(10);\n"
  )

  if (!is.null(output_path)) {
    writeLines(syntax, output_path)
    message("Mplus B-ESEM syntax written to: ", output_path)
  }

  invisible(syntax)
}


# -- Print / Summary for besem_fit ---------------------------------------------

#' Print Method for besem_fit
#' @param x A \code{besem_fit} object.
#' @param ... Ignored.
#' @return Invisibly returns \code{x}; called for the side effect of printing
#'   a formatted B-ESEM model overview to the console.
#' @export
print.besem_fit <- function(x, ...) {
  cat("\n======================================================\n")
  cat(" bifactory: Bifactor ESEM (B-ESEM)\n")
  cat("======================================================\n\n")

  cat("Call:\n "); print(x$call); cat("\n")
  cat("General factor  :", x$g_name, "(loads on all", length(x$indicators), "items)\n")
  cat("Specific factors:", paste(names(x$specific_factors), collapse = ", "), "\n")
  cat("Rotation        : Orthogonal target (uncorrelated factors)\n")

  lav <- x$lavaan_fit
  if (!is.null(lav)) {
    converged <- lavaan::lavInspect(lav, "converged")
    cat("Converged       :", if (converged) "Yes" else "NO -- check model!", "\n")
  }

  # From-scratch WLSMV stats take precedence (ordered B-ESEM path)
  ws <- x$wlsmv_stats
  if (!is.null(ws)) {
    cat("\nModel Fit (WLSMV -- Asparouhov & Muthen 2009):\n")
    cat(sprintf("  X2(%g) = %.3f, p = %.4f\n", ws$df, ws$chisq, ws$pvalue))
    cat(sprintf("  CFI = %.3f  |  TLI = %.3f\n", ws$cfi, ws$tli))
    cat(sprintf("  RMSEA = %.3f\n", ws$rmsea))
    cat(sprintf("  SRMR = %.3f\n", ws$srmr))
  } else if (!is.null(lav)) {
    fi <- .fit_indices_lav(lav)
    if (!is.null(fi)) {
      cat("\nModel Fit:\n")
      cat(sprintf("  X2(%g) = %.3f, p = %.3f\n", fi["df"], fi["chisq"], fi["pvalue"]))
      cat(sprintf("  CFI = %.3f  |  TLI = %.3f\n", fi["cfi"], fi["tli"]))
      cat(sprintf("  RMSEA = %.3f [%.3f, %.3f]\n",
                  fi["rmsea"], fi["rmsea.ci.lower"], fi["rmsea.ci.upper"]))
      cat(sprintf("  SRMR = %.3f\n", fi["srmr"]))
    }
  }

  cat("\nStandardised Loadings:\n")
  sl <- tryCatch(std_loadings(x), error = function(e) NULL)
  if (!is.null(sl)) {
    # Show G first, then specific
    col_order <- intersect(c(x$g_name, names(x$specific_factors)), colnames(sl))
    print(round(sl[, col_order, drop = FALSE], 3))
  }

  cat("\nNote: All factor correlations are fixed to 0 (orthogonal bifactor).\n")
  cat("Use summary() for parameter estimates and full output.\n\n")
  invisible(x)
}


# -- Internal syntax builder ---------------------------------------------------

.build_besem_syntax <- function(factor_names, indicators, g_name) {
  # Same efa() block structure as regular ESEM but with all factors
  # (G + specific) sharing one block and orthogonal rotation applied externally
  lhs <- paste(
    paste0('efa("besem")*', factor_names),
    collapse = " +\n    "
  )
  chunks   <- split(indicators, ceiling(seq_along(indicators) / 6))
  rhs_rows <- vapply(chunks, paste, character(1L), collapse = " + ")
  rhs      <- paste(rhs_rows, collapse = " +\n    ")

  paste0(
    "# bifactory: Bifactor ESEM (B-ESEM)\n",
    "# General factor '", g_name, "' + ", length(factor_names) - 1,
    " specific factors\n",
    "# Orthogonal target rotation = Mplus TARGET(orthogonal)\n\n",
    lhs, " =~\n    ", rhs, "\n"
  )
}
