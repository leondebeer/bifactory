# -- S3 Methods for esem_fit ---------------------------------------------------
# All standard lavaan generics are forwarded to the embedded lavaan_fit object
# so users rarely need to extract it manually.


# Internal: scaled-aware fit-index extractor used by print methods so the
# values they show line up with fit_indices() / .build_comparison_table().
# For WLSMV/MLR/MLM/DWLS lavaan reports both naive and scaled versions; the
# scaled ones are the Mplus-comparable, correctly-corrected statistics.
# Also applies the n_pairs/n_wls SRMR correction (Mplus denominator includes
# threshold residuals, which lavaan omits).
#
# Returns a named numeric vector with chisq, df, pvalue, cfi, tli, rmsea,
# rmsea.ci.lower, rmsea.ci.upper, srmr, npar -- or NULL if extraction fails.
.fit_indices_lav <- function(lav) {
  est <- tryCatch(lavaan::lavInspect(lav, "options")$estimator,
                  error = function(e) "")
  has_scaled <- grepl("DWLS|WLSMV|MLR|MLM|WLSM", est, ignore.case = TRUE)

  keys <- if (has_scaled)
    c(npar = "npar",
      chisq = "chisq.scaled", df = "df.scaled", pvalue = "pvalue.scaled",
      cfi = "cfi.scaled", tli = "tli.scaled", rmsea = "rmsea.scaled",
      rmsea.ci.lower = "rmsea.scaled.ci.lower",
      rmsea.ci.upper = "rmsea.scaled.ci.upper",
      srmr = "srmr")
  else
    c(npar = "npar",
      chisq = "chisq", df = "df", pvalue = "pvalue",
      cfi = "cfi", tli = "tli", rmsea = "rmsea",
      rmsea.ci.lower = "rmsea.ci.lower",
      rmsea.ci.upper = "rmsea.ci.upper",
      srmr = "srmr")

  fm <- tryCatch(lavaan::fitMeasures(lav), error = function(e) NULL)
  if (is.null(fm)) return(NULL)

  out <- vapply(keys, function(k)
    if (k %in% names(fm)) unname(fm[[k]]) else NA_real_, 0.0)
  names(out) <- names(keys)

  # Some lavaan versions emit `rmsea.scaled.ci.lower` whereas earlier ones
  # used `rmsea.ci.lower.scaled`; check the alternate spelling on miss.
  if (has_scaled) {
    if (is.na(out["rmsea.ci.lower"]) && "rmsea.ci.lower.scaled" %in% names(fm))
      out["rmsea.ci.lower"] <- unname(fm[["rmsea.ci.lower.scaled"]])
    if (is.na(out["rmsea.ci.upper"]) && "rmsea.ci.upper.scaled" %in% names(fm))
      out["rmsea.ci.upper"] <- unname(fm[["rmsea.ci.upper.scaled"]])
  }

  if (has_scaled && !is.na(out["srmr"])) {
    wls_obs <- tryCatch(lavaan::lavInspect(lav, "wls.obs"),
                        error = function(e) NULL)
    if (!is.null(wls_obs)) {
      if (is.list(wls_obs)) wls_obs <- wls_obs[[1L]]
      n_wls  <- length(wls_obs)
      cor_ov <- tryCatch(lavaan::lavInspect(lav, "cor.ov"),
                         error = function(e) NULL)
      if (!is.null(cor_ov)) {
        if (is.list(cor_ov)) cor_ov <- cor_ov[[1L]]
        n_items <- nrow(cor_ov)
        n_pairs <- n_items * (n_items - 1L) / 2L
        if (n_wls > n_pairs) out["srmr"] <- out["srmr"] * sqrt(n_pairs / n_wls)
      }
    }
  }
  out
}


#' Print Method for esem_fit
#'
#' @param x An \code{esem_fit} object.
#' @param ... Ignored.
#' @return Invisibly returns \code{x}; called for the side effect of printing
#'   a formatted model overview to the console.
#' @export
print.esem_fit <- function(x, ...) {
  cat("\n===============================================\n")
  cat(" bifactory: Exploratory Structural Equation Model\n")
  cat("===============================================\n\n")

  cat("Call:\n ")
  print(x$call)
  cat("\n")

  cat("Factors       :", x$nfactors, "\n")
  cat("Factor names  :", paste(x$factor_names, collapse = ", "), "\n")
  cat("Rotation      :", x$rotation, "\n")

  # Lavaan convergence
  lav <- x$lavaan_fit
  converged <- lavaan::lavInspect(lav, "converged")
  cat("Converged     :", if (converged) "Yes" else "NO (check model!)", "\n")

  # Key fit indices -- prefer custom WLSMV stats (B-ESEM DWLS from-scratch)
  ws <- x$wlsmv_stats
  if (!is.null(ws)) {
    cat("\nModel Fit (WLSMV -- Asparouhov & Muthen 2009):\n")
    cat(sprintf("  X2(%g) = %.3f, p = %.4f\n", ws$df, ws$chisq, ws$pvalue))
    cat(sprintf("  CFI = %.3f  |  TLI = %.3f\n", ws$cfi, ws$tli))
    cat(sprintf("  RMSEA = %.3f\n", ws$rmsea))
    cat(sprintf("  SRMR = %.3f\n", ws$srmr))
  } else {
    fi <- .fit_indices_lav(lav)
    if (!is.null(fi)) {
      cat("\nModel Fit:\n")
      cat(sprintf("  X2(%g) = %.3f, p = %.3f\n", fi["df"], fi["chisq"], fi["pvalue"]))
      cat(sprintf("  CFI = %.3f  |  TLI = %.3f\n", fi["cfi"], fi["tli"]))
      cat(sprintf("  RMSEA = %.3f [%.3f, %.3f]\n",
                  fi["rmsea"], fi["rmsea.ci.lower"], fi["rmsea.ci.upper"]))
      cat(sprintf("  SRMR = %.3f\n", fi["srmr"]))
      cat(sprintf("  Free parameters: %g\n", fi["npar"]))
    }
  }

  cat("\nStandardised Loadings (STDYX):\n")
  sl <- tryCatch(std_loadings(x, suppress = 0), error = function(e) NULL)
  if (!is.null(sl)) print(round(sl, 3))

  # Heywood correction summary
  hl <- x$heywood_log
  if (!is.null(hl)) {
    cat("\n-- Heywood correction applied ")
    cat(strrep("-", 38), "\n", sep = "")
    if (!is.null(hl$final_loadings)) {
      for (i in seq_len(nrow(hl$final_loadings))) {
        r      <- hl$final_loadings[i, ]
        tg     <- hl$final_targets[r$item]
        rd     <- if (!is.na(r$resolved_round))
                    sprintf("resolved in round %d", r$resolved_round) else "unresolved"
        if (is.na(tg)) {
          cat(sprintf("  %-12s (%-5s): %s -- final target (unknown), lam = %.3f\n",
                      r$item, r$factor, rd, r$lambda))
        } else {
          cat(sprintf("  %-12s (%-5s): %s -- final target %.2f, lam = %.3f\n",
                      r$item, r$factor, rd, tg, r$lambda))
        }
      }
    }
    if (!is.null(hl$delta_cfi) && !is.na(hl$delta_cfi))
      cat(sprintf("  dCFI (corrected - initial): %+.3f\n", hl$delta_cfi))
    if (!is.null(hl$original_rotation) &&
        !identical(hl$original_rotation, hl$final_rotation))
      cat(sprintf("  Note: rotation changed from '%s' to '%s' during correction.\n",
                  hl$original_rotation, hl$final_rotation))
    if (!hl$resolved)
      cat("  WARNING: one or more loadings could not be fully resolved.\n")
    cat(sprintf(
      "  Report: target rotation applied with item-specific targets as listed above (%d round(s)).\n",
      hl$rounds
    ))
    cat(strrep("-", 68), "\n", sep = "")
  }

  cat("\nUse summary() for parameter estimates and full output.\n")
  cat("Use fitMeasures(), modindices(), or lavaan_fit() for more.\n\n")

  invisible(x)
}


#' Summary Method for esem_fit
#'
#' Prints a full summary for an \code{esem_fit}.  For models fit through the
#' DWLS-from-scratch path (B-ESEM ordered, or \code{esem_ordered(method =
#' "rotation")}), the underlying \code{lavaan_fit} slot holds an auxiliary
#' 1-factor CFA used only for weight-matrix extraction; in that case we print
#' a self-contained bifactory summary built from the stored rotated
#' loadings, standard errors, and WLSMV statistics.  Otherwise the call is
#' forwarded to \code{lavaan::summary()}.
#'
#' @param object An \code{esem_fit} object.
#' @param fit.measures Logical. Include fit indices? Default \code{TRUE}.
#' @param standardized Logical. Include standardized solution? Default \code{TRUE}.
#' @param rsquare Logical. Include R-squared for endogenous variables? Default \code{FALSE}.
#' @param ... Additional arguments passed to \code{lavaan::summary()}.
#'
#' @return Invisibly returns \code{object}; called for the side effect of
#'   printing the model summary.
#' @export
summary.esem_fit <- function(object,
                              fit.measures = TRUE,
                              standardized = TRUE,
                              rsquare      = FALSE,
                              ...) {

  # DWLS from-scratch path: lavaan_fit is an auxiliary 1-factor CFA used only
  # for weight extraction.  Forwarding to lavaan::summary() would summarise
  # that aux model and produce no useful information about the actual ESEM /
  # B-ESEM solution -- so build the summary from our stored slots instead.
  if (!is.null(object$wlsmv_stats))
    return(.summary_from_scratch(object))

  cat("\n======================================================\n")
  cat(" bifactory ESEM Summary\n")
  cat("======================================================\n\n")

  cat("Rotation:", object$rotation, "| Factors:", object$nfactors, "\n\n")

  cat("--- ESEM Solution (lavaan native efa() block) -------- \n")
  lavaan::summary(object$lavaan_fit,
                  fit.measures = fit.measures,
                  standardized = standardized,
                  rsquare      = rsquare,
                  ...)

  invisible(object)
}


.summary_from_scratch <- function(object) {
  is_besem <- inherits(object, "besem_fit")
  banner   <- if (is_besem)
    "Bifactor ESEM (B-ESEM) Summary -- DWLS from-scratch path"
  else
    "ESEM Summary -- DWLS from-scratch path"

  cat("\n======================================================\n")
  cat(" bifactory ", banner, "\n", sep = "")
  cat("======================================================\n\n")
  cat("Estimator    : WLSMV (Asparouhov & Muthen 2009)\n")
  cat("Rotation     : ", object$rotation, "\n", sep = "")
  cat("Factors      : ", object$nfactors %||% length(object$factor_names), "\n", sep = "")
  if (is_besem) {
    cat("General      : ", object$g_name, " (loads on all ",
        length(object$indicators), " items)\n", sep = "")
    cat("Specific     : ",
        paste(names(object$specific_factors), collapse = ", "), "\n", sep = "")
  } else {
    cat("Factor names : ",
        paste(object$factor_names, collapse = ", "), "\n", sep = "")
  }

  ws <- object$wlsmv_stats
  cat("\nModel Fit (mean- and variance-adjusted chi-square):\n")
  cat(sprintf("  X2(%g) = %.3f, p = %.4f\n", ws$df, ws$chisq, ws$pvalue))
  cat(sprintf("  CFI = %.3f  |  TLI = %.3f\n", ws$cfi, ws$tli))
  cat(sprintf("  RMSEA = %.3f\n", ws$rmsea))
  cat(sprintf("  SRMR = %.3f\n", ws$srmr))

  cat("\nStandardised Loadings (STDYX) with rotation-corrected SEs:\n")
  parameters(object, type = "loadings")

  if (is_besem)
    cat("\nNote: All factor correlations fixed to 0 (orthogonal bifactor).\n")

  invisible(object)
}


#' Extract Fit Measures from an esem_fit
#'
#' Forwards to \code{lavaan::fitMeasures()} when the fit used lavaan's \code{efa()}
#' path. When \code{wlsmv_stats} is present (custom DWLS/WLSMV from
#' \code{\link{besem_ordered}(method = "rotation")} or
#' \code{\link{esem_ordered}(method = "rotation")}), returns those indices instead.
#'
#' @param object An \code{esem_fit} object.
#' @param fit.measures Character vector of fit index names. Default returns a
#'   standard set: CFI, TLI, RMSEA, SRMR, AIC, BIC.
#' @param ... Passed to \code{lavaan::fitMeasures()}.
#'
#' @return A named numeric vector of fit indices (by default CFI, TLI, RMSEA
#'   with its confidence interval, SRMR, and, on lavaan paths, AIC and BIC).
#' @exportS3Method lavaan::fitMeasures
fitMeasures.esem_fit <- function(object,
                                  fit.measures = c("cfi", "tli", "rmsea",
                                                   "rmsea.ci.lower",
                                                   "rmsea.ci.upper",
                                                   "srmr", "aic", "bic"),
                                  ...) {
  # B-ESEM DWLS from-scratch: return custom WLSMV stats
  if (!is.null(object$wlsmv_stats)) {
    ws  <- object$wlsmv_stats
    all <- c(cfi   = ws$cfi,   tli    = ws$tli,   rmsea  = ws$rmsea,
             srmr  = ws$srmr,  chisq  = ws$chisq, df     = ws$df,
             pvalue = ws$pvalue)
    avail <- intersect(fit.measures, names(all))
    return(all[avail])
  }
  lavaan::fitMeasures(object$lavaan_fit, fit.measures = fit.measures, ...)
}


#' Extract Modification Indices from an esem_fit
#'
#' Forwards to \code{lavaan::modindices()} on \code{lavaan_fit(x)}. For custom
#' WLSMV paths where \code{lavaan_fit} is an auxiliary one-factor CFA, indices
#' refer to that auxiliary model, not the rotated ESEM/B-ESEM solution; see
#' \link[=doc_estimator_paths]{The lavaan_fit slot on custom WLSMV fits}.
#'
#' @param object An \code{esem_fit} object.
#' @param sort. Logical. Sort by modification index value? Default \code{TRUE}.
#' @param maximum.number Integer. Maximum number of indices to return. Default 20.
#' @param ... Passed to \code{lavaan::modindices()}.
#'
#' @return A \code{data.frame} of modification indices, as returned by
#'   \code{lavaan::modindices()}.
#' @exportS3Method lavaan::modindices
modindices.esem_fit <- function(object,
                                 sort.           = TRUE,
                                 maximum.number  = 20,
                                 ...) {
  lavaan::modindices(object$lavaan_fit,
                     sort.           = sort.,
                     maximum.number  = maximum.number,
                     ...)
}


#' Extract Parameter Estimates from an esem_fit
#'
#' Forwards to lavaan on \code{lavaan_fit(x)}. On custom WLSMV rotation paths,
#' use \code{\link{parameters}} or \code{\link{std_loadings}} instead; see
#' \link[=doc_estimator_paths]{The lavaan_fit slot on custom WLSMV fits}.
#'
#' @param object An \code{esem_fit} object.
#' @param standardized Logical. Return standardized estimates? Default \code{FALSE}.
#' @param ... Passed to \code{lavaan::parameterEstimates()}.
#'
#' @return A \code{data.frame} of parameter estimates from
#'   \code{lavaan::parameterEstimates()}, or the standardized solution from
#'   \code{lavaan::standardizedsolution()} when \code{standardized = TRUE}.
#' @export
coef.esem_fit <- function(object, standardized = FALSE, ...) {
  if (standardized) {
    lavaan::standardizedsolution(object$lavaan_fit, ...)
  } else {
    lavaan::parameterEstimates(object$lavaan_fit, ...)
  }
}


#' Display Model Parameters
#'
#' Prints a formatted parameter table -- standardized loadings, standard errors,
#' z-values, and p-values -- for one fit object or all three models in a pipeline.
#' Analogous to the STDYX section of Mplus output.
#'
#' @param x An \code{esem_fit}, \code{besem_fit}, \code{ewc_fit}, or
#'   \code{esem_comparison_pipeline} object. To inspect parameters at a
#'   specific level of an \code{esem_invariance} result, pass the level's
#'   fit object directly (e.g. \code{parameters(inv$models$strict)}).
#' @param model Character. When \code{x} is a pipeline, which models to show:
#'   \code{"all"} (default), \code{"CFA"}, \code{"ESEM"}, or \code{"BESEM"}.
#' @param type Character. \code{"loadings"} (default) shows only factor loadings;
#'   \code{"all"} also shows residual variances and, for ordered models, thresholds.
#' @param suppress Numeric. Hide loadings with |z| below this value. Default 0
#'   (show all).
#' @param digits Integer. Decimal places. Default 3.
#' @param highlight_primary Logical. Colour the target (primary) loadings using
#'   ANSI codes. Supported in RStudio and most terminals; falls back silently to
#'   plain output when colour is unavailable. Default \code{TRUE}.
#'
#' @return Invisibly returns the parameter table as a data frame.
#' @export
parameters <- function(x,
                       model             = "all",
                       type              = "loadings",
                       suppress          = 0,
                       digits            = 3,
                       highlight_primary = TRUE) {

  .sig_stars <- function(p) {
    ifelse(is.na(p), "   ",
    ifelse(p < .001, "***",
    ifelse(p < .01,  "** ",
    ifelse(p < .05,  "*  ",
    ifelse(p < .10,  ".  ", "   ")))))
  }

  # Determine primary items per factor for a fit object.
  # factor_items: optional named list (factor -> character vector of items),
  #   passed in from spec$factors for CFA/ESEM models in a pipeline.
  # Returns NULL when no mapping is available (loading_type will be "cross"
  #   everywhere, i.e., no colour highlighting).
  .primary_map <- function(fit_obj, factor_items = NULL) {
    # Guard: lavaan S4 objects don't support $ -- use factor_items directly
    if (isS4(fit_obj)) {
      if (!is.null(factor_items))
        return(setNames(lapply(factor_items, function(i) tolower(i)),
                        tolower(names(factor_items))))
      return(NULL)
    }
    # B-ESEM bifactor structure (specific_factors + g_name present)
    if (!is.null(fit_obj$specific_factors) && !is.null(fit_obj$g_name)) {
      gname <- tolower(fit_obj$g_name)
      # Both keys AND values must be lowercase so tolower(fac) lookups work
      specs <- setNames(
        lapply(fit_obj$specific_factors, function(items) tolower(items)),
        tolower(names(fit_obj$specific_factors))
      )
      all_items <- tolower(unlist(fit_obj$specific_factors, use.names = FALSE))
      return(c(setNames(list(all_items), gname), specs))
    }
    # CFA / ESEM without bifactor structure: use factor_items if supplied
    if (!is.null(factor_items))
      return(setNames(lapply(factor_items, function(i) tolower(i)),
                      tolower(names(factor_items))))
    NULL   # no map available -- nothing will be highlighted
  }

  # Extract a clean parameter table from one fit object.
  # factor_items: optional named list for primary-loading colouring in CFA/ESEM.
  .params_one <- function(fit_obj, model_label, factor_items = NULL) {

    pmap <- .primary_map(fit_obj, factor_items)

    # B-ESEM rotation path: use rotation matrices directly
    if (inherits(fit_obj, "esem_fit") &&
        !is.null(fit_obj$std_rotated_loadings) &&
        !is.null(fit_obj$se_loadings)) {
      L  <- fit_obj$std_rotated_loadings
      SE <- fit_obj$se_loadings
      z  <- L / SE
      p  <- 2 * pnorm(-abs(z))
      out <- do.call(rbind, lapply(colnames(L), function(fac)
        data.frame(factor       = fac,
                   item         = rownames(L),
                   std          = round(L[, fac],  digits),
                   se           = round(SE[, fac], digits),
                   z            = round(z[, fac],  2),
                   p            = round(p[, fac],  4),
                   loading_type = if (!is.null(pmap[[tolower(fac)]])) {
                     ifelse(tolower(rownames(L)) %in% pmap[[tolower(fac)]],
                            "primary", "cross")
                   } else "cross",
                   stringsAsFactors = FALSE)))
      if (suppress > 0) out <- out[abs(out$z) >= suppress, , drop = FALSE]
      return(out)
    }

    # ESEM efa() path with Heywood correction: use corrected efa_loadings
    # with SEs from lavaan's standardizedSolution (rotation is post-hoc so
    # SEs apply to the original loadings, not the corrected ones -- but they
    # remain the best available estimate).
    if (inherits(fit_obj, "esem_fit") && !is.null(fit_obj$efa_loadings) &&
        !is.null(fit_obj$heywood_log) && type == "loadings") {
      L   <- fit_obj$efa_loadings
      lav <- fit_obj$lavaan_fit
      ss  <- lavaan::standardizedsolution(lav)
      ss  <- ss[ss$op == "=~", , drop = FALSE]
      se_lkp <- setNames(ss$se,     paste0(ss$lhs, ".", ss$rhs))
      z_lkp  <- setNames(ss$z,      paste0(ss$lhs, ".", ss$rhs))
      p_lkp  <- setNames(ss$pvalue, paste0(ss$lhs, ".", ss$rhs))
      out <- do.call(rbind, lapply(colnames(L), function(fac) {
        keys <- paste0(fac, ".", rownames(L))
        ltype <- if (!is.null(pmap[[tolower(fac)]])) {
          ifelse(tolower(rownames(L)) %in% pmap[[tolower(fac)]], "primary", "cross")
        } else "cross"
        data.frame(factor       = fac,
                   item         = rownames(L),
                   std          = round(L[, fac], digits),
                   se           = round(se_lkp[keys], digits),
                   z            = round(L[, fac] / se_lkp[keys], 2),
                   p            = round(2 * pnorm(-abs(L[, fac] / se_lkp[keys])), 4),
                   loading_type = ltype,
                   stringsAsFactors = FALSE)
      }))
      if (suppress > 0) out <- out[!is.na(out$z) & abs(out$z) >= suppress, , drop = FALSE]
      return(out)
    }

    # lavaan path: CFA, ESEM (no Heywood), continuous BESEM
    lav <- if (inherits(fit_obj, "esem_fit")) fit_obj$lavaan_fit else fit_obj
    ss  <- lavaan::standardizedsolution(lav)

    if (type == "loadings") {
      ss <- ss[ss$op == "=~", , drop = FALSE]
    } else {
      ss <- ss[ss$op %in% c("=~", "~~", "|"), , drop = FALSE]
    }

    fac_col  <- ifelse(ss$op == "=~", ss$lhs, paste0(ss$lhs, ss$op, ss$rhs))
    item_col <- ifelse(ss$op == "=~", ss$rhs, "")
    ltype    <- mapply(function(fac, item) {
      if (is.null(pmap)) return("cross")   # no map -> no highlighting
      pri <- pmap[[tolower(fac)]]
      if (is.null(pri)) return("cross")
      if (tolower(item) %in% pri) "primary" else "cross"
    }, fac_col, item_col)

    out <- data.frame(
      factor       = fac_col,
      item         = item_col,
      std          = round(ss$est.std, digits),
      se           = round(ss$se,      digits),
      z            = round(ss$z,       2),
      p            = round(ss$pvalue,  4),
      loading_type = ltype,
      stringsAsFactors = FALSE
    )
    if (suppress > 0) out <- out[!is.na(out$z) & abs(out$z) >= suppress, , drop = FALSE]
    # Within each factor: primary (target) loadings first, crosses after.
    # order() is stable so original item sequence is preserved inside each group.
    if ("loading_type" %in% names(out) && nrow(out) > 0L) {
      out <- out[order(match(out$factor, unique(out$factor)),
                       out$loading_type != "primary"), , drop = FALSE]
    }
    out
  }

  .print_block <- function(tbl, model_label, fit_obj) {
    cat(sprintf("\n======================================================\n"))
    cat(sprintf(" %s -- Standardized Parameters (STDYX)\n", model_label))
    cat(sprintf("======================================================\n"))

    # Print fit indices on one line
    ws <- if (inherits(fit_obj, "esem_fit")) fit_obj$wlsmv_stats else NULL
    if (!is.null(ws)) {
      cat(sprintf(" X2(%g)=%.3f  CFI=%.3f  TLI=%.3f  RMSEA=%.3f  SRMR=%.3f\n",
                  ws$df, ws$chisq, ws$cfi, ws$tli, ws$rmsea, ws$srmr))
    } else if (inherits(fit_obj, "esem_fit") || inherits(fit_obj, "lavaan")) {
      lav <- if (inherits(fit_obj, "esem_fit")) fit_obj$lavaan_fit else fit_obj
      # WLSMV/DWLS/MLR -> .scaled variants (Mplus-matched); plain ML -> plain.
      est_ok <- tryCatch(lavaan::lavInspect(lav, "options")$estimator,
                         error = function(e) "")
      has_scaled <- grepl("DWLS|WLSMV|MLR|MLM|WLSM", est_ok, ignore.case = TRUE)
      keys <- if (has_scaled)
        c("chisq.scaled","df.scaled","cfi.scaled","tli.scaled","rmsea.scaled","srmr")
      else
        c("chisq","df","cfi","tli","rmsea","srmr")
      fi <- tryCatch(lavaan::fitMeasures(lav, keys), error = function(e) NULL)
      if (!is.null(fi))
        cat(sprintf(" X2(%g)=%.3f  CFI=%.3f  TLI=%.3f  RMSEA=%.3f  SRMR=%.3f\n",
                    fi[2], fi[1], fi[3], fi[4], fi[5], fi[6]))
    }

    cat(sprintf("\n  %-10s  %-12s  %*s  %*s  %6s  %6s\n",
                "Factor", "Item",
                digits + 4, "Std.Est",
                digits + 3, "SE",
                "z", "p"))
    cat(sprintf("  %-10s  %-12s  %*s  %*s  %6s  %6s\n",
                "----------", "------------",
                digits + 4, "-------",
                digits + 3, "-----",
                "------", "------"))

    # Detect ANSI colour support: RStudio sets RSTUDIO env var; also check TERM
    use_colour <- highlight_primary && (
      nzchar(Sys.getenv("RSTUDIO")) ||
      (!identical(Sys.getenv("TERM"), "dumb") && nzchar(Sys.getenv("TERM"))) ||
      isTRUE(getOption("bifactory.ansi_colour"))
    )
    col_on  <- if (use_colour) "\033[34m" else ""  # blue for primary loadings
    col_off <- if (use_colour) "\033[0m"  else ""

    for (fac in unique(tbl$factor)) {
      rows <- tbl[tbl$factor == fac, , drop = FALSE]
      first <- TRUE
      for (i in seq_len(nrow(rows))) {
        r      <- rows[i, ]
        is_pri <- !is.null(r$loading_type) && !is.na(r$loading_type) &&
                  r$loading_type == "primary"
        c1 <- if (is_pri) col_on  else ""
        c0 <- if (is_pri) col_off else ""
        stars <- .sig_stars(r$p)
        p_str <- if (is.na(r$p)) "      " else
                 if (r$p < .001) " <.001" else sprintf(" %.4f", r$p)
        cat(sprintf("  %s%-10s  %-12s  %*.3f  %*.3f  %6.2f  %s %s%s\n",
                    c1,
                    if (first) fac else "",
                    r$item,
                    digits + 4, r$std,
                    digits + 3, r$se,
                    r$z,
                    p_str, stars,
                    c0))
        first <- FALSE
      }
      cat("\n")
    }
    invisible(tbl)
  }

  # -- Dispatch ----------------------------------------------------------------
  if (inherits(x, "esem_comparison_pipeline")) {
    show <- if (model == "all") c("CFA", "ESEM", "BESEM") else toupper(model)
    fits <- list(CFA   = x$fit_cfa,
                 ESEM  = x$fit_esem,
                 BESEM = x$fit_besem)
    # spec$factors provides the factor->items map for CFA and ESEM highlighting;
    # BESEM uses its own specific_factors field and ignores this argument.
    fi <- if (!is.null(x$spec)) x$spec$factors else NULL
    all_tbls <- lapply(show, function(m) {
      tbl <- .params_one(fits[[m]], m, factor_items = fi)
      .print_block(tbl, m, fits[[m]])
    })
    invisible(do.call(rbind, mapply(function(t, m) { t$model <- m; t },
                                   all_tbls, show, SIMPLIFY = FALSE)))
  } else if (inherits(x, "esem_fit")) {
    lbl <- if (!is.null(x$call)) deparse(x$call[[1]]) else "Model"
    # Spec attached by run_comparison() -> use for primary/cross highlighting.
    # B-ESEM ignores factor_items (uses its own specific_factors instead).
    fi  <- if (!is.null(x$spec)) x$spec$factors else NULL
    tbl <- .params_one(x, lbl, factor_items = fi)
    .print_block(tbl, lbl, x)
    invisible(tbl)
  } else if (inherits(x, "ewc_fit")) {
    # EWC: raw lavaan S4 inside; use spec$factors as the primary-item map so
    # target loadings get the same blue highlight as ESEM output.
    fi  <- if (!is.null(x$spec)) x$spec$factors else NULL
    ewc_label <- if (isTRUE(x$besem)) "B-EWC" else "EWC"
    tbl <- .params_one(x$lavaan_fit, ewc_label, factor_items = fi)
    .print_block(tbl, ewc_label, x$lavaan_fit)
    invisible(tbl)
  } else {
    stop("`x` must be an esem_fit, ewc_fit, or esem_comparison_pipeline ",
         "object.", call. = FALSE)
  }
}


#' Extract the Underlying lavaan Fit Object
#'
#' Returns the raw \code{lavaan} object stored in \code{x$lavaan_fit}. For
#' continuous ESEM/B-ESEM and \code{method = "lavaan"} ordered fits, this is the
#' fitted model. For \code{\link{besem_ordered}(method = "rotation")} (and
#' \code{\link{esem_ordered}(method = "rotation")}), it is often an **auxiliary**
#' one-factor CFA used for DWLS weights only---see
#' \link[=doc_estimator_paths]{The lavaan_fit slot on custom WLSMV fits}.
#' Use \code{lavaan::lavTestScore()}, \code{lavaan::lavInspect()}, etc. only when
#' that slot represents the model you intend to analyse.
#'
#' @param x An \code{esem_fit} object.
#'
#' @return A \code{lavaan} object.
#'
#' @examples
#' \dontrun{
#' fit <- esem(mydata, nfactors = 3)
#' lav <- lavaan_fit(fit)
#'
#' # Use semPlot
#' semPlot::semPaths(lav, whatLabels = "std", layout = "tree")
#'
#' # Use lavaan::lavInspect
#' lavaan::lavInspect(lav, "cor.lv")   # factor correlations
#' }
#'
#' @export
lavaan_fit <- function(x) {
  if (!inherits(x, "esem_fit"))
    stop("`x` must be an esem_fit object.", call. = FALSE)
  x$lavaan_fit
}


#' Extract the Generated lavaan Syntax
#'
#' @param x An \code{esem_fit} object.
#' @param cat Logical. If \code{TRUE} (default), also print to console.
#'
#' @return The lavaan model string (invisibly).
#' @export
get_syntax <- function(x, cat = TRUE) {
  if (!inherits(x, "esem_fit"))
    stop("`x` must be an esem_fit object.", call. = FALSE)
  if (cat) base::cat(x$syntax, "\n")
  invisible(x$syntax)
}


#' Extract Factor Correlations from an esem_fit
#'
#' Extracts the estimated correlation matrix among latent factors.
#'
#' @param x An \code{esem_fit} object.
#' @param digits Integer. Number of decimal places. Default 3.
#'
#' @return A symmetric matrix of factor correlations.
#' @export
factor_correlations <- function(x, digits = 3) {
  if (!inherits(x, "esem_fit"))
    stop("`x` must be an esem_fit object.", call. = FALSE)

  # DWLS path (esem_ordered_dwls / besem_ordered rotation): Phi stored directly.
  if (!is.null(x$factor_correlations) && is.matrix(x$factor_correlations))
    return(round(x$factor_correlations, digits))

  cor_mat <- lavaan::lavInspect(x$lavaan_fit, "cor.lv")
  if (is.list(cor_mat)) cor_mat <- cor_mat[[1L]]   # multi-group: use group 1
  # Restore factor names if lavInspect stripped them (e.g. WLSMV EFA)
  if (is.null(rownames(cor_mat)) && !is.null(x$factor_names)) {
    fn <- x$factor_names
    if (nrow(cor_mat) == length(fn)) dimnames(cor_mat) <- list(fn, fn)
  }
  round(cor_mat, digits)
}


#' Extract Standardized Loadings from an esem_fit
#'
#' Returns the standardized factor loading matrix. Uses \code{std_rotated_loadings}
#' when present (custom WLSMV rotation path or Heywood-corrected loadings);
#' otherwise extracts from \code{lavaan_fit(x)}. Analogous to Mplus STDYX loadings.
#'
#' @param x An \code{esem_fit} object.
#' @param digits Integer. Rounding digits. Default 3.
#' @param suppress Numeric. Loadings with absolute value below this threshold
#'   are replaced with \code{NA} for readability. Default \code{0} (show all).
#'
#' @return A matrix of standardized loadings (items x factors).
#' @export
std_loadings <- function(x, digits = 3, suppress = 0) {
  if (!inherits(x, "esem_fit"))
    stop("`x` must be an esem_fit object.", call. = FALSE)

  # Rotation method: return the post-hoc rotated + standardised loading matrix
  # (WLSMV unrestricted model with bifactor targetT rotation + STDYX scaling).
  # The lavaan_fit for this method holds the echelon-parameterised solution,
  # which is not interpretable directly -- use std_rotated_loadings instead.
  if (!is.null(x$std_rotated_loadings)) {
    lmat <- round(x$std_rotated_loadings, digits)
    if (suppress > 0) lmat[abs(lmat) < suppress] <- NA_real_
    return(lmat)
  }

  ss <- lavaan::standardizedsolution(x$lavaan_fit)
  # Filter to factor loading rows (op == "=~")
  ss_lv <- ss[ss$op == "=~", , drop = FALSE]

  # Pivot to matrix
  factors <- unique(ss_lv$lhs)
  items   <- unique(ss_lv$rhs)

  lmat <- matrix(NA_real_, nrow = length(items), ncol = length(factors),
                 dimnames = list(items, factors))

  for (i in seq_len(nrow(ss_lv))) {
    f <- ss_lv$lhs[i]
    v <- ss_lv$rhs[i]
    lmat[v, f] <- ss_lv$est.std[i]
  }

  lmat <- round(lmat, digits)

  if (suppress > 0)
    lmat[abs(lmat) < suppress] <- NA_real_

  lmat
}
