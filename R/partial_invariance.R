# -- Partial Measurement Invariance --------------------------------------------
#
# NOTE on parTable() lhs format for ESEM EFA blocks (verified Task 1):
# lhs uses BARE factor names (e.g. "EX"), NOT the compound 'efa("esem")*EX' form.
# The sub() prefix-strip below is a safe no-op but kept for defensive robustness.


.parse_score_labels <- function(lav_fit, uni_df) {

  pt <- lavaan::parTable(lav_fit)

  # Defensive: locate the score column by name (lavaan version differences)
  # Verified Task 1: column is "X2" in current lavaan version
  score_col <- intersect(c("X2", "x2", "Chisq", "chisq"), names(uni_df))[1L]
  if (is.na(score_col))
    stop("Cannot find score statistic column in lavTestScore()$uni. ",
         "Columns present: ", paste(names(uni_df), collapse=", "), call. = FALSE)

  grp_levels <- lavaan::lavInspect(lav_fit, "group.label")

  results <- lapply(seq_len(nrow(uni_df)), function(i) {
    plabel <- uni_df$lhs[i]   # e.g. ".p14."
    score  <- uni_df[[score_col]][i]

    # Skip rows where lavaan could not assign a plabel (e.g. factor mean
    # constraints appear as NA/NA in some lavaan versions)
    if (is.na(plabel)) return(NULL)

    # Find the actual parameter row (not the == row)
    param_row <- pt[!is.na(pt$plabel) & pt$plabel == plabel &
                    pt$op != "==" & pt$op != ":=", ]

    if (nrow(param_row) == 0L) return(NULL)

    # Prefer non-reference group row (group > 1) if available
    non_ref <- param_row[param_row$group > 1L, ]
    pr <- if (nrow(non_ref) > 0L) non_ref[1L, ] else param_row[1L, ]

    lhs   <- pr$lhs
    op    <- pr$op
    rhs   <- pr$rhs
    grp   <- pr$group

    # Strip efa("...") prefix from lhs if present (safe no-op for bare names)
    lhs_clean <- sub('^efa\\("[^"]*"\\)\\*', '', lhs)

    # Construct group.partial label based on operator
    gp_label <- switch(op,
      "=~" = paste0(lhs_clean, "=~", rhs),
      "~1" = paste0(lhs_clean, "~1"),
      "|"  = paste0(lhs_clean, "|",  rhs),
      paste0(lhs_clean, op, rhs)   # fallback
    )

    # Human-readable label
    human_label <- switch(op,
      "=~" = paste0(lhs_clean, " =~ ", rhs),
      "~1" = paste0(lhs_clean, " (intercept)"),
      "|"  = paste0(lhs_clean, " | ", rhs, " (threshold)"),
      paste0(lhs_clean, " ", op, " ", rhs)
    )

    # Group name (grp_levels computed once above the lapply)
    grp_name   <- if (!is.null(grp_levels) && !is.na(grp) && grp >= 1L &&
                      length(grp_levels) >= grp)
                    as.character(grp_levels[grp]) else paste0("G", grp)

    data.frame(
      label               = human_label,
      group_partial_label = gp_label,
      group_name          = grp_name,
      score               = score,
      stringsAsFactors    = FALSE
    )
  })

  results <- do.call(rbind, Filter(Negate(is.null), results))
  if (is.null(results) || nrow(results) == 0L)
    return(data.frame(label=character(), group_partial_label=character(),
                      group_name=character(), score=numeric()))

  results[order(-results$score), ]
}


#' Partial Measurement Invariance Testing
#'
#' Given a failed invariance level, uses a greedy score-test loop
#' (\code{lavTestScore} + \code{group.partial}) to identify and free the
#' minimum set of non-invariant equality constraints needed to restore
#' acceptable fit (DCFI >= -0.010). Downstream invariance levels are then
#' re-tested under the same partial constraints.
#'
#' @param inv An \code{esem_invariance} object from \code{\link{esem_invariance}}.
#' @param level Character. The invariance level to partially free:
#'   \code{"weak"}, \code{"strong"}, or \code{"strict"}.
#' @param max_free Integer. Maximum parameters to free before stopping.
#'   Default 10. Byrne et al. (1989) recommend freeing the minimum number
#'   -- in practice >5 rarely recovers invariance.
#' @param delta_cfi_cutoff Numeric. Stopping criterion: loop stops when
#'   DCFI >= this value. Default -0.010 (Cheung & Rensvold, 2002).
#' @param verbose Logical. Print progress messages. Default \code{TRUE}.
#'
#' @return An object of class \code{"esem_partial_invariance"} with:
#' \describe{
#'   \item{\code{$freed_params}}{Data frame: round, label, group_name,
#'     score (LM), delta_cfi (baseline-relative), converged.}
#'   \item{\code{$partial_fit}}{esem_fit at \code{level} with freed params.}
#'   \item{\code{$downstream}}{Named list of esem_fit for levels above target.}
#'   \item{\code{$downstream_lrt}}{lavTestLRT() comparisons among downstream.}
#'   \item{\code{$table}}{Combined fit table (original + partial + downstream).}
#'   \item{\code{$converged}}{Logical: did DCFI pass cutoff?}
#'   \item{\code{$group_partial}}{Character vector: final group.partial labels.}
#' }
#'
#' @section Limitations:
#' B-ESEM models are not supported. Use
#' \code{lavaan::lavTestScore(inv$models[[level]]$lavaan_fit)} directly.
#'
#' @seealso \code{\link{esem_invariance}}
#' @export
partial_invariance <- function(inv,
                               level,
                               max_free         = 10L,
                               delta_cfi_cutoff = -0.010,
                               verbose          = TRUE) {

  # -- 1. Validate inputs ------------------------------------------------------
  if (!inherits(inv, "esem_invariance"))
    stop("`inv` must be an esem_invariance object from esem_invariance().",
         call. = FALSE)

  if (!is.null(inv$model) && inv$model == "besem")
    stop("partial_invariance() is not supported for B-ESEM models. ",
         "Use lavaan::lavTestScore() on inv$models[[level]]$lavaan_fit directly.",
         call. = FALSE)

  level <- match.arg(level, c("weak", "strong", "strict"))

  if (is.null(inv$models[[level]]))
    stop(sprintf("inv$models$%s is NULL -- the original %s model failed to fit.",
                 level, level), call. = FALSE)

  if (!.inv_converged(inv$models[[level]]))
    warning(sprintf("inv$models$%s did not converge. Results may be unreliable.",
                    level), call. = FALSE)

  # -- 2. Identify baseline (less-constrained) model --------------------------
  baseline_map   <- c(weak = "configural", strong = "weak", strict = "strong")
  baseline_level <- baseline_map[[level]]
  baseline_fit   <- inv$models[[baseline_level]]

  if (is.null(baseline_fit))
    stop(sprintf("Baseline model (inv$models$%s) is NULL.", baseline_level),
         call. = FALSE)

  is_ordered <- !is.null(inv$spec$ordered) && length(inv$spec$ordered) > 0

  # Resolve missing -- match esem_invariance() convention (pairwise for WLSMV)
  missing_arg <- if (!is.null(inv$spec$missing)) {
    inv$spec$missing
  } else if (is_ordered) {
    "pairwise"
  } else {
    "listwise"
  }

  # Reconstruct constraint list (same logic as esem_invariance)
  constraints <- if (is_ordered) {
    list(
      weak   = "loadings",
      strong = c("loadings", "thresholds"),
      strict = c("loadings", "thresholds", "residuals")
    )
  } else {
    list(
      weak   = "loadings",
      strong = c("loadings", "intercepts"),
      strict = c("loadings", "intercepts", "residuals")
    )
  }

  baseline_cfi <- .extract_fit_inv(baseline_fit)[["cfi"]]

  # -- 3. Initialise greedy loop -----------------------------------------------
  current_fit   <- inv$models[[level]]
  group_partial <- character(0)
  freed_log     <- list()
  converged     <- FALSE

  if (verbose) {
    cat("======================================================\n")
    cat(sprintf(" Partial Invariance: %s -- level: %s\n", inv$spec$label, level))
    cat(sprintf(" Baseline: %s | dCFI cutoff: %.3f | max_free: %d\n",
                baseline_level, delta_cfi_cutoff, max_free))
    cat("======================================================\n\n")
  }

  # -- 4. Greedy loop ---------------------------------------------------------
  for (round_i in seq_len(max_free)) {

    if (verbose) cat(sprintf("  Round %d: running lavTestScore ...", round_i))

    sc <- tryCatch(
      lavaan::lavTestScore(current_fit$lavaan_fit),
      error = function(e) {
        warning(sprintf("lavTestScore() failed in round %d: %s",
                        round_i, conditionMessage(e)), call. = FALSE)
        NULL
      }
    )

    if (is.null(sc) || is.null(sc$uni) || nrow(sc$uni) == 0L) {
      if (verbose) cat(" FAILED\n")
      break
    }

    # Parse labels; exclude already-freed parameters
    lbl_df <- .parse_score_labels(current_fit$lavaan_fit, sc$uni)
    lbl_df <- lbl_df[!lbl_df$group_partial_label %in% group_partial, ]

    if (nrow(lbl_df) == 0L) {
      if (verbose) cat(" no more constraints to free\n")
      break
    }

    best <- lbl_df[1L, ]  # already sorted descending by score
    group_partial <- c(group_partial, best$group_partial_label)

    if (verbose) cat(sprintf(" freeing %s (score=%.2f) ...",
                             best$label, best$score))

    # Refit with updated group.partial
    new_fit <- tryCatch(
      .fit_invariance_model(
        spec          = inv$spec,
        group_equal   = constraints[[level]],
        is_ordered    = is_ordered,
        model         = "esem",
        missing       = missing_arg,
        verbose       = FALSE,
        group.partial = group_partial
      ),
      error = function(e) {
        warning(sprintf("Refit failed in round %d: %s",
                        round_i, conditionMessage(e)), call. = FALSE)
        NULL
      }
    )

    if (is.null(new_fit)) {
      if (verbose) cat(" FAILED\n")
      group_partial <- group_partial[-length(group_partial)]   # undo last addition
      break
    }

    new_cfi  <- .extract_fit_inv(new_fit)[["cfi"]]
    dcfi     <- round(new_cfi, 3) - round(baseline_cfi, 3)

    # Store dcfi as pre-rounded (difference of two individually-rounded values,
    # matching the display rounding in the table). Do NOT call round(dcfi, 3)
    # again -- that would double-round and could make freed_params$delta_cfi
    # disagree with the table by +/-0.001 at boundary values.
    freed_log[[round_i]] <- data.frame(
      round      = round_i,
      label      = best$label,
      group_name = best$group_name,
      score      = round(best$score, 2),
      delta_cfi  = dcfi,
      converged  = dcfi >= delta_cfi_cutoff,
      stringsAsFactors = FALSE
    )

    if (verbose) cat(sprintf(" dCFI = %+.3f\n", dcfi))

    if (dcfi >= delta_cfi_cutoff) {
      converged   <- TRUE
      current_fit <- new_fit
      break
    }

    current_fit <- new_fit
  }

  if (!converged && verbose)
    warning(sprintf(
      "partial_invariance: dCFI cutoff (%.3f) not reached after freeing %d parameter(s).",
      delta_cfi_cutoff, length(group_partial)
    ), call. = FALSE)

  partial_fit  <- current_fit
  freed_params <- if (length(freed_log) > 0L)
    do.call(rbind, freed_log)
  else
    data.frame(round=integer(), label=character(), group_name=character(),
               score=numeric(), delta_cfi=numeric(), converged=logical())

  # -- 5. Fit downstream levels ------------------------------------------------
  downstream_map <- list(
    weak   = c("strong", "strict"),
    strong = c("strict"),
    strict = character(0)
  )
  downstream_levels <- downstream_map[[level]]

  downstream     <- list()
  downstream_lrt <- list()
  # prev_ds_fit starts as the partial fit at the target level (not a downstream
  # fit), so the first LRT compares partial-target vs first-downstream. Subsequent
  # iterations compare adjacent downstream levels.
  prev_ds_fit    <- partial_fit

  for (ds_lv in downstream_levels) {
    if (verbose) cat(sprintf("  Fitting downstream: %s (partial) ...", ds_lv))

    ds_fit <- tryCatch(
      .fit_invariance_model(
        spec          = inv$spec,
        group_equal   = constraints[[ds_lv]],
        is_ordered    = is_ordered,
        model         = "esem",
        missing       = missing_arg,
        verbose       = FALSE,
        group.partial = group_partial
      ),
      error = function(e) {
        warning(sprintf("Downstream %s fit failed: %s", ds_lv,
                        conditionMessage(e)), call. = FALSE)
        NULL
      }
    )

    conv_msg <- if (is.null(ds_fit)) " FAILED\n" else " OK\n"
    if (verbose) cat(conv_msg)

    downstream[[ds_lv]] <- ds_fit

    if (!is.null(ds_fit) && !is.null(prev_ds_fit)) {
      downstream_lrt[[ds_lv]] <- tryCatch(
        lavaan::lavTestLRT(prev_ds_fit$lavaan_fit, ds_fit$lavaan_fit),
        error = function(e) NULL
      )
    }

    prev_ds_fit <- ds_fit
  }

  # -- 6. Build combined table -------------------------------------------------
  table_out <- .build_partial_table(
    inv                  = inv,
    level                = level,
    partial_fit          = partial_fit,
    downstream           = downstream,
    downstream_lrt       = downstream_lrt,
    group_partial_labels = group_partial,
    is_ordered           = is_ordered,
    delta_cfi_cutoff     = delta_cfi_cutoff
  )

  if (verbose) {
    cat("\n")
    .print_partial_table(table_out, is_ordered, freed_params)
  }

  structure(
    list(
      freed_params   = freed_params,
      partial_fit    = partial_fit,
      downstream     = downstream,
      downstream_lrt = downstream_lrt,
      table          = table_out,
      spec           = inv$spec,
      model          = inv$model,
      level          = level,
      converged      = converged,
      group_partial  = group_partial
    ),
    class = "esem_partial_invariance"
  )
}


.build_partial_table <- function(inv, level, partial_fit, downstream,
                                 downstream_lrt, group_partial_labels,
                                 is_ordered, delta_cfi_cutoff) {

  level_order  <- c("configural", "weak", "strong", "strict")
  target_idx   <- which(level_order == level)
  baseline_lv  <- c(weak="configural",strong="weak",strict="strong")[[level]]
  baseline_idx <- which(level_order == baseline_lv)

  # Level display labels
  level_labels <- c(
    configural = "1. Configural",
    weak       = "2. Weak (metric)",
    strong     = "3. Strong (scalar)",
    strict     = "4. Strict"
  )
  partial_prefix <- c(weak="2p",strong="3p",strict="4p")[[level]]
  # NOTE: must use a nested list, not c() -- c() flattens and mangles names.
  # Pre-subset by level so ds_prefix is a simple named vector for the loop.
  ds_prefix <- list(
    weak   = c(strong = "3p", strict = "4p"),
    strong = c(strict = "4p"),
    strict = character(0)
  )[[level]]

  # -- Original rows up to and including target level --------------------------
  orig_rows <- inv$table[seq_len(target_idx), , drop = FALSE]

  # Relabel target row: "scalar" -> "full"
  orig_rows$Model[target_idx] <- gsub("scalar", "full",
                                      orig_rows$Model[target_idx])
  orig_rows$Model[target_idx] <- gsub("\\(metric\\)", "(metric, full)",
                                      orig_rows$Model[target_idx])
  # Mark with x via a new column
  orig_rows$pass <- NA_character_
  orig_rows$pass[target_idx] <- "x"   # x

  # -- Partial row -------------------------------------------------------------
  # DCFI for partial row is relative to baseline (same as greedy loop)
  baseline_cfi  <- round(inv$table$CFI[baseline_idx], 3)
  partial_fi    <- .extract_fit_inv(partial_fit)
  partial_cfi_r <- round(partial_fi[["cfi"]], 3)
  partial_dcfi  <- partial_cfi_r - baseline_cfi

  partial_row <- data.frame(
    Model       = paste0(partial_prefix, ". ", sub("^[0-9]+\\. ", "",
                   level_labels[[level]]), " (partial)"),
    chisq       = unname(partial_fi["chisq"]),
    df          = unname(partial_fi["df"]),
    delta_chisq = NA_real_,
    delta_df    = NA_real_,
    pvalue      = NA_real_,
    CFI         = unname(partial_fi["cfi"]),
    dCFI        = partial_dcfi,
    TLI         = unname(partial_fi["tli"]),
    dTLI        = NA_real_,
    RMSEA       = unname(partial_fi["rmsea"]),
    dRMSEA      = NA_real_,
    SRMR        = unname(partial_fi["srmr"]),
    dSMR        = NA_real_,
    pass        = "v",   # v
    row.names   = NULL,
    stringsAsFactors = FALSE
  )

  # -- Downstream rows ----------------------------------------------------------
  ds_rows    <- list()
  prev_cfi_r <- partial_cfi_r

  for (ds_lv in names(downstream)) {
    # Guard: downstream refit may have failed (stored as NULL)
    if (is.null(downstream[[ds_lv]])) next

    ds_fi    <- .extract_fit_inv(downstream[[ds_lv]])
    ds_cfi_r <- round(ds_fi[["cfi"]], 3)
    ds_dcfi  <- ds_cfi_r - prev_cfi_r
    lrt_row  <- .extract_lrt_row(downstream_lrt[[ds_lv]])
    pfx      <- if (!is.null(ds_prefix[[ds_lv]]) && !is.na(ds_prefix[[ds_lv]]))
                  ds_prefix[[ds_lv]] else ds_lv

    ds_rows[[ds_lv]] <- data.frame(
      Model       = paste0(pfx, ". ",
                     sub("^[0-9]+\\. ", "", level_labels[[ds_lv]]),
                     " (partial)"),
      chisq       = unname(ds_fi["chisq"]),
      df          = unname(ds_fi["df"]),
      delta_chisq = unname(lrt_row["delta_chisq"]),
      delta_df    = unname(lrt_row["delta_df"]),
      pvalue      = unname(lrt_row["pvalue"]),
      CFI         = unname(ds_fi["cfi"]),
      dCFI        = ds_dcfi,
      TLI         = unname(ds_fi["tli"]),
      dTLI        = NA_real_,
      RMSEA       = unname(ds_fi["rmsea"]),
      dRMSEA      = NA_real_,
      SRMR        = unname(ds_fi["srmr"]),
      dSMR        = NA_real_,
      pass        = NA_character_,
      row.names   = NULL,
      stringsAsFactors = FALSE
    )
    prev_cfi_r <- ds_cfi_r
  }

  all_rows <- c(list(orig_rows), list(partial_row), ds_rows)
  do.call(rbind, all_rows)
}


.print_partial_table <- function(tbl, is_ordered, freed_params) {

  .fmt <- function(x, d = 3) {
    s <- ifelse(is.na(x), "", formatC(round(x, d), format = "f", digits = d))
    # formatC("-0.000") is confusing; collapse to "0.000"
    zero_pat <- paste0("^-0\\.", paste(rep("0", d), collapse = ""), "$")
    gsub(zero_pat, paste0("0.", paste(rep("0", d), collapse = "")), s)
  }
  .fmt_int <- function(x) ifelse(is.na(x), "", as.character(round(x, 0)))

  # Build display data frame
  display <- data.frame(
    "Model" = paste0(tbl$Model, ifelse(is.na(tbl$pass), "", paste0(" ", tbl$pass))),
    "chi2"  = .fmt(tbl$chisq),
    "df"    = .fmt_int(tbl$df),
    "CFI"   = .fmt(tbl$CFI),
    "dCFI"  = .fmt(tbl$dCFI),
    "RMSEA" = .fmt(tbl$RMSEA),
    "SRMR"  = .fmt(tbl$SRMR),
    stringsAsFactors = FALSE,
    check.names      = FALSE
  )
  names(display)[4:5] <- c("CFI", "dCFI")
  print(display, row.names = FALSE)

  # Freed parameters table
  if (nrow(freed_params) > 0L) {
    cat("\nFreed parameters (greedy):\n")
    cat(sprintf("  %-5s  %-30s  %-8s  %9s  %12s\n",
                "Round", "Parameter", "Group", "Score(LM)", "dCFI (vs baseline)"))
    cat(sprintf("  %-5s  %-30s  %-8s  %9s  %12s\n",
                "-----", "------------------------------",
                "--------", "---------", "------------"))
    for (i in seq_len(nrow(freed_params))) {
      r    <- freed_params[i, ]
      tick <- if (isTRUE(r$converged)) " v" else ""
      cat(sprintf("  %-5d  %-30s  %-8s  %9.2f  %+12.3f%s\n",
                  r$round, r$label, r$group_name, r$score, r$delta_cfi, tick))
    }
  }

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

  cat("\nAccess results:\n")
  cat("  pinv$partial_fit    -- esem_fit (partial)\n")
  cat("  pinv$downstream     -- list of esem_fit (downstream partial levels)\n")
  cat("  pinv$freed_params   -- data frame of freed parameters\n")
  cat("  pinv$table          -- combined fit table\n\n")

  invisible(NULL)
}


#' Print Method for esem_partial_invariance
#'
#' @param x An \code{esem_partial_invariance} object.
#' @param ... Ignored.
#' @return Invisibly returns \code{x}; called for the side effect of printing
#'   the partial-invariance summary.
#' @export
print.esem_partial_invariance <- function(x, ...) {
  is_ordered  <- !is.null(x$spec$ordered) && length(x$spec$ordered) > 0
  n_freed     <- nrow(x$freed_params)
  status      <- if (x$converged) "converged" else "did not reach cutoff"

  cat("======================================================\n")
  cat(sprintf(" ESEM Partial Measurement Invariance: %s\n", x$spec$label))
  cat(sprintf(" Level tested: %s  |  %d parameter(s) freed  |  %s\n",
              x$level, n_freed, status))
  cat("======================================================\n\n")

  .print_partial_table(x$table, is_ordered, x$freed_params)

  invisible(x)
}
