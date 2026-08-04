#' Compare R and Mplus Standardised Loadings
#'
#' Builds a long-format data frame with STDYX loadings, standard errors,
#' z-scores, p-values, and R \eqn{-} Mplus differences for every loading
#' (primary and cross) across CFA, ESEM, and B-ESEM.  Mplus columns are
#' \code{NA} when no Mplus results are present in the pipeline object.
#'
#' @param results An \code{esem_comparison_pipeline} object from
#'   \code{\link{run_comparison}}.
#'
#' @return A data frame with columns:
#' \describe{
#'   \item{\code{model}}{CFA, ESEM, or BESEM.}
#'   \item{\code{factor}}{Factor name (lowercase).}
#'   \item{\code{item}}{Item name (lowercase).}
#'   \item{\code{loading_type}}{\code{"primary"} or \code{"cross"}.}
#'   \item{\code{R_std}, \code{R_se}, \code{R_z}, \code{R_p}}{R estimates.}
#'   \item{\code{Mplus_std}, \code{Mplus_se}, \code{Mplus_z}, \code{Mplus_p}}{Mplus estimates (\code{NA} if unavailable).}
#'   \item{\code{diff_std}}{R \eqn{-} Mplus standardised loading difference (\code{NA} if unavailable).}
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
#' results <- run_comparison(spec, n_starts = 5L)
#' lc <- compare_loadings(results)
#' head(lc)
#'
#' # Primary loadings only (Mplus_* columns are NA without a Mplus run)
#' lc[lc$loading_type == "primary", ]
#' }
#'
#' @export
compare_loadings <- function(results) {

  if (!inherits(results, "esem_comparison_pipeline"))
    stop("`results` must be an esem_comparison_pipeline object from run_comparison().",
         call. = FALSE)

  spec         <- results$spec
  mplus_res    <- results$mplus_results

  r_cfa   <- .r_loadings_long(results$fit_cfa,   "CFA")
  r_esem  <- .r_loadings_long(results$fit_esem,  "ESEM")
  r_besem <- .r_loadings_long(results$fit_besem, "BESEM")

  mp_cfa   <- .mplus_loadings_long(if (!is.null(mplus_res)) mplus_res$cfa   else NULL, "CFA")
  mp_esem  <- .mplus_loadings_long(if (!is.null(mplus_res)) mplus_res$esem  else NULL, "ESEM")
  mp_besem <- .mplus_loadings_long(if (!is.null(mplus_res)) mplus_res$besem else NULL, "BESEM")

  rbind(
    .merge_loading_tables(r_cfa,   mp_cfa,   spec),
    .merge_loading_tables(r_esem,  mp_esem,  spec),
    .merge_loading_tables(r_besem, mp_besem, spec)
  )
}


#' Save Pipeline Results to CSV or xlsx
#'
#' Writes all standard output files to \code{output_folder}:
#' \itemize{
#'   \item \code{<label>_fit_indices.csv} -- CFI/TLI/RMSEA/SRMR for all models.
#'   \item \code{<label>_CFA_loadings.csv} -- standardised loading matrix (wide).
#'   \item \code{<label>_ESEM_loadings.csv} -- standardised loading matrix (wide).
#'   \item \code{<label>_BESEM_loadings.csv} -- standardised loading matrix (wide).
#'   \item \code{<label>_loadings_comparison.csv} -- R vs Mplus loadings, long
#'     format with SEs, z-scores, p-values, and \code{diff_std}.
#'   \item \code{<label>_omega.csv} -- Reliability indices (if
#'     \code{indices} is supplied).
#'   \item \code{<label>_results.xlsx} -- Single xlsx workbook with all results
#'     (when \code{xlsx = TRUE}).
#' }
#'
#' @param results An \code{esem_comparison_pipeline} object from
#'   \code{\link{run_comparison}}.
#' @param omega Optional. A \code{reliability_indices} object from
#'   \code{\link{compute_indices}}. If supplied, reliability indices are saved.
#' @param output_folder Character. Path to the folder where files are written.
#'   Created if it does not exist.
#' @param label Character. Prefix for output file names. Defaults to
#'   \code{results$spec$label}.
#' @param xlsx Logical. If \code{TRUE}, writes a single
#'   \code{<label>_results.xlsx} file instead of individual CSVs.  Requires
#'   the \code{openxlsx2} package (\code{install.packages("openxlsx2")}).
#'   Sheets: \code{Fit_Indices}, \code{CFA_Loadings},
#'   \code{ESEM_Loadings} (if ESEM fit present),
#'   \code{BESEM_Loadings} (if B-ESEM fit present),
#'   \code{Reliability} (if \code{omega} supplied),
#'   \code{Loadings_Comparison} (if Mplus results present and both
#'   ESEM and B-ESEM fits available).
#'   Primary loadings in ESEM/BESEM sheets are bolded; G-column bolded in
#'   BESEM sheet.  Default \code{FALSE}.
#'
#' @return Invisibly returns a character vector of file paths written.
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
#' results <- run_comparison(spec, n_starts = 5L)
#' indices <- compute_indices(results)
#' save_results(results, indices,
#'              output_folder = file.path(tempdir(), "esem_results"))
#' }
#'
#' @export
save_results <- function(results, omega = NULL,
                         output_folder,
                         label = NULL,
                         xlsx = FALSE) {

  if (!inherits(results, "esem_comparison_pipeline"))
    stop("`results` must be an esem_comparison_pipeline object from run_comparison().",
         call. = FALSE)

  if (isTRUE(xlsx) && !requireNamespace("openxlsx2", quietly = TRUE))
    stop("'xlsx = TRUE' requires the openxlsx2 package. Install it with:\n",
         "  install.packages(\"openxlsx2\")", call. = FALSE)

  if (is.null(label)) label <- results$spec$label
  dir.create(output_folder, showWarnings = FALSE, recursive = TRUE)

  written <- character(0)

  fp <- function(suffix) file.path(output_folder, paste0(label, suffix))

  lc <- NULL

  if (!isTRUE(xlsx)) {
    # Fit indices
    write.csv(results$comparison_table, fp("_fit_indices.csv"), row.names = FALSE)
    written <- c(written, fp("_fit_indices.csv"))
    if (isTRUE(attr(results$comparison_table, "wlsmv_note"))) {
      cat(
        paste0(
          "* WLSMV delta-chi2 computed via lavaan scaled difference test;",
          " may not match Mplus DIFFTEST.",
          " Mplus delta-p omitted for WLSMV (statistic is not chi2-distributed).",
          " B-ESEM vs ESEM delta-chi2 suppressed on WLSMV path",
          " (auxiliary lavaan fit not comparable).\n"
        ),
        file = fp("_fit_indices.csv"), append = TRUE
      )
    }

    # Wide loading matrices
    cfa_ss   <- lavaan::standardizedsolution(results$fit_cfa)
    cfa_lmat <- with(cfa_ss[cfa_ss$op == "=~", ],
                     tapply(est.std, list(rhs, lhs), identity))
    write.csv(cfa_lmat, fp("_CFA_loadings.csv"))
    written <- c(written, fp("_CFA_loadings.csv"))

    write.csv(std_loadings(results$fit_esem),  fp("_ESEM_loadings.csv"))
    written <- c(written, fp("_ESEM_loadings.csv"))

    write.csv(std_loadings(results$fit_besem), fp("_BESEM_loadings.csv"))
    written <- c(written, fp("_BESEM_loadings.csv"))

    # Loading comparison (R vs Mplus)
    lc <- compare_loadings(results)
    write.csv(lc, fp("_loadings_comparison.csv"), row.names = FALSE)
    written <- c(written, fp("_loadings_comparison.csv"))

    # Omega -- .build_omega_table() returns NULL if all slots are non-model (e.g. alpha-only)
    if (!is.null(omega)) {
      omega_tbl <- .build_omega_table(omega)
      if (!is.null(omega_tbl)) {
        write.csv(.pivot_omega_wide(omega_tbl), fp("_omega.csv"), row.names = FALSE, na = "-")
        written <- c(written, fp("_omega.csv"))
      }
    }
  } else {
    # -- xlsx path -------------------------------------------------------------
    wb_data   <- list()
    bold_info <- list()

    # Fit indices
    wb_data[["Fit_Indices"]] <- results$comparison_table

    # CFA loadings -- sorted by spec factor order; NAs blanked via na_strings
    cfa_ss   <- lavaan::standardizedsolution(results$fit_cfa)
    cfa_lmat <- with(cfa_ss[cfa_ss$op == "=~", ],
                     tapply(est.std, list(rhs, lhs), mean))
    spec_items_lower <- tolower(unlist(results$spec$factors, use.names = FALSE))
    mat_rows_lower   <- tolower(rownames(cfa_lmat))
    order_idx        <- match(spec_items_lower, mat_rows_lower)
    order_idx        <- order_idx[!is.na(order_idx)]
    if (length(order_idx) > 0L) cfa_lmat <- cfa_lmat[order_idx, , drop = FALSE]
    cfa_df   <- cbind(Item = rownames(as.data.frame(cfa_lmat)),
                      as.data.frame(cfa_lmat))
    wb_data[["CFA_Loadings"]] <- cfa_df

    # ESEM loadings
    if (!is.null(results$fit_esem)) {
      esem_mat <- std_loadings(results$fit_esem)
      wb_data[["ESEM_Loadings"]]   <- cbind(Item = rownames(esem_mat),
                                             as.data.frame(esem_mat))
      bold_info[["ESEM_Loadings"]] <- list(mat      = esem_mat,
                                            spec     = results$spec,
                                            fit_obj  = results$fit_esem,
                                            is_besem = FALSE)
    }

    # BESEM loadings
    if (!is.null(results$fit_besem)) {
      besem_mat <- std_loadings(results$fit_besem)
      wb_data[["BESEM_Loadings"]]   <- cbind(Item = rownames(besem_mat),
                                              as.data.frame(besem_mat))
      bold_info[["BESEM_Loadings"]] <- list(mat      = besem_mat,
                                             spec     = results$spec,
                                             fit_obj  = results$fit_besem,
                                             is_besem = TRUE)
    }

    # Reliability
    if (!is.null(omega)) {
      omega_tbl <- .build_omega_table(omega)
      if (!is.null(omega_tbl))
        wb_data[["Reliability"]] <- .pivot_omega_wide(omega_tbl)
    }

    # Loadings comparison -- only when Mplus results AND both fits available
    # (compare_loadings() has no NULL guard for fit_esem/fit_besem internally)
    if (!is.null(results$mplus_results) &&
        !is.null(results$fit_esem) &&
        !is.null(results$fit_besem)) {
      lc <- compare_loadings(results)
      wb_data[["Loadings_Comparison"]] <- lc
    }

    out_path <- fp("_results.xlsx")
    .write_xlsx(wb_data, bold_info, out_path)
    written <- c(written, out_path)
  }

  # Print summary
  message(paste("Results saved to:", output_folder))
  for (f in written) message(paste(" ", basename(f)))

  # Max |R - Mplus| summary if lc was populated
  if (!is.null(lc)) {
    safe_max <- function(x) { x <- abs(x[!is.na(x)]); if (length(x) == 0) NA_real_ else max(x) }
    message("\nMax |R - Mplus| STDYX:")
    for (mod in c("CFA", "ESEM", "BESEM")) {
      sl <- lc[lc$model == mod, ]
      mp <- safe_max(sl$diff_std[sl$loading_type == "primary"])
      mc <- safe_max(sl$diff_std[sl$loading_type == "cross"])
      message(sprintf("  %-6s  primary: %s   cross: %s",
                  mod,
                  if (is.na(mp)) "   --  " else sprintf("%.4f", mp),
                  if (is.na(mc)) "   --  " else sprintf("%.4f", mc)))
    }
    bl <- lc[lc$model == "BESEM", ]
    if (!all(is.na(bl$diff_std)) && safe_max(bl$diff_std) > 0.10)
      message(paste0(
          "\n  Note: BESEM loadings differ > 0.10 between R and Mplus.\n",
          "  This is normal when both converge to different rotation solutions.\n",
          "  Fit indices (CFI/TLI/RMSEA) are rotation-invariant and should match."))
  }

  invisible(written)
}


# -- Internal helpers ----------------------------------------------------------

.build_omega_table <- function(omega) {
  omega_rows <- lapply(omega, function(r) {
    if (is.null(r)) return(NULL)
    if (!is.list(r) || is.null(r$model)) return(NULL)
    rows <- list(
      data.frame(model = r$model, level = "total", factor = "--",
                 statistic = "omega_total", value = r$omega_total,
                 rotation_invariant = "yes", stringsAsFactors = FALSE)
    )
    if (!is.na(r$omega_h_g))
      rows <- c(rows, list(
        data.frame(model = r$model, level = "total", factor = "G",
                   statistic = "omega_H_G", value = r$omega_h_g,
                   rotation_invariant = "no", stringsAsFactors = FALSE),
        data.frame(model = r$model, level = "total", factor = "G",
                   statistic = "ECV", value = r$ecv,
                   rotation_invariant = "no", stringsAsFactors = FALSE)
      ))
    if (!is.null(r$H))
      rows <- c(rows, lapply(names(r$H), function(f)
        data.frame(model = r$model, level = "total", factor = f,
                   statistic = "H", value = r$H[[f]],
                   rotation_invariant = "no", stringsAsFactors = FALSE)))
    if (!is.null(r$subscales))
      rows <- c(rows, lapply(r$subscales, function(s) {
        if (is.null(s)) return(NULL)
        do.call(rbind, list(
          data.frame(model = r$model, level = "subscale", factor = s$factor,
                     statistic = "omega_S",   value = s$omega_specific,
                     rotation_invariant = "no", stringsAsFactors = FALSE),
          data.frame(model = r$model, level = "subscale", factor = s$factor,
                     statistic = "omega_sub", value = s$omega_sub,
                     rotation_invariant = "no", stringsAsFactors = FALSE),
          data.frame(model = r$model, level = "subscale", factor = s$factor,
                     statistic = "omega_H_G", value = s$omega_h_g,
                     rotation_invariant = "no", stringsAsFactors = FALSE),
          data.frame(model = r$model, level = "subscale", factor = s$factor,
                     statistic = "ECV_s",     value = s$ecv_sub,
                     rotation_invariant = "no", stringsAsFactors = FALSE),
          data.frame(model = r$model, level = "subscale", factor = s$factor,
                     statistic = "H_Fs",      value = s$H_specific,
                     rotation_invariant = "no", stringsAsFactors = FALSE)
        ))
      }))
    do.call(rbind, rows[!sapply(rows, is.null)])
  })
  tbl <- do.call(rbind, omega_rows[!sapply(omega_rows, is.null)])
  if (is.null(tbl) || nrow(tbl) == 0L) return(NULL)
  tbl
}

.pivot_omega_wide <- function(omega_tbl) {
  # Human-readable stat labels
  stat_labels <- c(
    omega_total = "omega_total",
    omega_H_G   = "omega_H_G",
    ECV         = "ECV",
    H           = "H",
    omega_S     = "omega_specific",
    omega_sub   = "omega_subscale",
    ECV_s       = "ECV_s",
    H_Fs        = "H_specific"
  )
  omega_tbl$statistic <- ifelse(
    omega_tbl$statistic %in% names(stat_labels),
    stat_labels[omega_tbl$statistic],
    omega_tbl$statistic
  )

  omega_tbl$factor[omega_tbl$factor == "--"] <- "Total"

  all_factors <- unique(omega_tbl$factor)
  factor_cols <- c(
    intersect("Total", all_factors),
    intersect("G",     all_factors),
    setdiff(all_factors, c("Total", "G"))
  )

  stat_order <- c("omega_total", "omega_subscale",
                  "omega_H_G", "omega_specific",
                  "ECV", "ECV_s", "H", "H_specific")
  model_order <- c("CFA", "ESEM", "BESEM")
  models <- c(intersect(model_order, omega_tbl$model),
              setdiff(omega_tbl$model, model_order))

  rows <- lapply(models, function(m) {
    sub       <- omega_tbl[omega_tbl$model == m, , drop = FALSE]
    all_stats <- unique(sub$statistic)
    stats     <- c(intersect(stat_order, all_stats), setdiff(all_stats, stat_order))

    lapply(stats, function(s) {
      sv  <- sub[sub$statistic == s, c("factor", "value"), drop = FALSE]
      row <- data.frame(Model = m, Statistic = s, stringsAsFactors = FALSE)
      for (fc in factor_cols) {
        idx       <- which(sv$factor == fc)
        row[[fc]] <- if (length(idx) == 1L) round(sv$value[idx], 3L) else NA_real_
      }
      row
    })
  })

  tbl <- do.call(rbind, unlist(rows, recursive = FALSE))
  rownames(tbl) <- NULL

  # Drop rows where every factor column is NA (bifactor-only metrics absent in CFA/ESEM)
  factor_cols_in_tbl <- intersect(factor_cols, colnames(tbl))
  all_na <- apply(tbl[, factor_cols_in_tbl, drop = FALSE], 1,
                  function(r) all(is.na(r)))
  tbl[!all_na, , drop = FALSE]
}

.style_sheet <- function(wb, sheet, nrow_data, ncol_data) {
  wb <- openxlsx2::wb_add_font(wb, sheet,
                                bold = TRUE,
                                dims = openxlsx2::wb_dims(rows = 1,
                                                           cols = seq_len(ncol_data)))
  wb <- openxlsx2::wb_freeze_pane(wb, sheet, first_row = TRUE)
  wb <- openxlsx2::wb_set_col_widths(wb, sheet,
                                      cols   = seq_len(ncol_data),
                                      widths = "auto")
  invisible(wb)
}

.bold_target_cells <- function(wb, sheet, mat, spec, fit_obj, is_besem) {
  primary_map <- setNames(
    tolower(rep(names(spec$factors), lengths(spec$factors))),
    tolower(unlist(spec$factors, use.names = FALSE))
  )
  mat_col_lower <- tolower(colnames(mat))

  for (i in seq_len(nrow(mat))) {
    pf <- primary_map[tolower(rownames(mat)[i])]
    if (is.na(pf)) next
    col_idx <- which(mat_col_lower == pf)
    if (length(col_idx) != 1L) next
    wb <- openxlsx2::wb_add_font(wb, sheet, bold = TRUE,
                                  dims = openxlsx2::wb_dims(rows = i + 1L,
                                                             cols = col_idx + 1L))
  }

  if (isTRUE(is_besem)) {
    g_name <- tolower(fit_obj$g_name %||% "G")
    g_col  <- which(mat_col_lower == g_name)
    if (length(g_col) == 1L) {
      for (i in seq_len(nrow(mat))) {
        wb <- openxlsx2::wb_add_font(wb, sheet, bold = TRUE,
                                      dims = openxlsx2::wb_dims(rows = i + 1L,
                                                                 cols = g_col + 1L))
      }
    }
  }

  invisible(wb)
}

.na_to_dash <- function(df) {
  for (j in seq_along(df)) {
    nas <- is.na(df[[j]])
    if (any(nas)) {
      df[[j]]       <- as.character(df[[j]])
      df[[j]][nas]  <- "-"
    }
  }
  df
}

.write_xlsx <- function(wb_data, bold_info = list(), path) {
  wb <- openxlsx2::wb_workbook()
  for (nm in names(wb_data)) {
    df <- .na_to_dash(wb_data[[nm]])
    wb <- openxlsx2::wb_add_worksheet(wb, nm)
    wb <- openxlsx2::wb_add_data(wb, nm, df)
    wb <- .style_sheet(wb, nm, nrow(df), ncol(df))
    if (!is.null(bold_info[[nm]])) {
      bi <- bold_info[[nm]]
      wb <- .bold_target_cells(wb, nm, bi$mat, bi$spec, bi$fit_obj, bi$is_besem)
    }
  }
  openxlsx2::wb_save(wb, path)
  invisible(path)
}

.r_loadings_long <- function(fit_obj, model_name) {
  # B-ESEM rotation -- esem_fit with rotation matrices
  if (inherits(fit_obj, "esem_fit") &&
      !is.null(fit_obj$std_rotated_loadings) &&
      !is.null(fit_obj$se_loadings)) {
    L  <- fit_obj$std_rotated_loadings
    SE <- fit_obj$se_loadings
    return(do.call(rbind, lapply(colnames(L), function(fac)
      data.frame(model  = model_name,
                 factor = tolower(fac),
                 item   = tolower(rownames(L)),
                 R_std  = round(L[, fac], 4),
                 R_se   = round(SE[, fac], 4),
                 R_z    = round(L[, fac] / SE[, fac], 3),
                 R_p    = round(2 * pnorm(-abs(L[, fac] / SE[, fac])), 4),
                 stringsAsFactors = FALSE))))
  }
  # efa() path with Heywood correction: use efa_loadings (corrected) with lavaan SEs
  if (inherits(fit_obj, "esem_fit") && !is.null(fit_obj$efa_loadings) &&
      !is.null(fit_obj$heywood_log)) {
    L   <- fit_obj$efa_loadings
    lav <- fit_obj$lavaan_fit
    ss  <- lavaan::standardizedsolution(lav)
    ss  <- ss[ss$op == "=~", , drop = FALSE]

    # Build SE/z/p lookup from lavaan
    se_lkp <- setNames(ss$se, paste0(ss$lhs, ".", ss$rhs))
    z_lkp  <- setNames(ss$z,  paste0(ss$lhs, ".", ss$rhs))
    p_lkp  <- setNames(ss$pvalue, paste0(ss$lhs, ".", ss$rhs))

    # NA out SEs for clamped items (boundary estimate -> Wald SE not meaningful)
    hw_keys <- paste0(fit_obj$heywood_log$detected$factor, ".",
                      fit_obj$heywood_log$detected$item)

    return(do.call(rbind, lapply(colnames(L), function(fac) {
      keys <- paste0(fac, ".", rownames(L))
      is_hw <- keys %in% hw_keys
      data.frame(model  = model_name,
                 factor = tolower(fac),
                 item   = tolower(rownames(L)),
                 R_std  = round(L[, fac], 4),
                 R_se   = round(ifelse(is_hw, NA_real_, se_lkp[keys]), 4),
                 R_z    = round(ifelse(is_hw, NA_real_, z_lkp[keys]),  3),
                 R_p    = round(ifelse(is_hw, NA_real_, p_lkp[keys]),  4),
                 stringsAsFactors = FALSE)
    })))
  }

  # ESEM (esem_fit wrapper) or plain lavaan CFA
  lav <- if (inherits(fit_obj, "esem_fit")) fit_obj$lavaan_fit else fit_obj
  ss  <- lavaan::standardizedsolution(lav)
  ss  <- ss[ss$op == "=~", , drop = FALSE]
  data.frame(model  = model_name,
             factor = tolower(ss$lhs),
             item   = tolower(ss$rhs),
             R_std  = round(ss$est.std, 4),
             R_se   = round(ss$se,      4),
             R_z    = round(ss$z,       3),
             R_p    = round(ss$pvalue,  4),
             stringsAsFactors = FALSE)
}

.mplus_loadings_long <- function(mplus_obj, model_name) {
  if (is.null(mplus_obj)) return(NULL)
  std  <- mplus_obj$parameters$stdyx.standardized
  if (is.null(std)) return(NULL)
  rows <- std[grepl("\\.BY$", std$paramHeader), , drop = FALSE]
  if (nrow(rows) == 0) return(NULL)
  data.frame(model     = model_name,
             factor    = tolower(sub("\\.BY$", "", rows$paramHeader)),
             item      = tolower(rows$param),
             Mplus_std = round(rows$est,    4),
             Mplus_se  = round(rows$se,     4),
             Mplus_z   = round(rows$est_se, 3),
             Mplus_p   = round(rows$pval,   4),
             stringsAsFactors = FALSE)
}

.merge_loading_tables <- function(r_df, mp_df, spec) {
  primary_map <- setNames(
    tolower(rep(names(spec$factors), lengths(spec$factors))),
    tolower(unlist(spec$factors, use.names = FALSE))
  )
  if (!is.null(mp_df)) {
    merged        <- merge(r_df, mp_df, by = c("model","factor","item"), all.x = TRUE)
    merged$diff_std <- round(merged$R_std - merged$Mplus_std, 4)
  } else {
    merged <- r_df
    merged$Mplus_std <- NA_real_; merged$Mplus_se <- NA_real_
    merged$Mplus_z   <- NA_real_; merged$Mplus_p  <- NA_real_
    merged$diff_std  <- NA_real_
  }
  g_col <- tolower(spec$g_name %||% "G")
  merged$loading_type <- ifelse(
    (!is.na(primary_map[merged$item]) & primary_map[merged$item] == merged$factor) |
      merged$factor == g_col,
    "primary", "cross"
  )
  merged[order(merged$factor, merged$loading_type, merged$item),
         c("model","factor","item","loading_type",
           "R_std","R_se","R_z","R_p",
           "Mplus_std","Mplus_se","Mplus_z","Mplus_p","diff_std")]
}
