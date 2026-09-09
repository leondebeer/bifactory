#' Alignment Ratio Check for ICM-CFA Specification
#'
#' Tests whether a standard CFA (ICM-CFA) is appropriately specified by
#' computing alignment ratios from the manifest correlation matrix. Based on
#' Mehrvarz & Rouder (2026), who prove that in a correctly specified ICM-CFA,
#' alignment ratios must be invariant across all admissible item quadruples,
#' with the common value equal to the squared latent correlation phi^2.
#'
#' @section What are alignment ratios?:
#' For any two items \eqn{\ell, \ell'} in cluster A and two items \eqn{k, k'}
#' in cluster B, the alignment ratio is:
#' \deqn{Q(\ell, \ell', k, k') = \frac{r_{\ell k} \cdot r_{\ell' k'}}{r_{\ell \ell'} \cdot r_{k k'}}}
#' Under a correctly specified ICM-CFA, all such ratios equal phi^2 and lie in
#' \eqn{[0, 1]} (Mehrvarz & Rouder, 2026, Eq. 5). Dispersion in the ratios --
#' measured by their log-scale standard deviation \eqn{s = \mathrm{sd}(\log Q)}
#' -- signals misspecification: either misassignment (items in the wrong
#' cluster) or misalignment (cross-loadings exist but are fixed to zero).
#'
#' @param data A \code{data.frame} or matrix of observed scores. May also be a
#'   correlation matrix (set \code{is_cor = TRUE}).
#' @param clusters A named character vector mapping item names to factor names.
#'   Example: \code{c(y1 = "F1", y2 = "F1", y3 = "F2", y4 = "F2")}.
#'   Alternatively, a named list: \code{list(F1 = c("y1","y2"), F2 = c("y3","y4"))}.
#' @param cfa_fit Optional. A fitted \code{lavaan} CFA object. When supplied,
#'   CFA-estimated latent correlations are compared against the alignment-implied
#'   phi to quantify inflation.
#' @param is_cor Logical. Is \code{data} already a correlation matrix?
#'   Default \code{FALSE}.
#' @param min_within_r Numeric. Alignment ratios whose denominator contains a
#'   within-cluster correlation below this value are excluded (near-zero
#'   within-cluster correlations make ratios numerically unstable). Default
#'   \code{0.05}.
#' @param log_sd_thresholds Named numeric vector of log-scale SD cutoffs for
#'   the slight / moderate / high misalignment classification. Defaults to the
#'   empirical terciles reported by Mehrvarz & Rouder (2026, p.24):
#'   \code{c(slight = 0.64, moderate = 1.6)}. Values of \eqn{s = \mathrm{sd}(
#'   \log Q)} below \code{slight} indicate an ICM-CFA approximately consistent
#'   with the data; above \code{moderate} indicate substantial misspecification.
#'
#' @return An object of class \code{"alignment_check"} (a list) containing:
#' \describe{
#'   \item{\code{pair_results}}{A \code{data.frame} with one row per factor
#'     pair, giving: geometric mean of alignment ratios, log-scale SD, implied phi,
#'     CFA-estimated phi (if supplied), inflation percentage, and verdict.}
#'   \item{\code{all_ratios}}{A named list of raw alignment ratio vectors,
#'     one element per factor pair.}
#'   \item{\code{recommendation}}{Character string: overall recommendation.}
#'   \item{\code{cor_matrix}}{The manifest correlation matrix used.}
#'   \item{\code{clusters}}{Resolved cluster assignments (named list).}
#'   \item{\code{call}}{The matched call.}
#' }
#'
#' @details
#' ## Interpreting the log-scale dispersion (sd of log Q)
#' \eqn{s = \mathrm{sd}(\log Q)} is the primary diagnostic. The log scale is
#' natural because alignment ratios are multiplicative objects. Mehrvarz &
#' Rouder (2026, p.24) report empirical terciles of \eqn{s} from their bifactor
#' misalignment simulation:
#' \itemize{
#'   \item \strong{s < 0.64} -- Slight misalignment. ICM-CFA is approximately
#'     consistent with the data; latent correlations are roughly trustworthy.
#'   \item \strong{s 0.64-1.6} -- Moderate misalignment. ICM-CFA is partially
#'     misspecified; interpret latent correlations with caution.
#'   \item \strong{s > 1.6} -- High misalignment (top tercile of the paper's
#'     simulation, up to ~4.4). ICM-CFA is substantially misspecified and
#'     latent correlations are likely inflated. Consider ESEM.
#' }
#' These cutpoints are the 33rd and 66th percentiles of \eqn{s} observed by
#' Mehrvarz & Rouder across a broad range of simulated loading configurations
#' (with \eqn{\mathrm{sd}(\log \kappa_j) \in [0, 3]}), not hard decision
#' boundaries; adjust \code{log_sd_thresholds} if your application calls for
#' stricter or looser cutoffs.
#'
#' ## Types of misspecification detected
#' \itemize{
#'   \item \strong{Misassignment}: items are assigned to the wrong cluster.
#'     Mehrvarz & Rouder prove this necessarily inflates phi. The alignment
#'     ratios split into three plateaus at phi^2, 1, and 1/phi^2.
#'   \item \strong{Misalignment}: cluster assignments are correct but
#'     cross-loadings exist (or the proportionality constraint is violated in
#'     a bifactor-like DGP). The ratios are dispersed rather than invariant.
#' }
#'
#' ## Number of alignment ratios per pair
#' Mehrvarz & Rouder (2026, p.14) count one alignment ratio per admissible
#' index quadruple: \eqn{n_Q = \binom{m_A}{2}\binom{m_B}{2}}. Because each
#' unordered quadruple \eqn{\{\ell, \ell'\} \times \{k, k'\}} admits two
#' between-cluster matchings (\eqn{(\ell,k)/(\ell',k')} and
#' \eqn{(\ell,k')/(\ell',k)}), both of which equal \eqn{\varphi^2} under
#' correctly specified ICM-CFA but diverge under misspecification, this
#' function records both matchings per quadruple -- yielding \eqn{2 \binom{m_A}
#' {2}\binom{m_B}{2}} ratios in \code{all_ratios} -- to maximise the
#' information available for the dispersion diagnostic.
#'
#' @references
#' Mehrvarz, M., & Rouder, J. N. (2026). The geometry and brittleness of latent
#' correlations in confirmatory factor analysis.
#'
#' @seealso \code{\link{esem}} for the recommended follow-up under moderate or
#'   high misalignment.
#'
#' @examples
#' data("HolzingerSwineford1939", package = "lavaan")
#' d <- HolzingerSwineford1939[, paste0("x", 1:9)]
#'
#' # Vector form of cluster assignment
#' clusters <- c(
#'   x1 = "Visual",  x2 = "Visual",  x3 = "Visual",
#'   x4 = "Textual", x5 = "Textual", x6 = "Textual",
#'   x7 = "Speed",   x8 = "Speed",   x9 = "Speed"
#' )
#'
#' \donttest{
#' # Run alignment check (data only)
#' check <- alignment_check(d, clusters)
#' print(check)
#'
#' # Run with a fitted CFA to quantify inflation
#' cfa_model <- "
#'   Visual  =~ x1 + x2 + x3
#'   Textual =~ x4 + x5 + x6
#'   Speed   =~ x7 + x8 + x9
#' "
#' cfa_fit <- lavaan::cfa(cfa_model, data = d, std.lv = TRUE)
#' check   <- alignment_check(d, clusters, cfa_fit = cfa_fit)
#' print(check)
#' }
#'
#' @export
alignment_check <- function(data,
                             clusters,
                             cfa_fit            = NULL,
                             is_cor             = FALSE,
                             min_within_r       = 0.05,
                             log_sd_thresholds  = c(slight = 0.64, moderate = 1.6)) {

  mc <- match.call()

  # -- 0. Parse cluster argument ---------------------------------------------
  clusters <- .parse_clusters(clusters)   # always returns named list

  factor_names <- names(clusters)
  nfactors     <- length(factor_names)
  if (nfactors < 2)
    stop("Need at least 2 factors to compute alignment ratios.", call. = FALSE)

  all_items <- unlist(clusters, use.names = FALSE)

  # -- 1. Compute manifest correlation matrix --------------------------------
  if (is_cor) {
    R <- as.matrix(data)
    if (!isSymmetric(R))
      stop("`data` is not symmetric. Set is_cor = FALSE if passing raw data.",
           call. = FALSE)
  } else {
    if (!is.data.frame(data) && !is.matrix(data))
      stop("`data` must be a data.frame, matrix, or correlation matrix.",
           call. = FALSE)
    missing_items <- setdiff(all_items, colnames(data))
    if (length(missing_items))
      stop("Items not found in `data`: ",
           paste(missing_items, collapse = ", "), call. = FALSE)
    n_obs <- colSums(!is.na(data[, all_items, drop = FALSE]))
    empty_items <- all_items[n_obs == 0L]
    if (length(empty_items))
      stop("Items with no non-missing values (cannot compute correlations): ",
           paste(empty_items, collapse = ", "), call. = FALSE)
    R <- cor(data[, all_items], use = "pairwise.complete.obs")
  }

  # -- 2. Extract CFA latent correlations (if provided) ---------------------
  cfa_phi <- NULL
  if (!is.null(cfa_fit)) {
    cfa_phi <- tryCatch(
      .extract_cfa_phi(cfa_fit, factor_names),
      error = function(e) {
        warning("Could not extract CFA latent correlations: ", conditionMessage(e),
                call. = FALSE)
        NULL
      }
    )
  }

  # -- 3. Compute alignment ratios for all factor pairs ---------------------
  pairs        <- utils::combn(factor_names, 2, simplify = FALSE)
  all_ratios   <- list()
  pair_results <- vector("list", length(pairs))

  for (p in seq_along(pairs)) {
    fA <- pairs[[p]][1]
    fB <- pairs[[p]][2]
    pair_key <- paste(fA, fB, sep = "-")

    itemsA <- clusters[[fA]]
    itemsB <- clusters[[fB]]

    ratios <- .compute_alignment_ratios(R, itemsA, itemsB, min_within_r)
    all_ratios[[pair_key]] <- ratios

    # Summary statistics
    if (length(ratios) == 0) {
      warning("No admissible alignment ratios for pair ", pair_key,
              " (too few items or near-zero within-cluster correlations).",
              call. = FALSE)
      pair_results[[p]] <- data.frame(
        pair          = pair_key,
        n_ratios      = 0L,
        geom_mean     = NA_real_,
        log_sd        = NA_real_,
        phi_alignment = NA_real_,
        phi_cfa       = NA_real_,
        inflation_pct = NA_real_,
        verdict       = "insufficient data",
        stringsAsFactors = FALSE
      )
      next
    }

    pos_ratios <- ratios[ratios > 0]             # log requires positive values
    gm         <- exp(mean(log(pos_ratios)))     # geometric mean (spec Sec.8.1)
    log_sd     <- if (length(pos_ratios) > 1)   # log-scale dispersion (spec Sec.8.2)
                    sd(log(pos_ratios)) else 0
    phi_align  <- sqrt(gm)                      # implied phi from alignment

    # CFA phi for this pair
    phi_cfa <- if (!is.null(cfa_phi)) cfa_phi[fA, fB] else NA_real_

    inflation_pct <- if (!is.na(phi_cfa) && phi_align > 0) {
      round((phi_cfa - phi_align) / phi_align * 100, 1)
    } else NA_real_

    verdict <- .alignment_verdict(log_sd, log_sd_thresholds)

    pair_results[[p]] <- data.frame(
      pair          = pair_key,
      n_ratios      = length(ratios),
      geom_mean     = round(gm, 4),
      log_sd        = round(log_sd, 4),
      phi_alignment = round(phi_align, 4),
      phi_cfa       = if (!is.na(phi_cfa)) round(phi_cfa, 4) else NA_real_,
      inflation_pct = inflation_pct,
      verdict       = verdict,
      stringsAsFactors = FALSE
    )
  }

  pair_results <- do.call(rbind, pair_results)

  # -- 4. Overall recommendation ---------------------------------------------
  verdicts       <- pair_results$verdict
  recommendation <- .overall_recommendation(verdicts, log_sd_thresholds)

  # -- 5. Return -------------------------------------------------------------
  structure(
    list(
      pair_results   = pair_results,
      all_ratios     = all_ratios,
      recommendation = recommendation,
      cor_matrix        = R,
      clusters          = clusters,
      log_sd_thresholds = log_sd_thresholds,
      call              = mc
    ),
    class = "alignment_check"
  )
}


# -- Print method --------------------------------------------------------------

#' Print Method for alignment_check
#'
#' @param x An \code{alignment_check} object.
#' @param ... Ignored.
#' @return Invisibly returns \code{x}; called for the side effect of printing
#'   the alignment-ratio summary.
#' @export
print.alignment_check <- function(x, ...) {
  cat("\nAlignment Ratio Check  (Mehrvarz & Rouder, 2026)\n")

  pr       <- x$pair_results
  has_cfa  <- !all(is.na(pr$phi_cfa))
  t1       <- x$log_sd_thresholds[1]
  t2       <- x$log_sd_thresholds[2]

  icon_map <- c(slight = "OK", moderate = " ~", high = " X",
                `insufficient data` = " ?")

  fmt_num <- function(v, d = 3) ifelse(is.na(v), "   -", formatC(v, digits = d, format = "f"))

  cols <- list(
    unname(ifelse(is.na(icon_map[pr$verdict]), " ?", icon_map[pr$verdict])),
    pr$pair,
    ifelse(is.na(pr$n_ratios), "-", format(pr$n_ratios)),
    fmt_num(pr$geom_mean),
    fmt_num(pr$log_sd),
    fmt_num(pr$phi_alignment)
  )
  headers <- c("", "Pair", "n", "gm", "sd(logQ)", "phi")
  if (has_cfa) {
    cols <- c(cols, list(
      fmt_num(pr$phi_cfa),
      ifelse(is.na(pr$inflation_pct), "   -", sprintf("%+5.1f", pr$inflation_pct))
    ))
    headers <- c(headers, "phi_cfa", "infl%")
  }
  cols    <- c(cols, list(toupper(ifelse(pr$verdict == "insufficient data",
                                         "n/a", pr$verdict))))
  headers <- c(headers, "Verdict")

  mat     <- do.call(cbind, cols)
  widths  <- pmax(nchar(headers), apply(mat, 2, function(col) max(nchar(col))))

  fmt_row <- function(row) {
    cells <- mapply(function(v, w) formatC(v, width = w, flag = "-"),
                    row, widths, USE.NAMES = FALSE)
    paste(cells, collapse = "  ")
  }
  cat("\n  ", fmt_row(headers), "\n", sep = "")
  cat("  ", strrep("-", sum(widths) + 2 * (length(widths) - 1)), "\n", sep = "")
  for (i in seq_len(nrow(mat))) cat("  ", fmt_row(mat[i, ]), "\n", sep = "")

  cat("\n  Key:  n         = # admissible alignment ratios Q\n")
  cat("        gm        = geometric mean of Q  (= phi^2 under correct ICM-CFA)\n")
  cat("        sd(logQ)  = log-scale dispersion of Q  (primary misalignment diagnostic)\n")
  cat("        phi       = sqrt(gm), alignment-implied latent correlation\n")
  if (has_cfa)
  cat("        phi_cfa   = latent correlation from fitted CFA\n        infl%     = (phi_cfa - phi) / phi, inflation of CFA vs alignment\n")
  cat(sprintf("\n  Thresholds (sd log Q): < %.2f slight,  < %.2f moderate\n", t1, t2))
  cat(paste(strwrap(paste0("Overall: ", x$recommendation), width = getOption("width", 80L),
                    indent = 2, exdent = 11), collapse = "\n"), "\n", sep = "")
  cat("\n  Reference:\n")
  cat("    Mehrvarz, M., & Rouder, J. N. (2026). The geometry and brittleness of\n")
  cat("    latent correlations in confirmatory factor analysis [Preprint].\n")
  cat("    OSF Preprints. https://doi.org/10.31234/osf.io/95enc_v3\n")

  invisible(x)
}


#' Plot Sorted Alignment Ratios
#'
#' Produces a sorted alignment ratio plot (as in Figure 2 of Mehrvarz & Rouder
#' 2026) for one or all factor pairs. Flat horizontal spread = invariant
#' (ICM-CFA tenable). Dispersion = misspecification.
#'
#' @param x An \code{alignment_check} object.
#' @param pair Character. Name of a specific factor pair (e.g., \code{"F1-F2"}).
#'   If \code{NULL} (default), plots all pairs in a grid.
#' @param ... Ignored.
#'
#' @return Called for its side effect of drawing the alignment-ratio plot;
#'   invisibly returns \code{NULL}.
#' @export
plot.alignment_check <- function(x, pair = NULL, ...) {
  pairs_to_plot <- if (!is.null(pair)) pair else names(x$all_ratios)

  n <- length(pairs_to_plot)
  old_par <- graphics::par(mfrow = c(ceiling(n / 2), min(n, 2)),
                           mar   = c(4, 4, 3, 1))
  on.exit(graphics::par(old_par))

  for (pname in pairs_to_plot) {
    ratios <- x$all_ratios[[pname]]
    pr_row <- x$pair_results[x$pair_results$pair == pname, ]

    if (is.null(ratios) || length(ratios) == 0) {
      graphics::plot.new()
      graphics::title(main = paste(pname, "(no ratios)"))
      next
    }

    sorted <- sort(ratios)
    ymax   <- max(sorted, 1.1, na.rm = TRUE) * 1.1

    graphics::plot(
      seq_along(sorted), sorted,
      pch  = 19, cex = 0.6,
      xlab = "Ordered Alignment Ratios",
      ylab = expression(Q[i]),
      main = paste("Pair:", pname),
      ylim = c(0, ymax),
      col  = "black"
    )

    # Geometric mean line (implied phi^2)
    graphics::abline(h = pr_row$geom_mean, col = "gray40", lty = 2, lwd = 1.5)

    # CFA phi^2 line (if available)
    if (!is.na(pr_row$phi_cfa)) {
      graphics::abline(h = pr_row$phi_cfa^2, col = "red", lty = 2, lwd = 1.5)
    }

    # Verdict annotation
    graphics::legend(
      "topleft",
      legend  = c(
        paste("Geom. mean =", round(pr_row$geom_mean, 3)),
        if (!is.na(pr_row$phi_cfa))
          paste("CFA phi^2 =", round(pr_row$phi_cfa^2, 3)) else NULL,
        paste("sd(log Q) =", round(pr_row$log_sd, 3)),
        paste("Verdict:", pr_row$verdict)
      ),
      bty     = "n",
      cex     = 0.75,
      text.col = c("gray40",
                   if (!is.na(pr_row$phi_cfa)) "red" else NULL,
                   "black", "black")
    )
  }
}


# -- Internal helpers ----------------------------------------------------------

.parse_clusters <- function(clusters) {
  if (is.list(clusters)) {
    if (is.null(names(clusters)))
      stop("`clusters` list must be named.", call. = FALSE)
    return(clusters)
  }
  if (is.character(clusters) && !is.null(names(clusters))) {
    # Named vector: names = items, values = factors
    factor_names <- unique(clusters)
    out <- lapply(factor_names, function(f) names(clusters)[clusters == f])
    names(out) <- factor_names
    return(out)
  }
  stop("`clusters` must be a named character vector or named list.",
       call. = FALSE)
}


.compute_alignment_ratios <- function(R, itemsA, itemsB, min_within_r) {
  mA <- length(itemsA)
  mB <- length(itemsB)

  if (mA < 2 || mB < 2) return(numeric(0))

  ratios <- numeric(0)

  # All unordered pairs within A and within B
  pairsA <- utils::combn(itemsA, 2, simplify = FALSE)
  pairsB <- utils::combn(itemsB, 2, simplify = FALSE)

  for (pA in pairsA) {
    l  <- pA[1]; lp <- pA[2]
    r_within_A <- R[l, lp]

    if (abs(r_within_A) < min_within_r) next  # skip near-zero

    for (pB in pairsB) {
      k  <- pB[1]; kp <- pB[2]
      r_within_B <- R[k, kp]

      if (abs(r_within_B) < min_within_r) next

      denom <- r_within_A * r_within_B
      if (abs(denom) < 1e-10) next

      # Two between-cluster matchings per unordered quadruple (Mehrvarz &
      # Rouder 2026, p.13-14: Q(l, l', k, k') vs Q(l, l', k', k)). Both equal
      # phi^2 under correctly-specified ICM-CFA (Eq. 5, 0 <= Q <= 1) but
      # diverge under misspecification, so we keep both.
      Q1 <- (R[l, k]  * R[lp, kp]) / denom   # matching 1: l-k, l'-k'
      Q2 <- (R[l, kp] * R[lp, k])  / denom   # matching 2: l-k', l'-k

      if (is.finite(Q1)) ratios <- c(ratios, Q1)
      if (is.finite(Q2)) ratios <- c(ratios, Q2)
    }
  }

  ratios
}


.extract_cfa_phi <- function(cfa_fit, factor_names) {
  phi <- lavaan::lavInspect(cfa_fit, "cor.lv")

  # Match factor names -- lavaan may use different ordering
  fn_in_model <- rownames(phi)
  matched     <- intersect(factor_names, fn_in_model)

  if (length(matched) < 2) {
    warning("Factor names in `cfa_fit` do not match those derived from `clusters`. ",
            "Provide matching factor names in both.", call. = FALSE)
    return(NULL)
  }

  phi[matched, matched]
}


.alignment_verdict <- function(log_sd, thresholds) {
  if (is.na(log_sd))              return("insufficient data")
  if (log_sd < thresholds[1])     return("slight")
  if (log_sd < thresholds[2])     return("moderate")
  return("high")
}


.overall_recommendation <- function(verdicts, thresholds) {
  t_slight   <- thresholds[1]
  t_moderate <- thresholds[2]
  n_high     <- sum(verdicts == "high",     na.rm = TRUE)
  n_mod      <- sum(verdicts == "moderate", na.rm = TRUE)
  n_slight   <- sum(verdicts == "slight",   na.rm = TRUE)
  n_tot      <- length(verdicts[!is.na(verdicts) & verdicts != "insufficient data"])

  if (n_tot == 0)
    return("Insufficient data to assess alignment.")

  if (n_high == 0 && n_mod == 0)
    return(sprintf(
      paste0("All factor pairs show SLIGHT misalignment (sd(log Q) < %.2f). ",
             "ICM-CFA is approximately consistent with the data; ",
             "latent correlations are roughly trustworthy."),
      t_slight
    ))

  if (n_high == 0 && n_mod > 0)
    return(sprintf(
      paste0("%d of %d factor pair(s) show MODERATE misalignment ",
             "(sd(log Q) between %.2f and %.2f). ",
             "ICM-CFA is partially misspecified; interpret latent correlations ",
             "with caution. Fit the ESEM model as a sensitivity check."),
      n_mod, n_tot, t_slight, t_moderate
    ))

  sprintf(
    paste0("%d of %d factor pair(s) show HIGH misalignment (sd(log Q) > %.2f). ",
           "ICM-CFA is substantially misspecified and latent correlations are ",
           "likely inflated. Fit the ESEM model, and re-examine cluster ",
           "assignments for misassignment."),
    n_high, n_tot, t_moderate
  )
}
