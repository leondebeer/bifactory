# Standardized Cronbach's alpha with automatic reversal detection.
# Delegates to psych::alpha(check.keys = TRUE), which reflects items whose
# first-principal-component loading is negative. This correctly handles
# mixed-keying scales (e.g. BFI) where a mean-inter-item-r heuristic mis-
# flags straight items as reversed. Returns NA_real_ for fewer than 2
# matched items or any psych::alpha failure.
.compute_alpha <- function(data, items) {
  matched <- colnames(data)[tolower(colnames(data)) %in% tolower(items)]
  if (length(matched) < 2) return(NA_real_)
  X <- apply(data[, matched, drop = FALSE], 2, as.numeric)
  a <- tryCatch(
    suppressWarnings(suppressMessages(
      psych::alpha(X, check.keys = TRUE, warnings = FALSE)
    )),
    error = function(e) NULL
  )
  if (is.null(a) || is.null(a$total$std.alpha)) return(NA_real_)
  round(a$total$std.alpha, 3)
}

#' Compute Reliability Indices for CFA, ESEM, and B-ESEM
#'
#' Computes McDonald's omega reliability indices for all three models in a
#' pipeline result, following the bifactor reporting framework recommended by
#' Morin, Arens & Marsh (2016) and Rodriguez, Reise & Haviland (2016).
#'
#' For the **total composite** (all items):
#' \itemize{
#'   \item \code{omega_total} -- total reliability (rotation-invariant for
#'     orthogonal models).
#'   \item \code{omega_H} -- hierarchical omega: G factor's contribution to
#'     total-score reliability (B-ESEM only).
#'   \item \code{ECV} -- explained common variance: G's share of all common
#'     variance (B-ESEM only).
#'   \item \code{H(G)} -- construct replicability of G (B-ESEM only).
#' }
#'
#' For each **subscale** (items of specific factor \eqn{F_s}):
#' \itemize{
#'   \item \code{omega_S} -- specific factor's unique contribution to subscale
#'     reliability: \eqn{(\sum_{i \in s} \lambda_{s,i})^2 /
#'     [(\sum_{i \in s} \lambda_{s,i})^2 + \sum_{i \in s} \psi_i]}.
#'     Uses target loadings and item residuals only; G does not enter the
#'     denominator because \eqn{\psi_i} already has G partialled out.
#'   \item \code{omega_sub} -- total reliability of the subscale sum score
#'     (G + specific factor combined).
#'   \item \code{omega_H_sub} -- G's contribution to subscale reliability.
#'   \item \code{ECV_s} -- G's share of common variance within the subscale.
#'   \item \code{H(Fs)} -- construct replicability of the specific factor.
#' }
#'
#' @param results An \code{esem_comparison_pipeline} object from
#'   \code{\link{run_comparison}}.
#'
#' @return An object of class \code{"reliability_indices"} -- a named list with
#'   elements \code{cfa}, \code{esem}, \code{besem}, each containing the
#'   computed indices for that model, plus \code{alpha} (Cronbach's alpha per
#'   subscale and for G). Pass to \code{print()} for a formatted table.
#'
#' @references
#' McDonald, R. P. (1999). \emph{Test theory: A unified treatment}. Erlbaum.
#'
#' Rodriguez, A., Reise, S. P., & Haviland, M. G. (2016). Evaluating bifactor
#' models: Calculating and interpreting statistical indices.
#' \emph{Psychological Methods}, \emph{21}(2), 137-150.
#'
#' Morin, A. J. S., Arens, A. K., & Marsh, H. W. (2016). A bifactor
#' exploratory structural equation modeling framework for the identification of
#' distinct sources of construct-relevant psychometric multidimensionality.
#' \emph{Structural Equation Modeling}, \emph{23}(1), 116-139.
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
#' print(indices)
#' }
#'
#' @export
compute_indices <- function(results) {

  if (!inherits(results, "esem_comparison_pipeline"))
    stop("`results` must be an esem_comparison_pipeline object from run_comparison().",
         call. = FALSE)

  factor_items <- results$spec$factors
  fit_cfa      <- results$fit_cfa
  fit_esem     <- results$fit_esem
  fit_besem    <- results$fit_besem

  # Cronbach's alpha -- computed once from raw data, model-agnostic.
  # Reverse-keyed items are auto-detected from the correlation matrix
  # (items whose mean inter-item r is negative are reflected), matching
  # psych::alpha(check.keys=TRUE). No model information is used.
  raw_data  <- results$spec$data
  all_items <- unlist(factor_items, use.names = FALSE)

  alpha_sub <- lapply(factor_items, function(items) .compute_alpha(raw_data, items))
  names(alpha_sub) <- names(factor_items)
  alpha_G   <- .compute_alpha(raw_data, all_items)
  alpha_vals <- c(list(G = alpha_G), alpha_sub)

  # PUC (Reise et al., 2013): proportion of off-diagonal correlations that are
  # between-subscale (uncontaminated by specific factors).
  # Purely structural -- depends only on subscale sizes, not on loadings.
  p             <- length(unlist(factor_items, use.names = FALSE))
  total_pairs   <- p * (p - 1L) / 2L
  within_pairs  <- sum(vapply(factor_items, function(items) {
    n <- length(items); n * (n - 1L) / 2L
  }, numeric(1L)))
  puc <- round((total_pairs - within_pairs) / total_pairs, 3)

  # Match items case-insensitively against matrix rownames
  .sub_rows <- function(item_vec, mat)
    which(tolower(rownames(mat)) %in% tolower(item_vec))

  # Extract Lambda, psi, Phi from lavaan or esem_fit
  .extract_lambda_psi_phi <- function(fit_obj) {
    lav_obj <- if (inherits(fit_obj, "esem_fit")) fit_obj$lavaan_fit else fit_obj
    ss      <- lavaan::standardizedsolution(lav_obj)
    lam     <- ss[ss$op == "=~", , drop = FALSE]
    items   <- unique(lam$rhs)
    factors <- unique(lam$lhs)
    Lambda  <- matrix(0, length(items), length(factors),
                      dimnames = list(items, factors))
    for (i in seq_len(nrow(lam)))
      Lambda[lam$rhs[i], lam$lhs[i]] <- lam$est.std[i]
    psi <- 1 - rowSums(Lambda^2)
    Phi <- tryCatch({
      m <- lavaan::lavInspect(lav_obj, "cor.lv")
      if (is.matrix(m) && nrow(m) == length(factors)) m else diag(length(factors))
    }, error = function(e) diag(length(factors)))
    list(Lambda = Lambda, psi = psi, Phi = Phi)
  }

  # Omega for correlated-factor models (CFA / ESEM)
  .omega_lavaan <- function(fit_obj, model_label) {
    e      <- .extract_lambda_psi_phi(fit_obj)
    Lambda <- e$Lambda; psi <- e$psi; Phi <- e$Phi
    c_vec  <- colSums(Lambda)
    numer  <- as.numeric(t(c_vec) %*% Phi %*% c_vec)
    H      <- sapply(colnames(Lambda), function(f) {
      r <- sum(Lambda[, f]^2 / psi); r / (1 + r)
    })
    sub_list <- lapply(names(factor_items), function(s) {
      idx <- .sub_rows(factor_items[[s]], Lambda)
      if (length(idx) < 2) return(NULL)
      Ls  <- Lambda[idx, , drop = FALSE]; ps <- psi[idx]
      cs  <- colSums(abs(Ls))
      n_s <- as.numeric(t(cs) %*% Phi %*% cs); v_s <- n_s + sum(ps)
      list(factor        = s,
           omega_sub     = round(n_s / v_s, 3),
           omega_h_g     = NA_real_,
           omega_specific = NA_real_,
           ecv_sub       = NA_real_,
           H_specific    = NA_real_)
    })
    names(sub_list) <- names(factor_items)
    list(model       = model_label,
         omega_total = round(numer / (numer + sum(psi)), 3),
         omega_h_g   = NA_real_,
         ecv         = NA_real_,
         H           = round(H, 3),
         subscales   = sub_list)
  }

  # Omega for B-ESEM rotation method (orthogonal bifactor)
  # Following Rodriguez et al. (2016) and Morin et al. (2016):
  #   omega_S = cs^2 / vs  (specific factor's share of total subscale variance)
  #   omega_H = cGs^2 / vs  (G's share of total subscale variance)
  #   omega_sub = (cGs^2 + cs^2) / vs  (total subscale reliability)
  #   omega_S + omega_H = omega_sub  (exact partition)
  .omega_besem_rot <- function(fit_b) {
    L     <- fit_b$std_rotated_loadings
    gname <- fit_b$g_name %||% colnames(L)[1]
    specs <- setdiff(colnames(L), gname)
    psi   <- pmax(1 - rowSums(L^2), 1e-6)
    c_vec <- colSums(L)
    total_var <- sum(c_vec^2) + sum(psi)
    H_vals <- sapply(colnames(L), function(f) { r <- sum(L[,f]^2/psi); r/(1+r) })

    sub_list <- lapply(specs, function(s) {
      idx  <- .sub_rows(factor_items[[s]], L)
      if (length(idx) < 2) return(NULL)
      Ls   <- L[idx, , drop = FALSE]; ps <- psi[idx]
      # Use abs() so reverse-keyed items (negative loadings) contribute positively
      # to the subscale sum score, matching the researcher's practice of
      # reverse-scoring before summing. Rodriguez et al. (2016) assume consistently-
      # keyed items; abs() generalises the formula to mixed-keyed subscales.
      cGs  <- sum(abs(Ls[, gname])); cs <- sum(abs(Ls[, s]))
      vs   <- cGs^2 + cs^2 + sum(ps)   # subscale composite variance
      list(factor         = s,
           omega_sub      = round((cGs^2 + cs^2) / vs, 3),
           omega_h_g      = round(cGs^2 / vs, 3),
           omega_specific = round(cs^2 / vs, 3),
           ecv_sub        = round(sum(Ls[,gname]^2) /
                                    (sum(Ls[,gname]^2) + sum(Ls[,s]^2)), 3),
           H_specific     = round({r <- sum(Ls[,s]^2/ps); r/(1+r)}, 3))
    })
    names(sub_list) <- specs

    list(model       = "BESEM",
         omega_total = round(sum(c_vec^2) / total_var, 3),
         omega_h_g   = round(c_vec[gname]^2 / total_var, 3),
         ecv         = if (sum(L^2) < 1e-14) NA_real_
                       else round(sum(L[, gname]^2) / sum(L^2), 3),
         puc         = puc,
         H           = round(H_vals, 3),
         subscales   = sub_list)
  }

  # Omega for Mplus BESEM -- build Lambda from STDYX loadings, psi from 1 - rowSums(L^2)
  .omega_besem_mplus <- function(mplus_besem) {
    std <- mplus_besem$parameters$stdyx.standardized
    if (is.null(std)) return(NULL)
    by_rows <- std[grepl("\\.BY$", std$paramHeader), , drop = FALSE]
    if (nrow(by_rows) == 0) return(NULL)

    all_factors <- unique(sub("\\.BY$", "", by_rows$paramHeader))
    all_items   <- unique(by_rows$param)
    L <- matrix(0, length(all_items), length(all_factors),
                dimnames = list(all_items, all_factors))
    for (i in seq_len(nrow(by_rows)))
      L[by_rows$param[i], sub("\\.BY$", "", by_rows$paramHeader[i])] <- by_rows$est[i]

    # Identify general factor: "G" if present, else first column
    gname <- if ("G" %in% colnames(L)) "G" else colnames(L)[1]
    specs <- setdiff(colnames(L), gname)
    psi   <- 1 - rowSums(L^2)
    psi   <- pmax(psi, 1e-6)   # guard against rounding below zero
    c_vec <- colSums(L)
    total_var <- sum(c_vec^2) + sum(psi)
    H_vals <- sapply(colnames(L), function(f) { r <- sum(L[,f]^2/psi); r/(1+r) })

    sub_list <- lapply(specs, function(s) {
      idx  <- .sub_rows(factor_items[[s]], L)
      if (length(idx) < 2) return(NULL)
      Ls   <- L[idx, , drop = FALSE]; ps <- psi[idx]
      cGs  <- sum(abs(Ls[, gname])); cs <- sum(abs(Ls[, s]))
      vs   <- cGs^2 + cs^2 + sum(ps)
      list(factor         = s,
           omega_sub      = round((cGs^2 + cs^2) / vs, 3),
           omega_h_g      = round(cGs^2 / vs, 3),
           omega_specific = round(cs^2 / vs, 3),
           ecv_sub        = round(sum(Ls[,gname]^2) /
                                    (sum(Ls[,gname]^2) + sum(Ls[,s]^2)), 3),
           H_specific     = round({r <- sum(Ls[,s]^2/ps); r/(1+r)}, 3))
    })
    names(sub_list) <- specs

    list(model       = "BESEM_Mplus",
         omega_total = round(sum(c_vec^2) / total_var, 3),
         omega_h_g   = round(c_vec[gname]^2 / total_var, 3),
         ecv         = if (sum(L^2) < 1e-14) NA_real_
                       else round(sum(L[, gname]^2) / sum(L^2), 3),
         puc         = puc,
         H           = round(H_vals, 3),
         subscales   = sub_list)
  }

  # Omega for continuous besem() -- loadings live in lavaan standardized solution
  .omega_besem_lavaan <- function(fit_b) {
    e     <- .extract_lambda_psi_phi(fit_b)
    L     <- e$Lambda
    psi   <- e$psi
    gname <- fit_b$g_name
    specs <- names(fit_b$specific_factors)
    col_order <- intersect(c(gname, specs), colnames(L))
    L <- L[, col_order, drop = FALSE]
    psi <- pmax(psi, 1e-6)
    c_vec     <- colSums(L)
    total_var <- sum(c_vec^2) + sum(psi)
    H_vals    <- sapply(colnames(L), function(f) { r <- sum(L[,f]^2/psi); r/(1+r) })
    sub_list <- lapply(specs, function(s) {
      idx <- .sub_rows(factor_items[[s]], L)
      if (length(idx) < 2) return(NULL)
      Ls  <- L[idx, , drop = FALSE]; ps <- psi[idx]
      cGs <- sum(abs(Ls[, gname])); cs <- sum(abs(Ls[, s]))
      vs   <- cGs^2 + cs^2 + sum(ps)
      list(factor         = s,
           omega_sub      = round((cGs^2 + cs^2) / vs, 3),
           omega_h_g      = round(cGs^2 / vs, 3),
           omega_specific = round(cs^2 / vs, 3),
           ecv_sub        = round(sum(Ls[,gname]^2) /
                                    (sum(Ls[,gname]^2) + sum(Ls[,s]^2)), 3),
           H_specific     = round({r <- sum(Ls[,s]^2/ps); r/(1+r)}, 3))
    })
    names(sub_list) <- specs
    list(model       = "BESEM",
         omega_total = round(sum(c_vec^2) / total_var, 3),
         omega_h_g   = round(c_vec[gname]^2 / total_var, 3),
         ecv         = if (sum(L^2) < 1e-14) NA_real_
                       else round(sum(L[, gname]^2) / sum(L^2), 3),
         puc         = puc,
         H           = round(H_vals, 3),
         subscales   = sub_list)
  }

  # Dispatch
  res_besem <- if (inherits(fit_besem, "esem_fit") &&
                   !is.null(fit_besem$std_rotated_loadings))
    .omega_besem_rot(fit_besem)                      # ordered path (WLSMV rotation)
  else if (inherits(fit_besem, "besem_fit") &&
           !is.null(fit_besem$g_name))
    .omega_besem_lavaan(fit_besem)                   # continuous path (MLR lavaan)
  else
    .omega_lavaan(fit_besem, "BESEM")                # fallback

  # Mplus BESEM omega (if available)
  res_besem_mplus <- if (!is.null(results$mplus_results) &&
                         !is.null(results$mplus_results$besem))
    tryCatch(.omega_besem_mplus(results$mplus_results$besem),
             error = function(e) NULL)
  else NULL

  result <- list(
    cfa          = .omega_lavaan(fit_cfa,  "CFA"),
    esem         = .omega_lavaan(fit_esem, "ESEM"),
    besem        = res_besem,
    besem_mplus  = res_besem_mplus
  )
  result$alpha <- alpha_vals
  class(result) <- "reliability_indices"
  result
}


#' Print a reliability_indices Object
#'
#' Displays McDonald's omega and Cronbach's alpha indices following the
#' Morin / Rodriguez et al. reporting framework.
#'
#' @param x A \code{reliability_indices} object from \code{\link{compute_indices}}.
#' @param ... Ignored.
#' @return Invisibly returns \code{x}; called for the side effect of printing
#'   the reliability table.
#' @export
print.reliability_indices <- function(x, ...) {
  W   <- 100
  r_b <- x$besem

  .wrap <- function(..., indent = 2, exdent = 4) {
    lines <- strwrap(paste0(...), width = W, indent = indent, exdent = exdent)
    cat(paste(lines, collapse = "\n"), "\n", sep = "")
  }
  .rule <- function(ch = "=") cat(strrep(ch, W), "\n", sep = "")
  .sec  <- function(title) {
    cat(sprintf("\n-- %s %s\n", title,
                strrep("-", max(1L, W - nchar(title, type = "chars") - 4L))))
  }

  .rule()
  cat(" McDonald's Omega and Bifactor Indices\n")
  cat(" McDonald (1999); Morin, Arens & Marsh (2016); Rodriguez, Reise & Haviland (2016)\n")
  .rule()

  has_mplus   <- !is.null(x$besem_mplus)
  show_models <- c("cfa", "esem", "besem")
  if (has_mplus) show_models <- c(show_models, "besem_mplus")

  # -- McDonald's omega -- subscales (CFA / ESEM) ---------------------------------
  cfa_subs  <- x$cfa$subscales
  esem_subs <- x$esem$subscales
  if (!is.null(cfa_subs) && length(cfa_subs) > 0) {
    .sec("McDonald's w (omega) -- subscales (CFA / ESEM)")
    fw  <- max(8L, max(nchar(names(cfa_subs))))
    fmt <- sprintf("  %%-%ds  %%-10s  %%-10s  %%-7s\n", fw)
    cat(sprintf(fmt, "Factor",          "w (CFA)",   "w (ESEM)",  "alpha"))
    cat(sprintf(fmt, strrep("-", fw),   "---------", "---------", "-------"))
    for (nm in names(cfa_subs)) {
      s_cfa  <- cfa_subs[[nm]]
      s_esem <- if (!is.null(esem_subs)) esem_subs[[nm]] else NULL
      alph   <- if (!is.null(x$alpha) && !is.na(x$alpha[[nm]]))
        sprintf("%.3f", x$alpha[[nm]]) else "--"
      cat(sprintf(fmt,
                  nm,
                  if (!is.null(s_cfa)  && !is.na(s_cfa$omega_sub))  sprintf("%.3f", s_cfa$omega_sub)  else "--",
                  if (!is.null(s_esem) && !is.na(s_esem$omega_sub)) sprintf("%.3f", s_esem$omega_sub) else "--",
                  alph))
    }
    cat("\n")
    .wrap("w = McDonald's omega for the subscale sum score: proportion of subscale score variance attributable to all common factors for those items.")
  }

  # -- McDonald's omega -- total composite ---------------------------------------
  .sec("McDonald's w (omega) -- total composite (all items)")
  cat(sprintf("  %-14s  %-11s  %-10s  %-8s  %-6s\n",
              "Model", "w (total)*", "w-H (G)+", "ECV(G)", "H(G)"))
  cat(sprintf("  %-14s  %-11s  %-10s  %-8s  %-6s\n",
              "--------------", "-----------", "----------", "--------", "------"))
  for (nm in show_models) {
    r <- x[[nm]]
    if (is.null(r)) next
    g_h <- if (!is.null(r$H) && "G" %in% names(r$H)) sprintf("%.3f", r$H[["G"]]) else "--"
    cat(sprintf("  %-14s  %-11s  %-10s  %-8s  %-6s\n",
                r$model,
                sprintf("%.3f", r$omega_total),
                if (!is.na(r$omega_h_g)) sprintf("%.3f", r$omega_h_g) else "--",
                if (!is.na(r$ecv))       sprintf("%.3f", r$ecv)       else "--",
                g_h))
  }
  cat("\n")
  .wrap("* w (total) = McDonald's omega for the total composite score -- the proportion of total score variance attributable to all common factors combined. Rotation-invariant.")
  .wrap("+ w-H (G) = omega-hierarchical: the proportion of total score variance attributable to G alone. Rotation-sensitive; B-ESEM only.")
  if (!is.null(x$alpha) && !is.na(x$alpha[["G"]]))
    .wrap(sprintf("Cronbach's alpha (all items, standardized) = %.3f", x$alpha[["G"]]))

  # -- Dimensionality assessment -------------------------------------------------
  if (!is.null(r_b$puc) && !is.na(r_b$omega_h_g) && !is.na(r_b$ecv)) {
    puc <- r_b$puc; ecv <- r_b$ecv; omH <- r_b$omega_h_g; omT <- r_b$omega_total
    H_G <- if (!is.null(r_b$H) && "G" %in% names(r_b$H)) r_b$H[["G"]] else NA

    .sec("Dimensionality assessment (Reise et al., 2013)")
    cat(sprintf("  PUC = %.3f  |  ECV(G) = %.3f  |  w-H = %.3f  |  w (total) = %.3f\n\n",
                puc, ecv, omH, omT))

    essentially_unidim <- (puc >= 0.80) || (puc < 0.80 && ecv >= 0.70)
    strong_G <- omH >= 0.75
    H_G_good <- !is.na(H_G) && H_G >= 0.80

    if (essentially_unidim && strong_G) {
      if (puc >= 0.80) {
        .wrap(sprintf("PUC = %.3f (>= .80): most item correlations are between-subscale, limiting specific-factor contamination of the correlation matrix.", puc))
      } else {
        .wrap(sprintf("PUC = %.3f (< .80) but ECV = %.3f (>= .70): although many correlations are within-subscale, G still explains the clear majority of common variance, overriding the low PUC concern.", puc, ecv))
      }
      .wrap(sprintf("ECV = %.3f: G accounts for %.0f%% of all common variance. w-H = %.3f means %.0f%% of total score reliability is attributable to G.", ecv, 100 * ecv, omH, 100 * omH / omT))
      if (H_G_good) .wrap(sprintf("H(G) = %.3f: G is highly replicable across samples.", H_G))
      cat("\n")
      .wrap("Verdict based on these data: The data are essentially unidimensional. Specific factors are present but secondary; treating total scores as a unidimensional composite involves minimal construct-relevant multidimensionality. Total scores are the primary interpretable unit; subscale scores may be reported for descriptive purposes but add limited incremental information.")

    } else if (!essentially_unidim && strong_G) {
      .wrap(sprintf("PUC = %.3f (< .80) and ECV = %.3f (< .70): the data are not essentially unidimensional -- specific factors carry meaningful variance beyond G. However, G is well-defined and dominant: w-H = %.3f means %.0f%% of total score reliability is attributable to G.", puc, ecv, omH, 100 * omH / omT))
      if (H_G_good) .wrap(sprintf("H(G) = %.3f: G is highly replicable across samples.", H_G))
      cat("\n")
      .wrap("Verdict based on these data: Multidimensional with a dominant G. Total scores reflect primarily G but carry construct-relevant multidimensionality; subscale scores add interpretive value.")

    } else {
      .wrap(sprintf("PUC = %.3f, ECV = %.3f, w-H = %.3f: G does not dominate the common variance structure. Specific factors account for a substantial share of item intercorrelations, and w-H = %.3f means only %.0f%% of total score reliability is attributable to a single general dimension.", puc, ecv, omH, omH, 100 * omH / omT))
      if (!is.na(H_G))
        .wrap(sprintf("H(G) = %.3f: G %s replicate well across samples.", H_G,
                      if (H_G >= 0.80) "would" else "may not"))
      cat("\n")
      .wrap("Verdict based on these data: Substantially multidimensional. Summing all items into a total score conflates distinct constructs and is not recommended. Subscale scores are the primary interpretable unit; if a total score is required, its multidimensional nature must be explicitly acknowledged.")
    }
    cat("\n")
    .wrap("Reference: Reise, S. P., Scheines, R., Widaman, K. F., & Haviland, M. G. (2013). Multidimensionality and structural coefficient bias in SEM: A bifactor perspective. Educational and Psychological Measurement, 73, 5-26.",
          indent = 2, exdent = 12)
  }

  # -- McDonald's omega -- subscales ---------------------------------------------
  if (!is.null(r_b$subscales) && length(r_b$subscales) > 0) {
    cat("\n")
    fnames <- vapply(r_b$subscales,
                     function(s) if (is.null(s)) "" else s$factor, character(1))
    fw     <- max(8L, max(nchar(fnames)))
    fdash  <- strrep("-", fw)
    if (has_mplus) {
      r_mp   <- x$besem_mplus
      mp_fmt <- sprintf("  %%-%ds  %%-9s  %%-9s  %%-7s  %%-10s  %%-10s  %%-7s  %%-7s  %%-8s  %%-9s\n", fw)
      cat(sprintf(mp_fmt,
                  "Factor",
                  "w-sub(R)", "w-spec(R)", "alpha",
                  "w-sub(Mp)", "w-spec(Mp)",
                  "w-H(R)", "w-H(Mp)", "ECV-s(R)", "ECV-s(Mp)"))
      cat(sprintf(mp_fmt,
                  fdash,
                  "---------", "---------", "-------",
                  "----------", "----------",
                  "-------", "-------", "--------", "---------"))
      for (s in r_b$subscales) {
        if (is.null(s)) next
        mp_s <- r_mp$subscales[[s$factor]]
        alph <- if (!is.null(x$alpha) && !is.na(x$alpha[[s$factor]]))
          sprintf("%.3f", x$alpha[[s$factor]]) else "--"
        cat(sprintf(mp_fmt,
                    s$factor,
                    sprintf("%.3f", s$omega_sub), sprintf("%.3f", s$omega_specific), alph,
                    if (!is.null(mp_s)) sprintf("%.3f", mp_s$omega_sub)      else "--",
                    if (!is.null(mp_s)) sprintf("%.3f", mp_s$omega_specific) else "--",
                    sprintf("%.3f", s$omega_h_g),
                    if (!is.null(mp_s)) sprintf("%.3f", mp_s$omega_h_g)      else "--",
                    sprintf("%.3f", s$ecv_sub),
                    if (!is.null(mp_s)) sprintf("%.3f", mp_s$ecv_sub)        else "--"))
      }
    } else {
      fmt <- sprintf("  %%-%ds  %%-14s  %%-14s  %%-7s  %%-12s  %%-7s  %%-6s\n", fw)
      cat(sprintf(fmt,
                  "Factor", "w-sub",          "w-spec",         "alpha",
                  "w-H(G|s)",     "ECV-s",   "H(Fs)"))
      cat(sprintf(fmt,
                  fdash,    "--------------", "--------------", "-------",
                  "------------", "-------", "------"))
      for (s in r_b$subscales) {
        if (is.null(s)) next
        alph <- if (!is.null(x$alpha) && !is.na(x$alpha[[s$factor]]))
          sprintf("%.3f", x$alpha[[s$factor]]) else "--"
        cat(sprintf(fmt,
                    s$factor,
                    sprintf("%.3f", s$omega_sub), sprintf("%.3f", s$omega_specific), alph,
                    sprintf("%.3f", s$omega_h_g), sprintf("%.3f", s$ecv_sub),
                    sprintf("%.3f", s$H_specific)))
      }
    }
    cat("\n")
    .wrap("w-sub   = McDonald's omega for the total subscale sum score combining G and the specific factor: (cGs^2 + cs^2) / (cGs^2 + cs^2 + S(psi)). Identity: w-spec + w-H = w-sub. [McDonald, 1970; also see Morin et al., 2016]", exdent = 16)
    .wrap("w-spec  = McDonald's omega for the specific factor's unique contribution to subscale reliability (G partialled out): cs^2 / (cGs^2 + cs^2 + S(psi)). [Rodriguez et al., 2016]", exdent = 16)
    .wrap("w-H(G|s) = G's proportion of subscale score variance (hierarchical omega at the subscale level): cGs^2 / (cGs^2 + cs^2 + S(psi)). [Rodriguez et al., 2016]", exdent = 16)
    .wrap("ECV-s       = G's share of common variance within the subscale: S(lam^2)_G / (S(lam^2)_G + S(lam^2)_s) for subscale items.", exdent = 16)
    .wrap("H(Fs)       = construct replicability of the specific factor: how well the subscale items jointly define their specific factor across samples. Values >= 0.80 indicate a well-defined factor. [Rodriguez et al., 2016]", exdent = 16)
    .wrap("Note: cs = S|lam_s| -- absolute loadings used to handle reverse-keyed items.", exdent = 7)
  }

  # -- H-index -------------------------------------------------------------------
  .sec("H-index (construct replicability -- all models)")
  for (nm in show_models) {
    r <- x[[nm]]
    if (is.null(r) || is.null(r$H) || length(r$H) == 0) next
    cat(sprintf("  %-14s  %s\n", r$model,
                paste(sprintf("%s=%.3f", names(r$H), r$H), collapse = "  ")))
  }
  cat("\n")
  invisible(x)
}


#' Refine B-ESEM Rotation Using Mplus Solution as Warm Start
#'
#' The bifactor target rotation criterion surface has many local minima. R and
#' Mplus sometimes converge to different ones, leading to different partitioning
#' of variance between G and specific factors (while fit indices remain
#' identical because they are rotation-invariant).
#'
#' This function uses orthogonal Procrustes rotation to compute the rotation
#' matrix \eqn{\mathbf{T}} that maps R's unrotated loading matrix toward the
#' Mplus STDYX solution, then evaluates whether the target criterion at that
#' \eqn{\mathbf{T}} is lower than R's best random-start criterion. If so, the
#' better rotation is adopted. This is methodologically sound because:
#' \enumerate{
#'   \item The same published criterion function is still being minimised.
#'   \item Procrustes only provides a well-informed starting point; the final
#'     solution is the converged \code{GPArotation::targetT} optimum from that
#'     start.
#'   \item The approach is analogous to Mansolf and Reise's (2016)
#'     recommendation to use Schmid-Leiman solutions as warm starts.
#' }
#'
#' @param results An \code{esem_comparison_pipeline} object from
#'   \code{\link{run_comparison}} that includes Mplus results
#'   (\code{results$mplus_results} must be non-NULL).
#'
#' @return The same \code{results} object with \code{fit_besem} updated if a
#'   better rotation was found, otherwise unchanged. A message reports whether
#'   the criterion improved. The \code{comparison_table} is also updated.
#'
#' @references
#' Mansolf, M., and Reise, S. P. (2016). Exploratory bifactor analysis: The
#' Schmid-Leiman orthogonalization and Jennrich-Bentler analytic rotations.
#' \emph{Multivariate Behavioral Research}, \emph{51}(5), 698--717.
#'
#' @examples
#' \dontrun{
#' # Requires run_comparison() results that include a Mplus rotation reference
#' # (results$mplus_results non-NULL), so a licensed Mplus install is needed.
#' results <- run_comparison(spec, mplus_folder = tempfile("mplus_"))
#' results <- refine_rotation(results)
#' omega   <- compute_omega(results)
#' print(omega)
#' }
#'
#' @export
refine_rotation <- function(results) {

  if (!inherits(results, "esem_comparison_pipeline"))
    stop("`results` must be an esem_comparison_pipeline object.", call. = FALSE)

  fit_b <- results$fit_besem

  if (is.null(fit_b$unrotated_loadings))
    stop("fit_besem does not contain unrotated_loadings. Re-run run_comparison() ",
         "with the current version of bifactory.", call. = FALSE)

  if (is.null(results$mplus_results) || is.null(results$mplus_results$besem))
    stop("No Mplus BESEM results found. Run with mplus_folder to enable refinement.",
         call. = FALSE)

  # -- Extract Mplus STDYX loading matrix ---------------------------------------
  std     <- results$mplus_results$besem$parameters$stdyx.standardized
  by_rows <- std[grepl("\\.BY$", std$paramHeader), , drop = FALSE]
  if (nrow(by_rows) == 0)
    stop("Could not extract Mplus STDYX loadings from besem output.", call. = FALSE)

  mp_factors <- unique(sub("\\.BY$", "", by_rows$paramHeader))
  mp_items   <- unique(by_rows$param)
  L_mp <- matrix(0, length(mp_items), length(mp_factors),
                 dimnames = list(mp_items, mp_factors))
  for (i in seq_len(nrow(by_rows)))
    L_mp[by_rows$param[i], sub("\\.BY$", "", by_rows$paramHeader[i])] <- by_rows$est[i]

  # Align column order to match R's factor order
  r_factors <- fit_b$factor_names
  shared    <- intersect(r_factors, colnames(L_mp))
  if (length(shared) < ncol(L_mp))
    stop("Factor names in Mplus and R do not match. Cannot compute Procrustes T.",
         call. = FALSE)
  L_mp <- L_mp[, r_factors, drop = FALSE]

  # Align row order to match R's item order (Mplus returns uppercase names)
  r_items    <- fit_b$indicators
  mp_rows_lc <- tolower(rownames(L_mp))
  r_items_lc <- tolower(r_items)
  row_idx    <- match(r_items_lc, mp_rows_lc)
  if (anyNA(row_idx))
    stop("Some R items not found in Mplus STDYX output: ",
         paste(r_items[is.na(row_idx)], collapse=", "), call. = FALSE)
  L_mp <- L_mp[row_idx, , drop = FALSE]
  rownames(L_mp) <- r_items   # restore R-style casing

  # -- Orthogonal Procrustes: find T minimising ||L_unrot @ T - L_mp||_F --------
  L_unrot <- fit_b$unrotated_loadings
  C       <- crossprod(L_mp, L_unrot)    # k  x  k
  sv      <- svd(C)
  T_proc  <- sv$v %*% t(sv$u)           # orthogonal Procrustes solution

  # -- Build target matrix (same convention as Stage 5) -------------------------
  btgt_rot <- results$spec$bifactor_target   # 0/1 matrix (cross=0, primary/G=1)
  btgt_na  <- btgt_rot
  btgt_na[btgt_rot == 1L] <- NA_real_
  cross_mask <- btgt_rot == 0L

  # -- Evaluate criterion at Procrustes T ---------------------------------------
  L_proc_start <- L_unrot %*% T_proc
  rot_proc <- tryCatch(
    suppressWarnings(
      GPArotation::targetT(L_proc_start, Target = btgt_na, maxit = 10000L, eps = 1e-8)
    ),
    error = function(e) NULL
  )

  if (is.null(rot_proc) || !all(is.finite(rot_proc$loadings))) {
    message("refine_rotation: Procrustes start failed to converge -- no change made.")
    return(invisible(results))
  }

  fn_proc <- 0.5 * sum(rot_proc$loadings[cross_mask]^2)
  fn_orig <- fit_b$rotation_criterion

  message(sprintf("Rotation criterion -- original: %.6f  |  Procrustes start: %.6f",
                  fn_orig, fn_proc))

  if (fn_proc >= fn_orig - 1e-8) {
    message("No improvement from Procrustes warm start. Original rotation retained.")
    return(invisible(results))
  }

  message(sprintf("Improvement found (delta = %.6f). Updating rotation.",
                  fn_orig - fn_proc))

  # -- Apply sign correction -----------------------------------------------------
  L_new <- rot_proc$loadings
  colnames(L_new) <- fit_b$factor_names
  rownames(L_new) <- fit_b$indicators
  k <- ncol(L_new)
  for (j in seq_len(k)) {
    prim_idx <- if (j == 1L) seq_len(nrow(L_new)) else which(btgt_rot[, j] == 1L)
    if (sum(L_new[prim_idx, j]) < 0) L_new[, j] <- -L_new[, j]
  }

  # -- Update fit_besem in place --------------------------------------------------
  results$fit_besem$rotated_loadings     <- L_new
  results$fit_besem$std_rotated_loadings <- L_new   # already standardised (polychoric fit)
  results$fit_besem$rotation_criterion   <- fn_proc
  results$fit_besem$rotation             <- paste0(
    fit_b$rotation, " [Procrustes-refined from Mplus]"
  )

  # SE recomputation would require re-running Stage 6 with the new rotation;
  # flag that SEs are from the original rotation.
  results$fit_besem$se_loadings_note <-
    "SEs computed at original rotation; refine_rotation() updated loadings only."

  invisible(results)
}

#' @rdname print.reliability_indices
#' @export
print.omega_result <- function(x, ...) print.reliability_indices(x, ...)


#' Deprecated: use \code{\link{compute_indices}} instead
#'
#' @param results An \code{esem_comparison_pipeline} object.
#' @return A \code{reliability_indices} object.
#' @export
compute_omega <- function(results) {
  .Deprecated("compute_indices")
  compute_indices(results)
}
