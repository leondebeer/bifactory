#' Per-Group Chi-Square Decomposition for Invariance Fits
#'
#' Extracts the per-group contribution to the WLSMV (or ML) chi-square at each
#' invariance level of an \code{esem_invariance} object. Optionally compares
#' against Mplus's \emph{Chi-Square Contribution From Each Group} table, which
#' is useful for pinpointing whether a basin difference between R and Mplus is
#' driven by a single group or distributed across all of them.
#'
#' lavaan's \code{stat.group} is the raw (unscaled) per-group chi-square; Mplus
#' prints the \emph{scaled} per-group contribution that sums to the scaled
#' total. To compare on the same scale, R per-group values are rescaled by the
#' global ratio \code{chisq.scaled / chisq} from \code{\link[lavaan]{fitMeasures}}.
#'
#' @param x An \code{esem_invariance} object.
#' @param mplus_dir Optional directory containing Mplus output files named
#'   \code{besem_inv_<level>.out}. When supplied, the result includes
#'   \code{chisq_M}, \code{mplus_group}, and \code{delta = chisq_R - chisq_M}.
#' @param levels Character vector of levels to include. Defaults to all levels
#'   present in \code{x$models}.
#'
#' @return A data frame of class \code{"chisq_decomp"} with columns
#'   \code{level}, \code{group}, \code{chisq_R} (rescaled to scaled scale), and
#'   when \code{mplus_dir} is supplied, \code{mplus_group}, \code{chisq_M},
#'   \code{delta}.
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
#' inv <- esem_invariance(spec)
#' chisq_decomp(inv)
#' }
#'
#' \dontrun{
#' # Supply a folder of Mplus .out files to add side-by-side deltas.
#' chisq_decomp(inv, mplus_dir = "validation/_bfi_g4_mplus_inv")
#' }
#' @export
chisq_decomp <- function(x, mplus_dir = NULL, levels = NULL) {
  if (!inherits(x, "esem_invariance"))
    stop("`x` must be an `esem_invariance` object.", call. = FALSE)

  available <- names(x$models)
  if (is.null(levels)) levels <- available
  levels <- intersect(levels, available)
  if (length(levels) == 0L)
    stop("No matching levels found in `x$models`.", call. = FALSE)

  rows <- list()
  for (lvl in levels) {
    fit <- x$models[[lvl]]$lavaan_fit
    if (is.null(fit) || !isS4(fit)) {
      warning("chisq_decomp: skipping level '", lvl,
              "' (no valid lavaan fit).", call. = FALSE)
      next
    }

    chi_R_raw <- .chisq_decomp_per_group(fit)
    fmR <- tryCatch(
      lavaan::fitMeasures(fit, c("chisq", "chisq.scaled", "df", "df.scaled")),
      error = function(e) NULL
    )
    if (is.null(fmR) || is.null(chi_R_raw)) {
      warning("chisq_decomp: skipping level '", lvl,
              "' (could not extract chi-square fit measures).", call. = FALSE)
      next
    }

    scale_factor <- if (!is.na(fmR["chisq.scaled"]) && fmR["chisq"] > 0)
      fmR["chisq.scaled"] / fmR["chisq"] else 1
    chi_R <- chi_R_raw * scale_factor

    df <- data.frame(
      level   = lvl,
      group   = names(chi_R),
      chisq_R = unname(as.numeric(chi_R)),
      stringsAsFactors = FALSE
    )

    if (!is.null(mplus_dir)) {
      out_path <- file.path(mplus_dir, sprintf("besem_inv_%s.out", lvl))
      chi_M <- if (file.exists(out_path)) .chisq_decomp_mplus_per_group(out_path) else NULL
      if (!is.null(chi_M)) {
        mapping <- .chisq_decomp_align(chi_R, chi_M)
        df$mplus_group <- vapply(df$group, function(g) {
          m <- mapping[[g]]; if (is.null(m)) NA_character_ else m
        }, character(1))
        df$chisq_M <- vapply(df$mplus_group, function(m) {
          if (is.na(m)) NA_real_ else unname(chi_M[m])
        }, numeric(1))
        df$delta <- df$chisq_R - df$chisq_M
      } else {
        df$mplus_group <- NA_character_
        df$chisq_M     <- NA_real_
        df$delta       <- NA_real_
      }
    }
    rows[[lvl]] <- df
  }

  if (length(rows) == 0L)
    stop("No usable lavaan fits found in `x$models`.", call. = FALSE)

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  attr(out, "has_mplus") <- !is.null(mplus_dir) && "chisq_M" %in% names(out)
  class(out) <- c("chisq_decomp", "data.frame")
  out
}

#' @export
print.chisq_decomp <- function(x, digits = 2L, ...) {
  has_mplus <- isTRUE(attr(x, "has_mplus"))
  cat("Per-group chi-square decomposition (WLSMV scaled scale)\n")
  if (has_mplus) cat("  R per-group rescaled by chisq.scaled / chisq.\n")
  cat("\n")
  for (lvl in unique(x$level)) {
    sub <- x[x$level == lvl, , drop = FALSE]
    cat(sprintf("--- %s ---\n", toupper(lvl)))
    if (has_mplus) {
      cat(sprintf("  %-15s %-15s %12s %12s %12s\n",
                  "R-group", "Mplus-group", "chisq_R", "chisq_M", "delta"))
      for (i in seq_len(nrow(sub))) {
        cat(sprintf("  %-15s %-15s %12.*f %12.*f %12.*f\n",
                    sub$group[i],
                    if (is.na(sub$mplus_group[i])) "??" else sub$mplus_group[i],
                    digits, sub$chisq_R[i],
                    digits, sub$chisq_M[i],
                    digits, sub$delta[i]))
      }
      cat(sprintf("  %-15s %-15s %12.*f %12.*f %12.*f\n",
                  "SUM", "",
                  digits, sum(sub$chisq_R, na.rm = TRUE),
                  digits, sum(sub$chisq_M, na.rm = TRUE),
                  digits, sum(sub$delta,   na.rm = TRUE)))
    } else {
      cat(sprintf("  %-15s %12s\n", "group", "chisq_R"))
      for (i in seq_len(nrow(sub))) {
        cat(sprintf("  %-15s %12.*f\n", sub$group[i], digits, sub$chisq_R[i]))
      }
      cat(sprintf("  %-15s %12.*f\n", "SUM",
                  digits, sum(sub$chisq_R, na.rm = TRUE)))
    }
    cat("\n")
  }
  invisible(x)
}

# ---- Internal helpers -------------------------------------------------------

.chisq_decomp_per_group <- function(fit) {
  tt <- fit@test
  std <- tt[[1]]
  if (is.null(std$stat.group)) return(NULL)
  g_lbls <- lavaan::lavInspect(fit, "group.label")
  if (length(g_lbls) == 0L) g_lbls <- as.character(seq_along(std$stat.group))
  setNames(as.numeric(std$stat.group), g_lbls)
}

.chisq_decomp_mplus_per_group <- function(out_path) {
  L <- readLines(out_path, warn = FALSE)
  i <- grep("Chi-Square Contribution From Each Group", L, fixed = TRUE)
  if (length(i) == 0L) return(NULL)
  i <- i[1]
  out <- list()
  k <- i + 2L
  while (k <= length(L)) {
    line <- trimws(L[k])
    if (line == "") { if (length(out) > 0L) break else { k <- k + 1L; next } }
    if (grepl("^Chi-Square", line)) break
    parts <- strsplit(line, "\\s+")[[1]]
    if (length(parts) >= 2L) {
      val <- suppressWarnings(as.numeric(parts[length(parts)]))
      lbl <- paste(parts[-length(parts)], collapse = "_")
      if (!is.na(val)) out[[lbl]] <- val
    }
    k <- k + 1L
  }
  unlist(out)
}

.chisq_decomp_align <- function(chi_R, chi_M) {
  R_remaining <- names(chi_R); M_remaining <- names(chi_M)
  pairs <- list()
  while (length(R_remaining) > 0L && length(M_remaining) > 0L) {
    best_d <- Inf; best <- NULL
    for (r in R_remaining) for (m in M_remaining) {
      d <- abs(chi_R[r] - chi_M[m])
      if (d < best_d) { best_d <- d; best <- c(r, m) }
    }
    pairs[[best[1]]] <- best[2]
    R_remaining <- setdiff(R_remaining, best[1])
    M_remaining <- setdiff(M_remaining, best[2])
  }
  pairs
}
