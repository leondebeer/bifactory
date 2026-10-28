#' Generate Mplus Syntax for ESEM with Target Rotation
#'
#' Generates a complete Mplus .inp file for ESEM with target rotation,
#' using the correct Mplus syntax: separate BY statements per factor,
#' each with its own (*1) label. Primary items always appear first,
#' followed by cross-loading items targeted to zero.
#'
#' Correct output format:
#' \preformatted{
#' EX BY
#'   batEX1 batEX2 ... batEX8       <- primary items first
#'   batMD1~0 ... batCI5~0 (*1);    <- cross-loadings after
#'
#' MD BY
#'   batMD1 batMD2 ... batMD5       <- primary items first
#'   batEX1~0 ... batCI5~0 (*1);    <- cross-loadings after
#' }
#'
#' @param factors Named list mapping factor names to their primary indicator
#'   names. Example: \code{list(EX = c("y1","y2"), MD = c("y3","y4"))}.
#' @param cfa_factors Optional named list of additional CFA factors.
#' @param regressions Optional character vector of regression statements.
#' @param covariances Optional character vector of covariance statements.
#' @param data_file Character. Name of the data file. Default \code{"mydata.dat"}.
#' @param missing_code Missing value code. Default \code{999}. Use \code{"."}
#'   for files written by \code{MplusAutomation::prepareMplusData()}, which
#'   writes \code{MISSING = .;}.
#' @param estimator Character. Default \code{"WLSMV"} when \code{ordered = TRUE},
#'   otherwise \code{"MLR"}.
#' @param output_path Character. Full path to write the .inp file. Returns
#'   syntax string invisibly if \code{NULL}.
#' @param ordered Logical. \code{TRUE} declares the factor indicators
#'   \code{CATEGORICAL} and uses the theta parameterization, as
#'   \code{\link{esem_ordered}} does. Default \code{FALSE}.
#'
#' @return The Mplus syntax as a character string (invisibly). Writes file
#'   if \code{output_path} is supplied.
#'
#' @examples
#' # Pure syntax generation: builds a Mplus .inp string, no Mplus install needed.
#' syntax <- generate_mplus_syntax(
#'   factors = list(
#'     Visual  = c("x1", "x2", "x3"),
#'     Textual = c("x4", "x5", "x6"),
#'     Speed   = c("x7", "x8", "x9")
#'   ),
#'   output_path = file.path(tempdir(), "esem_model.inp")
#' )
#' cat(syntax)
#'
#' @export
generate_mplus_syntax <- function(factors,
                                   cfa_factors  = NULL,
                                   regressions  = NULL,
                                   covariances  = NULL,
                                   data_file    = "mydata.dat",
                                   missing_code = 999,
                                   estimator    = if (ordered) "WLSMV" else "MLR",
                                   output_path  = NULL,
                                   ordered      = FALSE) {

  if (!is.list(factors) || is.null(names(factors)))
    stop("`factors` must be a named list.", call. = FALSE)

  factor_names   <- names(factors)
  all_esem_items <- unique(unlist(factors))

  cfa_items <- if (!is.null(cfa_factors)) unique(unlist(cfa_factors)) else character(0)

  # Extract covariate names from regression statements
  covariates <- character(0)
  if (!is.null(regressions)) {
    all_factor_names <- c(factor_names, names(cfa_factors))
    for (reg in regressions) {
      parts <- strsplit(reg, " ON ")[[1]]
      if (length(parts) == 2) {
        rhs_vars <- trimws(strsplit(parts[2], " ")[[1]])
        rhs_vars <- rhs_vars[rhs_vars != ""]
        covariates <- unique(c(covariates, setdiff(rhs_vars, all_factor_names)))
      }
    }
  }

  all_vars <- c(all_esem_items, cfa_items, covariates)

  # VARIABLE section
  names_line   <- paste0(.mplus_wrap("  NAMES = ",        all_vars), ";")
  usevars_line <- paste0(.mplus_wrap("  USEVARIABLES = ", all_vars), ";")
  cat_line     <- if (ordered)
    paste0(.mplus_wrap("  CATEGORICAL = ", c(all_esem_items, cfa_items)), ";\n") else ""

  # BY blocks: primary items first, then cross-loadings~0, then (*1)
  by_blocks <- vapply(factor_names, function(f)
    .mplus_esem_by(f, factors[[f]], setdiff(all_esem_items, factors[[f]])),
    FUN.VALUE = character(1))

  by_section <- paste(by_blocks, collapse = "\n\n")

  # CFA factors
  cfa_section <- ""
  if (!is.null(cfa_factors)) {
    cfa_blocks <- vapply(names(cfa_factors), function(f) {
      paste0(.mplus_wrap(paste0("  ", f, " BY "), cfa_factors[[f]]), ";")
    }, FUN.VALUE = character(1))
    cfa_section <- paste0(
      "\n  ! CFA factors\n",
      paste(cfa_blocks, collapse = "\n")
    )
  }

  # Regressions
  reg_section <- ""
  if (!is.null(regressions))
    reg_section <- paste0(
      "\n  ! Structural paths\n",
      paste(paste0("  ", regressions, ";"), collapse = "\n")
    )

  # Covariances
  cov_section <- ""
  if (!is.null(covariances))
    cov_section <- paste0(
      "\n  ! Residual covariances\n",
      paste(paste0("  ", covariances, ";"), collapse = "\n")
    )

  syntax <- paste0(
    "TITLE: ESEM Model (Target Rotation)\n\n",
    "DATA:\n  FILE = ", data_file, ";\n\n",
    "VARIABLE:\n",
    names_line, "\n",
    usevars_line, "\n",
    cat_line,
    .mplus_missing_line(missing_code), "\n\n",
    "ANALYSIS:\n",
    "  ESTIMATOR = ", estimator, ";\n",
    "  ROTATION = TARGET;",
    .mplus_conv_lines(), "\n",
    if (ordered) "  PARAMETERIZATION = THETA;\n" else "",
    "\n",
    "MODEL:\n",
    "  ! Separate BY statement per factor, each with (*1)\n",
    "  ! Primary items listed first; cross-loadings marked ~0 after\n\n",
    by_section,
    cfa_section,
    reg_section,
    cov_section,
    "\n\nOUTPUT:\n  STDYX;\n  MODINDICES(10);\n"
  )

  if (!is.null(output_path)) {
    writeLines(syntax, output_path)
    message("Mplus syntax written to: ", output_path)
  }

  invisible(syntax)
}


# -- Internal helper -----------------------------------------------------------

# One rotated BY statement: primary items, then cross-loadings targeted to
# zero, then the (*1) rotation label; wrapped so no line passes Mplus's
# 90-column limit.
.mplus_esem_by <- function(f, primary, cross = character(0)) {
  body <- .mplus_wrap("    ", primary)
  if (length(cross))
    body <- paste0(body, "\n", .mplus_wrap("    ", paste0(cross, "~0")))
  paste0("  ", f, " BY\n", body, " (*1);")
}

# Mplus reads "." as missing only through MISSING = .; MISSING ARE ALL (.)
# is an input error.  MplusAutomation::prepareMplusData() writes "." for NA.
.mplus_missing_line <- function(missing_code) {
  if (identical(as.character(missing_code), ".")) "  MISSING = .;"
  else paste0("  MISSING ARE ALL (", missing_code, ");")
}

.mplus_wrap <- function(prefix, vars, width = 70) {
  lines   <- character(0)
  indent  <- paste(rep(" ", nchar(prefix)), collapse = "")
  current <- prefix

  for (v in vars) {
    candidate <- paste0(current, v, " ")
    if (nchar(candidate) > width && current != prefix) {
      lines   <- c(lines, trimws(current, which = "right"))
      current <- paste0(indent, v, " ")
    } else {
      current <- candidate
    }
  }
  lines <- c(lines, trimws(current, which = "right"))
  paste(lines, collapse = "\n")
}
