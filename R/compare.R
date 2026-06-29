#' Compare ESEM Against a Standard CFA
#'
#' Fits a user-specified CFA model on the same data and compares it against an
#' \code{esem_fit} object using fit indices and, where applicable, a chi-square
#' difference test (when models are nested).
#'
#' @param esem_model An \code{esem_fit} object from \code{\link{esem}}.
#' @param cfa_model A lavaan model string for the comparison CFA. All items
#'   must be the same as in the ESEM.
#' @param data A \code{data.frame} used to fit the CFA. If \code{NULL}
#'   (default), the data are re-extracted from the \code{esem_model} lavaan
#'   object.
#' @param estimator Estimator for the CFA model. Defaults to the same estimator
#'   used in \code{esem_model}.
#' @param ... Additional arguments passed to \code{lavaan::cfa()}.
#'
#' @return A list of class \code{"esem_comparison"} with:
#' \describe{
#'   \item{\code{fit_table}}{A \code{data.frame} of fit indices for both models.}
#'   \item{\code{esem_fit}}{The original \code{esem_fit} object.}
#'   \item{\code{cfa_fit}}{The fitted \code{lavaan} CFA object.}
#'   \item{\code{lavtest}}{Output of \code{lavaan::lavTestLRT()} if models are nested, else \code{NULL}.}
#' }
#'
#' @seealso \code{\link{esem}}
#'
#' @examples
#' \dontrun{
#' # Fit ESEM
#' esem_result <- esem(mydata, nfactors = 3)
#'
#' # Define comparison CFA (no cross-loadings)
#' cfa_model <- "
#'   F1 =~ y1 + y2 + y3 + y4 + y5
#'   F2 =~ y6 + y7 + y8 + y9 + y10
#'   F3 =~ y11 + y12 + y13 + y14 + y15
#' "
#' comparison <- esem_compare(esem_result, cfa_model, data = mydata)
#' print(comparison)
#' }
#'
#' @importFrom lavaan cfa lavTestLRT fitMeasures
#' @export
esem_compare <- function(esem_model,
                          cfa_model,
                          data      = NULL,
                          estimator = NULL,
                          ...) {

  if (!inherits(esem_model, "esem_fit"))
    stop("`esem_model` must be an esem_fit object.", call. = FALSE)

  if (!is.character(cfa_model))
    stop("`cfa_model` must be a lavaan model string.", call. = FALSE)

  # Extract data and estimator from esem_model if not provided
  lav <- esem_model$lavaan_fit
  if (!isS4(lav))
    stop("`esem_model$lavaan_fit` is not a valid lavaan object; ",
         "the model may have failed to fit.", call. = FALSE)
  if (is.null(data))
    data <- lavaan::lavInspect(lav, "data")[[1]]  # works for single group

  if (is.null(estimator))
    estimator <- lav@Options$estimator

  # Fit CFA
  message("Fitting comparison CFA...")
  cfa_fit <- tryCatch(
    lavaan::cfa(cfa_model, data = data, estimator = estimator,
                std.lv = TRUE, ...),
    error = function(e) stop("CFA fitting failed: ", conditionMessage(e), call. = FALSE)
  )

  # Collect fit indices
  indices <- c("npar", "chisq", "df", "pvalue", "cfi", "tli",
                "rmsea", "rmsea.ci.lower", "rmsea.ci.upper", "srmr", "aic", "bic")

  esem_fi <- lavaan::fitMeasures(lav, fit.measures = indices)
  cfa_fi  <- lavaan::fitMeasures(cfa_fit, fit.measures = indices)

  fit_table <- data.frame(
    index = indices,
    ESEM  = round(as.numeric(esem_fi), 4),
    CFA   = round(as.numeric(cfa_fi), 4),
    stringsAsFactors = FALSE
  )

  # Interpret direction of differences
  fit_table$Better <- apply(fit_table, 1, function(row) {
    idx  <- row["index"]
    esem <- suppressWarnings(as.numeric(row["ESEM"]))
    cfa  <- suppressWarnings(as.numeric(row["CFA"]))
    if (is.na(esem) || is.na(cfa)) return("")
    # Higher is better: cfi, tli
    # Lower is better: chisq, rmsea, srmr, aic, bic
    higher_better <- c("cfi", "tli")
    lower_better  <- c("chisq", "rmsea", "rmsea.ci.lower", "rmsea.ci.upper",
                        "srmr", "aic", "bic")
    if (idx %in% higher_better) {
      if (esem > cfa) "ESEM" else if (cfa > esem) "CFA" else "="
    } else if (idx %in% lower_better) {
      if (esem < cfa) "ESEM" else if (cfa < esem) "CFA" else "="
    } else ""
  })

  # LRT if models might be nested (ESEM has more df used = fewer df residual)
  lavtest <- tryCatch(
    lavaan::lavTestLRT(cfa_fit, lav),
    error   = function(e) NULL,
    warning = function(w) NULL
  )

  out <- structure(
    list(
      fit_table = fit_table,
      esem_fit  = esem_model,
      cfa_fit   = cfa_fit,
      lavtest   = lavtest
    ),
    class = "esem_comparison"
  )

  invisible(out)
}


#' Print Method for esem_comparison
#'
#' @param x An \code{esem_comparison} object.
#' @param ... Ignored.
#' @return Invisibly returns \code{x}; called for the side effect of printing
#'   the comparison table.
#' @export
print.esem_comparison <- function(x, ...) {
  cat("\n=======================================\n")
  cat(" bifactory: ESEM vs CFA Comparison\n")
  cat("=======================================\n\n")

  cat("Fit Indices:\n")
  print(x$fit_table, row.names = FALSE)

  if (!is.null(x$lavtest)) {
    cat("\nChi-square Difference Test (CFA vs ESEM):\n")
    print(x$lavtest)
  } else {
    cat("\n(Models may not be nested -- chi-square difference test not computed.)\n")
  }

  invisible(x)
}
