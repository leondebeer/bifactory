#' bifactory: Bifactor ESEM, CFA, and Measurement Invariance with Mplus-Compatible Output
#'
#' Mplus-compatible CFA, ESEM, and bifactor ESEM (B-ESEM) in R, with
#' multi-group measurement invariance, ESEM-within-CFA conversion, and
#' McDonald's omega reliability. Built on \pkg{lavaan} and \pkg{GPArotation}.
#'
#' Estimator paths differ by data type and model; see \code{?esem}, \code{?besem},
#' \code{?esem_ordered}, and \code{?besem_ordered}. Validation summary:
#' \code{system.file("VALIDATION.md", package = "bifactory")}.
#'
#' @keywords internal
#' @importFrom stats complete.cases cor cov factanal na.omit nlminb pchisq
#'   pnorm rnorm runif sd setNames uniroot
#' @importFrom utils combn write.csv file.edit
#' @importFrom graphics abline legend par plot plot.new title
#' @importFrom methods getMethod isGeneric existsMethod
"_PACKAGE"
