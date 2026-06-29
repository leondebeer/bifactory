# Shared roxygen sections (@inheritSection doc_estimator_paths <Section title>).

#' Shared documentation: estimator paths and defaults
#'
#' @name doc_estimator_paths
#' @keywords internal
#'
#' @section Estimator paths:
#' \describe{
#'   \item{Continuous ESEM / B-ESEM}{
#'     \code{\link{esem}} and \code{\link{besem}} use \code{lavaan::cfa()} with a
#'     native \code{efa()} block (integrated MLR estimation and rotation).
#'   }
#'   \item{Ordered ESEM}{
#'     \code{esem(ordered = ...)} calls \code{\link{esem_ordered}}. Default
#'     \code{method = "lavaan"} (lavaan WLSMV, \code{efa()} block, post-hoc rotation).
#'     \code{method = "rotation"} uses a custom DWLS pipeline (polychoric correlations
#'     + \pkg{GPArotation}).
#'   }
#'   \item{Ordered B-ESEM}{
#'     \code{besem(ordered = ...)} calls \code{\link{besem_ordered}}. Default
#'     \code{method = "rotation"} (custom DWLS/WLSMV + orthogonal \code{targetT}; this
#'     is what \code{\link{run_comparison}} uses for ordered data). \code{method = "set-esem"}
#'     fits a lavaan WLSMV bifactor **CFA** with non-primary specific loadings fixed at
#'     zero; it does **not** match Mplus B-ESEM loadings.
#'   }
#'   \item{Multi-group invariance}{
#'     \code{\link{esem_invariance}} uses lavaan multi-group \code{efa()} models plus
#'     explicit syntax patches for ordered B-ESEM. See that help page for scope limits.
#'   }
#' }
#'
#' @section Missing-data defaults:
#' Defaults differ by entry point; pass \code{missing} explicitly when fitting ESEM and
#' B-ESEM separately on the same ordered dataset.
#' \itemize{
#'   \item \code{\link{specify_model}} / \code{\link{run_comparison}}: \code{"pairwise"}
#'     for ordered data, \code{"listwise"} for continuous.
#'   \item \code{\link{esem_ordered}}: \code{"pairwise"}.
#'   \item \code{\link{besem_ordered}}: \code{"listwise"} (pipeline still passes
#'     \code{missing} from the spec when used via \code{run_comparison}).
#'   \item \code{\link{esem}} / \code{\link{besem}} (continuous): \code{"listwise"}.
#' }
#'
#' @section The lavaan_fit slot on custom WLSMV fits:
#' For \code{\link{besem_ordered}(method = "rotation")} and
#' \code{\link{esem_ordered}(method = "rotation")}, \code{$lavaan_fit} is often an
#' **auxiliary** one-factor WLSMV CFA used only to extract DWLS weight matrices---not
#' the fitted ESEM/B-ESEM model. Use \code{\link{std_loadings}}, \code{\link{parameters}},
#' and \code{\link[lavaan:fitMeasures]{fitMeasures}(x)} on the \code{esem_fit} wrapper;
#' \code{\link[lavaan:fitMeasures]{fitMeasures}(x)} reads \code{wlsmv_stats} when present.
#' Do not interpret \code{summary(x$lavaan_fit)}, \code{modindices(x)}, or
#' \code{coef(x)} as the rotated solution unless you know the fit used the lavaan
#' \code{efa()} path (\code{method = "lavaan"} or \code{method = "set-esem"} for B-ESEM).
#'
#' @section Heywood fix and loadings:
#' When \code{heywood_fix = TRUE} (default on single-group \code{\link{esem}} and
#' \code{\link{esem_ordered}}), \code{\link{std_loadings}} may reflect a post-hoc
#' rotation correction while \code{lavaan_fit} still holds the pre-correction lavaan
#' solution. Prefer \code{std_loadings(x)} for reported loadings in that case.
NULL
