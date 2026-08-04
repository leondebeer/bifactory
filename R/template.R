# -- User-facing analysis template helper -------------------------------------

#' Open or Copy the bifactory Analysis Template
#'
#' Opens the shipped analysis template in your editor, or copies it to a
#' destination of your choice. The template walks through the full bifactory
#' pipeline (CFA / ESEM / B-ESEM comparison, reliability, ESEM-within-CFA,
#' factor scores, and optional multi-group invariance) using
#' \code{psych::bfi} as the demo dataset.
#'
#' @param to Optional file path. If supplied, the template is copied there
#'   (with overwrite protection). If \code{NULL} (default), the template is
#'   opened in the editor via \code{\link[utils]{file.edit}}.
#' @param overwrite Logical. Overwrite an existing file at \code{to}?
#'   Default \code{FALSE}.
#'
#' @return Invisibly returns the path to the template (or the destination
#'   when \code{to} is supplied).
#'
#' @examples
#' # Find the template path
#' system.file("templates", "template.R", package = "bifactory")
#'
#' # Copy the template to a file (here a temporary one)
#' dest <- file.path(tempdir(), "my_analysis.R")
#' bifactory_template(to = dest, overwrite = TRUE)
#'
#' \dontrun{
#' # Open the template directly in your editor (interactive session only)
#' bifactory_template()
#' }
#'
#' @export
bifactory_template <- function(to = NULL, overwrite = FALSE) {
  src <- system.file("templates", "template.R", package = "bifactory")
  if (!nzchar(src))
    stop("Template not found in installed package. Reinstall bifactory.",
         call. = FALSE)

  if (is.null(to)) {
    utils::file.edit(src)
    return(invisible(src))
  }

  if (file.exists(to) && !isTRUE(overwrite))
    stop("File already exists at '", to,
         "'. Pass overwrite = TRUE to replace it.", call. = FALSE)

  ok <- file.copy(src, to, overwrite = isTRUE(overwrite))
  if (!ok)
    stop("Could not copy template to '", to, "'.", call. = FALSE)

  message("Template written to: ", normalizePath(to, mustWork = FALSE))
  invisible(to)
}
