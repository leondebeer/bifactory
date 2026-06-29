#' Create a Target Loading Matrix for Target Rotation
#'
#' Constructs a target matrix suitable for use with \code{rotation = "target"}
#' in \code{\link{esem}}. Items listed in \code{keys} are assigned a target
#' value of \code{1} for their primary factor; all other cells are set to
#' \code{0} (penalised to be near zero) or \code{NA} (free, no penalty).
#'
#' @param keys A named list where each element is a vector of item indices
#'   (integers) or item names (characters) that are hypothesised to load
#'   primarily on that factor.  Names become factor names.
#'   Example: \code{list(F1 = 1:5, F2 = 6:10, F3 = 11:15)}.
#' @param nitems Integer. Total number of items. Required when \code{keys}
#'   uses integer indices and \code{item_names} is not supplied.
#' @param item_names Optional character vector of item names of length
#'   \code{nitems}. When supplied, \code{keys} may use either names or indices.
#'   Becomes the row names of the target matrix.
#' @param cross_loading_value Numeric or \code{NA}. Value assigned to cells
#'   that are \emph{not} the primary factor. Use \code{0} (default) to penalise
#'   cross-loadings toward zero, or \code{NA} to leave them completely free
#'   (soft target rotation).
#'
#' @return A numeric matrix (items  x  factors) with row names set to item names
#'   and column names set to factor names from \code{keys}.
#'
#' @details
#' ## Target vs. Soft Target Rotation
#'
#' - **Hard target** (\code{cross_loading_value = 0}): Cross-loadings are
#'   penalised toward zero. Use this when you have strong theory.
#' - **Soft target** (\code{cross_loading_value = NA}): Only primary loadings
#'   are targeted; cross-loadings are completely free. Use this when you are
#'   less certain about the zero pattern.
#'
#' ## Items Loading on Multiple Factors
#'
#' An item can appear in multiple \code{keys} entries if it is expected to have
#' meaningful loadings on more than one factor. In that case its target value
#' will be \code{1} for both factors and \code{cross_loading_value} elsewhere.
#'
#' @seealso \code{\link{esem}}
#'
#' @examples
#' # Simple 15-item, 3-factor target (integer keys)
#' tgt <- make_target(
#'   keys   = list(Extrav = 1:5, Agree = 6:10, Open = 11:15),
#'   nitems = 15
#' )
#'
#' # Named items
#' tgt2 <- make_target(
#'   keys       = list(F1 = c("y1", "y2", "y3"), F2 = c("y4", "y5", "y6")),
#'   item_names = paste0("y", 1:6)
#' )
#'
#' # Soft target (cross-loadings free)
#' tgt_soft <- make_target(
#'   keys                 = list(F1 = 1:5, F2 = 6:10),
#'   nitems               = 10,
#'   cross_loading_value  = NA
#' )
#'
#' @export
make_target <- function(keys,
                        nitems              = NULL,
                        item_names          = NULL,
                        cross_loading_value = 0) {

  if (!is.list(keys) || length(keys) == 0)
    stop("`keys` must be a non-empty named list.", call. = FALSE)

  nfactors     <- length(keys)
  factor_names <- names(keys)
  if (is.null(factor_names) || any(factor_names == ""))
    stop("All elements of `keys` must be named (use names() to assign factor names).",
         call. = FALSE)

  # Resolve item indices vs. names
  if (!is.null(item_names)) {
    nitems     <- length(item_names)
    # Convert any character keys to integer indices
    keys <- lapply(keys, function(k) {
      if (is.character(k)) {
        idx <- match(k, item_names)
        if (any(is.na(idx)))
          stop("Some item names in `keys` were not found in `item_names`: ",
               paste(k[is.na(idx)], collapse = ", "), call. = FALSE)
        idx
      } else {
        k
      }
    })
  } else {
    if (is.null(nitems)) {
      # Infer nitems from maximum index
      all_idx <- unlist(lapply(keys, function(k) {
        if (is.integer(k) || is.numeric(k)) k else NA
      }))
      if (any(is.na(all_idx)))
        stop("Supply `nitems` or `item_names` when `keys` uses character indices.",
             call. = FALSE)
      nitems <- max(all_idx)
      message("Inferring nitems = ", nitems, " from maximum index in `keys`.")
    }
    item_names <- paste0("y", seq_len(nitems))
  }

  # Validate indices
  for (f in seq_len(nfactors)) {
    idx <- keys[[f]]
    if (!is.numeric(idx))
      stop("Keys for factor '", factor_names[f], "' must be numeric indices.",
           call. = FALSE)
    if (any(idx < 1) || any(idx > nitems))
      stop("Index out of range for factor '", factor_names[f], "'. ",
           "Indices must be between 1 and ", nitems, ".", call. = FALSE)
  }

  # Build matrix
  tmat <- matrix(cross_loading_value, nrow = nitems, ncol = nfactors,
                 dimnames = list(item_names, factor_names))

  for (j in seq_len(nfactors)) {
    tmat[keys[[j]], j] <- 1
  }

  class(tmat) <- c("esem_target", "matrix", "array")
  tmat
}


#' Print a Target Rotation Matrix
#'
#' Displays a \code{make_target()} result in a compact, readable format,
#' marking primary loadings (\code{1}), penalised cells (\code{0}), and
#' free cells (\code{NA}).
#'
#' @param x An \code{esem_target} matrix from \code{\link{make_target}}.
#' @param ... Ignored.
#'
#' @return Invisibly returns \code{x}; called for the side effect of printing
#'   the target loading matrix.
#' @export
print.esem_target <- function(x, ...) {
  cat("Target loading matrix (", nrow(x), " items x ", ncol(x), " factors)\n",
      sep = "")
  cat("  1 = primary loading (targeted)\n")
  cat("  0 = penalised toward zero\n")
  cat(" NA = free (no penalty)\n\n")

  # Replace 0/1/NA with readable symbols for display
  display <- matrix(
    ifelse(is.na(x), " . ", ifelse(x == 1, " 1 ", " 0 ")),
    nrow = nrow(x),
    dimnames = dimnames(x)
  )
  print(display, quote = FALSE)
  invisible(x)
}


#' Modify an Existing Target Matrix
#'
#' Convenience function to set specific cells of a target matrix after
#' initial construction with \code{\link{make_target}}.
#'
#' @param target An \code{esem_target} matrix.
#' @param items Integer indices or character names of items to modify.
#' @param factors Integer indices or character names of factors to modify.
#' @param value New value: \code{1}, \code{0}, or \code{NA}.
#'
#' @return The modified target matrix.
#' @export
set_target <- function(target, items, factors, value) {
  if (!inherits(target, "esem_target"))
    stop("`target` must be an esem_target matrix from make_target().", call. = FALSE)

  target[items, factors] <- value
  target
}
