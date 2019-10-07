#' @title Compute the combination of all the elements of the given vectors
#' corresponding to the given sizes
#'
#' @param vec_x First input vector
#'
#' @param vec_y Second input vector
#'
#' @param size1 Number of elements to take from vec_x
#'
#' @param size2 Number of elements to take from vec_y
#'
#' @return A vector of column products made from the unique combinations
#' based on the params
#'
#' @importFrom utils combn
#' @importFrom assertthat are_equal
#'
#' @export
#'
Func_Cross <- function(vec_x, vec_y, size1, size2) {

  ## First, P-r small and r-1 caps:
  vec_x_combo <- utils::combn(vec_x, m = size1)
  vec_y_combo <- utils::combn(vec_y, m = size2)

  ### This needs reversal to counter the count
  if (size2 == 1) {
    vec_y_combo <- t(as.matrix(vec_y_combo[, ncol(vec_y_combo):1]))
  } else {
    vec_y_combo <- vec_y_combo[, ncol(vec_y_combo):1]
  }

  stopifnot(assertthat::are_equal(ncol(vec_x_combo), ncol(vec_y_combo)))

  ## Now that we have the matrix of combos,
  ## we multiple X and Y combo column wise
  return(vapply(c(1:ncol(vec_x_combo)),
    function(x) {
      prod(c(vec_x_combo[, x], vec_y_combo[, x]))
    },
    FUN.VALUE = 0
  ))
}
