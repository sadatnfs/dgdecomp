#' @title Compute the combination of all the elements of the given vectors
#' corresponding to the given sizes (using Func_Create_Combn)
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
#' of the *data*
#'
#'
#' @export
#'
Func_Cross <- compiler::cmpfun(function(vec_x, vec_y, size1, size2) {

  ## Simulate sequences of "positions" of the two vectors
  posit_seqs <- Func_Create_Combn(length(vec_x), size1, size2)

  ## Apply combo positioning to the input data vectors
  #### NOTE that this is traversing across columns of the combo set
  #### NOT the number of factors
  return(
    vapply(c(1:ncol(posit_seqs$vec_x_pos)),
      function(x) {
        prod(c(
          vec_x[posit_seqs$vec_x_pos[, x]],
          vec_y[posit_seqs$vec_y_pos[, x]]
        ))
      },
      FUN.VALUE = 0
    )
  )
})
