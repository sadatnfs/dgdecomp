#' @title Compute the combination of all the elements of the given vectors
#' corresponding to the given sizes (using Func_Create_Combn)
#'
#' @param Pfac Number of factors minus 1
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
Func_Cross <- function(Pfac, vec_x, vec_y, size1, size2) {
  
  ## Simulate sequences of "positions" of the two vectors
  posit_seqs <- Func_Create_Combn(Pfac, size1, size2)
  
  ## Apply combo positioning to the input data vectors
  return(vapply(c(1:ncol(vec_x_sim)),
         function(x) {
           prod(c(vec_x[vec_x_sim[,x]], vec_y[vec_y_sim[,x]] ))
         },
         FUN.VALUE = 0
  ))
  
  
}
