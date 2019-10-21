#' @title Compute the combination of all positions of the given factor
#' segmented into two pieces
#'
#' @param Pfac Number of factors minus 1
#'
#' @param size1 Number of elements to take from vec_x
#'
#' @param size2 Number of elements to take from vec_y
#'
#' @return A vector of positions made from the unique combinations
#' of size1 and size2
#'
#' @importFrom utils combn
#' @importFrom assertthat are_equal
#'
#' @export
#'
Func_Create_Combn <- function(Pfac, size1, size2) {
  
  ## Make sure that that the sizes are consistent with input
  stopifnot(assertthat::are_equal(
    Pfac, 
    eval(size1 + size2)))
  
  ## Simulate sequences and combn off of those
  vec_x_sim <- utils::combn(c(1:Pfac), m = size1)
  vec_y_sim <- utils::combn(c(1:Pfac), m = size2)
  
  ## And because we will have the same number of factors in the matrix,
  ## we can use the same ordering of the sequences to project on to the matrix!
  
  ## Reverse the combo order of y
  if (size2 == 1) {
    vec_y_sim <- t(as.matrix(vec_y_sim[, ncol(vec_y_sim):1]))
  } else {
    vec_y_sim <- vec_y_sim[, ncol(vec_y_sim):1]
  }  
  
  ## Make sure output dimensions are good
  stopifnot(assertthat::are_equal(
    ncol(vec_x_pos), 
    ncol(vec_y_pos)))
  
  ## Return the named list of position vectors
  return(list("vec_x_pos" = vec_x_sim,
              "vec_y_pos" = vec_y_sim))
  
}
