#' @rdname Func_Cross
#'
#' @useDynLib dgdecomp
#' @export
#'
Func_Cross_Matx <- function(vec_x, vec_y, size1, size2) {
  
  ## Simulate sequences of "positions" of the two vectors
  if (class(vec_x) == "matrix") {
    posit_seqs <- Func_Create_Combn(ncol(vec_x), size1, size2)
  } else {
    posit_seqs <- Func_Create_Combn(length(vec_x), size1, size2)
  }
  
  ## Apply combo positioning to the input data vectors
  #### NOTE that this is traversing across columns of the combo set
  #### NOT the number of factors
  return(
    .Call("ArmaCross", vec_x, vec_y, posit_seqs$vec_x_pos, posit_seqs$vec_y_pos)
  )
  
}