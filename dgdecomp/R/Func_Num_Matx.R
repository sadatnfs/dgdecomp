#' @rdname Func_Num
#'
#' @export
#'
Func_Num_Matx <- function(P, r, vec_x, vec_y) {


  ## Create positional vectors
  size1 <- P - r ## r < P and so will never be >= P-1
  size2 <- r - 1 ## Will never be <= 1 because r > 1
  posit_seqs <- Func_Create_Combn(ncol(vec_x), size1, size2)
  return(
    .Call(
      "ArmaNum", number_of_factors, r, vec_x,
      vec_y, posit_seqs$vec_x_pos, posit_seqs$vec_y_pos
    )
  )
}
