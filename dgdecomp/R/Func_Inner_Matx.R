#' @rdname Func_Inner
#' 
#' @export
#'
Func_Inner_Matx <- function(P, r, vec_x, vec_y) {
  posit_seqs <- Func_Create_Combn(ncol(vec_x), size1, size2)
  return(
    .Call("ArmaInner", number_of_factors, r, vec_x, 
          vec_y, posit_seqs$vec_x_pos, posit_seqs$vec_y_pos)
  )
}