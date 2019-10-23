#' @rdname Func_Inner
#'
#' @export
#'
Func_Inner_Matx <- function(P, r, vec_x, vec_y, vec_x_pos, vec_y_pos) {
  return(
    .Call(
      "ArmaNum", P, r, vec_x,
      vec_y, vec_x_pos, vec_y_pos
    ) / Func_Dem(P, r)
  )
}
