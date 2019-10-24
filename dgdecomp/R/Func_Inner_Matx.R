#' @rdname Func_Inner
#'
#' @export
#'
Func_Inner_Matx <- function(P, r, vec_x, vec_y) {
  return(
    .Call(
      "ArmaInner", P, r, vec_x, vec_y
    )
  )
}
