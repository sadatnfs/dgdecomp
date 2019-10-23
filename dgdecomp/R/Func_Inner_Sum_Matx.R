#' @rdname Func_Inner_Sum
#'
#' @export
#'
Func_Inner_Sum_Matx <- function(P, vec_x, vec_y) {
  P_upper <- ifelse(P %% 2 == 0, P / 2, (P + 1) * 0.5)

  sum_count <- rep(0, nrow(vec_x))

  for (r in 1:P_upper) {
    if (r > 1) {
      ## Create positional vectors
      size1 <- P - r ## r < P and so will never be >= P-1
      size2 <- r - 1 ## Will never be <= 1 because r > 1
      posit_seqs <- Func_Create_Combn(ncol(vec_x), size1, size2)
    } else {
      posit_seqs <- list(
        "vec_x_pos" = as.matrix(0),
        "vec_y_pos" = as.matrix(0)
      )
    }


    sum_count <- sum_count + Func_Inner_Matx(
      P, r, vec_x, vec_y,
      posit_seqs$vec_x_pos,
      posit_seqs$vec_y_pos
    )
  }
  return(sum_count)
}
