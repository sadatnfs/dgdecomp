#' @title Sum over all inner sums for the DG Decomposition
#'
#' @param P Number of factors
#'
#' @param vec_x First input vector
#'
#' @param vec_y Second input vector
#'
#' @return A numeric value with the full inner sum for the given effect
#'
#' @export
#'
Func_Inner_Sum <- compiler::cmpfun(function(P, vec_x, vec_y) {
  sum_count <- 0
  P_upper <- ifelse(P %% 2 == 0, P / 2, (P + 1) * 0.5)
  for (r in 1:P_upper) {
    sum_count <- sum_count + Func_Inner(P, r, vec_x, vec_y)
  }
  return(sum_count)
})
