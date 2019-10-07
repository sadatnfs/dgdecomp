#' @title Compute the fraction for the inner sum in the DG Decomposition
#' where all the params gets passed to Func_Num() and Func_Den()
#'
#' @param P Number of factors
#'
#' @param r The summing indicator
#'
#' @param vec_x First input vector
#'
#' @param vec_y Second input vector
#'
#' @return The fraction of the results of Func_Num() and Func_Den
#'
#' @export
#'
Func_Inner <- function(P, r, vec_x, vec_y) {
  Func_Num(P, r, vec_x, vec_y) / Func_Dem(P, r)
}
