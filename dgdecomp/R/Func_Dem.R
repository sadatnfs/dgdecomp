#' @title Compute the denomiator of the inner sums in the DG Decomposition
#'
#' @param P Number of factors
#'
#' @param r The summing indicator
#'
#' @return A numeric with value of \code{P * choose((P - 1), (r - 1))}
#'
#' @export
#'
Func_Dem <- function(P, r) {
  P * choose((P - 1), (r - 1))
}
