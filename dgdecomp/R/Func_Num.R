#' @title Compute the numerator of the inner sums in the DG Decomposition
#'
#' @param P Number of factors
#'
#' @param r The summing indicator
#'
#' @param vec_x First input vector
#'
#' @param vec_y Second input vector
#'
#' @return A single numeric from the sums of Func_Cross()
#'
#' @export
#'
Func_Num <- function(P, r, vec_x, vec_y) {
  if (r == 1) {
    this_count <- prod(vec_x) + prod(vec_y)
  } else {
    
    ## Pick : P-r from SMALL and mult with r-1 from CAPS
    ## with : P-r from CAPS and mult with r-1 from SMALL
    
    size1 <- P - r ## r < P and so will never be >= P-1
    size2 <- r - 1 ## Will never be <= 1 because r > 1
    
    ## First, P-r small and r-1 caps:
    this_count <- sum(Func_Cross(vec_x, vec_y, size1, size2))
    
    ## Next, P-r caps and r-1 small:
    ### ONLY applicable if we are not comparing identical sizes
    if (size1 != size2) {
      this_count <- this_count + sum(
        Func_Cross(vec_x, vec_y, size2, size1))
    }
  }
  
  return(this_count)
}
