#' @title Compute the marginal decomposition effects from given input
#'
#' @description \code{Decomp_Factor_Matx()} can be used if the input has
#' multiple rows of data to decompose, whereas \code{Decomp_Factor()} only
#' takes in vectors as inputs.
#'
#' @param vec_x First input vector (represents t-1)
#'
#' @param vec_y Second input vector (represents t)
#'
#' @param mat_x First input matrix (only for \code{Decomp_Factor_Matx()})
#'
#' @param mat_y Second input matrix (only for \code{Decomp_Factor_Matx()})
#'
#' @param return_dt A boolean on whether to return a data.table or a vector
#' 
#' @param use_cpp A boolean on whether to use the C++ compiled code or not.
#' Default: TRUE
#'
#' @param ... extra parameters to be passed to \code{all.equal()}, for e.g.
#' the tolerance.
#'
#' @return A data.table or vector of decomposed effects for each factors,
#' which is already multiplied by the change values
#'
#' @importFrom data.table data.table
#'
#' @export
#'
Decomp_Factors <- compiler::cmpfun(function(vec_x, vec_y, return_dt = TRUE, ...) {

  # Simple assertions
  # stopifnot(length(vec_x) == length(vec_y))

  # Gather the number of factors
  num_factors <- length(vec_x)

  # Compute each marginal effect
  effects_all <- vapply(
    c(1:num_factors),
    function(x) Func_Inner_Sum(
        num_factors,
        vec_x[-1 * x],
        vec_y[-1 * x]
      ) * (vec_y[x] - vec_x[x]),
    FUN.VALUE = 0
  )

  if (return_dt) {
    effects_all <- data.table(t(effects_all))
  }

  # Assertion on whether the decomp actually worked
  stopifnot(base::all.equal(
    prod(vec_y) - prod(vec_x),
    sum(effects_all),
    ...
  ))

  # Return the effects
  return(effects_all)
})
