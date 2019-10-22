#' @rdname Decomp_Factors
#'
#' @export
#'
Decomp_Factors_Matx <- function(vec_x, vec_y, return_dt = TRUE, ...) {
  
  # Simple assertions
  # stopifnot(length(vec_x) == length(vec_y))
  
  # Coerce vectors into matrices
  if (class(vec_x) == "numeric") {
    vec_x <- t(as.matrix(vec_x))
    vec_y <- t(as.matrix(vec_y))
  } 
  
  
  # Gather the number of factors
  num_factors <- ncol(vec_x)
  
  # Compute each marginal effect
  effects_all <- sapply(
    c(1:num_factors),
    function(x) Func_Inner_Sum_Matx(
      P = num_factors,
      vec_x = vec_x[,-1 * x],
      vec_y = vec_y[,-1 * x]
    ) * (vec_y[,x] - vec_x[,x])
  )
  
  if (ncol(effects_all) != num_factors) {
    effects_all <- t(effects_all)
  }
  
  if (return_dt) {
    effects_all <- data.table(effects_all)
  }
  
  # Assertion on whether the decomp actually worked
  stopifnot(base::all.equal(
    apply(vec_y, 1, prod) - apply(vec_x, 1, prod),
    apply(effects_all, 1, sum),
    ...
  ))
  
  # Return the effects
  return(effects_all)
}
