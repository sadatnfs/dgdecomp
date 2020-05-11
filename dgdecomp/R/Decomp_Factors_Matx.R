#' @rdname Decomp_Factors
#'
#' @export
#'
Decomp_Factors_Matx <- function(mat_x, mat_y, return_dt = TRUE, use_cpp = TRUE,
                                parallel = 1, cpplib = 'arma',
                                equality_check = TRUE,
                                ...) {

  # Coerce vectors into matrices
  if (!("matrix" %in% class(mat_x) )) {
    mat_x <- t(as.matrix(mat_x))
    mat_y <- t(as.matrix(mat_y))
  }


  # Gather the number of factors
  num_factors <- ncol(mat_x)


  # Redefine parallel if we're overallocating threads
  # that is, if threads > num_factors
  parallel <- min(num_factors, parallel)

  # Compute each marginal effect
  if (use_cpp) {
    if (cpplib == "arma") {
      effects_all <- .Call("ArmaDFInnerLoop", num_factors, mat_x, mat_y, parallel)
    } else if (cpplib == "eigen") {
      effects_all <- .Call("EigDFInnerLoop", num_factors, mat_x, mat_y, parallel)
    } else {
      stop("Need to specify either arma or eigen for cpplig")
    }
  } else {
    effects_all <- sapply(
      c(1:num_factors),
      function(x) {
        input_x <- mat_x[, -1 * x]
        input_y <- mat_y[, -1 * x]

        ## Edge case for inputting vectors
        if (class(input_x) != "matrix") {
          input_x <- t(as.matrix(input_x))
          input_y <- t(as.matrix(input_y))
        }

        Func_Inner_Sum_Matx(
          P = num_factors,
          vec_x = input_x,
          vec_y = input_y
        ) * as.matrix(mat_y[, x] - mat_x[, x])
      }
    )
  }


  ## Make sure that the output is a matrix
  if ("numeric" %in% class(effects_all)) {
    effects_all <- t(as.matrix(effects_all))
  }

  ## Fix in case the matrix is flipped
  if (ncol(effects_all) != num_factors) {
    effects_all <- t(effects_all)
  }

  if (return_dt) {
    effects_all <- data.table(effects_all)
  }

  # Assertion on whether the decomp actually worked
  if (equality_check) {
    stopifnot(base::all.equal(
      apply(mat_y, 1, prod) - apply(mat_x, 1, prod),
      apply(effects_all, 1, sum),
      ...
    ))
  }

  # Return the effects
  return(effects_all)
}
