#' @title Apply DG Decomposition to data.table columns
#'
#' @param input_data A data.table with the factors, which must already
#' be sorted
#'
#' @param factor_names A vector of column names for the factor
#'
#' @param time_col A string for the column name
#'
#' @param bycol The 'by' slicer which must make sure that the data is
#' reduced to just 2 rows per group after slicing
#' 
#' @param use_cpp A boolean on whether to use the C++ compiled code for the
#' factor for-loop or not (passes to Decomp_Factor_Matx()). Default: TRUE
#' 
#' @param ... extra parameters to be passed through \code{Decomp_Factors()}
#' to \code{all.equal()}
#'
#' @return A data.table of the same size as input, but instead with
#' the additive decomposition results (first row will be NA as being the
#' starting period)
#'
#' @importFrom data.table data.table rbindlist shift .SD setnames
#' @export
#'
Decomp_on_DT <- compiler::cmpfun(function(input_data,
                                          factor_names,
                                          time_col,
                                          bycol, 
                                          use_cpp = TRUE,
                                          ...) {

  ## Create the lag and current matrices from the DT
  lag_mat <- input_data[,
    as.list(as.matrix(.SD)[1, ]),
    .SDcols = factor_names, by = bycol
  ]
  curr_mat <- input_data[,
    as.list(as.matrix(.SD)[2, ]),
    .SDcols = factor_names, by = bycol
  ]

  ## Make sure that the two datasets' dimensions line up
  stopifnot((ncol(lag_mat) == ncol(curr_mat)) & (nrow(lag_mat) == nrow(curr_mat)))


  ## Apply decomp to with each of these new matrices
  decomp_output <- Decomp_Factors_Matx(
    mat_x = as.matrix(lag_mat[, .SD, .SDcols = factor_names]),
    mat_y = as.matrix(curr_mat[, .SD, .SDcols = factor_names]),
    use_cpp = use_cpp,
    ...
  )

  ## cbind the bycols to this output
  decomp_output <- cbind(lag_mat[, .SD, .SDcols = bycol], decomp_output)

  data.table::setnames(
    decomp_output,
    c(bycol, paste0("decomp_", factor_names))
  )

  return(decomp_output)
})
