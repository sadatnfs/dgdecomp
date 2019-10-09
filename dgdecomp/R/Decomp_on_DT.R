#' @title Apply DG Decomposition to data.table columns
#'
#' @param input_data A data.table with the factors, which must already
#' be sorted
#'
#' @param factor_names A vector of column names for the factor
#'
#' @param bycol The 'by' slicer which must make sure that the data is
#' reduced to just 2 rows per group after slicing
#'
#' @return A data.table of the same size as input, but instead with
#' the additive decomposition results (first row will be NA as being the
#' starting period)
#'
#' @importFrom data.table data.table rbindlist shift .SD setnames
#' @export
#'
Decomp_on_DT <- function(input_data, factor_names, bycol) {

  ## Apply decomp to each row of data
  decomp_output <- input_data[,
    as.list(Decomp_Factors(as.matrix(.SD)[1, ], as.matrix(.SD)[2, ])),
    by = bycol, .SDcols = factor_names
  ]

  data.table::setnames(
    decomp_output,
    c(bycol, paste0("decomp_", factor_names))
  )

  return(decomp_output)
}
