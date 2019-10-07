#' @title Apply DG Decomposition to data.table columns
#'
#' @param input_data A data.table with the factors, which must already
#' be sorted
#'
#' @param factor_names A vector of column names for the factor
#'
#' @return A data.table of the same size as input, but instead with
#' the additive decomposition results (first row will be NA as being the
#' starting period)
#'
#' @importFrom data.table data.table rbindlist shift .SD setnames
#' @export
#'
Decomp_on_DT <- function(input_data, factor_names) {
  dt_out <- rbindlist(

    ## Apply decomp to each row of data
    lapply(c(1:nrow(input_data)), function(x)
      Decomp_Factors(

        ## First input: value of X yesterday
        unlist(input_data[, data.table::shift(.SD),
          .SDcols = factor_names
        ][x, ]),

        ## Second input: value of X today
        unlist(input_data[, .SD, .SDcols = factor_names][x, ])
      ))
  )
  data.table::setnames(dt_out, paste0("decomp_", factor_names))
  return(dt_out)
}
