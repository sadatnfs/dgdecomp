#' @title Simulate grouped data for decomp analysis (T by IDI by factors)
#'
#' @param T_term Number of time periods
#'
#' @param num_factors Number of factors (gets slower for large numbers)
#'
#' @param id_grps Number of groups
#'
#' @return A data.table with T_term*id_grps rows and num_factors+1 columns
#' where X_1, ... , X_p are the factors, and Y is the rowwise product of the
#' factors
#'
#' @importFrom data.table data.table setnames setkeyv shift
#' 
#' @importFrom MASS mvrnorm
#' 
#' @importFrom matrixStats rowProds
#'
#' @export
#'
simulate_decomp_data_fullmat <- function(T_term, num_factors, id_grps) {

  sim_dt <- data.table(expand.grid(t = c(1:T_term), Id = c(1:id_grps)))

  factor_sims <- data.table(
    matrix(mvrnorm(
      T_term * num_factors * id_grps,
      c(-1, 1),
      matrix(c(10, 3, 3, 2), 2, 2)
    ),
    nrow = T_term * id_grps, ncol = num_factors
    ) / 1
  )
  setnames(factor_sims, paste0("X_", c(1:num_factors)))
  sim_dt <- cbind(sim_dt, factor_sims)
  sim_dt[, Y := rowProds(as.matrix(.SD)),
    .SDcols = paste0("X_", c(1:num_factors))
  ]
  setkeyv(sim_dt, c("Id", "t"))

  return(sim_dt)
}
