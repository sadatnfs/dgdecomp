

#### Just some testing sandbox code

rm(list = ls())
detach("package:dgdecomp")
pacman::p_load(data.table, dgdecomp, Rcpp, RcppArmadillo, microbenchmark)
number_of_factors <- 7



dyn.load("/opt/compiled_code_for_R/setthreads.so")
.C("setOMPthreads", 10)
.C("setMKLthreads", 10)

### Simulate 10 time periods; P factors; 3 groups
sim_dt <- simulate_decomp_data_fullmat(2, number_of_factors, 3)



### Let's prep the data to feed into Func_Cross()
factor_names <- paste0("X_", 1:number_of_factors)
time_name <- "t"
group_name <- "Id"
output_name <- "Y"

## Get marginal effect of one factor
## (and so we deal with P-1 values in the vector)
marginal_factor <- 1
x_lag <- as.matrix(sim_dt[t == 1, .SD, .SDcols = factor_names])
x_today <- as.matrix(sim_dt[t == 2, .SD, .SDcols = factor_names])

## Remove the Pth colummn
# vec_x <- x_lag[, -1 * marginal_factor]
# vec_y <- x_today[, -1 * marginal_factor]

r <- 3
size1 <- number_of_factors - r ## r < P and so will never be >= P-1
size2 <- r - 1 ## Will never be <= 1 because r > 1
Pfac <- number_of_factors - 1



run_decomp_sim <- function(Time = 2, P, G, use_cpp) {
  sim_dt <- simulate_decomp_data_fullmat(Time, P, G)

  # ## Prep lag and today matrices
  # lag_mat <- as.matrix(
  #   sim_dt[t == 1, .SD, .SDcols = paste0("X_", c(1:P))]
  #   )
  # today_mat <- as.matrix(
  #   sim_dt[t == 2, .SD, .SDcols = paste0("X_", c(1:P))]
  #   )
  #
  # decomp_out_DT <- Decomp_Factors_Matx(
  #   lag_mat,
  #   today_mat,
  #   tolerance = 1e-5
  # )
  #
  # data.table::setnames(
  #   decomp_out_DT,
  #   paste0("decomp_", paste0("X_", c(1:P)))
  # )

  ## Use the data.table method for doing decomp
  #### NOTE that there's an overhead of computing the lag columns
  decomp_out_DT <- Decomp_on_DT(
    input_data = sim_dt,
    factor_names = paste0("X_", c(1:P)),
    bycol = c("Id"),
    time_col = "t",
    use_cpp = use_cpp
  )


  true_delta <- sim_dt[,
    .(Ydelta = Y - shift(Y)),
    by = "Id"
  ][!is.na(Ydelta), .(Id, Ydelta)]
  decomp_delta <- decomp_out_DT[, .(decomp_delta = rowSums(.SD)),
    .SDcols = paste0("decomp_X_", 1:P)
  ]

  return(decomp_delta)
}


run_sim <- function(cppuse) {
  for (grou in c(5000, 10000, 15000)) {
    for (facto in 5:10) {
      tmp <- system.time(
        run_decomp_sim(P = facto, G = grou, use_cpp = cppuse))[3]
      print(paste0(
        "factors = ", facto,
        " groups = ", grou,
        " time_elapsed = ", format(tmp, digits = 4)
      ))
    }
  }
}



mbm <- microbenchmark(
  withcpp = run_sim(TRUE),
  without_cpp = run_sim(FALSE),
  times=5
)
mbm











