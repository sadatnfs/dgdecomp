


rm(list = ls())
pacman::p_load(data.table, dgdecomp, Rcpp, RcppArmadillo, microbenchmark)


## A simple function to simulate a dataset and run decomp
run_decomp_sim <- function(Time = 2, P, G, threads) {
  sim_dt <- simulate_decomp_data_fullmat(Time, P, G)

  ## Use the data.table method for doing decomp
  #### NOTE that there's an overhead of computing the lag columns
  decomp_out_DT <- Decomp_on_DT(
    input_data = sim_dt,
    factor_names = paste0("X_", c(1:P)),
    bycol = c("Id"),
    time_col = "t",
    use_cpp = TRUE,
    parallel = threads
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


## Drive for the above function with varying factors and groups
run_sim <- function(factoC = c(5, 10, 15),
                    groupo = c(2500, 5000, 10000, 20000),
                    threads = 1) {
  for (grou in groupo) {
    for (facto in factoC) {
      tmp <- system.time(
        run_decomp_sim(P = facto, G = grou, threads = threads)
      )[3]
      print(paste0(
        "factors = ", facto,
        " groups = ", grou,
        " time_elapsed = ", format(tmp, digits = 4)
      ))
    }
  }
}



## Run with single thread
run_sim(threads = 1, factoC = c(5, 10))

## Run with 10 threads
run_sim(threads = 10)
