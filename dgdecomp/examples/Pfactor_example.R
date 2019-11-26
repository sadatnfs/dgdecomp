
### Simple tests on a P factor example

rm(list = ls())
# try(detach("package:dgdecomp"))
pacman::p_load(data.table, matrixStats, MASS, foreach, doParallel, dgdecomp)
number_of_factors <- 8

##### Test on simple 1-D data ----

## Simulate a listed 1-D dataset
ffdecomp_simdata <- simulate_decomp_data(number_of_factors)

## Compute decomp with our glorious functions
decomp_out <- Decomp_Factors_Matx(
  ffdecomp_simdata$vec_X_lag,
  ffdecomp_simdata$vec_X_today,
  tolerance = 1e-4
)
(all.equal(sum(decomp_out), ffdecomp_simdata$y_today - ffdecomp_simdata$y_lag))


##### Testing on data.table input ----

### Simulate 2 time periods; P factors; 3 groups
sim_dt <- simulate_decomp_data_fullmat(2, number_of_factors, 3)

decomp_out_DT <- Decomp_on_DT(
  input_data = sim_dt,
  factor_names = paste0("X_", c(1:number_of_factors)),
  bycol = "Id",
  time_col = "t"
)

true_delta <- sim_dt[,
  .(Ydelta = Y - shift(Y)),
  by = "Id"
][!is.na(Ydelta), .(Id, Ydelta)]
decomp_delta <- decomp_out_DT[, .(Id,
  decomp_delta = rowSums(.SD)
),
.SDcols = paste0("decomp_X_", 1:number_of_factors)
]

(all.equal(true_delta$Ydelta, decomp_delta$decomp_delta))


##### Testing on multiple time points input ----

### Simulate 3 time periods; P factors; 3 groups
sim_dt_multiv <- simulate_decomp_data_fullmat(3, number_of_factors, 3)

decomp_t1_t2 <- Decomp_on_DT(
  input_data = sim_dt_multiv[t %in% c(1,2)],
  factor_names = paste0("X_", c(1:number_of_factors)),
  bycol = "Id",
  time_col = "t"
)

decomp_t2_t3 <- Decomp_on_DT(
  input_data = sim_dt_multiv[t %in% c(2,3)],
  factor_names = paste0("X_", c(1:number_of_factors)),
  bycol = "Id",
  time_col = "t"
)

decomp_t1_t3 <- Decomp_on_DT(
  input_data = sim_dt_multiv[t %in% c(1,3)],
  factor_names = paste0("X_", c(1:number_of_factors)),
  bycol = "Id",
  time_col = "t"
)






