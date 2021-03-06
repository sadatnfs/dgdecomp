
### Simple tests on a P factor example

rm(list = ls())
# try(detach("package:dgdecomp"))
pacman::p_load(data.table, matrixStats, MASS, foreach, doParallel, dgdecomp)
number_of_factors <- 8

##### Test on simple 1-D data ----

## Simulate a listed 1-D dataset
ffdecomp_simdata <- simulate_decomp_data(number_of_factors)

## Compute decomp with our armadillo functions
decomp_out_arma <- Decomp_Factors_Matx(
  ffdecomp_simdata$vec_X_lag,
  ffdecomp_simdata$vec_X_today,
  tolerance = 1e-4,
  cpplib = "arma"
)
(all.equal(sum(decomp_out_arma), ffdecomp_simdata$y_today - ffdecomp_simdata$y_lag))

decomp_out_eigen <- Decomp_Factors_Matx(
  ffdecomp_simdata$vec_X_lag,
  ffdecomp_simdata$vec_X_today,
  tolerance = 1e-4,
  cpplib = "eigen"
)
(all.equal(sum(decomp_out_eigen), ffdecomp_simdata$y_today - ffdecomp_simdata$y_lag))

decomp_out_blaze <- Decomp_Factors_Matx(
  ffdecomp_simdata$vec_X_lag,
  ffdecomp_simdata$vec_X_today,
  tolerance = 1e-4,
  use_cpp = F,
  cpplib = "blaze"
)
(all.equal(sum(decomp_out_blaze), ffdecomp_simdata$y_today - ffdecomp_simdata$y_lag))





##### Testing on data.table input ----

### Simulate 2 time periods; P factors; 3 groups
number_of_factors <- 15
sim_dt <- dgdecomp::simulate_decomp_data_fullmat(2, number_of_factors, 10)

## Create the lag and current matrices from the DT
lag_mat <- sim_dt[ ,
                      as.list(as.matrix(.SD)[1, ]),
                      .SDcols = paste0("X_", c(1:number_of_factors)), by = "Id"
                      ][Id %in% c(1:3)]
curr_mat <- sim_dt[ ,
                       as.list(as.matrix(.SD)[2, ]),
                       .SDcols = paste0("X_", c(1:number_of_factors)), by = "Id"
                       ][Id %in% c(1:3)]

## Apply decomp to with each of these new matrices
decomp_output_TestArma <- Decomp_Factors_Matx(
  mat_x = as.matrix(lag_mat[, .SD, .SDcols = paste0("X_", c(1:number_of_factors))]),
  mat_y = as.matrix(curr_mat[, .SD, .SDcols = paste0("X_", c(1:number_of_factors))]),
  use_cpp = TRUE,
  parallel = 1,
  cpplib = "arma"
)
decomp_output_TestBlaze <- Decomp_Factors_Matx(
  mat_x = as.matrix(lag_mat[, .SD, .SDcols = paste0("X_", c(1:number_of_factors))]),
  mat_y = as.matrix(curr_mat[, .SD, .SDcols = paste0("X_", c(1:number_of_factors))]),
  use_cpp = TRUE,
  parallel = 1,
  cpplib = "blaze", 
  equality_check = T
)



system.time(
decomp_out_DT <- dgdecomp::Decomp_on_DT(
  input_data = sim_dt,
  factor_names = paste0("X_", c(1:number_of_factors)),
  bycol = "Id",
  time_col = "t",
  parallel = 1,
  cpplib = "arma"
))


system.time(
decomp_out_DT_true_blaze <- dgdecomp::Decomp_on_DT(
  input_data = sim_dt,
  factor_names = paste0("X_", c(1:number_of_factors)),
  bycol = "Id",
  time_col = "t",
  parallel = 1,
  cpplib = "blaze",
  equality_check = F
))

decomp_out_DT_false_arma <- Decomp_on_DT(
  input_data = sim_dt,
  factor_names = paste0("X_", c(1:number_of_factors)),
  bycol = "Id",
  time_col = "t",
  use_cpp = FALSE,
  cpplib = 'arma'
)
decomp_out_DT_false_eigen <- Decomp_on_DT(
  input_data = sim_dt,
  factor_names = paste0("X_", c(1:number_of_factors)),
  bycol = "Id",
  time_col = "t",
  use_cpp = FALSE,
  cpplib = 'eigen'
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






