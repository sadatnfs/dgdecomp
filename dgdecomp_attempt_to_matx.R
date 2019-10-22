


rm(list = ls())
detach("package:dgdecomp")
pacman::p_load(data.table, matrixStats, MASS, foreach, 
               doParallel, dgdecomp, Rcpp, RcppArmadillo)
number_of_factors <- 7



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
vec_x <- x_lag[, -1 * marginal_factor]
vec_y <- x_today[, -1 * marginal_factor]

r <- 3
size1 <- number_of_factors - r ## r < P and so will never be >= P-1
size2 <- r - 1 ## Will never be <= 1 because r > 1
Pfac <- number_of_factors - 1


##### Func_Cross() ####

## Simulate sequences of "positions" of the two vectors

posit_seqs <- Func_Create_Combn(eval(number_of_factors - 1), size1, size2)

## Apply combo positioning to the input data vectors
#### NOTE that this is traversing across columns of the combo set
#### NOT the number of factors
# vapply(c(1:ncol(posit_seqs$vec_x_pos)),
#   function(x) {
#     prod(c(
#       vec_x[, posit_seqs$vec_x_pos[, x]],
#       vec_y[, posit_seqs$vec_y_pos[, x]]
#     ))
#   },
#   FUN.VALUE = 0
# )
# 
## Expected output ##
t(sapply(c(1:3), function(roo) {
  vapply(c(1:ncol(posit_seqs$vec_x_pos)),
    function(x) {
      prod(c(
        vec_x[roo, ][posit_seqs$vec_x_pos[, x]],
        vec_y[roo, ][posit_seqs$vec_y_pos[, x]]
      ))
    },
    FUN.VALUE = 0
  )
}))

try(dyn.unload("~/Documents/repos/dgdecomp_repo/dgdecomp/src/decomp_funks.so"))
dyn.load("~/Documents/repos/dgdecomp_repo/dgdecomp/src/decomp_funks.so")
is.loaded("ArmaCross")

.Call("ArmaCross", vec_x, vec_y, posit_seqs$vec_x_pos, posit_seqs$vec_y_pos)

### CROSS WORKS??
rbind(Func_Cross(vec_x[1,], vec_y[1,], size1, size2),
      Func_Cross(vec_x[2,], vec_y[2,], size1, size2),
      Func_Cross(vec_x[3,], vec_y[3,], size1, size2))




run_decomp_sim <- function(Time = 2, P, G) {
  sim_dt <- simulate_decomp_data_fullmat(Time, P, G)

  decomp_out_DT <- Decomp_on_DT(
    input_data = sim_dt,
    factor_names = paste0("X_", c(1:P)),
    bycol = "Id"
  )

  true_delta <- sim_dt[,
    .(Ydelta = Y - shift(Y)),
    by = "Id"
  ][!is.na(Ydelta), .(Id, Ydelta)]
  decomp_delta <- decomp_out_DT[, .(Id,
    decomp_delta = rowSums(.SD)
  ),
  .SDcols = paste0("decomp_X_", 1:P)
  ]

  return(decomp_delta)
}




### Test grid of outputs
### when we serialize everything
log_table <- rbindlist(
  lapply(c(5), function(grou) {
    rbindlist(
      lapply(c(5), function(facto) {
      tmp <- system.time(run_decomp_sim(P = facto, G = grou))[3]
      return(data.table(
        "factors" = facto,
        "groups" = grou,
        "time_elapsed" = tmp
      ))
    })
    )
  })
)















