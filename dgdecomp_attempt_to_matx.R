


rm(list = ls())
detach("package:dgdecomp")
pacman::p_load(data.table, matrixStats, MASS, foreach, 
               doParallel, dgdecomp, Rcpp, RcppArmadillo, RcppAlgos)
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

###### DOES NOT HAVE THE COMBN YET
dgdecomp:::Func_Cross_Matx(vec_x, vec_y, size1, size2)



try(dyn.unload("~/Documents/repos/dgdecomp_repo/dgdecomp/src/armacross.so"))
dyn.load("~/Documents/repos/dgdecomp_repo/dgdecomp/src/armacross.so")


lapply(c(1:3), function(xx) Func_Inner_Sum(number_of_factors, vec_x[xx,], vec_y[xx,]))



q("no")
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















