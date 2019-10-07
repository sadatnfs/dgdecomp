### P factor decomp formula: ----
## Let R_1 = A_1 * A_2 * ... A_p
## Let R_2 = a_1 * a_2 * ... a_p
### Then R_2 - R_1 = (alpha_1 effect) +
### (alpha_2 effect) + ... + (alpha_p effect)
### where:
## alpha_1 = Q * (a_1 - A_1)
### Q = Q(a_2, a_3, ..., a_p, A_2, A_3, ... , A_p)
###   = sum_{r=1)^{P} [ sum of all  (P-1) letter terms with
###          (P-r) small letters and (r-1) cap letter product +
###          (P-r) cap letters and (r-1) small letter product]
### S = P/2 if P is even; (P+1)/2 if P is odd
### Scheme:
### Let the inside sum be a function Inner():
### Inner(P,r) = Num(P,r) / Den(P,r)
### Den(P,r) = P*(choose(P-1, r-1))
### Num(P,r) = Cross(r-1 caps, P-r smalls) + \
###            Cross(P-r caps, r-1 smalls)
### Cross([x_i],[y_j]) = sum_{i=1}^{N1} sum_{j=1}^{N2} x_i * y_j
### for [x_i] = x_1, ... , x_N1;
###     [y_j] = y_1, ... , y_N2;
#' @title Simulate simple random decomp data for P factors such that the
#' product of P factors equal a measure for 2 time periods
#'
#' @param num_fac Number of factors to simulate
#'
#' @return A named list with the vector of P factors for 2 time periods,
#' and 2 numeric measures for each time period, which are just the
#' product of each of the two factor vectors
#'
#' @importFrom stats rnorm
#'
#' @export
#'
simulate_decomp_data <- function(num_fac) {
  vec_X_lag <- rnorm(n = num_fac) * 10
  vec_X_today <- (1:num_fac)
  y_lag <- prod(vec_X_lag)
  y_today <- prod(vec_X_today)

  return(list(
    "vec_X_lag" = vec_X_lag,
    "vec_X_today" = vec_X_today,
    "y_lag" = y_lag,
    "y_today" = y_today
  ))
}
