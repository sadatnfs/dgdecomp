#include <iostream>
// #include <Rcpp.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace std;
using namespace arma;


/// Compile with
//  PKG_CXXFLAGS="-I/home/nasadat/R/x86_64-pc-linux-gnu-library/3.6/RcppArmadillo/include -I/home/nasadat/R/x86_64-pc-linux-gnu-library/3.6/Rcpp/include" PKG_LIBS=`Rscript -e 'Rcpp:::LdFlags()'` R CMD SHLIB ~/Documents/repos/dgdecomp_repo/dgdecomp/src/decomp_funks.cpp


////////// We can write a C++ function as is
////////// But then we must pass it through SEXP
////////// if we are to call it from R

///// Here, I wrote a mostly-pure C++ function
///// and then an export function with SEXP inputs and outputs

//// Output the positioned matrix (Func_Cross)
arma::mat ArmaCross_(
    arma::mat vec_x, arma::mat vec_y,
    arma::mat vec_x_pos, arma::mat vec_y_pos){
  
  
  // Create output matrix
  arma::mat prod_inner(vec_x.n_rows, vec_x_pos.n_cols);
  
  // Temp store
  double tmpprod = 1.;
  
  for( int row1 = 0; row1 < prod_inner.n_rows; row1++){
    
    for( int col1 = 0; col1 < vec_x_pos.n_cols; col1++){ // invariant between vec_x_pos and vec_y_pos
      
      // Go over each row of vec_x_pos
      for( int row2 = 0; row2 < vec_x_pos.n_rows; row2++){ // vec_x_pos specific
        tmpprod = tmpprod * vec_x(row1, vec_x_pos(row2, col1)-1); // -1 correction for zero index
      }

      // Go over each row of vec_y_pos
      for( int row3 = 0; row3 < vec_y_pos.n_rows; row3++){ // vec_x_pos specific
        tmpprod = tmpprod * vec_y(row1, vec_y_pos(row3, col1)-1); // -1 correction for zero index
      }
      
      // Record in prod_inner
      prod_inner(row1, col1) = tmpprod;
      tmpprod = 1.;
      
    }
    
  }
    
  return prod_inner;
  
}

// Create numerator (Func_Num)
arma::mat ArmaNum_(
    int P, int r,
    arma::mat vec_x, arma::mat vec_y){
  
  if (r == 1) {
  } else {
    
  }
  
}

//// The returner
// [[Rcpp::Export]]
RcppExport SEXP ArmaCross( 
    SEXP vec_x_, SEXP vec_y_,
    SEXP vec_x_pos_, SEXP vec_y_pos_){
  
  // Must coerce vector to arma::vec type (or whatever type to feed into funk)
  // so that we convert from input type from R (SEXP) to whatever our funk needs
  
  return Rcpp::wrap(
    ArmaCross_(
      as<arma::mat>(vec_x_),
      as<arma::mat>(vec_y_),
      as<arma::mat>(vec_x_pos_),
      as<arma::mat>(vec_y_pos_)
    ));
  
}

