#include <iostream>
// #include "mkl.h"
#include <omp.h>

// #include <iterator>
#include <vector>
#include <algorithm>
#include <cstdlib>
// #include <discreture.hpp>
// #include <Rcpp.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace std;
using namespace arma;


/// Compile with
//  PKG_CXXFLAGS="-I/home/nasadat/R/x86_64-pc-linux-gnu-library/3.6/RcppArmadillo/include -I/home/nasadat/R/x86_64-pc-linux-gnu-library/3.6/Rcpp/include" PKG_LIBS=`Rscript -e 'Rcpp:::LdFlags()'` R CMD SHLIB ~/Documents/repos/dgdecomp_repo/dgdecomp/src/armacross.cpp



////////// We can write a C++ function as is
////////// But then we must pass it through SEXP
////////// if we are to call it from R

///// Here, I wrote a mostly-pure C++ function
///// and then an export function with SEXP inputs and outputs



// Create n-choose-k value
int nchook(const int N, const int K) {
  int nCk = 0;
  // if (K * 2 > N) K = N - K;
  if (K == 0) {
    nCk = 1;
  } else {
    nCk = N;
    for ( int i = 2; i <= K; ++i ) {
      nCk *= (N - i + 1);
      nCk /= i;
    }
  }

  return nCk;
}

int fact(int n);
double nCr(const int n, const int r)
{
  return fact(n) / (fact(r) * fact(n - r));
}

// Returns factorial of n
int fact(const int n)
{
  int res = 1;
  for (int i = 2; i <= n; i++)
    res = res * i;
  return res;
}


// Create cominations
arma::mat ArmaCombn_(int N, int K) {

  std::string bitmask(K, 1); // K leading 1s
  bitmask.resize(N, 0); // N-K trailing 0s


  // Create n-choose-k value
  int nCk = nCr(N, K);

  // Create output matrix
  int matrow = nCk;
  int matcol = K;

  // print integers and permute bitmask
  // and store in a VECTOR first
  arma::vec ech(matrow * matcol);
  int echh = 0;
  do {
    for (int i = 0; i < N; ++i) // [0..N-1] integers
    {
      if (bitmask[i]) {
        ech[echh] = i;
        echh++;
        // std::cout << " " << i;
      }
    }
    // std::cout << std::endl;
  } while (std::prev_permutation(bitmask.begin(), bitmask.end()));


  // Output matrix (from vector)
  arma::mat outmat(matrow, matcol);

  int counter = 0;
  for (int i = 0; i < matrow; i++) {
    for (int j = 0; j < matcol; j++) {
      outmat(i, j) = ech[counter] + 1; // +1 for usage in R space
      counter++;
    }
  }


  // Finally, transpose this matrix and return!
  return arma::trans(outmat);

}


// [[Rcpp::Export]]
RcppExport SEXP ArmaCombn(
  SEXP N_, SEXP K_) {

  return Rcpp::wrap(ArmaCombn_(
                      as<int>(N_), as<int>(K_)
                    ));

}

// [[Rcpp::Export]]
RcppExport SEXP ArmaCombnTEST(
  SEXP N_, SEXP K_) {

  arma::mat mat1 = ArmaCombn_(
                     as<int>(N_), as<int>(K_)
                   );
  arma::mat mat2 = arma::reverse(ArmaCombn_(
                                   as<int>(N_), as<int>(N_) - as<int>(K_)
                                 ), 1);

  return Rcpp::List::create(mat1, mat2);

}


//// Output the positioned matrix (Func_Cross)
arma::mat ArmaCross_(
  arma::mat vec_x, arma::mat vec_y,
  arma::mat vec_x_pos, arma::mat vec_y_pos) {


  // mkl_set_num_threads(10);
  // omp_set_num_threads(10);

  // Create output matrix
  arma::mat prod_inner(vec_x.n_rows, vec_x_pos.n_cols);

  // Temp store
  double tmpprod = 1.;

  for ( int row1 = 0; row1 < prod_inner.n_rows; row1++) {

    for ( int col1 = 0; col1 < vec_x_pos.n_cols; col1++) { // invariant between vec_x_pos and vec_y_pos

      // Go over each row of vec_x_pos
      for ( int row2 = 0; row2 < vec_x_pos.n_rows; row2++) { // vec_x_pos specific
        tmpprod = tmpprod * vec_x(row1, vec_x_pos(row2, col1) - 1); // -1 correction for zero index
      }

      // Go over each row of vec_y_pos
      for ( int row3 = 0; row3 < vec_y_pos.n_rows; row3++) { // vec_x_pos specific
        tmpprod = tmpprod * vec_y(row1, vec_y_pos(row3, col1) - 1); // -1 correction for zero index
      }

      // Record in prod_inner
      prod_inner(row1, col1) = tmpprod;
      tmpprod = 1.;

    }

  }

  return prod_inner;

}


//// The returner
// [[Rcpp::Export]]
RcppExport SEXP ArmaCross(
  SEXP vec_x_, SEXP vec_y_,
  SEXP vec_x_pos_, SEXP vec_y_pos_) {

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



// // [[Rcpp::export]]
// NumericMatrix callFunction(int x, int size1, int size2, Function f) {
//     NumericMatrix res = (x);
//     return res;
// }

// Create numerator (Func_Num)
// which will be a vector
arma::vec ArmaNum_(
  const int P, const int r,
  arma::mat vec_x, arma::mat vec_y) {

  // Create output placeholder
  arma::vec this_count(vec_x.n_rows);

  if (r == 1) {
    this_count = arma::prod(vec_x, 1) + arma::prod(vec_y, 1);
    return this_count;
  } else {

    // Pick : P-r from SMALL and mult with r-1 from CAPS
    // with : P-r from CAPS and mult with r-1 from SMALL

    int size1 = P - r;  // r < P and so will never be >= P-1
    int size2 = r - 1; //  Will never be <= 1 because r > 1

    // Create combinations of size1 and size2,
    // where size2 will be a reversal
    arma::mat vec_x_pos = ArmaCombn_(vec_x.n_cols, size1);
    arma::mat vec_y_pos = arma::reverse(ArmaCombn_(vec_x.n_cols, size2), 1);

    // printf("%4.3f", vec_x_pos.size());

    // First, P-r small and r-1 caps:
    this_count = arma::sum(ArmaCross_(vec_x, vec_y, vec_x_pos, vec_y_pos), 1);

    // Next, P-r caps and r-1 small:
    // ONLY applicable if we are not comparing identical sizes
    if (size1 != size2) {
      this_count = this_count + arma::sum(ArmaCross_(vec_x, vec_y, vec_y_pos, vec_x_pos), 1);
    }

    return this_count;

  }

}

// [[Rcpp::Export]]
RcppExport SEXP ArmaNum(
  const SEXP P_, const SEXP r_,
  SEXP vec_x_, SEXP vec_y_) {


  return Rcpp::wrap(ArmaNum_(
                      as<int>(P_), as<int>(r_),
                      as<arma::mat>(vec_x_),
                      as<arma::mat>(vec_y_)
                    ));

}


// Inner frac (Func_Inner)
arma::vec ArmaInner_(
  const int P, const int r,
  arma::mat vec_x, arma::mat vec_y) {

  if (r == 1) {
    return (arma::prod(vec_x, 1) + arma::prod(vec_y, 1)) / P;
  } else {
    return ArmaNum_(P, r, vec_x, vec_y) / (P * nchook(P - 1, r - 1));
  }

}


// [[Rcpp::Export]]
RcppExport SEXP ArmaInner(
  const SEXP P_, const SEXP r_,
  SEXP vec_x_, SEXP vec_y_) {

  return Rcpp::wrap(ArmaInner_(
                      as<int>(P_), as<int>(r_),
                      as<arma::mat>(vec_x_),
                      as<arma::mat>(vec_y_)
                    ));

}


// Inner sum (Func_Inner_Sum) [returns vector of size vec_x_.n_rows]
// this will loop over 1:P_upper and call ArmaInner, and add the values
arma::colvec ArmaInnerSum_ (
  const int P,
  const arma::mat vec_x, const arma::mat vec_y) {

  // Define the upper bound of the loop (P_upper)
  int P_upper = ( P % 2 == 0 ? (P / 2) : ((P + 1) * 0.5) );

  // Allocate output vector
  arma::vec sum_count = zeros<vec>(vec_x.n_rows);

  // Loop
  for (int Rx = 1; Rx <= P_upper; Rx++) {
    sum_count = sum_count + ArmaInner_(
                  P, Rx, vec_x, vec_y
                );
  }

  return sum_count;

}


// [[Rcpp::Export]]
RcppExport SEXP ArmaInnerSum(
  const SEXP P_,
  const SEXP vec_x_, const SEXP vec_y_) {

  return Rcpp::wrap(ArmaInnerSum_(
                      as<int>(P_),
                      as<arma::mat>(vec_x_),
                      as<arma::mat>(vec_y_)
                    ));

}


// The marginal effect loop inside Decomp_Factors_Matx
arma::mat ArmaDFInnerLoop_(
  const int num_facts,
  const arma::mat mat_x, const arma::mat mat_y) {

  // Allocate output matrix
  arma::mat outmat(mat_x.n_rows, num_facts);

  // We'll need to make copies of the inputs
  arma::mat mat_x_copy = mat_x;
  arma::mat mat_y_copy = mat_y;


  // We will loop over the factors 1 through num_facts
  for (int x = 1; x < num_facts + 1; x++) {

    // Remove the x'th factor from input matrix and
    // use over ArmaInnerSum_()
    arma::mat xshed = mat_x;
    arma::mat yshed = mat_y;

    // We need to refer from the original input copies everytime
    // because shed_col() is a self-inflicted void funk
    xshed.shed_col(x - 1);
    yshed.shed_col(x - 1);

    outmat.col(x - 1) = ArmaInnerSum_(num_facts, xshed, yshed);

    // Multiply with the sliced matrix
    arma::colvec tmpcol = mat_y.col(x - 1) - mat_x.col(x - 1);

    outmat.col(x - 1) = outmat.col(x - 1) % tmpcol;

  }

  return outmat;

}

// [[Rcpp::Export]]
RcppExport SEXP ArmaDFInnerLoop(
  const SEXP num_facts,
  SEXP mat_x_, SEXP mat_y_) {

  return Rcpp::wrap(ArmaDFInnerLoop_(
                      as<int>(num_facts),
                      as<arma::mat>(mat_x_),
                      as<arma::mat>(mat_y_)
                    ));

}
