#include <iostream>
#include "mkl.h"
#include <omp.h>

// #include <iterator>
#include <vector>
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



//// Output the positioned matrix (Func_Cross)
arma::mat ArmaCross_(
    arma::mat vec_x, arma::mat vec_y,
    arma::mat vec_x_pos, arma::mat vec_y_pos) {

  
  mkl_set_num_threads(10);
  omp_set_num_threads(10);
  
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
    int P, int r,
    arma::mat vec_x, arma::mat vec_y,
    arma::mat vec_x_pos, arma::mat vec_y_pos) {

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
    SEXP P_, SEXP r_,
    SEXP vec_x_, SEXP vec_y_,
    SEXP vec_x_pos_, SEXP vec_y_pos_) {

	return Rcpp::wrap(ArmaNum_(
	                      as<int>(P_), as<int>(r_),
	                      as<arma::mat>(vec_x_),
	                      as<arma::mat>(vec_y_),
	                      as<arma::mat>(vec_x_pos_),
	                      as<arma::mat>(vec_y_pos_)
	                  ));

}


int fact(int n);   
double nCr(int n, int r) 
{ 
    return fact(n) / (fact(r) * fact(n - r)); 
} 
  
// Returns factorial of n 
int fact(int n) 
{ 
    int res = 1; 
    for (int i = 2; i <= n; i++) 
        res = res * i; 
    return res; 
} 


// Inner frac (Func_Inner)
arma::vec ArmaInner_(
    int P, int r,
    arma::mat vec_x, arma::mat vec_y,
    arma::mat vec_x_pos, arma::mat vec_y_pos) {

	if (r == 1) {
		return (arma::prod(vec_x, 1) + arma::prod(vec_y, 1)) / P;
	} else {
		return ArmaNum_(P, r, vec_x, vec_y, vec_x_pos, vec_y_pos) / (P * nCr(P - 1, r - 1));	
	}
	
}


// [[Rcpp::Export]]
RcppExport SEXP ArmaInner(
    SEXP P_, SEXP r_,
    SEXP vec_x_, SEXP vec_y_,
    SEXP vec_x_pos_, SEXP vec_y_pos_) {

	return Rcpp::wrap(ArmaInner_(
	                      as<int>(P_), as<int>(r_),
	                      as<arma::mat>(vec_x_),
	                      as<arma::mat>(vec_y_),
	                      as<arma::mat>(vec_x_pos_),
	                      as<arma::mat>(vec_y_pos_)
	                  ));

}





// arma::mat CustomCombn_(int x_, int m_) {

// 	// x <- seq_len(x)
// 	vec x = regspace(0, x_);
// 	int n = x.size();

// 	int e = 0;
// 	int h = m_;
// 	vec a = regspace(0, m_);

// 	int len_r = a.size();
// 	int count = NCR(n, m_);

// 	vec r(a.size());
// 	for (int rr = 0; rr < r.size(); r++) {
// 		r[rr] = x[a[rr]];
// 	}
// 	arma::mat outmat(len_r, count);

// 	for (int rowz = 0; rowz < outmat.n_rows; r++) {
// 		int j = 0;
// 		for (int colz = 0; colz < outmat.n_cols; r++) {
// 			outmat(rowz, colz) = r[j];
// 		}
// 	}

// 	int i = 2;
// 	// int j = 0;
// 	int nmmp1 = n - m_ + 1;
// 	while ( a[0] < nmmp1 ) {
// 		if (e < n - h) {
// 			h = 0;
// 			e = a[m_];
// 			int j = 0;
// 			a[m_ - h + j] = e + j;
// 		} else {
// 			e = a[m_ - h];
// 			h = h + 1;
// 			vec j = regspace(0, h);
// 			for (int jer = 0; jer < j.size(); jer++) {
// 				a[m_ - h + j[jer]] = e + j[jer];
// 			}
// 		}


// 		for (int rr = 0; rr < r.size(); r++) {
// 			outmat(rr, i) = x[a[rr]];
// 		}

// 		i++;

// 	}

// 	return outmat;

// }


// // [[Rcpp::Export]]
// RcppExport SEXP CustomCombn(
//     SEXP Pfac, SEXP size1, SEXP size2) {

// 	auto vec_x_sim = discreture::combinations(as<int>(Pfac), as<int>(size1));
// 	auto vec_y_sim = reversed(discreture::combinations(as<int>(Pfac), as<int>(size2)));

// 	// arma::mat vec_x_sim_out(as<int>(size1), NCR(as<int>(Pfac), 2));
// 	arma::mat vec_x_sim_out = as<arma::mat>(vec_x_sim);

// 	// // int i = 0;
// 	// for (int i = 0; i < vec_x_sim.size(); ++i) {
// 	//     vec_x_sim_out.col(i) = as<arma::colvec>(vec_x_sim[i]);
// 	//     // i++;
// 	// }

// 	return Rcpp::wrap(vec_x_sim_out);

// }
