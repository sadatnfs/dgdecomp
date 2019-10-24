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
int nchook(int N, int K) {
  int nCk = 0;
  if (K * 2 > N) K = N-K;
  if (K == 0) {
    nCk = 1;
  } else {
    nCk = N;
    for( int i = 2; i <= K; ++i ) {
      nCk *= (N-i+1);
      nCk /= i;
    }
  }
  
  return nCk;
}

// Create cominations
arma::mat ArmaCombn_(int N, int K) {
  
  std::string bitmask(K, 1); // K leading 1s
  bitmask.resize(N, 0); // N-K trailing 0s
  
  
  // Create n-choose-k value
  int nCk = nchook(N, K);
  
  // Create output matrix
  int matrow = nCk;
  int matcol = K;
  
  // print integers and permute bitmask
  // and store in a VECTOR first
  arma::vec ech(matrow*matcol);
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
  for(int i = 0; i < matrow; i++) {
    for(int j = 0; j < matcol; j++) {
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
    SEXP P_, SEXP r_,
    SEXP vec_x_, SEXP vec_y_) {
  
  
  return Rcpp::wrap(ArmaNum_(
      as<int>(P_), as<int>(r_),
      as<arma::mat>(vec_x_),
      as<arma::mat>(vec_y_)
  ));
  
}



// Inner frac (Func_Inner)
arma::vec ArmaInner_(
    int P, int r,
    arma::mat vec_x, arma::mat vec_y) {
  
  if (r == 1) {
    return (arma::prod(vec_x, 1) + arma::prod(vec_y, 1)) / P;
  } else {
    return ArmaNum_(P, r, vec_x, vec_y) / (P * nchook(P - 1, r - 1));	
  }
  
}


// [[Rcpp::Export]]
RcppExport SEXP ArmaInner(
    SEXP P_, SEXP r_,
    SEXP vec_x_, SEXP vec_y_) {
  
  return Rcpp::wrap(ArmaInner_(
      as<int>(P_), as<int>(r_),
      as<arma::mat>(vec_x_),
      as<arma::mat>(vec_y_)
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
// 
// 	List outthis = List::create(discreture::combinations(as<int>(Pfac), as<int>(size1)));
// 	
// 	NumericMatrix outz(as<int>(Pfac), as<int>(size1));
// 	for(int row = 0; row < outz.rows(); row++) {
// 		for(int col = 0; col < outz.cols(); col++) {
// 			outz(row,col) = outthis[0](row)[col];
// 		}
// 	}
// 
// 
// 	return Rcpp::wrap(outz);

// }
