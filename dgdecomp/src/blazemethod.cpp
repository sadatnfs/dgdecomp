#include <iostream>
#include <omp.h>

#include <vector>
#include <thread>
#include <algorithm>
#include <cstdlib>
#include <RcppEigen.h>
#include <RcppBlaze3.h>


// [[Rcpp::depends(RcppBlaze)]]

using namespace Rcpp;
using namespace std;
using namespace blaze;
using Rcpp::as;
using blaze::DynamicMatrix;
using blaze::DynamicVector;
using blaze::columnVector;
using blaze::rowwise;
using blaze::prod;


//// Compile code
// g++  -std=gnu++17 -I"/opt/R/lib/R/include" -DNDEBUG  -I'/data/opt/R/lib/R/library/Rcpp/include' -I'/data/opt/R/lib/R/library/RcppArmadillo/include' -I'/data/opt/R/lib/R/library/RcppEigen/include' -I'/data/opt/R/lib/R/library/RcppBlaze/include' -w -fopenmp -fno-gnu-unique -fno-optimize-sibling-calls -DMKL_ILP64 -m64 -I/opt/intel/compilers_and_libraries/linux/mkl/include -I/usr/include -I/opt/conda/include                             -msse2 -msse3 -msse4.1 -msse4.2 -mfma -mavx -mavx2 -ffp-contract=fast  -mfpmath=sse -g -O3   -fopenmp -lgomp  -fPIC   -w -fopenmp -fno-gnu-unique -fno-optimize-sibling-calls -DMKL_ILP64 -m64 -I/opt/intel/compilers_and_libraries/linux/mkl/include -I/usr/include -I/opt/conda/include                             -msse2 -msse3 -msse4.1 -msse4.2 -mfma -mavx -mavx2 -ffp-contract=fast  -mfpmath=sse -g -O3   -c /home/nasadat/blazetest.cpp -o /home/nasadat/blazetest.o; g++ -std=gnu++17 -shared -L/opt/R/lib/R/lib -L/usr/local/lib -liconv -L/opt/conda/lib -licui18n -licuuc -licudata -L/opt/hadoop/lib/native:/opt/intel/compilers_and_libraries/linux/mkl/lib/intel64_lin:/opt/intel/lib/intel64_lin:/usr/lib:/usr/local/lib:/opt/intel/compilers_and_libraries_2020.1.217/linux/tbb/lib/intel64_lin/gcc4.7:/opt/intel/compilers_and_libraries_2020.1.217/linux/compiler/lib/intel64_lin:/opt/intel/compilers_and_libraries_2020.1.217/linux/mkl/lib/intel64_lin:/opt/hadoop/lib/native:/opt/intel/compilers_and_libraries/linux/mkl/lib/intel64_lin:/opt/intel/lib/intel64_lin:/usr/lib:/usr/local/lib:/opt/intel/compilers_and_libraries_2020.1.217/linux/tbb/lib/intel64_lin/gcc4.7:/opt/intel/compilers_and_libraries_2020.1.217/linux/compiler/lib/intel64_lin:/opt/intel/compilers_and_libraries_2020.1.217/linux/mkl/lib/intel64_lin:/opt/conda/lib:/opt/conda/lib -L/usr/lib -o /home/nasadat/blazetest.so /home/nasadat/blazetest.o -fopenmp -lgomp -L/opt/R/lib/R/lib -lR

// Create n-choose-k value
long int nchook3(const int N, const int K) {
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

// [[Rcpp::Export]]
RcppExport SEXP Blaze__nchook3(
    SEXP N_, SEXP K_) {
  
  return Rcpp::wrap(nchook3(
      as<int>(N_), as<int>(K_)
  ));
  
}

// Create cominations
blaze::DynamicMatrix<double> BlazeCombn_(int N, int K) {
  
  std::string bitmask(K, 1); // K leading 1s
  bitmask.resize(N, 0); // N-K trailing 0s
  
  
  // Create n-choose-k value
  long int nCk = nchook3(N, K);
  
  // Create output matrix
  int matrow = nCk;
  int matcol = K;
  
  // print integers and permute bitmask
  // and store in a VECTOR first
  blaze::DynamicVector<double> ech(matrow * matcol);
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
  blaze::DynamicMatrix<double> outmat(matrow, matcol);
  
  int counter = 0;
  for (int i = 0; i < matrow; i++) {
    for (int j = 0; j < matcol; j++) {
      outmat(i, j) = ech[counter] + 1; // +1 for usage in R space
      counter++;
    }
  }
  
  
  // Finally, transpose this matrix and return!
  transpose(outmat);
  return outmat;
  
}


// [[Rcpp::Export]]
RcppExport SEXP Blaze__Combn(
    SEXP N_, SEXP K_) {
  
  return Rcpp::wrap(BlazeCombn_(
      as<int>(N_), as<int>(K_)
  ));
  
}


//// Output the positioned matrix (Func_Cross)
blaze::DynamicMatrix<double> BlazeCross_(
    blaze::DynamicMatrix<double> vec_x, blaze::DynamicMatrix<double> vec_y,
    blaze::DynamicMatrix<double> vec_x_pos, blaze::DynamicMatrix<double> vec_y_pos) {
  
  
  // Create output matrix
  blaze::DynamicMatrix<double> prod_inner(vec_x.rows(), vec_x_pos.columns());
  
  // Temp store
  double tmpprod = 1.;
  
  for ( int row1 = 0; row1 < prod_inner.rows(); row1++) {
    
    for ( int col1 = 0; col1 < vec_x_pos.columns(); col1++) { 
      // invariant between vec_x_pos and vec_y_pos
      
      // Go over each row of vec_x_pos
      for ( int row2 = 0; row2 < vec_x_pos.rows(); row2++) { 
        // vec_x_pos specific
        tmpprod = tmpprod * vec_x(row1, vec_x_pos(row2, col1) - 1); 
        // -1 correction for zero index
      }
      
      // Go over each row of vec_y_pos
      for ( int row3 = 0; row3 < vec_y_pos.rows(); row3++) { 
        // vec_x_pos specific
        tmpprod = tmpprod * vec_y(row1, vec_y_pos(row3, col1) - 1); 
        // -1 correction for zero index
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
RcppExport SEXP Blaze__Cross(
    SEXP vec_x_, SEXP vec_y_,
    SEXP vec_x_pos_, SEXP vec_y_pos_) {
  
  // Must coerce vector to blaze::DynamicVector<double> type (or whatever type to feed into funk)
  // so that we convert from input type from R (SEXP) to whatever our funk needs
  
  return Rcpp::wrap(
    BlazeCross_(
      as<blaze::DynamicMatrix<double>>(vec_x_),
      as<blaze::DynamicMatrix<double>>(vec_y_),
      as<blaze::DynamicMatrix<double>>(vec_x_pos_),
      as<blaze::DynamicMatrix<double>>(vec_y_pos_)
    ));
  
}




// Create numerator (Func_Num)
// which will be a vector
blaze::DynamicVector<double> BlazeNum_(
    const int P, const int r,
    blaze::DynamicMatrix<double> vec_x, blaze::DynamicMatrix<double> vec_y) {
  
  // Create output placeholder
  blaze::DynamicVector<double> this_count(vec_x.rows());
  
  if (r == 1) {
    this_count = blaze::prod<rowwise>(vec_x) + blaze::prod<rowwise>(vec_y);
    return this_count;
  } else {
    
    // Pick : P-r from SMALL and mult with r-1 from CAPS
    // with : P-r from CAPS and mult with r-1 from SMALL
    
    int size1 = P - r;  // r < P and so will never be >= P-1
    int size2 = r - 1; //  Will never be <= 1 because r > 1
    
    // Create combinations of size1 and size2,
    // where size2 will be a reversal
    blaze::DynamicMatrix<double> vec_x_pos = BlazeCombn_(vec_x.columns(), size1);
    blaze::DynamicMatrix<double> vec_y_pos = blaze::reverse<rowwise>(BlazeCombn_(vec_x.columns(), size2));
    
    // printf("%4.3f", vec_x_pos.size());
    
    // First, P-r small and r-1 caps:
    this_count = blaze::sum<rowwise>(BlazeCross_(vec_x, vec_y, vec_x_pos, vec_y_pos));
    
    // Next, P-r caps and r-1 small:
    // ONLY applicable if we are not comparing identical sizes
    if (size1 != size2) {
      this_count = this_count + blaze::sum<rowwise>(
        BlazeCross_(vec_x, vec_y, vec_y_pos, vec_x_pos));
    }
    
    return this_count;
    
  }
  
}

// [[Rcpp::Export]]
RcppExport SEXP Blaze__Num(
    const SEXP P_, const SEXP r_,
    SEXP vec_x_, SEXP vec_y_) {
  
  
  return Rcpp::wrap(BlazeNum_(
      as<int>(P_), as<int>(r_),
      as<blaze::DynamicMatrix<double>>(vec_x_),
      as<blaze::DynamicMatrix<double>>(vec_y_)
  ));
  
}



// Inner frac (Func_Inner)
blaze::DynamicVector<double> BlazeInner_(
    const int P, const int r,
    blaze::DynamicMatrix<double> vec_x, blaze::DynamicMatrix<double> vec_y) {
  
  if (r == 1) {
    return (blaze::prod<rowwise>(vec_x) + blaze::prod<rowwise>(vec_y)) / P;
  } else {
    return BlazeNum_(P, r, vec_x, vec_y) / (P * nchook3(P - 1, r - 1));
  }
  
}


// [[Rcpp::Export]]
RcppExport SEXP Blaze__Inner(
    const SEXP P_, const SEXP r_,
    SEXP vec_x_, SEXP vec_y_) {
  
  return Rcpp::wrap(BlazeInner_(
      as<int>(P_), as<int>(r_),
      as<blaze::DynamicMatrix<double>>(vec_x_),
      as<blaze::DynamicMatrix<double>>(vec_y_)
  ));
  
}



// Inner sum (Func_Inner_Sum) [returns vector of size vec_x_.rows()]
// this will loop over 1:P_upper and call ArmaInner, and add the values
blaze::DynamicVector<double> BlazeInnerSum_ (
    const int P,
    const blaze::DynamicMatrix<double> vec_x, const blaze::DynamicMatrix<double> vec_y) {
  
  // Define the upper bound of the loop (P_upper)
  int P_upper = ( P % 2 == 0 ? (P / 2) : ((P + 1) * 0.5) );
  
  // Allocate output vector
  blaze::DynamicVector<double> sum_count(vec_x.rows());
  
  // Loop
  for (int Rx = 1; Rx <= P_upper; Rx++) {
    sum_count = sum_count + BlazeInner_(
      P, Rx, vec_x, vec_y
    );
  }
  
  return sum_count;
  
}


// [[Rcpp::Export]]
RcppExport SEXP Blaze__InnerSum(
    const SEXP P_,
    const SEXP vec_x_, const SEXP vec_y_) {
  
  return Rcpp::wrap(BlazeInnerSum_(
      as<int>(P_),
      as<blaze::DynamicMatrix<double>>(vec_x_),
      as<blaze::DynamicMatrix<double>>(vec_y_)
  ));
  
}




blaze::CompressedVector<double> BlazeSlice_(
    blaze::DynamicMatrix<double> mat_x
) {
  
  // Allocate output matrix
  blaze::CompressedVector<double> outmat;
  
  for(int row = 0; row < mat_x.rows(); row++) {
    outmat.insert(row, mat_x(row, 1));
  }
  
  return outmat;
}

// [[Rcpp::Export]]
RcppExport SEXP Blaze__Slice(
    const SEXP mat_x) {
  
  return Rcpp::wrap(BlazeSlice_(
      as<blaze::DynamicMatrix<double>>(mat_x)
  ));
  
}




// The marginal effect loop inside Decomp_Factors_Matx
blaze::DynamicMatrix<double> BlazeDFInnerLoop_(
  const int num_facts,
  const blaze::DynamicMatrix<double> mat_x, 
  const blaze::DynamicMatrix<double> mat_y,
  const int threads) {

  // Allocate output matrix
  blaze::DynamicMatrix<double> outmat(mat_x.rows(), num_facts);

  // // We will loop over the factors 1 through num_facts
  // #pragma omp parallel num_threads(threads)
  // {
  //   #pragma omp for
    for (int x = 1; x < num_facts + 1; x++) {
  
      // Remove the x'th factor from input matrix and
      // use over BlazeInnerSum_()
      blaze::DynamicMatrix<double> xshed = mat_x;
      blaze::DynamicMatrix<double> yshed = mat_y;
  
      // We need to refer from the original input copies everytime
      // because shed_col() is a self-inflicted void funk
      xshed.shed_col(x - 1);
      yshed.shed_col(x - 1);
  
      outmat.col(x - 1) = BlazeInnerSum_(num_facts, xshed, yshed);
  
      // Multiply with the sliced matrix
      Blaze::DynamicVector<double> tmpcol = mat_y.col(x - 1) - mat_x.col(x - 1);
  
      outmat.col(x - 1) = outmat.col(x - 1) % tmpcol;
  
    }
  // }

  return outmat;

}

// [[Rcpp::Export]]
RcppExport SEXP BlazeDFInnerLoop(
  const SEXP num_facts,
  SEXP mat_x_, SEXP mat_y_,
  const SEXP threads_ ) {

  return Rcpp::wrap(BlazeDFInnerLoop_(
                      as<int>(num_facts),
                      as<blaze::DynamicMatrix<double>>(mat_x_),
                      as<blaze::DynamicMatrix<double>>(mat_y_),
                      as<int>(threads_)
                    ));

}






