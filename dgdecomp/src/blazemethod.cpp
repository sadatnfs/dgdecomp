#include <iostream>
#include <omp.h>

#include <vector>
#include <thread>
#include <algorithm>
#include <cstdlib>
#include <RcppEigen.h>
#include <RcppBlaze.h>


// [[Rcpp::depends(RcppBlaze)]]

using namespace Rcpp;
using namespace std;
using namespace Eigen;
using Eigen::MatrixXd;
using Eigen::VectorXd;

using namespace blaze;
using blaze::DynamicMatrix;
using blaze::DynamicVector;

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

// long int fact(int n);
// long int nCr(const int n, const int r)
// {
//   return fact(n) / (fact(r) * fact(n - r));
// }
// 
// // Returns factorial of n
// long int fact(const int n)
// {
//   int res = 1;
//   for (int i = 2; i <= n; i++)
//     res = res * i;
//   return res;
// }







// Create cominations
MatrixXd BlazeCombn_(int N, int K) {
  
  std::string bitmask(K, 1); // K leading 1s
  bitmask.resize(N, 0); // N-K trailing 0s
  
  
  // Create n-choose-k value
  long int nCk = nchook3(N, K);
  
  // Create output matrix
  int matrow = nCk;
  int matcol = K;
  
  // print integers and permute bitmask
  // and store in a VECTOR first
  VectorXd ech(matrow * matcol);
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
  MatrixXd outmat(matrow, matcol);
  
  int counter = 0;
  for (int i = 0; i < matrow; i++) {
    for (int j = 0; j < matcol; j++) {
      outmat(i, j) = ech[counter] + 1; // +1 for usage in R space
      counter++;
    }
  }
  
  
  // Finally, transpose this matrix and return!
  outmat.transposeInPlace();
  return outmat;
  
}


// [[Rcpp::Export]]
RcppExport SEXP BlazeCombn(
    SEXP N_, SEXP K_) {
  
  return Rcpp::wrap(BlazeCombn_(
      as<int>(N_), as<int>(K_)
  ));
  
}



//// Output the positioned matrix (Func_Cross)
MatrixXd BlazeCross_(
    MatrixXd vec_x, MatrixXd vec_y,
    MatrixXd vec_x_pos, MatrixXd vec_y_pos) {
  
  
  // Create output matrix
  MatrixXd prod_inner(vec_x.rows(), vec_x_pos.cols());
  
  // Temp store
  double tmpprod = 1.;
  
  for ( int row1 = 0; row1 < prod_inner.rows(); row1++) {
    
    for ( int col1 = 0; col1 < vec_x_pos.cols(); col1++) { 
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
RcppExport SEXP BlazeCross(
    SEXP vec_x_, SEXP vec_y_,
    SEXP vec_x_pos_, SEXP vec_y_pos_) {
  
  // Must coerce vector to arma::vec type (or whatever type to feed into funk)
  // so that we convert from input type from R (SEXP) to whatever our funk needs
  
  return Rcpp::wrap(
    BlazeCross_(
      as<MatrixXd>(vec_x_),
      as<MatrixXd>(vec_y_),
      as<MatrixXd>(vec_x_pos_),
      as<MatrixXd>(vec_y_pos_)
    ));
  
}


// Create numerator (Func_Num)
// which will be a vector
VectorXd BlazeNum_(
    const int P, const int r,
    MatrixXd mat_x, MatrixXd mat_y) {
  
  // Create output placeholder
  VectorXd this_count(mat_x.rows());
  
  if (r == 1) {
    this_count = mat_x.rowwise().prod() + mat_y.rowwise().prod();
    return this_count;
  } else {
    
    // Pick : P-r from SMALL and mult with r-1 from CAPS
    // with : P-r from CAPS and mult with r-1 from SMALL
    
    int size1 = P - r;  // r < P and so will never be >= P-1
    int size2 = r - 1; //  Will never be <= 1 because r > 1
    
    // Create combinations of size1 and size2,
    // where size2 will be a reversal
    MatrixXd mat_x_pos = BlazeCombn_(mat_x.cols(), size1);
    MatrixXd mat_y_pos = BlazeCombn_(mat_x.cols(), size2).rowwise().reverse();
    
    // printf("%4.3f", mat_x_pos.size());
    
    // First, P-r small and r-1 caps:
    this_count = BlazeCross_(mat_x, mat_y, mat_x_pos, mat_y_pos).rowwise().sum();
    
    // Next, P-r caps and r-1 small:
    // ONLY applicable if we are not comparing identical sizes
    if (size1 != size2) {
      this_count = this_count + BlazeCross_(
        mat_x, mat_y, mat_y_pos, mat_x_pos).rowwise().sum();
    }
    
    return this_count;
    
  }
  
}



// Inner frac (Func_Inner)
VectorXd BlazeInner_(
    const int P, const int r,
    MatrixXd vec_x, MatrixXd vec_y) {
  
  if (r == 1) {
    return (vec_x.rowwise().prod() + vec_y.rowwise().prod()) / P;
  } else {
    return BlazeNum_(P, r, vec_x, vec_y) / (P * nchook3(P - 1, r - 1));
  }
  
}



// Inner sum (Func_Inner_Sum) [returns vector of size vec_x_.rows()]
// this will loop over 1:P_upper and call ArmaInner, and add the values
VectorXd BlazeInnerSum_ (
    const int P,
    const MatrixXd vec_x, const MatrixXd vec_y) {
  
  // Define the upper bound of the loop (P_upper)
  int P_upper = ( P % 2 == 0 ? (P / 2) : ((P + 1) * 0.5) );
  
  // Allocate output vector with zeroes
  VectorXd sum_count = vec_x * 0;
  
  // Loop
  for (int Rx = 1; Rx <= P_upper; Rx++) {
    sum_count = sum_count + BlazeInner_(
      P, Rx, vec_x, vec_y
    );
  }
  
  return sum_count;
  
}

// [[Rcpp::Export]]
RcppExport SEXP BlazeInnerSum(
    const SEXP P_,
    const SEXP vec_x_, const SEXP vec_y_) {
  
  return Rcpp::wrap(BlazeInnerSum_(
      as<int>(P_),
      as<MatrixXd>(vec_x_),
      as<MatrixXd>(vec_y_)
  ));
  
}


// void removeRow2(Eigen::MatrixXd& matrix, unsigned int rowToRemove)
// {
//   unsigned int numRows = matrix.rows()-1;
//   unsigned int numCols = matrix.cols();
//   
//   if( rowToRemove < numRows )
//     matrix.block(rowToRemove,0,numRows-rowToRemove,numCols) = matrix.bottomRows(numRows-rowToRemove);
//   
//   matrix.conservativeResize(numRows,numCols);
// }

void removeColumn2(Eigen::MatrixXd& matrix, unsigned int colToRemove)
{
  unsigned int numRows = matrix.rows();
  unsigned int numCols = matrix.cols()-1;
  
  if( colToRemove < numCols )
    matrix.block(0,colToRemove,numRows,numCols-colToRemove) = matrix.rightCols(numCols-colToRemove);
  
  matrix.conservativeResize(numRows,numCols);
}


// The marginal effect loop inside Decomp_Factors_Matx
MatrixXd BlazeDFInnerLoop_(
    const int num_facts,
    const MatrixXd mat_x, const MatrixXd mat_y,
    const int threads) {
  
  // Allocate output matrix
  MatrixXd outmat(mat_x.rows(), num_facts);
  
  // We will loop over the factors 1 through num_facts
#pragma omp parallel num_threads(threads)
{
#pragma omp for
  for (int x = 1; x < num_facts + 1; x++) {
    
    // Remove the x'th factor from input matrix and
    // use over ArmaInnerSum_()
    MatrixXd xshed = mat_x;
    MatrixXd yshed = mat_y;
    
    // We need to refer from the original input copies everytime
    // because shed_col() is a self-inflicted void funk
    // xshed.shed_col(x - 1);
    // yshed.shed_col(x - 1);
    
    removeColumn2(xshed, x-1);
    removeColumn2(yshed, x-1);
    
    outmat.col(x - 1) = BlazeInnerSum_(num_facts, xshed, yshed);
    
    // Multiply with the sliced matrix
    VectorXd tmpcol = mat_y.col(x - 1) - mat_x.col(x - 1);
    
    outmat.col(x - 1) = outmat.col(x - 1) * tmpcol;
    
  }
}

return outmat;

}

// [[Rcpp::Export]]
RcppExport SEXP BlazeDFInnerLoop(
    const SEXP num_facts,
    SEXP mat_x_, SEXP mat_y_,
    const SEXP threads_ ) {
  
  return Rcpp::wrap(BlazeDFInnerLoop_(
      as<int>(num_facts),
      as<MatrixXd>(mat_x_),
      as<MatrixXd>(mat_y_),
      as<int>(threads_)
  ));
  
}
