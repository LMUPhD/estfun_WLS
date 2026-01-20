// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
// [[Rcpp::interfaces(r, cpp)]]

using namespace Rcpp;
using namespace arma;



// [[Rcpp::export]]
List x2GLIST_cpp(NumericVector x,
                 List m_free_idx,
                 List x_free_idx,
                 List GLIST_template,
                 LogicalVector isSymmetric = LogicalVector()) {
  
  // Create a deep copy of the template
  List GLIST = clone(GLIST_template);
  int n_matrices = GLIST.size();
  
  // Check if symmetric handling is needed
  bool do_symmetric = (isSymmetric.size() > 0);
  
  for (int mm = 0; mm < n_matrices; mm++) {
    // Get current matrix
    NumericMatrix mat = GLIST[mm];
    
    // Get index vectors for this matrix
    IntegerVector m_idx = m_free_idx[mm];
    IntegerVector x_idx = x_free_idx[mm];
    
    // Skip if no free parameters
    if (m_idx.size() == 0) continue;
    
    // Convert R indices (1-based) to C++ indices (0-based)
    IntegerVector m_idx_0 = m_idx - 1;
    IntegerVector x_idx_0 = x_idx - 1;
    
    // Map matrix positions
    int n_free = m_idx_0.size();
    for (int i = 0; i < n_free; i++) {
      int matrix_pos = m_idx_0[i];
      int param_pos = x_idx_0[i];
      
      // Convert flat index to row/col for matrix
      int nrow = mat.nrow();
      int row = matrix_pos % nrow;
      int col = matrix_pos / nrow;
      
      // Assign value
      mat(row, col) = x[param_pos];
    }
    
    // Make symmetric if needed
    if (do_symmetric && isSymmetric[mm]) {
      int n = mat.nrow();
      int m = mat.ncol();
      
      if (n == m) {  // Must be square to be symmetric
        for (int i = 0; i < n; i++) {
          for (int j = i + 1; j < n; j++) {
            mat(i, j) = mat(j, i);
          }
        }
      }
    }
  }
  
  return GLIST;
}


// [[Rcpp::export]]
List x2GLIST_withtheta(NumericVector x,
                        List m_free_idx,
                        List x_free_idx,
                        List GLIST_template,
                        LogicalVector isSymmetric) {
  
  List GLIST = clone(GLIST_template);
  int n_matrices = GLIST.size();
  bool do_symmetric = (isSymmetric.size() > 0);
  
  // Step 1: Assign free parameters
  for (int mm = 0; mm < n_matrices; mm++) {
    NumericMatrix mat = GLIST[mm];
    IntegerVector m_idx = m_free_idx[mm];
    IntegerVector x_idx = x_free_idx[mm];
    
    if (m_idx.size() == 0) continue;
    
    for (int i = 0; i < m_idx.size(); i++) {
      int pos = m_idx[i] - 1;
      int nrow = mat.nrow();
      int row = pos % nrow;
      int col = pos / nrow;
      mat(row, col) = x[x_idx[i] - 1];
    }
    
    if (do_symmetric && isSymmetric[mm]) {
      int n = mat.nrow();
      if (n == mat.ncol()) {
        for (int i = 0; i < n; i++) {
          for (int j = i + 1; j < n; j++) {
            mat(i, j) = mat(j, i);
          }
        }
      }
    }
  }
  
  // Step 2: Compute theta for delta parameterization (categorical)
  // Simplified: theta = 1 - diag(lambda %*% psi %*% t(lambda))
  // since delta = 1 for all variables
  
  NumericMatrix lambda = GLIST[0];
  NumericMatrix psi = GLIST[2];
  NumericMatrix theta = GLIST[1];
  
  int nvar = lambda.nrow();
  int nfac = lambda.ncol();
  
  // Compute lambda %*% psi %*% t(lambda)
  for (int i = 0; i < nvar; i++) {
    double diag_sum = 0.0;
    for (int j = 0; j < nfac; j++) {
      for (int k = 0; k < nfac; k++) {
        diag_sum += lambda(i, j) * psi(j, k) * lambda(i, k);
      }
    }
    // delta = 1, so 1/delta^2 = 1
    theta(i, i) = 1.0 - diag_sum;
  }
  
  GLIST[1] = theta;
  
  return GLIST;
}



// [[Rcpp::export]]
arma::mat compute_SigmaHat(const arma::mat& lambda,
                               const arma::mat& psi,
                               const arma::mat& theta) {
  return lambda * psi * lambda.t() + theta;  
}

// [[Rcpp::export]]
arma::mat compute_scores(const arma::mat& Delta, 
                             const arma::mat& W_inv, 
                             const arma::mat& e){

  arma::mat Scores = (Delta.t() * W_inv * e.t()).t();

  return(Scores);
  
}