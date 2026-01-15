// [[Rcpp::depends(RcppArmadillo)]]
// #include <RcppArmadillo.h>
#include <Rcpp.h>

#include <vector>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <limits>


using namespace Rcpp;

// [[Rcpp::interfaces(r, cpp)]]


const double pi = 3.1415926535897;

///********************************************************************
///** pbv_rcpp_pnorm0
double pbv_rcpp_pnorm0( double z)
{
    double y = ::Rf_pnorm5(z, 0.0, 1.0, 1, 0);
    //--- OUTPUT
    return y;
}
///********************************************************************

///********************************************************************
///** pbv_rcpp_pnorm
Rcpp::NumericVector pbv_rcpp_pnorm( Rcpp::NumericVector x)
{
    int N = x.size();
    Rcpp::NumericVector y(N);
    for (int nn=0; nn<N; nn++){
        y[nn] = pbv_rcpp_pnorm0( x[nn] );
    }
    //--- OUTPUT
    return y;
}
///********************************************************************



///********************************************************************
///** Drezner & Wesolowksy, 1990, JCSC
///** pbv_rcpp_pbvnorm0
double pbv_rcpp_pbvnorm0( double h1, double hk, double r)
{
    int NX=5;
    Rcpp::NumericVector X(NX);
    Rcpp::NumericVector W(NX);
    // data
    X[0]=.04691008;
    X[1]=.23076534;
    X[2]=.5;
    X[3]=.76923466;
    X[4]=.95308992;
    W[0]=.018854042;
    W[1]=.038088059;
    W[2]=.0452707394;
    W[3]=.038088059;
    W[4]=.018854042;
    // declarations
    double bv = 0;
    double r1, r2, rr, rr2, r3, h3, h5, h6, h7, aa, ab, h11;
    double cor_max = 0.7;
    double bv_fac1 = 0.13298076;
    double bv_fac2 = 0.053051647;

    // computation
    double h2 = hk;
    double h12 = (h1*h1+h2*h2)/2;
    double r_abs = std::abs(r);
    if (r_abs > cor_max){
        r2 = 1.0 - r*r;
        r3 = std::sqrt(r2);
        if (r<0){
            h2 = -h2;
        }
        h3 = h1*h2;
        h7 = std::exp( -h3 / 2.0);
        if ( r_abs < 1){
            h6 = std::abs(h1-h2);
            h5 = h6*h6 / 2.0;
            h6 = h6 / r3;
            aa = 0.5 - h3 / 8.0;
            ab = 3.0 - 2.0 * aa * h5;
            bv = bv_fac1*h6*ab*(1-pbv_rcpp_pnorm0(h6))-std::exp(-h5/r2)*(ab + aa*r2)*bv_fac2;
            for (int ii=0; ii<NX; ii++){
                r1 = r3*X[ii];
                rr = r1*r1;
                r2 = std::sqrt( 1.0 - rr);
                bv += - W[ii]*std::exp(- h5/rr)*(std::exp(-h3/(1.0+r2))/r2/h7 - 1.0 - aa*rr);
            }
        }
        h11 = std::min(h1,h2);
        bv = bv*r3*h7 + pbv_rcpp_pnorm0(h11);
        if (r < 0){
            bv = pbv_rcpp_pnorm0(h1) - bv;
        }

    } else {
        h3=h1*h2;
        for (int ii=0; ii<NX; ii++){
            r1 = r*X[ii];
            rr2 = 1.0 - r1*r1;
            bv += W[ii] * std::exp(( r1*h3 - h12)/rr2)/ std::sqrt(rr2);
        }
        bv = pbv_rcpp_pnorm0(h1)*pbv_rcpp_pnorm0(h2) + r*bv;
    }
    //--- OUTPUT
    return bv;
}
///********************************************************************




///********************************************************************
///** pbv_rcpp_pbvnorm
Rcpp::NumericVector pbv_rcpp_pbvnorm( Rcpp::NumericVector x, Rcpp::NumericVector y,
        Rcpp::NumericVector rho)
{
    int N = x.size();
    Rcpp::NumericVector res(N);
    for (int ii=0; ii<N; ii++){
        res[ii] = pbv_rcpp_pbvnorm0(x[ii], y[ii], rho[ii]);
    }
    //--- OUTPUT
    return res;
}
///********************************************************************



////** own function pvivnorm_wls
double pbivnorm__rcpp_wls0(double x, double y, double rho) {
  // Define infinity constants
  const double pos_inf = std::numeric_limits<double>::infinity();
  const double neg_inf = -std::numeric_limits<double>::infinity();
  
  // Check for infinite cases
  if (x == pos_inf && y == pos_inf) {
    return 1.0;
  } else if (x == neg_inf || y == neg_inf) {
    return 0.0;
  } else if (x == pos_inf) {
    return pbv_rcpp_pnorm0(y);
  } else if (y == pos_inf) {
    return pbv_rcpp_pnorm0(x);
  } else {
    return pbv_rcpp_pbvnorm0(x, y, rho);
  }
}

Rcpp::NumericVector pbivnorm__rcpp_wls( Rcpp::NumericVector x, Rcpp::NumericVector y,
                                      Rcpp::NumericVector rho)
{
  int N = x.size();
  Rcpp::NumericVector res(N);
  for (int ii=0; ii<N; ii++){
    res[ii] = pbivnorm__rcpp_wls0(x[ii], y[ii], rho[ii]);
  }
  //--- OUTPUT
  return res;
}

///********************************************************************



///adjust so that it begins to count with with 0!
std::vector<std::vector<int>> getCols(const std::vector<int>& lv, int nvar) {
  std::vector<std::vector<int>> result(nvar, std::vector<int>(2));
  int maxcol = -1;
  
  for (int i = 0; i < nvar; i++) {
    int mincol = maxcol + 1;
    maxcol = maxcol + lv[i] - 1;
    result[i][0] = mincol;
    result[i][1] = maxcol;
  }
  
  return result;
}




// [[Rcpp::export]]
List get_unique_values_per_column(const NumericMatrix& X) {
  int ncol = X.ncol();
  List result(ncol);
  
  for (int col = 0; col < ncol; col++) {
    // Extract column
    NumericVector column = X(_, col);
    int nrow = column.size();
    
    // Use unordered_set for O(1) insertion
    std::unordered_set<double> unique_set;
    unique_set.reserve(nrow);
    
    for (int i = 0; i < nrow; i++) {
      unique_set.insert(column[i]);
    }
    
    // Convert to vector and sort
    std::vector<double> unique_vec(unique_set.begin(), unique_set.end());
    std::sort(unique_vec.begin(), unique_vec.end());
    
    // Convert to NumericVector
    result[col] = NumericVector(unique_vec.begin(), unique_vec.end());
  }
  
  return result;
}


// [[Rcpp::export]]
List create_combs(int nvar, const NumericMatrix& sigma_mat) {
  if (sigma_mat.nrow() != nvar || sigma_mat.ncol() != nvar) {
    stop("sigma_mat must be a square matrix of size nvar x nvar");
  }
  
  int n_combs = nvar * (nvar - 1) / 2;
  
  // Create vectors for each row
  IntegerVector var1(n_combs);
  IntegerVector var2(n_combs);
  NumericVector corr(n_combs);
  
  int idx = 0;
  for (int i = 0; i < nvar; i++) {
    for (int j = i + 1; j < nvar; j++) {
      var1[idx] = i;
      var2[idx] = j;
      corr[idx] = sigma_mat(i, j);
      idx++;
    }
  }
  
  // Return as a list (like a data frame)
  return List::create(
    Named("var1") = var1,
    Named("var2") = var2,
    Named("correlation") = corr
  );
}





// [[Rcpp::export]]
NumericVector apply_get_joint_exp(List combs,
                                  NumericVector th,
                                  IntegerVector lv,
                                  int nvar,
                                  List catvals) {
  
  IntegerVector var1_vec = combs["var1"];
  IntegerVector var2_vec = combs["var2"];
  NumericVector corr_vec = combs["correlation"];
  
  int ncols = var1_vec.size();
  
  if (var2_vec.size() != ncols || corr_vec.size() != ncols) {
    stop("All elements in combs list must have the same length");
  }
  
  // Precompute threshold vectors for each variable
  std::vector<std::vector<double>> th_vectors(nvar);
  std::vector<std::vector<int>> catvals_vec(nvar);
  
  // Precompute column ranges
  std::vector<std::pair<int, int>> col_ranges(nvar);
  int maxcol = -1;
  
  for (int var = 0; var < nvar; var++) {
    int mincol = maxcol + 1;
    maxcol = maxcol + lv[var] - 1;
    col_ranges[var] = std::make_pair(mincol, maxcol);
    
    // Convert catvals for this variable
    IntegerVector temp_catvals = catvals[var];
    catvals_vec[var] = as<std::vector<int>>(temp_catvals);
    
    // Create threshold vector for this variable with -Inf, actual thresholds, Inf
    std::vector<double> th_var;
    th_var.push_back(-std::numeric_limits<double>::infinity());
    
    int start = col_ranges[var].first;
    int end = col_ranges[var].second;
    for (int i = start; i <= end; i++) {
      th_var.push_back(th[i]);
    }
    
    th_var.push_back(std::numeric_limits<double>::infinity());
    th_vectors[var] = th_var;
  }
  
  // Convert threshold vector once
  std::vector<double> th_vec = as<std::vector<double>>(th);
  std::vector<int> lv_vec = as<std::vector<int>>(lv);
  
  NumericVector results(ncols);
  
  // Process each pair
  for (int col = 0; col < ncols; col++) {
    int var1_idx = var1_vec[col];
    int var2_idx = var2_vec[col];
    double correlation = corr_vec[col];
    
    // Get references to avoid copies
    const std::vector<double>& th_var1 = th_vectors[var1_idx];
    const std::vector<double>& th_var2 = th_vectors[var2_idx];
    const std::vector<int>& vals_var1 = catvals_vec[var1_idx];
    const std::vector<int>& vals_var2 = catvals_vec[var2_idx];
    
    int n_cats1 = lv_vec[var1_idx];
    int n_cats2 = lv_vec[var2_idx];
    
    double mu_joint = 0.0;
    
    // Process all category combinations
    for (int i = 1; i <= n_cats1; i++) {
      double th1_i = th_var1[i];
      double th1_im1 = th_var1[i-1];
      int val1 = vals_var1[i-1];
      
      for (int j = 1; j <= n_cats2; j++) {
        double th2_j = th_var2[j];
        double th2_jm1 = th_var2[j-1];
        int val2 = vals_var2[j-1];
        
        // Calculate probability for this category combination
        double term1 = pbivnorm__rcpp_wls0(th1_i, th2_j, correlation);
        double term2 = pbivnorm__rcpp_wls0(th1_im1, th2_j, correlation) * -1.0;
        double term3 = pbivnorm__rcpp_wls0(th1_i, th2_jm1, correlation) * -1.0;
        double term4 = pbivnorm__rcpp_wls0(th1_im1, th2_jm1, correlation);
        
        double p_katkat = term1 + term2 + term3 + term4;
        
        // Add weighted contribution
        mu_joint += val1 * val2 * p_katkat;
      }
    }
    
    results[col] = mu_joint;
  }
  
  return results;
}



// [[Rcpp::export]]
Rcpp::NumericVector apply_get_mus(const std::vector<double>& th, 
                                  const std::vector<int>& lv, 
                                  int nvar, 
                                  const std::vector<std::vector<int>>& catvals){
  
  
  Rcpp::NumericVector result;
  
  //get cols
  std::vector<std::vector<int>> selcols = getCols(lv, nvar);
  
  for(int var = 0; var < nvar; var++){
    
    // Get threshold indices
    int wth_start = selcols[var][0];
    int wth_end = selcols[var][1];
    
    //Get thresholds
    std::vector<double> p_item;
    for (int i = wth_start; i <= wth_end; i++){
      p_item.push_back(th[i]);
    }
    
    //get probas over categories
    // Get max category index
    int max_lv_idx = lv[var] - 1; 
    
    double mu = 0.0;
    
    for (int i = 0; i<= max_lv_idx; i++){
      double prob_cat;
      double cat_val = catvals[var][i];
      
      if (i==0){
        prob_cat = pbv_rcpp_pnorm0(p_item[0]);
      } else if (i == max_lv_idx){
        prob_cat = pbv_rcpp_pnorm0(p_item.back() * -1);
      } else {
        prob_cat = pbv_rcpp_pnorm0(p_item[i]) - pbv_rcpp_pnorm0(p_item[i-1]);
      }
      
      mu += cat_val * prob_cat; 
    }
    
    result.push_back(mu);
  }
  
  return result;
}





// [[Rcpp::export]]
NumericMatrix create_weight_matrix(const NumericVector& th_r,
                                   const IntegerVector& lv_r,
                                   int nvar,
                                   const NumericMatrix& polychors) {
  
  // Convert inputs
  std::vector<double> th = as<std::vector<double>>(th_r);
  std::vector<int> lv = as<std::vector<int>>(lv_r);
  int n_th = th.size();
  
  // Initialize matrix
  NumericMatrix weight_mat(n_th, n_th);
  
  // Pre-compute selcols
  std::vector<std::vector<int>> selcols = getCols(lv, nvar);
  
  // Create seqnc vector - NOTE: var is now 0-based
  std::vector<std::pair<int, int>> seqnc_pairs;  // Store (variable, category) pairs
  for (int var = 0; var < nvar; var++) {
    for (int cat = 1; cat <= lv[var] - 1; cat++) {
      seqnc_pairs.push_back(std::make_pair(var, cat));
    }
  }
  
  int n_pairs = seqnc_pairs.size();
  
  // Pre-compute threshold indices for each seqnc element
  std::vector<int> th_indices(n_pairs);
  for (int i = 0; i < n_pairs; i++) {
    int var = seqnc_pairs[i].first;
    int cat = seqnc_pairs[i].second;
    th_indices[i] = selcols[var][0] + (cat - 1);  // selcols returns 0-based indices
  }
  
  // Compute all pairwise sigmas and fill matrix
  for (int i = 0; i < n_pairs; i++) {
    int var_i = seqnc_pairs[i].first;
    int cat_i = seqnc_pairs[i].second;
    int th_idx_i = th_indices[i];
    
    // Get threshold for this category
    double th1_val;
    if (cat_i <= lv[var_i] - 1) {
      int th_pos_i = selcols[var_i][0] + (cat_i - 1);  // selcols returns 0-based indices
      th1_val = -th[th_pos_i];
    } else {
      th1_val = R_PosInf;
    }
    
    for (int j = i + 1; j < n_pairs; j++) {
      int var_j = seqnc_pairs[j].first;
      int cat_j = seqnc_pairs[j].second;
      int th_idx_j = th_indices[j];
      
      // Get correlation
      double p;
      if (var_i == var_j) {
        p = 1.0;
      } else {
        p = polychors(var_i, var_j);  // 0-based indexing
      }
      
      // Get second threshold
      double th2_val;
      if (cat_j <= lv[var_j] - 1) {
        int th_pos_j = selcols[var_j][0] + (cat_j - 1);  // selcols returns 0-based indices
        th2_val = -th[th_pos_j];
      } else {
        th2_val = R_PosInf;
      }
      
      // Compute sigma
      double p_katkat = pbv_rcpp_pbvnorm0(th1_val, th2_val, p);
      double pnorm_th1 = pbv_rcpp_pnorm0(th1_val);
      double pnorm_th2 = pbv_rcpp_pnorm0(th2_val);
      double sigma = p_katkat - pnorm_th1 * pnorm_th2;
      
      // Fill matrix (lower triangle)
      if (th_idx_i > th_idx_j) {
        weight_mat(th_idx_i, th_idx_j) = sigma;
      } else {
        weight_mat(th_idx_j, th_idx_i) = sigma;
      }
    }
  }
  
  // Make symmetric
  for (int i = 0; i < n_th; i++) {
    for (int j = i + 1; j < n_th; j++) {
      weight_mat(i, j) = weight_mat(j, i);
    }
  }
  
  // Fill diagonal
  for (int i = 0; i < n_th; i++) {
    double th_pr = pbv_rcpp_pnorm0(-th[i]);
    weight_mat(i, i) = th_pr * (1.0 - th_pr);
  }
  
  return weight_mat;
}


