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
// [[Rcpp::export]]
double pbv_rcpp_pnorm0( double z)
{
    double y = ::Rf_pnorm5(z, 0.0, 1.0, 1, 0);
    //--- OUTPUT
    return y;
}
///********************************************************************

///********************************************************************
///** pbv_rcpp_pnorm
// [[Rcpp::export]]
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
// [[Rcpp::export]]
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
// [[Rcpp::export]]
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
// [[Rcpp::export]]
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

// [[Rcpp::export]]
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

// [[Rcpp::export]]
std::vector<std::vector<int>> getCols(const std::vector<int>& lv, int nvar) {
  std::vector<std::vector<int>> result(nvar, std::vector<int>(2));
  int maxcol = 0;
  
  for (int i = 0; i < nvar; i++) {
    int mincol = maxcol + 1;
    maxcol = maxcol + lv[i] - 1;
    result[i][0] = mincol;
    result[i][1] = maxcol;
  }
  
  return result;
}


















// Alternative: Use a struct for clarity
struct VariablePair {
  int var1;
  int var2;
  double correlation;
};

double get_joint_exp(const VariablePair& pair,
                            const std::vector<double>& th, 
                            const std::vector<int>& lv, 
                            int nvar, 
                            const std::vector<std::vector<int>>& catvals) {
  
  std::vector<std::vector<int>> selcols = getCols(lv, nvar);
  
  // Create all combinations of category indices (1-indexed like R)
  std::vector<std::pair<int, int>> cat_combs;
  int var1_idx = pair.var1 - 1;  // Convert to 0-based
  int var2_idx = pair.var2 - 1;  // Convert to 0-based
  
  for (int i = 1; i <= lv[var1_idx]; i++) {
    for (int j = 1; j <= lv[var2_idx]; j++) {
      cat_combs.push_back(std::make_pair(i, j));
    }
  }
  
  const std::vector<int>& vals_var1 = catvals[var1_idx];
  const std::vector<int>& vals_var2 = catvals[var2_idx];
  
  // Get threshold indices (convert from 1-based to 0-based)
  int wth1_start = selcols[var1_idx][0] - 1;
  int wth1_end = selcols[var1_idx][1] - 1;
  int wth2_start = selcols[var2_idx][0] - 1;
  int wth2_end = selcols[var2_idx][1] - 1;
  
  // Create thresholds with -Inf and Inf at boundaries
  std::vector<double> th_var1, th_var2;
  
  // Add -Inf
  th_var1.push_back(-std::numeric_limits<double>::infinity());
  th_var2.push_back(-std::numeric_limits<double>::infinity());
  
  // Add actual thresholds
  for (int i = wth1_start; i <= wth1_end; i++) {
    th_var1.push_back(th[i]);
  }
  for (int i = wth2_start; i <= wth2_end; i++) {
    th_var2.push_back(th[i]);
  }
  
  // Add Inf
  th_var1.push_back(std::numeric_limits<double>::infinity());
  th_var2.push_back(std::numeric_limits<double>::infinity());
  
  double mu_joint = 0.0;
  
  for (const auto& comb : cat_combs) {
    int i = comb.first;  // 1-indexed category for var1
    int j = comb.second; // 1-indexed category for var2
    
    // Calculate probability for this category combination
    double term1 = pbivnorm__rcpp_wls0(th_var1[i], th_var2[j], pair.correlation);
    double term2 = pbivnorm__rcpp_wls0(th_var1[i-1], th_var2[j], pair.correlation) * -1.0;
    double term3 = pbivnorm__rcpp_wls0(th_var1[i], th_var2[j-1], pair.correlation) * -1.0;
    double term4 = pbivnorm__rcpp_wls0(th_var1[i-1], th_var2[j-1], pair.correlation);
    
    double p_katkat = term1 + term2 + term3 + term4;
    
    // Add weighted contribution
    mu_joint += vals_var1[i-1] * vals_var2[j-1] * p_katkat;
  }
  
  return mu_joint;
}


// [[Rcpp::export]]
NumericVector apply_get_joint_exp(NumericMatrix combs,
                                    NumericVector th,
                                    IntegerVector lv,
                                    int nvar,
                                    List catvals) {
  
  int ncols = combs.ncol();
  NumericVector results(ncols);
  
  // Convert to easier-to-use formats
  std::vector<double> th_vec = as<std::vector<double>>(th);
  std::vector<int> lv_vec = as<std::vector<int>>(lv);
  
  // Convert catvals
  std::vector<std::vector<int>> catvals_vec;
  for (int i = 0; i < catvals.size(); i++) {
    IntegerVector temp = catvals[i];
    catvals_vec.push_back(as<std::vector<int>>(temp));
  }
  
  // Pre-compute selcols once
  std::vector<std::vector<int>> selcols = getCols(lv_vec, nvar);
  
  for (int col = 0; col < ncols; col++) {
    // Extract column values
    double var1_idx = combs(0, col);
    double var2_idx = combs(1, col);
    double correlation = combs(2, col);
    
    VariablePair pair;
    pair.var1 = static_cast<int>(var1_idx);
    pair.var2 = static_cast<int>(var2_idx);
    pair.correlation = correlation;
    
    // Call get_joint_exp for this column
    results[col] = get_joint_exp(pair, th_vec, lv_vec, nvar, catvals_vec);
  }
  
  return results;
}