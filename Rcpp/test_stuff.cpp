// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
// [[Rcpp::interfaces(r, cpp)]]

using namespace Rcpp;
using namespace arma;


double pbv_rcpp_pnorm0( double z)
{
  double y = ::Rf_pnorm5(z, 0.0, 1.0, 1, 0);
  //--- OUTPUT
  return y;
}

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

NumericVector concat_two(NumericVector A, NumericVector B) {
  int nA = A.length();
  int nB = B.length();
  NumericVector C(nA + nB);
  
  // Copy A
  for(int i = 0; i < nA; i++) {
    C[i] = A[i];
  }
  
  // Copy B
  for(int j = 0; j < nB; j++) {
    C[nA + j] = B[j];
  }
  
  return C;
}

std::vector<std::vector<int>> getCols(const IntegerVector& lv, int nvar) {
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


Rcpp::NumericMatrix compute_sigma_hat(const Rcpp::NumericMatrix& lambda,
                            const Rcpp::NumericMatrix& psi,
                            const Rcpp::NumericMatrix& theta) {
  arma::mat lambda_arma = Rcpp::as<arma::mat>(lambda);
  arma::mat psi_arma = Rcpp::as<arma::mat>(psi);	
  arma::mat theta_arma = Rcpp::as<arma::mat>(theta);
  return wrap(lambda_arma * psi_arma * lambda_arma.t() + theta_arma);  
}


Rcpp::NumericVector apply_get_mus(const NumericVector& th, 
                                  const IntegerVector& lv, 
                                  int nvar, 
                                  const List& catvals){
  
  
  Rcpp::NumericVector result(nvar);
  
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
      
      IntegerVector catvals_var = catvals[var];
      double cat_val = catvals_var[i];
      
      if (i==0){
        prob_cat = pbv_rcpp_pnorm0(p_item[0]);
      } else if (i == max_lv_idx){
        prob_cat = pbv_rcpp_pnorm0(p_item.back() * -1);
      } else {
        prob_cat = pbv_rcpp_pnorm0(p_item[i]) - pbv_rcpp_pnorm0(p_item[i-1]);
      }
      
      mu += cat_val * prob_cat; 
    }
    
    result[var] = mu;
  }
  
  return result;
}


Rcpp::List create_combs(int nvar, const Rcpp::NumericMatrix& sigma_mat) {
  if (sigma_mat.nrow() != nvar || sigma_mat.ncol() != nvar) {
    stop("sigma_mat must be a square matrix of size nvar x nvar");
  }
  
  int n_combs = nvar * (nvar - 1) / 2;
  
  // Create vectors for each row
  Rcpp::IntegerVector var1(n_combs);
  Rcpp::IntegerVector var2(n_combs);
  Rcpp::NumericVector corr(n_combs);
  
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
  return Rcpp::List::create(
    Named("var1") = var1,
    Named("var2") = var2,
    Named("correlation") = corr
  );
}


Rcpp::NumericVector apply_get_joint_exp(Rcpp::List combs,
                                        Rcpp::NumericVector th,
                                        Rcpp::IntegerVector lv,
                                        int nvar,
                                        Rcpp::List catvals) {
  
  Rcpp::IntegerVector var1_vec = combs["var1"];
  Rcpp::IntegerVector var2_vec = combs["var2"];
  Rcpp::NumericVector corr_vec = combs["correlation"];
  
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
    Rcpp::IntegerVector temp_catvals = catvals[var];
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
  
  Rcpp::NumericVector results(ncols);
  
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


NumericVector compute_sigma(const NumericVector& joint_exps, 
                            const NumericVector& mus) {
  int n = mus.size();
  
  // Pre-allocate result
  NumericVector sigma = clone(joint_exps);  // Copy joint_exps
  
  // Subtract mu_i * mu_j from each element
  int idx = 0;
  for (int i = 0; i < n; i++) {
    double mu_i = mus[i];
    for (int j = i + 1; j < n; j++) {
      sigma[idx] -= mu_i * mus[j];
      idx++;
    }
  }
  
  return sigma;
}














//******************************************************************************

// [[Rcpp::export]]
NumericVector compute_moments(Rcpp::NumericVector params, Rcpp::List model_context){
  
  //GLIST
  Rcpp::List m_free_idx = model_context["m.free.idx"];
  Rcpp::List x_free_idx = model_context["x.free.idx"];
  Rcpp::List GLIST_template = model_context["GLIST.template"];
  Rcpp::LogicalVector isSymmetric = model_context["isSymmetric"];
  
  
  Rcpp::List GLIST = x2GLIST_withtheta(params,m_free_idx,x_free_idx,
                                       GLIST_template,isSymmetric);
  
  
  //simga_hat
  NumericMatrix lambda = GLIST["lambda"];	
  NumericMatrix psi = GLIST["psi"];
  NumericMatrix theta = GLIST["theta"];
  			   
  Rcpp::NumericMatrix sigma_hat = compute_sigma_hat(lambda,psi,theta);
  
  //th.pr
  NumericMatrix tau = GLIST["tau"];
  NumericVector th = tau;
  NumericVector th_pr = pbv_rcpp_pnorm(tau*-1);
  
  //mus
  IntegerVector lv = model_context["lv"];
  int nvar = model_context["nvar"];
  List catvals = model_context["catvals"];
  NumericVector mus = apply_get_mus(th, lv, nvar,  catvals);
  
  //sigma
  List combs = create_combs(nvar,sigma_hat);
  NumericVector joint_exps = apply_get_joint_exp(combs,th,lv,nvar,catvals);
  NumericVector sigma = compute_sigma(joint_exps,mus);
  
  //concat
  return(concat_two(th_pr,sigma));
}