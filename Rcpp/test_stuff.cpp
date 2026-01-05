#include <Rcpp.h>
using namespace Rcpp;

double pbv_rcpp_pnorm0( double z)
{
  double y = ::Rf_pnorm5(z, 0.0, 1.0, 1, 0);
  //--- OUTPUT
  return y;
}

//getCols with different containers...
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






// translate get_mus to C++
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
    int wth_start = selcols[var][0] - 1;
    int wth_end = selcols[var][1] - 1;
    
    //Get thresholds
    std::vector<double> p_item;
    for (int i = wth_start; i <= wth_end; i++){
      p_item.push_back(th[i]);
    }
    
    //get probas over categories
    // Get max category index from 1-based to 0-based
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
