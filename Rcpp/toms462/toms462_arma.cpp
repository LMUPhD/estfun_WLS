// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

//****************************************************************************80
// Helper functions
//****************************************************************************80

double gauss ( double t ) {
  // Returns the area of the lower tail of the normal curve
  return ( 1.0 + erf ( t / sqrt ( 2.0 ) ) ) / 2.0;
}

double r8_abs ( double x ) {
  // Returns the absolute value
  return ( x < 0.0 ) ? -x : x;
}

double r8_max ( double x, double y ) {
  // Returns the maximum of two values
  return ( y < x ) ? x : y;
}

double r8_min ( double x, double y ) {
  // Returns the minimum of two values
  return ( y < x ) ? y : x;
}

//****************************************************************************80
// Main bivariate normal CDF function
//****************************************************************************80

// [[Rcpp::export]]
double bivnor ( double ah, double ak, double r ) {
  // Purpose: BIVNOR computes the bivariate normal CDF.
  // Discussion: Computes the probability for two normal variates X and Y
  // whose correlation is R, that AH <= X and AK <= Y.
  
  double a2;
  double ap;
  double b;
  double cn;
  double con;
  double conex;
  double ex;
  double g2;
  double gh;
  double gk;
  double gw;
  double h2;
  double h4;
  int i;
  static int idig = 15;
  int is;
  double rr;
  double s1;
  double s2;
  double sgn;
  double sn;
  double sp;
  double sqr;
  double t;
  static double twopi = 6.283185307179587;
  double w2;
  double wh;
  double wk;
  
  b = 0.0;
  
  gh = gauss ( - ah ) / 2.0;
  gk = gauss ( - ak ) / 2.0;
  
  if ( r == 0.0 ) {
    b = 4.00 * gh * gk;
    b = r8_max ( b, 0.0 );
    b = r8_min ( b, 1.0 );
    return b;
  }
  
  rr = ( 1.0 + r ) * ( 1.0 - r );
  
  if ( rr < 0.0 ) {
    Rcpp::stop("BIVNOR - Fatal error! 1 < |R|.");
  }
  
  if ( rr == 0.0 ) {
    if ( r < 0.0 ) {
      if ( ah + ak < 0.0 ) {
        b = 2.0 * ( gh + gk ) - 1.0;
      }
    } else {
      if ( ah - ak < 0.0 ) {
        b = 2.0 * gk;
      } else {
        b = 2.0 * gh;
      }
    }
    b = r8_max ( b, 0.0 );
    b = r8_min ( b, 1.0 );
    return b;
  }
  
  sqr = sqrt ( rr );
  
  if ( idig == 15 ) {
    con = twopi * 1.0E-15 / 2.0;
  } else {
    con = twopi / 2.0;
    for ( i = 1; i <= idig; i++ ) {
      con = con / 10.0;
    }
  }
  
  // (0,0)
  if ( ah == 0.0 && ak == 0.0 ) {
    b = 0.25 + asin ( r ) / twopi;
    b = r8_max ( b, 0.0 );
    b = r8_min ( b, 1.0 );
    return b;
  }
  
  // (0,nonzero)
  if ( ah == 0.0 && ak != 0.0 ) {
    b = gk;
    wh = -ak;
    wk = ( ah / ak - r ) / sqr;
    gw = 2.0 * gk;
    is = 1;
  }
  // (nonzero,0)
  else if ( ah != 0.0 && ak == 0.0 ) {
    b = gh;
    wh = -ah;
    wk = ( ak / ah - r ) / sqr;
    gw = 2.0 * gh;
    is = -1;
  }
  // (nonzero,nonzero)
  else if ( ah != 0.0 && ak != 0.0 ) {
    b = gh + gk;
    if ( ah * ak < 0.0 ) {
      b = b - 0.5;
    }
    wh = - ah;
    wk = ( ak / ah - r ) / sqr;
    gw = 2.0 * gh;
    is = -1;
  }
  
  for ( ; ; ) {
    sgn = -1.0;
    t = 0.0;
    
    if ( wk != 0.0 ) {
      if ( r8_abs ( wk ) == 1.0 ) {
        t = wk * gw * ( 1.0 - gw ) / 2.0;
        b = b + sgn * t;
      } else {
        if ( 1.0 < r8_abs ( wk ) ) {
          sgn = -sgn;
          wh = wh * wk;
          g2 = gauss ( wh );
          wk = 1.0 / wk;
          
          if ( wk < 0.0 ) {
            b = b + 0.5;
          }
          b = b - ( gw + g2 ) / 2.0 + gw * g2;
        }
        h2 = wh * wh;
        a2 = wk * wk;
        h4 = h2 / 2.0;
        ex = exp ( - h4 );
        w2 = h4 * ex;
        ap = 1.0;
        s2 = ap - ex;
        sp = ap;
        s1 = 0.0;
        sn = s1;
        conex = r8_abs ( con / wk );
        
        for ( ; ; ) {
          cn = ap * s2 / ( sn + sp );
          s1 = s1 + cn;
          
          if ( r8_abs ( cn ) <= conex ) {
            break;
          }
          sn = sp;
          sp = sp + 1.0;
          s2 = s2 - w2;
          w2 = w2 * h4 / sp;
          ap = - ap * a2;
        }
        t = ( atan ( wk ) - wk * s1 ) / twopi;
        b = b + sgn * t;
      }
    }
    if ( 0 <= is ) {
      break;
    }
    if ( ak == 0.0 ) {
      break;
    }
    wh = -ak;
    wk = ( ah / ak - r ) / sqr;
    gw = 2.0 * gk;
    is = 1;
  }
  
  b = r8_max ( b, 0.0 );
  b = r8_min ( b, 1.0 );
  
  return b;
}

//****************************************************************************80
// Vectorized version for R compatibility
//****************************************************************************80

// [[Rcpp::export]]
NumericVector pbivnor(NumericVector x, NumericVector y, NumericVector r) {
  // Vectorized version that accepts vectors and recycles parameters
  int n = std::max(x.size(), std::max(y.size(), r.size()));
  NumericVector result(n);
  
  for(int i = 0; i < n; i++) {
    double xi = x[i % x.size()];
    double yi = y[i % y.size()];
    double ri = r[i % r.size()];
    
    result[i] = bivnor(xi, yi, ri);
  }
  
  return result;
}