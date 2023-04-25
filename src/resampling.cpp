#include <Rcpp.h>

using namespace Rcpp;
using namespace std;  

// [[Rcpp::export]]
IntegerVector multinomial_resampling(const NumericVector & weights, int ndraws){
  NumericVector cumsumw = cumsum(weights);
  IntegerVector ancestors(ndraws);
  NumericVector rand = runif(ndraws);
  int i;
  for (int n = 0; n < ndraws; n++){
    i = 0;
    while (cumsumw(i) < rand(n)){
      i = i + 1;
    }
    ancestors(n) = i + 1;
  }
  return ancestors;
}

// [[Rcpp::export]]
IntegerVector systematic_resampling(const NumericVector & weights, int ndraws){
  IntegerVector ancestors(ndraws);
  double rand = runif(1)[0];
  rand = rand / ndraws;
  int j = 0;
  double csw = weights(0);
  for(int k = 0; k < ndraws; k++){
    while(csw < rand){
      j++;
      csw += weights(j);
    }
    rand = rand + 1. / ndraws;
    ancestors(k) = j + 1;
  }
  return ancestors;
}