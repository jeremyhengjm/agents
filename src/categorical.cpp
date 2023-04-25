#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
NumericVector rcat_cpp(const NumericVector & support, 
                       const NumericVector & prob, 
                       const int & num_samples){
  NumericVector I = sample(support, num_samples, true, prob);
  return I;
}

// [[Rcpp::export]]
NumericVector rcat_vectorized_cpp(const NumericVector & support, 
                                  const NumericMatrix & prob){
  int M = prob.nrow();
  int P = prob.ncol();
  NumericVector current_prob(M);
  NumericVector I(P);
  for (int p = 0; p < P; p++){
    current_prob = prob(_, p);
    I[p] = sample(support, 1, true, current_prob)[0];
  }
  return I;
}
