#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
NumericVector dsumbin_cpp(const NumericVector & eval, 
                          const int & n1, 
                          const double & p1,
                          const int & n2, 
                          const double & p2){
  int M = eval.length();
  NumericVector log_pmf(M);
  double maxlogsummand;
  for (int m = 0; m < M; m++){
    IntegerVector eval_seq = seq(0, eval[m]);
    NumericVector log_binom1 = dbinom(eval_seq, n1, p1, true);
    NumericVector log_binom2 = dbinom(rev(eval_seq), n2, p2, true);
    NumericVector log_summand = log_binom1 + log_binom2;
    maxlogsummand = max(log_summand);
    NumericVector summand = exp(log_summand - maxlogsummand);
    log_pmf(m) = log(sum(summand)) + maxlogsummand;
  }

  return log_pmf;
}

// [[Rcpp::export]]
NumericMatrix dsumbin_vectorized_cpp(const NumericVector & eval, 
                                     const IntegerVector & n1, 
                                     const NumericVector & p1,
                                     const IntegerVector & n2, 
                                     const double & p2){
  
  int M = eval.length();
  int K = n1.length();
  NumericMatrix log_pmf(M, K);
  
  for (int k = 0; k < K; k++){
    log_pmf(_, k) = dsumbin_cpp(eval, n1[k], p1[k], n2[k], p2);
  }
  
  return log_pmf;
}
