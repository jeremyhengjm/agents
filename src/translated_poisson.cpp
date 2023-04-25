#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector dtranspoi_poibin_approx_cpp(const NumericVector & eval, const NumericVector & alpha){
  int M = eval.length();
  int N = alpha.length();
  NumericVector log_pmf(M);
  double mu = 0;
  double sigmasq = 0;
  
  for(int n = 0; n < N; n++){
    mu += alpha[n];
    sigmasq += alpha[n] * (1.0 - alpha[n]);
  }
  double difference = mu - sigmasq;
  double diff_floor = floor(difference);
  double diff_frac = difference - diff_floor;
  double lambda;
  
  for(int m = 0; m < M; m++){
    if ((eval[m] >= diff_floor) && (eval[m] <= N)){
      // approximate with Poisson probability 
      lambda = sigmasq + diff_frac;
      log_pmf(m) = R::dpois(eval[m] - diff_floor, lambda, true);
    } else{
      log_pmf(m) = R_NegInf;
    }
  }
  
  return log_pmf;
}

// [[Rcpp::export]]
NumericMatrix dtranspoi_poibin_approx_vectorized_cpp(const NumericVector & eval, const NumericMatrix & alpha){
  int M = eval.length();
  int P = alpha.ncol();
  NumericMatrix log_pmf(M, P);
  
  for(int p = 0; p < P; p++){
    log_pmf(_, p) = dtranspoi_poibin_approx_cpp(eval, alpha(_, p));
  }
  
  return log_pmf;
}  

// [[Rcpp::export]]
NumericVector dtranspoi_cpp(const NumericVector & eval, const double & mu, const double & sigmasq, const double & end_support){
  int M = eval.length();
  NumericVector log_pmf(M);

  double difference = mu - sigmasq;
  double diff_floor = floor(difference);
  double diff_frac = difference - diff_floor;
  double lambda = sigmasq + diff_frac;

  for(int m = 0; m < M; m++){
    if ((eval[m] >= diff_floor) && (eval[m] <= end_support)){
      // approximate with Poisson probability
      log_pmf(m) = R::dpois(eval[m] - diff_floor, lambda, true);
    } else{
      log_pmf(m) = R_NegInf;
    }
  }

  return log_pmf;
}

// [[Rcpp::export]]
NumericMatrix dtranspoi_vectorized_cpp(const NumericVector & eval, const NumericVector & mu, const NumericVector & sigmasq, const double & end_support){
  int M = eval.length();
  int P = mu.length();
  NumericMatrix log_pmf(M, P);

  for(int p = 0; p < P; p++){
    log_pmf(_, p) = dtranspoi_cpp(eval, mu[p], sigmasq[p], end_support);
  }

  return log_pmf;
}
