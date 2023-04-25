#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
NumericMatrix compute_qvalues_dpoibin(const NumericVector & alpha){
  int N = alpha.length();
  NumericMatrix log_qvalues(N+1, N);
  std::fill(log_qvalues.begin(), log_qvalues.end(), R_NegInf);
  
  // initialize the recursion with q(0,N) and q(1,N)
  log_qvalues(0, N-1) = log(1.0 - alpha[N-1]);
  log_qvalues(1, N-1) = log(alpha[N-1]);
  
  // initialize the recursion with q(0,n) for n = 1,...,N-1
  for(int n = N-1; n > 0; n--){
    log_qvalues(0, n-1) = log_qvalues(0, n) + log(1.0 - alpha[n-1]);
  }
  
  // iterate the recursion: q(i,n) = alpha[n] * q(i-1, n+1) + (1 - alpha[n]) * q(i, n+1) 
  double ls1, ls2, maxls;
  for(int n = N - 1; n > 0; n--){
    for (int i = 1; i <= N - n + 1 ; i++){
      ls1 = log(alpha[n-1]) + log_qvalues(i-1, n);
      ls2 = log(1.0 - alpha[n-1]) + log_qvalues(i, n);
      maxls = max(ls1, ls2);
      log_qvalues(i, n-1) = log(exp(ls1 - maxls) + exp(ls2 - maxls)) + maxls;
      if (Rcpp::traits::is_nan<REALSXP>(log_qvalues(i, n-1))){
        log_qvalues(i, n-1) = R_NegInf;
      }
    }
  }
  return log_qvalues;
}

// [[Rcpp::export]]
NumericVector dpoibin_cpp(const NumericVector & eval, const NumericVector & alpha){
  int M = eval.length();
  NumericVector log_pmf(M);
  
  // compute q-values with recursion
  NumericMatrix log_qvalues = compute_qvalues_dpoibin(alpha);
  for(int m = 0; m < M; m++){
    // get Poisson binomial probability 
    log_pmf(m) = log_qvalues(eval[m], 0);
  }
  
  return log_pmf;
}

// [[Rcpp::export]]
NumericMatrix dpoibin_vectorized_cpp(const NumericVector & eval, const NumericMatrix & alpha){
  int M = eval.length();
  int P = alpha.ncol();
  NumericMatrix log_pmf(M, P);
  
  for(int p = 0; p < P; p++){
    log_pmf(_, p) = dpoibin_cpp(eval, alpha(_, p));
  }
  
  return log_pmf;
}
