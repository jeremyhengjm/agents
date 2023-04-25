#include <Rcpp.h>

using namespace Rcpp;
using namespace std;  

// [[Rcpp::export]]
NumericMatrix sis_alpha_full(const LogicalMatrix & X, 
                             const NumericVector & I,
                             const NumericVector & lambda, 
                             const NumericVector & gamma){
  int N = X.nrow();
  int P = X.ncol();
  NumericMatrix alpha(N, P);
  for (int p = 0; p < P; p++){
    for (int n = 0; n < N; n++){
      alpha(n,p) = X(n,p) * (1.0 - gamma(n)) + (1.0 - X(n,p)) * lambda(n) * I(p) / N; 
    }
  }
  return alpha;
}

// [[Rcpp::export]]
NumericMatrix sis_alpha_net(const LogicalMatrix & X, 
                            const NumericVector & lambda, 
                            const NumericVector & gamma,
                            const NumericMatrix & adjacency,
                            const NumericVector & degree){
  int N = X.nrow();
  int P = X.ncol();
  double infected;
  NumericMatrix alpha(N, P);
  
  // first compute no. of infected neighbors
  NumericMatrix I(N, P); 
  for (int p = 0; p < P; p++){
    for (int n = 0; n < N; n++){
      infected = 0;
      for (int m = 0; m < N; m++){
        infected += X(m,p) * adjacency(n, m);
      }
      I(n,p) = infected;
    }
  }
  
  // then compute transition probabilities
  for (int p = 0; p < P; p++){
    for (int n = 0; n < N; n++){
      alpha(n,p) = X(n,p) * (1.0 - gamma(n)) + (1.0 - X(n,p)) * lambda(n) * I(n,p) / degree(n); 
    }
  }
  return alpha;
}

// [[Rcpp::export]]
List sis_compute_infected_prob(const NumericMatrix & log_poibin, const NumericVector & log_policy){
  int M = log_poibin.nrow();
  int P = log_poibin.ncol();
  
  NumericMatrix normprob(M, P);
  NumericVector log_normconstant(P);
  NumericVector log_prob(M);
  NumericVector prob(M);
  double max_log_prob;
  double sum_prob;
  
  for (int p = 0; p < P; p++){
    for (int m = 0; m < M; m++){
      log_prob[m] = log_poibin(m,p) + log_policy[m];
    }
    max_log_prob = max(log_prob);
    sum_prob = 0;
    for (int m = 0; m < M; m++){
      prob[m] = exp(log_prob[m] - max_log_prob);
      sum_prob += prob[m];
    }
    for (int m = 0; m < M; m++){
      normprob(m,p) = prob[m] / sum_prob;
    }
    log_normconstant[p] = log(sum_prob) + max_log_prob;
  }
  
  return List::create(Named("normprob") = normprob,
                      Named("log_normconstant") = log_normconstant);
}

// [[Rcpp::export]]
NumericVector sis_approximate_infected_prob(const NumericMatrix & log_sumbin, const NumericVector & log_policy){
  int M = log_sumbin.nrow();
  int P = log_sumbin.ncol();
  
  NumericVector log_normconstant(P);
  NumericVector log_prob(M);
  NumericVector prob(M);
  double max_log_prob;
  double sum_prob;
  
  for (int p = 0; p < P; p++){
    for (int m = 0; m < M; m++){
      log_prob[m] = log_sumbin(m,p) + log_policy[m];
    }
    max_log_prob = max(log_prob);
    sum_prob = 0;
    for (int m = 0; m < M; m++){
      prob[m] = exp(log_prob[m] - max_log_prob);
      sum_prob += prob[m];
    }
    log_normconstant[p] = log(sum_prob) + max_log_prob;
  }
  
  return log_normconstant;
}
