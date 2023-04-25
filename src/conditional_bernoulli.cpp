#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
NumericMatrix compute_qvalues_condber(const int & I, const NumericVector & alpha){
  int N = alpha.length();
  NumericMatrix log_qvalues(I+1, N);
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
    for(int i = 1; (i <= N - n + 1) && (i <= I); i++){ // i = 0,1,...,I
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
LogicalVector rcondber_cpp(const int & I,
                           const NumericVector & alpha){
  int N = alpha.length();
  LogicalVector X(N);
  
  if (I == 0){
    for(int n = 0; n < N; n++){
      X[n] = false;
    }
  } else if (I == N){
    for(int n = 0; n < N; n++){
      X[n] = true;
    }
  } else {
    int partial_I = 0; 
    double log_prob;
    
    // compute q-values with recursion
    NumericMatrix log_qvalues = compute_qvalues_condber(I, alpha); // (I+1) x N
    
    // sample from conditional Bernoulli distribution sequentially (see id-checking method in Chen and Liu (1997))
    NumericVector rand = runif(N-1);
    for(int n = 0; n < (N-1); n++){
      log_prob = log(alpha[n]) + log_qvalues(I-partial_I-1, n+1) - log_qvalues(I-partial_I, n);
      if(log(rand[n]) < log_prob){
        X[n] = true;
        partial_I++;
      } else{
        X[n] = false;
      }
      if(partial_I == I) break;
    }
    X[N-1] = (partial_I == (I - 1));
  }
  
  return X;
}

// [[Rcpp::export]]
LogicalMatrix rcondber_vectorized_cpp(const IntegerVector & I,
                                      const NumericMatrix & alpha){
  int N = alpha.nrow();
  int P = alpha.ncol();
  
  LogicalMatrix X(N, P);
  
  for(int p = 0; p < P; p++){
    X(_, p) = rcondber_cpp(I[p], alpha(_, p));
  }
  
  return X;
}

// sample uniformly in {0,...,n}
inline int randWrapper(const int n) { return floor(unif_rand()*n); }

// [[Rcpp::export]]
LogicalVector rcondber_mcmc_cpp(const int & I,
                                const NumericVector & alpha, 
                                const int & mcmc_iterations){
  int N = alpha.length();
  vector<int> s0; // set of indices n such that x_n = 0 
  vector<int> s1; // set of indices n such that x_n = 1 
  int i0, i1, tmp;
  double log_odds_i0, log_odds_i1, log_uniform, log_odds_diff;
  
  // initialization
  GetRNGstate();
  s0.clear(); s1.clear(); // reset s0, s1
  s0.reserve(N); s1.reserve(N);
  IntegerVector indices = seq(0, N-1);
  std::random_shuffle(indices.begin(), indices.end(), randWrapper); // random permutation of indices in 0,...,N-1
  for (int n = 0; n < N; n++){
    if (n < I){
      s1.push_back(indices(n)); // first I shuffled indices are set to one in vector x
    } else {
      s0.push_back(indices(n));
    }
  }
  PutRNGstate();
  
  // MCMC iterations
  GetRNGstate();
  for (int i = 0; i < mcmc_iterations; i++){
    // sample i0 in s0 and i1 in s1
    i0 = randWrapper(N-I);
    i1 = randWrapper(I);
    
    // log accept ratio
    log_odds_i0 = log(alpha[s0.at(i0)] / (1.0 - alpha[s0.at(i0)]));
    log_odds_i1 = log(alpha[s1.at(i1)] / (1.0 - alpha[s1.at(i1)]));
    log_odds_diff = log_odds_i0 - log_odds_i1;
    log_uniform = log(unif_rand());
    
    // comparison done on log scale
    if (log_uniform < log_odds_diff){ // accept and swap elements
      tmp = s0.at(i0);
      s0.at(i0) = s1.at(i1);
      s1.at(i1) = tmp;
    } 
  }
  PutRNGstate();
  
  // output vector X
  LogicalVector X = rep(false, N);
  for (int n = 0; n < s1.size(); n++){
    X[s1.at(n)] = true;
  }
  return X;
}

// [[Rcpp::export]]
LogicalMatrix rcondber_mcmc_vectorized_cpp(const IntegerVector & I,
                                           const NumericMatrix & alpha, 
                                           const int & mcmc_iterations){
  int N = alpha.nrow();
  int P = alpha.ncol();
  
  LogicalMatrix X(N, P);
  
  for(int p = 0; p < P; p++){
    X(_, p) = rcondber_mcmc_cpp(I[p], alpha(_, p), mcmc_iterations);
  }
  
  return X;
}
