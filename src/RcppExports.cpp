// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// rcat_cpp
NumericVector rcat_cpp(const NumericVector& support, const NumericVector& prob, const int& num_samples);
RcppExport SEXP _agents_rcat_cpp(SEXP supportSEXP, SEXP probSEXP, SEXP num_samplesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type support(supportSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type prob(probSEXP);
    Rcpp::traits::input_parameter< const int& >::type num_samples(num_samplesSEXP);
    rcpp_result_gen = Rcpp::wrap(rcat_cpp(support, prob, num_samples));
    return rcpp_result_gen;
END_RCPP
}
// rcat_vectorized_cpp
NumericVector rcat_vectorized_cpp(const NumericVector& support, const NumericMatrix& prob);
RcppExport SEXP _agents_rcat_vectorized_cpp(SEXP supportSEXP, SEXP probSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type support(supportSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type prob(probSEXP);
    rcpp_result_gen = Rcpp::wrap(rcat_vectorized_cpp(support, prob));
    return rcpp_result_gen;
END_RCPP
}
// compute_qvalues_condber
NumericMatrix compute_qvalues_condber(const int& I, const NumericVector& alpha);
RcppExport SEXP _agents_compute_qvalues_condber(SEXP ISEXP, SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int& >::type I(ISEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type alpha(alphaSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_qvalues_condber(I, alpha));
    return rcpp_result_gen;
END_RCPP
}
// rcondber_cpp
LogicalVector rcondber_cpp(const int& I, const NumericVector& alpha);
RcppExport SEXP _agents_rcondber_cpp(SEXP ISEXP, SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int& >::type I(ISEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type alpha(alphaSEXP);
    rcpp_result_gen = Rcpp::wrap(rcondber_cpp(I, alpha));
    return rcpp_result_gen;
END_RCPP
}
// rcondber_vectorized_cpp
LogicalMatrix rcondber_vectorized_cpp(const IntegerVector& I, const NumericMatrix& alpha);
RcppExport SEXP _agents_rcondber_vectorized_cpp(SEXP ISEXP, SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const IntegerVector& >::type I(ISEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type alpha(alphaSEXP);
    rcpp_result_gen = Rcpp::wrap(rcondber_vectorized_cpp(I, alpha));
    return rcpp_result_gen;
END_RCPP
}
// rcondber_mcmc_cpp
LogicalVector rcondber_mcmc_cpp(const int& I, const NumericVector& alpha, const int& mcmc_iterations);
RcppExport SEXP _agents_rcondber_mcmc_cpp(SEXP ISEXP, SEXP alphaSEXP, SEXP mcmc_iterationsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int& >::type I(ISEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const int& >::type mcmc_iterations(mcmc_iterationsSEXP);
    rcpp_result_gen = Rcpp::wrap(rcondber_mcmc_cpp(I, alpha, mcmc_iterations));
    return rcpp_result_gen;
END_RCPP
}
// rcondber_mcmc_vectorized_cpp
LogicalMatrix rcondber_mcmc_vectorized_cpp(const IntegerVector& I, const NumericMatrix& alpha, const int& mcmc_iterations);
RcppExport SEXP _agents_rcondber_mcmc_vectorized_cpp(SEXP ISEXP, SEXP alphaSEXP, SEXP mcmc_iterationsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const IntegerVector& >::type I(ISEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const int& >::type mcmc_iterations(mcmc_iterationsSEXP);
    rcpp_result_gen = Rcpp::wrap(rcondber_mcmc_vectorized_cpp(I, alpha, mcmc_iterations));
    return rcpp_result_gen;
END_RCPP
}
// compute_qvalues_dpoibin
NumericMatrix compute_qvalues_dpoibin(const NumericVector& alpha);
RcppExport SEXP _agents_compute_qvalues_dpoibin(SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type alpha(alphaSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_qvalues_dpoibin(alpha));
    return rcpp_result_gen;
END_RCPP
}
// dpoibin_cpp
NumericVector dpoibin_cpp(const NumericVector& eval, const NumericVector& alpha);
RcppExport SEXP _agents_dpoibin_cpp(SEXP evalSEXP, SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type eval(evalSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type alpha(alphaSEXP);
    rcpp_result_gen = Rcpp::wrap(dpoibin_cpp(eval, alpha));
    return rcpp_result_gen;
END_RCPP
}
// dpoibin_vectorized_cpp
NumericMatrix dpoibin_vectorized_cpp(const NumericVector& eval, const NumericMatrix& alpha);
RcppExport SEXP _agents_dpoibin_vectorized_cpp(SEXP evalSEXP, SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type eval(evalSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type alpha(alphaSEXP);
    rcpp_result_gen = Rcpp::wrap(dpoibin_vectorized_cpp(eval, alpha));
    return rcpp_result_gen;
END_RCPP
}
// multinomial_resampling
IntegerVector multinomial_resampling(const NumericVector& weights, int ndraws);
RcppExport SEXP _agents_multinomial_resampling(SEXP weightsSEXP, SEXP ndrawsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< int >::type ndraws(ndrawsSEXP);
    rcpp_result_gen = Rcpp::wrap(multinomial_resampling(weights, ndraws));
    return rcpp_result_gen;
END_RCPP
}
// systematic_resampling
IntegerVector systematic_resampling(const NumericVector& weights, int ndraws);
RcppExport SEXP _agents_systematic_resampling(SEXP weightsSEXP, SEXP ndrawsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< int >::type ndraws(ndrawsSEXP);
    rcpp_result_gen = Rcpp::wrap(systematic_resampling(weights, ndraws));
    return rcpp_result_gen;
END_RCPP
}
// sis_alpha_full
NumericMatrix sis_alpha_full(const LogicalMatrix& X, const NumericVector& I, const NumericVector& lambda, const NumericVector& gamma);
RcppExport SEXP _agents_sis_alpha_full(SEXP XSEXP, SEXP ISEXP, SEXP lambdaSEXP, SEXP gammaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const LogicalMatrix& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type I(ISEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type gamma(gammaSEXP);
    rcpp_result_gen = Rcpp::wrap(sis_alpha_full(X, I, lambda, gamma));
    return rcpp_result_gen;
END_RCPP
}
// sis_alpha_net
NumericMatrix sis_alpha_net(const LogicalMatrix& X, const NumericVector& lambda, const NumericVector& gamma, const NumericMatrix& adjacency, const NumericVector& degree);
RcppExport SEXP _agents_sis_alpha_net(SEXP XSEXP, SEXP lambdaSEXP, SEXP gammaSEXP, SEXP adjacencySEXP, SEXP degreeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const LogicalMatrix& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type adjacency(adjacencySEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type degree(degreeSEXP);
    rcpp_result_gen = Rcpp::wrap(sis_alpha_net(X, lambda, gamma, adjacency, degree));
    return rcpp_result_gen;
END_RCPP
}
// sis_compute_infected_prob
List sis_compute_infected_prob(const NumericMatrix& log_poibin, const NumericVector& log_policy);
RcppExport SEXP _agents_sis_compute_infected_prob(SEXP log_poibinSEXP, SEXP log_policySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type log_poibin(log_poibinSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type log_policy(log_policySEXP);
    rcpp_result_gen = Rcpp::wrap(sis_compute_infected_prob(log_poibin, log_policy));
    return rcpp_result_gen;
END_RCPP
}
// sis_approximate_infected_prob
NumericVector sis_approximate_infected_prob(const NumericMatrix& log_sumbin, const NumericVector& log_policy);
RcppExport SEXP _agents_sis_approximate_infected_prob(SEXP log_sumbinSEXP, SEXP log_policySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type log_sumbin(log_sumbinSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type log_policy(log_policySEXP);
    rcpp_result_gen = Rcpp::wrap(sis_approximate_infected_prob(log_sumbin, log_policy));
    return rcpp_result_gen;
END_RCPP
}
// dsumbin_cpp
NumericVector dsumbin_cpp(const NumericVector& eval, const int& n1, const double& p1, const int& n2, const double& p2);
RcppExport SEXP _agents_dsumbin_cpp(SEXP evalSEXP, SEXP n1SEXP, SEXP p1SEXP, SEXP n2SEXP, SEXP p2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type eval(evalSEXP);
    Rcpp::traits::input_parameter< const int& >::type n1(n1SEXP);
    Rcpp::traits::input_parameter< const double& >::type p1(p1SEXP);
    Rcpp::traits::input_parameter< const int& >::type n2(n2SEXP);
    Rcpp::traits::input_parameter< const double& >::type p2(p2SEXP);
    rcpp_result_gen = Rcpp::wrap(dsumbin_cpp(eval, n1, p1, n2, p2));
    return rcpp_result_gen;
END_RCPP
}
// dsumbin_vectorized_cpp
NumericMatrix dsumbin_vectorized_cpp(const NumericVector& eval, const IntegerVector& n1, const NumericVector& p1, const IntegerVector& n2, const double& p2);
RcppExport SEXP _agents_dsumbin_vectorized_cpp(SEXP evalSEXP, SEXP n1SEXP, SEXP p1SEXP, SEXP n2SEXP, SEXP p2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type eval(evalSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type n1(n1SEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type p1(p1SEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type n2(n2SEXP);
    Rcpp::traits::input_parameter< const double& >::type p2(p2SEXP);
    rcpp_result_gen = Rcpp::wrap(dsumbin_vectorized_cpp(eval, n1, p1, n2, p2));
    return rcpp_result_gen;
END_RCPP
}
// dtranspoi_poibin_approx_cpp
NumericVector dtranspoi_poibin_approx_cpp(const NumericVector& eval, const NumericVector& alpha);
RcppExport SEXP _agents_dtranspoi_poibin_approx_cpp(SEXP evalSEXP, SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type eval(evalSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type alpha(alphaSEXP);
    rcpp_result_gen = Rcpp::wrap(dtranspoi_poibin_approx_cpp(eval, alpha));
    return rcpp_result_gen;
END_RCPP
}
// dtranspoi_poibin_approx_vectorized_cpp
NumericMatrix dtranspoi_poibin_approx_vectorized_cpp(const NumericVector& eval, const NumericMatrix& alpha);
RcppExport SEXP _agents_dtranspoi_poibin_approx_vectorized_cpp(SEXP evalSEXP, SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type eval(evalSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type alpha(alphaSEXP);
    rcpp_result_gen = Rcpp::wrap(dtranspoi_poibin_approx_vectorized_cpp(eval, alpha));
    return rcpp_result_gen;
END_RCPP
}
// dtranspoi_cpp
NumericVector dtranspoi_cpp(const NumericVector& eval, const double& mu, const double& sigmasq, const double& end_support);
RcppExport SEXP _agents_dtranspoi_cpp(SEXP evalSEXP, SEXP muSEXP, SEXP sigmasqSEXP, SEXP end_supportSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type eval(evalSEXP);
    Rcpp::traits::input_parameter< const double& >::type mu(muSEXP);
    Rcpp::traits::input_parameter< const double& >::type sigmasq(sigmasqSEXP);
    Rcpp::traits::input_parameter< const double& >::type end_support(end_supportSEXP);
    rcpp_result_gen = Rcpp::wrap(dtranspoi_cpp(eval, mu, sigmasq, end_support));
    return rcpp_result_gen;
END_RCPP
}
// dtranspoi_vectorized_cpp
NumericMatrix dtranspoi_vectorized_cpp(const NumericVector& eval, const NumericVector& mu, const NumericVector& sigmasq, const double& end_support);
RcppExport SEXP _agents_dtranspoi_vectorized_cpp(SEXP evalSEXP, SEXP muSEXP, SEXP sigmasqSEXP, SEXP end_supportSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type eval(evalSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type mu(muSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type sigmasq(sigmasqSEXP);
    Rcpp::traits::input_parameter< const double& >::type end_support(end_supportSEXP);
    rcpp_result_gen = Rcpp::wrap(dtranspoi_vectorized_cpp(eval, mu, sigmasq, end_support));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_agents_rcat_cpp", (DL_FUNC) &_agents_rcat_cpp, 3},
    {"_agents_rcat_vectorized_cpp", (DL_FUNC) &_agents_rcat_vectorized_cpp, 2},
    {"_agents_compute_qvalues_condber", (DL_FUNC) &_agents_compute_qvalues_condber, 2},
    {"_agents_rcondber_cpp", (DL_FUNC) &_agents_rcondber_cpp, 2},
    {"_agents_rcondber_vectorized_cpp", (DL_FUNC) &_agents_rcondber_vectorized_cpp, 2},
    {"_agents_rcondber_mcmc_cpp", (DL_FUNC) &_agents_rcondber_mcmc_cpp, 3},
    {"_agents_rcondber_mcmc_vectorized_cpp", (DL_FUNC) &_agents_rcondber_mcmc_vectorized_cpp, 3},
    {"_agents_compute_qvalues_dpoibin", (DL_FUNC) &_agents_compute_qvalues_dpoibin, 1},
    {"_agents_dpoibin_cpp", (DL_FUNC) &_agents_dpoibin_cpp, 2},
    {"_agents_dpoibin_vectorized_cpp", (DL_FUNC) &_agents_dpoibin_vectorized_cpp, 2},
    {"_agents_multinomial_resampling", (DL_FUNC) &_agents_multinomial_resampling, 2},
    {"_agents_systematic_resampling", (DL_FUNC) &_agents_systematic_resampling, 2},
    {"_agents_sis_alpha_full", (DL_FUNC) &_agents_sis_alpha_full, 4},
    {"_agents_sis_alpha_net", (DL_FUNC) &_agents_sis_alpha_net, 5},
    {"_agents_sis_compute_infected_prob", (DL_FUNC) &_agents_sis_compute_infected_prob, 2},
    {"_agents_sis_approximate_infected_prob", (DL_FUNC) &_agents_sis_approximate_infected_prob, 2},
    {"_agents_dsumbin_cpp", (DL_FUNC) &_agents_dsumbin_cpp, 5},
    {"_agents_dsumbin_vectorized_cpp", (DL_FUNC) &_agents_dsumbin_vectorized_cpp, 5},
    {"_agents_dtranspoi_poibin_approx_cpp", (DL_FUNC) &_agents_dtranspoi_poibin_approx_cpp, 2},
    {"_agents_dtranspoi_poibin_approx_vectorized_cpp", (DL_FUNC) &_agents_dtranspoi_poibin_approx_vectorized_cpp, 2},
    {"_agents_dtranspoi_cpp", (DL_FUNC) &_agents_dtranspoi_cpp, 4},
    {"_agents_dtranspoi_vectorized_cpp", (DL_FUNC) &_agents_dtranspoi_vectorized_cpp, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_agents(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
