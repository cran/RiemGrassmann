// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// cpp_prangles
arma::vec cpp_prangles(arma::mat U, arma::mat V);
RcppExport SEXP _RiemGrassmann_cpp_prangles(SEXP USEXP, SEXP VSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type U(USEXP);
    Rcpp::traits::input_parameter< arma::mat >::type V(VSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_prangles(U, V));
    return rcpp_result_gen;
END_RCPP
}
// cpp_pairdist
double cpp_pairdist(arma::mat U, arma::mat V, std::string dname);
RcppExport SEXP _RiemGrassmann_cpp_pairdist(SEXP USEXP, SEXP VSEXP, SEXP dnameSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type U(USEXP);
    Rcpp::traits::input_parameter< arma::mat >::type V(VSEXP);
    Rcpp::traits::input_parameter< std::string >::type dname(dnameSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_pairdist(U, V, dname));
    return rcpp_result_gen;
END_RCPP
}
// cpp_pdist
arma::mat cpp_pdist(arma::field<arma::mat> INPUT1, std::string dname);
RcppExport SEXP _RiemGrassmann_cpp_pdist(SEXP INPUT1SEXP, SEXP dnameSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::field<arma::mat> >::type INPUT1(INPUT1SEXP);
    Rcpp::traits::input_parameter< std::string >::type dname(dnameSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_pdist(INPUT1, dname));
    return rcpp_result_gen;
END_RCPP
}
// cpp_pdist2
arma::mat cpp_pdist2(arma::field<arma::mat> INPUT1, arma::field<arma::mat> INPUT2, std::string dname);
RcppExport SEXP _RiemGrassmann_cpp_pdist2(SEXP INPUT1SEXP, SEXP INPUT2SEXP, SEXP dnameSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::field<arma::mat> >::type INPUT1(INPUT1SEXP);
    Rcpp::traits::input_parameter< arma::field<arma::mat> >::type INPUT2(INPUT2SEXP);
    Rcpp::traits::input_parameter< std::string >::type dname(dnameSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_pdist2(INPUT1, INPUT2, dname));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_RiemGrassmann_cpp_prangles", (DL_FUNC) &_RiemGrassmann_cpp_prangles, 2},
    {"_RiemGrassmann_cpp_pairdist", (DL_FUNC) &_RiemGrassmann_cpp_pairdist, 3},
    {"_RiemGrassmann_cpp_pdist", (DL_FUNC) &_RiemGrassmann_cpp_pdist, 2},
    {"_RiemGrassmann_cpp_pdist2", (DL_FUNC) &_RiemGrassmann_cpp_pdist2, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_RiemGrassmann(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
