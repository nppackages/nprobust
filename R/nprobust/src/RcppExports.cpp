// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// lpbwce
List lpbwce(arma::vec y, arma::vec x, arma::vec K, arma::vec L, arma::vec res, double c, int p, int q, double h, double b, int deriv, int fact);
RcppExport SEXP _nprobust_lpbwce(SEXP ySEXP, SEXP xSEXP, SEXP KSEXP, SEXP LSEXP, SEXP resSEXP, SEXP cSEXP, SEXP pSEXP, SEXP qSEXP, SEXP hSEXP, SEXP bSEXP, SEXP derivSEXP, SEXP factSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type K(KSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type L(LSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type res(resSEXP);
    Rcpp::traits::input_parameter< double >::type c(cSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type q(qSEXP);
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    Rcpp::traits::input_parameter< int >::type deriv(derivSEXP);
    Rcpp::traits::input_parameter< int >::type fact(factSEXP);
    rcpp_result_gen = Rcpp::wrap(lpbwce(y, x, K, L, res, c, p, q, h, b, deriv, fact));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_nprobust_lpbwce", (DL_FUNC) &_nprobust_lpbwce, 12},
    {NULL, NULL, 0}
};

RcppExport void R_init_nprobust(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
