// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// runLinearisedSISModel
List runLinearisedSISModel(double lambda, double gamma, double max_time);
RcppExport SEXP _tmi_runLinearisedSISModel(SEXP lambdaSEXP, SEXP gammaSEXP, SEXP max_timeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< double >::type max_time(max_timeSEXP);
    rcpp_result_gen = Rcpp::wrap(runLinearisedSISModel(lambda, gamma, max_time));
    return rcpp_result_gen;
END_RCPP
}
// runIndividualSISModel
List runIndividualSISModel(double beta, double gamma, int num_individuals, int num_infected, double max_time);
RcppExport SEXP _tmi_runIndividualSISModel(SEXP betaSEXP, SEXP gammaSEXP, SEXP num_individualsSEXP, SEXP num_infectedSEXP, SEXP max_timeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< double >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< int >::type num_individuals(num_individualsSEXP);
    Rcpp::traits::input_parameter< int >::type num_infected(num_infectedSEXP);
    Rcpp::traits::input_parameter< double >::type max_time(max_timeSEXP);
    rcpp_result_gen = Rcpp::wrap(runIndividualSISModel(beta, gamma, num_individuals, num_infected, max_time));
    return rcpp_result_gen;
END_RCPP
}
// runIndividualSIIModel
List runIndividualSIIModel(double beta, double gamma, int num_individuals, int num_infected, double max_time);
RcppExport SEXP _tmi_runIndividualSIIModel(SEXP betaSEXP, SEXP gammaSEXP, SEXP num_individualsSEXP, SEXP num_infectedSEXP, SEXP max_timeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< double >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< int >::type num_individuals(num_individualsSEXP);
    Rcpp::traits::input_parameter< int >::type num_infected(num_infectedSEXP);
    Rcpp::traits::input_parameter< double >::type max_time(max_timeSEXP);
    rcpp_result_gen = Rcpp::wrap(runIndividualSIIModel(beta, gamma, num_individuals, num_infected, max_time));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_tmi_runLinearisedSISModel", (DL_FUNC) &_tmi_runLinearisedSISModel, 3},
    {"_tmi_runIndividualSISModel", (DL_FUNC) &_tmi_runIndividualSISModel, 5},
    {"_tmi_runIndividualSIIModel", (DL_FUNC) &_tmi_runIndividualSIIModel, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_tmi(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
