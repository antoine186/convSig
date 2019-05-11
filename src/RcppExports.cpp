// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// feat2table3
std::vector<int> feat2table3(std::string b1, std::string b2, std::string b3);
RcppExport SEXP _convSig_feat2table3(SEXP b1SEXP, SEXP b2SEXP, SEXP b3SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type b1(b1SEXP);
    Rcpp::traits::input_parameter< std::string >::type b2(b2SEXP);
    Rcpp::traits::input_parameter< std::string >::type b3(b3SEXP);
    rcpp_result_gen = Rcpp::wrap(feat2table3(b1, b2, b3));
    return rcpp_result_gen;
END_RCPP
}
// reverse_transform
DataFrame reverse_transform(DataFrame alleles);
RcppExport SEXP _convSig_reverse_transform(SEXP allelesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame >::type alleles(allelesSEXP);
    rcpp_result_gen = Rcpp::wrap(reverse_transform(alleles));
    return rcpp_result_gen;
END_RCPP
}
// shallow_loop3
S4 shallow_loop3(S4 mat, DataFrame fasta, DataFrame mut_file, CharacterVector uniq_samples);
RcppExport SEXP _convSig_shallow_loop3(SEXP matSEXP, SEXP fastaSEXP, SEXP mut_fileSEXP, SEXP uniq_samplesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< S4 >::type mat(matSEXP);
    Rcpp::traits::input_parameter< DataFrame >::type fasta(fastaSEXP);
    Rcpp::traits::input_parameter< DataFrame >::type mut_file(mut_fileSEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type uniq_samples(uniq_samplesSEXP);
    rcpp_result_gen = Rcpp::wrap(shallow_loop3(mat, fasta, mut_file, uniq_samples));
    return rcpp_result_gen;
END_RCPP
}
// RM_nonSNP
LogicalVector RM_nonSNP(DataFrame startend, SEXP ar);
RcppExport SEXP _convSig_RM_nonSNP(SEXP startendSEXP, SEXP arSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame >::type startend(startendSEXP);
    Rcpp::traits::input_parameter< SEXP >::type ar(arSEXP);
    rcpp_result_gen = Rcpp::wrap(RM_nonSNP(startend, ar));
    return rcpp_result_gen;
END_RCPP
}
// timesTwo
RcppExport SEXP timesTwo(SEXP x);
RcppExport SEXP _convSig_timesTwo(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(timesTwo(x));
    return rcpp_result_gen;
END_RCPP
}
// timesTwoList
RcppExport SEXP timesTwoList(SEXP x);
RcppExport SEXP _convSig_timesTwoList(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(timesTwoList(x));
    return rcpp_result_gen;
END_RCPP
}
// WeirdVector
NumericVector WeirdVector(NumericVector x);
RcppExport SEXP _convSig_WeirdVector(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(WeirdVector(x));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_convSig_feat2table3", (DL_FUNC) &_convSig_feat2table3, 3},
    {"_convSig_reverse_transform", (DL_FUNC) &_convSig_reverse_transform, 1},
    {"_convSig_shallow_loop3", (DL_FUNC) &_convSig_shallow_loop3, 4},
    {"_convSig_RM_nonSNP", (DL_FUNC) &_convSig_RM_nonSNP, 2},
    {"_convSig_timesTwo", (DL_FUNC) &_convSig_timesTwo, 1},
    {"_convSig_timesTwoList", (DL_FUNC) &_convSig_timesTwoList, 1},
    {"_convSig_WeirdVector", (DL_FUNC) &_convSig_WeirdVector, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_convSig(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
