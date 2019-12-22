// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// simu_prep
StringMatrix simu_prep(CharacterMatrix ref, int tot_len, int omit);
RcppExport SEXP _convSig_simu_prep(SEXP refSEXP, SEXP tot_lenSEXP, SEXP omitSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterMatrix >::type ref(refSEXP);
    Rcpp::traits::input_parameter< int >::type tot_len(tot_lenSEXP);
    Rcpp::traits::input_parameter< int >::type omit(omitSEXP);
    rcpp_result_gen = Rcpp::wrap(simu_prep(ref, tot_len, omit));
    return rcpp_result_gen;
END_RCPP
}
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
// feat2table5
std::vector<int> feat2table5(std::string b1, std::string b2, std::string b3, std::string b4, std::string b5);
RcppExport SEXP _convSig_feat2table5(SEXP b1SEXP, SEXP b2SEXP, SEXP b3SEXP, SEXP b4SEXP, SEXP b5SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type b1(b1SEXP);
    Rcpp::traits::input_parameter< std::string >::type b2(b2SEXP);
    Rcpp::traits::input_parameter< std::string >::type b3(b3SEXP);
    Rcpp::traits::input_parameter< std::string >::type b4(b4SEXP);
    Rcpp::traits::input_parameter< std::string >::type b5(b5SEXP);
    rcpp_result_gen = Rcpp::wrap(feat2table5(b1, b2, b3, b4, b5));
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
// shallow_loop5
S4 shallow_loop5(S4 mat, DataFrame fasta, DataFrame mut_file, CharacterVector uniq_samples);
RcppExport SEXP _convSig_shallow_loop5(SEXP matSEXP, SEXP fastaSEXP, SEXP mut_fileSEXP, SEXP uniq_samplesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< S4 >::type mat(matSEXP);
    Rcpp::traits::input_parameter< DataFrame >::type fasta(fastaSEXP);
    Rcpp::traits::input_parameter< DataFrame >::type mut_file(mut_fileSEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type uniq_samples(uniq_samplesSEXP);
    rcpp_result_gen = Rcpp::wrap(shallow_loop5(mat, fasta, mut_file, uniq_samples));
    return rcpp_result_gen;
END_RCPP
}
// conv_bimap
NumericMatrix conv_bimap(int N, int numbase);
RcppExport SEXP _convSig_conv_bimap(SEXP NSEXP, SEXP numbaseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type numbase(numbaseSEXP);
    rcpp_result_gen = Rcpp::wrap(conv_bimap(N, numbase));
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
// icgc_creater
DataFrame icgc_creater(DataFrame vcf_data, NumericMatrix sample_data, CharacterVector sample_names, int height, int old_height);
RcppExport SEXP _convSig_icgc_creater(SEXP vcf_dataSEXP, SEXP sample_dataSEXP, SEXP sample_namesSEXP, SEXP heightSEXP, SEXP old_heightSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame >::type vcf_data(vcf_dataSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type sample_data(sample_dataSEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type sample_names(sample_namesSEXP);
    Rcpp::traits::input_parameter< int >::type height(heightSEXP);
    Rcpp::traits::input_parameter< int >::type old_height(old_heightSEXP);
    rcpp_result_gen = Rcpp::wrap(icgc_creater(vcf_data, sample_data, sample_names, height, old_height));
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
// string_comptest
bool string_comptest(CharacterVector x);
RcppExport SEXP _convSig_string_comptest(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(string_comptest(x));
    return rcpp_result_gen;
END_RCPP
}
// string_comptest2
int string_comptest2(CharacterVector x, CharacterVector y, int a);
RcppExport SEXP _convSig_string_comptest2(SEXP xSEXP, SEXP ySEXP, SEXP aSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< int >::type a(aSEXP);
    rcpp_result_gen = Rcpp::wrap(string_comptest2(x, y, a));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_convSig_simu_prep", (DL_FUNC) &_convSig_simu_prep, 3},
    {"_convSig_feat2table3", (DL_FUNC) &_convSig_feat2table3, 3},
    {"_convSig_feat2table5", (DL_FUNC) &_convSig_feat2table5, 5},
    {"_convSig_reverse_transform", (DL_FUNC) &_convSig_reverse_transform, 1},
    {"_convSig_shallow_loop3", (DL_FUNC) &_convSig_shallow_loop3, 4},
    {"_convSig_shallow_loop5", (DL_FUNC) &_convSig_shallow_loop5, 4},
    {"_convSig_conv_bimap", (DL_FUNC) &_convSig_conv_bimap, 2},
    {"_convSig_RM_nonSNP", (DL_FUNC) &_convSig_RM_nonSNP, 2},
    {"_convSig_icgc_creater", (DL_FUNC) &_convSig_icgc_creater, 5},
    {"_convSig_timesTwo", (DL_FUNC) &_convSig_timesTwo, 1},
    {"_convSig_timesTwoList", (DL_FUNC) &_convSig_timesTwoList, 1},
    {"_convSig_WeirdVector", (DL_FUNC) &_convSig_WeirdVector, 1},
    {"_convSig_string_comptest", (DL_FUNC) &_convSig_string_comptest, 1},
    {"_convSig_string_comptest2", (DL_FUNC) &_convSig_string_comptest2, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_convSig(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
