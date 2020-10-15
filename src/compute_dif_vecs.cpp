#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
List compute_dif_vecs(IntegerMatrix case_genetic_data, IntegerMatrix complement_genetic_data, LogicalMatrix case_comp_dif,
                      IntegerMatrix cases_minus_complements, LogicalMatrix both_one_mat,
                      IntegerVector weight_lookup, int n_different_snps_weight = 2, int n_both_one_weight = 1,
                      IntegerVector prev_informative_families = IntegerVector::create(NA_INTEGER)){

  // determine whether families are informative for the set of target_snps
  IntegerVector total_different_snps = colSums(case_comp_dif);
  LogicalVector informative_families_l = total_different_snps != 0;
  IntegerVector family_idx = seq_len(total_different_snps.length());
  IntegerVector informative_families = family_idx[informative_families_l];
  int n_informative_families = sum(informative_families_l);
  IntegerMatrix case_inf = subset_matrix_cols(case_genetic_data, informative_families);
  IntegerMatrix comp_inf = subset_matrix_cols(comp_genetic_data, informative_families);
  IntegerMatrix cases_minus_complements_inf = subset_matrix_cols(cases_minus_complements, informative_families);

  // compute weights
  LogicalMatrix both_one_inf = subset_matrix_cols_l(both_one_mat, informative_families);
  IntegerVector both_one = colSums(both_one_inf);
  IntegerVector total_different_inf = total_different_snps[informative_families_l];
  IntegerVector weighted_informativeness = n_both_one_weight * both_one + n_different_snps_weight * total_different_inf;
  IntegerVector family_weights = weight_lookup[weighted_informativeness - 1];
  double invsum_family_weights = 1 / sum(family_weights);

  // compute weighted difference vectors for cases vs complements
  NumericVector std_family_weights = family_weights * invsum_family_weights;
  arma::mat dif_vecs = mult_rows_by_scalars(transpose(cases_minus_complements_inf), std_family_weights);

  // take the sum of the case - complement difference vectors over families
  arma::rowvec sum_dif_vecs arma::sum(dif_vecs, 0);

  // determine how many cases and complements actually have the proposed risk set
  arma::rowvec risk_dirs = sign(sum_dif_vecs);
  LogicalVector pos_risk = pos_risk_dir(risk_dirs);
  LogicalVector neg_risk = neg_risk_dir(risk_dirs);
  int n_target = case_comp_dif.nrow();
  IntegerVector idx_vec = seq_len(n_target);

  // Use different cases to cull out useless operations
  int n_pos = sum(pos_risk);
  int n_neg = sum(neg_risk);
  if (n_neg > 0) {

    IntegerVector neg_risk_int = idx_vec[neg_risk];
    IntegerMatrix case_inf_neg subset_int_matrix_rows(case_inf, neg_risk_int);
    IntegerMatrix comp_inf_neg subset_int_matrix_rows(comp_inf, neg_risk_int);
    LogicalMatrix n1 = comp_mat_lt(case_inf_neg, 2);
    LogicalMatrix n2 = comp_mat_lt(comp_inf_neg, 2);

  } else {

    IntegerVector neg_risk_int = IntegerVector::create(NA_INTEGER);

  }
  if (n_pos > 0) {

    IntegerVector pos_risk_int = idx_vec[pos_risk];
    IntegerMatrix case_inf_pos subset_int_matrix_rows(case_inf, pos_risk_int);
    IntegerMatrix comp_inf_pos subset_int_matrix_rows(comp_inf, pos_risk_int);
    LogicalMatrix p1 = comp_mat_gt(case_inf_pos, 0);
    LogicalMatrix p2 = comp_mat_gt(comp_inf_pos, 0);

  } else {

    IntegerVector pos_risk_int = IntegerVector::create(NA_INTEGER);

  }
  if (n_pos == 0){

    IntegerVector n1_colsums = colSums(n1);
    IntegerVector n2_colsums = colSums(n2);
    LogicalVector case_high_risk = n1_colsums == n_target;
    LogicalVector comp_high_risk = n2_colsums == n_target;

  } else if (n_neg == 0){

    IntegerVector p1_colsums = colSums(p1);
    IntegerVector p2_colsums = colSums(p2);
    LogicalVector case_high_risk = p1_colsums == n_target;
    LogicalVector comp_high_risk = p2_colsums == n_target;

  } else {

    IntegerVector n1_colsums = colSums(n1);
    IntegerVector n2_colsums = colSums(n2);
    IntegerVector p1_colsums = colSums(p1);
    IntegerVector p2_colsums = colSums(p2);
    LogicalVector case_high_risk = (p1_colsums + n1_colsums) == n_target;
    LogicalVector comp_high_risk = (p2_colsums + n2_colsums) == n_target;

  }

  int n_case_high_risk = sum(case_high_risk);
  in n_comp_high_risk = sum(comp_high_risk);
  IntegerVector case_high_risk_int = family_idx[case_high_risk];
  IntegerVector comp_high_risk_int = family_idx[comp_high_risk];
  IntegerMatrix case_high_inf = subset_matrix_cols(case_inf, case_high_risk_int);
  IntegerMatrix comp_high_inf = subset_matrix_cols(comp_inf, comp_high_risk_int);

  // for epi test
  bool prev_inf = all(is_na(prev_informative_families));
  if (!prev_inf){

    informative_families = prev_informative_families[informative_families - 1];

  }

  //return list
  List res = List::create(Named("sum_dif_vecs") = sum_dif_vecs,
                          Named("family_weights") = family_weights,
                          Named("invsum_family_weights") = invsum_family_weights,
                          Named("n_case_high_risk") = n_case_high_risk,
                          Named("n_comp_high_risk") = n_comp_high_risk,
                          Named("case_high_inf") = case_high_inf,
                          Named("comp_high_inf") = comp_high_inf,
                          Named("n_pos") = n_pos,
                          Named("n_neg") = n_neg,
                          Named("pos_risk_int") = pos_risk_int,
                          Named("pos_risk") = pos_risk,
                          Named("neg_risk_int") = neg_risk_int,
                          Named("neg_risk") = neg_risk,
                          Named("case_inf") = case_inf,
                          Named("comp_inf") = comp_inf,
                          Named("cases_minus_complements_inf") = cases_minus_complements_inf,
                          Named("both_one_inf") = both_one_inf,
                          Named("n_informative_families") = n_informative_families,
                          Named("case_high_inf") = case_high_inf,
                          Named("comp_high_inf") = comp_high_inf,
                          Named("informative_families") = informative_families);
  return(res)

}
