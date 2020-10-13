#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

List chrom_fitness_score(IntegerMatrix case_genetic_data, IntegerMatrix complement_genetic_data, LogicalMatrix case_comp_differences,
                         IntegerVector target_snps, IntegerMatrix cases_minus_complements, LogicalMatrix both_one_mat,
                         LogicalMatrix block_ld_mat, IntegerVector weight_lookup,
                         int n_different_snps_weight = 2, int n_both_one_weight = 1,
                         int n_case_high_risk_thresh = 20, int outlier_sd = 2.5, bool epi_test = false) {

  // pick out the differences for the target snps
  LogicalMatrix case_comp_dif = transpose(subset_matrix_cols_l(case_comp_differences, target_snps));
  int n_target = target_snps.length();

  // also subset to target cols for the rest of the inputs
  case_genetic_data = transpose(subset_matrix_cols(case_genetic_data, target_snps));
  complement_genetic_data = transpose(subset_matrix_cols(complement_genetic_data, target_snps));
  cases_minus_complements = transpose(subset_matrix_cols(cases_minus_complements, target_snps));
  both_one_mat = transpose(subset_matrix_cols_l(both_one_mat, target_snps));

  // compute weighted difference vectors, determine risk related alleles
  List dif_vec_list = compute_dif_vecs(case_genetic_data, complement_genetic_data, case_comp_dif,
                                       cases_minus_complements, both_one_mat,
                                       weight_lookup, n_different_snps_weight, n_both_one_weight);

  // pick out the required pieces from function output
  arma::rowvec sum_dif_vecs = dif_vec_list["sum_dif_vecs"];
  IntegerVector family_weights = dif_vec_list["family_weights"];
  double invsum_family_weights = dif_vec_list["invsum_family_weights"];
  IntegerVector informative_families = dif_vec_list["informative_families"];
  int n_case_high_risk = dif_vec_list["n_case_high_risk"];
  int n_comp_high_risk = dif_vec_list["n_comp_high_risk"];
  IntegerMatrix case_high_inf = dif_vec_list["case_high_inf"];
  IntegerMatrix comp_high_inf = dif_vec_list["comp_high_inf"];
  int n_pos = dif_vec_list["n_pos"];
  int n_neg = dif_vec_list["n_neg"];
  IntegerVector pos_risk_int = dif_vec_list["pos_risk_int"];
  LogicalVector pos_risk = dif_vec_list["pos_risk"];
  IntegerVector neg_risk_int = dif_vec_list["neg_risk_int"];
  LogicalVector neg_risk = dif_vec_list["neg_risk"];
  IntegerMatrix case_inf = dif_vec_list["case_inf"];
  IntegerMatrix comp_inf = dif_vec_list["comp_inf"];
  IntegerMatrix both_one_inf = dif_vec_list["both_one_inf"];
  int n_informative_families = dif_vec_list["n_informative_families"];
  IntegerMatrix case_high_inf = dif_vec_list["case_high_inf"];
  IntegerMatrix comp_high_inf = dif_vec_list["comp_high_inf"];
  IntegerMatrix cases_minus_complements_inf = dif_vec_list["cases_minus_complements_inf"];

  // pick out misclassifications via outlier detection, indicating recessive mode of inheritance

  // initialize vector of pattern of inheritance
  CharacterVector risk_set_alleles(target_snps.length(), "1+");

  // only applies if we have at least 20 high risk case
  if (n_case_high_risk > n_case_high_risk_thresh) {

    NumericVector case_high_risk_means = rowMeans(case_high_inf);
    NumericVector case_high_risk_sd = rowSds(case_high_inf);
    double high_outlier_thresh = case_high_risk_means + outlier_sd * case_high_risk_sd;
    double low_outlier_thresh = case_high_risk_means - outlier_sd * case_high_risk_sd;
    IntegerMatrix all_high_risk = transpose(cbind(case_high_inf, comp_high_inf));
    int n_high_risk = all_high_risk.nrow();
    LogicalVector outliers(n_target, false);

    if (n_pos > 0){

      IntegerMatrix positive_high_risk = subset_matrix_cols(all_high_risk, pos_risk_int);
      LogicalMatrix positive_outlier_mat = check_pos_risk(positive_high_risk, low_outlier_thresh[pos_risk]);
      LogicalVector positive_outliers = colSums(positive_outlier_mat) > 0;
      outliers[pos_risk] = positive_outliers;

    }

    if (n_neg > 0){

      IntegerMatrix negative_high_risk = subset_matrix_cols(all_high_risk, neg_risk_int);
      LogicalMatrix negative_outlier_mat = check_neg_risk(negative_high_risk, high_outlier_thresh[neg_risk]);
      LogicalVector negative_outliers = colSums(negative_outlier_mat) > 0;
      outliers[neg_risk] = negative_outliers;

    }

    // if there are outliers, recompute the weights and associated statistics
    int n_outliers = sum(outliers);
    if (n_outliers > 0) {

      risk_set_alleles[outliers] = "2";
      LogicalVector pos_risk_outliers = outliers[pos_risk_int - 1];
      LogicalVector neg_risk_outliers = outliers[neg_risk_int - 1];
      IntegerVector pos_outlier_cols = pos_risk_int[pos_risk_outliers];
      IntegerVector neg_outlier_cols = neg_risk_int[neg_risk_outliers];

      // recode instances where the model appears to be recessive
      int n_pos_outliers = sum(pos_risk_outliers);
      int n_neg_outliers = sum(neg_risk_outliers);
      if (n_pos_outliers > 0){

        for (int k = 0; k < n_pos_outliers; k++){

          int outlier_idx = pos_outlier_cols[k] - 1;

          // first cases
          IntegerMatrix::Row case_target_row = case_inf(outlier_idx, _);
          IntegerVector case_new_vals = case_inf(outlier_idx, _);
          case_new_vals[case_new_vals == 1] = 0;
          case_target_row = case_new_vals;

          //then complements
          IntegerMatrix::Row comp_target_row = comp_inf(outlier_idx, _);
          IntegerVector comp_new_vals = comp_inf(outlier_idx, _);
          comp_new_vals[comp_new_vals == 1] = 0;
          comp_target_row = comp_new_vals;

          //then the both one matrix
          LogicalMatrix::Row both_one_target_row = both_one_inf(outlier_idx, _);
          LogicalVector both_one_new_vals(n_informative_families, false);
          both_one_target_row = both_one_new_vals;

        }

      }
      if (n_neg_outliers > 0){

        for (int k = 0; k < n_neg_outliers; k++){

          int outlier_idx = neg_outlier_cols[k] - 1;

          // first cases
          IntegerMatrix::Row case_target_row = case_inf(outlier_idx, _);
          IntegerVector case_new_vals = case_inf(outlier_idx, _);
          case_new_vals[case_new_vals == 1] = 2;
          case_target_row = case_new_vals;

          //then complements
          IntegerMatrix::Row comp_target_row = comp_inf(outlier_idx, _);
          IntegerVector comp_new_vals = comp_inf(outlier_idx, _);
          comp_new_vals[comp_new_vals == 1] = 2;
          comp_target_row = comp_new_vals;

          //then the both one matrix
          LogicalMatrix::Row both_one_target_row = both_one_inf(outlier_idx, _);
          LogicalVector both_one_new_vals(n_informative_families, false);
          both_one_target_row = both_one_new_vals;

        }

      }

      //recompute the number of informative families
      cases_minus_complements = sign_subtract_mat(case_inf, comp_inf);
      case_comp_dif = comp_mat_ne(cases_minus_complements, 0);
      dif_vec_list = compute_dif_vecs(case_inf, comp_inf, case_comp_dif, cases_minus_complements, both_one_inf,
                                      weight_lookup, n_different_snps_weight, n_both_one_weight, informative_families);

      //pick out required pieces
      sum_dif_vecs = dif_vec_list["sum_dif_vecs"];
      family_weights = dif_vec_list["family_weights"];
      invsum_family_weights = dif_vec_list["invsum_family_weights"];
      informative_families = dif_vec_list["informative_families"];
      n_case_high_risk = dif_vec_list["n_case_high_risk"];
      n_comp_high_risk = dif_vec_list["n_comp_high_risk"];
      case_high_inf = dif_vec_list["case_high_inf"];
      comp_high_inf = dif_vec_list["comp_high_inf"];
      n_pos = dif_vec_list["n_pos"];
      n_neg = dif_vec_list["n_neg"];
      pos_risk_int = dif_vec_list["pos_risk_int"];
      pos_risk = dif_vec_list["pos_risk"];
      neg_risk_int = dif_vec_list["neg_risk_int"];
      neg_risk = dif_vec_list["neg_risk"];
      case_inf = dif_vec_list["case_inf"];
      comp_inf = dif_vec_list["comp_inf"];
      both_one_inf = dif_vec_list["both_one_inf"];
      n_informative_families = dif_vec_list["n_informative_families"];
      case_high_inf = dif_vec_list["case_high_inf"];
      comp_high_inf = dif_vec_list["comp_high_inf"];
      cases_minus_complements_inf = dif_vec_list["cases_minus_complements_inf"];

    }
  }

  // count the number of risk alleles in those with the full risk set

  //cases
  IntegerMatrix case_high_inf_pos = subset_int_matrix_rows(case_high_inf, pos_risk_int);
  IntegerMatrix case_low_inf_pos = two_minus_mat(case_high_inf);
  case_low_inf_pos = subset_int_matrix_rows(case_low_inf_pos, neg_risk_int);
  IntegerVector case_high_risk_alleles = colSums(case_high_inf_pos) + colSums(case_low_inf_pos);
  NumericVector total_case_high_risk_alleles = sum(case_high_risk_alleles);

  //complements
  IntegerMatrix comp_high_inf_pos = subset_int_matrix_rows(comp_high_inf, pos_risk_int);
  IntegerMatrix comp_low_inf_pos = two_minus_mat(comp_high_inf);
  comp_low_inf_pos =  subset_int_matrix_rows(comp_low_inf_pos, neg_risk_int);
  IntegerVector comp_high_risk_alleles = colSums(comp_high_inf_pos) + colSums(comp_low_inf_pos);
  NumericVector total_comp_high_risk_alleles = sum(comp_high_risk_alleles);

  // compute scaling factor
  NumericVector rr = total_case_high_risk_alleles/(total_case_high_risk_alleles + total_comp_high_risk_alleles);
  if ((rr[0] <= 0) | all(is_na(rr)) | all(is_nan(rr))){
    rr = pow(10, -10);
  }

  // compute pseudo hotelling t2
  arma::rowvec mu_hat = sum_dif_vecs;
  arma::mat mu_hat_mat(n_informative_families, n_target);
  for (int i = 0; i < n_informative_families; i++){
    mu_hat_mat.row(i) = mu_hat;
  }
  arma::mat x = as<arma::mat>(transpose(cases_minus_complements_inf));
  arma::mat x_minus_mu_hat = x - mu_hat_mat;
  arma::mat weighted_x_minus_mu_hat = mult_rows_by_scalars_arma(x_minus_mu_hat, family_weights);
  arma::mat cov_mat = invsum_family_weights * trans(x) * y;

  // set cov elements to zero if SNPs are not in same ld block
  LogicalMatrix target_block_ld_mat = subset_matrix_l(block_ld_mat, target_snps, target_snps);
  for (int i = 0; i < n_target; i++){
    for (int j = 0; j < n_target; j++){
      if (!target_block_ld_mat(i, j)){

        cov_mat(i,j) = 0;

      }
    }
  }

  // get info for function output
  arma::vec elem_vars = sqrt(cov_mat.diag());
  sum_dif_vecs = sum_dif_vecs/trans(elem_vars);

  // if no variance, just make the element small
  arma::uvec zero_var = find(elem_vars == 0);
  double rep_val = pow(10, -10);
  sum_dif_vecs.elem(zero_var).fill(rep_val);

  // compute svd of cov_mat
  arma::mat U, V;
  arma::vec S;
  arma::svd(U, S, V, cov_mat, "dc");
  double lower_limit = std::numeric_limits<double>::epsilon();
  arma::uvec zero_sv = find(S < sqrt(lower_limit));
  double sv_rep_val = pow(10, 10);
  S.elem(zero_sv).fill(sv_rep_val);

  // compute fitness score
  double pseudo_t2 = (1/(1000*invsum_family_weights)) * sum(pow((mu_hat * U), 2) / trans(S));
  double fitness_score = pow(rr, 2) * pseudo_t2;

  // if desired, return the required information for the epistasis test
  if (epi_test){

    List res = List::create(Named("fitness_score") = fitness_score,
                            Named("sum_dif_vecs") = sum_dif_vecs,
                            Named("rr") = rr,
                            Named("pseudo_t2") = pseudo_t2,
                            Named("risk_set_alleles") = risk_set_alleles,
                            Named("high_risk_families") = informative_families);
    return(res);

  } else {

    List res = List::create(Named("fitness_score") = fitness_score,
                            Named("sum_dif_vecs") = sum_dif_vecs,
                            Named("rr") = rr,
                            Named("pseudo_t2") = pseudo_t2,
                            Named("risk_set_alleles") = risk_set_alleles);
    return(res);

  }

}
