#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
int sign_scalar(int x) {

  if (x > 0) {

    return 1;

  } else if (x == 0) {

    return 0;

  } else {

    return -1;

  }

}

// [[Rcpp::export]]
IntegerMatrix subset_matrix(IntegerMatrix in_matrix, IntegerVector rows, IntegerVector cols){

  int n_rows = rows.length();
  int n_cols = cols.length();
  IntegerMatrix out_matrix(n_rows, n_cols);

  for (int i = 0; i < n_cols; i++){

    IntegerMatrix::Column original_col = in_matrix(_, cols[i]-1);
    IntegerMatrix::Column new_col = out_matrix(_, i);

    for (int j = 0; j < n_rows; j++){

      new_col[j] = original_col[rows[j] - 1];

    }
  }
  return(out_matrix);
}

// [[Rcpp::export]]
NumericVector weighted_sub_colsums(IntegerMatrix in_mat, IntegerVector target_rows,
                                   IntegerVector target_cols, NumericVector row_weights){

  int n_cols = target_cols.length();
  int n_rows = target_rows.length();
  NumericVector out_vec(n_cols, 0.0);
  for (int i = 0; i < n_rows; i++){

    int this_row = target_rows[i] - 1;
    double row_weight = row_weights[i];

    for (int j = 0; j < n_cols; j ++){

      int this_col = target_cols[j] - 1;
      out_vec[j] += (row_weight*in_mat(this_row, this_col));

    }
  }
  return(out_vec);
}


// [[Rcpp::export]]
IntegerVector sub_rowsums_start(IntegerMatrix in_mat, IntegerVector target_cols){

  int n_cols = target_cols.length();
  int n_rows = in_mat.nrow();
  IntegerVector out_vec(n_rows, 0);
  for (int i = 0; i < n_rows; i++){

    for (int j = 0; j < n_cols; j ++){

      int this_col = target_cols[j] - 1;
      out_vec[i] += in_mat(i, this_col);

    }
  }
  return(out_vec);
}

// [[Rcpp::export]]
IntegerVector sub_colsums(IntegerMatrix in_mat, IntegerVector target_rows, IntegerVector target_cols){

  int n_cols = target_cols.length();
  int n_rows = target_rows.length();
  IntegerVector out_vec(n_cols, 0);
  for (int i = 0; i < n_rows; i++){

    int this_row = target_rows[i] - 1;

    for (int j = 0; j < n_cols; j ++){

      int this_col = target_cols[j] - 1;
      out_vec[j] += in_mat(this_row, this_col);

    }
  }
  return(out_vec);
}



// [[Rcpp::export]]
IntegerVector sub_colsums_2minus(IntegerMatrix in_mat, IntegerVector target_rows, IntegerVector target_cols){

  int n_cols = target_cols.length();
  int n_rows = target_rows.length();
  IntegerVector out_vec(n_cols, 0);
  for (int i = 0; i < n_rows; i++){

    int this_row = target_rows[i] - 1;

    for (int j = 0; j < n_cols; j ++){

      int this_col = target_cols[j] - 1;
      out_vec[j] += 2 - in_mat(this_row, this_col);

    }
  }
  return(out_vec);
}



// [[Rcpp::export]]
IntegerVector sub_rowsums(IntegerMatrix in_mat, IntegerVector target_rows, IntegerVector target_cols){

  int n_cols = target_cols.length();
  int n_rows = target_rows.length();
  IntegerVector out_vec(n_rows, 0);
  for (int i = 0; i < n_rows; i++){

    int this_row = target_rows[i] - 1;

    for (int j = 0; j < n_cols; j ++){

      int this_col = target_cols[j] - 1;
      out_vec[i] += in_mat(this_row, this_col);

    }
  }
  return(out_vec);
}


// [[Rcpp::export]]
IntegerVector n_pos_high_risk(IntegerMatrix in_mat, IntegerVector target_rows, IntegerVector target_cols){

  int n_cols = target_cols.length();
  int n_rows = target_rows.length();
  IntegerVector out_vec(n_rows, 0);
  for (int i = 0; i < n_rows; i++){

    int this_row = target_rows[i] - 1;

    for (int j = 0; j < n_cols; j ++){

      int this_col = target_cols[j] - 1;
      int obs_val = in_mat(this_row, this_col);
      if (obs_val > 0){

        out_vec[i] += 1;

      }
    }
  }
  return(out_vec);
}

// [[Rcpp::export]]
IntegerVector n_neg_high_risk(IntegerMatrix in_mat, IntegerVector target_rows, IntegerVector target_cols){

  int n_cols = target_cols.length();
  int n_rows = target_rows.length();
  IntegerVector out_vec(n_rows, 0);
  for (int i = 0; i < n_rows; i++){

    int this_row = target_rows[i] - 1;

    for (int j = 0; j < n_cols; j ++){

      int this_col = target_cols[j] - 1;
      int obs_val = in_mat(this_row, this_col);
      if (obs_val < 2){

        out_vec[i] += 1;

      }
    }
  }
  return(out_vec);
}


// [[Rcpp::export]]
NumericVector sub_colmeans(IntegerMatrix in_mat, IntegerVector target_rows, IntegerVector target_cols){

  IntegerVector target_colsums = sub_colsums(in_mat, target_rows, target_cols);
  NumericVector out_vec = as<NumericVector>(target_colsums);
  double n = target_rows.length();
  out_vec = out_vec / n;
  return(out_vec);

}

// [[Rcpp::export]]
NumericVector sub_colsds(IntegerMatrix in_mat, IntegerVector target_rows, IntegerVector target_cols){

  NumericVector col_means = sub_colmeans(in_mat, target_rows, target_cols);
  int n_cols = target_cols.length();
  int n_rows = target_rows.length();
  NumericVector out_vec(n_cols, 0.0);
  for (int i = 0; i < n_rows; i++){

    int this_row = target_rows[i] - 1;

    for (int j = 0; j < n_cols; j ++){

      int this_col = target_cols[j] - 1;
      double col_mean = col_means[j];
      double obs_val = in_mat(this_row, this_col) - col_mean;
      double fill_val = pow(obs_val, 2);
      out_vec[j] += fill_val;

    }
  }
  double n_minus_1 = n_rows - 1;
  out_vec = sqrt(out_vec / n_minus_1);
  return(out_vec);
}




// [[Rcpp::export]]
List find_high_risk(int n_target, int n_pos, int n_neg, IntegerVector neg_risk_int, IntegerVector pos_risk_int,
                    IntegerMatrix case_data, IntegerMatrix comp, IntegerVector informative_families, IntegerVector target_snps){

  if (n_neg > 0 & n_pos > 0) {

    IntegerVector neg_snp_cols = target_snps[neg_risk_int - 1];
    IntegerVector n1 = n_neg_high_risk(case_data, informative_families, neg_snp_cols);
    IntegerVector n2 = n_neg_high_risk(comp, informative_families, neg_snp_cols);

    IntegerVector pos_snp_cols = target_snps[pos_risk_int - 1];
    IntegerVector p1 = n_pos_high_risk(case_data, informative_families, pos_snp_cols);
    IntegerVector p2 = n_pos_high_risk(comp, informative_families, pos_snp_cols);

    LogicalVector case_high_risk = (p1 + n1) == n_target;
    LogicalVector comp_high_risk = (p2 + n2) == n_target;

    List res = List::create(Named("case_high_risk") = case_high_risk,
                            Named("comp_high_risk") = comp_high_risk);
    return(res);

  } else if (n_neg > 0 & n_pos == 0) {

    IntegerVector neg_snp_cols = target_snps[neg_risk_int - 1];
    IntegerVector n1 = n_neg_high_risk(case_data, informative_families, neg_snp_cols);
    IntegerVector n2 = n_neg_high_risk(comp, informative_families, neg_snp_cols);

    LogicalVector case_high_risk = n1 == n_target;
    LogicalVector comp_high_risk = n2 == n_target;

    List res = List::create(Named("case_high_risk") = case_high_risk,
                            Named("comp_high_risk") = comp_high_risk);
    return(res);

  } else {

    IntegerVector pos_snp_cols = target_snps[pos_risk_int - 1];
    IntegerVector p1 = n_pos_high_risk(case_data, informative_families, pos_snp_cols);
    IntegerVector p2 = n_pos_high_risk(comp, informative_families, pos_snp_cols);

    LogicalVector case_high_risk = p1 == n_target;
    LogicalVector comp_high_risk = p2 == n_target;

    List res = List::create(Named("case_high_risk") = case_high_risk,
                            Named("comp_high_risk") = comp_high_risk);
    return(res);

  }

}

// [[Rcpp::export]]
List compute_dif_vecs(IntegerMatrix case_genetic_data, IntegerMatrix comp_genetic_data, IntegerMatrix case_comp_dif,
                      IntegerVector target_snps, IntegerMatrix cases_minus_complements, IntegerMatrix both_one_mat,
                      NumericVector weight_lookup, int n_different_snps_weight = 2, int n_both_one_weight = 1,
                      IntegerVector prev_informative_families = IntegerVector::create(NA_INTEGER)){

  // determine whether families are informative for the set of target_snps
  IntegerVector total_different_snps = sub_rowsums_start(case_comp_dif, target_snps);
  LogicalVector informative_families_l = total_different_snps != 0;
  IntegerVector family_idx = seq_len(total_different_snps.length());
  IntegerVector informative_families = family_idx[informative_families_l];

  // check prev_informative_families and update inf families if needed
  bool prev_inf = all(is_na(prev_informative_families));
  if (!prev_inf){

    IntegerVector changed_families = setdiff(informative_families, prev_informative_families);
    if (changed_families.length() > 0){

      informative_families_l[changed_families - 1] = false;
      total_different_snps[changed_families - 1] = 0;
      informative_families = family_idx[informative_families_l];

    }

  }

  LogicalVector uninformative_families_l = total_different_snps == 0;
  total_different_snps = total_different_snps[total_different_snps > 0];

  // compute weights
  IntegerVector both_one = sub_rowsums(both_one_mat, informative_families, target_snps);
  IntegerVector weighted_informativeness = n_both_one_weight * both_one + n_different_snps_weight * total_different_snps;
  NumericVector family_weights(weighted_informativeness.length());
  for (int i = 0; i < weighted_informativeness.length(); i++){

    int weight_i = weighted_informativeness[i];
    family_weights[i] = weight_lookup[weight_i - 1];

  }
  double invsum_family_weights = 1 / sum(family_weights);

  // compute weighted difference vectors for cases vs complements
  NumericVector std_family_weights = family_weights * invsum_family_weights;
  //IntegerVector sum_dif_vecs_start = weighted_sub_colsums(cases_minus_complements, informative_families,
  //                                                        target_snps, family_weights);
  NumericVector sum_dif_vecs = weighted_sub_colsums(cases_minus_complements, informative_families,
                                                    target_snps, std_family_weights);

  // take the sum of the case - complement difference vectors over families
  //NumericVector sum_dif_vecs = as<NumericVector>(sum_dif_vecs_start);

  //for (int i = 0; i < sum_dif_vecs.length(); i++){

  //  sum_dif_vecs[i] *= invsum_family_weights;

  //}

  // determine how many cases and complements actually have the proposed risk set
  IntegerVector risk_dirs = sign(sum_dif_vecs);

  LogicalVector pos_risk = risk_dirs > 0;
  int n_pos = sum(pos_risk);

  LogicalVector neg_risk = risk_dirs <= 0;
  int n_neg = sum(neg_risk);

  int n_target = target_snps.length();
  IntegerVector idx_vec = seq_len(n_target);

  IntegerVector pos_risk_int = idx_vec[pos_risk];
  IntegerVector neg_risk_int = idx_vec[neg_risk];

  //using helper function to pick out high risk families
  List high_risk = find_high_risk(n_target, n_pos, n_neg, neg_risk_int, pos_risk_int, case_genetic_data,
                                  comp_genetic_data, informative_families, target_snps);
  LogicalVector case_high_risk = high_risk["case_high_risk"];
  LogicalVector comp_high_risk = high_risk["comp_high_risk"];

  // return list
  List res = List::create(Named("sum_dif_vecs") = sum_dif_vecs,
                          Named("family_weights") = family_weights,
                          Named("invsum_family_weights") = invsum_family_weights,
                          Named("n_pos") = n_pos,
                          Named("n_neg") = n_neg,
                          Named("pos_risk") = pos_risk,
                          Named("neg_risk") = neg_risk,
                          Named("case_high_risk") = case_high_risk,
                          Named("comp_high_risk") = comp_high_risk,
                          Named("informative_families") = informative_families);
  return(res);
}


// [[Rcpp::export]]
List chrom_fitness_score(IntegerMatrix case_genetic_data_in, IntegerMatrix complement_genetic_data_in, IntegerMatrix case_comp_differences_in,
                         IntegerVector target_snps, IntegerMatrix cases_minus_complements_in, IntegerMatrix both_one_mat_in,
                         LogicalMatrix block_ld_mat, NumericVector weight_lookup,
                         int n_different_snps_weight = 2, int n_both_one_weight = 1,
                         int n_case_high_risk_thresh = 20, double outlier_sd = 2.5, bool epi_test = false) {

  // need to deep copy inputs to avoid overwriting
  IntegerMatrix case_genetic_data = clone(case_genetic_data_in);
  IntegerMatrix complement_genetic_data = clone(complement_genetic_data_in);
  IntegerMatrix case_comp_differences = clone(case_comp_differences_in);
  IntegerMatrix cases_minus_complements = clone(cases_minus_complements_in);
  IntegerMatrix both_one_mat = both_one_mat_in;

  // n target snps
  int n_target = target_snps.length();

  // compute weighted difference vectors, determine risk related alleles
  List dif_vec_list = compute_dif_vecs(case_genetic_data, complement_genetic_data, case_comp_differences,
                                       target_snps, cases_minus_complements, both_one_mat,
                                       weight_lookup, n_different_snps_weight, n_both_one_weight);

  // pick out the required pieces from function output
  NumericVector sum_dif_vecs = dif_vec_list["sum_dif_vecs"];
  NumericVector family_weights = dif_vec_list["family_weights"];
  double invsum_family_weights = dif_vec_list["invsum_family_weights"];
  IntegerVector informative_families = dif_vec_list["informative_families"];
  int n_pos = dif_vec_list["n_pos"];
  int n_neg = dif_vec_list["n_neg"];
  LogicalVector pos_risk = dif_vec_list["pos_risk"];
  LogicalVector neg_risk = dif_vec_list["neg_risk"];
  LogicalVector case_high_risk = dif_vec_list["case_high_risk"];
  LogicalVector comp_high_risk = dif_vec_list["comp_high_risk"];
  int n_informative_families = informative_families.length();
  int n_case_high_risk = sum(case_high_risk);

  // pick out misclassifications via outlier detection, indicating recessive mode of inheritance

  // initialize vector of pattern of inheritance
  CharacterVector risk_set_alleles(target_snps.length(), "1+");

  // count outliers
  int outlier_count = 0;

  // only applies if we have at least 20 high risk case
  if (n_case_high_risk > n_case_high_risk_thresh) {

    IntegerVector case_high_inf_rows = informative_families[case_high_risk];
    IntegerVector comp_high_inf_rows = informative_families[comp_high_risk];

    NumericVector case_high_risk_means = sub_colmeans(case_genetic_data, case_high_inf_rows, target_snps);
    NumericVector case_high_risk_sd = sub_colsds(case_genetic_data, case_high_inf_rows, target_snps);
    NumericVector high_outlier_thresh = case_high_risk_means + outlier_sd * case_high_risk_sd;
    NumericVector low_outlier_thresh = case_high_risk_means - outlier_sd * case_high_risk_sd;

    LogicalVector outliers(n_target, false);

    if (n_pos > 0){

      IntegerVector p_tmp = seq_len(n_target);
      IntegerVector pos_risk_idx = p_tmp[pos_risk];
      IntegerVector pos_cols = target_snps[pos_risk];

      for (int i = 0; i < n_pos; i++){

        int this_col = pos_cols[i] - 1;
        int this_thresh = pos_risk_idx[i] - 1;
        double thresh_val = low_outlier_thresh[this_thresh];

        if (thresh_val > 1){

          // check cases for outliers
          int pos_outliers = 0;
          for (int m = 0; m < case_high_inf_rows.length(); m++){

            int this_row = case_high_inf_rows[m] - 1;
            if (case_genetic_data(this_row, this_col) < thresh_val){

              pos_outliers += 1;
              outlier_count += 1;
              break;

            }

          }

          // if we have found an outlier, move to recoding
          // otherwise check complements for outliers
          if (pos_outliers == 0){

            for (int m = 0; m < comp_high_inf_rows.length(); m++){

              int this_row = comp_high_inf_rows[m] - 1;
              if (complement_genetic_data(this_row, this_col) < thresh_val){

                pos_outliers += 1;
                outlier_count += 1;
                break;

              }

            }

          }

          // if we found any outliers, recode
          if (pos_outliers != 0){

            // indicate we need two risk alleles
            int this_index = pos_risk_idx[i] - 1;
            risk_set_alleles[this_index] = "2";

            // loop over informative families
            for (int n = 0; n < informative_families.length(); n++){

              int family_row = informative_families[n] - 1;

              // recode cases
              if (case_genetic_data(family_row, this_col) == 1){

                case_genetic_data(family_row, this_col) = 0;

              }
              //recode complements
              if (complement_genetic_data(family_row, this_col) == 1){

                complement_genetic_data(family_row, this_col) = 0;

              }
              // recode both one
              if (both_one_mat(family_row, this_col) == 1){

                both_one_mat(family_row, this_col) = 0;

              }

              //recode cases minus complements
              if ((case_genetic_data(family_row, this_col) == 1) | (complement_genetic_data(family_row, this_col) == 1)){

                cases_minus_complements(family_row, this_col) = sign_scalar(case_genetic_data(family_row, this_col) - complement_genetic_data(family_row, this_col));
                case_comp_differences(family_row, this_col) = abs(cases_minus_complements(family_row, this_col));

              }

            }

          }

        }

      }

    }

    if (n_neg > 0){

      IntegerVector n_tmp = seq_len(n_target);
      IntegerVector neg_risk_idx = n_tmp[neg_risk];
      IntegerVector neg_cols = target_snps[neg_risk];

      for (int i = 0; i < n_neg; i++){

        int this_col = neg_cols[i] - 1;
        int this_thresh = neg_risk_idx[i] - 1;
        double thresh_val = high_outlier_thresh[this_thresh];

        if (thresh_val < 1){

          // check cases for outliers
          int neg_outliers = 0;
          for (int m = 0; m < case_high_inf_rows.length(); m++){

            int this_row = case_high_inf_rows[m] - 1;
            if (case_genetic_data(this_row, this_col) > thresh_val){

              neg_outliers += 1;
              outlier_count += 1;
              break;

            }

          }

          // if we have found an outlier, move to recoding
          // otherwise check complements for outliers
          if (neg_outliers == 0){

            for (int m = 0; m < comp_high_inf_rows.length(); m++){

              int this_row = comp_high_inf_rows[m] - 1;
              if (complement_genetic_data(this_row, this_col) > thresh_val){

                neg_outliers += 1;
                outlier_count += 1;
                break;

              }

            }

          }

          // if we found any outliers, recode
          if (neg_outliers != 0){

            // indicate we need two risk alleles
            int this_index = neg_risk_idx[i] - 1;
            risk_set_alleles[this_index] = "2";

            // loop over informative families
            for (int n = 0; n < informative_families.length(); n++){

              int family_row = informative_families[n] - 1;

              // recode cases
              if (case_genetic_data(family_row, this_col) == 1){

                case_genetic_data(family_row, this_col) = 2;

              }
              //recode complements
              if (complement_genetic_data(family_row, this_col) == 1){

                complement_genetic_data(family_row, this_col) = 2;

              }
              // recode both one
              if (both_one_mat(family_row, this_col) == 1){

                both_one_mat(family_row, this_col) = 0;

              }

              //recode cases minus complements
              if ((case_genetic_data(family_row, this_col) == 1) | (complement_genetic_data(family_row, this_col) == 1)){

                cases_minus_complements(family_row, this_col) = sign_scalar(case_genetic_data(family_row, this_col) - complement_genetic_data(family_row, this_col));
                case_comp_differences(family_row, this_col) = abs(cases_minus_complements(family_row, this_col));

              }

            }

          }

        }

      }

    }

    // if there are outliers, recompute the weights and associated statistics
    if (outlier_count > 0) {

      //recompute the number of informative families
      List dif_vec_list = compute_dif_vecs(case_genetic_data, complement_genetic_data, case_comp_differences, target_snps,
                                      cases_minus_complements, both_one_mat,
                                      weight_lookup, n_different_snps_weight, n_both_one_weight, informative_families);

      //pick out required pieces
      // the tmp var versions are used because direct reassignment
      // causes compile errors
      NumericVector sum_dif_vecs_tmp = dif_vec_list["sum_dif_vecs"];
      sum_dif_vecs = sum_dif_vecs_tmp;
      NumericVector family_weights_tmp = dif_vec_list["family_weights"];
      family_weights = family_weights_tmp;
      double invsum_family_weights_tmp = dif_vec_list["invsum_family_weights"];
      invsum_family_weights = invsum_family_weights_tmp;
      IntegerVector informative_families_tmp = dif_vec_list["informative_families"];
      informative_families = informative_families_tmp;
      int n_pos_tmp = dif_vec_list["n_pos"];
      n_pos = n_pos_tmp;
      int n_neg_tmp = dif_vec_list["n_neg"];
      n_neg = n_neg_tmp;
      LogicalVector pos_risk_tmp = dif_vec_list["pos_risk"];
      pos_risk = pos_risk_tmp;
      LogicalVector neg_risk_tmp = dif_vec_list["neg_risk"];
      neg_risk = neg_risk_tmp;
      LogicalVector case_high_risk_tmp = dif_vec_list["case_high_risk"];
      case_high_risk = case_high_risk_tmp;
      LogicalVector comp_high_risk_tmp = dif_vec_list["comp_high_risk"];
      comp_high_risk = comp_high_risk_tmp;

    }
  }

  // count the number of risk alleles in those with the full risk set
  double total_case_high_risk_alleles = 0;
  double total_comp_high_risk_alleles = 0;
  IntegerVector case_high_risk_int = informative_families[case_high_risk];
  IntegerVector comp_high_risk_int = informative_families[comp_high_risk];

  if (n_pos > 0){

    IntegerVector pos_cols = target_snps[pos_risk];

    // high risk case alleles
    IntegerVector case_pos_alleles = sub_colsums(case_genetic_data, case_high_risk_int, pos_cols);
    total_case_high_risk_alleles += sum(case_pos_alleles);

    // high risk comp alleles
    IntegerVector comp_pos_alleles = sub_colsums(complement_genetic_data, comp_high_risk_int, pos_cols);
    total_comp_high_risk_alleles += sum(comp_pos_alleles);

  }
  if (n_neg > 0){

    IntegerVector neg_cols = target_snps[neg_risk];

    // high risk case alleles
    IntegerVector case_neg_alleles = sub_colsums_2minus(case_genetic_data, case_high_risk_int, neg_cols);
    total_case_high_risk_alleles += sum(case_neg_alleles);

    // high risk comp alleles
    IntegerVector comp_neg_alleles = sub_colsums_2minus(complement_genetic_data, comp_high_risk_int, neg_cols);
    total_comp_high_risk_alleles += sum(comp_neg_alleles);

  }

  // compute scaling factor
  double rr;
  if ( (total_case_high_risk_alleles == 0 & total_comp_high_risk_alleles == 0) |
       R_isnancpp(total_case_high_risk_alleles) | R_isnancpp(total_comp_high_risk_alleles)){
    rr = pow(10, -10);
  } else {
    rr = total_case_high_risk_alleles/(total_case_high_risk_alleles + total_comp_high_risk_alleles);
  }

  // compute pseudo hotelling t2
 // NumericVector mu_hat = sum_dif_vecs;
//  NumericMatrix mu_hat_mat(n_informative_families, n_target);
  //for (int i = 0; i < n_informative_families; i++){

  //  NumericMatrix::Row this_row = mu_hat_mat(i, _);
//    this_row = mu_hat;
//  }

  n_informative_families = informative_families.length();
  arma::rowvec mu_hat = as<arma::rowvec>(sum_dif_vecs);
  arma::mat mu_hat_mat(n_informative_families, n_target);
  for (int i = 0; i < n_informative_families; i++){
    mu_hat_mat.row(i) = mu_hat;
  }
  arma::mat x = as<arma::mat>(subset_matrix(cases_minus_complements, informative_families, target_snps));
  arma::vec inf_family_weights = as<arma::vec>(family_weights);
  //NumericVector inf_family_weights = family_weights[family_weights > 0];
  //NumericMatrix x_minus_mu_hat = subtract_mat_n(as<NumericMatrix>(x), mu_hat_mat);
  //NumericMatrix weighted_x_minus_mu_hat = mult_rows_by_scalars(x_minus_mu_hat, inf_family_weights);
  arma::mat x_minus_mu_hat = x - mu_hat_mat;
  arma::mat weighted_x_minus_mu_hat = x_minus_mu_hat;
  weighted_x_minus_mu_hat.each_col() %= inf_family_weights;
  arma::mat cov_mat = invsum_family_weights * trans(weighted_x_minus_mu_hat) * x_minus_mu_hat;

  // set cov elements to zero if SNPs are not in same ld block
  for (int i = 0; i < n_target; i++){
    int this_row = target_snps[i] - 1;
    for (int j = 0; j < n_target; j++){
      int this_col = target_snps[j] - 1;
      if (block_ld_mat(this_row, this_col) != 1){

        cov_mat(i,j) = 0;

      }
    }
  }

  // get info for function output
  NumericVector elem_vars = wrap(sqrt(cov_mat.diag()));
  sum_dif_vecs = sum_dif_vecs/elem_vars;

  // if no variance, just make the element small
  sum_dif_vecs[elem_vars == 0] = pow(10, -10);

  // compute svd of cov_mat
  arma::mat U, V;
  arma::vec S;
  arma::svd(U, S, V, cov_mat, "dc");
  double lower_limit = std::numeric_limits<double>::epsilon();
  arma::uvec zero_sv = find(S < sqrt(lower_limit));
  double sv_rep_val = pow(10, 10);
  S.elem(zero_sv).fill(sv_rep_val);

  // compute fitness score
  double pseudo_t2 = (1/(1000*invsum_family_weights)) * arma::sum(pow((mu_hat * U), 2) / trans(S));
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


// [[Rcpp::export]]
List chrom_fitness_list(IntegerMatrix case_genetic_data, IntegerMatrix complement_genetic_data, IntegerMatrix case_comp_differences,
                        List chromosome_list, IntegerMatrix cases_minus_complements, IntegerMatrix both_one_mat,
                        LogicalMatrix block_ld_mat, NumericVector weight_lookup,
                        int n_different_snps_weight = 2, int n_both_one_weight = 1,
                        int n_case_high_risk_thresh = 20, double outlier_sd = 2.5, bool epi_test = false){

  List scores = chromosome_list.length();
  for (int i = 0; i < chromosome_list.length(); i++){

    IntegerVector target_snps = chromosome_list[i];
    scores[i] = chrom_fitness_score(case_genetic_data, complement_genetic_data, case_comp_differences,
                                    target_snps, cases_minus_complements, both_one_mat,
                                    block_ld_mat, weight_lookup, n_different_snps_weight,
                                    n_both_one_weight, n_case_high_risk_thresh, outlier_sd, epi_test);

  }
  return(scores);

}




