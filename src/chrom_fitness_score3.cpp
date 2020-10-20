#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


// [[Rcpp::export]]
List unique_chrom_list(List& chromosome_list, int chrom_size) {

  int in_size = chromosome_list.size();

  List out_list(chromosome_list);
  int s = 1;
  for (int i = 1; i < in_size; i++){

    NumericVector xi = chromosome_list[i];
    int l = 0;

    for (int j = 0; j < s; j++){

      NumericVector xj = chromosome_list[j];

      if ((sum(xi == xj) == chrom_size)){

        break;

      } else {

        l++;

      }
    }

    if (l == s){

      out_list[s] = xi;
      s++;

    }
  }
  return head(out_list, s);
}

// [[Rcpp::export]]
Rcpp::IntegerMatrix subset_matrix_cols(Rcpp::IntegerMatrix in_matrix, Rcpp::IntegerVector cols){

  int n_rows = in_matrix.nrow();
  int n_cols = cols.length();
  Rcpp::IntegerMatrix out_matrix(n_rows, n_cols);

  for (int i = 0; i < n_cols; i++){
    Rcpp::IntegerMatrix::Column original_col = in_matrix(_, cols[i]-1);
    Rcpp::IntegerMatrix::Column new_col = out_matrix(_, i);
    new_col = original_col;
  }
  return(out_matrix);
}


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
  NumericVector sum_dif_vecs = weighted_sub_colsums(cases_minus_complements, informative_families,
                                                    target_snps, family_weights);
  sum_dif_vecs = sum_dif_vecs * invsum_family_weights;

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
                         IntegerVector target_snps_in, IntegerMatrix cases_minus_complements_in, IntegerMatrix both_one_mat_in,
                         LogicalMatrix block_ld_mat, NumericVector weight_lookup,
                         int n_different_snps_weight = 2, int n_both_one_weight = 1,
                         int n_case_high_risk_thresh = 20, double outlier_sd = 2.5, bool epi_test = false) {

  // need to deep copy inputs to avoid overwriting
  IntegerMatrix case_genetic_data = subset_matrix_cols(case_genetic_data_in, target_snps_in);
  IntegerMatrix complement_genetic_data = subset_matrix_cols(complement_genetic_data_in, target_snps_in);
  IntegerMatrix case_comp_differences = subset_matrix_cols(case_comp_differences_in, target_snps_in);
  IntegerMatrix cases_minus_complements = subset_matrix_cols(cases_minus_complements_in, target_snps_in);
  IntegerMatrix both_one_mat = subset_matrix_cols(both_one_mat_in, target_snps_in);

  // now redefining target snps
  int n_target = target_snps_in.length();
  IntegerVector target_snps = seq_len(n_target);

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
              bool change_difs = false;

              // recode cases
              if (case_genetic_data(family_row, this_col) == 1){

                case_genetic_data(family_row, this_col) = 0;
                change_difs = true;

              }
              //recode complements
              if (complement_genetic_data(family_row, this_col) == 1){

                complement_genetic_data(family_row, this_col) = 0;
                if (!change_difs){
                  change_difs = true;
                }

              }
              // recode both one
              if (both_one_mat(family_row, this_col) == 1){

                both_one_mat(family_row, this_col) = 0;

              }

              //recode cases minus complements
              if (change_difs){

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
              bool change_difs = false;

              // recode cases
              if (case_genetic_data(family_row, this_col) == 1){

                case_genetic_data(family_row, this_col) = 2;
                change_difs = true;

              }
              //recode complements
              if (complement_genetic_data(family_row, this_col) == 1){

                complement_genetic_data(family_row, this_col) = 2;
                if (!change_difs){
                  change_difs = true;
                }

              }
              // recode both one
              if (both_one_mat(family_row, this_col) == 1){

                both_one_mat(family_row, this_col) = 0;

              }

              //recode cases minus complements
              if (change_difs){

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
  if ( (total_case_high_risk_alleles == 0) |
       R_isnancpp(total_case_high_risk_alleles) | R_isnancpp(total_comp_high_risk_alleles)){
    rr = pow(10, -10);
  } else {
    rr = total_case_high_risk_alleles/(total_case_high_risk_alleles + total_comp_high_risk_alleles);
  }

  // compute pseudo hotelling t2
  n_informative_families = informative_families.length();
  arma::rowvec mu_hat = as<arma::rowvec>(sum_dif_vecs);
  arma::mat mu_hat_mat(n_informative_families, n_target);
  for (int i = 0; i < n_informative_families; i++){
    mu_hat_mat.row(i) = mu_hat;
  }
  arma::mat x = as<arma::mat>(subset_matrix(cases_minus_complements, informative_families, target_snps));
  arma::vec inf_family_weights = as<arma::vec>(family_weights);
  arma::mat x_minus_mu_hat = x - mu_hat_mat;
  arma::mat weighted_x_minus_mu_hat = x_minus_mu_hat;
  weighted_x_minus_mu_hat.each_col() %= inf_family_weights;
  arma::mat cov_mat = invsum_family_weights * trans(weighted_x_minus_mu_hat) * x_minus_mu_hat;

  // set cov elements to zero if SNPs are not in same ld block
  for (int i = 0; i < n_target; i++){
    int this_row = target_snps_in[i] - 1;
    for (int j = 0; j < n_target; j++){
      int this_col = target_snps_in[j] - 1;
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
                            Named("risk_set_alleles") = risk_set_alleles,
                            Named("family_weights") = family_weights);
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

// [[Rcpp::export]]
List evolve_island(int n_migrations = 20, IntegerMatrix case_genetic_data, IntegerMatrix complement_genetic_data,
                   IntegerMatrix case_comp_different, IntegerMatrix case_minus_comp, IntegerMatrix both_one_mat,
                   IntegerMatrix block_ld_mat, int n_chromosomes, int n_candidate_snps, int chromosome_size,
                   int start_generation, NumericVector snp_chisq, IntegerVector original_col_numbers,
                   bool all_converged = false, int n_different_snps_weight = 2, int n_both_one_weight = 1,
                   IntegerVector weight_lookup, int migration_interval = 50, int gen_same_fitness = 50,
                   int max_generations = 500, double tol = 10^-6, int n_top_chroms = 100,
                   bool initial_sample_duplicates = false, string snp_sampling_type = "chisq",
                   double crossover_prop = 0.8, List chromosome_list = NULL, NumericMatrix fitness_score_mat = NULL,
                   NumericVector top_fitness = NULL, bool last_gens_equal = NULL,
                   CharacterVector top_generation_chromosome = NULL, List chromosome_mat_list = NULL,
                   List sum_dif_vec_list = NULL, List risk_allele_vec_list = NULL, int n_case_high_risk_thresh = 20,
                   double outlier_sd = 2.5){

  // initialize groups of candidate solutions if generation 1
  int generation = start_generation;
  int generations = min(NumericVector::create(start_generation + migration_interval - 1, max_generations));

  // initiate chromosome list if this is the first generation
  if (generation == 1) {

    int n_possible_unique_combn = n_chromosomes * chromosome_size;
    if ((case_genetic_data.ncol() < n_possible_unique_combn) & !initial_sample_duplicates) {

      Rcout << "Not enough SNPs present to allow for no initial sample duplicate SNPs, now allowing initial sample duplicate snps.\n";
      initial_sample_duplicates = true;

    }
    List chromosome_list(n_chromosomes);
    IntegerVector all_snps_idx = seq_len(case_genetic_data.ncol());

    for (int i = 0; i < n_chromosomes; i++) {

      IntegerVector snp_idx = sort(sample(all_snps_idx, chromosome_size, false));
      if (!initial_sample_duplicates) {

        LogicalVector not_these = !in(all_snps_idx, snp_idx);
        all_snps_idx = all_snps_idx[not_these];

      }
      chromosome_list[i] = snp_idx;

    }

    NumericMatrix fitness_score_mat(max_generations, n_chromosomes);
    std::fill(fitness_score_mat(), fitness_score_mat.end(), NumericVector::get_na());
    NumericVector top_fitness(max_generations);
    bool last_gens_equal = false;
    List top_generation_chromosome(max_generations);
    List chromosome_mat_list(max_generations);
    List sum_dif_vec_list(max_generations);
    List risk_allele_vec_list(max_generations);

  }

  // iterate over generations
  while (generation <= generations & !all_converged) {

    // 1. compute the fitness score for each set of candidate snps
    List fitness_score_list = chrom_fitness_list(case_genetic_data, complement_genetic_data, case_comp_different,
                                                 chromosome_list, case_minus_comp, both_one_mat, block_ld_mat, weight_lookup,
                                                 n_different_snps_weight, n_both_one_weight, n_case_high_risk_thresh, outlier_sd);
    NumericVector fitness_scores(n_chromosomes);
    NumericMatrix sum_dif_vecs(n_chromosomes, chromosome_size);
    IntegerMatrix gen_original_cols(n_chromosomes, chromosome_size);
    CharacterMatrix risk_allele_vecs(n_chromosomes, chromosome_size);

    for (int i = 0; i < n_chromosomes; i++){

      fitness_scores[i] = fitness_score_list[i]["fitness_score"];
      sum_dif_vecs.row(i) = fitness_score_list[i]["sum_dif_vecs"];
      risk_allele_vecs.row(i) = fitness_score_list[i]["risk_set_alleles"]
      IntegerVector chrom_col_idx = chromosome_list[i];
      gen_original_cols.row(i) = original_col_numbers[chrom_col_idx - 1];

    }

    // store the fitness scores, elements (snps) of the chromosomes, sum of the difference vectors
    fitness_score_mat.row(generation) = fitness_scores;
    chromosome_mat_list[generation] = gen_original_cols;
    sum_dif_vec_list[generation] = sum_dif_vecs;
    risk_allele_vec_list[generation] = risk_allele_vecs;

    // 2. identify the top scoring candidate solution(s) and fitness score
    double max_fitness = max(fitness_scores);
    LogicalVector top_chromosome_idx = fitness_scores == max_fitness;
    LogicalVector lower_chromosome_idx = fitness_scores == max_fitness;
    IntegerVector top_chromosomes = chromosome_list[top_chromosome_idx];
    IntegerVector lower_chromosomes = chromosome_list[lower_chromosome_idx];

    // if we have the same chromosome multiple times, just take one of them and put the duplicates in
    // the pool to be resampled
    int n_top_chroms = top_chromosomes.length();
    if (n_top_chroms  > 1) {

      int top_chromosome = top_chromosomes[0];
      IntegerVector duplicate_idx = seq(1, n_top_chroms - 1);
      IntegerVector duplicate_top_chromosomes = top_chromosomes[duplicate_idx];
      for (int i = 0, i < duplicate_top_chromosomes.length(); i++){

        int dup_top_chrom = duplicate_top_chromosomes[i];
        lower_chromosomes.push_back(dup_top_chrom);

      }

    } else {

      int top_chromosome = top.chromosomes[0];

    }


    // 3. Sample with replacement from the existing chromosomes, allowing the top
    // scoring chromosome to be sampled, but only sample from the unique chromosomes available
    Function d("duplicated");
    LogicalVector dont_sample_these = d(chromosome_list);
    LogicalVector sample_these = !(dont_sample_these);
    IntegerVector chrom_list_idx = seq_len(n_chromosomes);
    IntegerVector which_sample_these = chrom_list_idx[sample_these];
    NumericVector sample_these_scores = fitness_scores[sample_these];
    IntegerVector sampled_lower_idx = sample(which_sample_these, lower_chromosomes.length(), true, sample_these_scores);
    List sampled_lower_chromosomes = chromosome_list[sampled_lower_idx];
    NumericMatrix sampled_lower_dif_vecs = subset_matrix_rows(sum_dif_vecs, sampled_lower_idx);
    NumericVector sampled_lower_fitness_scores = fitness_scores[sampled_lower_idx];

    // 4. Determine whether each lower chromosome will be subject to mutation or crossing over

    // only allowing cross-overs between distinct chromosomes (i.e, if a chromosome was sampled twice,
    // it can't cross over with itself)

    unique.lower.idx <- unique(sampled.lower.idx)
      cross.overs <- rep(FALSE, length(unique.lower.idx))

      if (round(length(unique.lower.idx) * crossover.prop)%%2 == 0 | length(unique.lower.idx) ==
          1) {

        cross.overs[sample(seq_len(length(unique.lower.idx)), size = round(length(unique.lower.idx) *
          crossover.prop))] <- TRUE

      } else {

        cross.overs[sample(seq_len(length(unique.lower.idx)), size = (round(length(unique.lower.idx) *
          crossover.prop) + 1))] <- TRUE

      }

# those not getting crossover will be mutated
      mutations <- !cross.overs

### 6. Execute crossing over for the relevant chromosomes ### print('Step 6/9')
        if (sum(cross.overs) > 1) {

          cross.over.positions <- match(unique.lower.idx[cross.overs], sampled.lower.idx)
          cross.over.starts <- seq(1, length(cross.over.positions), by = 2)
          for (i in cross.over.starts) {

# grab pair of chromosomes to cross over
            chrom1 <- sampled.lower.chromosomes[[cross.over.positions[i]]]
            chrom1.dif.vecs <- sampled.lower.dif.vecs[cross.over.positions[i], ]
            chrom1.fitness.score <- sampled.lower.fitness.scores[cross.over.positions[i]]

            chrom2 <- sampled.lower.chromosomes[[cross.over.positions[i + 1]]]
            chrom2.dif.vecs <- sampled.lower.dif.vecs[cross.over.positions[i + 1], ]
            chrom2.fitness.score <- sampled.lower.fitness.scores[cross.over.positions[i + 1]]

# check for overlapping snps
            c1.c2.matching.snp.positions <- match(chrom1, chrom2)
              c1.c2.matching.snp.positions <- c1.c2.matching.snp.positions[!is.na(c1.c2.matching.snp.positions)]
            c1.c2.not.matching.snp.positions <- setdiff(seq_len(chromosome.size), c1.c2.matching.snp.positions)

              c2.c1.matching.snp.positions <- match(chrom2, chrom1)
              c2.c1.matching.snp.positions <- c2.c1.matching.snp.positions[!is.na(c2.c1.matching.snp.positions)]
            c2.c1.not.matching.snp.positions <- setdiff(seq_len(chromosome.size), c2.c1.matching.snp.positions)

# for the non-matching snps, order by the decreasing magnitude of the difference vector in the
# chromsome with the higher fitness score and order by the increasing magntidue of the difference
# vector in the chromosome with the lower fitness score **ultimately will be used to substitute the
# higher magnitude elements in the lower scoring chromosome for the lower magnitude elements in the
# higher scoring chromosome**
              if (chrom1.fitness.score >= chrom2.fitness.score) {

                c1.c2.not.matching.snp.positions <- c1.c2.not.matching.snp.positions[order(abs(chrom2.dif.vecs[c1.c2.not.matching.snp.positions]))]
                c2.c1.not.matching.snp.positions <- c2.c1.not.matching.snp.positions[order(abs(chrom1.dif.vecs[c2.c1.not.matching.snp.positions]),
                                                                                           decreasing = TRUE)]

              } else {

                c1.c2.not.matching.snp.positions <- c1.c2.not.matching.snp.positions[order(abs(chrom2.dif.vecs[c1.c2.not.matching.snp.positions]),
                                                                                           decreasing = TRUE)]
                c2.c1.not.matching.snp.positions <- c2.c1.not.matching.snp.positions[order(abs(chrom1.dif.vecs[c2.c1.not.matching.snp.positions]))]

              }

# order the chromosomes, first by the overlapping snps and the non-overlapping
              chrom2 <- chrom2[c(c1.c2.matching.snp.positions, c1.c2.not.matching.snp.positions)]
              chrom1 <- chrom1[c(c2.c1.matching.snp.positions, c2.c1.not.matching.snp.positions)]

# determine how many snps could be crossed over and also make sure we don't simply swap chromosomes
              possible.cut.points <- sort(which(chrom1 != chrom2), decreasing = TRUE)
                if (length(possible.cut.points) == chromosome.size) {

                  n.possible.crosses <- chromosome.size - 1

                } else {

                  n.possible.crosses <- length(possible.cut.points)

                }

# determine how many snps will actually be crossed over
                n.crosses <- sample.int(n.possible.crosses, 1)

# pick out their positions
                  cross.points <- c(seq_len(chromosome.size))[(chromosome.size - n.crosses + 1):chromosome.size]

# exchange the high magnitude elements from the lower scoring chromosome with the low magnitude
# elements from the high scoring chromsosome
                  chrom1.cross <- chrom1
                    chrom1.cross[cross.points] <- chrom2[cross.points]
                  chrom2.cross <- chrom2
                    chrom2.cross[cross.points] <- chrom1[cross.points]

# replace in the chromosome list
                  sampled.lower.chromosomes[[cross.over.positions[i]]] <- chrom1.cross
                    sampled.lower.chromosomes[[cross.over.positions[i + 1]]] <- chrom2.cross

          }

        } else {

          cross.over.positions <- 0

        }

### 7. Mutate the chromosomes that were not crossed over ### print('Step 7/9')
        mutation.positions <- (seq_len(length(sampled.lower.chromosomes)))
          snps.for.mutation <- sample(seq_len(n.candidate.snps), n.candidate.snps, prob = snp.chisq,
                                      replace = TRUE)
            ulen <- length(unique(snps.for.mutation))

            for (i in mutation.positions) {

# grab the chromosome and its difference vector
              target.chrom <- sampled.lower.chromosomes[[i]]
              target.dif.vec <- as.vector(t(sampled.lower.dif.vecs[i, ]))

# sort the chromosome elements from lowest absolute difference vector to highest
                target.chrom <- target.chrom[order(abs(target.dif.vec))]

# determine which snps to mutate
                total.mutations <- min(max(1, rbinom(1, chromosome.size, 0.5)), ulen)

# get random mutation set
                  mutated.snps <- sample(snps.for.mutation, total.mutations)

# check for duplicates versus chromosome or within sample. If so, re-sample with exclusion checks (slow)
                    if (length(mutated.snps) != length(unique(mutated.snps)) || any(mutated.snps %in% target.chrom)) {

                      possible.snps.for.mutation <- snps.for.mutation[!snps.for.mutation %in% target.chrom]
                      ulen.psm <- length(unique(possible.snps.for.mutation))
                        if (ulen.psm < total.mutations){

                          total.mutations <- ulen.psm

                        }
                        mutated.snps <- rep(NA, total.mutations)

# execute mutations
                          for (j in seq_len(total.mutations)) {
                            sampled.snp <- sample(possible.snps.for.mutation, 1)
                            mutated.snps[j] <- sampled.snp
                            possible.snps.for.mutation <- possible.snps.for.mutation[possible.snps.for.mutation !=
                              sampled.snp]
                          }
                    }

# substitute in mutations
                    target.chrom[1:total.mutations] <- mutated.snps
                      sampled.lower.chromosomes[[i]] <- target.chrom

            }

### 8. Combine into new population (i.e., the final collection of chromosomes for the next
### generation) print('Step 8/9')
            chromosome.list <- lapply(c(top.chromosome, sampled.lower.chromosomes), sort, method = "radix")

### 9.Increment Iterators print('Step 9/9')
              top.fitness[generation] <- max.fitness
                top.generation.chromosome[[generation]] <- original.col.numbers[top.chromosome[[1]]]
# print(paste0('Max fitness score:', max.fitness)) print('Top Chromosome(s):')
# print(original.col.numbers[top.chromosome[[1]]])
              if (generation >= gen.same.fitness) {

# check to see if enough of the last generations have had the same top chromosome to terminate
                last.gens <- top.fitness[(generation - (gen.same.fitness - 1)):generation]
                last.gens.equal <- abs(max(last.gens) - min(last.gens)) < tol
                  if (is.null(n.migrations) & last.gens.equal) {

                    all.converged <- TRUE

                  }

              }

              generation <- generation + 1

  }
### If the algorithm hasn't hit the max number of generations or converged, return a partial list of
### the results ###
  if (generation < max.generations & !all.converged & !is.null(n.migrations)) {

# pick out the top fitness scores, and bottom fitness scores
    fitness.score.list <- lapply(seq_len(length(chromosome.list)), function(x) {

      chrom.fitness.score(case.genetic.data, complement.genetic.data, case.comp.different, chromosome.list[[x]],
                          case.minus.comp, both.one.mat, block.ld.mat, weight.lookup, n.different.snps.weight, n.both.one.weight,
                          n.case.high.risk.thresh, outlier.sd)

    })
      fitness.scores <- vapply(fitness.score.list, function(x) x$fitness.score, 1.0)
      chromosome.list <- chromosome.list[order(fitness.scores, decreasing = TRUE)]

### identify the chromosomes that will migrate to other islands ###
    migrations <- chromosome.list[seq_len(n.migrations)]

### remove the lowest scoring chromosomes in preparation for migration ###
    chromosome.list <- chromosome.list[seq_len((length(chromosome.list) - n.migrations))]

### return list of results ###
    return(list(migrations = migrations, chromosome.list = chromosome.list, fitness.score.mat = fitness.score.mat,
                top.fitness = top.fitness, last.gens.equal = last.gens.equal, top.generation.chromosome = top.generation.chromosome,
                chromosome.mat.list = chromosome.mat.list, sum.dif.vec.list = sum.dif.vec.list, risk.allele.vec.list = risk.allele.vec.list,
                generation = generation))

  } else {

### Otherwise Return the best chromosomes, their fitness scores, difference vectors and the number of
### generations ###
    last.generation <- generation - 1
    all.chrom.dt <- rbindlist(chromosome.mat.list[seq_len(last.generation)], use.names = FALSE)
      all.chrom.dif.vec.dt <- rbindlist(sum.dif.vec.list[seq_len(last.generation)], use.names = FALSE)
      all.chrom.risk.allele.vec.dt <- rbindlist(risk.allele.vec.list[seq_len(last.generation)], use.names = FALSE)
      unique.chromosome.dt <- unique(all.chrom.dt)
      colnames(unique.chromosome.dt) <- paste0("snp", seq_len(ncol(unique.chromosome.dt)))
      unique.chrom.dif.vec.dt <- all.chrom.dif.vec.dt[!duplicated(all.chrom.dt), ]
    unique.chrom.risk.allele.vec.dt <- all.chrom.risk.allele.vec.dt[!duplicated(all.chrom.dt), ]
    colnames(unique.chrom.dif.vec.dt) <- paste0("snp", seq_len(ncol(unique.chrom.dif.vec.dt)),
             ".diff.vec")
      colnames(unique.chrom.risk.allele.vec.dt) <- paste0("snp", seq_len(ncol(unique.chrom.dif.vec.dt)),
               ".allele.copies")
      unique.fitness.score.vec <- as.vector(t(fitness.score.mat[seq_len(last.generation), ]))[!duplicated(all.chrom.dt)]
    unique.results <- cbind(unique.chromosome.dt, unique.chrom.dif.vec.dt, unique.chrom.risk.allele.vec.dt)
      unique.results[, `:=`(raw.fitness.score, unique.fitness.score.vec)]
    unique.results[, `:=`(min.elem, min(abs(.SD))), by = seq_len(nrow(unique.results)), .SDcols = (1 +
      chromosome.size):(2 * chromosome.size)]
    unique.results[, `:=`(fitness.score, min.elem * raw.fitness.score)]
    setorder(unique.results, -fitness.score)
      final.result <- unique.results[seq_len(n.top.chroms), ]

# print(paste('Algorithm terminated after', last.generation, 'generations.'))
    res.list <- list(top.chromosome.results = final.result, n.generations = last.generation)
      return(res.list)

  }
}



