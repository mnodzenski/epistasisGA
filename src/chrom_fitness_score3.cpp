#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// HELPER FUNCTIONS

// [[Rcpp::export]]
int scalar_min(int x, int y){

  if ( x <= y){

    return(x);

  } else {

    return(y);

  }

}

// [[Rcpp::export]]
int scalar_max(int x, int y){

  if ( x >= y){

    return(x);

  } else {

    return(y);

  }

}

// [[Rcpp::export]]
IntegerVector concat(IntegerVector x, IntegerVector y){

  IntegerVector out_vec(x.length() + y.length());
  IntegerVector x_pos = seq_len(x.length());
  out_vec[x_pos - 1] = x;
  IntegerVector y_pos = seq(x.length() + 1, out_vec.length());
  out_vec[y_pos - 1] = y;
  return(out_vec);

}

// [[Rcpp::export]]
IntegerVector sort_by_order(IntegerVector x, NumericVector y, int sort_type) {

  //sort_type 1 for ascending, 2 for descending
  IntegerVector idx = seq_along(x) - 1;
  if (sort_type == 1){

    std::sort(idx.begin(), idx.end(), [&](int i, int j){return y[i] < y[j];});

  } else if (sort_type == 2){

    std::sort(idx.begin(), idx.end(), [&](int i, int j){return y[i] > y[j];});

  }

  return x[idx];
}

// [[Rcpp::export]]
IntegerVector seq_by2(int l){

  IntegerVector out_vec(l);
  out_vec[0] = 1;
  for (int i = 1; i < l; i++){

    out_vec[i] = out_vec[i - 1] + 2;

  }
  return(out_vec);

}

// [[Rcpp::export]]
LogicalVector unique_chrom_list(List chromosome_list, int chrom_size) {

  int in_size = chromosome_list.size();
  LogicalVector out_vec(in_size, false);
  out_vec[0] = true;
  List comp_these = chromosome_list[out_vec];

  int s = 1;
  for (int i = 1; i < in_size; i++){

    NumericVector xi = chromosome_list[i];
    int l = 0;

    comp_these = chromosome_list[out_vec];

    for (int j = 0; j < comp_these.length(); j++){

      NumericVector xj = comp_these[j];

      if ((sum(xi == xj) == chrom_size)){

        break;

      } else {

        l++;

      }
    }

    if (l == s){

      out_vec[i] = true;
      s++;

    }
  }
  return out_vec;
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

// function to pick out the families with case or complement with the full risk set

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

// function to compute the difference vectors in computing the fitness score

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

// actually computes the fitness score

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


// helper to apply the fitness score function to a list of chromosomes

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

// function to initiate island populations
List initiate_population(IntegerMatrix case_genetic_data, int n_chromosomes, int chromosome_size,
                         int generation, NumericVector snp_chisq,
                         int max_generations = 500,
                         bool initial_sample_duplicates = false, String snp_sampling_type = "chisq",
                         Nullable<List> chromosome_list = R_NilValue,
                         Nullable<List> fitness_score_list = R_NilValue, Nullable<NumericVector> top_fitness = R_NilValue,
                         Nullable<List> top_generation_chromosome = R_NilValue,
                         Nullable<List> gen_chromosome_list = R_NilValue, Nullable<List> sum_dif_vec_list = R_NilValue,
                         Nullable<List> risk_allele_vec_list = R_NilValue){

  if (generation == 1) {

    int n_possible_unique_combn = n_chromosomes * chromosome_size;
    if ((case_genetic_data.ncol() < n_possible_unique_combn) & !initial_sample_duplicates) {

      Rcout << "Not enough SNPs present to allow for no initial sample duplicate SNPs, now allowing initial sample duplicate snps.\n";
      initial_sample_duplicates = true;

    }
    List chromosome_list(n_chromosomes);
    IntegerVector all_snps_idx = seq_len(case_genetic_data.ncol());

    for (int i = 0; i < n_chromosomes; i++) {

      IntegerVector snp_idx = sample(all_snps_idx, chromosome_size, false);
      std::sort(snp_idx.begin(), snp_idx.end());
      if (!initial_sample_duplicates) {

        LogicalVector not_these = !in(all_snps_idx, snp_idx);
        all_snps_idx = all_snps_idx[not_these];

      }
      chromosome_list[i] = snp_idx;

    }

    List fitness_score_list(max_generations);
    NumericVector top_fitness(max_generations);
    List top_generation_chromosome(max_generations);
    List gen_chromosome_list(max_generations);
    List sum_dif_vec_list(max_generations);
    List risk_allele_vec_list(max_generations);
    List res = List::create(Named("chromosome_list") = chromosome_list,
                            Named("fitness_score_list") = fitness_score_list,
                            Named("top_generation_chromosome") = top_generation_chromosome,
                            Named("gen_chromosome_list") = gen_chromosome_list,
                            Named("sum_dif_vec_list") = sum_dif_vec_list,
                            Named("risk_allele_vec_list") = risk_allele_vec_list,
                            Named("top_fitness") = top_fitness);
    return(res);

  } else {

    //otherwise initiate the objects passed as optional arguments
    List chromosome_list_out(chromosome_list);
    List fitness_score_list_out(fitness_score_list);
    List top_generation_chromosome_out(top_generation_chromosome);
    List gen_chromosome_list_out(gen_chromosome_list);
    List sum_dif_vec_list_out(sum_dif_vec_list);
    List risk_allele_vec_list_out(risk_allele_vec_list);
    NumericVector top_fitness_out(top_fitness);
    List res = List::create(Named("chromosome_list") = chromosome_list_out,
                            Named("fitness_score_list") = fitness_score_list_out,
                            Named("top_generation_chromosome") = top_generation_chromosome_out,
                            Named("gen_chromosome_list") = gen_chromosome_list_out,
                            Named("sum_dif_vec_list") = sum_dif_vec_list_out,
                            Named("risk_allele_vec_list") = risk_allele_vec_list_out,
                            Named("top_fitness") = top_fitness_out);
      return(res);
  }
}



// function to evolve island populations

// [[Rcpp::export]]
List evolve_island(int n_migrations, IntegerMatrix case_genetic_data, IntegerMatrix complement_genetic_data,
                   IntegerMatrix case_comp_different, IntegerMatrix case_minus_comp, IntegerMatrix both_one_mat,
                   LogicalMatrix block_ld_mat, int n_chromosomes, int chromosome_size, NumericVector weight_lookup,
                   int start_generation, NumericVector snp_chisq, IntegerVector original_col_numbers,
                   bool all_converged = false, int n_different_snps_weight = 2, int n_both_one_weight = 1,
                   int migration_interval = 50, int gen_same_fitness = 50,
                   int max_generations = 500, double tol = 10^-6, int n_top_chroms = 100,
                   bool initial_sample_duplicates = false, String snp_sampling_type = "chisq",
                   double crossover_prop = 0.8, Nullable<List> chromosome_list_in = R_NilValue,
                   Nullable<List> fitness_score_list_in = R_NilValue, Nullable<NumericVector> top_fitness_in = R_NilValue,
                   bool last_gens_equal = false, Nullable<List> top_generation_chromosome_in = R_NilValue,
                   Nullable<List> gen_chromosome_list_in = R_NilValue, Nullable<List> sum_dif_vec_list_in = R_NilValue,
                   Nullable<List> risk_allele_vec_list_in = R_NilValue, int n_case_high_risk_thresh = 20,
                   double outlier_sd = 2.5){

  // initialize groups of candidate solutions if generation 1
  int generation = start_generation;
  int generations = min(NumericVector::create(start_generation + migration_interval - 1, max_generations));

  // initialize populations, or grab input population data
  List population = initiate_population(case_genetic_data, n_chromosomes, chromosome_size, generation, snp_chisq,
                                        max_generations, initial_sample_duplicates, snp_sampling_type, chromosome_list_in,
                                        fitness_score_list_in, top_fitness_in, top_generation_chromosome_in, gen_chromosome_list_in,
                                        sum_dif_vec_list_in, risk_allele_vec_list_in);
  List chromosome_list = population["chromosome_list"];
  List fitness_score_list = population["fitness_score_list"];
  List top_generation_chromosome = population["top_generation_chromosome"];
  List gen_chromosome_list = population["gen_chromosome_list"];
  List sum_dif_vec_list = population["sum_dif_vec_list"];
  List risk_allele_vec_list = population["risk_allele_vec_list"];
  NumericVector top_fitness = population["top_fitness"];

  // iterate over generations
  while ((generation <= generations) & !all_converged) {

    // 1. compute the fitness score for each set of candidate snps
    List chrom_fitness_score_list = chrom_fitness_list(case_genetic_data, complement_genetic_data, case_comp_different,
                                                 chromosome_list, case_minus_comp, both_one_mat, block_ld_mat, weight_lookup,
                                                 n_different_snps_weight, n_both_one_weight, n_case_high_risk_thresh, outlier_sd);
    NumericVector fitness_scores(n_chromosomes);
    List sum_dif_vecs(n_chromosomes);
    List gen_original_cols(n_chromosomes);
    List risk_allele_vecs(n_chromosomes);

    for (int i = 0; i < n_chromosomes; i++){

      fitness_scores[i] = chrom_fitness_score_list[i]["fitness_score"];
      sum_dif_vecs[i] = chrom_fitness_score_list[i]["sum_dif_vecs"];
      risk_allele_vecs[i] = chrom_fitness_score_list[i]["risk_set_alleles"];
      IntegerVector chrom_col_idx = chromosome_list[i];
      gen_original_cols[i] = original_col_numbers[chrom_col_idx - 1];

    }

    // store the fitness scores, elements (snps) of the chromosomes, sum of the difference vectors
    fitness_score_list[generation - 1] = fitness_scores;
    gen_chromosome_list[generation - 1] = gen_original_cols;
    sum_dif_vec_list[generation - 1] = sum_dif_vecs;
    risk_allele_vec_list[generation - 1] = risk_allele_vecs;

    // 2. identify the top scoring candidate solution(s) and fitness score
    double max_fitness = max(fitness_scores);
    LogicalVector top_chromosome_idx = fitness_scores == max_fitness;
    LogicalVector lower_chromosome_idx = fitness_scores != max_fitness;
    List top_chromosomes = chromosome_list[top_chromosome_idx];
    List lower_chromosomes = chromosome_list[lower_chromosome_idx];

    // if we have the same chromosome multiple times, just take one of them and put the duplicates in
    // the pool to be resampled
    int n_top_chroms = top_chromosomes.length();
    IntegerVector top_chromosome(chromosome_size);
    if (n_top_chroms  > 1) {

      top_chromosome = top_chromosomes[0];
      IntegerVector duplicate_idx = seq(1, n_top_chroms - 1);
      List duplicate_top_chromosomes = top_chromosomes[duplicate_idx];
      List new_lower_chromosomes(lower_chromosomes.length() + duplicate_top_chromosomes.length());
      for (int i = 0; i < lower_chromosomes.length(); i++){

        new_lower_chromosomes[i] = lower_chromosomes[i];

      }
      for (int i = 0; i < duplicate_top_chromosomes.length(); i++){

        new_lower_chromosomes[i + lower_chromosomes.length()] = duplicate_top_chromosomes[i];

      }
      lower_chromosomes = new_lower_chromosomes;

    } else {

      top_chromosome = top_chromosomes[0];

    }

    // 3. Sample with replacement from the existing chromosomes, allowing the top
    // scoring chromosome to be sampled, but only sample from the unique chromosomes available
    LogicalVector sample_these = unique_chrom_list(chromosome_list, chromosome_size);
    IntegerVector chrom_list_idx = seq_len(n_chromosomes);
    IntegerVector which_sample_these = chrom_list_idx[sample_these];
    NumericVector sample_these_scores = fitness_scores[sample_these];
    IntegerVector sampled_lower_idx = sample(which_sample_these, lower_chromosomes.length(), true, sample_these_scores);
    List sampled_lower_chromosomes = chromosome_list[sampled_lower_idx - 1];
    List sampled_lower_dif_vecs = sum_dif_vecs[sampled_lower_idx - 1];
    NumericVector sampled_lower_fitness_scores = fitness_scores[sampled_lower_idx - 1];

    // 4. Determine whether each lower chromosome will be subject to mutation or crossing over

    // only allowing cross-overs between distinct chromosomes (i.e, if a chromosome was sampled twice,
    // it can't cross over with itself)
    IntegerVector unique_lower_idx = unique(sampled_lower_idx);
    LogicalVector cross_overs(unique_lower_idx.length(), false);
    int rounded_crosses = round(unique_lower_idx.length() * crossover_prop);
    IntegerVector possible_crosses = seq_len(unique_lower_idx.length());
    int n_crosses = 0;
    if ((rounded_crosses % 2 == 0) | (unique_lower_idx.length() == 1)){

      n_crosses = round(unique_lower_idx.length() * crossover_prop);

    } else {

      n_crosses = round(unique_lower_idx.length() * crossover_prop) + 1;

    }
    IntegerVector cross_these = sample(possible_crosses, n_crosses, false);
    cross_overs[cross_these - 1] = true;

    // 5. Execute crossing over for the relevant chromosomes
    IntegerVector cross_over_positions(1, 0);
    if (sum(cross_overs) > 1){

      IntegerVector lower_idx_cross_overs = unique_lower_idx[cross_overs];
      cross_over_positions = match(lower_idx_cross_overs, sampled_lower_idx);
      int this_length = cross_over_positions.length() / 2;
      IntegerVector crossover_starts = seq_by2(this_length);

      for (int k = 0; k < crossover_starts.length(); k++){

        // grab pair of chromosomes to cross over
        int chrom1_sub_idx = crossover_starts[k];
        int chrom1_idx = cross_over_positions[chrom1_sub_idx - 1];
        IntegerVector chrom1 = sampled_lower_chromosomes[chrom1_idx - 1];
        NumericVector chrom1_dif_vecs = sampled_lower_dif_vecs[chrom1_idx - 1];
        double chrom1_fitness_score = sampled_lower_fitness_scores[chrom1_idx - 1];

        int chrom2_idx = cross_over_positions[chrom1_sub_idx];
        IntegerVector chrom2 = sampled_lower_chromosomes[chrom2_idx - 1];
        NumericVector chrom2_dif_vecs = sampled_lower_dif_vecs[chrom2_idx - 1];
        double chrom2_fitness_score = sampled_lower_fitness_scores[chrom2_idx - 1];

        //check for overlapping snps
        IntegerVector chrom_positions = seq_len(chromosome_size);
        IntegerVector c1_c2_matching_snp_positions = match(chrom1, chrom2);
        c1_c2_matching_snp_positions = c1_c2_matching_snp_positions[!is_na(c1_c2_matching_snp_positions)];
        IntegerVector c1_c2_not_matching_snp_positions = setdiff(chrom_positions, c1_c2_matching_snp_positions);

        IntegerVector c2_c1_matching_snp_positions = match(chrom2, chrom1);
        c2_c1_matching_snp_positions = c2_c1_matching_snp_positions[!is_na(c2_c1_matching_snp_positions)];
        IntegerVector c2_c1_not_matching_snp_positions = setdiff(chrom_positions, c2_c1_matching_snp_positions);

        // for the non-matching snps, order by the decreasing magnitude of the difference vector in the
        // chromsome with the higher fitness score and order by the increasing magntidue of the difference
        // vector in the chromosome with the lower fitness score **ultimately will be used to substitute the
        // higher magnitude elements in the lower scoring chromosome for the lower magnitude elements in the
        // higher scoring chromosome**
        if (chrom1_fitness_score >= chrom2_fitness_score){

          NumericVector c2_not_matching_dif_vecs = chrom2_dif_vecs[c1_c2_not_matching_snp_positions - 1];
          IntegerVector c1_c2_order = sort_by_order(chrom_positions, abs(c2_not_matching_dif_vecs), 1);
          c1_c2_not_matching_snp_positions = c1_c2_not_matching_snp_positions[c1_c2_order - 1];

          NumericVector c1_not_matching_dif_vecs = chrom1_dif_vecs[c2_c1_not_matching_snp_positions - 1];
          IntegerVector c2_c1_order = sort_by_order(chrom_positions, abs(c1_not_matching_dif_vecs), 2);
          c2_c1_not_matching_snp_positions = c2_c1_not_matching_snp_positions[c2_c1_order - 1];

        } else {

          NumericVector c2_not_matching_dif_vecs = chrom2_dif_vecs[c1_c2_not_matching_snp_positions - 1];
          IntegerVector c1_c2_order = sort_by_order(chrom_positions, abs(c2_not_matching_dif_vecs), 2);
          c1_c2_not_matching_snp_positions = c1_c2_not_matching_snp_positions[c1_c2_order - 1];

          NumericVector c1_not_matching_dif_vecs = chrom1_dif_vecs[c2_c1_not_matching_snp_positions - 1];
          IntegerVector c2_c1_order = sort_by_order(chrom_positions, abs(c1_not_matching_dif_vecs), 1);
          c2_c1_not_matching_snp_positions = c2_c1_not_matching_snp_positions[c2_c1_order - 1];

        }

        // order the chromosomes, first by overlapping snps and then non-overlapping
        IntegerVector c2_order = concat(c1_c2_matching_snp_positions, c1_c2_not_matching_snp_positions);
        chrom2 = chrom2[c2_order - 1];

        IntegerVector c1_order = concat(c2_c1_matching_snp_positions, c2_c1_not_matching_snp_positions);
        chrom1 = chrom1[c1_order - 1];

        // determine how many snps could be crossed over and also make sure we don't simply swap chromosomes
        IntegerVector possible_cut_points = chrom_positions[chrom1 != chrom2];
        std::sort(possible_cut_points.begin(), possible_cut_points.end());
        std::reverse(possible_cut_points.begin(), possible_cut_points.end());
        int n_possible_crosses = 0;
        if (possible_cut_points.length() == chromosome_size){

          n_possible_crosses = chromosome_size - 1;

        } else {

          n_possible_crosses = possible_cut_points.length();

        }
        // decide how many snps will actually be crossed
        IntegerVector sample_from = seq_len(n_possible_crosses);
        int n_crosses = sample(sample_from, 1)[0];

        // pick out their positions
        int start_here = chromosome_size - n_crosses;
        IntegerVector cross_points = seq(start_here, chromosome_size - 1);

        // exchange the high magnitude elements from the lower scoring chromosome with the low magnitude
        // elements from the high scoring chromsosome
        IntegerVector chrom1_cross = chrom1;
        chrom1_cross[cross_points] = chrom2[cross_points];
        IntegerVector chrom2_cross = chrom2;
        chrom2_cross[cross_points] = chrom1[cross_points];

        // replace in chromosome list
        sampled_lower_chromosomes[chrom1_idx - 1] = chrom1_cross;
        sampled_lower_chromosomes[chrom2_idx - 1] = chrom2_cross;

      }

    }
    // 6. mutate chromosomes
    IntegerVector mutation_positions = seq_len(sampled_lower_chromosomes.length());
    IntegerVector candiate_snp_idx = seq_len(case_genetic_data.ncol());
    IntegerVector snps_for_mutation = sample(candiate_snp_idx, candiate_snp_idx.length(),
                                             true, snp_chisq);
    int ulen = unique(snps_for_mutation).length();
    for (int l = 0; l < mutation_positions.length(); l++){

      int mutation_position = mutation_positions[l];

      // grab chromosome and its difference vector
      IntegerVector target_chrom = sampled_lower_chromosomes[mutation_position - 1];
      NumericVector target_dif_vec = sampled_lower_dif_vecs[mutation_position - 1];

      // sort the chromosome elements from lowest absolute difference vector to highest
      target_chrom = sort_by_order(target_chrom, abs(target_dif_vec), 1);

      // determine which snps to mutate
      int binom_sample = rbinom(1, chromosome_size, 0.5)[0];
      int part1 = scalar_max(1, binom_sample);
      int total_mutations = scalar_min(part1, ulen);

      // grab random mutation set
      IntegerVector mutated_snps = sample(snps_for_mutation, total_mutations, false);

      // check for duplicates and, if so, resample with exclusion checks
      if ( (mutated_snps.length() != unique(mutated_snps).length()) | (sum(in(mutated_snps, target_chrom)) > 0) ){

        IntegerVector possible_snps_for_mutation = snps_for_mutation[!in(snps_for_mutation, target_chrom)];
        int ulen_psm = unique(possible_snps_for_mutation).length();
        if (ulen_psm < total_mutations){

          total_mutations = ulen_psm;

        }
        IntegerVector new_mutated_snps(total_mutations);
        for (int p = 0; p < total_mutations; p++){

          int sampled_snp = sample(possible_snps_for_mutation, 1)[0];
          new_mutated_snps[p] = sampled_snp;
          possible_snps_for_mutation = possible_snps_for_mutation[possible_snps_for_mutation != sampled_snp];

        }
        mutated_snps = new_mutated_snps;

      }
      // substitute in mutations
      IntegerVector mutate_here = seq(0, total_mutations - 1);
      target_chrom[mutate_here] = mutated_snps;
      sampled_lower_chromosomes[mutation_position - 1] = target_chrom;
    }
    // 7. Combined into new population (i.e., the final collection of chromosomes for the next generation)
    std::sort(top_chromosome.begin(), top_chromosome.end());
    chromosome_list[0] = top_chromosome;
    for (int i = 1; i < chromosome_list.length(); i++){

      IntegerVector sampled_lower_chromosomes_i = sampled_lower_chromosomes[i - 1];
      std::sort(sampled_lower_chromosomes_i.begin(), sampled_lower_chromosomes_i.end());
      chromosome_list[i] = sampled_lower_chromosomes_i;

    }

    // 8. Increment iterators
    top_fitness[generation - 1] = max_fitness;
    top_generation_chromosome[generation - 1] = original_col_numbers[top_chromosome - 1];

    // check to see if we can terminate
    if (generation >= gen_same_fitness){

      int start_gen = generation - gen_same_fitness;
      IntegerVector check_these_gens = seq(start_gen, generation - 1);
      NumericVector last_gens = top_fitness[check_these_gens];
      double abs_diff = abs(max(last_gens) - min(last_gens));
      last_gens_equal = abs_diff < tol;
      if ( (n_migrations == 0) & last_gens_equal){

        all_converged = true;

      }

    }
    generation += 1;

  }

  // if the algorithm hasn't hit the max number of generations or converged, return a partial list of
  // the results
  if (generation < max_generations & !all_converged & n_migrations > 0){

    // pick out the highest and lowest fitness scores of the new chromosomes
    List chrom_fitness_score_list = chrom_fitness_list(case_genetic_data, complement_genetic_data, case_comp_different,
                                            chromosome_list, case_minus_comp, both_one_mat, block_ld_mat, weight_lookup,
                                            n_different_snps_weight, n_both_one_weight, n_case_high_risk_thresh, outlier_sd);
    NumericVector fitness_scores(n_chromosomes);
    for (int i = 0; i < n_chromosomes; i++){

      fitness_scores[i] = chrom_fitness_score_list[i]["fitness_score"];

    }
    IntegerVector chrom_list_order = sort_by_order(seq_len(chromosome_list), fitness_scores, 2);
    chromosome_list = chromosome_list[chrom_list_order - 1];

    // identify chromosomes that will migrate to other islands
    IntegerVector migration_idx = seq_len(n_migrations);
    IntegerVector migrations = chromosome_list[ migration_idx - 1];

    // remove the lowest scoring chromosomes in preparation for migration
    IntegerVector keep_these = seq_len(n_chromosomes - n_migrations);
    chromosome_list = chromosome_list[keep_these - 1];

    List res = List::create(Named("migrations") = migrations,
                            Named("chromosome_list") = chromosome_list,
                            Named("fitness_score_list") = fitness_score_list,
                            Named("top_fitness") = top_fitness,
                            Named("last_gens_equal") = last_gens_equal,
                            Named("top_generation_chromosome") = top_generation_chromosome,
                            Named("gen_chromosome_list") = gen_chromosome_list,
                            Named("sum_dif_vec_list") = sum_dif_vec_list,
                            Named("risk_allele_vec_list") = risk_allele_vec_list,
                            Named("generation") = generation);
    return(res);

  } else {

    // Otherwise return the best chromosomes, their fitness scores, difference vectors and the number of
    // generations
    int last_generation = generation - 1;

    // first picking out the unique results over islands/generations
    int ll = last_generation*n_chromosomes;
    List u_chrom_list(ll);
    List u_dif_vec_list(ll);
    List u_risk_allele_list(ll);
    NumericVector u_fitness_score_vec(ll);
    int counter = 0;
    for (int i = 0; i < last_generation; i++){

      List gen_chromosome_list_i = gen_chromosome_list[i];
      List sum_dif_vec_list_i = sum_dif_vec_list[i];
      List risk_allele_vec_list_i = risk_allele_vec_list[i];
      NumericVector fitness_score_vec_i = fitness_score_list[i];

      for(int j = 0; j < n_chromosomes; j++){

        u_chrom_list[counter] = gen_chromosome_list_i[j];
        u_dif_vec_list[counter] = sum_dif_vec_list_i[j];
        u_risk_allele_list[counter] = risk_allele_vec_list_i[j];
        u_fitness_score_vec[counter] = fitness_score_vec_i[j];
        counter += 1;

      }
    }

    LogicalVector unique_res_idx = unique_chrom_list(u_chrom_list, chromosome_size);
    u_chrom_list = u_chrom_list[unique_res_idx];
    u_dif_vec_list = u_dif_vec_list[unique_res_idx];
    u_risk_allele_list = u_risk_allele_list[unique_res_idx];
    u_fitness_score_vec = u_fitness_score_vec[unique_res_idx];

    // now pick out the highest scoring chromosomes
    int n_unique_chroms = sum(unique_res_idx);
    IntegerVector tmp_idx = seq_len(n_unique_chroms);
    IntegerVector top_res_order = sort_by_order(tmp_idx, u_fitness_score_vec, 2);
    int n_out_res = scalar_min(n_top_chroms, n_unique_chroms);
    IntegerVector top_res_idx = top_res_order[seq_len(n_out_res) - 1];

    u_chrom_list = u_chrom_list[top_res_idx - 1];
    u_dif_vec_list = u_dif_vec_list[top_res_idx - 1];
    u_risk_allele_list = u_risk_allele_list[top_res_idx - 1];
    out_fitness_score_vec = u_fitness_score_vec[top_res_idx - 1];

    // put together pieces of results
    IntegerMatrix out_chrom_mat(n_out_res, chromosome_size);
    NumericMatrix out_dif_vec_mat(n_out_res, chromosome_size);
    CharacterMatrix out_risk_allele_mat(n_out_res, chromosome_size);
    for (int i = 0; i < n_out_res; i++){

      out_chrom_mat(i, _) = u_chrom_list[i];
      out_dif_vec_mat(i, _) = u_dif_vec_list[i];
      out_risk_allele_mat(i, _) = u_risk_allele_list[i];

    }
    List res = List::create(Named("top_chromosomes") = out_chrom_mat,
                            Named("sum_dif_vecs") = out_dif_vec_mat,
                            Named("risk_alleles") = out_risk_allele_mat,
                            Named("raw_fitness_scores") = out_fitness_score_vec);
    return(res);

  }
}






