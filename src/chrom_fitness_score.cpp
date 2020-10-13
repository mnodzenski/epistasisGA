#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
IntegerMatrix sign_subtract_mat(IntegerMatrix x, IntegerMatrix y){

  int nrow = x.nrow();
  int ncol = x.ncol();
  IntegerMatrix out_mat(nrow, ncol);
  for (int i = 0; i < ncol; i++){

    IntegerMatrix::Column x_col = x(_, i);
    IntegerMatrix::Column y_col = y(_, i);
    IntegerMatrix::Column out_col = out_mat(_, i);
    out_col = sign(x_col - y_col);

  }
  return(out_mat);

}

// [[Rcpp::export]]
LogicalMatrix subset_matrix_cols_l(LogicalMatrix in_matrix, IntegerVector cols){

  int n_rows = in_matrix.nrow();
  int n_cols = cols.length();
  LogicalMatrix out_matrix(n_rows, n_cols);

  for (int i = 0; i < n_cols; i++){
    LogicalMatrix::Column original_col = in_matrix(_, cols[i]-1);
    LogicalMatrix::Column new_col = out_matrix(_, i);
    new_col = original_col;
  }
  return(out_matrix);
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
LogicalMatrix subset_matrix_l(LogicalMatrix in_matrix, IntegerVector rows, IntegerVector cols){

  int n_rows = rows.length();
  int n_cols = cols.length();
  LogicalMatrix out_matrix(n_rows, n_cols);

  for (int i = 0; i < n_cols; i++){

    LogicalMatrix::Column original_col = in_matrix(_, cols[i]-1);
    LogicalMatrix::Column new_col = out_matrix(_, i);

    for (int j = 0; j < n_rows; j++){

      new_col[j] = original_col[rows[j] - 1];

    }
  }
  return(out_matrix);
}


// [[Rcpp::export]]
Rcpp::NumericMatrix subset_matrix_rows(Rcpp::NumericMatrix in_matrix, Rcpp::IntegerVector rows){

  int n_cols = in_matrix.ncol();
  int n_rows = rows.length();
  Rcpp::NumericMatrix out_matrix(n_rows, n_cols);

  for (int i = 0; i < n_rows; i++){

    Rcpp::NumericMatrix::Row original_row = in_matrix(rows[i]-1, _);
    Rcpp::NumericMatrix::Row new_row = out_matrix(i, _);
    new_row = original_row;
  }
  return(out_matrix);
}

// [[Rcpp::export]]
IntegerMatrix subset_int_matrix_rows(IntegerMatrix in_matrix, IntegerVector rows){

  int n_cols = in_matrix.ncol();
  int n_rows = rows.length();
  IntegerMatrix out_matrix(n_rows, n_cols);

  for (int i = 0; i < n_rows; i++){

    IntegerMatrix::Row original_row = in_matrix(rows[i] - 1, _);
    IntegerMatrix::Row new_row = out_matrix(i, _);
    new_row = original_row;
  }
  return(out_matrix);
}


// [[Rcpp::export]]
NumericVector rowSds(IntegerMatrix x){

  int nrows = x.nrow();
  NumericVector out_vec(nrows);
  for (int i = 0; i < nrows; i++){

    IntegerVector this_row = x(i, _);
    double row_var = var(this_row);
    if (R_isnancpp(row_var)){
      row_var = 0;
    }
    out_vec[i] = sqrt(row_var);

  }
  return(out_vec);

}

// [[Rcpp::export]]
LogicalMatrix comp_mat_gt(IntegerMatrix in_mat, int comp_val){

  int nrows = in_mat.nrow();
  int ncols = in_mat.ncol();
  LogicalMatrix out_mat(nrows, ncols);
  for ( int i = 0; i < nrows; i++ ) {

    for ( int j = 0; j < ncols; j++ ) {

      out_mat(i, j) = in_mat(i, j) > comp_val;

    }

  }
  return out_mat;
}

// [[Rcpp::export]]
LogicalMatrix check_pos_risk(IntegerMatrix in_mat, NumericVector comp_vals){

  int nrows = in_mat.nrow();
  int ncols = in_mat.ncol();
  LogicalMatrix out_mat(nrows, ncols);
  for ( int i = 0; i < nrows; i++ ) {

    for ( int j = 0; j < ncols; j++ ) {

      out_mat(i, j) = in_mat(i, j) < comp_vals[j];

    }

  }
  return out_mat;
}


// [[Rcpp::export]]
LogicalMatrix check_neg_risk(IntegerMatrix in_mat, NumericVector comp_vals){

  int nrows = in_mat.nrow();
  int ncols = in_mat.ncol();
  LogicalMatrix out_mat(nrows, ncols);
  for ( int i = 0; i < nrows; i++ ) {

    for ( int j = 0; j < ncols; j++ ) {

      out_mat(i, j) = in_mat(i, j) > comp_vals[j];

    }

  }
  return out_mat;
}


// [[Rcpp::export]]
IntegerMatrix two_minus_mat(IntegerMatrix in_mat){

  int nrows = in_mat.nrow();
  int ncols = in_mat.ncol();
  IntegerMatrix out_mat(nrows, ncols);
  for ( int i = 0; i < ncols; i++ ) {

    IntegerMatrix::Column in_col = in_mat(_, i);
    IntegerMatrix::Column out_col = out_mat(_, i);
    out_col = 2 - in_col;

  }
  return out_mat;
}


// [[Rcpp::export]]
LogicalMatrix comp_mat_lt(IntegerMatrix in_mat, int comp_val){

  int nrows = in_mat.nrow();
  int ncols = in_mat.ncol();
  LogicalMatrix out_mat(nrows, ncols);
  for ( int i = 0; i < nrows; i++ ) {

    for ( int j = 0; j < ncols; j++ ) {

      out_mat(i, j) = in_mat(i, j) < comp_val;

    }

  }
  return out_mat;
}


// [[Rcpp::export]]
IntegerMatrix comp_mat_ne(IntegerMatrix in_mat, int comp_val){

  int nrows = in_mat.nrow();
  int ncols = in_mat.ncol();
  IntegerMatrix out_mat(nrows, ncols);
  for ( int i = 0; i < nrows; i++ ) {

    for ( int j = 0; j < ncols; j++ ) {

      out_mat(i, j) = in_mat(i, j) != comp_val;

    }

  }
  return out_mat;
}

// [[Rcpp::export]]
LogicalVector pos_risk_dir(arma::rowvec x){

  IntegerVector x_int = wrap(x);
  LogicalVector y = x_int > 0;
  return(y);

}


// [[Rcpp::export]]
LogicalVector neg_risk_dir(arma::rowvec x){

  IntegerVector x_int = wrap(x);
  LogicalVector y = x_int <= 0;
  return(y);

}

// [[Rcpp::export]]
arma::mat mult_rows_by_scalars(IntegerMatrix in_mat, NumericVector scalars){

  //convert to arma versions
  arma::mat arma_in_mat = as<arma::mat>(in_mat);
  arma::vec arma_scalars = as<arma::vec>(scalars);
  arma_in_mat.each_col() %= arma_scalars;
  return(arma_in_mat);
}

// [[Rcpp::export]]
arma::mat mult_rows_by_scalars_arma(arma::mat arma_in_mat, IntegerVector scalars){

  //convert to arma versions
  arma::vec arma_scalars = as<arma::vec>(scalars);
  arma_in_mat.each_col() %= arma_scalars;
  return(arma_in_mat);

}

// [[Rcpp::export]]
List find_high_risk(int n_target, int n_pos, int n_neg, LogicalVector neg_risk, LogicalVector pos_risk,
                    IntegerMatrix case_inf, IntegerMatrix comp_inf, IntegerVector idx_vec){

  if (n_neg > 0 & n_pos > 0) {

    IntegerVector neg_risk_int = idx_vec[neg_risk];
    IntegerMatrix case_inf_neg = subset_int_matrix_rows(case_inf, neg_risk_int);
    IntegerMatrix comp_inf_neg = subset_int_matrix_rows(comp_inf, neg_risk_int);
    LogicalMatrix n1 = comp_mat_lt(case_inf_neg, 2);
    LogicalMatrix n2 = comp_mat_lt(comp_inf_neg, 2);

    IntegerVector pos_risk_int = idx_vec[pos_risk];
    IntegerMatrix case_inf_pos = subset_int_matrix_rows(case_inf, pos_risk_int);
    IntegerMatrix comp_inf_pos = subset_int_matrix_rows(comp_inf, pos_risk_int);
    LogicalMatrix p1 = comp_mat_gt(case_inf_pos, 0);
    LogicalMatrix p2 = comp_mat_gt(comp_inf_pos, 0);

    IntegerVector n1_colsums = colSums(n1);
    IntegerVector n2_colsums = colSums(n2);
    IntegerVector p1_colsums = colSums(p1);
    IntegerVector p2_colsums = colSums(p2);
    LogicalVector case_high_risk = (p1_colsums + n1_colsums) == n_target;
    LogicalVector comp_high_risk = (p2_colsums + n2_colsums) == n_target;

    List res = List::create(Named("case_high_risk") = case_high_risk,
                            Named("comp_high_risk") = comp_high_risk,
                            Named("pos_risk_int") = pos_risk_int,
                            Named("neg_risk_int") = neg_risk_int);
    return(res);

  } else if (n_neg > 0 & n_pos == 0) {

    IntegerVector neg_risk_int = idx_vec[neg_risk];
    IntegerMatrix case_inf_neg = subset_int_matrix_rows(case_inf, neg_risk_int);
    IntegerMatrix comp_inf_neg = subset_int_matrix_rows(comp_inf, neg_risk_int);
    LogicalMatrix n1 = comp_mat_lt(case_inf_neg, 2);
    LogicalMatrix n2 = comp_mat_lt(comp_inf_neg, 2);

    IntegerVector n1_colsums = colSums(n1);
    IntegerVector n2_colsums = colSums(n2);
    LogicalVector case_high_risk = n1_colsums == n_target;
    LogicalVector comp_high_risk = n2_colsums == n_target;

    IntegerVector pos_risk_int = IntegerVector::create(NA_INTEGER);

    List res = List::create(Named("case_high_risk") = case_high_risk,
                            Named("comp_high_risk") = comp_high_risk,
                            Named("pos_risk_int") = pos_risk_int,
                            Named("neg_risk_int") = neg_risk_int);
    return(res);

  } else {

    IntegerVector pos_risk_int = idx_vec[pos_risk];
    IntegerMatrix case_inf_pos = subset_int_matrix_rows(case_inf, pos_risk_int);
    IntegerMatrix comp_inf_pos = subset_int_matrix_rows(comp_inf, pos_risk_int);
    LogicalMatrix p1 = comp_mat_gt(case_inf_pos, 0);
    LogicalMatrix p2 = comp_mat_gt(comp_inf_pos, 0);

    IntegerVector p1_colsums = colSums(p1);
    IntegerVector p2_colsums = colSums(p2);
    LogicalVector case_high_risk = p1_colsums == n_target;
    LogicalVector comp_high_risk = p2_colsums == n_target;

    IntegerVector neg_risk_int = IntegerVector::create(NA_INTEGER);

    List res = List::create(Named("case_high_risk") = case_high_risk,
                            Named("comp_high_risk") = comp_high_risk,
                            Named("pos_risk_int") = pos_risk_int,
                            Named("neg_risk_int") = neg_risk_int);
    return(res);

  }

}

// [[Rcpp::export]]
List compute_dif_vecs(IntegerMatrix case_genetic_data, IntegerMatrix comp_genetic_data, IntegerMatrix case_comp_dif,
                      IntegerMatrix cases_minus_complements, IntegerMatrix both_one_mat,
                      NumericVector weight_lookup, int n_different_snps_weight = 2, int n_both_one_weight = 1,
                      IntegerVector prev_informative_families = IntegerVector::create(NA_INTEGER)){

  // determine whether families are informative for the set of target_snps
  IntegerVector total_different_snps = colSums(case_comp_dif);
  LogicalVector informative_families_l = total_different_snps != 0;
  IntegerVector family_idx = seq_len(total_different_snps.length());
  IntegerVector informative_families = family_idx[informative_families_l];
  IntegerMatrix case_inf = subset_matrix_cols(case_genetic_data, informative_families);
  IntegerMatrix comp_inf = subset_matrix_cols(comp_genetic_data, informative_families);
  IntegerMatrix cases_minus_complements_inf = subset_matrix_cols(cases_minus_complements, informative_families);

  // compute weights
  IntegerMatrix both_one_inf = subset_matrix_cols(both_one_mat, informative_families);
  IntegerVector both_one = colSums(both_one_inf);
  IntegerVector total_different_inf = total_different_snps[informative_families_l];
  IntegerVector weighted_informativeness = n_both_one_weight * both_one + n_different_snps_weight * total_different_inf;
  NumericVector family_weights(weighted_informativeness.length());
  for (int i = 0; i < weighted_informativeness.length(); i++){

    int weight_i = weighted_informativeness[i];
    family_weights[i] = weight_lookup[weight_i - 1];

  }
  double invsum_family_weights = 1 / sum(family_weights);

  // compute weighted difference vectors for cases vs complements
  NumericVector std_family_weights = family_weights * invsum_family_weights;
  arma::mat dif_vecs = mult_rows_by_scalars(transpose(cases_minus_complements_inf), std_family_weights);

  // take the sum of the case - complement difference vectors over families
  arma::rowvec sum_dif_vecs = arma::sum(dif_vecs, 0);

  // determine how many cases and complements actually have the proposed risk set
  arma::rowvec risk_dirs = sign(sum_dif_vecs);
  LogicalVector pos_risk = pos_risk_dir(risk_dirs);
  int n_pos = sum(pos_risk);
  LogicalVector neg_risk = neg_risk_dir(risk_dirs);
  int n_neg = sum(neg_risk);
  int n_target = case_comp_dif.nrow();
  IntegerVector idx_vec = seq_len(n_target);

  //using helper function to pick out high risk families
  List high_risk = find_high_risk(n_target, n_pos, n_neg, neg_risk, pos_risk, case_inf,
                                  comp_inf, idx_vec);
  LogicalVector case_high_risk = high_risk["case_high_risk"];
  LogicalVector comp_high_risk = high_risk["comp_high_risk"];
  IntegerVector pos_risk_int = high_risk["pos_risk_int"];
  IntegerVector neg_risk_int = high_risk["neg_risk_int"];

  int n_case_high_risk = sum(case_high_risk);
  int n_comp_high_risk = sum(comp_high_risk);

  IntegerVector inf_family_idx = seq_len(case_high_risk.length());
  IntegerVector case_high_risk_int = inf_family_idx[case_high_risk];
  IntegerVector comp_high_risk_int = inf_family_idx[comp_high_risk];
  IntegerMatrix case_high_inf = subset_matrix_cols(case_inf, case_high_risk_int);
  IntegerMatrix comp_high_inf = subset_matrix_cols(comp_inf, comp_high_risk_int);

  // for epi test
  bool prev_inf = all(is_na(prev_informative_families));
  if (!prev_inf){

    informative_families = prev_informative_families[informative_families - 1];
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
                            Named("informative_families") = informative_families);
    return(res);

  } else {

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
                            Named("informative_families") = informative_families);
    return(res);

  }

}


// [[Rcpp::export]]
List chrom_fitness_score(IntegerMatrix case_genetic_data, IntegerMatrix complement_genetic_data, IntegerMatrix case_comp_differences,
                         IntegerVector target_snps, IntegerMatrix cases_minus_complements, IntegerMatrix both_one_mat,
                         LogicalMatrix block_ld_mat, NumericVector weight_lookup,
                         int n_different_snps_weight = 2, int n_both_one_weight = 1,
                         int n_case_high_risk_thresh = 20, double outlier_sd = 2.5, bool epi_test = false) {

  // pick out the differences for the target snps

  IntegerMatrix case_comp_dif = subset_matrix_cols(case_comp_differences, target_snps);
  case_comp_dif = transpose(case_comp_dif);
  int n_target = target_snps.length();

  // also subset to target cols for the rest of the inputs
  case_genetic_data = transpose(subset_matrix_cols(case_genetic_data, target_snps));
  complement_genetic_data = transpose(subset_matrix_cols(complement_genetic_data, target_snps));
  cases_minus_complements = transpose(subset_matrix_cols(cases_minus_complements, target_snps));
  both_one_mat = transpose(subset_matrix_cols(both_one_mat, target_snps));

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
  int n_informative_families = informative_families.length();
  IntegerMatrix cases_minus_complements_inf = dif_vec_list["cases_minus_complements_inf"];

  // pick out misclassifications via outlier detection, indicating recessive mode of inheritance

  // initialize vector of pattern of inheritance
  CharacterVector risk_set_alleles(target_snps.length(), "1+");

  // only applies if we have at least 20 high risk case
  if (n_case_high_risk > n_case_high_risk_thresh) {

    NumericVector case_high_risk_means = rowMeans(case_high_inf);
    NumericVector case_high_risk_sd = rowSds(case_high_inf);
    NumericVector high_outlier_thresh = case_high_risk_means + outlier_sd * case_high_risk_sd;
    NumericVector low_outlier_thresh = case_high_risk_means - outlier_sd * case_high_risk_sd;
    IntegerMatrix all_high_risk = transpose(cbind(case_high_inf, comp_high_inf));
    LogicalVector outliers(n_target, false);

    if (n_pos > 0){

      IntegerMatrix positive_high_risk = subset_matrix_cols(all_high_risk, pos_risk_int);
      LogicalMatrix positive_outlier_mat = check_pos_risk(positive_high_risk, low_outlier_thresh[pos_risk]);
      LogicalVector positive_outliers = colSums(positive_outlier_mat) > 0;
      outliers[pos_risk_int - 1] = positive_outliers;

    }

    if (n_neg > 0){

      IntegerMatrix negative_high_risk = subset_matrix_cols(all_high_risk, neg_risk_int);
      LogicalMatrix negative_outlier_mat = check_neg_risk(negative_high_risk, high_outlier_thresh[neg_risk]);
      LogicalVector negative_outliers = colSums(negative_outlier_mat) > 0;
      outliers[neg_risk_int - 1] = negative_outliers;

    }

    // if there are outliers, recompute the weights and associated statistics
    int n_outliers = sum(outliers);
    if (n_outliers > 0) {

      risk_set_alleles[outliers] = "2";
      int n_pos_outliers = 0;
      IntegerVector pos_outlier_cols = IntegerVector::create(NA_INTEGER);
      if (sum(pos_risk) > 0){

        LogicalVector pos_risk_outliers = outliers[pos_risk_int - 1];
        pos_outlier_cols = pos_risk_int[pos_risk_outliers];
        n_pos_outliers = sum(pos_risk_outliers);

      }

      int n_neg_outliers = 0;
      IntegerVector neg_outlier_cols = IntegerVector::create(NA_INTEGER);
      if(sum(neg_risk) > 0){

        LogicalVector neg_risk_outliers = outliers[neg_risk_int - 1];
        neg_outlier_cols = neg_risk_int[neg_risk_outliers];
        n_neg_outliers = sum(neg_risk_outliers);

      }

      // recode instances where the model appears to be recessive
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
          IntegerMatrix::Row both_one_target_row = both_one_inf(outlier_idx, _);
          IntegerVector both_one_new_vals(n_informative_families, 0);
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
          IntegerMatrix::Row both_one_target_row = both_one_inf(outlier_idx, _);
          IntegerVector both_one_new_vals(n_informative_families, 0);
          both_one_target_row = both_one_new_vals;

        }

      }

      //recompute the number of informative families
      cases_minus_complements = sign_subtract_mat(case_inf, comp_inf);
      case_comp_dif = comp_mat_ne(cases_minus_complements, 0);
      dif_vec_list = compute_dif_vecs(case_inf, comp_inf, case_comp_dif, cases_minus_complements, both_one_inf,
                                      weight_lookup, n_different_snps_weight, n_both_one_weight, informative_families);

      //pick out required pieces
      arma::rowvec sum_dif_vecs_tmp = dif_vec_list["sum_dif_vecs"];
      sum_dif_vecs = sum_dif_vecs_tmp;
      family_weights = dif_vec_list["family_weights"];
      invsum_family_weights = dif_vec_list["invsum_family_weights"];
      informative_families = dif_vec_list["informative_families"];
      n_case_high_risk = dif_vec_list["n_case_high_risk"];
      n_comp_high_risk = dif_vec_list["n_comp_high_risk"];
      IntegerMatrix case_high_inf_tmp = dif_vec_list["case_high_inf"];
      case_high_inf = case_high_inf_tmp;
      IntegerMatrix comp_high_inf_tmp = dif_vec_list["comp_high_inf"];
      comp_high_inf = comp_high_inf_tmp;
      n_pos = dif_vec_list["n_pos"];
      n_neg = dif_vec_list["n_neg"];
      pos_risk_int = dif_vec_list["pos_risk_int"];
      pos_risk = dif_vec_list["pos_risk"];
      neg_risk_int = dif_vec_list["neg_risk_int"];
      neg_risk = dif_vec_list["neg_risk"];
      IntegerMatrix case_inf_tmp = dif_vec_list["case_inf"];
      case_inf = case_inf_tmp;
      IntegerMatrix comp_inf_tmp = dif_vec_list["comp_inf"];
      comp_inf = comp_inf_tmp;
      IntegerMatrix both_one_inf_tmp = dif_vec_list["both_one_inf"];
      both_one_inf = both_one_inf_tmp;
      n_informative_families = informative_families.length();
      IntegerMatrix cases_minus_complements_inf_tmp = dif_vec_list["cases_minus_complements_inf"];
      cases_minus_complements_inf = cases_minus_complements_inf_tmp;

    }
  }

  // count the number of risk alleles in those with the full risk set

  double total_case_high_risk_alleles = 0;
  double total_comp_high_risk_alleles = 0;
  if (sum(pos_risk) > 0){

    IntegerMatrix case_high_inf_pos = subset_int_matrix_rows(case_high_inf, pos_risk_int);
    IntegerVector case_pos_alleles = colSums(case_high_inf_pos);
    total_case_high_risk_alleles += sum(case_pos_alleles);

    IntegerMatrix comp_high_inf_pos = subset_int_matrix_rows(comp_high_inf, pos_risk_int);
    IntegerVector comp_pos_alleles = colSums(comp_high_inf_pos);
    total_comp_high_risk_alleles += sum(comp_pos_alleles);

  }
  if (sum(neg_risk) > 0){

    IntegerMatrix case_low_inf_pos = two_minus_mat(case_high_inf);
    case_low_inf_pos = subset_int_matrix_rows(case_low_inf_pos, neg_risk_int);
    IntegerVector case_neg_alleles = colSums(case_low_inf_pos);
    total_case_high_risk_alleles += sum(case_neg_alleles);

    IntegerMatrix comp_low_inf_pos = two_minus_mat(comp_high_inf);
    comp_low_inf_pos =  subset_int_matrix_rows(comp_low_inf_pos, neg_risk_int);
    IntegerVector comp_neg_alleles = colSums(comp_low_inf_pos);
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
  arma::rowvec mu_hat = sum_dif_vecs;
  arma::mat mu_hat_mat(n_informative_families, n_target);
  for (int i = 0; i < n_informative_families; i++){
    mu_hat_mat.row(i) = mu_hat;
  }
  arma::mat x = as<arma::mat>(transpose(cases_minus_complements_inf));
  arma::mat x_minus_mu_hat = x - mu_hat_mat;
  arma::mat weighted_x_minus_mu_hat = mult_rows_by_scalars_arma(x_minus_mu_hat, family_weights);
  arma::mat cov_mat = invsum_family_weights * trans(weighted_x_minus_mu_hat) * x_minus_mu_hat;

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



