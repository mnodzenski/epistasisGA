#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

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

IntegerMatrix subtract_mat(IntegerMatrix x, IntegerMatrix y){

  int nrow = x.nrow();
  int ncol = x.ncol();
  IntegerMatrix out_mat(nrow, ncol);
  for (int i = 0; i < ncol; i++){

    IntegerMatrix::Column x_col = x(_, i);
    IntegerMatrix::Column y_col = y(_, i);
    IntegerMatrix::Column out_col = out_mat(_, i);
    out_col = x_col - y_col;

  }
  return(out_mat);

}

NumericMatrix subtract_mat_n(NumericMatrix x, NumericMatrix y){

  int nrow = x.nrow();
  int ncol = x.ncol();
  NumericMatrix out_mat(nrow, ncol);
  for (int i = 0; i < ncol; i++){

    NumericMatrix::Column x_col = x(_, i);
    NumericMatrix::Column y_col = y(_, i);
    NumericMatrix::Column out_col = out_mat(_, i);
    out_col = x_col - y_col;

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
  for ( int i = 0; i < ncols; i++ ) {

    IntegerMatrix::Column original_col = in_mat(_, i);
    LogicalMatrix::Column new_col = out_mat(_, i);
    new_col = original_col > comp_val;

  }
  return out_mat;
}

// [[Rcpp::export]]
LogicalMatrix comp_mat_lt(IntegerMatrix in_mat, int comp_val){

  int nrows = in_mat.nrow();
  int ncols = in_mat.ncol();
  LogicalMatrix out_mat(nrows, ncols);
  for ( int i = 0; i < ncols; i++ ) {

    IntegerMatrix::Column original_col = in_mat(_, i);
    LogicalMatrix::Column new_col = out_mat(_, i);
    new_col = original_col <= comp_val;

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
IntegerMatrix abs_matrix(IntegerMatrix in_mat){

  int nrows = in_mat.nrow();
  int ncols = in_mat.ncol();
  IntegerMatrix out_mat(nrows, ncols);
  for ( int i = 0; i < ncols; i++ ) {

    IntegerMatrix::Column original_col = in_mat(_, i);
    IntegerMatrix::Column new_col = out_mat(_, i);
    new_col = abs(original_col);

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
IntegerMatrix mult_cols_by_scalars(IntegerMatrix in_matrix, NumericVector scalars){

  int n_rows = in_matrix.nrow();
  int n_cols = in_matrix.ncol();
  IntegerMatrix out_matrix(n_rows, n_cols);
  //NumericMatrix in_matrix_n = as<NumericMatrix>(in_matrix);

  for (int i = 0; i < n_cols; i++){
    IntegerMatrix::Column original_col = in_matrix(_, i);
    double scale_val = scalars[i];
    IntegerMatrix::Column new_col = out_matrix(_, i);
    new_col = scale_val*original_col;
  }
  return(out_matrix);

}

// [[Rcpp::export]]
NumericMatrix mult_rows_by_scalars(NumericMatrix in_matrix, NumericVector scalars){

  int n_rows = in_matrix.nrow();
  int n_cols = in_matrix.ncol();
  NumericMatrix out_matrix(n_rows, n_cols);

  for (int i = 0; i < n_rows; i++){
    NumericMatrix::Row original_row = in_matrix(i, _);
    double scale_val = scalars[i];
    NumericMatrix::Row new_row = out_matrix(i, _);
    new_row = scale_val*original_row;
  }
  return(out_matrix);

}

// [[Rcpp::export]]
IntegerVector subset_colsums(LogicalMatrix in_mat, IntegerVector row_idx){

  int n_cols = in_mat.ncol();
  IntegerVector out_vec(n_cols, 0);
  for (int i = 0; i < n_cols; i++){
    for (int j = 0; j < row_idx.length(); j++){

      int this_row = row_idx[j] - 1;
      out_vec[i] += in_mat(this_row, i);

    }
  }
  return(out_vec);
}

// [[Rcpp::export]]
IntegerVector subset_colsums_int(IntegerMatrix in_mat, IntegerVector row_idx){

  int n_cols = in_mat.ncol();
  IntegerVector out_vec(n_cols, 0);
  for (int i = 0; i < n_cols; i++){
    for (int j = 0; j < row_idx.length(); j++){

      int this_row = row_idx[j] - 1;
      out_vec[i] += in_mat(this_row, i);

    }
  }
  return(out_vec);
}


// [[Rcpp::export]]
IntegerVector subset_rowsums(IntegerMatrix in_mat, IntegerVector row_idx){

  int nrows = row_idx.length();
  IntegerVector out_vec(nrows);
  for (int i = 0; i < nrows; i++){

    int this_row = row_idx[i] - 1;
    IntegerMatrix::Row target = in_mat(this_row, _);
    out_vec[i] = sum(target);

  }
  return(out_vec);

}



// [[Rcpp::export]]
List find_high_risk(int n_target, int n_pos, int n_neg, IntegerVector neg_risk_int, IntegerVector pos_risk_int,
                    IntegerMatrix case_data, IntegerMatrix comp, LogicalVector uninformative_families_l){

  if (n_neg > 0 & n_pos > 0) {

    LogicalMatrix n1 = comp_mat_lt(case_data, 2);
    LogicalMatrix n2 = comp_mat_lt(comp, 2);

    LogicalMatrix p1 = comp_mat_gt(case_data, 0);
    LogicalMatrix p2 = comp_mat_gt(comp, 0);

    IntegerVector n1_colsums = subset_colsums(n1, neg_risk_int);
    IntegerVector n2_colsums = subset_colsums(n2, neg_risk_int);

    IntegerVector p1_colsums = subset_colsums(p1, pos_risk_int);
    IntegerVector p2_colsums = subset_colsums(p2, pos_risk_int);

    LogicalVector case_high_risk = (p1_colsums + n1_colsums) == n_target;
    case_high_risk[uninformative_families_l] = false;

    LogicalVector comp_high_risk = (p2_colsums + n2_colsums) == n_target;
    comp_high_risk[uninformative_families_l] = false;

    List res = List::create(Named("case_high_risk") = case_high_risk,
                            Named("comp_high_risk") = comp_high_risk);
    return(res);

  } else if (n_neg > 0 & n_pos == 0) {

    LogicalMatrix n1 = comp_mat_lt(case_data, 2);
    LogicalMatrix n2 = comp_mat_lt(comp, 2);

    IntegerVector n1_colsums = subset_colsums(n1, neg_risk_int);
    IntegerVector n2_colsums = subset_colsums(n2, neg_risk_int);

    LogicalVector case_high_risk = n1_colsums == n_target;
    case_high_risk[uninformative_families_l] = false;

    LogicalVector comp_high_risk = n2_colsums == n_target;
    comp_high_risk[uninformative_families_l] = false;

    List res = List::create(Named("case_high_risk") = case_high_risk,
                            Named("comp_high_risk") = comp_high_risk);
    return(res);

  } else {

    LogicalMatrix p1 = comp_mat_gt(case_data, 0);
    LogicalMatrix p2 = comp_mat_gt(comp, 0);

    IntegerVector p1_colsums = subset_colsums(p1, pos_risk_int);
    IntegerVector p2_colsums = subset_colsums(p2, pos_risk_int);

    LogicalVector case_high_risk = p1_colsums == n_target;
    case_high_risk[uninformative_families_l] = false;

    LogicalVector comp_high_risk = p2_colsums == n_target;
    comp_high_risk[uninformative_families_l] = false;

    List res = List::create(Named("case_high_risk") = case_high_risk,
                            Named("comp_high_risk") = comp_high_risk);
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

  // compute weights
  IntegerVector both_one = colSums(both_one_mat);
  IntegerVector weighted_informativeness = n_both_one_weight * both_one + n_different_snps_weight * total_different_snps;
  //NumericVector family_weights = weight_lookup[weighted_informativeness - 1];
  NumericVector family_weights(weighted_informativeness.length());
  for (int i = 0; i < weighted_informativeness.length(); i++){

    int weight_i = weighted_informativeness[i];
    family_weights[i] = weight_lookup[weight_i - 1];

  }
  family_weights[uninformative_families_l] = 0;
  double invsum_family_weights = 1 / sum(family_weights);

  // compute weighted difference vectors for cases vs complements
  //NumericVector std_family_weights = family_weights * invsum_family_weights;
  IntegerMatrix dif_vecs = mult_cols_by_scalars(cases_minus_complements, family_weights);

  // take the sum of the case - complement difference vectors over families
  NumericVector sum_dif_vecs = as<NumericVector>(rowSums(dif_vecs));
  for (int i = 0; i < sum_dif_vecs.length(); i++){

    sum_dif_vecs[i] *= invsum_family_weights;

  }

  // determine how many cases and complements actually have the proposed risk set
  IntegerVector risk_dirs = sign(sum_dif_vecs);

  LogicalVector pos_risk = risk_dirs > 0;
  int n_pos = sum(pos_risk);

  LogicalVector neg_risk = risk_dirs <= 0;
  int n_neg = sum(neg_risk);

  int n_target = case_comp_dif.nrow();
  IntegerVector idx_vec = seq_len(n_target);

  IntegerVector pos_risk_int = idx_vec[pos_risk];
  IntegerVector neg_risk_int = idx_vec[neg_risk];

  //using helper function to pick out high risk families
  List high_risk = find_high_risk(n_target, n_pos, n_neg, neg_risk_int, pos_risk_int, case_genetic_data,
                                  comp_genetic_data, uninformative_families_l);
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
  NumericVector sum_dif_vecs = dif_vec_list["sum_dif_vecs"];
  IntegerVector family_weights = dif_vec_list["family_weights"];
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

  IntegerVector tmp_idx = seq_len(case_genetic_data.ncol());
  IntegerVector case_high_risk_int = tmp_idx[case_high_risk];
  IntegerVector comp_high_risk_int = tmp_idx[comp_high_risk];

  // pick out misclassifications via outlier detection, indicating recessive mode of inheritance

  // initialize vector of pattern of inheritance
  CharacterVector risk_set_alleles(target_snps.length(), "1+");

  // only applies if we have at least 20 high risk case
  if (n_case_high_risk > n_case_high_risk_thresh) {

    IntegerMatrix case_high_inf = subset_matrix_cols(case_genetic_data, case_high_risk_int);
    IntegerMatrix comp_high_inf = subset_matrix_cols(complement_genetic_data, comp_high_risk_int);

    NumericVector case_high_risk_means = rowMeans(case_high_inf);
    NumericVector case_high_risk_sd = rowSds(case_high_inf);
    NumericVector high_outlier_thresh = case_high_risk_means + outlier_sd * case_high_risk_sd;
    NumericVector low_outlier_thresh = case_high_risk_means - outlier_sd * case_high_risk_sd;

    int outlier_count = 0;
    LogicalVector outliers(n_target, false);

    if (n_pos > 0){

      IntegerVector p_tmp = seq_len(n_target);
      IntegerVector pos_risk_int = p_tmp[pos_risk];

      for (int i = 0; i < n_pos; i++){

        int this_row = pos_risk_int[i] - 1;
        double thresh_val = low_outlier_thresh[this_row];

        if (thresh_val > 1){

          // check cases for outliers
          IntegerMatrix::Row case_row = case_high_inf(this_row, _);
          int case_outliers = sum(case_row < thresh_val);

          //check comps for outliers
          IntegerMatrix::Row comp_row = comp_high_inf(this_row, _);
          int comp_outliers = sum(comp_row < thresh_val);

          // if there are outliers, recode
          if ((case_outliers + comp_outliers) > 0){

            // increment outlier count
            outlier_count += 1;

            // indicate we need two risk alleles
            risk_set_alleles[this_row] = "2";

            // recode cases
            IntegerMatrix::Row case_target_row = case_genetic_data(this_row, _);
            IntegerVector case_new_vals = case_target_row;
            case_new_vals[case_new_vals == 1] = 0;
            case_target_row = case_new_vals;

            //then complements
            IntegerMatrix::Row comp_target_row = complement_genetic_data(this_row, _);
            IntegerVector comp_new_vals = comp_target_row;
            comp_new_vals[comp_new_vals == 1] = 0;
            comp_target_row = comp_new_vals;

            //then the both one matrix
            IntegerMatrix::Row both_one_target_row = both_one_mat(this_row, _);
            IntegerVector both_one_new_vals = both_one_target_row;
            both_one_new_vals[both_one_new_vals == 1] = 0;
            both_one_target_row = both_one_new_vals;

          }

        }

      }

    }

    if (n_neg > 0){

      IntegerVector n_tmp = seq_len(n_target);
      IntegerVector neg_risk_int = n_tmp[neg_risk];

      for (int i = 0; i < n_neg; i++){

        int this_row = neg_risk_int[i] - 1;
        double thresh_val = high_outlier_thresh[this_row];

        if (thresh_val < 1){

          // check cases for outliers
          IntegerMatrix::Row case_row = case_high_inf(this_row, _);
          int case_outliers = sum(case_row > thresh_val);

          //check comps for outliers
          IntegerMatrix::Row comp_row = comp_high_inf(this_row, _);
          int comp_outliers = sum(comp_row > thresh_val);

          // if there are outliers, recode
          if ((case_outliers + comp_outliers) > 0){

            // increment outlier count
            outlier_count += 1;

            // indicate we need two risk alleles
            risk_set_alleles[this_row] = "2";

            // recode cases
            IntegerMatrix::Row case_target_row = case_genetic_data(this_row, _);
            IntegerVector case_new_vals = case_target_row;
            case_new_vals[case_new_vals == 1] = 2;
            case_target_row = case_new_vals;

            //then complements
            IntegerMatrix::Row comp_target_row = complement_genetic_data(this_row, _);
            IntegerVector comp_new_vals = comp_target_row;
            comp_new_vals[comp_new_vals == 1] = 2;
            comp_target_row = comp_new_vals;

            //then the both one matrix
            IntegerMatrix::Row both_one_target_row = both_one_mat(this_row, _);
            IntegerVector both_one_new_vals = both_one_target_row;
            both_one_new_vals[both_one_new_vals == 1] = 0;
            both_one_target_row = both_one_new_vals;

          }

        }

      }

    }

    // if there are outliers, recompute the weights and associated statistics
    if (outlier_count > 0) {

      //recompute the number of informative families
      cases_minus_complements = sign_subtract_mat(case_genetic_data, complement_genetic_data);
      case_comp_dif = abs_matrix(cases_minus_complements);
      dif_vec_list = compute_dif_vecs(case_genetic_data, complement_genetic_data, case_comp_dif, cases_minus_complements, both_one_mat,
                                      weight_lookup, n_different_snps_weight, n_both_one_weight, informative_families);

      //pick out required pieces
      // the tmp var versions are used because direct reassignment
      // causes compile errors
      NumericVector sum_dif_vecs_tmp = dif_vec_list["sum_dif_vecs"];
      sum_dif_vecs = sum_dif_vecs_tmp;
      IntegerVector family_weights_tmp = dif_vec_list["family_weights"];
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

      IntegerVector tmp_idx = seq_len(case_genetic_data.ncol());
      IntegerVector case_high_risk_int = tmp_idx[case_high_risk];
      IntegerVector comp_high_risk_int = tmp_idx[comp_high_risk];

    }
  }

  // count the number of risk alleles in those with the full risk set

  double total_case_high_risk_alleles = 0;
  double total_comp_high_risk_alleles = 0;
  IntegerMatrix case_high_inf = subset_matrix_cols(case_genetic_data, case_high_risk_int);
  IntegerMatrix comp_high_inf = subset_matrix_cols(complement_genetic_data, comp_high_risk_int);

  if (n_pos > 0){

    IntegerVector p_tmp = seq_len(n_target);
    IntegerVector pos_risk_int = p_tmp[pos_risk];

    // high risk case alleles
    IntegerVector case_pos_alleles = subset_colsums_int(case_high_inf, pos_risk_int);
    total_case_high_risk_alleles += sum(case_pos_alleles);

    // high risk comp alleles
    IntegerVector comp_pos_alleles = subset_colsums_int(comp_high_inf, pos_risk_int);
    total_comp_high_risk_alleles += sum(comp_pos_alleles);

  }
  if (n_neg > 0){

    IntegerVector n_tmp = seq_len(n_target);
    IntegerVector neg_risk_int = n_tmp[neg_risk];

    // high risk case alleles
    IntegerVector case_neg_alleles = subset_colsums_int(case_high_inf, neg_risk_int);
    total_case_high_risk_alleles += sum(case_neg_alleles);

    // high risk comp alleles
    IntegerVector comp_neg_alleles = subset_colsums_int(comp_high_inf, neg_risk_int);
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
  arma::mat x = as<arma::mat>(transpose(subset_matrix_cols(cases_minus_complements, informative_families)));
  arma::vec inf_family_weights = as<arma::vec>(family_weights[family_weights > 0]);
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




