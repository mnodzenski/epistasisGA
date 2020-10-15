#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
arma::mat mult_rows_by_scalars_arma(arma::mat arma_in_mat, arma::vec arma_scalars){

  //convert to arma versions
  arma_in_mat.each_col() %= arma_scalars;
  return(arma_in_mat);

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
arma::mat mult_rows_by_scalars_orig(IntegerMatrix in_mat, NumericVector scalars){

  //convert to arma versions
  arma::mat arma_in_mat = as<arma::mat>(in_mat);
  arma::vec arma_scalars = as<arma::vec>(scalars);
  arma_in_mat.each_col() %= arma_scalars;
  return(arma_in_mat);
}



arma::mat mult_rows_by_scalars(IntegerMatrix in_mat, NumericVector scalars){

  //convert to arma versions
  arma::mat arma_in_mat = as<arma::mat>(in_mat);
  arma::vec arma_scalars = as<arma::vec>(scalars);
  arma_in_mat.each_row() %= arma_scalars;
  return(arma_in_mat);
}

// [[Rcpp::export]]
List find_high_risk(int n_target, int n_pos, int n_neg, arma::uvec neg_risk, arma::uvec pos_risk,
                    arma::mat case_data, arma::mat comp){

  if (n_neg > 0 & n_pos > 0) {

    arma::vec n1_colsums(case_data.n_cols, arma::fill::zeros);
    for (int i = 0; i < neg_risk.n_elem; i++){

      int this_row = neg_risk[i];

      for (int j = 0; j < case_data.n_cols; j++){

        if (case_data(this_row, j) < 2) n1_colsums(j) += 1;

      }

    }

    arma::vec n2_colsums(comp.n_cols, arma::fill::zeros);
    for (int i = 0; i < neg_risk.n_elem; i++){

      int this_row = neg_risk[i];

      for (int j = 0; j < comp.n_cols; j++){

        if (comp(this_row, j) < 2) n2_colsums(j) += 1;

      }

    }

    arma::vec p1_colsums(case_data.n_cols, arma::fill::zeros);
    for (int i = 0; i < pos_risk.n_elem; i++){

      int this_row = pos_risk[i];

      for (int j = 0; j < case_data.n_cols; j++){

        if (case_data(this_row, j) > 0) p1_colsums(j) += 1;

      }

    }

    arma::vec p2_colsums(comp.n_cols, arma::fill::zeros);
    for (int i = 0; i < pos_risk.n_elem; i++){

      int this_row = pos_risk[i];

      for (int j = 0; j < comp.n_cols; j++){

        if (comp(this_row, j) > 0) p2_colsums(j) += 1;

      }

    }

    arma::vec case_colsums = p1_colsums + n1_colsums;
    arma::uvec case_high_risk = find(case_colsums == n_target);

    arma::vec comp_colsums = p2_colsums + n2_colsums;
    arma::uvec comp_high_risk = find(comp_colsums == n_target);

    List res = List::create(Named("case_high_risk") = case_high_risk,
                            Named("comp_high_risk") = comp_high_risk);
    return(res);

  } else if (n_neg > 0 & n_pos == 0) {

    arma::vec n1_colsums(case_data.n_cols, arma::fill::zeros);
    for (int i = 0; i < neg_risk.n_elem; i++){

      int this_row = neg_risk[i];

      for (int j = 0; j < case_data.n_cols; j++){

        if (case_data(this_row, j) < 2) n1_colsums(j) += 1;

      }

    }

    arma::vec n2_colsums(comp.n_cols, arma::fill::zeros);
    for (int i = 0; i < neg_risk.n_elem; i++){

      int this_row = neg_risk[i];

      for (int j = 0; j < comp.n_cols; j++){

        if (comp(this_row, j) < 2) n2_colsums(j) += 1;

      }

    }

    arma::uvec case_high_risk = find(n1_colsums == n_target);
    arma::uvec comp_high_risk = find(n2_colsums == n_target);

    List res = List::create(Named("case_high_risk") = case_high_risk,
                            Named("comp_high_risk") = comp_high_risk);
    return(res);

  } else {

    arma::vec p1_colsums(case_data.n_cols, arma::fill::zeros);
    for (int i = 0; i < pos_risk.n_elem; i++){

      int this_row = pos_risk[i];

      for (int j = 0; j < case_data.n_cols; j++){

        if (case_data(this_row, j) > 0) p1_colsums(j) += 1;

      }

    }

    arma::vec p2_colsums(comp.n_cols, arma::fill::zeros);
    for (int i = 0; i < pos_risk.n_elem; i++){

      int this_row = pos_risk[i];

      for (int j = 0; j < comp.n_cols; j++){

        if (comp(this_row, j) > 0) p2_colsums(j) += 1;

      }

    }

    arma::uvec case_high_risk = find(p1_colsums == n_target);
    arma::uvec comp_high_risk = find(p2_colsums == n_target);

    List res = List::create(Named("case_high_risk") = case_high_risk,
                            Named("comp_high_risk") = comp_high_risk);
    return(res);

  }

}

// [[Rcpp::export]]
List compute_dif_vecs(arma::mat case_genetic_data, arma::mat comp_genetic_data, arma::mat case_comp_dif,
                      arma::mat cases_minus_complements, arma::mat both_one_mat,
                      arma::vec weight_lookup, double n_different_snps_weight = 2, double n_both_one_weight = 1,
                      IntegerVector prev_informative_families = IntegerVector::create(NA_INTEGER)){


  // determine whether families are informative for the set of target_snps
  arma::rowvec total_different_snps = arma::sum(case_comp_dif, 0);
  arma::uvec informative_families = find(total_different_snps > 0);

  // check prev_informative_families and update inf families if needed
  bool prev_inf = all(is_na(prev_informative_families));
  if (!prev_inf){

    IntegerVector informative_families_rcpp = wrap(informative_families);
    IntegerVector changed_families = setdiff(informative_families_rcpp, prev_informative_families);
    if (changed_families.length() > 0){

      arma::uvec change_these = as<arma::uvec>(changed_families);
      total_different_snps.elem(change_these).fill(0);

    }

    arma::uvec prev_inf_families = as<arma::uvec>(prev_informative_families);
    informative_families = intersect(informative_families, prev_inf_families);

  }
  arma::uvec uninformative_families = find(total_different_snps == 0);

  // compute weights
  arma::rowvec both_one = sum(both_one_mat, 0);
  IntegerVector weighted_informativeness = wrap(n_both_one_weight * both_one + n_different_snps_weight * total_different_snps);
  arma::uvec weighted_informativeness_idx = as<arma::uvec>(weighted_informativeness);
  arma::vec family_weights = weight_lookup.elem(weighted_informativeness_idx - 1);
  family_weights.elem(uninformative_families).fill(0);
  double invsum_family_weights = 1 / arma::sum(family_weights);

  // compute weighted difference vectors for cases vs complements
  arma::vec std_family_weights = family_weights * invsum_family_weights;
  //arma::mat dif_vecs = cases_minus_complements.each_col() %= std_family_weights;
  arma::mat dif_vecs = cases_minus_complements;

  // take the sum of the case - complement difference vectors over families
  arma::rowvec sum_dif_vecs = arma::sum(dif_vecs, 1).t();

  // determine how many cases and complements actually have the proposed risk set
  arma::rowvec risk_dirs = sign(sum_dif_vecs);
  arma::uvec pos_risk = find(risk_dirs > 0);
  arma::uvec neg_risk = find(risk_dirs <= 0);
  int n_pos = pos_risk.n_elem;
  int n_neg = neg_risk.n_elem;
  int n_target = case_comp_dif.n_rows;

  IntegerVector idx_vec = seq_len(n_target);

  //using helper function to pick out high risk families
  List high_risk = find_high_risk(n_target, n_pos, n_neg, neg_risk, pos_risk, case_genetic_data,
                                  comp_genetic_data);
  arma::uvec case_high_risk = high_risk["case_high_risk"];
  arma::uvec comp_high_risk = high_risk["comp_high_risk"];

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
List chrom_fitness_score(arma::mat case_genetic_data, arma::mat complement_genetic_data, arma::mat case_comp_dif,
                         IntegerVector snp_cols, arma::mat cases_minus_complements, arma::mat both_one_mat,
                         arma::mat block_ld_mat, arma::vec weight_lookup,
                         int n_different_snps_weight = 2, int n_both_one_weight = 1,
                         int n_case_high_risk_thresh = 20, double outlier_sd = 2.5, bool epi_test = false) {

  // pick out the differences for the target snps
  arma::uvec target_snps = as<arma::uvec>(snp_cols) - 1;
  case_comp_dif = case_comp_dif.cols(target_snps).t();
  int n_target = target_snps.n_elem;

  // also subset to target cols for the rest of the inputs
  case_genetic_data = case_genetic_data.cols(target_snps).t();
  complement_genetic_data = complement_genetic_data.cols(target_snps).t();
  cases_minus_complements = cases_minus_complements.cols(target_snps).t();
  both_one_mat = both_one_mat.cols(target_snps).t();

  // compute weighted difference vectors, determine risk related alleles
  List dif_vec_list = compute_dif_vecs(case_genetic_data, complement_genetic_data, case_comp_dif,
                                       cases_minus_complements, both_one_mat,
                                       weight_lookup, n_different_snps_weight, n_both_one_weight);

  // pick out the required pieces from function output
  arma::rowvec sum_dif_vecs = dif_vec_list["sum_dif_vecs"];
  arma::vec family_weights = dif_vec_list["family_weights"];
  double invsum_family_weights = dif_vec_list["invsum_family_weights"];
  arma::uvec informative_families = dif_vec_list["informative_families"];
  int n_pos = dif_vec_list["n_pos"];
  int n_neg = dif_vec_list["n_neg"];
  arma::uvec pos_risk = dif_vec_list["pos_risk"];
  arma::uvec neg_risk = dif_vec_list["neg_risk"];
  arma::uvec case_high_risk = dif_vec_list["case_high_risk"];
  arma::uvec comp_high_risk = dif_vec_list["comp_high_risk"];

  int n_informative_families = informative_families.n_elem;
  arma::uvec case_high_inf = intersect(informative_families, case_high_risk);
  arma::uvec comp_high_inf = intersect(informative_families, comp_high_risk);

  int n_case_high_risk = case_high_inf.n_elem;
  int n_comp_high_risk = comp_high_inf.n_elem;

  // pick out misclassifications via outlier detection, indicating recessive mode of inheritance

  // initialize vector of pattern of inheritance
  CharacterVector risk_set_alleles(n_target, "1+");

  // only applies if we have at least 20 high risk case
  if (n_case_high_risk > n_case_high_risk_thresh) {

    int outlier_count = 0;
    arma::mat case_high_inf_mat = case_genetic_data.cols(case_high_inf);
    arma::mat comp_high_inf_mat = complement_genetic_data.cols(comp_high_inf);
    arma::vec case_high_risk_means = arma::mean(case_high_inf_mat, 1);
    arma::vec case_high_risk_sd = arma::stddev(case_high_inf_mat, 0, 1);
    arma::vec high_outlier_thresh = case_high_risk_means + outlier_sd * case_high_risk_sd;
    arma::vec low_outlier_thresh = case_high_risk_means - outlier_sd * case_high_risk_sd;
    LogicalVector outliers(n_target, false);

    if (n_pos > 0){

      for (int i = 0; i < n_pos; i++){

        arma::uvec this_row(1);
        this_row.fill(pos_risk(i));
        double thresh_val = low_outlier_thresh(i);

        // check cases
        arma::rowvec case_outlier_target = case_high_inf_mat.rows(this_row);
        arma::uvec case_pos_row_outlier = find(case_outlier_target < thresh_val);

        // check comps
        arma::rowvec comp_outlier_target = comp_high_inf_mat.rows(this_row);
        arma::uvec comp_pos_row_outlier = find(comp_outlier_target < thresh_val);

        // recode if outliers
        int n_outliers = case_pos_row_outlier.n_elem + comp_pos_row_outlier.n_elem;
        if (n_outliers > 0){

            risk_set_alleles[pos_risk(i)] = "2";
            outlier_count += 1;

            // cases
            arma::rowvec case_row = case_genetic_data.rows(this_row);
            arma::uvec case_recode_these = find(case_row == 1);
            IntegerVector row_start(case_genetic_data.n_cols, this_row(0));
            arma::uvec row_uvec = as<arma::uvec>(row_start);
            case_genetic_data.submat(row_uvec, case_recode_these).fill(0);

            // complements
            arma::rowvec comp_row = complement_genetic_data.rows(this_row);
            arma::uvec comp_recode_these = find(comp_row == 1);
            complement_genetic_data.submat(row_uvec, comp_recode_these).fill(0);

            // both one mat
            arma::rowvec both_one_row = both_one_mat.rows(this_row);
            arma::uvec both_one_recode_these = find(both_one_row == 1);
            both_one_mat.submat(row_uvec, both_one_recode_these).fill(0);

        }
      }
    }

    if (n_neg > 0){

      for (int i = 0; i < n_neg; i++){

        arma::uvec this_row(1);
        this_row.fill(neg_risk(i));
        double thresh_val = low_outlier_thresh(i);

        // check cases
        arma::rowvec case_outlier_target = case_high_inf_mat.rows(this_row);
        arma::uvec case_neg_row_outlier = find(case_outlier_target > thresh_val);

        // check comps
        arma::rowvec comp_outlier_target = comp_high_inf_mat.rows(this_row);
        arma::uvec comp_neg_row_outlier = find(comp_outlier_target > thresh_val);

        // recode if outliers
        int n_outliers = case_neg_row_outlier.n_elem + comp_neg_row_outlier.n_elem;
        if (n_outliers > 0){

          outlier_count += 1;
          risk_set_alleles[neg_risk(i)] = "2";

          // cases
          arma::rowvec case_row = case_genetic_data.rows(this_row);
          arma::uvec case_recode_these = find(case_row == 1);
          IntegerVector row_start(case_genetic_data.n_cols, this_row(0));
          arma::uvec row_uvec = as<arma::uvec>(row_start);
          case_genetic_data.submat(row_uvec, case_recode_these).fill(2);

          // complements
          arma::rowvec comp_row = complement_genetic_data.rows(this_row);
          arma::uvec comp_recode_these = find(comp_row == 1);
          complement_genetic_data.submat(row_uvec, comp_recode_these).fill(2);

          // both one mat
          arma::rowvec both_one_row = both_one_mat.rows(this_row);
          arma::uvec both_one_recode_these = find(both_one_row == 1);
          both_one_mat.submat(row_uvec, both_one_recode_these).fill(0);

        }
      }
    }

    if (outlier_count > 0){

      //recompute the number of informative families
      cases_minus_complements = sign(case_genetic_data - complement_genetic_data);
      case_comp_dif = abs(cases_minus_complements);
      IntegerVector informative_families_iv = wrap(informative_families);
      dif_vec_list = compute_dif_vecs(case_genetic_data, complement_genetic_data, case_comp_dif, cases_minus_complements, both_one_mat,
                                      weight_lookup, n_different_snps_weight, n_both_one_weight, informative_families_iv);

      // pick out the required pieces
      arma::rowvec sum_dif_vecs_tmp = dif_vec_list["sum_dif_vecs"];
      sum_dif_vecs = sum_dif_vecs_tmp;
      arma::vec family_weights_tmp = dif_vec_list["family_weights"];
      family_weights = family_weights_tmp;
      double invsum_family_weights_tmp = dif_vec_list["invsum_family_weights"];
      invsum_family_weights = invsum_family_weights_tmp;
      arma::uvec informative_families_tmp = dif_vec_list["informative_families"];
      informative_families = informative_families_tmp;
      int n_pos_tmp = dif_vec_list["n_pos"];
      n_pos = n_pos_tmp;
      int n_neg_tmp = dif_vec_list["n_neg"];
      n_neg = n_neg_tmp;
      arma::uvec pos_risk_tmp = dif_vec_list["pos_risk"];
      pos_risk = pos_risk_tmp;
      arma::uvec neg_risk_tmp = dif_vec_list["neg_risk"];
      neg_risk = neg_risk_tmp;
      arma::uvec case_high_risk_tmp = dif_vec_list["case_high_risk"];
      case_high_risk = case_high_risk_tmp;
      arma::uvec comp_high_risk_tmp = dif_vec_list["comp_high_risk"];
      comp_high_risk = comp_high_risk_tmp;

      n_informative_families = informative_families.n_elem;
      case_high_inf = intersect(informative_families, case_high_risk);
      comp_high_inf = intersect(informative_families, comp_high_risk);

      n_case_high_risk = case_high_inf.n_elem;
      n_comp_high_risk = comp_high_inf.n_elem;

    }
  }

  // count the number of risk alleles in those with the full risk set

  double total_case_high_risk_alleles = 0;
  double total_comp_high_risk_alleles = 0;
  arma::mat case_high_inf_mat = case_genetic_data.cols(case_high_inf);
  arma::mat comp_high_inf_mat = complement_genetic_data.cols(comp_high_inf);

  if (n_pos > 0){

    arma::rowvec case_pos_alleles = sum(case_high_inf_mat.rows(pos_risk), 0);
    total_case_high_risk_alleles += sum(case_pos_alleles);

    arma::rowvec comp_pos_alleles = sum(comp_high_inf_mat.rows(pos_risk), 0);
    total_comp_high_risk_alleles += sum(comp_pos_alleles);

  }
  if (n_neg > 0){

    arma::mat case_all2(n_neg, case_high_inf_mat.n_cols);
    case_all2.fill(2);
    arma::rowvec case_neg_alleles = sum(case_all2 - case_high_inf_mat.rows(neg_risk), 0);
    total_case_high_risk_alleles += sum(case_neg_alleles);

    arma::mat comp_all2(n_neg, comp_high_inf_mat.n_cols);
    comp_all2.fill(2);
    arma::rowvec comp_neg_alleles = sum(comp_all2 - comp_high_inf_mat.rows(neg_risk), 0);
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
  arma::mat x = cases_minus_complements.cols(informative_families).t();
  arma::mat x_minus_mu_hat = x - mu_hat_mat;
  arma::mat weighted_x_minus_mu_hat = x_minus_mu_hat.each_col() %= family_weights;
  arma::mat cov_mat = invsum_family_weights * weighted_x_minus_mu_hat.t() * x_minus_mu_hat;

  // set cov elements to zero if SNPs are not in same ld block
  for (int i = 0; i < n_target; i++){
    int this_row = target_snps(i);
    for (int j = 0; j < n_target; j++){
      int this_col = target_snps(j);
      if (block_ld_mat(this_row, this_col) != 1){

        cov_mat(i,j) = 0;

      }
    }
  }

  // get info for function output
  arma::vec elem_vars = sqrt(cov_mat.diag());
  sum_dif_vecs = sum_dif_vecs/elem_vars.t();

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
List chrom_fitness_list(arma::mat case_genetic_data, arma::mat complement_genetic_data, arma::mat case_comp_differences,
                        List chromosome_list, arma::mat cases_minus_complements, arma::mat both_one_mat,
                        arma::mat block_ld_mat, arma::vec weight_lookup,
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



