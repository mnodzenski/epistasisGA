#include <RcppArmadillo.h>
#include <boost/lexical_cast.hpp>
// [[Rcpp::depends(RcppArmadillo, BH, bigmemory)]]
using namespace Rcpp;
using namespace arma;
using boost::lexical_cast;
#include <bigmemory/BigMatrix.h>
#include <bigmemory/MatrixAccessor.hpp>

//////////////////////
// HELPER FUNCTIONS
/////////////////////

// [[Rcpp::export]]
IntegerVector as_integer(CharacterVector x){

  int out_length = x.length();
  IntegerVector out_vec(out_length);
  for (int i = 0; i < out_length; i++){

    out_vec[i] = lexical_cast<int>(x[i]);

  }
  return(out_vec);

}

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
List concat_list(List x, List y){

  List out_list(x.length() + y.length());
  IntegerVector x_pos = seq_len(x.length());
  out_list[x_pos - 1] = x;
  IntegerVector y_pos = seq(x.length() + 1, out_list.length());
  out_list[y_pos - 1] = y;
  return(out_list);

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

    IntegerVector xi = chromosome_list[i];
    int l = 0;

    for (int j = 0; j < comp_these.length(); j++){

      IntegerVector xj = comp_these[j];

      if ((sum(xi == xj) == chrom_size)){

        break;

      } else {

        l++;

      }
    }

    if (l == s){

      out_vec[i] = true;
      comp_these = chromosome_list[out_vec];
      s++;

    }

  }
  return out_vec;
}

// [[Rcpp::export]]
IntegerMatrix subset_matrix_cols(IntegerMatrix in_matrix, IntegerVector cols){

  int n_rows = in_matrix.nrow();
  int n_cols = cols.length();
  IntegerMatrix out_matrix(n_rows, n_cols);

  for (int i = 0; i < n_cols; i++){
    IntegerMatrix::Column original_col = in_matrix(_, cols[i]-1);
    IntegerMatrix::Column new_col = out_matrix(_, i);
    new_col = original_col;
  }
  return(out_matrix);
}

// [[Rcpp::export]]
LogicalMatrix subset_lmatrix_cols(LogicalMatrix in_matrix, IntegerVector cols){

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
Rcpp::IntegerMatrix subset_matrix_rows(Rcpp::IntegerMatrix in_matrix, Rcpp::IntegerVector rows){

  int n_rows = rows.length();
  int n_cols = in_matrix.ncol();
  Rcpp::IntegerMatrix out_matrix(n_rows, n_cols);

  for (int i = 0; i < n_rows; i++){
    Rcpp::IntegerMatrix::Row original_row = in_matrix(rows[i]-1, _);
    Rcpp::IntegerMatrix::Row new_row = out_matrix(i, _);
    new_row = original_row;
  }
  return(out_matrix);
}

// [[Rcpp::export]]
Rcpp::LogicalMatrix subset_lmatrix_rows(Rcpp::LogicalMatrix in_matrix, Rcpp::IntegerVector rows){

  int n_rows = rows.length();
  int n_cols = in_matrix.ncol();
  Rcpp::LogicalMatrix out_matrix(n_rows, n_cols);

  for (int i = 0; i < n_rows; i++){
    Rcpp::LogicalMatrix::Row original_row = in_matrix(rows[i]-1, _);
    Rcpp::LogicalMatrix::Row new_row = out_matrix(i, _);
    new_row = original_row;
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
LogicalMatrix subset_lmatrix(LogicalMatrix in_matrix, IntegerVector rows, IntegerVector cols){

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
NumericVector weighted_sub_colsums(IntegerMatrix x, IntegerMatrix y, IntegerVector target_rows,
                                   IntegerVector target_cols, IntegerVector row_weights){

  int n_cols = target_cols.length();
  int n_rows = target_rows.length();
  NumericVector out_vec(n_cols, 0.0);
  for (int i = 0; i < n_rows; i++){

    int this_row = target_rows[i] - 1;
    int row_weight = row_weights[i];

    for (int j = 0; j < n_cols; j ++){

      int this_col = target_cols[j] - 1;
      int this_diff = x(this_row, this_col) - y(this_row, this_col);
      out_vec[j] += (row_weight*this_diff);

    }
  }
  return(out_vec);
}


// [[Rcpp::export]]
IntegerVector sub_rowsums_start(IntegerMatrix x, IntegerMatrix y, IntegerVector target_cols){

  int n_cols = target_cols.length();
  int n_rows = x.nrow();
  IntegerVector out_vec(n_rows, 0);
  for (int i = 0; i < n_rows; i++){

    for (int j = 0; j < n_cols; j ++){

      int this_col = target_cols[j] - 1;
      if (x(i , this_col) != y(i, this_col)){

        out_vec[i] += 1;

      }

    }
  }
  return(out_vec);
}

// [[Rcpp::export]]
List sub_rowsums_parent_weights(IntegerMatrix x, IntegerMatrix y, IntegerVector target_cols){

  int n_cols = target_cols.length();
  int n_rows = x.nrow();
  IntegerVector out_vec(n_rows, 0);
  IntegerMatrix out_mat(n_rows, n_cols);
  for (int i = 0; i < n_rows; i++){

    for (int j = 0; j < n_cols; j ++){

      int this_col = target_cols[j] - 1;
      if (x(i , this_col) != y(i, this_col)){

        out_vec[i] += 1;
        out_mat(i, j) = 1;

      }

    }
  }
  List out_list = List::create(Named("n_differences_vec") = out_vec,
                               Named("difference_mat") = out_mat);
  return(out_list);
}

// [[Rcpp::export]]
IntegerVector sub_colsums(IntegerMatrix in_mat, IntegerVector target_rows,
                          IntegerVector target_cols, int target_val){

  int n_cols = target_cols.length();
  int n_rows = target_rows.length();
  IntegerVector out_vec(n_cols, 0);
  for (int i = 0; i < n_rows; i++){

    int this_row = target_rows[i] - 1;

    for (int j = 0; j < n_cols; j ++){

      int this_col = target_cols[j] - 1;
      int this_val = in_mat(this_row, this_col);
      if (this_val == target_val){

        out_vec[j] += 1;

      }

    }
  }
  return(out_vec);
}

// [[Rcpp::export]]
IntegerVector sub_rowsums_both_one(IntegerMatrix x, IntegerMatrix y, IntegerVector target_rows,
                                   IntegerVector target_cols){

  int n_cols = target_cols.length();
  int n_rows = target_rows.length();
  IntegerVector out_vec(n_rows, 0);
  for (int i = 0; i < n_rows; i++){

    int this_row = target_rows[i] - 1;

    for (int j = 0; j < n_cols; j ++){

      int this_col = target_cols[j] - 1;
      if ((x(this_row , this_col) == 1) & (y(this_row, this_col) == 1)){

        out_vec[i] += 1;

      }

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
NumericVector sub_colmeans(IntegerMatrix in_mat, IntegerVector target_rows,
                           IntegerVector target_cols, int target_val){

  IntegerVector target_colsums = sub_colsums(in_mat, target_rows, target_cols, target_val);
  NumericVector out_vec = as<NumericVector>(target_colsums);
  double n = target_rows.length();
  out_vec = out_vec / n;
  return(out_vec);

}

// [[Rcpp::export]]
ListOf<IntegerMatrix> split_int_mat(IntegerMatrix in_mat, IntegerVector in_vec, IntegerVector uni_in_vec){

  // initiate list to be output
  int n_levels = uni_in_vec.length();
  List out_list(n_levels);
  IntegerVector all_rows = seq_along(in_vec);

  // loop over values of the in_vec
  for (int i = 0; i < n_levels; i++){

    int target_level = uni_in_vec[i];
    LogicalVector these_rows_l = in_vec == target_level;
    IntegerVector these_rows = all_rows[these_rows_l];
    IntegerMatrix target_level_mat = subset_matrix_rows(in_mat, these_rows);
    out_list[i] = target_level_mat;

  }

  // return final list
  return(out_list);

}

// [[Rcpp::export]]
ListOf<LogicalMatrix> split_logical_mat(LogicalMatrix in_mat, IntegerVector in_vec, IntegerVector uni_in_vec){

  // initiate list to be output
  int n_levels = uni_in_vec.length();
  List out_list(n_levels);
  IntegerVector all_rows = seq_along(in_vec);

  // loop over values of the in_vec
  for (int i = 0; i < n_levels; i++){

    int target_level = uni_in_vec[i];
    LogicalVector these_rows_l = in_vec == target_level;
    IntegerVector these_rows = all_rows[these_rows_l];
    LogicalMatrix target_level_mat = subset_lmatrix_rows(in_mat, these_rows);
    out_list[i] = target_level_mat;

  }

  // return final list
  return(out_list);

}

// [[Rcpp::export]]
IntegerVector get_target_snps_ld_blocks(IntegerVector target_snps_in,
                                        IntegerVector ld_block_vec){

  int n_target = target_snps_in.length();
  IntegerVector target_snps_block(n_target);

  for (int i = 0; i < n_target; i++){

    int target_snp = target_snps_in[i];

    for (int j = 0; j < ld_block_vec.length(); j++){

      int block_upper_limit = ld_block_vec[j];

      if (target_snp <= block_upper_limit){

        target_snps_block[i] = j;
        break;

      }

    }

  }
  return(target_snps_block);

}



////////////////////////////////////////////////////////////////////
// The following functions actually implement the GADGETS method
///////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////
// function to pick out the families with case or complement with the full risk set
////////////////////////////////////////////////////////////////////////////////////

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

/////////////////////////////////////////////////////////////////////////////
// function to compute the difference vectors in computing the fitness score
/////////////////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
List compute_dif_vecs(IntegerMatrix case_genetic_data, IntegerMatrix comp_genetic_data, IntegerVector target_snps,
                      IntegerVector both_one_snps, IntegerVector weight_lookup, int n_different_snps_weight = 2,
                      int n_both_one_weight = 1){

  // determine whether families are informative for the set of target_snps
  IntegerVector total_different_snps = sub_rowsums_start(case_genetic_data, comp_genetic_data, target_snps);
  LogicalVector informative_families_l = total_different_snps != 0;
  if (sum(informative_families_l) == 0){

    List res = List::create(Named("no_informative_families") = true);
    return(res);

  } else {

    IntegerVector family_idx = seq_len(total_different_snps.length());
    IntegerVector informative_families = family_idx[informative_families_l];

    //LogicalVector uninformative_families_l = total_different_snps == 0;
    total_different_snps = total_different_snps[total_different_snps > 0];

    // compute weights
    IntegerVector both_one = sub_rowsums_both_one(case_genetic_data, comp_genetic_data, informative_families, both_one_snps);
    IntegerVector weighted_informativeness = n_both_one_weight * both_one + n_different_snps_weight * total_different_snps;
    IntegerVector family_weights(informative_families_l.length(), 0);
    for (int i = 0; i < informative_families.length(); i++){

      int inf_family_i = informative_families[i];
      int weight_i = weighted_informativeness[i];
      family_weights[inf_family_i - 1] = weight_lookup[weight_i - 1];

    }
    double sum_family_weights = sum(family_weights);
    double invsum_family_weights = 1 / sum_family_weights;
    IntegerVector inf_family_weights = family_weights[informative_families_l];
    NumericVector sum_dif_vecs = weighted_sub_colsums(case_genetic_data, comp_genetic_data, informative_families,
                                                      target_snps, inf_family_weights);

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
    LogicalVector both_high_risk = case_high_risk & comp_high_risk;
    IntegerVector case_high_inf_rows = informative_families[case_high_risk & !comp_high_risk];
    IntegerVector comp_high_inf_rows = informative_families[comp_high_risk & ! case_high_risk];

    // remove families where both case and complement have high risk genotype
    case_high_risk = case_high_risk[!both_high_risk];
    comp_high_risk = comp_high_risk[!both_high_risk];

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
                            Named("case_high_inf_rows") = case_high_inf_rows,
                            Named("comp_high_inf_rows") = comp_high_inf_rows,
                            Named("informative_families") = informative_families,
                            Named("no_informative_families") = false);
    return(res);


  }
}


/////////////////////////////////////////////////////////
// Function to compute the weights for parents only GxE
/////////////////////////////////////////////////////////

// [[Rcpp::export]]
List compute_parent_weights(IntegerMatrix case_genetic_data, IntegerMatrix comp_genetic_data, IntegerVector target_snps,
                      IntegerVector weight_lookup, int n_different_snps_weight = 1, int n_both_one_weight = 0){

  // determine whether families are informative for the set of target_snps
  List total_different_snps_list = sub_rowsums_parent_weights(case_genetic_data, comp_genetic_data, target_snps);
  IntegerVector total_different_snps = total_different_snps_list["n_differences_vec"];
  int n_fams = total_different_snps.length();
  LogicalVector informative_families_l = total_different_snps != 0;
  if (sum(informative_families_l) == 0){

    IntegerVector family_weights(n_fams, 0);
    total_different_snps_list["family_weights"] = family_weights;
    return(total_different_snps_list);

  } else {

    IntegerVector family_idx = seq_len(n_fams);
    IntegerVector informative_families = family_idx[informative_families_l];

    total_different_snps = total_different_snps[total_different_snps > 0];

    // compute weights
    IntegerVector family_weights(n_fams, 0);
    IntegerVector weighted_informativeness(total_different_snps.length());
    if (n_both_one_weight > 0){

      IntegerVector both_one = sub_rowsums_both_one(case_genetic_data, comp_genetic_data, informative_families, target_snps);
      weighted_informativeness = n_both_one_weight * both_one + n_different_snps_weight * total_different_snps;

    } else {

      weighted_informativeness = n_different_snps_weight * total_different_snps;

    }

    for (int i = 0; i < informative_families.length(); i++){

      int inf_family_i = informative_families[i];
      int weight_i = weighted_informativeness[i];
      family_weights[inf_family_i - 1] = weight_lookup[weight_i - 1];

    }

    // return vector of family weights
    total_different_snps_list["family_weights"] = family_weights;
    return(total_different_snps_list);

  }
}

////////////////////////////////////////////////////////////////
// Function to compute the fitness score
///////////////////////////////////////////////////////////////

// [[Rcpp::export]]
List chrom_fitness_score(IntegerMatrix case_genetic_data_in, IntegerMatrix complement_genetic_data_in, IntegerVector target_snps_in,
                         IntegerVector ld_block_vec, IntegerVector weight_lookup, int n_different_snps_weight = 2, int n_both_one_weight = 1,
                         double recessive_ref_prop = 0.75, double recode_test_stat = 1.64, bool epi_test = false, bool GxE = false,
                         Nullable<IntegerVector> both_one_snps_in = R_NilValue) {

  // need to deep copy inputs to avoid overwriting if we recode recessives
  IntegerMatrix case_genetic_data = subset_matrix_cols(case_genetic_data_in, target_snps_in);
  IntegerMatrix complement_genetic_data = subset_matrix_cols(complement_genetic_data_in, target_snps_in);

  // now redefining target snps
  int n_target = target_snps_in.length();
  IntegerVector target_snps = seq_len(n_target);

  // determine which SNPs need to be checked for both one
  IntegerVector both_one_snps;
  if (both_one_snps_in.isNotNull()){

    both_one_snps = both_one_snps_in;

  } else {

    both_one_snps = target_snps;

  }

  // compute weighted difference vectors, determine risk related alleles
  List dif_vec_list = compute_dif_vecs(case_genetic_data, complement_genetic_data, target_snps, both_one_snps,
                                       weight_lookup, n_different_snps_weight, n_both_one_weight);

  // pick out the required pieces from function output
  bool no_informative_families = dif_vec_list["no_informative_families"];

  // require at least one informative family
  if (no_informative_families){

    if (GxE){

      List res = List::create(Named("no_informative_families") = no_informative_families);
      return(res);

    } else {

      double fitness_score = pow(10, -10);
      NumericVector sum_dif_vecs(n_target, 1.0);
      double q = pow(10, -10);
      CharacterVector risk_set_alleles(n_target, "1+");
      double n_case_high_risk = 0.0;
      double n_comp_high_risk = 0.0;
      List res = List::create(Named("fitness_score") = fitness_score,
                              Named("sum_dif_vecs") = sum_dif_vecs,
                              Named("q") = q,
                              Named("risk_set_alleles") = risk_set_alleles,
                              Named("n_case_risk_geno") = n_case_high_risk,
                              Named("n_comp_risk_geno") = n_comp_high_risk);
      return(res);

    }

  } else {

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
    IntegerVector case_high_inf_rows = dif_vec_list["case_high_inf_rows"];
    IntegerVector comp_high_inf_rows = dif_vec_list["comp_high_inf_rows"];
    double n_case_high_risk = sum(case_high_risk);
    double n_comp_high_risk = sum(comp_high_risk);

    // initialize vector of pattern of inheritance
    CharacterVector risk_set_alleles(target_snps.length(), "1+");
    LogicalVector recoded(target_snps.length(), false);

    // examine recessive SNPs
    int recessive_count = 0;

    if (n_case_high_risk > 0){

      if (n_pos > 0){

        IntegerVector pos_cols = target_snps[pos_risk];
        NumericVector prop2_case = sub_colmeans(case_genetic_data, case_high_inf_rows, pos_cols, 2);
        NumericVector prop2_comp = sub_colmeans(complement_genetic_data, comp_high_inf_rows, pos_cols, 2);

        for (int i = 0; i < n_pos; i++){

          int this_col = pos_cols[i] - 1;
          double phat_case = prop2_case[i];
          double phat_comp = prop2_comp[i];
          double xi1 = (phat_case - recessive_ref_prop)/sqrt((recessive_ref_prop*(1 - recessive_ref_prop)/n_case_high_risk));
          double xi2 = recode_test_stat;

          // check expected sample sizes, and compare to comps where possible
          double case1 = n_case_high_risk*phat_case;
          double case0 = n_case_high_risk*(1 - phat_case);
          double comp1 = n_comp_high_risk*phat_comp;
          double comp0 = n_comp_high_risk*(1 - phat_comp);
          if ((case1 >= 5) & (case0 >= 5) & (comp1 >= 5) & (comp0 >= 5)){

            double phat_pooled = (n_case_high_risk*phat_case + n_comp_high_risk*phat_comp)/(n_case_high_risk + n_comp_high_risk);
            xi2 = (phat_case - phat_comp)/sqrt((phat_pooled*(1 - phat_pooled)*(1/n_case_high_risk + 1/n_comp_high_risk)));

          }

          // if test stat exceeds threshold, recode
          if ((xi1 >= recode_test_stat) & (xi2 >= recode_test_stat)){

            // increment the counter
            recessive_count += 1;

            // indicate we need two risk alleles
            int this_index = pos_cols[i] - 1;
            risk_set_alleles[this_index] = "2";
            recoded[this_index] = true;

            // loop over informative families
            for (int n = 0; n < informative_families.length(); n++){

              int family_row = informative_families[n] - 1;

              // recode cases
              int case_val = case_genetic_data(family_row, this_col);
              if (case_val == 1){

                case_genetic_data(family_row, this_col) = 0;

              } else if (case_val == 2){

                case_genetic_data(family_row, this_col) = 1;

              }

              //recode complements
              int comp_val = complement_genetic_data(family_row, this_col);
              if (comp_val == 1){

                complement_genetic_data(family_row, this_col) = 0;

              } else if (comp_val == 2){

                complement_genetic_data(family_row, this_col) = 1;

              }

            }

          }

        }

      }

      if (n_neg > 0){

        IntegerVector neg_cols = target_snps[neg_risk];
        NumericVector prop0_case = sub_colmeans(case_genetic_data, case_high_inf_rows, neg_cols, 0);
        NumericVector prop0_comp = sub_colmeans(complement_genetic_data, comp_high_inf_rows, neg_cols, 0);

        for (int i = 0; i < n_neg; i++){

          int this_col = neg_cols[i] - 1;
          double phat_case = prop0_case[i];
          double phat_comp = prop0_comp[i];
          double xi1 = (phat_case - recessive_ref_prop)/sqrt((recessive_ref_prop*(1 - recessive_ref_prop)/n_case_high_risk));
          double xi2 = recode_test_stat;

          // check expected sample sizes, and compare to comps where possible
          double case1 = n_case_high_risk*phat_case;
          double case0 = n_case_high_risk*(1 - phat_case);
          double comp1 = n_comp_high_risk*phat_comp;
          double comp0 = n_comp_high_risk*(1 - phat_comp);

          if ((case1 >= 5) & (case0 >= 5) & (comp1 >= 5) & (comp0 >= 5)){

            double phat_pooled = (n_case_high_risk*phat_case + n_comp_high_risk*phat_comp)/(n_case_high_risk + n_comp_high_risk);
            xi2 = (phat_case - phat_comp)/sqrt((phat_pooled*(1 - phat_pooled)*(1/n_case_high_risk + 1/n_comp_high_risk)));

          }

          // if test stat exceeds threshold, recode
          if ((xi1 >= recode_test_stat) & (xi2 >= recode_test_stat)){

            // increment the counter
            recessive_count += 1;

            // indicate we need two risk alleles
            int this_index = neg_cols[i] - 1;
            risk_set_alleles[this_index] = "2";
            recoded[this_index] = true;

            // loop over informative families
            for (int n = 0; n < informative_families.length(); n++){

              int family_row = informative_families[n] - 1;

              // recode cases
              int case_val = case_genetic_data(family_row, this_col);
              if (case_val == 1){

                case_genetic_data(family_row, this_col) = 2;

              } else if (case_val == 0){

                case_genetic_data(family_row, this_col) = 1;

              }
              //recode complements
              int comp_val = complement_genetic_data(family_row, this_col);
              if (comp_val == 1){

                complement_genetic_data(family_row, this_col) = 2;

              } else if (comp_val == 0){

                complement_genetic_data(family_row, this_col) = 1;

              }

            }

          }

        }

      }

    }

    // if there are recessive snps, recompute the weights and associated statistics
    if (recessive_count > 0) {

      both_one_snps = target_snps[!recoded];

      //recompute the number of informative families
      List dif_vec_list = compute_dif_vecs(case_genetic_data, complement_genetic_data, target_snps, both_one_snps,
                                           weight_lookup, n_different_snps_weight, n_both_one_weight);

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
      n_case_high_risk = sum(case_high_risk);
      n_comp_high_risk = sum(comp_high_risk);

    }

    // compute scaling factor
    double q;
    double total_high_risk = n_case_high_risk + n_comp_high_risk;
    if ( (total_high_risk == 0) | (n_case_high_risk == 0) |
         R_isnancpp(n_case_high_risk) | R_isnancpp(n_comp_high_risk) |
         !arma::is_finite(n_case_high_risk) | !arma::is_finite(n_comp_high_risk)){

         q = pow(10, -10);

    } else {

      q = n_case_high_risk/total_high_risk;
      if (q <= 0.5){

        q = pow(10, -10);

      }

    }

    // compute t2
    sum_dif_vecs = q*sum_dif_vecs;
    arma::rowvec mu_hat = as<arma::rowvec>(sum_dif_vecs);
    arma::mat mu_hat_mat(case_genetic_data.nrow(), n_target);
    for (int i = 0; i < case_genetic_data.nrow(); i++){
      mu_hat_mat.row(i) = mu_hat;
    }
    arma::mat x = as<arma::mat>(case_genetic_data) - as<arma::mat>(complement_genetic_data);
    arma::vec inf_family_weights = as<arma::vec>(family_weights);
    arma::mat x_minus_mu_hat = x - mu_hat_mat;
    arma::mat weighted_x_minus_mu_hat = x_minus_mu_hat;
    weighted_x_minus_mu_hat.each_col() %= inf_family_weights;
    arma::mat cov_mat = invsum_family_weights * trans(weighted_x_minus_mu_hat) * x_minus_mu_hat;

    // set cov elements to zero if SNPs are not in same ld block
    if (ld_block_vec.length() > 1){

      IntegerVector target_snps_block = get_target_snps_ld_blocks(target_snps_in, ld_block_vec);
      IntegerVector uni_target_blocks = unique(target_snps_block);

      if (uni_target_blocks.length() > 1){

        IntegerVector target_block_pos = seq_along(target_snps_block);

        for (unsigned int i = 0; i < uni_target_blocks.length(); i++){

          unsigned int block_i = uni_target_blocks[i];
          LogicalVector these_pos_l = target_snps_block == block_i;
          IntegerVector these_pos = target_block_pos[these_pos_l];
          IntegerVector not_these_pos = setdiff(target_block_pos, these_pos);
          for (unsigned int j = 0; j < these_pos.length(); j++){

            unsigned int these_pos_j = these_pos[j];

            for (unsigned int k = 0; k < not_these_pos.length(); k++){

              unsigned int not_these_pos_k = not_these_pos[k];
              cov_mat(not_these_pos_k - 1, these_pos_j - 1) = 0;

            }

          }

        }

      }

    }

    // get info for function output
    NumericVector elem_vars = wrap(sqrt(cov_mat.diag()));
    sum_dif_vecs = sum_dif_vecs/elem_vars;

    // if no variance, just make the element small
    sum_dif_vecs[elem_vars == 0] = pow(10, -10);

    // compute fitness score
    double fitness_score = (1/(1000*invsum_family_weights)) * as_scalar(mu_hat * arma::pinv(cov_mat) * mu_hat.t());

    // if the fitness score is zero or undefined (either due to zero variance or mean), reset to small number
    if ( (fitness_score <= 0) | R_isnancpp(fitness_score) | !arma::is_finite(fitness_score) ){
      fitness_score = pow(10, -10);
    }

    // if desired, return the required information for the epistasis test
    if (epi_test & !GxE){

      List res = List::create(Named("fitness_score") = fitness_score,
                              Named("sum_dif_vecs") = sum_dif_vecs,
                              Named("q") = q,
                              Named("risk_set_alleles") = risk_set_alleles,
                              Named("n_case_risk_geno") = n_case_high_risk,
                              Named("n_comp_risk_geno") = n_comp_high_risk,
                              Named("inf_families") = informative_families);
      return(res);

    } else if (GxE & !epi_test){

        List res = List::create(Named("xbar") = mu_hat,
                                Named("sum_dif_vecs") = sum_dif_vecs,
                                Named("sigma") = cov_mat,
                                Named("w") = 1/invsum_family_weights,
                                Named("family_weights") = family_weights,
                                Named("risk_set_alleles") = risk_set_alleles,
                                Named("n_case_risk_geno") = n_case_high_risk,
                                Named("n_comp_risk_geno") = n_comp_high_risk);
        return(res);

    } else if (GxE & epi_test){

      List res = List::create(Named("xbar") = mu_hat,
                              Named("sum_dif_vecs") = sum_dif_vecs,
                              Named("sigma") = cov_mat,
                              Named("w") = 1/invsum_family_weights,
                              Named("family_weights") = family_weights,
                              Named("risk_set_alleles") = risk_set_alleles,
                              Named("n_case_risk_geno") = n_case_high_risk,
                              Named("n_comp_risk_geno") = n_comp_high_risk,
                              Named("inf_families") = informative_families);
      return(res);

    } else {

      List res = List::create(Named("fitness_score") = fitness_score,
                              Named("sum_dif_vecs") = sum_dif_vecs,
                              Named("q") = q,
                              Named("risk_set_alleles") = risk_set_alleles,
                              Named("n_case_risk_geno") = n_case_high_risk,
                              Named("n_comp_risk_geno") = n_comp_high_risk);
      return(res);

    }

  }

}

///////////////////////////////////////////////////////////
// fitness score for gene by environment interactions
/////////////////////////////////////////////////////////

// [[Rcpp::export]]
List GxE_fitness_score_parents_only(ListOf<IntegerMatrix> case_genetic_data_list, ListOf<IntegerMatrix> complement_genetic_data_list,
                                    IntegerVector target_snps, IntegerVector weight_lookup,
                                    int n_different_snps_weight = 1, int n_both_one_weight = 0){

  // divide the input data based on exposure and get components for fitness score //
  int n_levels = case_genetic_data_list.size();
  List family_weights_list(n_levels);

  for (int i = 0; i < n_levels; i++){

    IntegerMatrix case_genetic_data = case_genetic_data_list[i];
    IntegerMatrix complement_genetic_data = complement_genetic_data_list[i];
    List family_weights_list_i = compute_parent_weights(case_genetic_data, complement_genetic_data, target_snps,
                                                        weight_lookup, n_different_snps_weight, n_both_one_weight);
    family_weights_list[i] = family_weights_list_i;

  }

  // linear models for each pair of exposure levels  //
  int n_pairs = Rf_choose(n_levels, 2);
  List pair_scores_list(n_pairs);
  int pos = 0;
  int max_score_pos = 0;
  double max_score;

  for (int i = 0; i < n_levels - 1; i++){

    for (int j = i + 1; j < n_levels; j++){

      double s;

      // results for specific exposure levels
      List family_weights_list_i = family_weights_list[i];
      List family_weights_list_j = family_weights_list[j];
      arma::vec weights_exp1 = family_weights_list_i["family_weights"];
      arma::vec weights_exp2 = family_weights_list_j["family_weights"];

      if (arma::sum(weights_exp1) == 0 | arma::sum(weights_exp2) == 0){

        s = pow(10, -10);
        int n_target = target_snps.length();
        arma::vec abs_tstats(n_target, fill::ones);
        pair_scores_list[pos] = List::create(Named("fitness_score") = s,
                                             Named("sum_dif_vecs") = abs_tstats.t());

      } else {

        arma::mat informativeness_mat_exp1 = family_weights_list_i["difference_mat"];
        arma::mat informativeness_mat_exp2 = family_weights_list_j["difference_mat"];
        arma::vec y_exp1(weights_exp1.n_elem, fill::zeros);
        arma::vec y_exp2(weights_exp2.n_elem, fill::ones);

        // prepare for linear model
        arma::vec weights = join_cols(weights_exp1, weights_exp2);
        arma::mat informativeness_mat = join_cols(informativeness_mat_exp1, informativeness_mat_exp2);

        // 'outcome' vector
        arma::vec y= join_cols(y_exp1, y_exp2);
        int n_fam = y.n_elem;

        // intercept only mat
        arma::vec x0(n_fam, fill::ones);

        // full model mat
        arma::mat x = join_rows(x0, informativeness_mat);

        // now transform wls to ols
        arma::vec sqrt_weights = arma::sqrt(weights);
        x.each_col() %= sqrt_weights;
        y.each_col() %= sqrt_weights;
        double dfe = x.n_rows - x.n_cols;
        double dfm = x.n_cols - 1;

        // get beta/se, f-stat
        arma::vec betas = solve(x, y, solve_opts::fast);
        arma::vec resid = y - x*betas;
        double mse = arma::as_scalar(arma::trans(resid)*resid/dfe);
        arma::mat cov_mat = mse * arma::pinv(arma::trans(x)*x);
        arma::vec se = arma::sqrt(arma::diagvec(cov_mat));
        arma::vec abs_tstats = arma::abs(betas/se);

        // compute overall test stat
        betas.head(1).fill(0);
        cov_mat.col(0).fill(0);
        cov_mat.row(0).fill(0);
        double s = arma::as_scalar(betas.t() * arma::pinv(cov_mat) * betas);

        // if the fitness score is zero or undefined (either due to zero variance or mean), reset to small number
        if ( (s <= 0) | R_isnancpp(s) | !arma::is_finite(s) ){
          s = pow(10, -10);
        }

        // make sure the t-stats are positive real numbers
        double fill_val = pow(10, -10);
        arma::uvec se_zero = arma::find(se == 0);
        abs_tstats.elem(se_zero).fill(fill_val);
        arma::uvec beta_zero = arma::find(betas == 0);
        abs_tstats.elem(beta_zero).fill(fill_val);

        // get rid of intercept abs t-stat
        abs_tstats = abs_tstats.tail(dfm);

        // store in list
        pair_scores_list[pos] = List::create(Named("fitness_score") = s,
                                             Named("sum_dif_vecs") = abs_tstats.t());

      }

      if (pos == 0){

        max_score = s;

      }

      if (pos > 0){

        if (s > max_score){

          max_score_pos = pos;

        }

      }
      pos += 1;

    }

  }

  List res = pair_scores_list[max_score_pos];
  return(res);

}

// [[Rcpp::export]]
List GxE_fitness_score(ListOf<IntegerMatrix> case_genetic_data_list, ListOf<IntegerMatrix> complement_genetic_data_list,
                       IntegerVector target_snps, IntegerVector ld_block_vec, IntegerVector weight_lookup,
                       IntegerVector exposure_levels, IntegerVector exposure_risk_levels, int n_different_snps_weight = 2,
                       int n_both_one_weight = 1, double recessive_ref_prop = 0.75, double recode_test_stat = 1.64,
                       bool epi_test = false, bool check_risk = true, int use_parents = 1){

  // divide the input data based on exposure and get components for fitness score //
  List score_by_exposure(exposure_levels.length());

  for (int i = 0; i < score_by_exposure.length(); i++){

    IntegerMatrix case_genetic_data = case_genetic_data_list[i];
    IntegerMatrix complement_genetic_data = complement_genetic_data_list[i];
    score_by_exposure[i] = chrom_fitness_score(case_genetic_data, complement_genetic_data, target_snps,
                                               ld_block_vec, weight_lookup, n_different_snps_weight, n_both_one_weight,
                                               recessive_ref_prop, recode_test_stat, epi_test, true);

  }

  // check lengths of difference vectors
  NumericVector xbar_length_by_exposure(score_by_exposure.length());

  for (int i = 0; i < exposure_levels.length(); i++){

    // mean difference vector length //
    List exposure_i_res = score_by_exposure[i];
    bool no_inf_families = exposure_i_res.containsElementNamed("no_informative_families");
    if (no_inf_families){

      xbar_length_by_exposure[i] = 0.0;

    } else {

      arma::rowvec xbar = exposure_i_res["xbar"];
      double xbar_length_sq = as_scalar(xbar * xbar.t());
      double xbar_length = pow(xbar_length_sq, 0.5);
      xbar_length_by_exposure[i] = xbar_length;

    }

  }

  // decide on risk order
  NumericVector sorted_xbar_length_by_exposure = clone(xbar_length_by_exposure).sort();
  IntegerVector xbar_length_order = match(xbar_length_by_exposure, sorted_xbar_length_by_exposure);

  // if needed, check the risk model
  bool correct_risk_model = true;
  if (check_risk){

    for (int i = 0; i < exposure_levels.length() - 1; i++){

      for (int j = i + 1; j < exposure_levels.length(); j++){

        // risk levels
        int risk_level_i = exposure_risk_levels[i];
        int risk_level_j = exposure_risk_levels[j];

        // vector lengths
        double xbar_length_i = xbar_length_by_exposure[i];
        double xbar_length_j = xbar_length_by_exposure[j];

        // make sure the fitness scores correspond to the hypothesized
        // risk levels (recall risk level 1 is the lowest risk level)
        if ((risk_level_i > risk_level_j) & (xbar_length_j > xbar_length_i)){

          correct_risk_model = false;
          break;

        } else if ((risk_level_j > risk_level_i) & (xbar_length_i > xbar_length_j)){

          correct_risk_model = false;
          break;

        }

      }

      if (!correct_risk_model){

        break;

      }

    }

  }

  // initialize results object
  List res(4);
  if (!correct_risk_model){

    double s = pow(10, -10);
    int n_target = target_snps.length();
    NumericVector std_diff_vecs(n_target, 1.0);
    res = List::create(Named("fitness_score") = s,
                            Named("sum_dif_vecs") = std_diff_vecs,
                            Named("xbar_length_order") = xbar_length_order,
                            Named("score_by_exposure") = score_by_exposure);
    return(res);

  } else {

    // compute two sample hotelling for each pairwise comparison //
    int n_pairs = Rf_choose(exposure_levels.length(), 2);
    List pair_scores_list(n_pairs);
    int pos = 0;
    int max_score_pos = 0;
    double max_score;

    for (int i = 0; i < exposure_levels.length() - 1; i++){

      for (int j = i + 1; j < exposure_levels.length(); j++){

        // results for specific exposure levels
        List exp1_list = score_by_exposure[i];
        List exp2_list = score_by_exposure[j];

        // make sure we have at least one informative
        // family for each exposure level
        bool exp1_no_inf_families = exp1_list.containsElementNamed("no_informative_families");
        bool exp2_no_inf_families = exp2_list.containsElementNamed("no_informative_families");

        // initialize output info
        double s = pow(10, -10);;
        int n_target = target_snps.length();
        NumericVector std_diff_vecs(n_target, 1.0);

        // return defaults if at least one level has no informative families
        // otherwise two sample hotelling t2
        if (!(exp1_no_inf_families | exp2_no_inf_families)){

          // mean difference vector //
          arma::rowvec xbar1 = exp1_list["xbar"];
          arma::rowvec xbar2 = exp2_list["xbar"];
          arma::rowvec xbar_diff = xbar1 - xbar2;

          // cov mat //
          double w1 = exp1_list["w"];
          arma::vec family_weights1 = exp1_list["family_weights"];
          arma::mat sigma1 = exp1_list["sigma"];
          double w2 = exp2_list["w"];
          arma::vec family_weights2 = exp2_list["family_weights"];
          arma::mat sigma2 = exp2_list["sigma"];
          arma::mat sigma_hat = (1/w1)*sigma1 + (1/w2)*sigma2;
          //arma::mat sigma_hat = (w1*sigma1 + w2*sigma2)/(w1 + w2);

          // two sample modified hotelling stat //
          //double weight_scalar = (w1*w2)/(w1 + w2);
          //s = (weight_scalar / 1000) * as_scalar(xbar_diff * arma::pinv(sigma_hat) * xbar_diff.t());
          arma::mat sigma_hat_inv = arma::pinv(sigma_hat);
          s = as_scalar(xbar_diff * sigma_hat_inv * xbar_diff.t());
          s = s / 1000;

          // if parents are to be used, compute the t-stat and multiply
          if (use_parents == 1){

            arma::colvec family_weights_vec = arma::join_cols(family_weights1, family_weights2);
            int n_exp1 = family_weights1.size();
            int n_exp2 = family_weights2.size();
            int n = n_exp1 + n_exp2;
            arma::colvec exp1_vals = arma::colvec(n_exp1, fill::zeros);
            arma::colvec exp2_vals = arma::colvec(n_exp2, fill::ones);
            arma::colvec y = join_cols(exp1_vals, exp2_vals);
            arma::colvec intercept = arma::colvec(n, fill::ones);
            arma::mat X = arma::join_horiz(intercept, family_weights_vec);
            arma::vec coef = solve(X, y, solve_opts::fast);
            double beta_hat = arma::as_scalar(coef[1]);
            arma::vec resid = y - X*coef;
            double sig2 = arma::as_scalar(arma::trans(resid)*resid/(n - 2));
            arma::vec se = arma::sqrt(sig2 * arma::diagvec(arma::inv(arma::trans(X)*X)));
            double se_hat = arma::as_scalar(se[1]);
            double abs_t = std::abs(beta_hat / se_hat);
            s = abs_t * s;

          }

          // if the fitness score is zero or undefined (either due to zero variance or mean), reset to small number
          if ( (s <= 0) | R_isnancpp(s) | !arma::is_finite(s) ){
            s = pow(10, -10);
          }

          // prepare results
          NumericVector se = wrap(sqrt(sigma_hat.diag()));
          NumericVector xbar_diff_n = wrap(xbar_diff);
          NumericVector std_diff_vecs_tmp = xbar_diff_n/se;
          std_diff_vecs = std_diff_vecs_tmp;
          std_diff_vecs[se == 0] = pow(10, -10);

        }
        pair_scores_list[pos] = List::create(Named("fitness_score") = s,
                                             Named("sum_dif_vecs") = std_diff_vecs,
                                             Named("xbar_length_order") = xbar_length_order,
                                             Named("score_by_exposure") = score_by_exposure);

        if (pos == 0){

          max_score = s;

        }

        if (pos > 0){

          if (s > max_score){

            max_score_pos = pos;

          }

        }
        pos += 1;

      }

    }

    res = pair_scores_list[max_score_pos];

  }

  return(res);

}

////////////////////////////////////////////////////////////////////////////////
// helper function to apply the fitness score function to a list of chromosomes
////////////////////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
List chrom_fitness_list(IntegerMatrix case_genetic_data, IntegerMatrix complement_genetic_data, List chromosome_list,
                        IntegerVector ld_block_vec, IntegerVector weight_lookup, int n_different_snps_weight = 2,
                        int n_both_one_weight = 1, double recessive_ref_prop = 0.75, double recode_test_stat = 1.64,
                        bool epi_test = false){

  List scores = chromosome_list.length();
  for (int i = 0; i < chromosome_list.length(); i++){

    IntegerVector target_snps = chromosome_list[i];
    scores[i] = chrom_fitness_score(case_genetic_data, complement_genetic_data, target_snps, ld_block_vec,
                                    weight_lookup, n_different_snps_weight, n_both_one_weight, recessive_ref_prop,
                                    recode_test_stat, epi_test);

  }
  return(scores);

}

///////////////////////////////////////////////////////////////////////////
// Helper function to apply the GxE fitness score to a list of chromosomes
///////////////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
List GxE_fitness_list(List case_genetic_data_list, List complement_genetic_data_list, List chromosome_list,
                        IntegerVector ld_block_vec, IntegerVector weight_lookup, IntegerVector exposure_levels,
                        IntegerVector exposure_risk_levels, int n_different_snps_weight = 2, int n_both_one_weight = 1,
                        double recessive_ref_prop = 0.75, double recode_test_stat = 1.64, bool epi_test = false,
                        bool check_risk = true, int use_parents = 1){

  List scores = chromosome_list.length();
  for (int i = 0; i < chromosome_list.length(); i++){

    IntegerVector target_snps = chromosome_list[i];
    if (use_parents != 2){

      scores[i] = GxE_fitness_score(case_genetic_data_list, complement_genetic_data_list, target_snps,
                                    ld_block_vec, weight_lookup, exposure_levels, exposure_risk_levels, n_different_snps_weight,
                                    n_both_one_weight, recessive_ref_prop, recode_test_stat, epi_test, check_risk, use_parents);
    } else {

      scores[i] = GxE_fitness_score_parents_only(case_genetic_data_list, complement_genetic_data_list, target_snps,
                                                 weight_lookup, n_different_snps_weight, n_both_one_weight);

    }

  }
  return(scores);

}

//////////////////////////////////////////////////////////////
//function to compute fitness score of list of chromosomes
// and parse the output
/////////////////////////////////////////////////////////////

// [[Rcpp::export]]
List compute_population_fitness(IntegerMatrix case_genetic_data, IntegerMatrix complement_genetic_data,
                                IntegerVector ld_block_vec, List chromosome_list,
                                IntegerVector weight_lookup,
                                ListOf<IntegerMatrix> case_genetic_data_list,
                                ListOf<IntegerMatrix> complement_genetic_data_list,
                                IntegerVector exposure_levels,
                                IntegerVector exposure_risk_levels,
                                int n_different_snps_weight = 2,
                                int n_both_one_weight = 1, double recessive_ref_prop = 0.75, double recode_test_stat = 1.64,
                                bool GxE = false, bool check_risk = true, int use_parents = 1){

  // initiate storage object for fitness scores
  List chrom_fitness_score_list;

  // get correct fitness score depending on whether GxE is desired
  if (GxE){

    chrom_fitness_score_list = GxE_fitness_list(case_genetic_data_list, complement_genetic_data_list, chromosome_list,
                                                ld_block_vec, weight_lookup, exposure_levels, exposure_risk_levels,
                                                n_different_snps_weight, n_both_one_weight, recessive_ref_prop, recode_test_stat,
                                                false, check_risk, use_parents);

    // for storing generation information

    // overall fitness score information
    int n_chromosomes = chromosome_list.length();
    NumericVector fitness_scores(n_chromosomes);
    List sum_dif_vecs(n_chromosomes);
    List gen_original_cols(n_chromosomes);
    List exposure_level_info(n_chromosomes);

    // loop over chromosomes and pick out the information
    for (int i = 0; i < n_chromosomes; i++){

      // general info
      List chrom_fitness_score_list_i = chrom_fitness_score_list[i];
      double fs = chrom_fitness_score_list_i["fitness_score"];
      NumericVector sdv = chrom_fitness_score_list_i["sum_dif_vecs"];
      IntegerVector chrom_col_idx = chromosome_list[i];

      fitness_scores[i] = fs;
      sum_dif_vecs[i] = sdv;
      gen_original_cols[i] = chrom_col_idx;
      if (use_parents != 2){

        List exposure_level_res_i = chrom_fitness_score_list_i["score_by_exposure"];
        IntegerVector xbar_length_order_i = chrom_fitness_score_list_i["xbar_length_order"];
        List exposure_info_i =  List::create(Named("score_by_exposure") = exposure_level_res_i,
                                             Named("xbar_length_order") = xbar_length_order_i);
        exposure_level_info[i] = exposure_info_i;

      }

    }

    if (use_parents != 2){

      List res = List::create(Named("chromosome_list") = chromosome_list,
                              Named("fitness_scores") = fitness_scores,
                              Named("sum_dif_vecs") = sum_dif_vecs,
                              Named("gen_original_cols") = gen_original_cols,
                              Named("exposure_level_info") = exposure_level_info);
      return(res);

    } else {

      List res = List::create(Named("chromosome_list") = chromosome_list,
                              Named("fitness_scores") = fitness_scores,
                              Named("sum_dif_vecs") = sum_dif_vecs,
                              Named("gen_original_cols") = gen_original_cols);
      return(res);

    }

  } else {

    chrom_fitness_score_list = chrom_fitness_list(case_genetic_data, complement_genetic_data, chromosome_list, ld_block_vec, weight_lookup,
                                                  n_different_snps_weight, n_both_one_weight, recessive_ref_prop, recode_test_stat, false);

    // for storing generation information
    int n_chromosomes = chromosome_list.length();
    NumericVector fitness_scores(n_chromosomes);
    List sum_dif_vecs(n_chromosomes);
    List gen_original_cols(n_chromosomes);
    List risk_allele_vecs(n_chromosomes);
    IntegerVector n_case_risk_geno_vec(n_chromosomes);
    IntegerVector n_comp_risk_geno_vec(n_chromosomes);

    for (int i = 0; i < n_chromosomes; i++){

      List chrom_fitness_score_list_i = chrom_fitness_score_list[i];
      double fs = chrom_fitness_score_list_i["fitness_score"];
      NumericVector sdv = chrom_fitness_score_list_i["sum_dif_vecs"];
      CharacterVector rav = chrom_fitness_score_list_i["risk_set_alleles"];
      int n_case_risk_geno = chrom_fitness_score_list_i["n_case_risk_geno"];
      int n_comp_risk_geno = chrom_fitness_score_list_i["n_comp_risk_geno"];

      fitness_scores[i] = fs;
      sum_dif_vecs[i] = sdv;
      risk_allele_vecs[i] = rav;
      IntegerVector chrom_col_idx = chromosome_list[i];
      gen_original_cols[i] = chrom_col_idx;
      n_case_risk_geno_vec[i] = n_case_risk_geno;
      n_comp_risk_geno_vec[i] = n_comp_risk_geno;

    }

    List res = List::create(Named("chromosome_list") = chromosome_list,
                            Named("fitness_scores") = fitness_scores,
                            Named("sum_dif_vecs") = sum_dif_vecs,
                            Named("gen_original_cols") = gen_original_cols,
                            Named("risk_allele_vecs") = risk_allele_vecs,
                            Named("n_case_risk_geno_vec") = n_case_risk_geno_vec,
                            Named("n_comp_risk_geno_vec") = n_comp_risk_geno_vec);
    return(res);

  }

}

/////////////////////////////////////////////////
// function to pick out the top chromosome
// and generate a list of unique chroms
// from a given population
/////////////////////////////////////////////////

// [[Rcpp::export]]
List find_top_chrom(NumericVector fitness_scores, List chromosome_list, int chromosome_size){

  // find the top fitness score
  double max_fitness = max(fitness_scores);

  // determine how many chromosomes have the top score
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

  List res = List::create(Named("top_chromosome") = top_chromosome,
                          Named("max_fitness") = max_fitness,
                          Named("lower_chromosomes") = lower_chromosomes);
  return(res);

}

///////////////////////////////////////////
// function to initiate island populations
///////////////////////////////////////////

// [[Rcpp::export]]
List initiate_population(int n_candidate_snps, IntegerMatrix case_genetic_data, IntegerMatrix complement_genetic_data,
                         IntegerVector ld_block_vec, int n_chromosomes, int chromosome_size,
                         IntegerVector weight_lookup,
                         ListOf<IntegerMatrix> case_genetic_data_list,
                         ListOf<IntegerMatrix> complement_genetic_data_list,
                         IntegerVector exposure_levels,
                         IntegerVector exposure_risk_levels,
                         int n_different_snps_weight = 2, int n_both_one_weight = 1,
                         double recessive_ref_prop = 0.75, double recode_test_stat = 1.64,
                         int max_generations = 500, bool initial_sample_duplicates = false,
                         bool GxE = false, bool check_risk = true, int use_parents = 1){

  int n_possible_unique_combn = n_chromosomes * chromosome_size;
  if ((n_candidate_snps < n_possible_unique_combn) & !initial_sample_duplicates) {

    Rcout << "Not enough SNPs present to allow for no initial sample duplicate SNPs, now allowing initial sample duplicate snps.\n";
    initial_sample_duplicates = true;

  }

  // make initial chromosome list
  List chromosome_list(n_chromosomes);
  IntegerVector all_snps_idx = seq_len(n_candidate_snps);

  for (int i = 0; i < n_chromosomes; i++) {

    IntegerVector snp_idx = sample(all_snps_idx, chromosome_size, false);
    std::sort(snp_idx.begin(), snp_idx.end());
    if (!initial_sample_duplicates) {

      LogicalVector not_these = !in(all_snps_idx, snp_idx);
      all_snps_idx = all_snps_idx[not_these];

    }
    chromosome_list[i] = snp_idx;

  }


  // initiate storage objects
  int generation = 1;
  bool last_gens_equal = false;
  bool top_chrom_migrated = false;

  // compute fitness scores for the chromosome list
  List current_fitness_list = compute_population_fitness(case_genetic_data, complement_genetic_data, ld_block_vec, chromosome_list,
                                                         weight_lookup, case_genetic_data_list, complement_genetic_data_list,
                                                         exposure_levels, exposure_risk_levels, n_different_snps_weight,
                                                         n_both_one_weight, recessive_ref_prop, recode_test_stat, GxE,
                                                         check_risk, use_parents);

  // pick out the top and bottom scores
  chromosome_list = current_fitness_list["chromosome_list"];
  NumericVector fitness_scores = current_fitness_list["fitness_scores"];
  List gen_top_chrom_list = find_top_chrom(fitness_scores, chromosome_list, chromosome_size);
  IntegerVector top_chromosome = gen_top_chrom_list["top_chromosome"];
  List lower_chromosomes = gen_top_chrom_list["lower_chromosomes"];

  // counter for same top chrom
  int same_top_chrom = 1;

  List res = List::create(Named("current_top_generation_chromosome") = top_chromosome,
                          Named("current_lower_chromosomes") = lower_chromosomes,
                          Named("top_chrom_migrated") = top_chrom_migrated,
                          Named("current_fitness") = current_fitness_list,
                          Named("generation") = generation,
                          Named("last_gens_equal") = last_gens_equal,
                          Named("n_gens_same_top_chrom") = same_top_chrom);
  return(res);

}

///////////////////////////////////////////
// function to evolve island populations
//////////////////////////////////////////

List evolve_island(int n_migrations, IntegerMatrix case_genetic_data, IntegerMatrix complement_genetic_data,
                   IntegerVector ld_block_vec, int n_chromosomes, int chromosome_size, IntegerVector weight_lookup,
                   ListOf<IntegerMatrix> case_genetic_data_list,
                   ListOf<IntegerMatrix> complement_genetic_data_list,
                   IntegerVector exposure_levels,
                   IntegerVector exposure_risk_levels,
                   NumericVector snp_chisq, List population,
                   int n_different_snps_weight = 2, int n_both_one_weight = 1,
                   int migration_interval = 50, int gen_same_fitness = 50,
                   int max_generations = 500,
                   bool initial_sample_duplicates = false,
                   double crossover_prop = 0.8, double recessive_ref_prop = 0.75, double recode_test_stat = 1.64,
                   bool GxE = false, bool check_risk = true, int use_parents = 1){

  // initialize groups of candidate solutions if generation 1
  int generation = population["generation"];
  int generations = 0;
  if (generation == 1){

    generations = migration_interval;

  } else {

    generations = generation + migration_interval;

  }

  // grab input population data
  IntegerVector top_chromosome = population["current_top_generation_chromosome"];
  List lower_chromosomes = population["current_lower_chromosomes"];
  int n_gens_same_top_chrom = population["n_gens_same_top_chrom"];
  bool same_last_gens = population["last_gens_equal"];
  List current_fitness_list = population["current_fitness"];
  bool top_chrom_migrated = population["top_chrom_migrated"];

  //parse the current population fitness scores
  List chromosome_list = current_fitness_list["chromosome_list"];
  NumericVector fitness_scores = current_fitness_list["fitness_scores"];
  List sum_dif_vecs = current_fitness_list["sum_dif_vecs"];
  List gen_original_cols = current_fitness_list["gen_original_cols"];

  // iterate over generations
  bool all_converged = false;
  while ((generation < generations) & !all_converged) {

    // 1. Sample with replacement from the existing chromosomes, allowing the top
    // scoring chromosome to be sampled, but only sample from the unique chromosomes available
    LogicalVector sample_these = unique_chrom_list(chromosome_list, chromosome_size);
    IntegerVector chrom_list_idx = seq_len(n_chromosomes);
    IntegerVector which_sample_these = chrom_list_idx[sample_these];
    NumericVector sample_these_scores = fitness_scores[sample_these];
    IntegerVector sampled_lower_idx = sample(which_sample_these, lower_chromosomes.length(), true, sample_these_scores);
    List sampled_lower_chromosomes = chromosome_list[sampled_lower_idx - 1];
    List sampled_lower_dif_vecs = sum_dif_vecs[sampled_lower_idx - 1];
    NumericVector sampled_lower_fitness_scores = fitness_scores[sampled_lower_idx - 1];

    // 2. Determine whether each lower chromosome will be subject to mutation or crossing over

    // only allowing cross-overs between distinct chromosomes (i.e, if a chromosome was sampled twice,
    // it can't cross over with itself)
    IntegerVector unique_lower_idx = unique(sampled_lower_idx);
    LogicalVector cross_overs(unique_lower_idx.length(), false);
    int rounded_crosses = round(unique_lower_idx.length() * crossover_prop);
    IntegerVector possible_crosses = seq_len(unique_lower_idx.length());
    int n_crosses = 0;
    if ((rounded_crosses % 2 == 0) | (unique_lower_idx.length() == 1)){

      n_crosses = rounded_crosses;

    } else {

      n_crosses = rounded_crosses + 1;

    }
    IntegerVector cross_these = sample(possible_crosses, n_crosses, false);
    cross_overs[cross_these - 1] = true;

    // 3. Execute crossing over for the relevant chromosomes
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
          c1_c2_not_matching_snp_positions = sort_by_order(c1_c2_not_matching_snp_positions, abs(c2_not_matching_dif_vecs), 1);

          NumericVector c1_not_matching_dif_vecs = chrom1_dif_vecs[c2_c1_not_matching_snp_positions - 1];
          c2_c1_not_matching_snp_positions = sort_by_order(c2_c1_not_matching_snp_positions, abs(c1_not_matching_dif_vecs), 2);

        } else {

          NumericVector c2_not_matching_dif_vecs = chrom2_dif_vecs[c1_c2_not_matching_snp_positions - 1];
          c1_c2_not_matching_snp_positions = sort_by_order(c1_c2_not_matching_snp_positions, abs(c2_not_matching_dif_vecs), 2);

          NumericVector c1_not_matching_dif_vecs = chrom1_dif_vecs[c2_c1_not_matching_snp_positions - 1];
          c2_c1_not_matching_snp_positions = sort_by_order(c2_c1_not_matching_snp_positions, abs(c1_not_matching_dif_vecs), 1);

        }

        // order the chromosomes, first by overlapping snps and then non-overlapping
        IntegerVector c2_order = concat(c1_c2_matching_snp_positions, c1_c2_not_matching_snp_positions);
        chrom2 = chrom2[c2_order - 1];

        IntegerVector c1_order = concat(c2_c1_matching_snp_positions, c2_c1_not_matching_snp_positions);
        chrom1 = chrom1[c1_order - 1];

        // determine how many snps could be crossed over and also make sure we don't simply swap chromosomes
        IntegerVector possible_cut_points = chrom_positions[chrom1 != chrom2];
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

    // 4. mutate chromosomes
    IntegerVector tmp_lower_idx = seq_len(sampled_lower_chromosomes.length());
    IntegerVector mutation_positions = setdiff(tmp_lower_idx, cross_over_positions);
    if (mutation_positions.length() > 0){

      IntegerVector candiate_snp_idx = seq_along(snp_chisq);

      IntegerVector snps_for_mutation = sample(candiate_snp_idx, candiate_snp_idx.length(),
                                               true, snp_chisq);

      int ulen = unique(snps_for_mutation).length();
      IntegerVector n_possible_mutations = seq_len(chromosome_size);
      for (int l = 0; l < mutation_positions.length(); l++){

        int mutation_position = mutation_positions[l];

        // grab chromosome and its difference vector
        IntegerVector target_chrom = sampled_lower_chromosomes[mutation_position - 1];
        NumericVector target_dif_vec = sampled_lower_dif_vecs[mutation_position - 1];

        // sort the chromosome elements from lowest absolute difference vector to highest
        target_chrom = sort_by_order(target_chrom, abs(target_dif_vec), 1);

        // determine which snps to mutate
        int sampled_n_mutations = sample(n_possible_mutations, 1)[0];

        // note that we account for the very unlikely case where the pool
        // of unique possible mutations is smaller than the sampled number of mutations
        int total_mutations = scalar_min(sampled_n_mutations, ulen);

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
        if (total_mutations > 0){

          IntegerVector mutate_here = seq(0, total_mutations - 1);
          target_chrom[mutate_here] = mutated_snps;
          sampled_lower_chromosomes[mutation_position - 1] = target_chrom;

        }

      }

    }

    // 5. Combined into new population (i.e., the final collection of chromosomes for the next generation)
    std::sort(top_chromosome.begin(), top_chromosome.end());
    chromosome_list[0] = top_chromosome;
    for (int i = 1; i < chromosome_list.size(); i++){

      IntegerVector sampled_lower_chromosomes_i = sampled_lower_chromosomes[i - 1];
      std::sort(sampled_lower_chromosomes_i.begin(), sampled_lower_chromosomes_i.end());
      chromosome_list[i] = sampled_lower_chromosomes_i;

    }

    // 6. Compute new population fitness
    current_fitness_list = compute_population_fitness(case_genetic_data, complement_genetic_data, ld_block_vec, chromosome_list,
                                                      weight_lookup, case_genetic_data_list, complement_genetic_data_list,
                                                      exposure_levels, exposure_risk_levels, n_different_snps_weight,
                                                      n_both_one_weight, recessive_ref_prop, recode_test_stat, GxE,
                                                      check_risk, use_parents);

    //7. Increment iterators
    generation += 1;

    // 8. pick out relevant pieces of the output
    NumericVector fitness_scores_tmp = current_fitness_list["fitness_scores"];
    fitness_scores = fitness_scores_tmp;
    sum_dif_vecs = current_fitness_list["sum_dif_vecs"];
    gen_original_cols = current_fitness_list["gen_original_cols"];

    // 9. identify the top scoring candidate solution(s) and fitness score
    List gen_top_chrom_list = find_top_chrom(fitness_scores, chromosome_list, chromosome_size);
    IntegerVector new_top_chromosome = gen_top_chrom_list["top_chromosome"];
    lower_chromosomes = gen_top_chrom_list["lower_chromosomes"];

    // check if the top chromosome is the same as the previous generation
    if (top_chrom_migrated){

      n_gens_same_top_chrom = 1;
      top_chrom_migrated = false;

    } else {

      if (is_true(all(new_top_chromosome == top_chromosome))){

        n_gens_same_top_chrom += 1;

      } else {

        n_gens_same_top_chrom = 1;

      }

    }

    // update the current top generation chromosome
    top_chromosome = new_top_chromosome;

    // check to see if we can stop iterating over generations
    if ((generation >= gen_same_fitness)) {

      same_last_gens = n_gens_same_top_chrom >= gen_same_fitness;

      if (n_migrations == 0){

        if (same_last_gens){

          all_converged = true;

        }

      }

    }

  }

  // if the algorithm hasn't hit the max number of generations prepare for potential migration
  if ((generation < max_generations) & (n_migrations > 0)){

    // order by fitness score
    IntegerVector chrom_list_order = sort_by_order(seq_len(chromosome_list.size()), fitness_scores, 2);
    chromosome_list = chromosome_list[chrom_list_order - 1];
    NumericVector fitness_scores_tmp = fitness_scores[chrom_list_order - 1];
    fitness_scores = fitness_scores_tmp;
    sum_dif_vecs = sum_dif_vecs[chrom_list_order - 1];
    gen_original_cols = gen_original_cols[chrom_list_order - 1];

    // store in population
    if (GxE){

      List pop_current_fitness_list = population["current_fitness"];
      pop_current_fitness_list["chromosome_list"] = chromosome_list;
      pop_current_fitness_list["fitness_scores"] = fitness_scores;
      pop_current_fitness_list["sum_dif_vecs"] = sum_dif_vecs;
      pop_current_fitness_list["gen_original_cols"] = gen_original_cols;
      if (use_parents != 2){

        List exposure_info = current_fitness_list["exposure_level_info"];
        exposure_info = exposure_info[chrom_list_order - 1];
        pop_current_fitness_list["exposure_level_info"] = exposure_info;

      }

      population["last_gens_equal"] = same_last_gens;
      population["current_top_generation_chromosome"] = top_chromosome;
      population["current_lower_chromosomes"] = lower_chromosomes;
      population["n_gens_same_top_chrom"] = n_gens_same_top_chrom;
      population["generation"] = generation;
      population["top_chrom_migrated"] = top_chrom_migrated;

      // return the population
      return(population);

    } else {

      List risk_allele_vecs = current_fitness_list["risk_allele_vecs"];
      risk_allele_vecs = risk_allele_vecs[chrom_list_order - 1];

      IntegerVector n_case_risk_geno_vec_tmp = current_fitness_list["n_case_risk_geno_vec"];
      IntegerVector n_case_risk_geno_vec = n_case_risk_geno_vec_tmp[chrom_list_order - 1];

      IntegerVector n_comp_risk_geno_vec_tmp = current_fitness_list["n_comp_risk_geno_vec"];
      IntegerVector n_comp_risk_geno_vec = n_comp_risk_geno_vec_tmp[chrom_list_order - 1];

      List pop_current_fitness_list = population["current_fitness"];
      pop_current_fitness_list["chromosome_list"] = chromosome_list;
      pop_current_fitness_list["fitness_scores"] = fitness_scores;
      pop_current_fitness_list["sum_dif_vecs"] = sum_dif_vecs;
      pop_current_fitness_list["gen_original_cols"] = gen_original_cols;
      pop_current_fitness_list["risk_allele_vecs"] = risk_allele_vecs;
      pop_current_fitness_list["n_case_risk_geno_vec"] = n_case_risk_geno_vec;
      pop_current_fitness_list["n_comp_risk_geno_vec"] = n_comp_risk_geno_vec;

      population["last_gens_equal"] = same_last_gens;
      population["current_top_generation_chromosome"] = top_chromosome;
      population["current_lower_chromosomes"] = lower_chromosomes;
      population["n_gens_same_top_chrom"] = n_gens_same_top_chrom;
      population["generation"] = generation;
      population["top_chrom_migrated"] = top_chrom_migrated;

      // return the population
      return(population);

    }

  // otherwise just store the current generation and return the population
  } else {

    if (GxE){

      List pop_current_fitness_list = population["current_fitness"];
      pop_current_fitness_list["chromosome_list"] = chromosome_list;
      pop_current_fitness_list["fitness_scores"] = fitness_scores;
      pop_current_fitness_list["sum_dif_vecs"] = sum_dif_vecs;
      pop_current_fitness_list["gen_original_cols"] = gen_original_cols;
      if (use_parents != 2){

        List exposure_info = current_fitness_list["exposure_level_info"];
        pop_current_fitness_list["exposure_level_info"] = exposure_info;

      }

      population["last_gens_equal"] = same_last_gens;
      population["current_top_generation_chromosome"] = top_chromosome;
      population["current_lower_chromosomes"] = lower_chromosomes;
      population["n_gens_same_top_chrom"] = n_gens_same_top_chrom;
      population["generation"] = generation;
      population["top_chrom_migrated"] = top_chrom_migrated;

      // return the population
      return(population);

    } else {

      List risk_allele_vecs = current_fitness_list["risk_allele_vecs"];
      IntegerVector n_case_risk_geno_vec = current_fitness_list["n_case_risk_geno_vec"];
      IntegerVector n_comp_risk_geno_vec = current_fitness_list["n_comp_risk_geno_vec"];

      List pop_current_fitness_list = population["current_fitness"];
      pop_current_fitness_list["chromosome_list"] = chromosome_list;
      pop_current_fitness_list["fitness_scores"] = fitness_scores;
      pop_current_fitness_list["sum_dif_vecs"] = sum_dif_vecs;
      pop_current_fitness_list["gen_original_cols"] = gen_original_cols;
      pop_current_fitness_list["risk_allele_vecs"] = risk_allele_vecs;
      pop_current_fitness_list["n_case_risk_geno_vec"] = n_case_risk_geno_vec;
      pop_current_fitness_list["n_comp_risk_geno_vec"] = n_comp_risk_geno_vec;

      population["last_gens_equal"] = same_last_gens;
      population["current_top_generation_chromosome"] = top_chromosome;
      population["current_lower_chromosomes"] = lower_chromosomes;
      population["n_gens_same_top_chrom"] = n_gens_same_top_chrom;
      population["generation"] = generation;
      population["top_chrom_migrated"] = top_chrom_migrated;

      return(population);

    }

  }
}

////////////////////////////////////////////
// function to check island convergence
//////////////////////////////////////////

// [[Rcpp::export]]
bool check_convergence(int island_cluster_size, List island_populations){

  // check convergence
  int n_converged = 0;
  bool all_converged = false;
  for (int i = 0; i < island_cluster_size; i++){

    List island_i = island_populations[i];
    bool last_gens_equal_i = island_i["last_gens_equal"];
    if (last_gens_equal_i){

      n_converged += 1;

    }

  }
  if (n_converged == island_cluster_size){

    all_converged = true;

  }
  return(all_converged);

}

////////////////////////////////////////////
// function to check maximum generations
//////////////////////////////////////////

// [[Rcpp::export]]
bool check_max_gens(List island_populations, int max_generations){

  // just pick one island and check the generation
  List this_island = island_populations[0];
  int generation = this_island["generation"];
  if (generation == max_generations){

    return(true);

  } else {

    return(false);

  }
}

////////////////////////////////////////////////////////////
// Function to put all the pieces together and run GADGETS
///////////////////////////////////////////////////////////

// [[Rcpp::export]]
List run_GADGETS(int island_cluster_size, int n_migrations,
                 IntegerVector ld_block_vec, int n_chromosomes, int chromosome_size, IntegerVector weight_lookup,
                 NumericVector snp_chisq, IntegerMatrix case_genetic_data_in, IntegerMatrix complement_genetic_data_in,
                 Nullable<IntegerVector> exposure_levels_in = R_NilValue,
                 Nullable<IntegerVector> exposure_risk_levels_in = R_NilValue,
                 Nullable<IntegerVector> exposure_in = R_NilValue,
                 int n_different_snps_weight = 2, int n_both_one_weight = 1, int migration_interval = 50,
                 int gen_same_fitness = 50, int max_generations = 500,
                 bool initial_sample_duplicates = false, double crossover_prop = 0.8,
                 double recessive_ref_prop = 0.75, double recode_test_stat = 1.64, int use_parents = 1){

  // instantiate input objects
  IntegerMatrix case_genetic_data;
  IntegerMatrix complement_genetic_data;

  ListOf<IntegerMatrix> case_genetic_data_list;
  ListOf<IntegerMatrix> complement_genetic_data_list;
  IntegerVector exposure_levels;
  IntegerVector exposure_risk_levels;
  IntegerVector exposure;
  bool GxE = false;
  bool check_risk = true;
  int n_candidate_snps = 0;

  if (exposure_in.isNotNull()){

    // grab exposures
    exposure_levels = exposure_levels_in;
    exposure_risk_levels = exposure_risk_levels_in;
    exposure = exposure_in;

    // grab input data
    IntegerMatrix case_genetic_data_tmp = case_genetic_data_in;
    IntegerMatrix complement_genetic_data_tmp = complement_genetic_data_in;

    // split the input data into a list by exposure
    case_genetic_data_list = split_int_mat(case_genetic_data_tmp, exposure, exposure_levels);
    complement_genetic_data_list = split_int_mat(complement_genetic_data_tmp, exposure, exposure_levels);
    GxE = true;

    // decide if risk model needs to be checked
    IntegerVector unique_exposure_risk_levels = unique(exposure_risk_levels);
    check_risk = unique_exposure_risk_levels.length() > 1;
    n_candidate_snps = case_genetic_data_tmp.ncol();

  } else {

    case_genetic_data = case_genetic_data_in;
    complement_genetic_data = complement_genetic_data_in;
    n_candidate_snps = case_genetic_data.ncol();

  }

  // go through first round of island evolution
  List island_populations(island_cluster_size);

  for (int i = 0; i < island_cluster_size; i++){

    List island_population_i = initiate_population(n_candidate_snps, case_genetic_data, complement_genetic_data,
                                                   ld_block_vec, n_chromosomes, chromosome_size,
                                                   weight_lookup, case_genetic_data_list,
                                                   complement_genetic_data_list,
                                                   exposure_levels, exposure_risk_levels,
                                                   n_different_snps_weight, n_both_one_weight,
                                                   recessive_ref_prop, recode_test_stat,
                                                   max_generations, initial_sample_duplicates, GxE,
                                                   check_risk, use_parents);

    island_populations[i] = evolve_island(n_migrations, case_genetic_data, complement_genetic_data,
                                          ld_block_vec, n_chromosomes, chromosome_size, weight_lookup,
                                          case_genetic_data_list, complement_genetic_data_list,
                                          exposure_levels, exposure_risk_levels, snp_chisq,
                                          island_population_i, n_different_snps_weight, n_both_one_weight,
                                          migration_interval, gen_same_fitness, max_generations,
                                          initial_sample_duplicates, crossover_prop, recessive_ref_prop, recode_test_stat,
                                          GxE, check_risk, use_parents);

  }

  // check convergence
  bool all_converged = check_convergence(island_cluster_size, island_populations);

  // check max gens
  bool max_gens = check_max_gens(island_populations, max_generations);

  // if no convergence/max gens, migrate chromosomes and continue evolution until
  // convergence/max gens
  while (!all_converged & !max_gens){

    IntegerVector donor_idx = seq_len(n_migrations) - 1;
    int start_here = n_chromosomes - n_migrations;
    IntegerVector receiver_idx = seq(start_here, n_chromosomes - 1);

    // migrate top chromosomes among islands
    for (int i = 0; i < island_cluster_size; i++){

      List first_island = island_populations[i];
      List first_island_current = first_island["current_fitness"];

      List second_island;
      if (i == 0){

        second_island = island_populations[island_cluster_size - 1];

      } else {

        second_island = island_populations[i-1];

      }
      List second_island_current = second_island["current_fitness"];

      // migrate chromosomes
      List second_island_chroms = second_island_current["chromosome_list"];
      List first_island_chroms = first_island_current["chromosome_list"];
      List chrom_migrations = second_island_chroms[donor_idx];
      first_island_chroms[receiver_idx] = chrom_migrations;

      // migrate corresponding fitness scores
      NumericVector first_island_fitness_scores = first_island_current["fitness_scores"];
      NumericVector second_island_fitness_scores = second_island_current["fitness_scores"];
      NumericVector fitness_migrations = second_island_fitness_scores[donor_idx];
      first_island_fitness_scores[receiver_idx] = fitness_migrations;

      // migrate corresponding dif vecs
      List second_island_dif_vecs = second_island_current["sum_dif_vecs"];
      List first_island_dif_vecs = first_island_current["sum_dif_vecs"];
      List dif_vec_migrations = second_island_dif_vecs[donor_idx];
      first_island_dif_vecs[receiver_idx] = dif_vec_migrations;

      // migrate original chromosome numbers
      List second_island_orig_chroms = second_island_current["gen_original_cols"];
      List first_island_orig_chroms = first_island_current["gen_original_cols"];
      List original_chrom_migrations = second_island_orig_chroms[donor_idx];
      first_island_orig_chroms[receiver_idx] = original_chrom_migrations;

      if (GxE){

        if (use_parents == 1){

          // migrate exposure info
          List first_island_exposure_info = first_island_current["exposure_level_info"];
          List second_island_exposure_info = second_island_current["exposure_level_info"];
          List exposure_info_migrations = second_island_exposure_info[donor_idx];
          first_island_exposure_info[receiver_idx] = exposure_info_migrations;

        }

      } else {

        // migrate risk allele vecs
        List second_island_risk_alleles = second_island_current["risk_allele_vecs"];
        List first_island_risk_alleles = first_island_current["risk_allele_vecs"];
        List risk_allele_migrations = second_island_risk_alleles[donor_idx];
        first_island_risk_alleles[receiver_idx] = risk_allele_migrations;

        // migrate n case high risk geno
        IntegerVector first_island_case_geno = first_island_current["n_case_risk_geno_vec"];
        IntegerVector second_island_case_geno = second_island_current["n_case_risk_geno_vec"];
        IntegerVector case_geno_migrations = second_island_case_geno[donor_idx];
        first_island_case_geno[receiver_idx] = case_geno_migrations;

        // migrate n comp high risk geno
        IntegerVector first_island_comp_geno = first_island_current["n_comp_risk_geno_vec"];
        IntegerVector second_island_comp_geno = second_island_current["n_comp_risk_geno_vec"];
        IntegerVector comp_geno_migrations = second_island_comp_geno[donor_idx];
        first_island_comp_geno[receiver_idx] = comp_geno_migrations;

      }

      // update the top/lower chromosomes, and the n_same_last_gens
      List gen_top_chrom_list = find_top_chrom(first_island_fitness_scores, first_island_chroms, chromosome_size);
      IntegerVector new_top_chromosome = gen_top_chrom_list["top_chromosome"];
      List lower_chromosomes = gen_top_chrom_list["lower_chromosomes"];
      first_island["current_lower_chromosomes"] = lower_chromosomes;
      IntegerVector top_chromosome = first_island["current_top_generation_chromosome"];

      if (!is_true(all(new_top_chromosome == top_chromosome))){

        first_island["n_gens_same_top_chrom"] = 1;
        first_island["last_gens_equal"] = false;
        first_island["current_top_generation_chromosome"] = new_top_chromosome;
        first_island["top_chrom_migrated"]  = true;

      }

    }

    // evolve islands using the new populations

    for (int i = 0; i < island_cluster_size; i++){

      List island_population_i = island_populations[i];
      island_populations[i] = evolve_island(n_migrations, case_genetic_data, complement_genetic_data,
                                            ld_block_vec, n_chromosomes, chromosome_size, weight_lookup,
                                            case_genetic_data_list, complement_genetic_data_list,
                                            exposure_levels, exposure_risk_levels, snp_chisq,
                                            island_population_i, n_different_snps_weight, n_both_one_weight,
                                            migration_interval, gen_same_fitness, max_generations,
                                            initial_sample_duplicates, crossover_prop, recessive_ref_prop, recode_test_stat,
                                            GxE, check_risk, use_parents);

    }

    // check convergence
    all_converged = check_convergence(island_cluster_size, island_populations);

    // check max gens
    max_gens = check_max_gens(island_populations, max_generations);

  }

  //return evolved island list
  return(island_populations);

}

//////////////////////////////////////////////////////////////////////////////
// Function shuffle rows of input case/comp data and recompute fitness score
/////////////////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
double epistasis_test_permute(arma::mat case_inf, arma::mat comp_inf, IntegerVector target_snps_ld_blocks,
                              IntegerVector uni_ld_blocks, int n_families,
                              IntegerVector weight_lookup, int n_different_snps_weight = 2, int n_both_one_weight = 1,
                              double recessive_ref_prop = 0.75, double recode_test_stat = 1.64){

  // make copies of the input data
  uint nrows = case_inf.n_rows;
  uint ncols = case_inf.n_cols;
  arma::mat case_permuted(nrows, ncols);
  arma::mat comp_permuted(nrows, ncols);
  IntegerVector target_snps_idx = seq_along(target_snps_ld_blocks);

  // loop over SNP LD blocks and shuffle rows
  for (int i = 0; i < uni_ld_blocks.length(); i++){

    IntegerVector family_idx = seq_len(n_families);
    arma::uvec row_order = arma::randperm(n_families, n_families);
    int ld_block_val = uni_ld_blocks[i];
    LogicalVector ld_block_snps_pos = target_snps_ld_blocks == ld_block_val;
    IntegerVector these_snps_rcpp = target_snps_idx[ld_block_snps_pos];
    arma::vec these_snps = as<arma::vec>(these_snps_rcpp);

    for (arma::uword j = 0; j < row_order.size(); j++){

      arma::uword in_row = row_order(j);

      for (arma::uword k = 0; k < these_snps.size(); k++){

        arma::uword snp_col = these_snps(k) - 1;
        int case_val_k = case_inf(in_row, snp_col);
        case_permuted(j, snp_col) = case_val_k;
        int comp_val_k = comp_inf(in_row, snp_col);
        comp_permuted(j, snp_col) = comp_val_k;

      }

    }

  }

  // convert back to regular rcpp
  IntegerMatrix case_rcpp = wrap(case_permuted);
  IntegerMatrix comp_rcpp = wrap(comp_permuted);

  // compute fitness score for permuted dataset
  IntegerVector target_snps = seq_len(case_rcpp.ncol());
  List fitness_list = chrom_fitness_score(case_rcpp, comp_rcpp, target_snps,
                                          target_snps_ld_blocks, weight_lookup,
                                          n_different_snps_weight, n_both_one_weight,
                                          recessive_ref_prop, recode_test_stat, false);

  double fitness = fitness_list["fitness_score"];
  return(fitness);

}

/////////////////////////////////////////////////////////////////
// function to generate n pemutes and compute the fitness score,
// to generate the null distribution for the epistasis test
/////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
NumericVector epistasis_test_null_scores(int n_permutes, arma::mat case_inf, arma::mat comp_inf,
                                         IntegerVector target_snps_ld_blocks, IntegerVector uni_ld_blocks,
                                         int n_families, IntegerVector weight_lookup,
                                         int n_different_snps_weight = 2, int n_both_one_weight = 1,
                                         double recessive_ref_prop = 0.75, double recode_test_stat = 1.64){

  // loop over number of permutes and output vector of null fitness scores
  NumericVector res(n_permutes);
  for (int i = 0; i < n_permutes; i++){

    double permute_score_i = epistasis_test_permute(case_inf, comp_inf, target_snps_ld_blocks, uni_ld_blocks, n_families,
                                                    weight_lookup, n_different_snps_weight, n_both_one_weight,
                                                    recessive_ref_prop, recode_test_stat);
    res[i] = permute_score_i;

  }
  return(res);

}

//////////////////////////////////
// function to run epistasis test
/////////////////////////////////

// [[Rcpp::export]]
List epistasis_test(IntegerVector snp_cols, List preprocessed_list, int n_permutes = 10000,
                    int n_different_snps_weight = 2, int n_both_one_weight = 1, int weight_function_int = 2,
                    double recessive_ref_prop = 0.75, double recode_test_stat = 1.64, bool warn = true){

  // pick out target columns in the preprocessed data
  IntegerVector target_snps = snp_cols;

  // grab LD info
  IntegerVector ld_block_vec = preprocessed_list["ld.block.vec"];
  IntegerVector target_snps_ld_blocks = get_target_snps_ld_blocks(snp_cols, ld_block_vec);
  IntegerVector uni_ld_blocks = unique(target_snps_ld_blocks);

  // require more than 1 ld block
  int n_ld_blocks = uni_ld_blocks.length();
  if (n_ld_blocks == 1){

    if (warn){

      Rcout << "All chromosome SNPs in linkage, returning NA for p-value \n";

    }
    List res = List::create(Named("pval") = NA_REAL,
                            Named("obs_fitness_score") = NA_REAL,
                            Named("perm_fitness_scores") = NumericVector::get_na());
    return(res);

  }

  // compute weight lookup table
  int max_weight = n_different_snps_weight;
  if (n_both_one_weight > n_different_snps_weight){

    max_weight = n_both_one_weight;

  }

  int max_sum = max_weight*snp_cols.length();
  IntegerVector exponents = seq_len(max_sum);
  IntegerVector weight_lookup(max_sum);
  for (int i = 0; i < max_sum; i++){

    int exponent_i = exponents[i];
    int lookup_i = pow(weight_function_int, exponent_i);
    weight_lookup[i] = lookup_i;

  }

  // get input genetic data
  IntegerMatrix case_genetic_data_in = preprocessed_list["case.genetic.data"];
  IntegerMatrix complement_genetic_data_in = preprocessed_list["complement.genetic.data"];
  IntegerMatrix case_genetic_data = subset_matrix_cols(case_genetic_data_in, snp_cols);
  IntegerMatrix complement_genetic_data = subset_matrix_cols(complement_genetic_data_in, snp_cols);

  // compute fitness score for observed data
  IntegerVector chrom_snps = seq_len(case_genetic_data.ncol());
  List obs_fitness_list = chrom_fitness_score(case_genetic_data, complement_genetic_data, chrom_snps,
                                              target_snps_ld_blocks, weight_lookup, n_different_snps_weight,
                                              n_both_one_weight, recessive_ref_prop, recode_test_stat, true);
  double obs_fitness_score = obs_fitness_list["fitness_score"];

  // restrict to informative families
  IntegerVector informative_families = obs_fitness_list["inf_families"];
  IntegerMatrix case_inf = subset_matrix_rows(case_genetic_data, informative_families);
  arma::mat case_inf_arma = as<arma::mat>(case_inf);
  IntegerMatrix comp_inf = subset_matrix_rows(complement_genetic_data, informative_families);
  arma::mat comp_inf_arma = as<arma::mat>(comp_inf);
  int n_families = informative_families.length();

  // loop over permuted datasets and compute fitness scores
  NumericVector perm_fitness_scores = epistasis_test_null_scores(n_permutes, case_inf_arma, comp_inf_arma, target_snps_ld_blocks,
                                                                 uni_ld_blocks, n_families, weight_lookup,
                                                                 n_different_snps_weight, n_both_one_weight,
                                                                 recessive_ref_prop, recode_test_stat);

  // compute p-value
  double N = n_permutes + 1;
  double B = sum(perm_fitness_scores >= obs_fitness_score);
  double pval = (B + 1)/N;

  //return result
  List res = List::create(Named("pval") = pval,
                          Named("obs_fitness_score") = obs_fitness_score,
                          Named("perm_fitness_scores") = perm_fitness_scores);
  return(res);

}

/////////////////////////////////////////////
// function to run interaction permutation test
////////////////////////////////////////////

// [[Rcpp::export]]
List GxE_test(IntegerVector snp_cols, List preprocessed_list, int n_permutes = 10000,
              int n_different_snps_weight = 2, int n_both_one_weight = 1, int weight_function_int = 2,
              double recessive_ref_prop = 0.75, double recode_test_stat = 1.64, int use_parents = 1){

  // pick out target columns in the preprocessed data
  IntegerVector target_snps = snp_cols;

  // grab LD mat
  IntegerVector ld_block_vec = preprocessed_list["ld.block.vec"];
  IntegerVector target_snps_ld_blocks = get_target_snps_ld_blocks(snp_cols, ld_block_vec);

  // Grab exposure info
  IntegerVector exposure = preprocessed_list["exposure"];
  IntegerVector exposure_levels = preprocessed_list["exposure.levels"];
  IntegerVector exposure_risk_levels = preprocessed_list["exposure.risk.levels"];

  // compute weight lookup table
  int max_weight = n_different_snps_weight;
  if (n_both_one_weight > n_different_snps_weight){

    max_weight = n_both_one_weight;

  }

  int max_sum = max_weight*snp_cols.length();
  IntegerVector exponents = seq_len(max_sum);
  IntegerVector weight_lookup(max_sum);
  for (int i = 0; i < max_sum; i++){

    int exponent_i = exponents[i];
    int lookup_i = pow(weight_function_int, exponent_i);
    weight_lookup[i] = lookup_i;

  }

  // get input genetic data
  IntegerMatrix case_genetic_data = preprocessed_list["case.genetic.data"];
  IntegerMatrix complement_genetic_data = preprocessed_list["complement.genetic.data"];

  // split by exposure
  ListOf<IntegerMatrix> case_genetic_data_list = split_int_mat(case_genetic_data, exposure, exposure_levels);
  ListOf<IntegerMatrix> complement_genetic_data_list = split_int_mat(complement_genetic_data, exposure, exposure_levels);

  int chrom_size = target_snps.length();
  bool check_risk = true;

  // see if there is more than 1 risk level
  IntegerVector unique_risk_levels = unique(exposure_risk_levels);
  if (unique_risk_levels.length() == 1){

    check_risk = false;

  }

  // compute fitness score for observed data
  IntegerVector chrom_snps = seq_len(chrom_size);
  List obs_fitness_list =  GxE_fitness_score(case_genetic_data_list, complement_genetic_data_list,
                                             chrom_snps, target_snps_ld_blocks, weight_lookup, exposure_levels,
                                             exposure_risk_levels, n_different_snps_weight, n_both_one_weight,
                                             recessive_ref_prop, recode_test_stat, true, check_risk, use_parents);

  double obs_fitness_score = obs_fitness_list["fitness_score"];

  // loop over number of permutations, randomizing exposure and recomputing fitness score
  NumericVector perm_fitness_scores(n_permutes);
  for (int i = 0; i < n_permutes; i++){

    // shuffle the observed exposures
    IntegerVector exposure_perm = sample(exposure, exposure.length(), false);
    IntegerVector exposure_levels_perm = sort_unique(exposure_perm);

    // make sure the risk levels vector matches with the shuffled exposures
    IntegerVector exposure_risk_levels_perm = exposure_risk_levels;
    if (check_risk){

      IntegerVector perm_matches = match(exposure_levels, exposure_levels_perm);
      exposure_risk_levels_perm = exposure_risk_levels[perm_matches - 1];

    }

    // split by exposure
    ListOf<IntegerMatrix> case_genetic_data_list_perm = split_int_mat(case_genetic_data, exposure_perm, exposure_levels_perm);
    ListOf<IntegerMatrix> complement_genetic_data_list_perm = split_int_mat(complement_genetic_data, exposure_perm, exposure_levels_perm);

    // compute fitness
    List perm_fitness_list =  GxE_fitness_score(case_genetic_data_list_perm, complement_genetic_data_list_perm,
                                                chrom_snps, target_snps_ld_blocks, weight_lookup, exposure_levels_perm,
                                                exposure_risk_levels_perm, n_different_snps_weight, n_both_one_weight,
                                                recessive_ref_prop, recode_test_stat, true, check_risk, use_parents);
    double perm_fitness = perm_fitness_list["fitness_score"];
    perm_fitness_scores[i] = perm_fitness;

  }

  // compute p-value
  double N = n_permutes + 1;
  double B = sum(perm_fitness_scores >= obs_fitness_score);
  double pval = (B + 1)/N;

  //return result
  List res = List::create(Named("pval") = pval,
                          Named("obs_fitness_score") = obs_fitness_score,
                          Named("perm_fitness_scores") = perm_fitness_scores);
  return(res);

}

/////////////////////////////////////////////////////////
// function to loop over a list of chromosomes and
// compute epistasis or GxE interaction p-value, and reuturn the -2log.
// This is used in the graphical scoring procedure.
/////////////////////////////////////////////////////////


// [[Rcpp::export]]
NumericVector n2log_epistasis_pvals(ListOf<IntegerVector> chromosome_list, List in_data_list, int n_permutes = 10000,
                                    int n_different_snps_weight = 2, int n_both_one_weight = 1, int weight_function_int = 2,
                                    double recessive_ref_prop = 0.75, double recode_test_stat = 1.64, bool GxE = false,
                                    int use_parents = 1){

  NumericVector n2log_epi_pvals(chromosome_list.size());
  double N = n_permutes + 1;
  for (int i = 0; i < chromosome_list.size(); i++){

    IntegerVector chromosome = chromosome_list[i];
    List perm_res;

    if (GxE){

      perm_res = GxE_test(chromosome, in_data_list, n_permutes, n_different_snps_weight, n_both_one_weight,
                          weight_function_int, recessive_ref_prop, recode_test_stat, use_parents);

    } else {

      perm_res = epistasis_test(chromosome, in_data_list, n_permutes,
                                n_different_snps_weight, n_both_one_weight, weight_function_int,
                                recessive_ref_prop, recode_test_stat, false);

    }

    NumericVector pval_vec = perm_res["pval"];
    double pval = pval_vec[0];
    LogicalVector pval_na_vec = NumericVector::is_na(pval);
    bool pval_na = pval_na_vec[0];
    if (pval_na){

      pval = 0.5;

    }

    if (pval == 1){

      pval = 1 - 1/N;

    }
    pval = -2*log(pval);
    n2log_epi_pvals[i] = pval;

  }
  return(n2log_epi_pvals);

}
