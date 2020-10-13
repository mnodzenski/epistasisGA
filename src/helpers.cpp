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

      out_mat(i, j) = in_mat(i, j) > comp_val[j];

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
LogicalMatrix comp_mat_ne(IntegerMatrix in_mat, int comp_val){

  int nrows = in_mat.nrow();
  int ncols = in_mat.ncol();
  LogicalMatrix out_mat(nrows, ncols);
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
arma::mat mult_rows_by_scalars(IntegerMatrix in_mat, IntegerVector scalars){

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


