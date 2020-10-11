#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

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


List chrom_fitness_score(IntegerMatrix case_genetic_data, IntegerMatrix complement_genetic_data, LogicalMatrix case_comp_differences,
                         IntegerVector target_snps, IntegerMatrix cases_minus_complements, LogicalMatrix both_one_mat,
                         LogicalMatrix block_ld_mat, IntegerVector weight_lookup,
                         int n_different_snps_weight = 2, int n_both_one_weight = 1,
                         int n_case_high_risk_thresh = 20, int outlier_sd = 2.5, bool epi_test = false) {

  // pick out the differences for the target snps
  IntegerMatrix case_comp_diff = subset_matrix_cols(case_comp_diff, target_snps);
  case_comp_diff = transpose(case_comp_diff);

  // determine whether families are informative for the set of target_snps
  IntegerVector total_different_snps = colSums(case_comp_diff);
  LogicalVector informative_families_l = total_different_snps != 0;
  IntegerVector family_idx = seq_len(total_different_snps.length());
  IntegerVector informative_families = family_idx[informative_families_l];
  int n_informative_families = sum(informative_families_l);
  IntegerMatrix case_inf = transpose(subset_matrix(case_genetic_data, informative_families, target_snps));
  IntegerMatrix comp_inf = transpose(subset_matrix(comp_genetic_data, informative_families, target_snps));
  IntegerMatrix cases_minus_complements_inf = transpose(subset_matrix(cases_minus_complements, informative_families, target_snps));

  // compute weights
  LogicalMatrix both_one_inf = transpose(subset_matrix(both_one_mat, informative_families, target_snps));
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
  int n_target = target_snps.length();
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

  }
  if (n_pos > 0) {

    IntegerVector pos_risk_int = idx_vec[pos_risk];
    IntegerMatrix case_inf_pos subset_int_matrix_rows(case_inf, pos_risk_int);
    IntegerMatrix comp_inf_pos subset_int_matrix_rows(comp_inf, pos_risk_int);
    LogicalMatrix p1 = comp_mat_gt(case_inf_pos, 0);
    LogicalMatrix p2 = comp_mat_gt(comp_inf_pos, 0);

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

  // pick out misclassifications via outlier detection, indicating recessive mode of inheritance

  // initialize vector of pattern of inheritance
  CharacterVector risk_set_alleles(target_snps.length(), "1+");

  // only applies if we have at least 20 high risk case
  if (n.case.high.risk > n.case.high.risk.thresh) {

    NumericVector case_high_risk_means = rowMeans(case_high_inf);
    NumericVector case_high_risk_sd = rowSds(case_high_inf);
    double high_outlier_thresh = case_high_risk_means + outlier_sd * case_high_risk_sd;
    double low_outlier_thresh = case_high_risk_means - outlier_sd * case_high_risk_sd;
    IntegerMatrix all_high_risk = transpose(cbind(case_high_inf, comp_high_inf));
    int n_high_risk = all_high_risk.nrow();
    LogicalVector outliers(n_target, false);

    if (n_pos > 0){

      IntegerMatrix positive_high_risk =


    }


    if (any(pos.risk)) {

      positive.high.risk <- all.high.risk[, pos.risk, drop = FALSE]
      low.outlier.mat <- matrix(rep(low.outlier.thresh[pos.risk], n.high.risk), nrow = n.high.risk,
                                byrow = TRUE)
      positive.outliers <- colSums(positive.high.risk < low.outlier.mat) > 0
      outliers[pos.risk] <- positive.outliers

    }

    if (any(neg.risk)) {

      negative.high.risk <- all.high.risk[, neg.risk, drop = FALSE]
      high.outlier.mat <- matrix(rep(high.outlier.thresh[neg.risk], n.high.risk), nrow = n.high.risk,
                                 byrow = TRUE)
      negative.outliers <- colSums(negative.high.risk > high.outlier.mat) > 0
      outliers[neg.risk] <- negative.outliers

    }

    // if there are outliers, recompute the weights and associated statistics ###
    if (any(outliers)) {

      risk.set.alleles[outliers] <- "2"
      pos.outlier.cols <- pos.risk[outliers[pos.risk]]
      neg.outlier.cols <- neg.risk[outliers[neg.risk]]

      // recode instances where the model appears to be recessive
      if (length(pos.outlier.cols) > 0) {

        case.inf[pos.outlier.cols, ][case.inf[pos.outlier.cols, ] == 1] <- 0
        comp.inf[pos.outlier.cols, ][comp.inf[pos.outlier.cols, ] == 1] <- 0
        both.one.inf[pos.outlier.cols, ] <- FALSE

      }

      if (length(neg.outlier.cols) > 0) {

        case.inf[neg.outlier.cols, ][case.inf[neg.outlier.cols, ] == 1] <- 2
        comp.inf[neg.outlier.cols, ][comp.inf[neg.outlier.cols, ] == 1] <- 2
        both.one.inf[neg.outlier.cols, ] <- FALSE

      }

      //recompute the number of informative families
      cases.minus.complements <- sign(case.inf - comp.inf)
      case.comp.diff <- cases.minus.complements != 0
      total.different.snps <- colSums(case.comp.diff)
      informative.families <- total.different.snps != 0
      inf.family.rows <- inf.family.rows[informative.families]
      n.informative.families <- sum(informative.families)
      case.inf <- case.inf[, informative.families]
      comp.inf <- comp.inf[, informative.families]
      cases.minus.complements.inf <- cases.minus.complements[, informative.families]

      // re-compute weights ###
      both.one.inf <- both.one.inf[, informative.families]
      both.one <- colSums(both.one.inf)
      weighted.informativeness <- n.both.one.weight * both.one + n.different.snps.weight * total.different.snps[informative.families]
      family.weights <- weight.lookup[weighted.informativeness]
      invsum.family.weights <- 1/sum(family.weights)

      // compute weighted difference vectors for cases vs complements ###
      dif.vecs <- as.matrix((family.weights *invsum.family.weights)*t(cases.minus.complements.inf))
      sum.dif.vecs <- colSums(dif.vecs)

      // re-compute proposed risk set ###
      risk.dirs <- sign(sum.dif.vecs)
      pos.risk <- which(risk.dirs > 0)
      neg.risk <- which(risk.dirs <= 0)

      n.pos <- length(pos.risk)
      n.neg <- length(neg.risk)
      if (n.neg > 0) {
        n1 <- case.inf[neg.risk, , drop = FALSE] < 2
        n2 <- comp.inf[neg.risk, , drop = FALSE] < 2
      }
      if (n.pos > 0) {
        p1 <- case.inf[pos.risk, , drop = FALSE] > 0
        p2 <- comp.inf[pos.risk, , drop = FALSE] > 0
      }
      if (n.pos == 0) {
        case.high.risk <- colSums(n1) == n.target
        comp.high.risk <- colSums(n2) == n.target
      } else if (n.neg == 0) {
        case.high.risk <- colSums(p1) == n.target
        comp.high.risk <- colSums(p2) == n.target
      } else {
        case.high.risk <- (colSums(p1) + colSums(n1)) == n.target
        comp.high.risk <- (colSums(p2) + colSums(n2)) == n.target
      }
      n.case.high.risk <- sum(case.high.risk)
        n.comp.high.risk <- sum(comp.high.risk)
        case.high.inf <- case.inf[, case.high.risk, drop = FALSE]
      comp.high.inf <- comp.inf[, comp.high.risk, drop = FALSE]

    }
  }

### count the number of risk alleles in those with the full risk set ###
                  case.high.risk.alleles <- colSums(case.high.inf[pos.risk, , drop = FALSE]) +
                    (colSums(2 - case.high.inf[neg.risk, , drop = FALSE]))
                    total.case.high.risk.alleles <- sum(case.high.risk.alleles)
                    comp.high.risk.alleles <- colSums(comp.high.inf[pos.risk, , drop = FALSE]) +
                      (colSums(2 - comp.high.inf[neg.risk, , drop = FALSE]))
                    total.comp.high.risk.alleles <- sum(comp.high.risk.alleles)

### compute scaling factor ###
                    rr <- total.case.high.risk.alleles/(total.case.high.risk.alleles + total.comp.high.risk.alleles)
                      rr <- ifelse(rr <= 0 | is.na(rr), 10^-10, rr)

### compute pseudo hotelling t2 ###
                      mu.hat <- sum.dif.vecs
                        mu.hat.mat <- matrix(rep(mu.hat, n.informative.families), nrow = n.informative.families, byrow = TRUE)
                        x <- t(as.matrix(cases.minus.complements.inf))
                        x.minus.mu.hat <- x - mu.hat.mat
                        weighted.x.minus.mu.hat <- family.weights * x.minus.mu.hat
                        cov.mat <- (invsum.family.weights) * crossprod(weighted.x.minus.mu.hat, x.minus.mu.hat)

                        target.block.ld.mat <- block.ld.mat[target.snps, target.snps]
                      cov.mat[!target.block.ld.mat] <- 0
                      elem.vars <- sqrt(diag(cov.mat))
                        sum.dif.vecs <- sum.dif.vecs/elem.vars
                        zero.var <- elem.vars == 0
                      sum.dif.vecs[zero.var] <- 10^-10

# compute svd of dif.vec.cov.mat
                      cov.mat.svd <- svd(cov.mat)
                        cov.mat.svd$d[cov.mat.svd$d < sqrt(.Machine$double.eps)] <- 10^10

# compute fitness score
                      pseudo.t2 <- (1/(invsum.family.weights*1000)) * rowSums((t(mu.hat) %*% cov.mat.svd$u)^2/cov.mat.svd$d)
                        fitness.score <- (rr^2) * pseudo.t2

# if desired, return the required information for the epistasis test
                        if (epi.test){

                          high.risk.families <- inf.family.rows
                          return(list(fitness.score = fitness.score, sum.dif.vecs = sum.dif.vecs, rr = rr, pseudo.t2 = pseudo.t2,
                                      risk.set.alleles = risk.set.alleles, high.risk.families = high.risk.families))

                        } else {

                          return(list(fitness.score = fitness.score, sum.dif.vecs = sum.dif.vecs, rr = rr, pseudo.t2 = pseudo.t2,
                                      risk.set.alleles = risk.set.alleles))

                        }

}
