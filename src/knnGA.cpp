#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector knn_prop_cases(NumericMatrix case_data, int k, double threshold, NumericVector case_status) {

  double n_informative_fams = case_data.nrow();
  // initialize vector of proportion of nearest neighbor cases
  NumericVector prop_cases(n_informative_fams);

  // loop over each individual and actually compute the proportion of knn cases
  for(int i = 0; i < n_informative_fams; ++i) {

    NumericVector distances = case_data(_, i);

    // make sure the subject, and the subject's complement, are not neighbors
    distances[i] = 99999;
    distances[i + n_informative_fams - 1] = 99999;

    // find the indices of nearest neighbors
    NumericVector ones_vec(n_informative_fams, 1);
    LogicalVector zeroes = distances == 0;
    LogicalVector ones = distances == 1;
    LogicalVector twos = distances == 2;
    NumericVector zeroes_n = ones_vec[zeroes];
    NumericVector ones_n = ones_vec[ones];
    NumericVector twos_n = ones_vec[twos];
    int n_zeroes = Rcpp::sum(zeroes_n);
    int n_ones = Rcpp::sum(ones_n);
    int n_twos = Rcpp::sum(twos_n);
    int n_zeroes_ones = n_zeroes + n_ones;
    int n_zeroes_ones_twos = n_zeroes + n_ones + n_twos;
    NumericVector zero_cases = case_status[zeroes];
    NumericVector one_cases = case_status[ones];
    NumericVector two_cases = case_status[twos];

    //compute proportion of nearest neighbors that are also cases
    double prop;
    if (n_zeroes >= k){
      prop = Rcpp::sum(zero_cases)/n_zeroes;
    } else if (n_zeroes_ones >= k){
      prop = (Rcpp::sum(zero_cases) + Rcpp::sum(one_cases))/n_zeroes_ones;
    } else if (n_zeroes_ones_twos >= k){
      prop = (Rcpp::sum(zero_cases) + Rcpp::sum(one_cases) + Rcpp::sum(two_cases))/n_zeroes_ones_twos;
    } else {
      prop = 0;
    }
    if (prop < threshold){
      prop = 0;
    }
    prop_cases[i] = prop;
  }
  return prop_cases;
}



