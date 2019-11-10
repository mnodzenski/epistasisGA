#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector knn_prop_cases(NumericMatrix case_data, int k, double threshold, NumericVector case_status) {
  // identify the total number of informative families
  int n_informative_fams = case_data.nrow();
  // initialize vector of proportion of nearest neighbor cases
  NumericVector prop_cases(n_informative_fams);
  // loop over each individual and actually compute the proportion of knn cases
  for(int i = 0; i < n_informative_fams; ++i) {
    NumericVector distances = case_data(_, i);
    // make sure the subject, and the subject's complement, are not neighbors
    distances[i] = 99999;
    distances[i + n_informative_fams - 1] = 99999;
    // find the indices of nearest neighbors
    LogicalVector zeroes = distances == 0;
    LogicalVector ones = distances == 1;
    LogicalVector twos = distances == 2;
    int n_zeroes = sum(zeroes);
    int n_ones = sum(ones);
    int n_twos = sum(twos);
    int n_zeroes_ones = n_zeroes + n_ones;
    int n_zeroes_ones_twos = n_zeroes + n_ones + n_twos;
    //compute proportion of nearest neighbors that are also cases
    double prop;
    if (n_zeroes >= k){
      prop = sum(case_status[zeroes])/n_zeroes;
    } else if (n_zeroes_ones >= k){
      prop = (sum(case_status[zeroes]) + sum(case_status[ones]))/n_zeroes_ones;
    } else if (n_zeroes_ones_twos >= k){
      prop = (sum(case_status[zeroes]) + sum(case_status[ones]) + sum(case_status[twos]))/n_zeroes_ones_twos;
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



