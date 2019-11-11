#include <Rcpp.h>
using namespace Rcpp;
#include <cmath>
#include <algorithm>

// [[Rcpp::export]]
NumericVector knn_prop_cases(int n_informative_fams, NumericMatrix case_distances, int k, double threshold, NumericVector case_status, NumericVector ones_vec) {

  // initialize vector of proportion of nearest neighbor cases
  NumericVector prop_cases(n_informative_fams);

  // loop over each individual and actually compute the proportion of knn cases
  for(int i = 0; i < n_informative_fams; ++i) {

    NumericVector knn_dist = case_distances(_, i);

    // make sure the subject, and the subject's complement, are not neighbors
    knn_dist[i] = 3;
    unsigned int complement_idx = i + n_informative_fams;
    knn_dist[complement_idx] = 3;

    // find the indices of nearest neighbors
    LogicalVector zeroes = knn_dist == 0;
    LogicalVector ones = knn_dist == 1;
    LogicalVector twos = knn_dist == 2;

    NumericVector zeroes_n = ones_vec[zeroes];
    NumericVector ones_n = ones_vec[ones];
    NumericVector twos_n = ones_vec[twos];
    unsigned int n_zeroes = sum(zeroes_n);
    unsigned int n_ones = sum(ones_n);
    unsigned int n_twos = sum(twos_n);
    unsigned int n_zeroes_ones = n_zeroes + n_ones;
    unsigned int n_zeroes_ones_twos = n_zeroes + n_ones + n_twos;

    //compute proportion of nearest neighbors that are also cases
    double prop;
    if (n_zeroes >= k){
      NumericVector zero_cases = case_status[zeroes];
      double sum_zero_cases = sum(zero_cases);
      prop = sum_zero_cases/n_zeroes;
    } else if (n_zeroes_ones >= k){
      NumericVector zero_cases = case_status[zeroes];
      double sum_zero_cases = sum(zero_cases);
      NumericVector one_cases = case_status[ones];
      double sum_one_cases = sum(one_cases);
      double sum_zero_one_cases = sum_zero_cases + sum_one_cases;
      prop = sum_zero_one_cases/n_zeroes_ones;
    } else if (n_zeroes_ones_twos >= k){
      NumericVector zero_cases = case_status[zeroes];
      double sum_zero_cases = sum(zero_cases);
      NumericVector one_cases = case_status[ones];
      double sum_one_cases = sum(one_cases);
      NumericVector two_cases = case_status[twos];
      double sum_two_cases = sum(two_cases);
      double sum_zero_one_two_cases = sum_zero_cases + sum_one_cases + sum_two_cases;
      prop = sum_zero_one_two_cases/n_zeroes_ones_twos;
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


// [[Rcpp::export]]
// generic function for city_block distance
template <typename InputIterator1, typename InputIterator2>
inline double city_block(InputIterator1 begin1, InputIterator1 end1,
                            InputIterator2 begin2) {

  // value to return
  double rval = 0;

  // set iterators to beginning of ranges
  InputIterator1 it1 = begin1;
  InputIterator2 it2 = begin2;

  // for each input item
  while (it1 != end1) {

    // take the value and increment the iterator
    double d1 = *it1++;
    double d2 = *it2++;

    // accumulate if appropirate
    if (d1 > 0 && d2 > 0)
      rval += std::abs(d1 - d2);
  }
  return rval;
}

// helper function for taking the average of two numbers
inline double average(double val1, double val2) {
  return (val1 + val2) / 2;
}

// [[Rcpp::export]]
NumericMatrix rcpp_city_block(NumericMatrix mat) {

  // allocate the matrix we will return
  NumericMatrix rmat(mat.nrow(), mat.nrow());

  for (int i = 0; i < rmat.nrow(); i++) {
    for (int j = 0; j < i; j++) {

      // rows we will operate on
      NumericMatrix::Row row1 = mat.row(i);
      NumericMatrix::Row row2 = mat.row(j);

      // compute the average using std::tranform from the STL
      std::vector<double> avg(row1.size());
      std::transform(row1.begin(), row1.end(), // input range 1
                     row2.begin(),             // input range 2
                     avg.begin(),              // output range
                     average);                 // function to apply

      // calculate divergences
      double d1 = kl_divergence(row1.begin(), row1.end(), avg.begin());
      double d2 = kl_divergence(row2.begin(), row2.end(), avg.begin());

      // write to output matrix
      rmat(i,j) = std::sqrt(.5 * (d1 + d2));
    }
  }

  return rmat;
}



