#include <Rcpp.h>
using namespace Rcpp;

#include <cmath>
#include <algorithm>

// generic function for manhattan distance
template <typename InputIterator1, typename InputIterator2>
inline double manhattan_distance(InputIterator1 begin1, InputIterator1 end1,
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
    rval += std::abs(d1 - d2);
  }
  return rval;
}

// [[Rcpp::export]]
NumericMatrix rcpp_manhattan_distance(NumericMatrix mat) {

  // allocate the matrix we will return
  NumericMatrix rmat(mat.nrow(), mat.nrow());

  for (int i = 0; i < rmat.nrow(); i++) {
    for (int j = 0; j < rmat.nrow(); j++) {

      // rows we will operate on
      NumericMatrix::Row row1 = mat.row(i);
      NumericMatrix::Row row2 = mat.row(j);

      // calculate divergences
      double d1 = manhattan_distance(row1.begin(), row1.end(), row2.begin());

      // write to output matrix
      rmat(i,j) = d1;
    }
  }

  return rmat;
}

