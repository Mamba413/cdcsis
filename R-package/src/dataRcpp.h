#ifndef DATARCPP_H_
#define DATARCPP_H_

#include <vector>
#include <Rcpp.h>
#include "global.h"

template<typename T>
std::vector<std::vector<T>> rcpp_matrix_to_vector2d(Rcpp::NumericMatrix& x) {
  
  std::vector<std::vector<T>> vector2d(x.nrow(), std::vector<T>(x.ncol()));
  
  for (uint i = 0; i < (uint) x.nrow(); ++i) {
    for (uint j = 0; j < (uint) x.ncol(); ++j) {
      vector2d[i][j] = x(i, j);
    }
  }
  
  return vector2d;
}


template<typename T>
std::vector<std::vector<T>> rcpp_distance_matrix_to_vector2d(Rcpp::NumericMatrix& x) {
  
  uint num_row = (uint) x.nrow();
  std::vector<std::vector<T> > vector2d(num_row, std::vector<T>(num_row));
  
  for (uint i = 0; i < num_row; ++i) {
    vector2d[i][i] = 0.0;
    for (uint j = 0; j < i; ++j) {
      vector2d[i][j] = x(i, j);
      vector2d[j][i] = vector2d[i][j];
    }
  }
  
  return vector2d;
}

#endif /* DATARCPP_H_ */