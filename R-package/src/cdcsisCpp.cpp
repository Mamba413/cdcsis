#include "global.h"
#include "utility.h"
#include "kde.h"
#include "cdc_statstics_method.h"
#include "cdc_stats.h"
#include "dataRcpp.h"

#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
Rcpp::List cdcsisCpp(uint stats_method, Rcpp::NumericMatrix& x, std::vector<uint> &variable_index, 
                     Rcpp::NumericMatrix& y, Rcpp::NumericMatrix& z, 
                     std::vector<double> &bandwidth, double distance_index, uint num_threads, 
                     uint num_bootstrap, uint seed, uint stats_type) {
  Rcpp::List result;
  
  try {
    StatsMethod statsMethod = StatsMethod(stats_method);
    
    std::vector<std::vector<double> > kernel = rcpp_matrix_to_vector2d<double>(z);
    KernelDensityEstimation kernelDensityEstimation = KernelDensityEstimation(kernel, bandwidth, 1);
    kernelDensityEstimation.compute_kernel_density_estimate();
    kernel = kernelDensityEstimation.get_kernel_density_estimate();
    
    std::vector<std::vector<double> > distance_y = rcpp_distance_matrix_to_vector2d<double>(y);
    
    CDCStatsticsMethod cDCStatsticsMethod = CDCStatsticsMethod();
    
    if (statsMethod == HYPOTHESIS_TESTS) {
      
      std::vector<std::vector<double> > distance_x = rcpp_distance_matrix_to_vector2d<double>(x);
      
      if (stats_type == 1) {
          cDCStatsticsMethod.conduct_cdc_test(distance_x, distance_y, kernel, num_bootstrap, seed, num_threads);
      } else if (stats_type == 3) {
          cDCStatsticsMethod.conduct_cbc_test(distance_x, distance_y, kernel, num_bootstrap, seed, num_threads);
      }
      
      result.push_back(cDCStatsticsMethod.getCdc_statistic(), "statistic");
      result.push_back(cDCStatsticsMethod.getP_value(), "p.value");
      
    } else if (statsMethod == SURE_INDEPENDENCE_SCREENING) {
      
      std::vector<std::vector<double> > x_matrix = rcpp_matrix_to_vector2d<double>(x);
      
      if (variable_index.size() == 0) {
        cDCStatsticsMethod.conduct_cdc_screening(x_matrix, distance_y, kernel, distance_index);
      } else {
        cDCStatsticsMethod.conduct_cdc_screening(x_matrix, variable_index, distance_y, kernel, distance_index);
      }
      
      result.push_back(cDCStatsticsMethod.getCdc_statistic(), "statistic");
    } else if (statsMethod == STATISTICS_VALUE) {
      std::vector<std::vector<double> > distance_x = rcpp_distance_matrix_to_vector2d<double>(x);
      
      ConditionDistanceCovarianceStats conditionDistanceCovarianceStats = ConditionDistanceCovarianceStats(distance_x, 
                                                                                                           distance_y, 
                                                                                                           kernel, 
                                                                                                           stats_type);
      conditionDistanceCovarianceStats.compute_stats();
      result.push_back(conditionDistanceCovarianceStats.getCondition_distance_covariance_stats(), "statistic");
    }
  } catch (std::exception& e) {
    if (strcmp(e.what(), "User interrupt.") != 0) {
      Rcpp::Rcerr << "Error: " << e.what() << " CDCSIS will EXIT now.\n";
    }
    return result;
  }
  
  return result;
}
