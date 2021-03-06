#include "global.h"
#include "utility.h"
#include "kde.h"
#include "cdc_statstics_method.h"
#include "cdc_stats.h"
#include "dataRcpp.h"

#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
Rcpp::List cdcsisCpp(unsigned int stats_method, 
                     Rcpp::NumericMatrix& x, std::vector<unsigned int> &variable_index, 
                     Rcpp::NumericMatrix& y, Rcpp::NumericMatrix& z, 
                     std::vector<double> &bandwidth, double distance_index, unsigned int num_threads, 
                     unsigned int num_bootstrap, unsigned int seed, unsigned int stats_type, 
                     unsigned int kernel_type, unsigned int distance_z) {
  Rcpp::List result;
  
  try {
    StatsMethod statsMethod = StatsMethod(stats_method);
    bool kde_distance = false;
    if (distance_z) {
      kde_distance = true;
    }
    std::vector<std::vector<double> > kernel = rcpp_matrix_to_vector2d<double>(z);
    KernelDensityEstimation kernelDensityEstimation = KernelDensityEstimation(kernel, bandwidth, kernel_type, kde_distance);
    kernelDensityEstimation.compute_kernel_density_estimate();
    kernel = kernelDensityEstimation.get_kernel_density_estimate();
    
    std::vector<std::vector<double> > distance_y = rcpp_distance_matrix_to_vector2d<double>(y);
    
    CDCStatisticsMethod cDCStatisticsMethod = CDCStatisticsMethod();
    
    if (statsMethod == HYPOTHESIS_TESTS) {
      
      std::vector<std::vector<double> > distance_x = rcpp_distance_matrix_to_vector2d<double>(x);
      
      if (stats_type == 1) {
        cDCStatisticsMethod.conduct_cdc_test(distance_x, distance_y, kernel, num_bootstrap, seed, num_threads);
      }
      
      result.push_back(cDCStatisticsMethod.getCdc_statistic(), "statistic");
      result.push_back(cDCStatisticsMethod.getP_value(), "p.value");
      
    } else if (statsMethod == SURE_INDEPENDENCE_SCREENING) {
      
      std::vector<std::vector<double> > x_matrix = rcpp_matrix_to_vector2d<double>(x);
      
      if (variable_index.size() == 0) {
        cDCStatisticsMethod.conduct_cdc_screening(x_matrix, distance_y, kernel, distance_index);
      } else {
        cDCStatisticsMethod.conduct_cdc_screening(x_matrix, variable_index, distance_y, kernel, distance_index);
      }
      
      result.push_back(cDCStatisticsMethod.getCdc_statistic(), "statistic");
    } else if (statsMethod == STATISTICS_VALUE) {
      std::vector<std::vector<double> > distance_x = rcpp_distance_matrix_to_vector2d<double>(x);
      ConditionDistanceCovarianceStats conditionCovarianceStats = ConditionDistanceCovarianceStats(distance_x, 
                                                                                                            distance_y, 
                                                                                                            kernel, 
                                                                                                            stats_type);
      conditionCovarianceStats.compute_stats();
      result.push_back(conditionCovarianceStats.getCondition_distance_covariance_stats(), "statistic");
      result.push_back(conditionCovarianceStats.getCondition_distance_covariance(), "cdc");
    }
  } catch (std::exception& e) {
    if (strcmp(e.what(), "User interrupt.") != 0) {
      Rcpp::Rcerr << "Error: " << e.what() << " CDCSIS will EXIT now.\n";
    }
    return result;
  }
  
  return result;
}
