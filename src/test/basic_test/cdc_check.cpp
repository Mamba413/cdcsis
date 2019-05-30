//
// Created by JinZhu on 2018/12/27.
//

#include "gtest/gtest.h"
#include "utility.h"
#include "cdc_stats.h"
#include "cdc_statstics_method.h"
#include "global.h"
#include "kde.h"
#include "test_setting.h"

TEST(cdc_stats, compute_condition_distance_covariance) {
    std::vector<double> x = {-0.62645381, 0.18364332, -0.83562861, 1.59528080, 0.32950777, -0.82046838,
                             0.48742905, 0.73832471, 0.57578135, -0.30538839, 1.51178117, 0.38984324,
                             -0.62124058, -2.21469989, 1.12493092, -0.04493361, -0.01619026, 0.94383621,
                             0.82122120, 0.59390132, 0.91897737, 0.78213630, 0.07456498, -1.98935170, 0.61982575};
    std::vector<double> y = {-0.05612874, -0.15579551, -1.47075238, -0.47815006, 0.41794156, 1.35867955,
                             -0.10278773, 0.38767161, -0.05380504, -1.37705956, -0.41499456, -0.39428995,
                             -0.05931340, 1.10002537, 0.76317575, -0.16452360, -0.25336168, 0.69696338,
                             0.55666320, -0.68875569, -0.70749516, 0.36458196, 0.76853292, -0.11234621, 0.88110773};
    std::vector<double> z = {0.39810588, -0.61202639, 0.34111969, -1.12936310, 1.43302370, 1.98039990,
                             -0.36722148, -1.04413463, 0.56971963, -0.13505460, 2.40161776, -0.03924000,
                             0.68973936, 0.02800216, -0.74327321, 0.18879230, -1.80495863, 1.46555486,
                             0.15325334, 2.17261167, 0.47550953, -0.70994643, 0.61072635, -0.93409763, -1.25363340};

    uint num_row = (uint) x.size();
    std::vector<std::vector<double>> x_2d(num_row, std::vector<double>(1));
    std::vector<std::vector<double>> y_2d(num_row, std::vector<double>(1));
    std::vector<std::vector<double>> z_2d(num_row, std::vector<double>(1));

    x_2d = vector_to_matrix(x, num_row, 1);
    y_2d = vector_to_matrix(y, num_row, 1);
    z_2d = vector_to_matrix(z, num_row, 1);

    std::vector<std::vector<double>> distance_x(num_row, std::vector<double>(num_row));
    std::vector<std::vector<double>> distance_y(num_row, std::vector<double>(num_row));
    distance_x = Euclidean_distance(x_2d, 1.0);
    distance_y = Euclidean_distance(y_2d, 1.0);

    std::vector<std::vector<double>> kde(num_row, std::vector<double>(num_row));
    std::vector<std::vector<double>> bandwidth(1, std::vector<double>(1));
    bandwidth[0][0] = 1;

    KernelDensityEstimation kernelDensityEstimation = KernelDensityEstimation(z_2d, bandwidth, 1);
    kernelDensityEstimation.compute_kernel_density_estimate();
    kde = kernelDensityEstimation.get_kernel_density_estimate();

    ConditionDistanceCovarianceStats conditionDistanceCovarianceStats = ConditionDistanceCovarianceStats(distance_x,
                                                                                                         distance_y,
                                                                                                         kde, 1);
    conditionDistanceCovarianceStats.compute_stats();

    std::vector<double> condition_distance_covariance = conditionDistanceCovarianceStats.getCondition_distance_covariance();
    double condition_distance_covariance_stats = conditionDistanceCovarianceStats.getCondition_distance_covariance_stats();

    std::vector<double> condition_distance_covariance_true = {0.04175982, 0.04676984, 0.04230362, 0.04873413,
                                                              0.04071863, 0.07817091, 0.04625996, 0.04821924,
                                                              0.04000199, 0.04551928, 0.13606705, 0.04507440,
                                                              0.03875126, 0.04470311, 0.04706571, 0.04360601,
                                                              0.05651247, 0.04173558, 0.04387401, 0.10217581,
                                                              0.04098363, 0.04698362, 0.03957068, 0.04769151,
                                                              0.04967046};

    double abs_error = 0.000001;
    double condition_distance_covariance_stats_crude = compute_condition_distance_covariance_crude(distance_x,
                                                                                                   distance_y, kde);
    EXPECT_NEAR(condition_distance_covariance_stats_crude, condition_distance_covariance_stats, abs_error);

    ConditionDistanceCovarianceStats conditionDistanceCorrelationStats = ConditionDistanceCovarianceStats(distance_x,
                                                                                                          distance_y,
                                                                                                          kde, 2);
    conditionDistanceCorrelationStats.compute_stats();
    double condition_distance_correlation_stats = conditionDistanceCorrelationStats.getCondition_distance_covariance_stats();
    EXPECT_NEAR(0.1810813, condition_distance_correlation_stats, abs_error);
}

TEST(cdc_statistics_method, conduct_cdc_screening) {
    std::vector<double> x = {-0.62645381, 0.18364332, -0.83562861, 1.59528080, 0.32950777, -0.82046838,
                             0.48742905, 0.73832471, 0.57578135, -0.30538839, 1.51178117, 0.38984324,
                             -0.62124058, -2.21469989, 1.12493092, -0.04493361, -0.01619026, 0.94383621,
                             0.82122120, 0.59390132, 0.91897737, 0.78213630, 0.07456498, -1.98935170, 0.61982575};
    std::vector<double> y = {-0.05612874, -0.15579551, -1.47075238, -0.47815006, 0.41794156, 1.35867955,
                             -0.10278773, 0.38767161, -0.05380504, -1.37705956, -0.41499456, -0.39428995,
                             -0.05931340, 1.10002537, 0.76317575, -0.16452360, -0.25336168, 0.69696338,
                             0.55666320, -0.68875569, -0.70749516, 0.36458196, 0.76853292, -0.11234621, 0.88110773};
    std::vector<double> z = {0.39810588, -0.61202639, 0.34111969, -1.12936310, 1.43302370, 1.98039990,
                             -0.36722148, -1.04413463, 0.56971963, -0.13505460, 2.40161776, -0.03924000,
                             0.68973936, 0.02800216, -0.74327321, 0.18879230, -1.80495863, 1.46555486,
                             0.15325334, 2.17261167, 0.47550953, -0.70994643, 0.61072635, -0.93409763, -1.25363340};

    uint num_row = (uint) x.size();
    std::vector<std::vector<double>> z_2d(num_row, std::vector<double>(1));
    z_2d = vector_to_matrix(z, num_row, 1);

    std::vector<std::vector<double>> x_2d;
    x_2d.push_back(x);
    x_2d.push_back(x);
    x_2d.push_back(x);

    std::vector<std::vector<double>> distance_y = Euclidean_distance(y, 1.0);

    std::vector<std::vector<double>> bandwidth(1, std::vector<double>(1));
    bandwidth[0][0] = 1;

    KernelDensityEstimation kernelDensityEstimation = KernelDensityEstimation(z_2d, bandwidth, 1);
    kernelDensityEstimation.compute_kernel_density_estimate();
    std::vector<std::vector<double>> kde = kernelDensityEstimation.get_kernel_density_estimate();

    std::vector<uint> variable_index = {1, 2};

    double abs_error = 0.000001;
    CDCStatsticsMethod cdcStatsticsMethod = CDCStatsticsMethod();
    cdcStatsticsMethod.conduct_cdc_screening(x_2d, variable_index, distance_y, kde, 1.0);
    std::vector<double> result = cdcStatsticsMethod.getCdc_statistic();
    EXPECT_NEAR(0.1810813, result[1], abs_error);
}

TEST(cbc_stats, compute_condition_ball_covariance) {
    uint num_row;
    num_row = (uint) get_array_length(CONTINUOUS_X);
    std::vector<std::vector<double>> x_2d(num_row, std::vector<double>(1));
    std::vector<std::vector<double>> y_2d(num_row, std::vector<double>(1));
    std::vector<std::vector<double>> distance_x(num_row, std::vector<double>(num_row));
    std::vector<std::vector<double>> distance_y(num_row, std::vector<double>(num_row));
    std::vector<std::vector<double>> z_2d(num_row, std::vector<double>(1));
    std::vector<std::vector<double>> kde(num_row, std::vector<double>(num_row));
    double bandwidth = 1.0;
    double condition_ball_covariance_stats, condition_ball_covariance_stats_crude;

    // Continuous case:
    x_2d = vector_to_matrix(CONTINUOUS_X, num_row, 1);
    y_2d = vector_to_matrix(CONTINUOUS_Y, num_row, 1);
    distance_x = Euclidean_distance(x_2d, 1.0);
    distance_y = Euclidean_distance(y_2d, 1.0);

    z_2d = vector_to_matrix(CONTINUOUS_Z, num_row, 1);
    KernelDensityEstimation kernelDensityEstimation1 = KernelDensityEstimation(z_2d, bandwidth, 1);
    kernelDensityEstimation1.compute_kernel_density_estimate();
    kde = kernelDensityEstimation1.get_kernel_density_estimate();

    ConditionBallCovarianceStats conditionBallCovarianceStats1 = ConditionBallCovarianceStats(distance_x,
                                                                                              distance_y,
                                                                                              kde, 3);
    conditionBallCovarianceStats1.compute_stats();
    condition_ball_covariance_stats = conditionBallCovarianceStats1.getCondition_ball_covariance_stats();
    condition_ball_covariance_stats_crude = compute_condition_ball_covariance_crude(distance_x, distance_y, kde);
    EXPECT_NEAR(condition_ball_covariance_stats_crude, condition_ball_covariance_stats, ABSOLUATE_ERROR);

    // Discrete case:
    num_row = (uint) get_array_length(DISCRETE_X);
    x_2d = vector_to_matrix(DISCRETE_X, num_row, 1);
    y_2d = vector_to_matrix(DISCRETE_Y, num_row, 1);
    distance_x = Euclidean_distance(x_2d, 1.0);
    distance_y = Euclidean_distance(y_2d, 1.0);

    z_2d = vector_to_matrix(DISCRETE_Z, num_row, 1);
    KernelDensityEstimation kernelDensityEstimation2 = KernelDensityEstimation(z_2d, bandwidth, 1);
    kernelDensityEstimation2.compute_kernel_density_estimate();
    kde = kernelDensityEstimation2.get_kernel_density_estimate();

    ConditionBallCovarianceStats conditionBallCovarianceStats2 = ConditionBallCovarianceStats(distance_x,
                                                                                              distance_y,
                                                                                              kde, 3);
    conditionBallCovarianceStats2.compute_stats();
    condition_ball_covariance_stats = conditionBallCovarianceStats2.getCondition_ball_covariance_stats();
    condition_ball_covariance_stats_crude = compute_condition_ball_covariance_crude(distance_x, distance_y, kde);
    EXPECT_NEAR(condition_ball_covariance_stats_crude, condition_ball_covariance_stats, ABSOLUATE_ERROR);
}
