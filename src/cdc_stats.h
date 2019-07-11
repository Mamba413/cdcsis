//
// Created by JinZhu on 2018/12/27.
//

#ifndef CDC_STATS_H
#define CDC_STATS_H

#include <vector>
#include "utility.h"
#include "global.h"

class ConditionDistanceCovarianceStats {

public:

    ConditionDistanceCovarianceStats(std::vector<std::vector<double>> &distance_x,
                                     std::vector<std::vector<double>> &distance_y,
                                     std::vector<std::vector<double>> &kernel_density_estimation, int stats_type) {
        this->distance_x = distance_x;
        this->distance_y = distance_y;
        this->kernel_density_estimation = kernel_density_estimation;
        this->statsType = StatsType(stats_type);
    }

    void compute_stats();

    double getCondition_distance_covariance_stats() const;

    const std::vector<double> &getCondition_distance_covariance() const;

    void setDistance_x(const std::vector<std::vector<double>> &distance_x);

    void setKernel_density_estimation(const std::vector<std::vector<double>> &kernel_density_estimation);

private:
    std::vector<std::vector<double>> distance_x;
    std::vector<std::vector<double>> distance_y;
    std::vector<std::vector<double>> kernel_density_estimation;
    StatsType statsType;
    double condition_distance_covariance_stats = 0.0;
    std::vector<double> condition_distance_covariance;

    std::vector<double> compute_condition_distance_covariance(std::vector<std::vector<double>> &distance_x,
                                                              std::vector<std::vector<double>> &distance_y,
                                                              std::vector<std::vector<double>> &kernel_density_estimation);

    std::vector<double> compute_condition_distance_correlation(std::vector<std::vector<double>> &distance_x,
                                                               std::vector<std::vector<double>> &distance_y,
                                                               std::vector<std::vector<double>> &kernel_density_estimation);

    double compute_condition_distance_covariance_stats(std::vector<std::vector<double>> &distance_x,
                                                       std::vector<std::vector<double>> &distance_y,
                                                       std::vector<std::vector<double>> &kernel_density_estimation) {
        std::vector<double> condition_distance_covariance(distance_x.size());
        condition_distance_covariance = compute_condition_distance_covariance(distance_x, distance_y,
                                                                              kernel_density_estimation);
        this->condition_distance_covariance = condition_distance_covariance;
        uint num = condition_distance_covariance.size();
        double kernel_sum_square, condition_distance_covariance_stats = 0.0;
        for (uint i = 0; i < num; ++i) {
            kernel_sum_square = vector_sum(kernel_density_estimation[i]);
            condition_distance_covariance_stats += 12 * pow(kernel_sum_square / num, 4.0) *
                                                   condition_distance_covariance[i];
        }
        condition_distance_covariance_stats /= num;
        return condition_distance_covariance_stats;
    }

    double compute_condition_distance_correlation_stats(std::vector<std::vector<double>> &distance_x,
                                                        std::vector<std::vector<double>> &distance_y,
                                                        std::vector<std::vector<double>> &kernel_density_estimation) {
        std::vector<double> condition_distance_correlation(distance_x.size());
        condition_distance_correlation = compute_condition_distance_correlation(distance_x, distance_y,
                                                                                kernel_density_estimation);
        this->condition_distance_covariance = condition_distance_correlation;
        return vector_mean(condition_distance_correlation);
    }
};


class ConditionBallCovarianceStats {

public:
    ConditionBallCovarianceStats() = default;

    ConditionBallCovarianceStats(std::vector<std::vector<double>> &distance_x,
                                 std::vector<std::vector<double>> &distance_y,
                                 std::vector<std::vector<double>> &kernel_density_estimation, int stats_type) {
        this->distance_x = distance_x;
        this->distance_y = distance_y;
        this->kernel_density_estimation = kernel_density_estimation;
        this->statsType = StatsType(stats_type);
    };

private:
    friend class CDCStatsticsMethod;

    std::vector<std::vector<double >> compute_weight_delta_x_matrix(std::vector<std::vector<double >> &distance_x,
                                                                    std::vector<double> &weight) {
        std::vector<std::vector<double >> delta_x_matrix;
        for (std::vector<double> distance_vector : distance_x) {
            delta_x_matrix.push_back(compute_weight_delta_x_vector(distance_vector, weight));
        }
        return delta_x_matrix;
    };

private:

    std::vector<std::vector<double>> distance_x;
    std::vector<std::vector<double>> distance_y;
    std::vector<std::vector<double>> kernel_density_estimation;
    StatsType statsType;

    double condition_ball_covariance_stats = 0.0;
    std::vector<double> condition_ball_covariance;
};


#endif //CDC_STATS_H
