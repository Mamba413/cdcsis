//
// Created by JinZhu on 2018/12/27.
//

#ifndef CDC_STATS_H
#define CDC_STATS_H

#include <vector>
#include "global.h"

class ConditionDistanceCovarianceStats {

public:
    double getCondition_distance_covariance_stats() const;

    const std::vector<double> &getCondition_distance_covariance() const;

    ConditionDistanceCovarianceStats(std::vector<std::vector<double>> &distance_x,
                                     std::vector<std::vector<double>> &distance_y,
                                     std::vector<std::vector<double>> &kernel_density_estimation, int stats_type) {
        this->distance_x = distance_x;
        this->distance_y = distance_y;
        this->kernel_density_estimation = kernel_density_estimation;
        this->statsType = StatsType(stats_type);
    }

    void compute_stats();

private:
    StatsType statsType;
    std::vector<std::vector<double>> distance_x;
    std::vector<std::vector<double>> distance_y;
    std::vector<std::vector<double>> kernel_density_estimation;
    double condition_distance_covariance_stats;
    std::vector<double> condition_distance_covariance;

    std::vector<double> compute_condition_distance_covariance(std::vector<std::vector<double>> &distance_x,
                                                              std::vector<std::vector<double>> &distance_y,
                                                              std::vector<std::vector<double>> &kernel_density_estimation);

    double compute_condition_distance_covariance_stats(std::vector<std::vector<double>> &distance_x,
                                                       std::vector<std::vector<double>> &distance_y,
                                                       std::vector<std::vector<double>> &kernel_density_estimation) {
        uint num = (uint) distance_x.size();
        std::vector<double> condition_distance_covariance(num);
        condition_distance_covariance = compute_condition_distance_covariance(distance_x, distance_y,
                                                                              kernel_density_estimation);

        double condition_distance_covariance_stats = 0.0;
        for (double cdc : condition_distance_covariance) {
            condition_distance_covariance_stats += cdc;
        }

        return condition_distance_covariance_stats / num;
    }

    double compute_condition_distance_correlation_stats(std::vector<std::vector<double>> &distance_x,
                                                        std::vector<std::vector<double>> &distance_y,
                                                        std::vector<std::vector<double>> &kernel_density_estimation) {
        return compute_condition_distance_covariance_stats(distance_x, distance_y, kernel_density_estimation);
    }
};


#endif //CDC_STATS_H
