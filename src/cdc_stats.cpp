//
// Created by JinZhu on 2018/12/27.
//

#include "cdc_stats.h"
#include "utility.h"

double ConditionDistanceCovarianceStats::getCondition_distance_covariance_stats() const {
    return condition_distance_covariance_stats;
}

const std::vector<double> &ConditionDistanceCovarianceStats::getCondition_distance_covariance() const {
    return condition_distance_covariance;
}

std::vector<double> ConditionDistanceCovarianceStats::compute_condition_distance_covariance(
        std::vector<std::vector<double>> &distance_x, std::vector<std::vector<double>> &distance_y,
        std::vector<std::vector<double>> &kernel_density_estimation) {

    uint num = (uint) distance_x.size();
    std::vector<std::vector<double>> anova_x(num, std::vector<double>(num));
    std::vector<std::vector<double>> anova_y(num, std::vector<double>(num));
    std::vector<double> condition_distance_covariance(num);

    double kernel_sum_square;
    for (uint i = 0; i < num; i++) {
        anova_x = weight_distance_anova(distance_x, kernel_density_estimation[i]);
        anova_y = weight_distance_anova(distance_y, kernel_density_estimation[i]);

        kernel_sum_square = vector_sum(kernel_density_estimation[i]);
        kernel_sum_square = kernel_sum_square * kernel_sum_square;

        for (uint k = 0; k < num; k++) {
            for (uint j = 0; j < num; j++) {
                condition_distance_covariance[i] += anova_x[k][j] * anova_y[k][j] * kernel_density_estimation[i][k] *
                                                    kernel_density_estimation[i][j];
            }
        }
        condition_distance_covariance[i] /= kernel_sum_square;
    }

    return condition_distance_covariance;
}


void ConditionDistanceCovarianceStats::compute_stats() {
    switch (this->statsType) {
        case CONDITION_DISTANCE_COVARIANCE:
            this->condition_distance_covariance_stats = compute_condition_distance_covariance_stats(this->distance_x,
                                                                                                   this->distance_y,
                                                                                                   this->kernel_density_estimation);
            break;
        case CONDITION_DISTANCE_CORRELATION:
            this->condition_distance_covariance_stats = compute_condition_distance_correlation_stats(this->distance_x,
                                                                                                     this->distance_y,
                                                                                                     this->kernel_density_estimation);
            break;
    }
}
