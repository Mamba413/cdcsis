//
// Created by JinZhu on 2018/12/27.
//

#include "cdc_stats.h"
#include "utility.h"
#include <tuple>

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


std::vector<double> ConditionDistanceCovarianceStats::compute_condition_distance_correlation(
        std::vector<std::vector<double>> &distance_x, std::vector<std::vector<double>> &distance_y,
        std::vector<std::vector<double>> &kernel_density_estimation) {

    uint num = (uint) distance_x.size();
    std::vector<std::vector<double>> anova_x(num, std::vector<double>(num));
    std::vector<std::vector<double>> anova_y(num, std::vector<double>(num));
    std::vector<double> condition_distance_covariance_xy(num);
    std::vector<double> condition_distance_covariance_xx(num);
    std::vector<double> condition_distance_covariance_yy(num);

    for (uint i = 0; i < num; i++) {
        anova_x = weight_distance_anova(distance_x, kernel_density_estimation[i]);
        anova_y = weight_distance_anova(distance_y, kernel_density_estimation[i]);

        for (uint k = 0; k < num; k++) {
            for (uint j = 0; j < num; j++) {
                condition_distance_covariance_xy[i] += anova_x[k][j] * anova_y[k][j] * kernel_density_estimation[i][k] *
                                                       kernel_density_estimation[i][j];
                condition_distance_covariance_xx[i] += anova_x[k][j] * anova_x[k][j] * kernel_density_estimation[i][k] *
                                                       kernel_density_estimation[i][j];
                condition_distance_covariance_yy[i] += anova_y[k][j] * anova_y[k][j] * kernel_density_estimation[i][k] *
                                                       kernel_density_estimation[i][j];
            }
        }
    }

    for (uint i = 0; i < num; i++) {
        condition_distance_covariance_xy[i] /= sqrt(
                condition_distance_covariance_xx[i] * condition_distance_covariance_yy[i]);
    }

    return condition_distance_covariance_xy;
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
        default:
            this->condition_distance_covariance_stats = compute_condition_distance_covariance_stats(this->distance_x,
                                                                                                    this->distance_y,
                                                                                                    this->kernel_density_estimation);
    }
}

void ConditionDistanceCovarianceStats::setDistance_x(const std::vector<std::vector<double>> &distance_x) {
    ConditionDistanceCovarianceStats::distance_x = distance_x;
}

void ConditionDistanceCovarianceStats::setKernel_density_estimation(
        const std::vector<std::vector<double>> &kernel_density_estimation) {
    ConditionDistanceCovarianceStats::kernel_density_estimation = kernel_density_estimation;
}

std::vector<std::vector<double >> ConditionBallCovarianceStats::compute_weight_delta_xy_matrix(
        std::vector<std::vector<double >> &delta_y_matrix, std::vector<std::vector<double >> &distance_x,
        std::vector<std::vector<double >> &distance_y, std::vector<double> &weight) {

    std::vector<std::vector<double >> weight_delta_xy_matrix;
    for (uint i = 0; i < delta_y_matrix.size(); ++i) {
        weight_delta_xy_matrix.push_back(compute_weight_delta_xy_vector(delta_y_matrix[i], distance_x[i], distance_y[i],
                                                                        weight));
    }

    return weight_delta_xy_matrix;
}

std::vector<double> ConditionBallCovarianceStats::compute_condition_ball_covariance(
        std::vector<std::vector<double>> &distance_x,
        std::vector<std::vector<double>> &distance_y,
        std::vector<std::vector<double>> &kernel_density_estimation) {
    int num = (int) kernel_density_estimation.size();
    std::vector<double> condition_ball_covariance(kernel_density_estimation.size(), 0.0);

    std::vector<std::vector<double >> weight_delta_x_matrix, weight_delta_y_matrix, weight_delta_xy_matrix;
    double kernel_sum;
    for (int l = 0; l < num; ++l) {
        weight_delta_x_matrix = compute_weight_delta_x_matrix(distance_x, kernel_density_estimation[l]);
        weight_delta_y_matrix = compute_weight_delta_x_matrix(distance_y, kernel_density_estimation[l]);
        weight_delta_xy_matrix = compute_weight_delta_xy_matrix(weight_delta_y_matrix, distance_x, distance_y,
                                                                kernel_density_estimation[l]);
        kernel_sum = vector_sum(kernel_density_estimation[l]);
        for (int i = 0; i < num; ++i) {
            for (int j = 0; j < num; ++j) {
                condition_ball_covariance[l] += pow(
                        weight_delta_xy_matrix[i][j] -
                        weight_delta_x_matrix[i][j] * weight_delta_y_matrix[i][j] / kernel_sum, 2);
            }
        }
        condition_ball_covariance[l] /= (kernel_sum * kernel_sum * num * num);
    }

    return condition_ball_covariance;
}

std::vector<double> ConditionBallCovarianceStats::compute_condition_ball_correlation(
        std::vector<std::vector<double>> &distance_x, std::vector<std::vector<double>> &distance_y,
        std::vector<std::vector<double>> &kernel_density_estimation) {
    std::vector<double> condition_ball_correlation(kernel_density_estimation.size(), 0.0);
    return condition_ball_correlation;
}

void ConditionBallCovarianceStats::compute_stats() {
    switch (this->statsType) {
        case CONDITION_BALL_COVARIANCE:
            this->condition_ball_covariance_stats = compute_condition_ball_covariance_stats(this->distance_x,
                                                                                            this->distance_y,
                                                                                            this->kernel_density_estimation);
            break;
        case CONDITION_BALL_CORRELATION:
            this->condition_ball_covariance_stats = compute_condition_ball_correlation_stats(this->distance_x,
                                                                                             this->distance_y,
                                                                                             this->kernel_density_estimation);
            break;
        default:
            this->condition_ball_covariance_stats = compute_condition_ball_covariance_stats(this->distance_x,
                                                                                            this->distance_y,
                                                                                            this->kernel_density_estimation);
    }
}

double ConditionBallCovarianceStats::getCondition_ball_covariance_stats() const {
    return condition_ball_covariance_stats;
}

const std::vector<double> &ConditionBallCovarianceStats::getCondition_ball_covariance() const {
    return condition_ball_covariance;
}

void ConditionBallCovarianceStats::setDistance_x(const std::vector<std::vector<double>> &distance_x) {
    ConditionBallCovarianceStats::distance_x = distance_x;
}

void ConditionBallCovarianceStats::setKernel_density_estimation(
        const std::vector<std::vector<double>> &kernel_density_estimation) {
    ConditionBallCovarianceStats::kernel_density_estimation = kernel_density_estimation;
}
