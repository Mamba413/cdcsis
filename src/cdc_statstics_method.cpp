//
// Created by JinZhu on 2018/12/28.
//

#include "cdc_statstics_method.h"
#include "cdc_stats.h"
#include "utility.h"

#include <algorithm>

const double &CDCStatsticsMethod::getP_value() const {
    return p_value;
}

const std::vector<double> &CDCStatsticsMethod::getCdc_statistic() const {
    return cdc_statistic;
}


void CDCStatsticsMethod::conduct_cdc_test(std::vector<std::vector<double>> &distance_x,
                                          std::vector<std::vector<double>> &distance_y,
                                          std::vector<std::vector<double>> &kernel, uint num_bootstrap, uint seed,
                                          uint num_threads) {

    ConditionDistanceCovarianceStats conditionDistanceCovarianceStats = ConditionDistanceCovarianceStats(distance_x,
                                                                                                         distance_y,
                                                                                                         kernel, 1);
    conditionDistanceCovarianceStats.compute_stats();
    this->cdc_statistic.push_back(conditionDistanceCovarianceStats.getCondition_distance_covariance_stats());

    if (num_bootstrap != 0) {
        // bootstrap sample index:
        if (seed == 0) {
            std::random_device random_device;
            random_number_generator.seed(random_device());
        } else {
            random_number_generator.seed(seed);
        }
        std::vector<std::vector<uint>> random_sample_index = generate_random_sample_index(num_bootstrap, kernel,
                                                                                          random_number_generator);

        // bootstrap test:
        std::vector<std::vector<double >> bootstrap_distance_x;
        for (uint i = 0; i < num_bootstrap; ++i) {
            bootstrap_distance_x = rearrange_matrix(distance_x, random_sample_index[i]);
            conditionDistanceCovarianceStats.setDistance_x(bootstrap_distance_x);
            conditionDistanceCovarianceStats.compute_stats();
            this->permuted_cdc_statistic.push_back(
                    conditionDistanceCovarianceStats.getCondition_distance_covariance_stats());
        }

        // compute p-value:
        this->p_value = compute_p_value(this->permuted_cdc_statistic, this->cdc_statistic[0]);
    }
}

void CDCStatsticsMethod::conduct_cdc_screening(std::vector<std::vector<double>> &x,
                                               std::vector<std::vector<double>> &distance_y,
                                               std::vector<std::vector<double>> &kernel, double index) {
    std::vector<std::vector<double>> distance_x;
    ConditionDistanceCovarianceStats conditionDistanceCovarianceStats = ConditionDistanceCovarianceStats(distance_x,
                                                                                                         distance_y,
                                                                                                         kernel, 2);
    for (std::vector<double> univariate_x : x) {
        distance_x = Euclidean_distance(univariate_x, index);
        conditionDistanceCovarianceStats.setDistance_x(distance_x);
        conditionDistanceCovarianceStats.compute_stats();
        this->cdc_statistic.push_back(conditionDistanceCovarianceStats.getCondition_distance_covariance_stats());
    }
}

void CDCStatsticsMethod::conduct_cdc_screening(std::vector<std::vector<double>> &x,
                                               std::vector<uint> &x_variable_index,
                                               std::vector<std::vector<double>> &distance_y,
                                               std::vector<std::vector<double>> &kernel, double index) {
    std::vector<uint> variable_sequence;
    std::vector<std::vector<double>> submatrix;
    std::vector<std::vector<double>> distance_x;
    ConditionDistanceCovarianceStats conditionDistanceCovarianceStats = ConditionDistanceCovarianceStats(distance_x,
                                                                                                         distance_y,
                                                                                                         kernel, 2);

    for (int i = (int) (x_variable_index.size() - 1); i >= 0; --i) {
        // extract variable index
        if (i != 0) {
            variable_sequence = generate_sequence(x_variable_index[i - 1] + 1, x_variable_index[i]);
        } else {
            variable_sequence = generate_sequence(0, x_variable_index[0]);
        }
        submatrix = extract_submatrix_byrow(x, variable_sequence);

        // compute statistics
        distance_x = Euclidean_distance_byrow(submatrix, index);
        conditionDistanceCovarianceStats.setDistance_x(distance_x);
        conditionDistanceCovarianceStats.compute_stats();
        this->cdc_statistic.push_back(conditionDistanceCovarianceStats.getCondition_distance_covariance_stats());
    }
    reverse(this->cdc_statistic.begin(), this->cdc_statistic.end());
}