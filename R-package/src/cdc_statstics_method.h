//
// Created by JinZhu on 2018/12/28.
//

#ifndef CDC_STATISTICS_METHOD_H
#define CDC_STATISTICS_METHOD_H

#include <vector>
#include <random>
#include "global.h"

class CDCStatisticsMethod {

public:
    const double &getP_value() const;

    const std::vector<double> &getCdc_statistic() const;

    void conduct_cdc_test(std::vector<std::vector<double>> &distance_x, std::vector<std::vector<double>> &distance_y,
                          std::vector<std::vector<double>> &kernel, uint num_bootstrap, uint seed, uint num_threads);

    /**
     * Conduct conditional distance correlation based screening for multivariate explanation variable
     * @param x : p * N matrix, p is variable numbers and N is sample size
     * @param x_variable_index :
     * @param distance_y : N * N distance matrix of response variable
     * @param kernel : N * N kernel density density matrix
     * @param index : range within (0, 2]
     */
    void conduct_cdc_screening(std::vector<std::vector<double>> &x, std::vector<uint> &x_variable_index,
                               std::vector<std::vector<double>> &distance_y,
                               std::vector<std::vector<double>> &kernel, double index);

    /**
     * Conduct conditional distance correlation based screening for each univariate explanation variable
     * @param x : p * N matrix, p is variable numbers and N is sample size
     * @param distance_y : N * N distance matrix of response variable
     * @param kernel : N * N kernel density density matrix
     * @param index : range within (0, 2]
     */
    void conduct_cdc_screening(std::vector<std::vector<double>> &x,
                               std::vector<std::vector<double>> &distance_y,
                               std::vector<std::vector<double>> &kernel, double index);

private:
    double p_value;
    std::vector<double> cdc_statistic;
    std::vector<double> permuted_cdc_statistic;

    double compute_p_value(std::vector<double> permuted_statistic, double statistic) {
        double larger_num = 0.0;
        for (double value : permuted_statistic) {
            larger_num += value >= statistic;
        }
        // p-value:
        return (1.0 + larger_num) / (1.0 + (double) permuted_statistic.size());
    }

protected:
    std::mt19937_64 random_number_generator;
};


#endif //CDC_STATISTICS_METHOD_H
