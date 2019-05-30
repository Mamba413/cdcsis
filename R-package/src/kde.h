//
// Created by JinZhu on 2018/12/26.
//

#ifndef KernelDensityEstimation_H
#define KernelDensityEstimation_H

#include <vector>
#include <iostream>
#include <math.h>
#include "global.h"

using namespace std;

class KernelDensityEstimation {
public:
    void setBandwidth_matrix(const vector<vector<double>> &bandwidth_matrix);

    void compute_kernel_density_estimate();

    std::vector<std::vector<double>> get_kernel_density_estimate() {
        return this->kernel_density_estimate;
    }

    KernelDensityEstimation(std::vector<std::vector<double>> &condition_variable,
                            std::vector<std::vector<double>> &bandwidth,
                            int kernel_type,
                            bool distance = false) {
        this->condition_variable = condition_variable;
        this->bandwidth_matrix = bandwidth;
        this->num = (uint) condition_variable.size();
        this->kernel_type = KernelType(kernel_type);
        this->distance = distance;
    }

    KernelDensityEstimation(std::vector<std::vector<double>> &condition_variable,
                            std::vector<double> &bandwidth,
                            int kernel_type,
                            bool distance = false) {
        if (bandwidth.size() == 1) {
            this->bandwidth_value = bandwidth[0];
        } else if (bandwidth.size() == pow(condition_variable[0].size(), 2)) {

            uint d = (uint) condition_variable[0].size();
            std::vector<std::vector<double>> bandwidth_matrix(d);
            for (uint i = 0; i < d; i++) {
                bandwidth_matrix[i].resize(d);
            }
            this->bandwidth_matrix = bandwidth_matrix;

            uint k = 0;
            for (uint i = 0; i < condition_variable[0].size(); ++i) {
                for (uint j = 0; j < condition_variable[0].size(); ++j) {
                    this->bandwidth_matrix[i][j] = bandwidth[k++];
                }
            }
        } else {
            this->bandwidth_vector = bandwidth;
        }
        this->condition_variable = condition_variable;
        this->num = (uint) condition_variable.size();
        this->kernel_type = KernelType(kernel_type);
        this->distance = distance;
    }

    KernelDensityEstimation(std::vector<std::vector<double>> &condition_variable,
                            double bandwidth,
                            int kernel_type,
                            bool distance = false) {
        this->condition_variable = condition_variable;
        this->bandwidth_value = bandwidth;
        this->num = (uint) condition_variable.size();
        this->kernel_type = KernelType(kernel_type);
        this->distance = distance;
    }

    ~KernelDensityEstimation() = default;

private:
    uint num;
    bool distance;
    KernelType kernel_type;
    std::vector<std::vector<double>> condition_variable;
    std::vector<std::vector<double>> bandwidth_matrix;
    std::vector<double> bandwidth_vector;
    double bandwidth_value;
    std::vector<std::vector<double>> kernel_density_estimate;

    std::vector<std::vector<double>>
    compute_gaussian_kernel_estimate(std::vector<std::vector<double>> &condition_variable,
                                     std::vector<std::vector<double>> &bandwidth);

    std::vector<std::vector<double>>
    compute_gaussian_kernel_estimate(std::vector<std::vector<double>> &condition_variable,
                                     std::vector<double> &bandwidth);

    std::vector<std::vector<double>>
    compute_gaussian_kernel_estimate(std::vector<std::vector<double>> &condition_variable,
                                     double bandwidth);

    std::vector<std::vector<double>>
    compute_gaussian_kernel_estimate_based_distance(std::vector<std::vector<double>> &distance_matrix,
                                                    std::vector<std::vector<double>> &bandwidth);

    std::vector<std::vector<double>>
    compute_gaussian_kernel_estimate_based_distance(std::vector<std::vector<double>> &distance_matrix,
                                                    std::vector<double> &bandwidth);

    std::vector<std::vector<double>>
    compute_gaussian_kernel_estimate_based_distance(std::vector<std::vector<double>> &distance_matrix,
                                                    double bandwidth);

    std::vector<std::vector<double>>
    compute_rectangle_kernel_estimate_based_distance(std::vector<std::vector<double>> &distance_matrix,
                                                     std::vector<std::vector<double>> &bandwidth);

    std::vector<std::vector<double>>
    compute_rectangle_kernel_estimate_based_distance(std::vector<std::vector<double>> &distance_matrix,
                                                     std::vector<double> &bandwidth);

    std::vector<std::vector<double>>
    compute_rectangle_kernel_estimate_based_distance(std::vector<std::vector<double>> &distance_matrix,
                                                     double bandwidth);

    std::vector<std::vector<double>>
    compute_rectangle_kernel_estimate(std::vector<std::vector<double>> &condition_variable,
                                      std::vector<std::vector<double>> &bandwidth);

    std::vector<std::vector<double >>
    compute_rectangle_kernel_estimate(std::vector<std::vector<double>> &condition_variable,
                                      std::vector<double> &bandwidth);

    std::vector<std::vector<double >>
    compute_rectangle_kernel_estimate(std::vector<std::vector<double>> &condition_variable,
                                      double bandwidth);

    // I think I find a valid method to auto selection strategy for bandwidth even when the distance matrix are computed by Non-Euclidean metric..
    // I briefly describe below:
    // 1. For i = 1, ..., n find the K nearest-neighbor and its associated distance d_{i_K}
    // 2. Take the average of d_{1_K}, ..., d_{n_K}} as the bandwidth
};

#endif //KernelDensityEstimation_H
