//
// Created by JinZhu on 2018/12/26.
//

#ifndef KernelDensityEstimation_H
#define KernelDensityEstimation_H

#include <vector>
#include <iostream>
#include "global.h"

using namespace std;

class KernelDensityEstimation {
public:
    void compute_kernel_density_estimate();

    std::vector<std::vector<double>> get_kernel_density_estimate() {
        return this->kernel_density_estimate;
    }

    KernelDensityEstimation(std::vector<std::vector<double>> &condition_variable,
                            std::vector<std::vector<double>> &bandwidth,
                            int kernel_type) {
        this->condition_variable = condition_variable;
        this->bandwidth = bandwidth;
        this->num = (uint) condition_variable.size();
        this->kernel_type = KernelType(kernel_type);
    }

    KernelDensityEstimation(std::vector<std::vector<double>> &condition_variable,
                            std::vector<double> &bandwidth,
                            int kernel_type) {
        this->condition_variable = condition_variable;
        this->bandwidth_vector = bandwidth;
        this->num = (uint) condition_variable.size();
        this->kernel_type = KernelType(kernel_type);
    }

    KernelDensityEstimation(std::vector<std::vector<double>> &condition_variable,
                            double bandwidth,
                            int kernel_type) {
        this->condition_variable = condition_variable;
        this->bandwidth_value = bandwidth;
        this->num = (uint) condition_variable.size();
        this->kernel_type = KernelType(kernel_type);
    }

    ~KernelDensityEstimation() = default;

private:
    uint num;
    KernelType kernel_type;
    std::vector<std::vector<double>> condition_variable;
    std::vector<std::vector<double>> bandwidth;
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
    compute_rectangle_kernel_estimate(std::vector<std::vector<double>> &condition_variable,
                                      std::vector<std::vector<double>> &bandwidth);
};

#endif //KernelDensityEstimation_H
