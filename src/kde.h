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

    vector<vector<double>> get_kernel_density_estimate() {
        return this->kernel_density_estimate;
    }

    KernelDensityEstimation(vector<vector<double>>& condition_variable, vector<vector<double>>& bandwidth, int kernel_type) {
        this->condition_variable = condition_variable;
        this->bandwidth = bandwidth;
        this->num = (uint) condition_variable.size();
        this->kernel_type = KernelType(kernel_type);
    }

private:
    uint num;
    KernelType kernel_type;
    vector<vector<double>> condition_variable;
    vector<vector<double>> bandwidth;
    vector<vector<double>> kernel_density_estimate;

    vector<vector<double>> compute_gaussian_kernel_estimate(vector<vector<double>> &condition_variable,
                                                            vector<vector<double>> &bandwidth);

    vector<vector<double>> compute_rectangle_kernel_estimate(vector<vector<double>> &condition_variable,
                                                             vector<vector<double>> &bandwidth);
};

#endif //KernelDensityEstimation_H
