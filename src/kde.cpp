//
// Created by JinZhu on 2018/12/26.
//

#include "kde.h"
#include "utility.h"

void KernelDensityEstimation::compute_kernel_density_estimate() {
    this->kernel_density_estimate.resize(this->num);
    for (int i = 0; i < this->num; ++i) {
        this->kernel_density_estimate[i].resize(this->num);
    }

    switch (this->kernel_type) {
        case KERNEL_GAUSSIAN:
            this->kernel_density_estimate = KernelDensityEstimation::compute_gaussian_kernel_estimate(
                    this->condition_variable,
                    this->bandwidth);
            break;
        case KERNEL_RECTANGLE:
            this->kernel_density_estimate = KernelDensityEstimation::compute_rectangle_kernel_estimate(
                    this->condition_variable,
                    this->bandwidth);
            break;
    }
}

vector<vector<double>>
KernelDensityEstimation::compute_gaussian_kernel_estimate(vector<vector<double>> &condition_variable,
                                                          vector<vector<double>> &bandwidth) {
    int i, j, k;
    size_t d = bandwidth.size();
    size_t num = condition_variable.size();
    std::vector<std::vector<double>> sigma(d, std::vector<double>(d));
    std::vector<std::vector<double>> kernel_density_estimate(num, std::vector<double>(num));
    std::vector<std::vector<double>> expect1(1, std::vector<double>(d));
    std::vector<std::vector<double>> expect2(d, std::vector<double>(1));

    double density, det, quadric_value;

    if (d == 1) {
        det = (bandwidth[0][0]) * (bandwidth[0][0]);
        sigma[0][0] = 1.0 / bandwidth[0][0];
    } else {
        for (i = 0; i < d; i++) {
            sigma[i] = bandwidth[i];
        }
        det = compute_matrix_determinant(bandwidth);
        compute_matrix_inversion(sigma);
    }

    density = 1.0 / (pow(2 * CDC_PI, d / 2.0) * pow(det, 0.5));

    for (i = 0; i < num; i++) {
        kernel_density_estimate[i][i] = density;
        for (j = 0; j < i; j++) {
            for (k = 0; k < d; k++) {
                expect1[0][k] = condition_variable[i][k] - condition_variable[j][k];
                expect2[k][0] = expect1[0][k];
            }
            quadric_value = quadratic_matrix_multiplication(expect1, sigma, expect2);
            kernel_density_estimate[i][j] = kernel_density_estimate[j][i] = exp(-(quadric_value) / 2.0) * density;
        }
    }

    return (kernel_density_estimate);
}

vector<vector<double>>
KernelDensityEstimation::compute_rectangle_kernel_estimate(vector<vector<double>> &condition_variable,
                                                           vector<vector<double>> &bandwidth) {
    size_t num = condition_variable.size();
    std::vector<std::vector<double>> kernel_density_estimate(num, std::vector<double>(num));
    // TODO:
    return (kernel_density_estimate);
}