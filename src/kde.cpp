//
// Created by JinZhu on 2018/12/26.
//

#include "kde.h"
#include "utility.h"

void KernelDensityEstimation::compute_kernel_density_estimate() {
    this->kernel_density_estimate.resize(this->num);
    for (uint i = 0; i < this->num; ++i) {
        this->kernel_density_estimate[i].resize(this->num);
    }

    if (this->distance) {
        switch (this->kernel_type) {
            case KERNEL_GAUSSIAN:
                if (!this->bandwidth_matrix.empty()) {
                    this->kernel_density_estimate = KernelDensityEstimation::compute_gaussian_kernel_estimate_based_distance(
                            this->condition_variable,
                            this->bandwidth_matrix);
                } else if (!this->bandwidth_vector.empty()) {
                    this->kernel_density_estimate = KernelDensityEstimation::compute_gaussian_kernel_estimate_based_distance(
                            this->condition_variable,
                            this->bandwidth_vector);
                } else {
                    this->kernel_density_estimate = KernelDensityEstimation::compute_gaussian_kernel_estimate_based_distance(
                            this->condition_variable,
                            this->bandwidth_value);
                }
                break;
            case KERNEL_RECTANGLE:
                if (!this->bandwidth_matrix.empty()) {
                    this->kernel_density_estimate = KernelDensityEstimation::compute_rectangle_kernel_estimate_based_distance(
                            this->condition_variable,
                            this->bandwidth_matrix);
                } else if (!this->bandwidth_vector.empty()) {
                    this->kernel_density_estimate = KernelDensityEstimation::compute_rectangle_kernel_estimate_based_distance(
                            this->condition_variable,
                            this->bandwidth_vector);
                } else {
                    this->kernel_density_estimate = KernelDensityEstimation::compute_rectangle_kernel_estimate_based_distance(
                            this->condition_variable,
                            this->bandwidth_value);
                }
                break;
        }
    } else {
        switch (this->kernel_type) {
            case KERNEL_GAUSSIAN:
                if (!this->bandwidth_matrix.empty()) {
                    this->kernel_density_estimate = KernelDensityEstimation::compute_gaussian_kernel_estimate(
                            this->condition_variable,
                            this->bandwidth_matrix);
                } else if (!this->bandwidth_vector.empty()) {
                    this->kernel_density_estimate = KernelDensityEstimation::compute_gaussian_kernel_estimate(
                            this->condition_variable,
                            this->bandwidth_vector);
                } else {
                    this->kernel_density_estimate = KernelDensityEstimation::compute_gaussian_kernel_estimate(
                            this->condition_variable,
                            this->bandwidth_value);
                }
                break;
            case KERNEL_RECTANGLE:
                if (!this->bandwidth_matrix.empty()) {
                    this->kernel_density_estimate = KernelDensityEstimation::compute_rectangle_kernel_estimate(
                            this->condition_variable,
                            this->bandwidth_matrix);
                } else if (!this->bandwidth_vector.empty()) {
                    this->kernel_density_estimate = KernelDensityEstimation::compute_rectangle_kernel_estimate(
                            this->condition_variable,
                            this->bandwidth_vector);
                } else {
                    this->kernel_density_estimate = KernelDensityEstimation::compute_rectangle_kernel_estimate(
                            this->condition_variable,
                            this->bandwidth_value);
                }
                break;
        }
    }
}

/**
 * The kernel estimation function is the same as the equation in the section 4.2 The Kernel Function and Bandwidth Selection
 * (Conditional Distance Correlation, 2015, JASA).
 *
 * I think the code provided by original author is wrong and I add follow two line:
 *    sigma = compute_matrix_multiplication(sigma, sigma);
 *
 */
std::vector<std::vector<double>>
KernelDensityEstimation::compute_gaussian_kernel_estimate(std::vector<std::vector<double>> &condition_variable,
                                                          std::vector<std::vector<double>> &bandwidth) {
    uint d = (uint) bandwidth.size();
    uint num = (uint) condition_variable.size();
    std::vector<std::vector<double>> sigma(d, std::vector<double>(d));
    std::vector<std::vector<double>> kernel_density_estimate(num, std::vector<double>(num));
    std::vector<std::vector<double>> expect1(1, std::vector<double>(d));
    std::vector<std::vector<double>> expect2(d, std::vector<double>(1));

    double density, det, quadric_value;
    if (d == 1) {
        det = (bandwidth[0][0]) * (bandwidth[0][0]);
        sigma[0][0] = 1.0 / det;
    } else {
        for (uint i = 0; i < d; i++) {
            sigma[i] = bandwidth[i];
        }
        det = compute_matrix_determinant(bandwidth);
        compute_matrix_inversion(sigma);
        sigma = compute_matrix_multiplication(sigma, sigma);
    }

    density = 1.0 / (pow(2 * CDC_PI, d / 2.0) * pow(det, 0.5));
    for (uint i = 0; i < num; i++) {
        kernel_density_estimate[i][i] = density;
        for (uint j = 0; j < i; j++) {
            for (uint k = 0; k < d; k++) {
                expect1[0][k] = condition_variable[i][k] - condition_variable[j][k];
                expect2[k][0] = expect1[0][k];
            }
            quadric_value = quadratic_matrix_multiplication(expect1, sigma, expect2);
            kernel_density_estimate[i][j] = kernel_density_estimate[j][i] = exp(-(quadric_value) / 2.0) * density;
        }
    }

    return kernel_density_estimate;
}

std::vector<std::vector<double>> KernelDensityEstimation::compute_gaussian_kernel_estimate(
        std::vector<std::vector<double>> &condition_variable, std::vector<double> &bandwidth) {
    uint num = (uint) condition_variable.size();
    uint d = (uint) condition_variable[0].size();
    double det = vector_prod(bandwidth);
    double density = 1.0 / (pow(2 * CDC_PI, d / 2.0) * det);

    std::vector<double> weight;
    for (double width : bandwidth) {
        weight.push_back(1.0 / pow(width, 2));
    }

    std::vector<std::vector<double>> kernel_density_estimate(num, std::vector<double>(num));
    for (uint i = 0; i < num; i++) {
        kernel_density_estimate[i][i] = density;
        for (uint j = 0; j < i; j++) {
            kernel_density_estimate[j][i] = exp(
                    -0.5 * weight_square_Euclidean_distance(condition_variable[i], condition_variable[j], weight));
            kernel_density_estimate[j][i] *= density;
            kernel_density_estimate[i][j] = kernel_density_estimate[j][i];
        }
    }

    return kernel_density_estimate;
}

std::vector<std::vector<double>> KernelDensityEstimation::compute_gaussian_kernel_estimate(
        std::vector<std::vector<double>> &condition_variable, double bandwidth) {
    uint num = (uint) condition_variable.size();
    uint d = (uint) condition_variable[0].size();
    double det = pow(bandwidth, d);
    double density = 1.0 / (pow(2 * CDC_PI, d / 2.0) * sqrt(det));

    std::vector<std::vector<double>> kernel_density_estimate(num, std::vector<double>(num));
    for (uint i = 0; i < num; i++) {
        kernel_density_estimate[i][i] = density;
        for (uint j = 0; j < i; j++) {
            kernel_density_estimate[j][i] = exp(
                    -0.5 * square_Euclidean_distance(condition_variable[i], condition_variable[j]) / pow(bandwidth, 2));
            kernel_density_estimate[j][i] *= density;
            kernel_density_estimate[i][j] = kernel_density_estimate[j][i];
        }
    }

    return kernel_density_estimate;
}

std::vector<std::vector<double>>
KernelDensityEstimation::compute_gaussian_kernel_estimate_based_distance(
        std::vector<std::vector<double>> &distance_matrix,
        std::vector<std::vector<double>> &bandwidth) {
    size_t num = distance_matrix.size();
    std::vector<std::vector<double>> kernel_density_estimate(num, std::vector<double>(num));
    // TODO:
    return (kernel_density_estimate);
}

std::vector<std::vector<double>>
KernelDensityEstimation::compute_gaussian_kernel_estimate_based_distance(
        std::vector<std::vector<double>> &distance_matrix,
        std::vector<double> &bandwidth) {
    size_t num = distance_matrix.size();
    std::vector<std::vector<double>> kernel_density_estimate(num, std::vector<double>(num));
    // TODO:
    return (kernel_density_estimate);
}

std::vector<std::vector<double>>
KernelDensityEstimation::compute_gaussian_kernel_estimate_based_distance(
        std::vector<std::vector<double>> &distance_matrix,
        double bandwidth) {
    size_t num = distance_matrix.size();
    std::vector<std::vector<double>> kernel_density_estimate(num, std::vector<double>(num));
    // TODO:
    return (kernel_density_estimate);
}

std::vector<std::vector<double>>
KernelDensityEstimation::compute_rectangle_kernel_estimate_based_distance(
        std::vector<std::vector<double>> &distance_matrix,
        std::vector<std::vector<double>> &bandwidth) {
    size_t num = distance_matrix.size();
    std::vector<std::vector<double>> kernel_density_estimate(num, std::vector<double>(num));
    // TODO:
    return (kernel_density_estimate);
}

std::vector<std::vector<double>>
KernelDensityEstimation::compute_rectangle_kernel_estimate_based_distance(
        std::vector<std::vector<double>> &distance_matrix,
        std::vector<double> &bandwidth) {
    size_t num = distance_matrix.size();
    std::vector<std::vector<double>> kernel_density_estimate(num, std::vector<double>(num));

    for (size_t i = 0; i < num; i++) {
        for (size_t j = 0; j < num; j++) {
            if (distance_matrix[i][j] <= bandwidth[i]) {
                kernel_density_estimate[i][j] = 1.0;
            } else {
                kernel_density_estimate[i][j] = 0.0;
            }
        }
    }
    return (kernel_density_estimate);
}

std::vector<std::vector<double>>
KernelDensityEstimation::compute_rectangle_kernel_estimate_based_distance(
        std::vector<std::vector<double>> &distance_matrix,
        double bandwidth) {
    size_t num = distance_matrix.size();
    std::vector<std::vector<double>> kernel_density_estimate(num, std::vector<double>(num));
    // TODO:
    return (kernel_density_estimate);
}

std::vector<std::vector<double>>
KernelDensityEstimation::compute_rectangle_kernel_estimate(std::vector<std::vector<double> > &condition_variable,
                                                           std::vector<std::vector<double> > &bandwidth) {
    size_t num = condition_variable.size();
    size_t z_dim = condition_variable[0].size();
    std::vector<std::vector<double> > kernel_density_estimate(num, std::vector<double>(num));

    for (size_t i = 0; i < num; i++) {
        for (size_t j = 0; j < num; j++) {
            kernel_density_estimate[i][j] = 1;
            for (size_t k = 0; k < z_dim; k++) {
                if (abs(condition_variable[i][k] - condition_variable[j][k]) < bandwidth[k][k])
                    kernel_density_estimate[i][j] *= 1.0;
                else
                    kernel_density_estimate[i][j] *= 0;

            }

        }
    }

    return (kernel_density_estimate);
}

std::vector<std::vector<double>>
KernelDensityEstimation::compute_rectangle_kernel_estimate(std::vector<std::vector<double>> &condition_variable,
                                                           std::vector<double> &bandwidth) {
    size_t num = condition_variable.size();
    size_t z_dim = condition_variable[0].size();
    std::vector<std::vector<double>> kernel_density_estimate(num, std::vector<double>(num));

    for (size_t i = 0; i < num; i++) {
        for (size_t j = 0; j < num; j++) {
            kernel_density_estimate[i][j] = 1;
            for (size_t k = 0; k < z_dim; k++) {
                if (abs(condition_variable[i][k] - condition_variable[j][k]) < bandwidth[k]) {
                    kernel_density_estimate[i][j] *= 1.0;
                } else {
                    kernel_density_estimate[i][j] *= 0;
                }
            }
        }
    }
    return (kernel_density_estimate);
}

std::vector<std::vector<double>>
KernelDensityEstimation::compute_rectangle_kernel_estimate(std::vector<std::vector<double> > &condition_variable,
                                                           double bandwidth) {
    size_t num = condition_variable.size();
    size_t z_dim = condition_variable[0].size();
    std::vector<std::vector<double> > kernel_density_estimate(num, std::vector<double>(num));

    for (size_t i = 0; i < num; i++) {
        for (size_t j = 0; j < num; j++) {
            kernel_density_estimate[i][j] = 1;
            for (size_t k = 0; k < z_dim; k++) {
                if (abs(condition_variable[i][k] - condition_variable[j][k]) < bandwidth)
                    kernel_density_estimate[i][j] *= 1.0;
                else
                    kernel_density_estimate[i][j] *= 0;

            }

        }
    }
    return (kernel_density_estimate);
}

void KernelDensityEstimation::setBandwidth_matrix(const vector<vector<double>> &bandwidth_matrix) {
    KernelDensityEstimation::bandwidth_matrix = bandwidth_matrix;
}
