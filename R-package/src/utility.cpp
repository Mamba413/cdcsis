//
// Created by JinZhu on 2018/12/26.
//

#include <cmath>
#include <vector>
#include <iostream>
#include <random>
#include <array>
#include "global.h"
#include "utility.h"


/**
 * Compute determinant of a matrix.
 * @param matrix: a N * N matrix
 */
double compute_matrix_determinant(std::vector<std::vector<double>> matrix) {
    uint matrix_dim = (uint) matrix.size();
    uint i, j, k, is, js;
    double f, det, q, d;
    is = 0;
    js = 0;
    f = 1.0;
    det = 1.0;
    for (k = 0; k <= matrix_dim - 2; k++) {
        q = 0.0;
        for (i = k; i <= matrix_dim - 1; i++) {
            for (j = k; j <= matrix_dim - 1; j++) {
                d = fabs(matrix[i][j]);
                if (d > q) {
                    q = d;
                    is = i;
                    js = j;
                }
            }
        }
        if (q + 1.0 == 1.0) {
            det = 0.0;
            return (det);
        }
        if (is != k) {
            f = -f;
            for (j = k; j <= matrix_dim - 1; j++) {
                d = matrix[k][j];
                matrix[k][j] = matrix[is][j];
                matrix[is][j] = d;
            }
        }
        if (js != k) {
            f = -f;
            for (i = k; i <= matrix_dim - 1; i++) {
                d = matrix[i][js];
                matrix[i][js] = matrix[i][k];
                matrix[i][k] = d;
            }
        }
        det = det * matrix[k][k];
        for (i = k + 1; i <= matrix_dim - 1; i++) {
            d = matrix[i][k] / matrix[k][k];
            for (j = k + 1; j <= matrix_dim - 1; j++) {
                matrix[i][j] = matrix[i][j] - d * matrix[k][j];
            }
        }
    }
    det = f * det * matrix[matrix_dim - 1][matrix_dim - 1];
    return (det);
}

/**
 * Computes Inverse of a matrix
 * @param matrix: a N * N matrix
 */
void compute_matrix_inversion(std::vector<std::vector<double>> &matrix) {
    int matrix_dim = (int) matrix.size();
    int *is, *js;
    int i, j, k = 0;
    double d, p;
    is = new int[matrix_dim];
    js = new int[matrix_dim];
    for (k = 0; k <= matrix_dim - 1; k++) {
        d = 0.0;
        for (i = k; i <= matrix_dim - 1; i++) {
            for (j = k; j <= matrix_dim - 1; j++) {
                p = fabs(matrix[i][j]);
                if (p > d) {
                    d = p;
                    is[k] = i;
                    js[k] = j;
                }
            }
        }
        if (d + 1.0 == 1.0) {
            delete[]is;
            delete[]js;
            throw "Error during estimate kernel density";
        }
        if (is[k] != k) {
            for (j = 0; j <= matrix_dim - 1; j++) {
                p = matrix[k][j];
                matrix[k][j] = matrix[is[k]][j];
                matrix[is[k]][j] = p;
            }
        }
        if (js[k] != k) {
            for (i = 0; i <= matrix_dim - 1; i++) {
                p = matrix[i][k];
                matrix[i][k] = matrix[i][js[k]];
                matrix[i][js[k]] = p;
            }
        }
        matrix[k][k] = 1.0 / matrix[k][k];
        for (j = 0; j <= matrix_dim - 1; j++) {
            if (j != k) {
                matrix[k][j] = matrix[k][j] * matrix[k][k];
            }
        }
        for (i = 0; i <= matrix_dim - 1; i++) {
            if (i != k) {
                for (j = 0; j <= matrix_dim - 1; j++) {
                    if (j != k) {
                        matrix[i][j] = matrix[i][j] - matrix[i][k] * matrix[k][j];
                    }
                }
            }
        }
        for (i = 0; i <= matrix_dim - 1; i++) {
            if (i != k) {
                matrix[i][k] = -matrix[i][k] * matrix[k][k];
            }
        }
    }
    for (k = matrix_dim - 1; k >= 0; k--) {
        if (js[k] != k) {
            for (j = 0; j <= matrix_dim - 1; j++) {
                p = matrix[k][j];
                matrix[k][j] = matrix[js[k]][j];
                matrix[js[k]][j] = p;
            }
        }
        if (is[k] != k) {
            for (i = 0; i <= matrix_dim - 1; i++) {
                p = matrix[i][k];
                matrix[i][k] = matrix[i][is[k]];
                matrix[i][is[k]] = p;
            }
        }
    }
    delete[]is;
    delete[]js;
}

/**
 * Matrix multiplication,
 * @param matrix_A: interpret as m rows and n columns matrix
 * @param matrix_B: interpret as n rows and k columns matrix
 */
std::vector<std::vector<double>> compute_matrix_multiplication(std::vector<std::vector<double>> &matrix_A,
                                                               std::vector<std::vector<double>> &matrix_B) {
    uint m = (uint) matrix_A.size();
    uint n = (uint) matrix_B.size();
    uint k = (uint) matrix_B[0].size();
    std::vector<std::vector<double>> matrix_c(m, std::vector<double>(k));

    uint i, j, l;
    for (i = 0; i <= m - 1; i++) {
        for (j = 0; j <= k - 1; j++) {
            matrix_c[i][j] = 0.0;
            for (l = 0; l <= n - 1; l++) {
                matrix_c[i][j] += matrix_A[i][l] * matrix_B[l][j];
            }
        }
    }
    return matrix_c;
}

/**
 * Compute quadratic form product, i.e., x*Matrix*y
 * @param x : 1 * d 2D matrix
 * @param matrix : d * d 2D matrix
 * @param y : d * 1 2D matrix
 */
double quadratic_matrix_multiplication(std::vector<std::vector<double>> &x, std::vector<std::vector<double>> &matrix,
                                       std::vector<std::vector<double>> &y) {
    std::vector<std::vector<double>> quadric_value_matrix(1, std::vector<double>(1));
    std::vector<std::vector<double>> tmp(1, std::vector<double>(x[0].size()));

    tmp = compute_matrix_multiplication(x, matrix);
    quadric_value_matrix = compute_matrix_multiplication(tmp, y);

    return quadric_value_matrix[0][0];
}

/**
 * Compute the summation of 1D std::vector
 * @param vector1 : 1D vector
 * @return summation value
 */
double vector_sum(std::vector<double> &vector1) {
    double sum_value = 0.0;
    for (double value : vector1) {
        sum_value += value;
    }
    return sum_value;
}

double vector_mean(std::vector<double> &vector1) {
    double vector_sum_value = vector_sum(vector1);
    return vector_sum_value / vector1.size();
}

double vector_weight_sum(std::vector<double> &vector1, std::vector<double> &weight) {
    double sum_value = 0.0;
    for (size_t i = 0; i < vector1.size(); ++i) {
        sum_value += vector1[i] * weight[i];
    }
    return sum_value;
}

/**
 * Compute the Euclidean distance matrix
 * @param matrix: interpret x as an matrix
 * @param index: see energy distance for reference
*/
std::vector<std::vector<double>> Euclidean_distance(std::vector<std::vector<double>> &matrix, double index) {
    uint n = (uint) matrix.size();
    uint d = (uint) matrix[0].size();
    std::vector<std::vector<double>> distance_matrix(n, std::vector<double>(n));

    double diff_sum, diff;
    for (uint i = 1; i < n; i++) {
        distance_matrix[i][i] = 0.0;
        for (uint j = 0; j < i; j++) {
            diff_sum = 0.0;
            for (uint k = 0; k < d; k++) {
                diff = matrix[i][k] - matrix[j][k];
                diff_sum += diff * diff;
            }
            distance_matrix[i][j] = distance_matrix[j][i] = pow(sqrt(diff_sum), index);
        }
    }
    return distance_matrix;
}

std::vector<std::vector<double>> Euclidean_distance_byrow(std::vector<std::vector<double>> &matrix, double index) {
    uint n = (uint) matrix[0].size();
    uint d = (uint) matrix.size();
    std::vector<std::vector<double>> distance_matrix(n, std::vector<double>(n));

    double diff_sum, diff;
    for (uint i = 1; i < n; i++) {
        distance_matrix[i][i] = 0.0;
        for (uint j = 0; j < i; j++) {
            diff_sum = 0.0;
            for (uint k = 0; k < d; k++) {
                diff = matrix[k][i] - matrix[k][j];
                diff_sum += diff * diff;
            }
            distance_matrix[i][j] = distance_matrix[j][i] = pow(sqrt(diff_sum), index);
        }
    }
    return distance_matrix;
}

/**
 * Compute Euclidean distance for univariate vector
 * @param matrix
 * @param index
 * @return
 */
std::vector<std::vector<double>> Euclidean_distance(std::vector<double> &matrix, double index) {
    uint n = (uint) matrix.size();
    std::vector<std::vector<double>> distance_matrix(n, std::vector<double>(n));

    double diff_sum;
    for (uint i = 1; i < n; i++) {
        distance_matrix[i][i] = 0.0;
        for (uint j = 0; j < i; j++) {
            diff_sum = fabs(matrix[i] - matrix[j]);
            distance_matrix[i][j] = distance_matrix[j][i] = pow(diff_sum, index);
        }
    }
    return distance_matrix;
}

std::vector<std::vector<double>> vector_to_matrix(std::vector<double> &vector, uint num_row, uint num_col) {

    std::vector<std::vector<double>> matrix(num_row, std::vector<double>(num_col));

    uint k = 0;
    for (uint i = 0; i < num_row; i++) {
        for (uint j = 0; j < num_col; j++) {
            matrix[i][j] = vector[k++];
        }
    }

    return matrix;
}

std::vector<std::vector<double>> weight_distance_anova(std::vector<std::vector<double>> &distance_matrix,
                                                       std::vector<double> &weight) {
    double weight_sum = vector_sum(weight);
    uint num = (uint) distance_matrix.size();

    std::vector<double> marginal_weight_distance(num);
    for (uint i = 0; i < num; ++i) {
        marginal_weight_distance[i] = vector_weight_sum(distance_matrix[i], weight);
    }

    double weight_distance_sum = 0.0;
    for (uint i = 0; i < num; ++i) {
        weight_distance_sum = vector_weight_sum(marginal_weight_distance, weight);
    }
    weight_distance_sum /= weight_sum * weight_sum;

    for (uint i = 0; i < num; ++i) {
        marginal_weight_distance[i] /= weight_sum;
    }

    std::vector<std::vector<double>> weight_distance_anova_table(num, std::vector<double>(num));
    for (uint k = 0; k < num; k++) {
        for (uint j = k; j < num; j++) {
            weight_distance_anova_table[k][j] =
                    distance_matrix[k][j] - marginal_weight_distance[k] - marginal_weight_distance[j] +
                    weight_distance_sum;
            weight_distance_anova_table[j][k] = weight_distance_anova_table[k][j];
        }
    }

    return weight_distance_anova_table;
}

std::vector<double> compositional_transform(std::vector<double> &vector) {
    double vector_sum_value = vector_sum(vector);
    std::vector<double> transformed_vector(vector.size());

    for (size_t i = 0; i < vector.size(); i++) {
        transformed_vector[i] = vector[i] / vector_sum_value;
    }

    return transformed_vector;
}

uint sample_multinomial_distribution(std::vector<double> &probability, std::mt19937_64 &random_number_generator) {
    std::discrete_distribution<uint> multinomial_sampler(probability.begin(), probability.end());
    return multinomial_sampler(random_number_generator);
}

/**
 * Generate unsigned sequence
 * @param start
 * @param end
 * @return <start, start + 1, ..., end>
 */
std::vector<uint> generate_sequence(uint start, uint end) {
    std::vector<uint> sequence;
    for (uint i = start; i <= end; ++i) {
        sequence.push_back(i);
    }

    return sequence;
}