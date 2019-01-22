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

double vector_prod(std::vector<double> &vector1) {
    double sum_value = 1.0;
    for (double value : vector1) {
        sum_value *= value;
    }
    return sum_value;
}

double vector_mean(std::vector<double> &vector1) {
    double vector_sum_value = vector_sum(vector1);
    return vector_sum_value / vector1.size();
}

double vector_sd(std::vector<double> &vector1) {
    double vec_mean = vector_mean(vector1);
    double vec_sd = 0.0;
    for (double vector1_value:vector1) {
        vec_sd += pow(vector1_value - vec_mean, 2);
    }
    vec_sd = vec_sd / (vector1.size() - 1);
    vec_sd = sqrt(vec_sd);
    return vec_sd;
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
 * @param matrix: interpret x as an matrix with size n*d
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

double square_Euclidean_distance(std::vector<double> &vector1, std::vector<double> &vector2) {
    double distance_value = 0.0;
    for (uint i = 0; i < vector1.size(); ++i) {
        distance_value += pow(vector1[i] - vector2[i], 2);
    }
    return distance_value;
}

double weight_square_Euclidean_distance(std::vector<double> &vector1, std::vector<double> &vector2,
                                        std::vector<double> &weight) {
    double distance_value = 0.0;
    for (uint i = 0; i < vector1.size(); ++i) {
        distance_value += weight[i] * pow(vector1[i] - vector2[i], 2);
    }
    return distance_value;
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

/**
 * Compute Delta_{ij}^{X}, j = 1,..., n
 * @param distance_vector a distance vector contains n element
 * @param weight kernel density estimation omega_{k}(Z_l), k = 1,..., n
 *
 * @return sum_{k=1}^{n} omega_{k}(Z_l) delta_{ij, l, n}^{X}, j = 1, ..., n
 *
 * @example: see test case
 */
std::vector<double> compute_weight_delta_x_vector(std::vector<double> &distance_vector,
                                                  std::vector<double> &weight) {
    int num = (int) distance_vector.size();
    std::vector<std::tuple<int, double, double>> dataset;
    for (int i = 0; i < num; ++i) {
        dataset.emplace_back(std::make_tuple(i, distance_vector[i], weight[i]));
    }
    quick_sort_dataset(dataset, 0, (int) (dataset.size() - 1));

    int last_number = 1, update_number = 0;
    double now_value, last_value = std::get<1>(dataset[0]);
    double cumulative_value = std::get<2>(dataset[0]);
    std::vector<double> delta_x_vector(distance_vector.size(), 0.0);
    delta_x_vector[std::get<0>(dataset[0])] = cumulative_value;
    for (int j = 1; j < num; ++j) {
        now_value = std::get<1>(dataset[j]);
        if (now_value == last_value) {
            last_number++;
            cumulative_value += std::get<2>(dataset[j]);
        } else {
            for (int i = 1; i <= last_number; ++i) {
                delta_x_vector[std::get<0>(dataset[j - i])] = cumulative_value;
                update_number++;
            }
            cumulative_value += std::get<2>(dataset[j]);
            last_value = now_value;
            last_number = 1;
        }
    }
    if (update_number < num) {
        for (int i = 1; i <= last_number; ++i) {
            delta_x_vector[std::get<0>(dataset[num - i])] = cumulative_value;
        }
    }

    return delta_x_vector;
}

/**
 * Compute Delta_{ij}^{XY}, j = 1,..., n
 * @param dataset a vector contains n tuple with 4 element.
 * The 1st element: index, 2nd element: distance_x, 3rd element: distance_y, 4th element: kernel density estimation
 * @param delta_y_vector Delta_{ij}^{Y}, j = 1,..., n
 *
 * @return sum_{k=1}^{n} omega_{k}(Z_l) delta_{ij, l, n}^{XY}, j = 1, ..., n
 *
 * @example see test case
 *
 * @note Only if distance_x without ties, this function is right.
 * @todo tackle the case when distance_x is free of ties
 */
std::vector<double> compute_weight_delta_xy_vector(std::vector<double> &delta_y_vector,
                                                   std::vector<double> &distance_x,
                                                   std::vector<double> &distance_y,
                                                   std::vector<double> &weight) {
    int num = (int) distance_x.size(), index;
    std::vector<std::tuple<int, double, double, double >> dataset;
    for (int j = 0; j < num; ++j) {
        dataset.emplace_back(std::make_tuple(j, distance_x[j], distance_y[j], weight[j]));
    }
    quick_sort_dataset(dataset, 0, num - 1);

    std::vector<double> distance_y_vector, weight_vector, weight_sum_count_smaller_number_after_self_vec;
    std::vector<double> weight_delta_xy_vector(dataset.size(), 0.0);

    for (int j = 0; j < num; ++j) {
        distance_y_vector.push_back(std::get<2>(dataset[j]));
        weight_vector.push_back(std::get<3>(dataset[j]));
    }
    weight_sum_count_smaller_number_after_self_vec = weight_sum_count_smaller_number_after_self(distance_y_vector,
                                                                                                weight_vector);

    // Alternative formula, see https://arxiv.org/pdf/1811.03750.pdf, equation 1:
    for (int j = 0; j < num; ++j) {
        index = std::get<0>(dataset[j]);
        weight_delta_xy_vector[index] =
                delta_y_vector[index] - weight_sum_count_smaller_number_after_self_vec[j];
    }

    return weight_delta_xy_vector;
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

std::vector<std::vector<uint>> generate_random_sample_index(uint replication_number,
                                                            std::vector<std::vector<double>> probability_matrix,
                                                            std::mt19937_64 &random_number_generator) {
    std::vector<std::vector<uint>> random_sample_index(replication_number,
                                                       std::vector<uint>(probability_matrix.size()));
    for (uint i = 0; i < probability_matrix.size(); ++i) {
        for (uint j = 0; j < replication_number; ++j) {
            random_sample_index[j][i] = sample_multinomial_distribution(probability_matrix[i], random_number_generator);
        }
    }

    return random_sample_index;
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

/**
 * Solution of Count of Smaller Numbers After Self Problem
 * @param a vector
 *
 * @refitem https://segmentfault.com/a/1190000012866241
 * @example
 * Input : [1, 2, 5, 4, 4, 3]
 * Output : [0, 0, 3, 1, 1, 0]
 */
std::vector<int> countSmaller(std::vector<int> &vector) {
    std::vector<int> right_smaller(vector.size(), 0);
    std::vector<std::pair<int, int>> vec(vector.size());
    for (uint i = 0; i < vector.size(); i++) {
        vec[i] = std::make_pair(i, vector[i]);
    }
    merge_sort(vec, 0, (int) vec.size(), right_smaller);
    return right_smaller;
}

void merge_sort(std::vector<std::pair<int, int>> &vec, int start, int end,
                std::vector<int> &right_smaller) {
    if (end - start <= 1) return;
    int mid = (start + end) >> 1;
    merge_sort(vec, start, mid, right_smaller);
    merge_sort(vec, mid, end, right_smaller);
    merge(vec, start, mid, end, right_smaller);
}

void merge(std::vector<std::pair<int, int>> &vec, int start, int mid, int end,
           std::vector<int> &right_smaller) {
    auto it_start = vec.begin() + start;
    auto it_mid = vec.begin() + mid;
    auto it_end = vec.begin() + end;
    std::vector<std::pair<int, int>> left(it_start, it_mid), right(it_mid, it_end);
    uint left_merged = 0, right_merged = 0, total_merged = 0;
    while (left_merged < left.size() && right_merged < right.size()) {
        if (left[left_merged].second < right[right_merged].second) {
            right_smaller[left[left_merged].first] += right_merged;
            vec[start + total_merged] = left[left_merged];
            ++left_merged;
            ++total_merged;
        } else {
            vec[start + total_merged] = right[right_merged];
            ++right_merged;
            ++total_merged;
        }
    }
    while (left_merged < left.size()) {
        right_smaller[left[left_merged].first] += right_merged;
        vec[start + total_merged] = left[left_merged];
        ++left_merged;
        ++total_merged;
    }
    while (right_merged < right.size()) {
        vec[start + total_merged] = right[right_merged];
        ++right_merged;
        ++total_merged;
    }
}

void quick_sort_dataset(std::vector<std::tuple<int, double, double>> &dataset, int start, int end) {
    if (start >= end)
        return;
    double mid = std::get<1>(dataset[end]);
    int left = start, right = end - 1;
    while (left < right) {
        while (std::get<1>(dataset[left]) < mid && left < right)
            left++;
        while (std::get<1>(dataset[right]) >= mid && left < right)
            right--;
        dataset[left].swap(dataset[right]);
    }
    if (std::get<1>(dataset[left]) >= std::get<1>(dataset[end])) {
        dataset[left].swap(dataset[end]);
    } else {
        left++;
    }
    if (left) {
        quick_sort_dataset(dataset, start, left - 1);
    }
    quick_sort_dataset(dataset, left + 1, end);
}

void quick_sort_dataset(std::vector<std::tuple<int, double, double, double>> &dataset, int start, int end) {
    if (start >= end)
        return;
    double mid = std::get<1>(dataset[end]);
    int left = start, right = end - 1;
    while (left < right) {
        while (std::get<1>(dataset[left]) < mid && left < right)
            left++;
        while (std::get<1>(dataset[right]) >= mid && left < right)
            right--;
        dataset[left].swap(dataset[right]);
    }
    if (std::get<1>(dataset[left]) >= std::get<1>(dataset[end])) {
        dataset[left].swap(dataset[end]);
    } else {
        left++;
    }
    if (left) {
        quick_sort_dataset(dataset, start, left - 1);
    }
    quick_sort_dataset(dataset, left + 1, end);
}

double compute_condition_ball_covariance_crude(std::vector<std::vector<double>> &distance_x,
                                               std::vector<std::vector<double>> &distance_y,
                                               std::vector<std::vector<double>> &kernel_density_estimation) {
    int num = (int) kernel_density_estimation.size();
    std::vector<double> condition_ball_covariance(kernel_density_estimation.size(), 0.0);
    double condition_ball_covariance_stats;

    std::vector<double> kernel_rowsum(kernel_density_estimation.size(), 0.0);
    for (int m = 0; m < num; ++m) {
        kernel_rowsum[m] += vector_sum(kernel_density_estimation[m]);
    }
    std::vector<std::vector<double >> weight_delta_x_matrix(kernel_rowsum.size(),
                                                            std::vector<double>(kernel_rowsum.size()));
    std::vector<std::vector<double >> weight_delta_y_matrix(kernel_rowsum.size(),
                                                            std::vector<double>(kernel_rowsum.size()));
    std::vector<std::vector<double >> weight_delta_xy_matrix(kernel_rowsum.size(),
                                                             std::vector<double>(kernel_rowsum.size()));
    double delta_ijl_xy, delta_ijl_x, delta_ijl_y;
    for (int l = 0; l < num; ++l) {
        for (int i = 0; i < num; ++i) {
            for (int j = 0; j < num; ++j) {
                delta_ijl_xy = delta_ijl_x = delta_ijl_y = 0.0;
                for (int k = 0; k < num; ++k) {
                    delta_ijl_xy += distance_x[i][j] >= distance_x[i][k] && distance_y[i][j] >= distance_y[i][k]
                                    ? kernel_density_estimation[k][l] : 0;
                    delta_ijl_x += distance_x[i][j] >= distance_x[i][k] ? kernel_density_estimation[k][l] : 0;
                    delta_ijl_y += distance_y[i][j] >= distance_y[i][k] ? kernel_density_estimation[k][l] : 0;
                }
                weight_delta_xy_matrix[i][j] = delta_ijl_xy;
                weight_delta_x_matrix[i][j] = delta_ijl_x;
                weight_delta_y_matrix[i][j] = delta_ijl_y;
            }
        }
        for (int i = 0; i < num; ++i) {
            for (int j = 0; j < num; ++j) {
                condition_ball_covariance[l] += pow(
                        weight_delta_xy_matrix[i][j] -
                        weight_delta_x_matrix[i][j] * weight_delta_y_matrix[i][j] / kernel_rowsum[l], 2);
            }
        }
        condition_ball_covariance[l] /= (kernel_rowsum[l] * kernel_rowsum[l] * num * num);
    }
    condition_ball_covariance_stats = vector_mean(condition_ball_covariance);

    return condition_ball_covariance_stats;
}