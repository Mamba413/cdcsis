//
// Created by JinZhu on 2018/12/27.
//

#ifndef SRC_UTILITY_H
#define SRC_UTILITY_H

#include <cmath>
#include <vector>
#include <array>
#include <random>
#include <iostream>
#include <algorithm>
#include "global.h"

double compute_matrix_determinant(std::vector<std::vector<double>> matrix);

void compute_matrix_inversion(std::vector<std::vector<double>> &matrix);

std::vector<std::vector<double>> compute_matrix_multiplication(std::vector<std::vector<double>> &matrix_A,
                                                               std::vector<std::vector<double>> &matrix_B);

double quadratic_matrix_multiplication(std::vector<std::vector<double>> &x, std::vector<std::vector<double>> &matrix,
                                       std::vector<std::vector<double>> &y);

template<typename T>
std::vector<std::vector<T>> extract_submatrix_byrow(std::vector<std::vector<T>> &matrix,
                                                    std::vector<uint> &index) {
    std::vector<std::vector<T>> submatrix(index.size(), std::vector<T>(matrix[0].size()));
    uint k = 0;
    for (uint index_value : index) {
        submatrix[k++] = matrix[index_value];
    }

    return submatrix;
}

template<typename T>
std::vector<std::vector<T>> extract_submatrix_bycol(std::vector<std::vector<T>> &matrix,
                                                    std::vector<uint> &index) {
    std::vector<std::vector<T>> submatrix(matrix.size(), std::vector<T>(index.size()));
    uint k = 0;
    for (uint index_value : index) {
        for (uint i = 0; i < matrix.size(); ++i) {
            submatrix[i][k] = matrix[i][index_value];
        }
        k++;
    }

    return submatrix;
}

double vector_sum(std::vector<double> &vector1);

double vector_mean(std::vector<double> &vector1);

double vector_sd(std::vector<double> &vector1);

double vector_weight_sum(std::vector<double> &vector1, std::vector<double> &weight);

std::vector<std::vector<double>> Euclidean_distance(std::vector<std::vector<double>> &matrix, double index);

std::vector<std::vector<double>> Euclidean_distance_byrow(std::vector<std::vector<double>> &matrix, double index);

std::vector<std::vector<double>> Euclidean_distance(std::vector<double> &matrix, double index);

std::vector<std::vector<double>> vector_to_matrix(std::vector<double> &vector, uint num_row, uint num_col);

std::vector<std::vector<double>> weight_distance_anova(std::vector<std::vector<double>> &distance_matrix,
                                                       std::vector<double> &weight);

/**
 * Rearrange N*N matrix
 * @param matrix
 * @param rearrange_index
 * @return
 */
template<typename T>
std::vector<std::vector<T>> rearrange_matrix(std::vector<std::vector<T>> &matrix,
                                             std::vector<uint> &rearrange_index) {
    std::vector<std::vector<T>> new_matrix1(rearrange_index.size(), std::vector<T>(rearrange_index.size()));
    std::vector<std::vector<T>> new_matrix2(rearrange_index.size(), std::vector<T>(rearrange_index.size()));
    uint k = 0;
    for (uint index : rearrange_index) {
        new_matrix1[k] = matrix[index];
        new_matrix2[k] = matrix[index];
        k++;
    }
    k = 0;
    for (uint index : rearrange_index) {
        for (uint i = 0; i < matrix.size(); ++i) {
            new_matrix1[i][k] = new_matrix2[i][index];
        }
        k++;
    }
    return new_matrix1;
}

std::vector<double> compositional_transform(std::vector<double> &vector);

uint sample_multinomial_distribution(std::vector<double> &probability, std::mt19937_64 &random_number_generator);

std::vector<uint> generate_sequence(uint start, uint end);

template<typename T>
T quartile_value(std::vector<T> vector, double quartile) {
    int q_index = (int) (vector.size() * quartile + 0.5) - 1;
    std::nth_element(vector.begin(), vector.begin() + q_index, vector.end());
    return vector[q_index];
}

void merge(std::vector<std::pair<int, int>> &vec, int start, int mid, int end,
           std::vector<int> &right_smaller);

void merge_sort(std::vector<std::pair<int, int>> &vec, int start, int end,
                std::vector<int> &right_smaller);

template<typename T>
void weight_sum_merge(std::vector<std::pair<int, T>> &vec, std::vector<std::pair<int, double>> &weight_vec,
                      int start, int mid, int end,
                      std::vector<double> &right_smaller_weight_sum) {
    std::vector<std::pair<int, int>> left(vec.begin() + start, vec.begin() + mid);
    std::vector<std::pair<int, int>> right(vec.begin() + mid, vec.begin() + end);
    std::vector<std::pair<int, int>> weight_left(weight_vec.begin() + start, weight_vec.begin() + mid);
    std::vector<std::pair<int, int>> weight_right(weight_vec.begin() + mid, weight_vec.begin() + end);
    int left_merged = 0, right_merged = 0, right_merged_tmp = 0, total_merged = 0;
    while (left_merged < left.size() && right_merged < right.size()) {
        if (left[left_merged].second < right[right_merged].second) {
            vec[start + total_merged] = left[left_merged];
            weight_vec[start + total_merged] = weight_left[left_merged];

            right_merged_tmp = right_merged;
            while (right_merged_tmp > 0) {
                right_smaller_weight_sum[left[left_merged].first] += weight_right[right_merged_tmp - 1].second;
                right_merged_tmp--;
            }

            ++left_merged;
            ++total_merged;
        } else {
            vec[start + total_merged] = right[right_merged];
            weight_vec[start + total_merged] = weight_right[right_merged];
            ++right_merged;
            ++total_merged;
        }
    }
    while (left_merged < left.size()) {
        vec[start + total_merged] = left[left_merged];
        weight_vec[start + total_merged] = weight_left[left_merged];

        right_merged_tmp = right_merged;
        while (right_merged_tmp > 0) {
            right_smaller_weight_sum[left[left_merged].first] += weight_right[right_merged_tmp - 1].second;
            right_merged_tmp--;
        }

        ++left_merged;
        ++total_merged;
    }
    while (right_merged < right.size()) {
        vec[start + total_merged] = right[right_merged];
        weight_vec[start + total_merged] = weight_right[right_merged];
        ++right_merged;
        ++total_merged;
    }
}

template<typename T>
void
weight_sum_merge_sort(std::vector<std::pair<int, T>> &vec, std::vector<std::pair<int, double>> &weight_vec,
                      int start, int end,
                      std::vector<double> &right_smaller_weight_sum) {
    if (end - start <= 1) return;
    int mid = (start + end) >> 1;
    weight_sum_merge_sort(vec, weight_vec, start, mid, right_smaller_weight_sum);
    weight_sum_merge_sort(vec, weight_vec, mid, end, right_smaller_weight_sum);
    weight_sum_merge(vec, weight_vec, start, mid, end, right_smaller_weight_sum);
}

/**
 * A summation version of "a typical “count of smaller numbers after self" problem
 * @refitem: https://leetcode.com/problems/count-of-smaller-numbers-after-self/
 *
 * @example:
 * Input:
 * vector = [1.0, 3.0, 5,0, 4.0, 2.0, 2.0]
 * weight = [3.0, 2.0, 4.0, 1.0, 2.0, 5.0]
 * Output:
 * [0.0, 7.0, 8.0, 7.0, 5.0, 0.0]
 */
template<typename T>
std::vector<T> weight_sum_count_smaller_number_after_self(std::vector<T> &vector, std::vector<double> &weight) {
    std::vector<double> right_smaller_weight_sum(vector.size(), 0);
    std::vector<std::pair<int, T>> vec(vector.size());
    std::vector<std::pair<int, double>> weight_vec(weight.size());
    for (int i = 0; i < vector.size(); i++) {
        vec[i] = std::make_pair(i, vector[i]);
        weight_vec[i] = std::make_pair(i, weight[i]);
    }
    weight_sum_merge_sort(vec, weight_vec, 0, (int) vec.size(), right_smaller_weight_sum);
    return right_smaller_weight_sum;
}

std::vector<int> countSmaller(std::vector<int> &nums);

#endif //SRC_UTILITY_H
