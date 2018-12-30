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

#endif //SRC_UTILITY_H
