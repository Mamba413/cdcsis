//
// Created by JinZhu on 2018/12/27.
//

#ifndef SRC_UTILITY_H
#define SRC_UTILITY_H

#include <cmath>
#include <vector>
#include <iostream>
#include "global.h"

double compute_matrix_determinant(std::vector<std::vector<double>> matrix);

void compute_matrix_inversion(std::vector<std::vector<double>>& matrix);

std::vector<std::vector<double>> compute_matrix_multiplication(std::vector<std::vector<double>> &matrix_A,
                                                               std::vector<std::vector<double>> &matrix_B);

double quadratic_matrix_multiplication(std::vector<std::vector<double>> &x, std::vector<std::vector<double>> &matrix,
                                       std::vector<std::vector<double>> &y);

double vector_sum(std::vector<double> &vector1);

double vector_weight_sum(std::vector<double> &vector1, std::vector<double> &weight);

std::vector<std::vector<double>> Euclidean_distance(std::vector<std::vector<double>> &matrix, double index);

std::vector<std::vector<double>> weight_distance_anova(std::vector<std::vector<double>> &distance_matrix,
                                                       std::vector<double> &weight);

#endif //SRC_UTILITY_H
