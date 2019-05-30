//
// Created by JinZhu on 2018/12/27.
//

#include <cmath>
#include <vector>
#include <iostream>
#include <random>
#include "gtest/gtest.h"
#include "utility.h"
#include "kde.h"
#include "global.h"

TEST(utility, compute_matrix_determinant) {
    std::vector<std::vector<double>> test_vector(2, std::vector<double>(2));
    test_vector[0][0] = 1;
    test_vector[0][1] = 0;
    test_vector[1][0] = 0;
    test_vector[1][1] = 1;
    EXPECT_EQ(compute_matrix_determinant(test_vector), 1.0);

    test_vector[0][0] = 4;
    test_vector[0][1] = 2;
    test_vector[1][0] = 1;
    test_vector[1][1] = 4;
    EXPECT_EQ(compute_matrix_determinant(test_vector), 14.0);
}

TEST(utility, compute_matrix_inversion) {
    std::vector<std::vector<double>> test_vector(2, std::vector<double>(2));
    test_vector[0][0] = 4;
    test_vector[0][1] = 2;
    test_vector[1][0] = 1;
    test_vector[1][1] = 4;

    compute_matrix_inversion(test_vector);
    double abs_error = 0.000001;
    EXPECT_NEAR(test_vector[0][0], 0.28571429, abs_error);
    EXPECT_NEAR(test_vector[0][1], -0.1428571, abs_error);
    EXPECT_NEAR(test_vector[1][0], -0.07142857, abs_error);
    EXPECT_NEAR(test_vector[1][1], 0.2857143, abs_error);
}

TEST(utility, compute_matrix_multiplication) {
    std::vector<std::vector<double>> matrix1(2, std::vector<double>(2));
    matrix1[0][0] = 0.68740;
    matrix1[0][1] = 0.56718;
    matrix1[1][0] = 0.30923;
    matrix1[1][1] = 1.67582;

    std::vector<std::vector<double>> matrix2(2, std::vector<double>(3));
    matrix2[0][0] = -0.53381;
    matrix2[0][1] = -1.53642;
    matrix2[0][2] = -1.60320;
    matrix2[1][0] = 1.06933;
    matrix2[1][1] = -1.26774;
    matrix2[1][2] = -1.41094;

    std::vector<std::vector<double>> matrix3(2, std::vector<double>(3));
    matrix3 = compute_matrix_multiplication(matrix1, matrix2);

    double abs_error = 0.000001;
    EXPECT_NEAR(0.2395616, matrix3[0][0], abs_error);
    EXPECT_NEAR(-1.775172, matrix3[0][1], abs_error);
    EXPECT_NEAR(-1.902297, matrix3[0][2], abs_error);
    EXPECT_NEAR(1.6269345, matrix3[1][0], abs_error);
    EXPECT_NEAR(-2.599611, matrix3[1][1], abs_error);
    EXPECT_NEAR(-2.860239, matrix3[1][2], abs_error);
}

TEST(utility, rearrange_matrix) {
    std::vector<std::vector<double>> matrix(3, std::vector<double>(3));
    matrix[0][0] = 0.0;
    matrix[0][1] = 1.0;
    matrix[0][2] = 2.0;
    matrix[1][0] = 1.0;
    matrix[1][1] = 0.0;
    matrix[1][2] = 3.0;
    matrix[2][0] = 2.0;
    matrix[2][1] = 3.0;
    matrix[2][2] = 0.0;

    std::vector<uint> rearrange_index = {1, 2, 0};

    std::vector<std::vector<double>> new_matrix(3, std::vector<double>(3));

    new_matrix = rearrange_matrix(matrix, rearrange_index);

    ASSERT_EQ(new_matrix[0][0], 0.0);
    ASSERT_EQ(new_matrix[1][1], 0.0);
    ASSERT_EQ(new_matrix[2][2], 0.0);
}

TEST(utility, compute_weight_delta_x_vector) {
    std::vector<double> distance_vector = {3.0, 2.0, 4.0, 5.0, 1.0};
    std::vector<double> weight = {1.0, 2.0, 3.0, 4.0, 5.0};
    std::vector<double> true_delta_x = {8.0, 7.0, 11.0, 15.0, 5.0};
    std::vector<double> delta_x = compute_weight_delta_x_vector(distance_vector, weight);
    for (int i = 0; i < true_delta_x.size(); ++i) {
        EXPECT_EQ(true_delta_x[i], delta_x[i]);
    }

    distance_vector = {2.0, 2.0, 2.0, 2.0, 2.0};
    weight = {1.0, 2.0, 3.0, 4.0, 5.0};
    true_delta_x = {15.0, 15.0, 15.0, 15.0, 15.0};
    delta_x = compute_weight_delta_x_vector(distance_vector, weight);
    for (int i = 0; i < true_delta_x.size(); ++i) {
        EXPECT_EQ(true_delta_x[i], delta_x[i]);
    }
}

TEST(utility, compute_weight_delta_xy_vector) {
    std::vector<double> distance_x = {1.0, 2.0, 3.0, 4.0, 5.0};
    std::vector<double> distance_y = {5.0, 4.0, 3.0, 2.0, 1.0};
    std::vector<double> weight = {2.0, 1.0, 2.0, 1.0, 3.0};
    std::vector<double> delta_y = compute_weight_delta_x_vector(distance_y, weight);
    std::vector<double> delta_xy = compute_weight_delta_xy_vector(delta_y, distance_x, distance_y, weight);

    std::vector<double> delta_xy_true(distance_x.size(), 0.0);
    for (int j = 0; j < distance_x.size(); ++j) {
        for (int k = 0; k < distance_x.size(); ++k) {
            delta_xy_true[j] += distance_x[j] >= distance_x[k] && distance_y[j] >= distance_y[k] ? weight[k] : 0;
        }
    }

    for (int i = 0; i < distance_x.size(); ++i) {
        EXPECT_EQ(delta_xy_true[i], delta_xy[i]);
    }
}

TEST(utility, compute_weight_delta_xy_vector_ties) {
    std::vector<double> distance_x = {1.0, 3.0, 3.0, 4.0, 5.0};
    std::vector<double> distance_y = {5.0, 4.0, 3.0, 2.0, 2.0};
    std::vector<double> weight = {2.0, 1.0, 2.0, 1.0, 3.0};

    std::vector<double> delta_y, delta_xy, delta_xy_true;
    delta_y = compute_weight_delta_x_vector(distance_y, weight);
    delta_xy = compute_weight_delta_xy_vector_ties(delta_y, distance_x, distance_y, weight);
    delta_xy_true = compute_weight_delta_xy_vector_crude(distance_x, distance_y, weight);
    for (int i = 0; i < distance_x.size(); ++i) {
        EXPECT_EQ(delta_xy_true[i], delta_xy[i]);
    }

    distance_x[2] = 4.0;
    delta_y = compute_weight_delta_x_vector(distance_y, weight);
    delta_xy = compute_weight_delta_xy_vector_ties(delta_y, distance_x, distance_y, weight);
    delta_xy_true = compute_weight_delta_xy_vector_crude(distance_x, distance_y, weight);
    for (int i = 0; i < distance_x.size(); ++i) {
        EXPECT_EQ(delta_xy_true[i], delta_xy[i]);
    }

    std::vector<double> distance_x1 = {0, 0, 0, 0, 1, 1, 1, 1};
    std::vector<double> distance_y1 = {0, 0, 1, 1, 0, 0, 1, 1};
    std::vector<double> weight_1 = {1, 2, 3, 4, 5, 6, 7, 8};
    delta_y.clear();
    delta_xy.clear();
    delta_xy_true.clear();
    delta_y = compute_weight_delta_x_vector(distance_y1, weight_1);
    delta_xy = compute_weight_delta_xy_vector_ties(delta_y, distance_x1, distance_y1, weight_1);
    delta_xy_true = compute_weight_delta_xy_vector_crude(distance_x1, distance_y1, weight_1);
    for (int i = 0; i < distance_x.size(); ++i) {
        EXPECT_EQ(delta_xy_true[i], delta_xy[i]);
    }
}

TEST(utility, sample_multinomial_distribution) {
    std::vector<double> probability = {0.1, 0.2, 0.3, 0.4};
    std::mt19937_64 random_number_generator;
    std::random_device random_device;
    random_number_generator.seed(random_device());

    uint my_sample = sample_multinomial_distribution(probability, random_number_generator);
    int get_target = 0;
    for (int i = 0; i < probability.size(); ++i) {
        get_target += my_sample == i;
    }
    EXPECT_EQ(1, get_target);

    probability.push_back(0.5);
    probability.push_back(0.6);
    probability.push_back(0.7);
    probability.push_back(0.8);
    random_number_generator.seed(1);
    uint my_sample1, my_sample2;
    for (int j = 0; j < 50; ++j) {
        my_sample1 = sample_multinomial_distribution(probability, random_number_generator);
        random_number_generator.seed(1);
        my_sample2 = sample_multinomial_distribution(probability, random_number_generator);
        EXPECT_EQ(my_sample1, my_sample2);
    }
}

TEST(utility, generate_sequence) {
    std::vector<uint> seq1 = generate_sequence(0, 0);
    EXPECT_EQ(seq1.size(), 1);
    EXPECT_EQ(seq1[0], 0);

    std::vector<uint> seq2 = generate_sequence(1, 3);
    EXPECT_EQ(seq2.size(), 3);
    EXPECT_EQ(seq2[0], 1);
    EXPECT_EQ(seq2[1], 2);
    EXPECT_EQ(seq2[2], 3);
}

TEST(utility, quartile_value) {
    std::vector<int> vector;
    for (int i = 1; i <= 10; i++) vector.push_back(i);
    int q_value = quartile_value(vector, 0.6);
    EXPECT_EQ(q_value, 6);
    vector.clear();

    std::vector<double> vector1;
    for (int i = 1; i <= 10; i++) vector1.push_back(static_cast<double>(i));
    double q_value1 = quartile_value(vector1, 0.4);
    EXPECT_EQ(q_value1, 4.0);
    vector.clear();
}

TEST(utility, countSmaller) {
    std::vector<int> data = {1, 2, 5, 6, 4, 3};
    std::vector<int> true_number = {0, 0, 2, 2, 1, 0};
    std::vector<int> number = countSmaller(data);
    for (int i = 0; i < true_number.size(); ++i) {
        EXPECT_EQ(number[i], true_number[i]);
    }

    data = {2, 2, 2, 2, 2, 2};
    true_number = {5, 4, 3, 2, 1, 0};
    number = countSmaller(data);
    for (int i = 0; i < true_number.size(); ++i) {
        EXPECT_EQ(number[i], true_number[i]);
    }
}

TEST(utility, weight_sum_count_smaller_number_after_self) {
    std::vector<double> data = {1.0, 3.0, 5.0, 4.0, 6.0, 2.0};
    std::vector<double> weight = {3.0, 2.0, 4.0, 1.0, 2.0, 5.0};
    std::vector<double> true_weight_sum = {0.0, 5.0, 6.0, 5.0, 5.0, 0.0};
    std::vector<double> weight_sum = weight_sum_count_smaller_number_after_self(data, weight);
    for (int i = 0; i < true_weight_sum.size(); ++i) {
        EXPECT_EQ(weight_sum[i], true_weight_sum[i]);
    }

    data = {2.0, 2.0, 2.0, 2.0, 2.0, 2.0};
    weight = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
    true_weight_sum = {20.0, 18.0, 15.0, 11.0, 6.0, 0.0};
    weight_sum = weight_sum_count_smaller_number_after_self(data, weight);
    for (int i = 0; i < true_weight_sum.size(); ++i) {
        EXPECT_EQ(weight_sum[i], true_weight_sum[i]);
    }
}

TEST(utility, quick_sort_dataset) {
    std::vector<double> data = {0.2, 0.4, 0.2, 2.0, 3.0, 0.2, 0.6};
    std::vector<double> weight = {1.0, 2.0, 2.0, 3.0, 4.0, 5.0, 0.2};
    std::vector<int> index = {1, 2, 3, 4, 5, 6, 7};
    std::vector<std::tuple<int, double, double >> dataset;
    for (int i = 0; i < index.size(); ++i) {
        dataset.emplace_back(make_tuple(index[i], data[i], weight[i]));
    }
    quick_sort_dataset(dataset);

    std::vector<double> true_value = {0.2, 0.2, 0.2, 0.4, 0.6, 2.0, 3.0};
    for (int j = 0; j < dataset.size(); ++j) {
        EXPECT_EQ(true_value[j], std::get<1>(dataset[j]));
    }
}

TEST(utility, quick_sort_dataset2) {
    std::vector<double> v1 = {0.2, 0.2, 0.2, 2.0, 3.0, 0.6};
    std::vector<double> v2 = {0.2, 0.4, 0.3, 2.0, 3.0, 0.2};
    std::vector<double> weight = {1.0, 2.0, 2.0, 3.0, 4.0, 5.0};
    std::vector<int> index = {1, 2, 3, 4, 5, 6};
    std::vector<std::tuple<int, double, double, double>> dataset;
    for (int i = 0; i < index.size(); ++i) {
        dataset.emplace_back(make_tuple(index[i], v1[i], v2[i], weight[i]));
    }
    quick_sort_dataset2(dataset);

    std::vector<double> true_v1 = {0.2, 0.2, 0.2, 0.6, 2.0, 3.0};
    std::vector<double> true_v2 = {0.2, 0.3, 0.4, 0.2, 2.0, 3.0};
    for (int j = 0; j < dataset.size(); ++j) {
        EXPECT_EQ(true_v1[j], std::get<1>(dataset[j]));
        EXPECT_EQ(true_v2[j], std::get<2>(dataset[j]));
    }
}

TEST(kde, compute_gaussian_kernel_estimate) {
    std::vector<std::vector<double>> conditional_variable(5, std::vector<double>(2));
    conditional_variable[0][0] = -0.6264538;
    conditional_variable[0][1] = -0.8204684;
    conditional_variable[1][0] = 0.1836433;
    conditional_variable[1][1] = 0.4874291;
    conditional_variable[2][0] = -0.8356286;
    conditional_variable[2][1] = 0.7383247;
    conditional_variable[3][0] = 1.5952808;
    conditional_variable[3][1] = 0.5757814;
    conditional_variable[4][0] = 0.3295078;
    conditional_variable[4][1] = -0.3053884;

    std::vector<std::vector<double>> bandwidth(2, std::vector<double>(2));
    bandwidth[0][0] = 1.0;
    bandwidth[0][1] = 0.0;
    bandwidth[1][0] = 0.0;
    bandwidth[1][1] = 1.0;

    KernelDensityEstimation kernelDensityEstimation = KernelDensityEstimation(conditional_variable, bandwidth, 1);
    kernelDensityEstimation.compute_kernel_density_estimate();

    double abs_error = 0.000001;
    std::vector<std::vector<double>> kernel(5, std::vector<double>(5));
    // Note: the kernel true can be obtain by running following R code:
    // set.seed(1)
    // xx <- cbind(rnorm(5), rnorm(5))
    // set.seed(1)
    // x <- rnorm(5)
    // y <- rnorm(5)
    // debug(cdcor)
    // cdcor(x, y, z = xx, width = diag(1, 2))
    // re[["k"]]
    std::vector<double> kernel_true = {0.159154943, 0.048737642, 0.046204778, 0.005088938, 0.088260550,
                                       0.048737642, 0.159154943, 0.091738878, 0.058534420, 0.115003234,
                                       0.046204778, 0.091738878, 0.159154943, 0.008182582, 0.046825986,
                                       0.005088938, 0.058534420, 0.008182582, 0.159154943, 0.048451218,
                                       0.088260550, 0.115003234, 0.046825986, 0.048451218, 0.159154943};
    kernel = kernelDensityEstimation.get_kernel_density_estimate();
    int k = 0;
    for (int i = 0; i < kernel.size(); ++i) {
        for (int j = 0; j < kernel[0].size(); ++j) {
            EXPECT_NEAR(kernel[i][j], kernel_true[k++], abs_error);
        }
    }
}
