//
// Created by JinZhu on 2018/12/27.
//

#include "gtest/gtest.h"
#include "utility.h"

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