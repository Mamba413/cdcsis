//
// Created by JinZhu on 2018/12/27.
//

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
