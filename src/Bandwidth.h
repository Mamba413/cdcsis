//
// Created by JinZhu on 2019/1/13.
//

#ifndef CDC_BANDWIDTH_H
#define CDC_BANDWIDTH_H

#include <vector>
#include "utility.h"

class Bandwidth {
public:

    Bandwidth(uint var_nun) {
        for (uint i = 0; i < var_nun; ++i) {
            this->bandwidth_matrix.emplace_back(std::vector<double > (var_nun));
        }
    }

    std::vector<std::vector<double>> bandwidth_selection(std::vector<double>& conditional_vec) {
        double bandwidth = bandwidth_selection_vector_stats_bwnrd0(conditional_vec);
        bandwidth_matrix[0][0] = bandwidth > 0 ? bandwidth : 1.0;
        return bandwidth_matrix;
    };

    std::vector<std::vector<double>> bandwidth_selection(std::vector<std::vector<double>>& conditional_mat) {
        std::vector<double> conditional_vec;
        for (int j = 0; j < conditional_mat[0].size(); ++j) {
            conditional_vec.clear();
            for (std::vector<double> conditional_record : conditional_mat) {
                conditional_vec.push_back(conditional_record[j]);
            }
            bandwidth_matrix[j][j] = bandwidth_selection_vector_stats_bwnrd0(conditional_vec);
        }

        return bandwidth_matrix;
    };

private:
    double bandwidth_selection_vector_stats_bwnrd0(std::vector<double>& conditional_vec) {
        double iqr = quartile_value(conditional_vec, 0.75) - quartile_value(conditional_vec, 0.25);
        double sd = vector_sd(conditional_vec);
        double lo = std::min(sd, iqr/1.34);
        return 0.9 * lo * pow(conditional_vec.size(), -0.2);
    };

    std::vector<std::vector<double>> bandwidth_matrix;
};


#endif //CDC_BANDWIDTH_H
