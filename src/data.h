//
// Created by JinZhu on 2019/1/12.
//

#ifndef CDC_DATA_H
#define CDC_DATA_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "global.h"
#include "utility.h"

class Data {
public:
    Data(std::string &input_file, std::vector<uint> x_index, std::vector<uint> y_index, std::vector<uint> z_index,
         bool distance) {
        this->num = count_num(input_file);
        this->dataset = read_dataset(input_file);
        if (distance) {

        } else {
            this->var_num = (uint) dataset.size() / num;
            std::vector<std::vector<double>> data_tmp;
            data_tmp = extract_data(this->dataset, x_index, this->num, this->var_num);
            this->distance_x = Euclidean_distance(data_tmp, 1);
            data_tmp = extract_data(this->dataset, y_index, this->num, this->var_num);
            this->distance_y = Euclidean_distance(data_tmp, 1);
            data_tmp.clear();
            this->data_z = extract_data(this->dataset, z_index, this->num, this->var_num);
            this->conditional_var_num = (uint) z_index.size();
        }
        dataset.clear();
    };

    uint getVar_num() const {
        return var_num;
    }

    uint getConditional_var_num() const {
        return conditional_var_num;
    }

    const std::vector<std::vector<double>> &getData_z() const {
        return data_z;
    }

    const std::vector<std::vector<double>> &getDistance_x() const {
        return distance_x;
    }

    const std::vector<std::vector<double>> &getDistance_y() const {
        return distance_y;
    }


private:
    uint num;
    uint var_num;
    uint conditional_var_num;
    std::vector<double> dataset;
    std::vector<std::vector<double>> data_z;
    std::vector<std::vector<double>> distance_x;
    std::vector<std::vector<double>> distance_y;

    uint count_num(std::string &input_file) {
        std::ifstream read_file;
        uint num = 0;
        std::string tmp;
        read_file.open(input_file, std::ios::in);
        if (read_file.fail()) {
            throw "Cannot read file!";
        } else {
            while (getline(read_file, tmp, '\n')) {
                num++;
            }
            read_file.close();
            return num;
        }
    };

    std::vector<double> read_dataset(std::string &input_file) {
        std::vector<double> dataset;
        double value;
        std::ifstream read_file;
        read_file.open(input_file, std::ios::in);
        while (!read_file.eof()) {
            read_file >> value;
            dataset.push_back(value);
        }
        read_file.close();
        dataset.pop_back();
        return dataset;
    };

    std::vector<std::vector<double>>
    extract_data(std::vector<double> &dataset, std::vector<uint> &index_vec, uint num, uint var_num) {
        std::vector<std::vector<double>> sub_data;
        std::vector<double> sub_vec (var_num);
        for (uint j = 0; j < num; ++j) {
            sub_vec.clear();
            for (uint column_index : index_vec) {
                sub_vec.push_back(dataset[j * var_num + column_index]);
            }
            sub_data.push_back(sub_vec);
        }

        return sub_data;
    }

    std::vector<double> extract_data(std::vector<double> &dataset, uint index, uint num, uint var_num) {
        std::vector<double> sub_vec;
        for (uint j = 0; j < num; ++j) {
            sub_vec.push_back(dataset[j * var_num + index]);
        }

        return sub_vec;
    }
};


#endif //CDC_DATA_H
