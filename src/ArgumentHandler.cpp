//
// Created by JinZhu on 2019/1/12.
//

#include "ArgumentHandler.h"

#include <sstream>
#include <iostream>
#include <stdexcept>
#include <map>

ArgumentHandler::ArgumentHandler(int argc, char **argv) :
        testtype(1), num_permutation(99), distance(false) {
    this->argc = argc;
    this->argv = argv;
}

void ArgumentHandler::displayHelp() {
    std::cout << "Usage: " << std::endl;
    std::cout << "    " << argv[0] << " [options]" << std::endl;
    std::cout << std::endl;

    std::cout << "Options:" << std::endl;
    std::cout << "    " << "--file FILE                   Filename of input data. Only numerical values are supported."
              << std::endl;
    std::cout << "    " << "--xindex                      index of x" << std::endl;
    std::cout << "    " << "--yindex                      index of y" << std::endl;
    std::cout << "    " << "--zindex                      index of z" << std::endl;
    std::cout << "    " << "--testtype TYPE               Set test type to:" << std::endl;
    std::cout << "    " << "                              TYPE = 1: Unbiased." << std::endl;
    std::cout << "    " << "                              TYPE = 2: Biased." << std::endl;
    std::cout << "    " << "--num_permutation             Local Bootstrap replication." << std::endl;
    std::cout << "    " << "--distance                    Input file is distance matrix." << std::endl;
    std::cout << std::endl;

    std::cout << "See README file for details and examples." << std::endl;
}

void ArgumentHandler::checkArguments() {

    if (filename.empty()) {
        throw std::runtime_error("Please specify an input filename with '--file'. See '--help' for details.");
    }

    if (x_index.empty()) {
        throw std::runtime_error("Please specify x-index with '--xindex'");
    }

    if (y_index.empty()) {
        throw std::runtime_error("Please specify y-index with '--yindex'");
    }

    if (z_index.empty()) {
        throw std::runtime_error("Please specify z-index with '--zindex'");
    }
}

uint ArgumentHandler::processArguments() {
    int option;
    int optional_index = 1;

    std::map<std::string, int> option_map;
    option_map.insert(std::pair<std::string, int>("file", 1));
    option_map.insert(std::pair<std::string, int>("xindex", 2));
    option_map.insert(std::pair<std::string, int>("yindex", 3));
    option_map.insert(std::pair<std::string, int>("zindex", 4));
    option_map.insert(std::pair<std::string, int>("testtype", 5));
    option_map.insert(std::pair<std::string, int>("num_permutation", 6));
    option_map.insert(std::pair<std::string, int>("distance", 7));

    std::string option_name;
    while (optional_index < argc) {

        option_name = argv[optional_index++];
        option_name.erase(0, 2);

        auto i = option_map.find(option_name);
        if (i == option_map.end()) {
            throw std::runtime_error("invalid input arguments!");
        } else {
            option = i->second;
        }

        switch (option) {
            case 1:
                filename = argv[optional_index];
                break;
            case 2:
                x_index = extract_index(argv[optional_index]);
                break;
            case 3:
                y_index = extract_index(argv[optional_index]);
                break;
            case 4:
                z_index = extract_index(argv[optional_index]);
                break;
            case 5:
                testtype = (uint) *argv[optional_index];
                break;
            case 6:
                num_permutation = (uint) *argv[optional_index];
                break;
            case 7:
                distance = (bool) *argv[optional_index];
            default:
                break;
        }

        optional_index++;
    }
    return 0;
}

std::vector<uint> ArgumentHandler::extract_index(std::string arg) {
    std::vector<uint> index;

    std::stringstream ss(arg);
    uint i;
    while (ss >> i) {
        // Start from 0 as traditional C++ format
        index.push_back(i - 1);

        if (ss.peek() == ',')
            ss.ignore();
    }
    return index;
}
