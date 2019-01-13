//
// Created by JinZhu on 2019/1/12.
//

#ifndef CDC_ARGUMENTSCHECKER_H
#define CDC_ARGUMENTSCHECKER_H

#include <string>
#include <vector>
#include "global.h"
#include <map>

class ArgumentHandler {
public:
    void checkArguments();
    ArgumentHandler(int argc, char **argv);
    uint processArguments();
    std::string filename;
    std::vector<uint> x_index;
    std::vector<uint> y_index;
    std::vector<uint> z_index;
    uint testtype;
    uint num_permutation;
    bool distance;

private:
    void displayHelp();
    char** argv;
    int argc;
    std::vector<uint > extract_index(std::string arg);

};


#endif //CDC_ARGUMENTSCHECKER_H
