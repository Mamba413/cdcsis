//
// Created by JinZhu on 2018/12/26.
//

#ifndef GLOBAL_H
#define GLOBAL_H

#ifndef CDC_PI
#define CDC_PI 3.14159265358979323846
#endif

typedef unsigned int uint;

enum KernelType {
    KERNEL_GAUSSIAN = 1,
    KERNEL_RECTANGLE = 3
};

enum StatsType {
    CONDITION_DISTANCE_COVARIANCE = 1,
    CONDITION_DISTANCE_CORRELATION = 2
};

const bool DISTANCE = false;

#endif //GLOBAL_H
