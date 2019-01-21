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
    CONDITION_DISTANCE_CORRELATION = 2,
    CONDITION_BALL_COVARIANCE = 3,
    CONDITION_BALL_CORRELATION = 4
};

enum StatsMethod {
    HYPOTHESIS_TESTS = 1,
    SURE_INDEPENDENCE_SCREENING = 2,
    STATISTICS_VALUE = 3
};

const uint DEFAULT_LOCAL_BOOTSTRAP_NUM = 99;

const uint DEFAULT_NUM_THREADS = 1;

const double DEFAULT_DISTANCE_INDEX = 1.0;

#endif //GLOBAL_H
