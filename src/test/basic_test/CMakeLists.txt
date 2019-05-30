include_directories(../lib/gtest-1.7.0/include ../lib/gtest-1.7.0 ../../../src)

add_executable(runBasicCDCTest
        basic_check.cpp cdc_check.cpp ../../utility.h
        ../../utility.cpp ../../kde.cpp
        ../../cdc_stats.cpp ../../cdc_statstics_method.cpp test_setting.h)

target_link_libraries(runBasicCDCTest gtest gtest_main)
