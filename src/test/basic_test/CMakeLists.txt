include_directories(../lib/gtest-1.7.0/include ../lib/gtest-1.7.0 ../../../src)

add_executable(runBasicCDCTest
        basic_check.cpp ../../utility.cpp ../../kde.cpp)

target_link_libraries(runBasicCDCTest gtest gtest_main)
