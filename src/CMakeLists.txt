cmake_minimum_required(VERSION 3.12)
project(cdc)

set(CMAKE_CXX_STANDARD 11)
set(SOURCE_FILES main.cpp cdc_stats.cpp cdc_stats.h utility.cpp utility.h cdc_statstics_method.cpp cdc_statstics_method.h data.cpp data.h ArgumentHandler.cpp ArgumentHandler.h Bandwidth.cpp Bandwidth.h)

add_executable(cdc ${SOURCE_FILES} kde.cpp kde.h global.h)

add_subdirectory(test)
include_directories(test/lib/gtest/include test/lib/gtest)
