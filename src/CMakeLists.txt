cmake_minimum_required(VERSION 3.12)
project(cdc)

set(CMAKE_CXX_STANDARD 11)
set(SOURCE_FILES main.cpp)

add_executable(cdc ${SOURCE_FILES} kde.cpp kde.h global.h utility.cpp utility.h)

add_subdirectory(test)
include_directories(test/lib/gtest/include test/lib/gtest)
