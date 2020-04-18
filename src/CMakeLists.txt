cmake_minimum_required(VERSION 3.10)

project(FIC VERSION 0.1 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

add_library(lib_commons
        common.h
        filehandler.h
        filehandler.cpp
        tsc_x86.h)
add_library(lib_naive naive.cpp)

add_executable(exec exec/compress_decompress.cpp)
add_executable(benchmark exec/benchmark.cpp)
add_executable(test_rotate exec/rotate.cpp)

target_link_libraries(exec PUBLIC lib_commons lib_naive)
target_link_libraries(benchmark PUBLIC lib_commons lib_naive)
