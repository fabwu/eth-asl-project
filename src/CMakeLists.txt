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
add_library(lib_rotate rotate.h rotate.cpp)

add_executable(compress_decompress exec/compress_decompress.cpp)
add_executable(benchmark exec/benchmark.cpp)
add_executable(test_rotate exec/test_rotate.cpp)

target_link_libraries(compress_decompress PUBLIC lib_commons lib_naive)
target_link_libraries(benchmark PUBLIC lib_commons lib_naive)
target_link_libraries(test_rotate PUBLIC lib_rotate)
