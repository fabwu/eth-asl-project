cmake_minimum_required(VERSION 3.10)

project(FIC VERSION 0.1 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 14)

add_executable(main main.cpp)
target_link_libraries(main)