cmake_minimum_required(VERSION 3.10)

project(FIC VERSION 0.1 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}) # place binaries in build root

#
# default flags (they get extended by options)
#
set(DEFAULT_FLAGS "-Wall -lboost_program_options")
set(RELEASE_FLAGS "-O3 -DNDEBUG")
set(DEBUG_FLAGS "-O0")

add_subdirectory(src)
