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
set(DEBUG_FLAGS "-g -O0")

option(OPT_ARCH "optimize for arch" ON)

if (OPT_ARCH)
    set(RELEASE_FLAGS "${RELEASE_FLAGS} -march=native")
endif ()

#
# override cmake flags (gives us full control)
#
set(CMAKE_C_FLAGS ${DEFAULT_FLAGS})
set(CMAKE_CXX_FLAGS ${DEFAULT_FLAGS})

set(CMAKE_C_FLAGS_RELEASE ${RELEASE_FLAGS})
set(CMAKE_CXX_FLAGS_RELEASE ${RELEASE_FLAGS})

set(CMAKE_C_FLAGS_DEBUG ${DEBUG_FLAGS})
set(CMAKE_CXX_FLAGS_DEBUG ${DEBUG_FLAGS})

add_subdirectory(src)
