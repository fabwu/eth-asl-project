cmake_minimum_required(VERSION 3.10)

project(FIC VERSION 0.1 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}) # place binaries in build root

#
# default flags (they get extended by options)
#
set(DEFAULT_FLAGS "-Wall -lboost_program_options -mavx2 -mfma")
set(RELEASE_FLAGS "-DNDEBUG")
set(DEBUG_FLAGS "-g -O0")

option(OPT_ARCH "optimize for arch" ON)
option(OPT_O1 "O1" OFF)
option(OPT_O2 "O2" OFF)
option(OPT_O3 "O3" ON)
option(OPT_OFAST "OFast" OFF)
option(OPT_VEC "vecorize" ON)
option(OPT_UNSAFE_MATH "unsafe math" OFF)
option(OPT_UNROLL_LOOPS "unroll loops" ON)
option(OPT_FAST_MATH "fast math" OFF)

if (OPT_ARCH)
    set(RELEASE_FLAGS "${RELEASE_FLAGS} -march=native")
endif ()

if (OPT_O1)
    set(RELEASE_FLAGS "${RELEASE_FLAGS} -O1")
endif()

if (OPT_O2)
    set(RELEASE_FLAGS "${RELEASE_FLAGS} -O2")
endif()

if (OPT_O3)
    set(RELEASE_FLAGS "${RELEASE_FLAGS} -O3")
endif()

if (OPT_OFAST)
    set(RELEASE_FLAGS "${RELEASE_FLAGS} -Ofast")
endif()

if (OPT_VEC)
    set(RELEASE_FLAGS "${RELEASE_FLAGS} -ftree-vectorize")
endif()

if (OPT_UNSAFE_MATH)
    set(RELEASE_FLAGS "${RELEASE_FLAGS} -funsafe-math-optimizations")
endif()

if (OPT_UNROLL_LOOPS)
    set(RELEASE_FLAGS "${RELEASE_FLAGS} -funroll-all-loops")
endif()

if (OPT_FAST_MATH)
    set(RELEASE_FLAGS "${RELEASE_FLAGS} -ffast-math")
endif()


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
