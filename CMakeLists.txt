cmake_minimum_required(VERSION 3.10)

project(FIC VERSION 0.1 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}) # place binaries in build root

#
# default flags (they get extended by options)
#
set(DEFAULT_FLAGS "-Wall -lboost_program_options -mavx2 -mfma")
set(RELEASE_FLAGS "-DNDEBUG -march=native")
set(DEBUG_FLAGS "-g -O0")

# see here: https://gcc.gnu.org/onlinedocs/gcc/Optimize-Options.html

option(OPT_ARCH "optimize for arch" ON)
option(OPT_O1 "O1" OFF)
option(OPT_O2 "O2" OFF)
option(OPT_O3 "O3" ON)
option(OPT_OFAST "OFast" OFF)
option(OPT_FMA "fma" OFF)
option(OPT_UNROLL_LOOPS "unroll loops" OFF)


if(ICC)
  set(CMAKE_CXX_COMPILER "/opt/intel/bin/icpc")
  set(CMAKE_C_COMPILER "/opt/intel/bin/icc")
  set(CMAKE_C_FLAGS "${INTEL_DEFAULT_FLAGS} -std=c99")
  set(CMAKE_CXX_FLAGS "${INTEL_DEFAULT_FLAGS} -std=c++11")
  set(CMAKE_C_FLAGS_RELEASE   ${INTEL_OPTIMIZATION_FLAGS})
  set(CMAKE_CXX_FLAGS_RELEASE ${INTEL_OPTIMIZATION_FLAGS})
endif(ICC)


if (OPT_O1)
  set(RELEASE_FLAGS "${RELEASE_FLAGS} -O1")
endif(OPT_O1)

if (OPT_O2)
  set(RELEASE_FLAGS "${RELEASE_FLAGS} -O2")
endif(OPT_O2)

if (OPT_O3)
  set(RELEASE_FLAGS "${RELEASE_FLAGS} -O3")
endif(OPT_O3)

if (OPT_OFAST)
  set(RELEASE_FLAGS "${RELEASE_FLAGS} -Ofast")
endif(OPT_OFAST)

if (OPT_FMA)
  if(ICC)
    set(RELEASE_FLAGS "${RELEASE_FLAGS} -fma")
  else(ICC)
    set(RELEASE_FLAGS "${RELEASE_FLAGS} -mfma")
  endif(ICC)
else(OPT_FMA)
  if(ICC)
    set(RELEASE_FLAGS "${RELEASE_FLAGS} -no-fma")
  else(ICC)
    set(RELEASE_FLAGS "${RELEASE_FLAGS} -mno-fma")
  endif(ICC)
endif(OPT_FMA)

if (OPT_UNROLL_LOOPS)
  set(RELEASE_FLAGS "${RELEASE_FLAGS} -funroll-loops")
else(OPT_UNROLL_LOOPS)
  set(RELEASE_FLAGS "${RELEASE_FLAGS} -fno-unroll-loops")
endif(OPT_UNROLL_LOOPS)


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
