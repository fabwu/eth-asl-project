FIND_PACKAGE(Boost COMPONENTS program_options REQUIRED)
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})

add_library(lib_commons
        lib/common.hpp
        lib/filehandler.hpp
        lib/filehandler.cpp
        lib/tsc_x86.h
        lib/performance.h
        lib/performance.c
        )
add_library(lib_types lib/types.h lib/types.c)
add_library(lib_rotate lib/rotate.h lib/rotate.c)
add_library(lib_queue lib/queue.h lib/queue.c)
add_library(lib_naive naive.c) # baseline
add_library(lib_fast fast.c) # working file
add_library(lib_01_precompute_indices 01_precompute_indices.c)
add_library(lib_02_precomputations 02_precomputations.c)

add_executable(test_rotate exec/test_rotate.cpp)
add_executable(baseline exec/fic.cpp)
add_executable(fic exec/fic.cpp)
add_executable(01_precompute_indices exec/fic.cpp)
add_executable(02_precomputations exec/fic.cpp)

target_link_libraries(lib_commons PRIVATE lib_types lib_queue)
target_link_libraries(lib_rotate PRIVATE lib_types)
target_link_libraries(lib_naive PRIVATE lib_types lib_commons lib_queue lib_rotate)
target_link_libraries(lib_fast PRIVATE lib_types lib_commons lib_queue lib_rotate)
target_link_libraries(lib_01_precompute_indices PRIVATE lib_types lib_commons lib_queue lib_rotate)
target_link_libraries(lib_02_precomputations PRIVATE lib_types lib_commons lib_queue lib_rotate)

target_link_libraries(test_rotate PUBLIC lib_rotate)
target_link_libraries(baseline PUBLIC lib_types lib_queue lib_commons lib_rotate lib_naive ${Boost_LIBRARIES})
target_link_libraries(fic PUBLIC lib_types lib_queue lib_commons lib_rotate lib_fast ${Boost_LIBRARIES})
target_link_libraries(01_precompute_indices PUBLIC lib_types lib_queue lib_commons lib_rotate lib_01_precompute_indices ${Boost_LIBRARIES})
target_link_libraries(02_precomputations PUBLIC lib_types lib_queue lib_commons lib_rotate lib_02_precomputations ${Boost_LIBRARIES})
