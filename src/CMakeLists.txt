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
# add_library(lib_03_jonas 03_jonas.c)
# add_library(lib_04_jonas 04_jonas_one_rtd_function_per_degree_no_optimizations.c)
add_library(lib_05_jonas 05_rtd_in_main_loop.c)
# add_library(lib_06_jonas 06_optimize_compute_and_brightness.c)
# add_library(lib_07_jopa  07_jopa.c)
add_library(lib_08_jopa  08_jopa.c)

##################################################################################

add_executable(test_rotate exec/test_rotate.cpp)
add_executable(bench_rtd exec/bench_rtd.c)
add_executable(baseline exec/fic.cpp)
add_executable(fic exec/fic.cpp)
add_executable(01_precompute_indices exec/fic.cpp)
# 2 has to be fixed for large images (see 5)
add_executable(02_precomputations exec/fic.cpp)
# delete 3
# delete 4
add_executable(05_rtd_no_precompute exec/fic.cpp)
# delete 6
# delete 7
add_executable(08_remove_function_and_add_fma exec/fic.cpp)


##################################################################################

target_link_libraries(lib_commons PRIVATE lib_types lib_queue)
target_link_libraries(lib_rotate PRIVATE lib_types)
target_link_libraries(lib_naive PRIVATE lib_types lib_commons lib_queue lib_rotate)
target_link_libraries(lib_fast PRIVATE lib_types lib_commons lib_queue lib_rotate)
target_link_libraries(lib_01_precompute_indices PRIVATE lib_types lib_commons lib_queue lib_rotate)
target_link_libraries(lib_02_precomputations PRIVATE lib_types lib_commons lib_queue lib_rotate)
# target_link_libraries(lib_03_jonas PRIVATE lib_types lib_commons lib_queue lib_rotate)
# target_link_libraries(lib_04_jonas PRIVATE lib_types lib_commons lib_queue lib_rotate)
target_link_libraries(lib_05_jonas PRIVATE lib_types lib_commons lib_queue lib_rotate)
# target_link_libraries(lib_06_jonas PRIVATE lib_types lib_commons lib_queue lib_rotate)
# target_link_libraries(lib_07_jopa PRIVATE lib_types lib_commons lib_queue lib_rotate)
target_link_libraries(lib_08_jopa PRIVATE lib_types lib_commons lib_queue lib_rotate)

##################################################################################

target_link_libraries(test_rotate PUBLIC lib_rotate)
target_link_libraries(bench_rtd PUBLIC -lm lib_types lib_queue)
target_link_libraries(baseline PUBLIC lib_types lib_queue lib_commons lib_rotate lib_naive ${Boost_LIBRARIES})
target_link_libraries(fic PUBLIC lib_types lib_queue lib_commons lib_rotate lib_fast ${Boost_LIBRARIES})
target_link_libraries(01_precompute_indices PUBLIC lib_types lib_queue lib_commons lib_rotate lib_01_precompute_indices ${Boost_LIBRARIES})
target_link_libraries(02_precomputations PUBLIC lib_types lib_queue lib_commons lib_rotate lib_02_precomputations ${Boost_LIBRARIES})
# target_link_libraries(03_jonas PUBLIC lib_types lib_queue lib_commons lib_rotate lib_03_jonas ${Boost_LIBRARIES})
# target_link_libraries(04_jonas PUBLIC lib_types lib_queue lib_commons lib_rotate lib_04_jonas ${Boost_LIBRARIES})
target_link_libraries(05_rtd_no_precompute PUBLIC lib_types lib_queue lib_commons lib_rotate lib_05_jonas ${Boost_LIBRARIES})
# target_link_libraries(06_jonas PUBLIC lib_types lib_queue lib_commons lib_rotate lib_06_jonas ${Boost_LIBRARIES})
# target_link_libraries(07_jopa PUBLIC lib_types lib_queue lib_commons lib_rotate lib_07_jopa ${Boost_LIBRARIES})
target_link_libraries(08_remove_function_and_add_fma PUBLIC lib_types lib_queue lib_commons lib_rotate lib_08_jopa ${Boost_LIBRARIES})
