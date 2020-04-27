add_library(lib_commons
        common.h
        filehandler.h
        filehandler.cpp
        tsc_x86.h
        performance.h)
add_library(lib_naive naive.cpp)
add_library(lib_rotate rotate.h rotate.cpp)

add_executable(fic exec/fic.cpp)
add_executable(test_rotate exec/test_rotate.cpp)

target_link_libraries(lib_naive PRIVATE lib_rotate)
target_link_libraries(fic PUBLIC lib_commons lib_naive)
target_link_libraries(test_rotate PUBLIC lib_rotate)
