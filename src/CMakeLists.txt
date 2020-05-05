FIND_PACKAGE(Boost COMPONENTS program_options REQUIRED)
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})

add_library(lib_commons
        lib/types.h
        lib/common.h
        lib/filehandler.h
        lib/filehandler.cpp
        lib/tsc_x86.h
        lib/performance.h
        lib/queue.h
        lib/queue.c
)
add_library(lib_naive naive.c)
add_library(lib_rotate lib/rotate.h lib/rotate.c)

add_executable(fic exec/fic.cpp)
add_executable(test_rotate exec/test_rotate.cpp)

target_link_libraries(lib_naive PRIVATE lib_rotate)
target_link_libraries(fic PUBLIC lib_commons lib_naive ${Boost_LIBRARIES})
target_link_libraries(test_rotate PUBLIC lib_rotate)
