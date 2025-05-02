cmake_minimum_required(VERSION 3.10)
project(bellman_ford_mpi C)

set(CMAKE_C_STANDARD 17)
set(CMAKE_C_STANDARD_REQUIRED ON)

if(NOT DEFINED ENV{MPICH_TARGET_DIR})
    message(FATAL_ERROR "MPICH_TARGET_DIR is not set. Please export it before running CMake.")
endif()

set(MPI_HOME $ENV{MPICH_TARGET_DIR})

add_executable(bellman_ford
    src/bellman-ford.c
    src/graph.c
)

target_include_directories(bellman_ford PRIVATE
    "${CMAKE_SOURCE_DIR}/include"
    "${MPI_HOME}/include"
)

target_link_directories(bellman_ford PRIVATE
    "${MPI_HOME}/lib"
)

target_link_libraries(bellman_ford PRIVATE mpi)