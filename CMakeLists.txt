cmake_minimum_required(VERSION 3.16)
project(Opt)


find_package(Eigen3 REQUIRED)
include_directories(include)
include_directories(${EIGEN3_INCLUDE_DIRS})

set(source wolfe_simplex/wolfe_simplex.cpp
        active_set_dual/QuadProg.cpp
        interior_point_qp/QuadraticProgramming.cpp)
add_executable(main main.cpp ${source})
