cmake_minimum_required(VERSION 3.16)
project(Opt)


SET (CMAKE_CXX_STANDARD 17)

find_package(Eigen3 REQUIRED)
include_directories(include)
include_directories(${EIGEN3_INCLUDE_DIRS})

#set(source qp/wolfe_simplex/wolfe_simplex.cpp
#        qp/active_set_dual/QuadProg.cpp
#        opt_global/SA/SA.cpp
#        qp/interior_point_qp/QuadraticProgramming.cpp
#        opt_global/GA/GA.cpp
#        opt_global/GA/GA.h
#        qp/test/test_qp.cpp
#        opt_unconstr/test/test_unconstr.cpp
#        line_search/test/test_line_search.cpp
#        opt_constr/test/test_constr.cpp)
#set(source
#        opt_global/GA/GA.cpp
#        opt_global/GA/GA.h
#        opt_global/test/test_global.cpp)

#set(source
#        opt_global/BB/bb.cpp
#        opt_global/BB/bb.h)
set(source
        lp/simplex.cpp
        lp/simplex.h)
#add_executable(main main.cpp ${source})
#add_executable(testGA ${source})
add_executable(simplex ${source})