cmake_minimum_required(VERSION 3.16)
project(Opt)

include_directories(include)
set(source wolfe_simplex.cpp)
add_executable(main main.cpp ${source})
