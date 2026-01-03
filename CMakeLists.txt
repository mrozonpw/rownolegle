cmake_minimum_required(VERSION 3.20)
project(simulated_annealing)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "-O2 -pthread")

add_executable(simulated_annealing
        main.cpp
        algorithm.cpp
        functions.cpp
        functions.cpp
        algorithm.h
        algorithm.cpp)
