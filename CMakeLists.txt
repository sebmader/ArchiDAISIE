cmake_minimum_required(VERSION 3.12)
project(ArchiDAISIE)

set(CMAKE_CXX_STANDARD 14)

include_directories(.)

add_executable(ArchiDAISIE
        ArchiMainSource.cpp
        Island.cpp
        Species.cpp)
