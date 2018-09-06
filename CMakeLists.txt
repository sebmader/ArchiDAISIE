cmake_minimum_required(VERSION 3.12)
project(ArchiDAISIE)

set(CMAKE_CXX_STANDARD 14)

include_directories(.)

add_executable(ArchiDAISIE
        ArchiMainSource.cpp
        DrawEvents.cpp
        Species.cpp
        Island.cpp
        Archipelago.cpp)
