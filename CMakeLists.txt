cmake_minimum_required(VERSION 3.12)
project(ArchiDAISIE)

set(CMAKE_CXX_STANDARD 14)

if(MSVC)
    add_compile_options(/W4)
else()
    add_compile_options(-Wall -Wextra -Wpedantic)
endif()

include_directories(.)

add_executable(ArchiDAISIE
        ArchiMainSource.cpp
        DrawEvents.cpp
        event_type.cpp
        SpeciesID.cpp
        Species.cpp
        Island.cpp
        Archipelago.cpp
        tests.cpp)
