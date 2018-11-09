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
        event_type.cpp
        SpeciesID.cpp
        DrawEvents.cpp
        Species.cpp
        Island.cpp
        Archipelago.cpp
        STT.cpp
        STTtable.cpp
        free_functions.cpp
        ArchiDAISIE.cpp
        tests.cpp
        ArchiMainSource.cpp)
