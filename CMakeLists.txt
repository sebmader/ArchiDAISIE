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
        free_functions.cpp
        DrawEvents.cpp
        event_type.cpp
        SpeciesID.cpp
        Species.cpp
        Island.cpp
        Archipelago.cpp
        ArchiDAISIE.cpp
        ArchiMainSource.cpp
        STT.cpp
        STTtable.cpp
        tests.cpp)

target_link_libraries(${PROJECT_NAME} stdc++fs)
set(CMAKE_CXX_FLAGS "-lstdc++fs")