cmake_minimum_required(VERSION 3.10)

project (Celestial)

set (CMAKE_CXX_STANDARD 14)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=c++17 -lm")

include_directories(src/include)
set (NBODY_SRCS
        src/main.cpp
        src/integration.cpp
        src/vector.cpp
        src/menu.cpp
        src/parser.cpp src/include/matrix.h src/matrix.cpp)

add_executable(Celestial ${NBODY_SRCS})
