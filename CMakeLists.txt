cmake_minimum_required(VERSION 3.10)
project(Celestial)

set(CMAKE_CXX_STANDARD 14)

include_directories(src/include)

add_executable(Celestial
        src/include/benchmark.h
        src/include/exit_codes.h
        src/include/integration.h
        src/include/menu.h
        src/include/parser.h
        src/include/planet_data.h
        src/include/structures.h
        src/include/astro_times.h
        src/main.cpp
        src/parser.cpp)
