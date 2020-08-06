cmake_minimum_required(VERSION 3.10)

project (Celestial)

set (CMAKE_CXX_STANDARD 14)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Werror -std=c++17 -lm") 

include_directories(src/include)
set (source_dir "${PROJECT_SOURCE_DIR}/src/")

file (GLOB source_files "${source_dir}/main.cpp")

add_executable(Celestial ${source_files} src/include/astro_constants.h src/integration.cpp)
