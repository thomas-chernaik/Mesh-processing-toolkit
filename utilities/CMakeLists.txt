# utilities/CMakeLists.txt
cmake_minimum_required(VERSION 3.10)
project(utilities)


# Add the library
add_library(utilities vertex.cpp face.cpp
        edge.h)

# Specify include directories for the library
target_include_directories(utilities PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})