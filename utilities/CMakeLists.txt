# utilities/CMakeLists.txt
cmake_minimum_required(VERSION 3.10)
project(utilities)

# Specify C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Add the library
add_library(utilities vertex.cpp face.cpp)

# Specify include directories for the library
target_include_directories(utilities PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})