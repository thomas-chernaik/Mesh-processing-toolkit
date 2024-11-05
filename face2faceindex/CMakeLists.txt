cmake_minimum_required(VERSION 3.10)
project(face2faceindex)

# Specify C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Add the executable
add_executable(face2faceindex main.cpp face2faceindex.cpp)