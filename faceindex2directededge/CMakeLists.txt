cmake_minimum_required(VERSION 3.10)
project(face2faceindex)

# Specify C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)


# Add the executable
add_executable(faceindex2directededge main.cpp faceindex2directededge.cpp)



# link the utilities library to the executable
target_link_libraries(faceindex2directededge PUBLIC utilities)