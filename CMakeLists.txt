cmake_minimum_required(VERSION 3.10)

# Set the project name
project(modellingcwk1)

# if no build type is specified, default to Release
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()
# set up the release flags
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(CMAKE_C_FLAGS_RELEASE "-O3")


# Add subdirectories for each executable
add_subdirectory(utilities)
add_subdirectory(face2faceindex)
add_subdirectory(faceindex2directededge)
add_subdirectory(manifoldtesting)
add_subdirectory(meshrepair)
add_subdirectory(simplification)