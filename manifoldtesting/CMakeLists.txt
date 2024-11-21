cmake_minimum_required(VERSION 3.10)
project(manifoldtesting)

# Specify C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Add the executable
add_executable(manifoldtesting main.cpp ../src/ManifoldTester.cpp ../src/DirectedEdge.cpp ../src/FaceIndex.cpp ../triangle_renderer/Cartesian3.cpp
)


# link the utilities library to the executable
target_link_libraries(manifoldtesting PUBLIC utilities)