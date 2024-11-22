cmake_minimum_required(VERSION 3.10)
project(simplification)

# Specify C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Add the executable
add_executable(simplification main.cpp ../src/Simplification.cpp ../src/FaceIndex.cpp ../src/DirectedEdge.cpp ../src/ManifoldTester.cpp ../triangle_renderer/Cartesian3.cpp ../src/MeshRepair.cpp)


# link the utilities library to the executable
target_link_libraries(simplification PUBLIC utilities)