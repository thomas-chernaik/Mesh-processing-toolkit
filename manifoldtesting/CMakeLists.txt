cmake_minimum_required(VERSION 3.10)
project(manifoldtesting)


# Add the executable
add_executable(manifoldtesting main.cpp ../src/ManifoldTester.cpp ../src/DirectedEdge.cpp ../src/FaceIndex.cpp ../triangle_renderer/Cartesian3.cpp
)


# link the utilities library to the executable
target_link_libraries(manifoldtesting PUBLIC utilities)