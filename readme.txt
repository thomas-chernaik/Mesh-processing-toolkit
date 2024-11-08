This project is built using CMake.

All aspects of the project will be built in one go.

Build instructions:
1. Create a build directory where you want the build to go. I recommend a build directory in the root of the project.
    e.g. `mkdir build && cd build` from the root of the project.
2. Run cmake. This should be run from inside the build directory, and given the argument of the path to the root of the project.
    e.g. `cmake ..`
3. Run make. This should be run from inside the build directory.
    e.g. `make`

The project will then be built.

Running instructions:

The face to face index converter program (Question 1(a) can be run as follows:
./<path_to_build>/face2faceindex/face2faceindex [path_to_.tri_file]
if no path is given, the program will prompt for the path.
A .face file will be created at the same path as the .tri file, with the same name.

The face index to directed edge converter program (Question 1(b)) can be run as follows:
./<path_to_build>/faceindex2directededge/faceindex2directededge [path_to_.face_file]
if no path is given, the program will prompt for the path.
A .diredge file will be created at the same path as the .face file, with the same name.
