The face to face index converter program can be built with the following steps:

step 1:
    Create the build directory and enter it
        e.g. mkdir build && cd build
step 2:
    run cmake
        cmake ..
step 3:
    run make
        make
step 4:
    run the program
    option 1:
        ./face2faceindex <path_to_.tri_file>
    option2:
        ./face2faceindex
        the file will then be prompted

the program will write to the same path but replace the .tri with .face