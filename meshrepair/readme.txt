The program will discard non-connected components except for the largest, where a connected comoponent is a set of faces that are connected by shared edges. This is what blender "clean up -> delete loose" does.

The program will write the repaired file to the same location, with the same extension, but add _repaired to the filename. If the mesh can't be repaired then no file will be written, and an error should be printed to the console.

see the readme in the root directory for information on compiling the program.

run with:
./<path to meshrepair> [path to file]
if no file is given, the program will prompt the user for a file.