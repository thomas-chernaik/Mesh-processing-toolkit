The program will write the simplified file to the same location, as a .tri file with the same name as the input file with "_simplified_<iterations>_iterations" appended to the filename.

If the mesh can't be simplified as many times as requested, then the program will output a file simplified as much as possible.

If you ask the mesh to simplify so much it won't have enough vertices left to form a mesh (4), the program will output an error.

See the readme in the root directory for information on compiling the program.

Run with:
./<path to simplification> [<path to file> [number of iterations]]

If no file is given, the program will prompt the user for a file.
If no number of iterations is given, the program will prompt the user for a number of iterations.
