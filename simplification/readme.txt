The simplification code is not guaranteed to output a manifold mesh. The triangulation can result in non-manifold edges. A cleanup function at the end of the simplification process is run to try and eliminate these edges, but it is not guaranteed to catch all cases where a large number of iterations have been run.

The program will write the simplified file to the same location, as a .tri file with the same name as the input file with "_simplified" appended to the filename.

If the mesh can't be simplified as many times as requested, then the program will output a file simplified as much as possible.

See the readme in the root directory for information on compiling the program.

Run with:
./<path to simplification> [<path to file> [number of iterations]]

If no file is given, the program will prompt the user for a file.
If no number of iterations is given, the program will prompt the user for a number of iterations.
