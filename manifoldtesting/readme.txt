The triangle intersection does not return an error and end the program, but instead returns a list of the triangles it thinks intersects. Whether the triangles actually intersect is up to the user to determine.
If the mesh contains multiple components, the program will output the number of components and the genus of each, assuming each individual component is a manifold.

Otherwise, the program will output the first encountered edge or vertex that causes the mesh to be non-manifold.

see the readme in the root directory for information on compiling the program.

run with:
./<path to manifoldtesting> [path to file]
if no file is given, the program will prompt the user for a file.