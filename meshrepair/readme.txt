The program deviates from the specification as follows:
When a hole is being repaired and contains a reflex angle on the boundary, the program will not repair the hole in the standard way of adding a new vertex, as this could create a self-intersection.
Instead, the program will triangulate the hole in a way that will not create a self-intersection.

The program will discard non-connected components except for the largest, where a connected comoponent is a set of faces that are connected by shared edges. This is what blender "clean up -> delete loose" does.
