## Mesh Conversion Utilities:
- face2faceindex: Converts raw triangle soup files (.tri) to a face-indexed format (.face).
- faceindex2directededge: Constructs a directed-edge (half-edge) structure from face-indexed files, outputting .diredge files for efficient traversal.
-  Manifold Validation: Detects non-manifold edges and vertices, flagging issues in a results file.
- Topological Analysis (part of the manifold testing): Computes surface genus using Euler’s formula for valid manifold meshes.
- Automated Repair: Closes holes in meshes by walking boundary edges and averaging vertex positions.
- Curvature-Driven Simplification: Reduces mesh complexity using Gaussian curvature metrics while preserving topological integrity.

## Build instructions
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

see the readme in the directories of the individual programs for instructions on running them, and more information on what they do.
