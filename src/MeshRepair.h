//
// Created by thomas on 15/11/24.
//

#ifndef MODELLINGCWK1_MESHREPAIR_H
#define MODELLINGCWK1_MESHREPAIR_H

#include "ManifoldTester.h"

class MeshRepair : public ManifoldTester
{
public:
    // set runwitherrors to true
    MeshRepair() : ManifoldTester(true)
    {};

    // method to repair the mesh
    void repairMesh();

    // method to write the repaired mesh to a .tri file
    void writeRepairedMeshTri(const std::string& filename);

    // method to write the repaired mesh to a .face file
    void writeRepairedMeshFace(const std::string& filename);

    // it will write to a .diredge file in the default inherited writeFile method


private:
    // method to return the first unpaired edge
    int findLooseEdge();
    // method to walk around an unpaired set of edges and return the boundary
    std::vector<edge> walkAroundEdge(int edgeIndex);
    // method to fill a hole in the mesh by adding a point at the center of gravity of the boundary and connecting it to all the vertices in the boundary
    void fillHole(std::vector<edge> boundary);

    void removeAllButLargestComponent();



};


#endif //MODELLINGCWK1_MESHREPAIR_H
