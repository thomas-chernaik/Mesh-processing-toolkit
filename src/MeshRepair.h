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
    int findLooseEdge();
    std::vector<edge> walkAroundEdge(int edgeIndex);
    void fillHole(std::vector<edge> boundary);



};


#endif //MODELLINGCWK1_MESHREPAIR_H
