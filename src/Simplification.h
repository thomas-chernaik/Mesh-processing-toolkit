//
// Created by thomas on 18/11/24.
//

#ifndef MODELLINGCWK1_SIMPLIFICATION_H
#define MODELLINGCWK1_SIMPLIFICATION_H

#include "MeshRepair.h"

class Simplification : public MeshRepair
{
public:
    // method to simplify the mesh
    void simplifyMesh();
private:
    // method to find the vertex with the nth smallest curvature
    int findSmallestCurvature(int n);

    // method to find the mean/gaussian curvature of a vertex
    float findCurvature(int vertexIndex);

    // method to remove a vertex from the mesh and triangulate the hole
    // returns false if it does not maintain the eulerian condition
    bool removeVertex(int vertexIndex);




    bool isEulerian();

    // stuff we need to be able to backtrack if we lose the eulerian condition
    // the number of faces we just added
    int facesAdded;
    // the edges of the hole we just filled
    std::vector<Edge> holeEdges;
    // the removed vertex
    Vertex removedVertex;
    void backtrack();

};


#endif //MODELLINGCWK1_SIMPLIFICATION_H
