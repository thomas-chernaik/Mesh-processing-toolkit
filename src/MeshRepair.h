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
    // method to walk around an unpaired set of edges and return the boundary as a set of edges
    std::vector<Edge> walkAroundEdge(int edgeIndex);

    // method to get all unpaired edges starting at a vertex
    std::vector<Edge> getUnpairedEdges(int vertexIndex);

    // method to choose the smallest angled edge from a set of edges
    Edge chooseSmallestAngledEdge(std::vector<Edge> edges, Edge currentEdge);

    // method to fill a hole in the mesh by adding a point at the center of gravity of the boundary and connecting it to all the vertices in the boundary
    void fillHole(std::vector<Edge> boundary);

    void removeAllButLargestComponent();

    // method to get the centre of gravity of a set of vertices
    Vertex getCentreOfGravity(const std::vector<int>& vertices);

    // method to remove any items in vector2 from vector1
    void removeItems(std::vector<Edge>& vector1, const std::vector<Edge>& vector2);


protected:
// method to get the angle between two edges
float getAngleBetweenEdges(Edge edge1, Edge edge2);
};


#endif //MODELLINGCWK1_MESHREPAIR_H
