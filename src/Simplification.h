//
// Created by thomas on 18/11/24.
//

#ifndef MODELLINGCWK1_SIMPLIFICATION_H
#define MODELLINGCWK1_SIMPLIFICATION_H

#include "MeshRepair.h"
#include "../triangle_renderer/Cartesian3.h"

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

    // method to find the mean curvature of a vertex
    float findMeanCurvature(int vertexIndex);

    // method to find the gaussian curvature of a vertex
    float findGaussianCurvature(int vertexIndex);

    // method to compute the discrete laplace beltrami operator of a vertex
    Cartesian3 computeLaplaceBeltrami(int vertexIndex);

    // method to get the vertices in the one ring of a vertex
    std::vector<int> getOneRingVertices(int vertexIndex);

    // method to get the angle between two vectors
    float getAngleBetweenVectors(Cartesian3 vector1, Cartesian3 vector2);

    // method to remove a vertex from the mesh and triangulate the hole
    // returns false if it does not maintain the eulerian condition
    bool removeVertex(int vertexIndex);

    // method to find the smallest angle between two adjacent edges on a boundary
    int findSmallestAngle(const std::vector<Edge> &boundary);

    // method to triangulate a hole
    void triangulateHole(std::vector<Edge> &boundary);




    // method to test if the eulerian condition is maintained
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
