//
// Created by thomas on 08/11/24.
//

#ifndef MODELLINGCWK1_MANIFOLDTESTER_H
#define MODELLINGCWK1_MANIFOLDTESTER_H

#include "DirectedEdge.h"
#include "edge.h"
#include "../triangle_renderer/Cartesian3.h"



// class manifold tester extends DirectedEdge, also containing methods to test for manifoldness
class ManifoldTester : public DirectedEdge
{
public:
    // method to test if the mesh is manifold
    void testManifold();

private:
    // method to test if the mesh has any pinch points
    void testPinchPoints();

    // method to test if the mesh has any self intersections
    void testSelfIntersections();

    // method to test the intersection of two triangles
    bool testTriangleIntersection(int f1, int f2);

    // method to test if the mesh contains multiple components
    void testMultipleComponents();

    // method to get the one ring (as a list of edge indices) of a vertex
    std::vector<Edge> getOneRing(int vertexIndex);

    // method to check if a list of edges is a single cycle
    bool isSingleCycle(std::vector<Edge> edges);

    float dotProduct(const Cartesian3 &a, const Cartesian3 &b)
    {
        return a.x * b.x + a.y * b.y + a.z * b.z;
    }
    float dotProduct(const Cartesian3 &a, const Vertex &b)
    {
        return a.x * b.x + a.y * b.y + a.z * b.z;
    }

    Cartesian3 VertexToCartesian3(const Vertex &v)
    {
        return {v.x, v.y, v.z};
    }

};


#endif //MODELLINGCWK1_MANIFOLDTESTER_H
