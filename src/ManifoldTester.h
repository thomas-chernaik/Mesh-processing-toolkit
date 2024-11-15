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
    // constructor with readwitherrors as false
    explicit ManifoldTester(bool readWithErrors = false) : DirectedEdge(readWithErrors)
    {}

    // method to test if the mesh is manifold
    void testManifold();

    // method to calculate the genus of the mesh using the euler formula
    std::vector<int> CalculateGenus();

    // a method to read in a .diredge file
    void readFileDiredge(const std::string &filename);

protected:
    // method to test if the mesh contains multiple components
    void testMultipleComponents();

    // method to get the one ring (as a list of edge indices) of a vertex
    std::vector<Edge> getOneRing(int vertexIndex);

    // method to get the one ring (as a list of vertex indices) of a vertex
    std::vector<int> getOneRingVertices(int vertexIndex);

private:
    // method to test if the mesh has any pinch points
    void testPinchPoints();

    // method to test if the mesh has any self intersections
    void testSelfIntersections();

    // method to test the intersection of two triangles
    bool testTriangleIntersection(int f1, int f2);


    // method to check if a list of edges is a single cycle
    static bool isSingleCycle(std::vector<Edge> edges);

    // returns the dot product of two Cartesian3 vectors
    static float dotProduct(const Cartesian3 &a, const Cartesian3 &b)
    {
        return a.x * b.x + a.y * b.y + a.z * b.z;
    }

    // returns the dot product of a Cartesian3 vector and a Vertex data type
    static float dotProduct(const Cartesian3 &a, const Vertex &b)
    {
        return a.x * b.x + a.y * b.y + a.z * b.z;
    }

    // converts a Vertex to a Cartesian3
    static Cartesian3 VertexToCartesian3(const Vertex &v)
    {
        return {v.x, v.y, v.z};
    }

    // returns if two line segments intersect
    static inline bool EdgesIntersect(Cartesian3 v1, Cartesian3 v2, Cartesian3 v3, Cartesian3 v4);

    // returns if a point is inside a triangle
    static inline bool
    TriangleContainsVertex(Cartesian3 &p, Cartesian3 &q, Cartesian3 &r, Cartesian3 &point);

    // contains the face indexes of each component
    std::vector<std::vector<int>> components;

};


#endif //MODELLINGCWK1_MANIFOLDTESTER_H
