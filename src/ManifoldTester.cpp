//
// Created by thomas on 08/11/24.
//

#include <algorithm>
#include <chrono>
#include "ManifoldTester.h"

void ManifoldTester::testManifold()
{
    // A surface mesh is manifold if:
    // 1. Each edge is incident to at most two faces. This is already tested for when the file is read in
    // 2. There are no pinch points
    // 3. There are no self intersections
    // These tests are assuming that it is in fact a surface mesh.
    // There is nothing stopping the mesh from being multiple components, in which case it isn't a manifold mesh, but multiple

    // test for pinch points
    // start a timer
    auto start = std::chrono::high_resolution_clock::now();
    testPinchPoints();
    // stop the timer
    auto stop = std::chrono::high_resolution_clock::now();
    // output the time taken
    std::cout << "Pinch points took "
              << std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count()
              << " milliseconds" << std::endl;
    // test for self intersections
    start = std::chrono::high_resolution_clock::now();
    testSelfIntersections();
    stop = std::chrono::high_resolution_clock::now();
    std::cout << "Self intersections took "
              << std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count()
              << " milliseconds" << std::endl;
    // test for multiple components
    testMultipleComponents();

}

void ManifoldTester::testPinchPoints()
{
    // we can test for pinch points by checking that the one ring of every vertex contains a single cycle

    // for each vertex
    for (int v = 0; v < vertices.size(); v++)
    {
        // get the one ring of the vertex
        std::vector<Edge> oneRing = getOneRing(v);
        // check that the one ring is a single cycle
        if (!isSingleCycle(oneRing))
        {
            std::cerr << "Mesh is not manifold as vertex " << v << " has a pinch point" << std::endl;
            exit(-3);
        }
    }

}

void ManifoldTester::testSelfIntersections()
{
    // we can test for self intersections by checking that no two triangles intersect
    // TODO: optimise with an octree
    for (int f1 = 0; f1 < faces.size(); f1++)
    {
        for (int f2 = f1 + 1; f2 < faces.size(); f2++)
        {
            if (testTriangleIntersection(f1, f2))
            {
                // output the edges and vertices of the intersecting faces
                std::cerr << "Mesh is not manifold as faces " << f1 << " and " << f2 << " intersect" << std::endl;
                exit(-4);
            }
        }
    }
}

bool ManifoldTester::testTriangleIntersection(int f1, int f2)
{
    // calculate the plane for each face
    // an early out is after the first plane calculation if the second face is on the same side of the plane
    Cartesian3 face1Normal = VertexToCartesian3(vertices[faces[f1][1]] - vertices[faces[f1][0]]).cross(
            VertexToCartesian3(
                    vertices[faces[f1][2]] - vertices[faces[f1][0]]));
    float face1D = dotProduct(face1Normal, vertices[faces[f1][0]]);
    // get the 3 signed distances of the second face to the plane of the first face
    float distances1[3];
    for (int v = 0; v < 3; v++)
    {
        distances1[v] = dotProduct(face1Normal, vertices[faces[f2][v]]) - face1D;
    }
    float e = std::numeric_limits<float>::epsilon();
    // if the distances1 are all zero then the faces are coplanar, and we check for intersection separately
    if (std::abs(distances1[0]) < e &&
        std::abs(distances1[1]) < e &&
        std::abs(distances1[2]) < e)
    {
        // TODO: test 2D
    }

    // if all the distances1 are the same sign then the faces don't intersect
    if (distances1[0] > -e && distances1[1] > -e && distances1[2] > -e)
    {
        return false;
    }
    if (distances1[0] < e && distances1[1] < e && distances1[2] < e)
    {
        return false;
    }

    // calculate the 2 edges of the first face that will intersect the plane of the second face
    // one vertex will have a different sign to the other two
    int oddVertex1 = 2;
    if ((distances1[0] < e && distances1[1] > -e && distances1[2] > -e) || (distances1[0] > -e && distances1[1] < e && distances1[2] < e))
    {
        oddVertex1 = 0;
    }
    else if ((distances1[1] < e && distances1[0] > -e && distances1[2] > -e) || (distances1[1] > -e && distances1[0] < e && distances1[2] < e))
    {
        oddVertex1 = 1;
    }

    // calculate the plane for the second face
    Cartesian3 face2Normal = VertexToCartesian3(vertices[faces[f2][1]] - vertices[faces[f2][0]]).cross(
            VertexToCartesian3(
                    vertices[faces[f2][2]] - vertices[faces[f2][0]]));
    float face2D = dotProduct(face2Normal, vertices[faces[f2][0]]);
    // get the 3 signed distances of the first face to the plane of the second face
    float distances2[3];
    for (int v = 0; v < 3; v++)
    {
        distances2[v] = dotProduct(face2Normal, vertices[faces[f1][v]]) - face2D;
    }
    // if all the distances are the same sign then the faces don't intersect
    if (distances2[0] > 0 && distances2[1] > 0 && distances2[2] > 0)
    {
        return false;
    }
    if (distances2[0] < 0 && distances2[1] < 0 && distances2[2] < 0)
    {
        return false;
    }

    // calculate the 2 edges of the first face that will intersect the plane of the second face
    // one vertex will have a different sign to the other two
    int oddVertex2 = 2;
    if ((distances2[0] < e && distances2[1] > -e && distances2[2] > -e) || (distances2[0] > -e && distances2[1] < e && distances2[2] < e))
    {
        oddVertex2 = 0;
    }
    else if ((distances2[1] < e && distances2[0] > -e && distances2[2] > -e) || (distances2[1] > -e && distances2[0] < e && distances2[2] < e))
    {
        oddVertex2 = 1;
    }


    // check if the two faces are parallel
    if (face1Normal.cross(face2Normal).length() < 1e-6)
    {
        // we have an earlier test for them being coplanar
        // if they are parallel then they don't intersect
        return false;
    }
    // calculate the intersection line of the two planes
    Cartesian3 intersectionLine = face1Normal.cross(face2Normal);
    float maxComponent = std::max(std::abs(intersectionLine.x), std::max(std::abs(intersectionLine.y), std::abs(intersectionLine.z)));
    float face1Projected[3];
    float face2Projected[3];

    // using equation 6 from https://fileadmin.cs.lth.se/cs/Personal/Tomas_Akenine-Moller/pubs/tritri.pdf to project the vertices onto the intersection line
    if (std::abs(intersectionLine.x) == maxComponent)
    {
        for(int i = 0; i < 3; i++)
        {
            face1Projected[i] = vertices[faces[f1][i]].x;
            face2Projected[i] = vertices[faces[f2][i]].x;
        }
    }
    else if (std::abs(intersectionLine.y) == maxComponent)
    {
        for(int i = 0; i < 3; i++)
        {
            face1Projected[i] = vertices[faces[f1][i]].y;
            face2Projected[i] = vertices[faces[f2][i]].y;
        }
    }
    else
    {
        for(int i = 0; i < 3; i++)
        {
            face1Projected[i] = vertices[faces[f1][i]].z;
            face2Projected[i] = vertices[faces[f2][i]].z;
        }
    }

    // calculate the t values for the intersection line
    float f1t1, f1t2, f2t1, f2t2;
    f1t1 = face1Projected[oddVertex1] + (face1Projected[(oddVertex1 + 1) % 3] - face1Projected[oddVertex1]) * distances1[oddVertex1] / (distances1[oddVertex1] - distances1[(oddVertex1 + 1) % 3]);
    f1t2 = face1Projected[oddVertex1] + (face1Projected[(oddVertex1 + 2) % 3] - face1Projected[oddVertex1]) * distances1[oddVertex1] / (distances1[oddVertex1] - distances1[(oddVertex1 + 2) % 3]);
    f2t1 = face2Projected[oddVertex2] + (face2Projected[(oddVertex2 + 1) % 3] - face2Projected[oddVertex2]) * distances2[oddVertex2] / (distances2[oddVertex2] - distances2[(oddVertex2 + 1) % 3]);
    f2t2 = face2Projected[oddVertex2] + (face2Projected[(oddVertex2 + 2) % 3] - face2Projected[oddVertex2]) * distances2[oddVertex2] / (distances2[oddVertex2] - distances2[(oddVertex2 + 2) % 3]);
    // make sure t1 < t2
    if (f1t1 > f1t2)
    {
        std::swap(f1t1, f1t2);
    }
    if (f2t1 > f2t2)
    {
        std::swap(f2t1, f2t2);
    }
    // see if the intervals overlap (f1t1 is inside f2t or f2t1 is inside f1t)
    if ((f1t1 > f2t1 && f1t1 < f2t2) || (f2t1 > f1t1 && f2t1 < f1t2))
    {
        return true;
    }
    return false;
}

void ManifoldTester::testMultipleComponents()
{
    // we have one component if we can reach every vertex from every other vertex i.e. the graph is connected
    // we can do this by traversing the graph from one starting vertex and seeing if we can reach every other vertex
    std::vector<bool> visited(vertices.size(), false);

    // TODO: traverse the graph


    // if any vertex is still unvisited then we have multiple components
    if (std::any_of(visited.begin(), visited.end(), [](bool v)
    { return !v; }))
    {
        std::cout << "Mesh has multiple components, but each component is manifold" << std::endl;
    }
}

std::vector<Edge> ManifoldTester::getOneRing(int vertexIndex)
{
    // get the one ring of a vertex by finding all the edges that are incident to the vertex
    std::vector<Edge> oneRing;
    // for each face
    for (auto & face : faces)
    {
        // check if the face contains the vertex
        for (int v = 0; v < 3; v++)
        {
            if (face[v] == vertexIndex)
            {
                // add the edge to the one ring
                oneRing.push_back(Edge{face[(v + 1) % 3], face[(v + 2) % 3]});
            }
        }
    }
    return oneRing;
}

bool ManifoldTester::isSingleCycle(std::vector<Edge> edges)
{
    // follow the cycle
    int numTraversalsLeft;
    Edge currentEdge = edges[0];
    for (numTraversalsLeft = (int)edges.size() - 1; numTraversalsLeft > 0; numTraversalsLeft--)
    {
        // find the next edge - this is the edge that has the same start as the current edge's end
        bool foundNextEdge = false;
        for (auto & edge : edges)
        {
            if (edge.start == currentEdge.end)
            {
                currentEdge = edge;
                foundNextEdge = true;
                break;
            }
        }
        // if we didn't find the next edge then this isn't a single cycle
        if (!foundNextEdge)
        {
            return false;
        }
        // check we aren't currently at edge 0
        // if we are then we completed the cycle early so there must be a pinch point
        if (currentEdge == edges[0])
        {
            return false;
        }
    }
    // if we this edge connects back to the start then it is a single cycle
    return currentEdge.end == edges[0].start;
}