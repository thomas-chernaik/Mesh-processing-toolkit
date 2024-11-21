//
// Created by thomas on 08/11/24.
//

#include <algorithm>
#include <chrono>
#include <unordered_set>
#include <sstream>
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
#ifdef TIMING
    // start a timer
    auto start = std::chrono::high_resolution_clock::now();
#endif
    testPinchPoints();
#ifdef TIMING
    // stop the timer
    auto stop = std::chrono::high_resolution_clock::now();
    // output the time taken
    std::cout << "Pinch points took "
              << std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count()
              << " milliseconds" << std::endl;

    start = std::chrono::high_resolution_clock::now();
#endif
    // test for self intersections
    testSelfIntersections();
#ifdef TIMING
    stop = std::chrono::high_resolution_clock::now();
    std::cout << "Self intersections took "
              << std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count()
              << " milliseconds" << std::endl;
    start = std::chrono::high_resolution_clock::now();
#endif
    // test for multiple components
    testMultipleComponents();
#ifdef TIMING
    stop = std::chrono::high_resolution_clock::now();
    std::cout << "Multiple components took "
              << std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count()
              << " milliseconds" << std::endl;
#endif
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
#ifndef DEBUG
                // output the vertices of the faces
                std::cerr << "Face " << f1 << " vertices: " << vertices[faces[f1][0]] << " " << vertices[faces[f1][1]]
                          << " " << vertices[faces[f1][2]] << std::endl;
                std::cerr << "Face " << f2 << " vertices: " << vertices[faces[f2][0]] << " " << vertices[faces[f2][1]]
                          << " " << vertices[faces[f2][2]] << std::endl;
#endif
                //exit(-4);
            }
        }
    }
    std::cout << "Tested " << numTestedIntersections << " intersections" << std::endl;
}

bool ManifoldTester::testTriangleIntersection(int f1, int f2)
{
    // calculate the plane for each face
    // an early out is after the first plane calculation if the second face is on the same side of the plane
    Cartesian3 face1Normal = (vertices[faces[f1][1]] - vertices[faces[f1][0]]).cross(
            (vertices[faces[f1][2]] - vertices[faces[f1][0]]));
    face1Normal = face1Normal.normalise();
    float face1D = dotProduct(face1Normal, vertices[faces[f1][0]]);
    // get the 3 signed distances of the second face to the plane of the first face
    float distances1[3];
    for (int v = 0; v < 3; v++)
    {
        distances1[v] = dotProduct(face1Normal, vertices[faces[f2][v]]) - face1D;
    }
    float e = 1e-5;
    // if the distances1 are all zero then the faces are coplanar, and we check for intersection separately
    if (std::abs(distances1[0]) < e &&
        std::abs(distances1[1]) < e &&
        std::abs(distances1[2]) < e)
    {
        numTestedIntersections++;
        // project the vertices of the faces onto the plane of the first face
        // if the 2D triangles intersect then the 3D triangles intersect
        // create basis vectors for the plane from the triangle edges
        Cartesian3 p = VertexToCartesian3(vertices[faces[f1][0]]);
        Cartesian3 q = VertexToCartesian3(vertices[faces[f1][1]]);
        Cartesian3 r = VertexToCartesian3(vertices[faces[f1][2]]);
        Cartesian3 u = q - p;
        u = u.normalise();
        Cartesian3 v = r - p;
        Cartesian3 n = u.cross(v);
        n = n.normalise();
        Cartesian3 w = n.cross(u);
        w = w.normalise();

        Cartesian3 pPrime = {dotProduct(p - p, u), dotProduct(p - p, w), dotProduct(p - p, n)};
        Cartesian3 qPrime = {dotProduct(q - p, u), dotProduct(q - p, w), dotProduct(q - p, n)};
        Cartesian3 rPrime = {dotProduct(r - p, u), dotProduct(r - p, w), dotProduct(r - p, n)};
        Cartesian3 triangle1[3] = {pPrime, qPrime, rPrime};
        Cartesian3 p2 = VertexToCartesian3(vertices[faces[f2][0]]);
        Cartesian3 q2 = VertexToCartesian3(vertices[faces[f2][1]]);
        Cartesian3 r2 = VertexToCartesian3(vertices[faces[f2][2]]);
        Cartesian3 pPrime2 = {dotProduct(p2 - p, u), dotProduct(p2 - p, w), dotProduct(p2 - p, n)};
        Cartesian3 qPrime2 = {dotProduct(q2 - p, u), dotProduct(q2 - p, w), dotProduct(q2 - p, n)};
        Cartesian3 rPrime2 = {dotProduct(r2 - p, u), dotProduct(r2 - p, w), dotProduct(r2 - p, n)};
        Cartesian3 triangle2[3] = {pPrime2, qPrime2, rPrime2};

        // if they share an edge, we will handle this differently, checking for point in triangle
        int sharedVertices = 0;
        int verticesShared[2];
        int verticesShared2[2];
        for (int v1 = 0; v1 < 3; v1++)
        {
            for (int v2 = 0; v2 < 3; v2++)
            {
                if (faces[f1][v1] == faces[f2][v2])
                {
                    verticesShared[sharedVertices] = v1;
                    verticesShared2[sharedVertices] = v2;
                    sharedVertices++;
                }
            }
        }
        if (sharedVertices >= 2)
        {
            // we only need to test that either of the non-shared vertices is in the other triangle
            for (int i = 0; i < 3; i++)
            {
                if (i != verticesShared[0] && i != verticesShared[1])
                {
                    if (TriangleContainsVertex(triangle2[0], triangle2[1], triangle2[2], triangle1[i]))
                    {
                        return true;
                    }
                }
                if (i != verticesShared2[0] && i != verticesShared2[1])
                {
                    if (TriangleContainsVertex(triangle1[0], triangle1[1], triangle1[2], triangle2[i]))
                    {
                        return true;
                    }
                }
            }
            return false;
        }

        for (int v1 = 0; v1 < 3; v1++)
        {
            for (int v2 = 0; v2 < 3; v2++)
            {
                if (EdgesIntersect(triangle1[v1], triangle1[(v1 + 1) % 3], triangle2[v2], triangle2[(v2 + 1) % 3]))
                {
                    // print out the intersecting edges
                    std::cout << "Edges (" << vertices[faces[f1][v1]].x << ", " << vertices[faces[f1][v1]].y << ", "
                              << vertices[faces[f1][v1]].z << ") to ("
                              << vertices[faces[f1][(v1 + 1) % 3]].x << ", " << vertices[faces[f1][(v1 + 1) % 3]].y
                              << ", "
                              << vertices[faces[f1][(v1 + 1) % 3]].z << ") and ("
                              << vertices[faces[f2][v2]].x << ", " << vertices[faces[f2][v2]].y << ", "
                              << vertices[faces[f2][v2]].z << ") to ("
                              << vertices[faces[f2][(v2 + 1) % 3]].x << ", " << vertices[faces[f2][(v2 + 1) % 3]].y
                              << ", "
                              << vertices[faces[f2][(v2 + 1) % 3]].z << ") intersect" << std::endl;

                    return true;
                }
            }
        }
        // if the edges don't intersect than either one triangle contains the other, or they don't intersect
        for (int v1 = 0; v1 < 3; v1++)
        {
            if (TriangleContainsVertex(triangle2[0], triangle2[1], triangle2[2], triangle1[v1]))
            {
                std::cout << "Triangle 1 contains triangle 2" << std::endl;
                return true;
            }
            if (TriangleContainsVertex(triangle1[0], triangle1[1], triangle1[2], triangle2[v1]))
            {
                std::cout << "Triangle 2 contains triangle 1" << std::endl;
                return true;
            }
        }
        return false;
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
    if ((distances1[0] < e && distances1[1] > -e && distances1[2] > -e) ||
        (distances1[0] > -e && distances1[1] < e && distances1[2] < e))
    {
        oddVertex1 = 0;
    } else if ((distances1[1] < e && distances1[0] > -e && distances1[2] > -e) ||
               (distances1[1] > -e && distances1[0] < e && distances1[2] < e))
    {
        oddVertex1 = 1;
    }

    // calculate the plane for the second face
    Cartesian3 face2Normal = (vertices[faces[f2][1]] - vertices[faces[f2][0]]).cross(
            (vertices[faces[f2][2]] - vertices[faces[f2][0]]));
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

    numTestedIntersections++;

    // calculate the 2 edges of the first face that will intersect the plane of the second face
    // one vertex will have a different sign to the other two
    int oddVertex2 = 2;
    if ((distances2[0] < e && distances2[1] > -e && distances2[2] > -e) ||
        (distances2[0] > -e && distances2[1] < e && distances2[2] < e))
    {
        oddVertex2 = 0;
    } else if ((distances2[1] < e && distances2[0] > -e && distances2[2] > -e) ||
               (distances2[1] > -e && distances2[0] < e && distances2[2] < e))
    {
        oddVertex2 = 1;
    }


    // check if the two faces are parallel
    if (face1Normal.cross(face2Normal).length() < e)
    {
        // we have an earlier test for them being coplanar
        // if they are parallel then they don't intersect
        return false;
    }
    // calculate the intersection line of the two planes
    Cartesian3 intersectionLine = face1Normal.cross(face2Normal);
    float maxComponent = std::max(std::abs(intersectionLine.x),
                                  std::max(std::abs(intersectionLine.y), std::abs(intersectionLine.z)));
    float face1Projected[3];
    float face2Projected[3];

//    // using equation 6 from https://fileadmin.cs.lth.se/cs/Personal/Tomas_Akenine-Moller/pubs/tritri.pdf to project the vertices onto the intersection line
//    if (std::abs(intersectionLine.x) == maxComponent)
//    {
//        for (int i = 0; i < 3; i++)
//        {
//            face1Projected[i] = vertices[faces[f1][i]].x;
//            face2Projected[i] = vertices[faces[f2][i]].x;
//        }
//    } else if (std::abs(intersectionLine.y) == maxComponent)
//    {
//        for (int i = 0; i < 3; i++)
//        {
//            face1Projected[i] = vertices[faces[f1][i]].y;
//            face2Projected[i] = vertices[faces[f2][i]].y;
//        }
//    } else
//    {
//        for (int i = 0; i < 3; i++)
//        {
//            face1Projected[i] = vertices[faces[f1][i]].z;
//            face2Projected[i] = vertices[faces[f2][i]].z;
//        }
//    }


    // using equation 3 from https://fileadmin.cs.lth.se/cs/Personal/Tomas_Akenine-Moller/pubs/tritri.pdf to project the vertices onto the intersection line
    // get O, a point on the intersection line
    Cartesian3 O = {0, 0, 0};
    // o = p1 + d * ((p1 - p2) dot n2) / (d dot n2)
    Cartesian3 p1 = VertexToCartesian3(vertices[faces[f1][0]]);
    Cartesian3 p2 = VertexToCartesian3(vertices[faces[f2][0]]);
    Cartesian3 d = intersectionLine;
    Cartesian3 n2 = face2Normal;
    O = p1 + d * dotProduct(p1 - p2, n2) / dotProduct(d, n2);
    // project the vertices of the first face onto the intersection line
    for (int i = 0; i < 3; i++)
    {
        face1Projected[i] = dotProduct(intersectionLine, VertexToCartesian3(vertices[faces[f1][i]]) - O);
    }
    // project the vertices of the second face onto the intersection line
    for (int i = 0; i < 3; i++)
    {
        face2Projected[i] = dotProduct(intersectionLine, VertexToCartesian3(vertices[faces[f2][i]]) - O);
    }
    // calculate the t values for the intersection line
    float f1t1, f1t2, f2t1, f2t2;
    f1t1 = face1Projected[oddVertex1] +
           (face1Projected[(oddVertex1 + 1) % 3] - face1Projected[oddVertex1]) * distances1[oddVertex1] /
           (distances1[oddVertex1] - distances1[(oddVertex1 + 1) % 3]);
    f1t2 = face1Projected[oddVertex1] +
           (face1Projected[(oddVertex1 + 2) % 3] - face1Projected[oddVertex1]) * distances1[oddVertex1] /
           (distances1[oddVertex1] - distances1[(oddVertex1 + 2) % 3]);
    f2t1 = face2Projected[oddVertex2] +
           (face2Projected[(oddVertex2 + 1) % 3] - face2Projected[oddVertex2]) * distances2[oddVertex2] /
           (distances2[oddVertex2] - distances2[(oddVertex2 + 1) % 3]);
    f2t2 = face2Projected[oddVertex2] +
           (face2Projected[(oddVertex2 + 2) % 3] - face2Projected[oddVertex2]) * distances2[oddVertex2] /
           (distances2[oddVertex2] - distances2[(oddVertex2 + 2) % 3]);
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
    if ((f1t1 > f2t1 + e && f1t1 < f2t2 - e) || (f2t1 > f1t1 + e && f2t1 < f1t2 - e))
    {
        return true;
    }
    return false;
}

void ManifoldTester::testMultipleComponents()
{
    // we have one component if we can reach every face from every other face via a paired edge
    // we can do this by traversing the graph from one starting vertex and seeing if we can reach every other vertex
    // create a vector of visited faces
    std::vector<bool> visited(faces.size(), false);

    int numComponents = 0;
    components.clear();

    while (std::any_of(visited.begin(), visited.end(), [](bool v)
    { return !v; }))
    {
        numComponents++;
        // create a stack of faces to visit
        components.emplace_back();
        std::vector<int> toVisit;
        // add the first unvisited face to the stack
        toVisit.push_back(std::find(visited.begin(), visited.end(), false) - visited.begin());
        // while there are still vertices to visit
        while (!toVisit.empty())
        {
            // pop the last face from the stack
            int currentFace = toVisit.back();
            toVisit.pop_back();
            // if we have already visited this face then skip it
            if (visited[currentFace])
            {
                continue;
            }
            // add the face to the component
            components[numComponents - 1].push_back(currentFace);
            // mark the face as visited
            visited[currentFace] = true;
            // get the three edges of the face
            edge edges[3] = {currentFace * 3, currentFace * 3 + 1, currentFace * 3 + 2};
            // get the three other half edges
            int otherHalfEdges[3] = {otherHalf[edges[0]], otherHalf[edges[1]], otherHalf[edges[2]]};
            // for each other half edge
            for (int i = 0; i < 3; i++)
            {
                // if the other half edge is unpaired then skip it
                if (otherHalfEdges[i] == -1)
                {
                    continue;
                }
                // get the face that the other half edge is paired to
                int nextFace = otherHalfEdges[i] / 3;
                // if we haven't visited the next face then add it to the stack
                if (!visited[nextFace])
                {
                    toVisit.push_back(nextFace);
                }
            }
        }

    }
    // if there is more than one component then output a warning
    if (numComponents > 1)
    {
        std::cerr << "Mesh has " << numComponents << " components" << std::endl;
    }
}

std::vector<Edge> ManifoldTester::getOneRing(int vertexIndex)
{
    // get the one ring of a vertex by finding all the edges that are incident to the vertex
    std::vector<Edge> oneRing;
    // for each face
    for (auto &face: faces)
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
    for (numTraversalsLeft = (int) edges.size() - 1; numTraversalsLeft > 0; numTraversalsLeft--)
    {
        // find the next edge - this is the edge that has the same start as the current edge's end
        bool foundNextEdge = false;
        for (auto &edge: edges)
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

bool ManifoldTester::EdgesIntersect(Cartesian3 v1, Cartesian3 v2, Cartesian3 v3, Cartesian3 v4)
{
    // check if the two line segments share a point
    if (v1 == v3 || v1 == v4 || v2 == v3 || v2 == v4)
    {
        return false;
    }
    // calculate the direction vectors of the two line segments
    Cartesian3 direction1 = v2 - v1;
    Cartesian3 direction2 = v4 - v3;

    // calculate the denominator of the t values
    float denominator = 1.f / (direction1.x * direction2.y - direction1.y * direction2.x);
    // if the denominator is zero then the lines are parallel
    if (std::abs(denominator) < std::numeric_limits<float>::epsilon() * 5.0f)
    {
        return false;
    }

    // calculate the numerators
    Cartesian3 u = v1 - v3;
    float numerator1 = u.x * direction2.y - u.y * direction2.x;
    float numerator2 = u.x * direction1.y - u.y * direction1.x;

    float t1 = numerator1 * denominator;
    float t2 = numerator2 * denominator;

    float e = std::numeric_limits<float>::epsilon() * 5.0f;
    // if the t values are between 0 and 1 then the line segments intersect
    return t1 > e && t1 < 1 - e && t2 > e && t2 < 1 - e;

}

// takes a triangle in 2D and a point in 2D and returns if the point is inside the triangle
// Cartesian3s are used but the z component is ignored
bool ManifoldTester::TriangleContainsVertex(Cartesian3 &p, Cartesian3 &q, Cartesian3 &r,
                                            Cartesian3 &point)
{
    // do the half plane test
    Cartesian3 pr = r - p;
    Cartesian3 rq = q - r;
    Cartesian3 qp = p - q;
    Cartesian3 pPoint = point - p;
    Cartesian3 qPoint = point - q;
    Cartesian3 rPoint = point - r;
    // find the normals of the lines
    Cartesian3 prNormal = Cartesian3(-pr.y, pr.x, 0);
    Cartesian3 rqNormal = Cartesian3(-rq.y, rq.x, 0);
    Cartesian3 qpNormal = Cartesian3(-qp.y, qp.x, 0);

    float epsilon = -std::numeric_limits<float>::epsilon() * 5.0f;

    if (dotProduct(prNormal, pPoint) >= epsilon && dotProduct(rqNormal, rPoint) >= epsilon &&
        dotProduct(qpNormal, qPoint) >= epsilon)
    {
        return true;
    }
    return false;
}

std::vector<int> ManifoldTester::CalculateGenus()
{
    // in any orientable mesh/polyhedron, the Euler characteristic is given by V - E + F = 2 - 2g
    // therefore, g = (2 - V + E - F) / 2
    // V can be found as the vertices.size()
    // E can be found as the size of otherHalf divided by 2
    // F can be found as the faces.size()

    // we will calculate the genus separately for each component
    std::vector<int> genuses(components.size());
    int currentComponent = 0;
    for (auto &component: components)
    {
        // the number of faces is the size of the component
        int iFacesInComponent = component.size();
        // the number of vertices is the number of unique vertices in the component
        std::unordered_set<int> verticesInComponent;
        for (int face: component)
        {
            verticesInComponent.insert(faces[face][0]);
            verticesInComponent.insert(faces[face][1]);
            verticesInComponent.insert(faces[face][2]);
        }
        int iVerticesInComponent = verticesInComponent.size();
        // the number of edges is the number of faces * 3 / 2
        int edgesInComponent = iFacesInComponent * 3 / 2;
        int genusesInComponent = (2 - iVerticesInComponent + edgesInComponent - iFacesInComponent) / 2;
        genuses[currentComponent++] = genusesInComponent;
    }
    return genuses;
}

std::vector<int> ManifoldTester::getOneRingVertices(int vertexIndex)
{
    std::vector<int> oneRingVertices;

    // get the one ring of the vertex
    for (auto &face: faces)
    {
        for (int v = 0; v < 3; v++)
        {
            if (face[v] == vertexIndex)
            {
                oneRingVertices.push_back(face[(v + 1) % 3]);
                oneRingVertices.push_back(face[(v + 2) % 3]);
            }
        }
    }
    return oneRingVertices;
}

void ManifoldTester::readFileDiredge(const std::string &filename)
{
    // clear the vertices, faces, directedEdges and otherHalf vectors
    vertices.clear();
    faces.clear();
    directedEdges.clear();
    otherHalf.clear();
    // open the file to read
    std::ifstream fileReader(filename);
    // check the file opened
    if (!fileReader.is_open())
    {
        std::cerr << "File failed to open" << std::endl;
        exit(-2);
    }
    // read out the comments - these are lines that start with a #
    std::string line;
    while (std::getline(fileReader, line))
    {
        if (line[0] != '#')
        {
            break;
        }
    }
    int currentIndex;
    // we are reading vertices while the line's first word is Vertex
    while (line.substr(0, 6) == "Vertex")
    {
        // create a new vertex
        Vertex vertex{};
        // split up the line by whitespace
        // [0] = "Vertex", [1] = index, [2] = x, [3] = y, [4] = z
        std::istringstream lineStream(line);
        lineStream >> line; // "Vertex"
        lineStream >> currentIndex;
        lineStream >> vertex.x;
        lineStream >> vertex.y;
        lineStream >> vertex.z;
        // resize the vertices to be the maximum of the currentIndex + 1 or the current size, with a default value of null vertex
        vertices.resize(std::max(currentIndex + 1, (int) vertices.size()), nullVertex);
        // add the vertex to the list
        vertices[currentIndex] = vertex;
        // read the next line
        std::getline(fileReader, line);
    }
    // check we don't have any null vertices remaining
    for (Vertex vertex: vertices)
    {
        if (vertex == nullVertex)
        {
            std::cerr << "File is missing vertices" << std::endl;
            exit(-3);
        }
    }
    // we are now reading first directed edges while the line's first word is FirstDirectedEdge
    while (line.substr(0, 17) == "FirstDirectedEdge")
    {
        int thisEdge;
        // split up the line by whitespace
        // [0] = "FirstDirectedEdge", [1] = index, [2] = edge index
        std::istringstream lineStream(line);
        lineStream >> line; // "FirstDirectedEdge"
        lineStream >> currentIndex;
        lineStream >> thisEdge;
        // resize the directedEdges to be the maximum of the currentIndex + 1 or the current size, with a default value of -1
        directedEdges.resize(std::max(currentIndex + 1, (int) directedEdges.size()), -1);
        // add the edge to the list
        directedEdges[currentIndex] = thisEdge;
        // read the next line
        std::getline(fileReader, line);
    }
    // check we don't have any null directed edges remaining
    for (int edge: directedEdges)
    {
        if (edge == -1)
        {
            std::cerr << "File is missing directed edges" << std::endl;
            exit(-3);
        }
    }
    // check we have as many directed edges as we have vertices
    if (directedEdges.size() != vertices.size())
    {
        std::cerr << "File has the wrong number of directed edges" << std::endl;
        exit(-3);
    }
    // now we read the faces
    // we are reading faces while the line's first word is Face
    while (line.substr(0, 4) == "Face")
    {
        // create a new face
        Face face{};
        // split up the line by whitespace
        // [0] = "Face", [1] = index, [2] = vertex1, [3] = vertex2, [4] = vertex3
        std::istringstream lineStream(line);
        lineStream >> line; // "Face"
        lineStream >> currentIndex;
        lineStream >> face[0];
        lineStream >> face[1];
        lineStream >> face[2];
        // check that the face indices are valid
        if (face[0] < 0 || face[0] > vertices.size() || face[1] < 0 || face[1] > vertices.size() || face[2] < 0 ||
            face[2] > vertices.size())
        {
            std::cerr << "Face " << currentIndex << " has an invalid vertex index" << std::endl;
            exit(-3);
        }
        // resize the faces to be the maximum of the currentIndex + 1 or the current size, with a default value of null face
        faces.resize(std::max(currentIndex + 1, (int) faces.size()), nullFace);
        // add the face to the list
        faces[currentIndex] = face;
        // read the next line
        std::getline(fileReader, line);
    }
    // check we don't have any null faces remaining
    for (Face face: faces)
    {
        if (face == nullFace)
        {
            std::cerr << "File is missing faces" << std::endl;
            exit(-3);
        }
    }
    // now we read the other halfs
    // we are reading other halfs while the line's first word is OtherHalf
    while (line.substr(0, 9) == "OtherHalf")
    {
        int thisHalf;
        // split up the line by whitespace
        // [0] = "OtherHalf", [1] = index, [2] = half index
        std::istringstream lineStream(line);
        lineStream >> line; // "OtherHalf"
        lineStream >> currentIndex;
        lineStream >> thisHalf;
        // resize the otherHalf to be the maximum of the currentIndex + 1 or the current size, with a default value of -2 - a -1 value means no other half but is valid
        otherHalf.resize(std::max(currentIndex + 1, (int) otherHalf.size()), -2);
        // add the half to the list
        otherHalf[currentIndex] = thisHalf;
        // read the next line
        std::getline(fileReader, line);
    }
    // check we don't have any null other halfs remaining
    for (int half: otherHalf)
    {
        if (half == -2)
        {
            std::cerr << "File is missing other halfs" << std::endl;
            exit(-3);
        }
    }
    // check we don't have anything else left but white space
    char whiteSpaceCheckerChar;
    while (fileReader.get(whiteSpaceCheckerChar))
    {
        if (!isspace(whiteSpaceCheckerChar))
        {
            std::cerr << "Finished reading the file but there are still non-whitespace characters left" << std::endl;
            exit(-3);
        }
    }
    // close the file
    fileReader.close();


}
