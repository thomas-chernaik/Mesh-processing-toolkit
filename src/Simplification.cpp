//
// Created by thomas on 18/11/24.
//

#include <algorithm>
#include <cmath>
#include "Simplification.h"



void Simplification::simplifyMesh(int maxIterations)
{
    // TODO: work out exit conditions
    //  for now, if there is no vertex that can be removed without breaking the eulerian condition, then we exit
    bool removed;
    curvatures.clear();
    // add the whole mesh to the component - every face
    components.emplace_back();
    components[0].resize(faces.size());
    for (int i = 0; i < faces.size(); i++)
    {
        components[0][i] = i;
    }
    startingGenus = CalculateGenus()[0];
//    std::cout << "Starting genus: " << startingGenus << std::endl;
    // generate the set of curvatures
    generateCurvatures();
    if(maxIterations > vertices.size() - 4)
    {
        std::cout << "There are only " << vertices.size() << " so it wouldn't make sense to simplify more than that" << std::endl;
        exit(-1);
    }
    std::cout << "Simplifying the mesh" << std::endl;
    for (int i = 0; i < maxIterations; i++)
    {
        //cleanUpNonManifoldEdges();
        // -debug - write the mesh to a file with the iteration number
//        // if its a multiple of 100
//        if (i % 1000 == 0)
//            writeObjFile("iteration" + std::to_string(i));
//        if (i == 11536)
//            writeRepairedMeshTri("iteration" + std::to_string(i));

        // progress bar
        int percentage = (i * 100) / maxIterations;
        printProgress(percentage);

        //writeRepairedMeshTri("iteration" + std::to_string(i));
        removed = false;
        int currentVertex = 0;
        int smallestCurvature = findSmallestCurvature(0);
        // print out the vertex with the smallest curvature
        while (!removed && currentVertex < vertices.size() - i)
        {
            //std::cout << "Smallest curvature: " << vertices[smallestCurvature] << std::endl;
            try
            {
                if (!removeVertex(smallestCurvature))
                {
                    //std::cout << "backtracking" << std::endl;
                    // if we can't remove the vertex then we need to backtrack
                    backtrack();
                    currentVertex++;
                    smallestCurvature = findSmallestCurvature(currentVertex);
                } else
                {
                    removed = true;
                    // update the curvatures
                    updateCurvatures();
                    // repair the mesh
                    //cleanUpNonManifoldEdges();
                }
            }
            catch(std::exception &e)
            {
                std::cout << std::endl;
                std::cout << "Cannot continue simplification because of error: " << e.what() << std::endl;
                std::cout << "Simplified mesh for a total vertex reduction of " << i << std::endl;
                //cleanUpNonManifoldEdges();
                writeObjFile("final");
                return;
            }
        }
        if(!removed)
        {
            //std::cout << "No vertex can be removed without breaking the eulerian condition" << std::endl;
            break;
        }
//        // update the component
//        // remove the faces in holeFaces from the component
//        int numRemoved = 0;
//        for(int j = 0; j < components.size(); j++)
//        {
//            // if the component is in hole faces, swap it with the last - numRemoved component
//            if (holeFaces.find(j - numRemoved) != holeFaces.end())
//            {
//                std::swap(components[j - numRemoved], components[components.size() - numRemoved - 1]);
//                numRemoved++;
//            }
//
//        }
//        // remove the last numRemoved components
//        components.resize(components.size() - numRemoved);
//        // now add in the new faces
//        // these are faces indexed at the end of the faces vector, so faces.size() - facesAdded to faces.size()
//        for(int j = faces.size() - facesAdded; j < faces.size(); j++)
//        {
//            components[0].push_back(j);
//        }


    }
    printProgress(100);
    std::cout << "Simplified mesh for a total vertex reduction of " << maxIterations << std::endl;
    // -debug - write an obj file
    //cleanUpNonManifoldEdges();
    writeObjFile("final");
}

int Simplification::findSmallestCurvature(int n)
{
    if(n == 0)
    {
        // return the smallest curvature
        // get the key of the minimum value in curvatures
        int minElement = -1;
        float minValue = std::numeric_limits<float>::max();
        // iterate across the curvatures dictionary
        for (auto &curvature: curvatures)
        {
            // if the curvature is less than the current minimum value
            if (curvature.second < minValue)
            {
                // set the minimum value to the curvature
                minValue = curvature.second;
                // set the minimum element to the key
                minElement = curvature.first;
            }
        }
        return minElement;
    }
    // return the nth smallest curvature
    // get the curvatures in a vector

    std::vector<float> curvaturesVector(vertices.size(), std::numeric_limits<float>::max());
    // set the indices to 0 - vertices.size() - 1
    std::vector<int> curvaturesIndices;
    curvaturesIndices.reserve(vertices.size());
    for(int i = 0; i < vertices.size(); i++)
    {
        curvaturesIndices.push_back(i);
    }
    // iterate across the curvatures dictionary and update the vectors where needed
    int i = 0;
    for (auto &curvature: curvatures)
    {
        curvaturesVector[curvature.first] = curvature.second;
    }
    // sort the indices by the curvatures
    std::sort(curvaturesIndices.begin(), curvaturesIndices.end(), [&curvaturesVector](int i1, int i2)
    {
        return curvaturesVector[i1] < curvaturesVector[i2];
    });
    // return the nth smallest curvature
    return curvaturesIndices[n];

}

void Simplification::backtrack()
{
    // remove the last added faces
    faces.resize(faces.size() - facesAdded);
    // add the removed faces back
    for (auto &edge: holeEdges)
    {
        faces.push_back({edge.start, edge.end, removedVertexIndex});
    }

}

bool Simplification::removeVertex(int vertexIndex)
{
    // store the vertex that is being removed
    removedVertexIndex = vertexIndex;
    holeEdges.clear();
    //holeFaces.clear();
    int faceIndex = 0;
    // remove any faces that contain the vertex, and add the edges of the hole to the holeEdges vector
    for (int i = 0; i < faces.size(); i++)
    {
        Face face = faces[i];
        if (face[0] == vertexIndex)
        {
            //holeFaces.insert(faceIndex);
            holeEdges.push_back({face[1], face[2]});
            // remove the face
            faces[i] = faces.back();
            faces.pop_back();
            i--;
        } else if (face[1] == vertexIndex)
        {
            //holeFaces.insert(faceIndex);
            holeEdges.push_back({face[2], face[0]});
            // remove the face
            faces[i] = faces.back();
            faces.pop_back();
            i--;
        } else if (face[2] == vertexIndex)
        {
            //holeFaces.insert(faceIndex);
            holeEdges.push_back({face[0], face[1]});
            // remove the face
            faces[i] = faces.back();
            faces.pop_back();
            i--;
        }
        faceIndex++;
    }
    int previousFaces = (int) faces.size();
    //std::cout << "filling hole with " << holeEdges.size() << " edges" << std::endl;
    // triangulate the hole
    triangulateHole(holeEdges);
    facesAdded = (int) faces.size() - previousFaces;

    // check that the mesh is still eulerian
    return testNonManifoldEdges();
}


float Simplification::findCurvature(int vertexIndex)
{
    return findGaussianCurvature(vertexIndex);
    // get the mean and gaussian curvature
    float meanCurvature = findMeanCurvature(vertexIndex);
    float gaussianCurvature = findGaussianCurvature(vertexIndex);
    // plug them into the quadratic equation for k1 and k2
    float k1 = meanCurvature + std::sqrt(meanCurvature * meanCurvature - gaussianCurvature);
    float k2 = meanCurvature - std::sqrt(meanCurvature * meanCurvature - gaussianCurvature);
    // return the larger of the two
    if (k1 > k2)
    {
        return k1;
    }
    return k2;
}

float Simplification::findMeanCurvature(int vertexIndex)
{
    // get the discrete laplace beltrami operator
    Cartesian3 laplaceBeltrami = computeLaplaceBeltrami(vertexIndex);
    // the mean curvature is half the magnitude of the laplace beltrami operator
    return laplaceBeltrami.lengthSqrt() / 2;
}

float Simplification::findGaussianCurvature(int vertexIndex)
{
    // this is the "angular deficit"
    // its equal to 1 / Ai (what is Ai?) * (2pi - sum of angles around the vertex)
    // get the one ring vertices
    //std::vector<int> oneRingVertices = getOneRingVertices(vertexIndex);
    // -debug test
    // get the one ring edges
    std::vector<Edge> oneRingEdges = getOneRing(vertexIndex);
    // order the one ring edges
    std::vector<Edge> orderedOneRingEdges;
    orderedOneRingEdges.push_back(oneRingEdges[0]);
    while (orderedOneRingEdges.size() < oneRingEdges.size())
    {
        // get the last vertex of the last edge in the ordered boundary
        int lastVertex = orderedOneRingEdges.back().end;
        // find the next edge that has the last vertex as the start
        for (auto &edge: oneRingEdges)
        {
            if (edge.start == lastVertex)
            {
                orderedOneRingEdges.push_back(edge);
                break;
            }
        }
    }
    // get the first vertex of each edge
    std::vector<int> oneRingVertices;
    oneRingVertices.reserve(orderedOneRingEdges.size());
for (auto &edge: orderedOneRingEdges)
    {
        oneRingVertices.push_back(edge.start);
    }
    // get the angles around the vertex
    float sumOfAngles = 0.f;
    for (int i = 0; i < oneRingVertices.size(); i++)
    {
        // get the two edge vectors, both coming from the vertex
        Cartesian3 edge1 = vertices[oneRingVertices[i]] - vertices[vertexIndex];
        Cartesian3 edge2 = vertices[oneRingVertices[(i + 1) % oneRingVertices.size()]] - vertices[vertexIndex];
        // get the angle between the two edges
        sumOfAngles += getAngleBetweenVectors(edge1, edge2);
    }
//    // if the sum of the angles is greater than 2pi then smth went wrong
//    if (sumOfAngles >= 2.f * M_PI + 0.01f)
//    {
//        std::cerr << "Sum of angles is greater than 2pi around vertex " << vertexIndex << " : " << vertices[vertexIndex]
//                  << std::endl;
//        // output the one ring vertices
//        for (auto &vertex: oneRingVertices)
//        {
//            std::cout << vertices[vertex] << std::endl;
//        }
//        exit(-3);
//    }
    // return the gaussian curvature
    return std::abs(2.f * M_PI - sumOfAngles);
}


void Simplification::triangulateHole(std::vector<Edge> &boundary)
{
    // order the boundary edges so that they form a single cycle
    std::vector<Edge> orderedBoundary;
    orderedBoundary.push_back(boundary[0]);
    while (orderedBoundary.size() < boundary.size())
    {
        // get the last vertex of the last edge in the ordered boundary
        int lastVertex = orderedBoundary.back().end;
        // find the next edge that has the last vertex as the start
        for (auto &edge: boundary)
        {
            if (edge.start == lastVertex)
            {
                orderedBoundary.push_back(edge);
                break;
            }
        }
    }
//    // do a fan triangulation
//    for(int i = 1; i < orderedBoundary.size() - 1; i++)
//    {
//        faces.push_back({orderedBoundary[0].start, orderedBoundary[i].start, orderedBoundary[i + 1].start});
//    }
//    return;
    // compute the normal of the hole
    Cartesian3 holeNormal = getHoleNormal(orderedBoundary);
    // while the boundary has more than 3 edges
    while (orderedBoundary.size() > 3)
    {
        // get two edges that aren't reflex angles
        int goodEdgeIndex = getConcaveEdge(orderedBoundary, holeNormal);
        // get the two edges
        Edge edge1 = orderedBoundary[goodEdgeIndex];
        Edge edge2 = orderedBoundary[(goodEdgeIndex + 1) % orderedBoundary.size()];
        // add the triangle to the faces
        faces.push_back({edge1.start, edge1.end, edge2.end});
        // remove the edges from the boundary and add the new edge
        // we need to maintain the order of the boundary, so we remove one edge and swap the other with the new edge
        // swap the new edge in with the second edge
        orderedBoundary[(goodEdgeIndex + 1) % orderedBoundary.size()] = {edge1.start, edge2.end};
        // remove the first edge
        orderedBoundary.erase(orderedBoundary.begin() + goodEdgeIndex);
    }
    // make sure the last triangle is in fact a triangle
    if (orderedBoundary[0].end != orderedBoundary[1].start || orderedBoundary[1].end != orderedBoundary[2].start ||
        orderedBoundary[2].end != orderedBoundary[0].start)
    {
        throw std::runtime_error("Last triangle is not a triangle");
    }
    // make sure the last triangle isn't degenerate
    // if any of the start vertices are the same then we have a degenerate triangle
    if (orderedBoundary[0].start == orderedBoundary[1].start || orderedBoundary[1].start == orderedBoundary[2].start ||
        orderedBoundary[2].start == orderedBoundary[0].start)
    {
        throw std::runtime_error("Last triangle is degenerate");
    }
    // make sure the last triangle isn't a duplicate of an existing face
    for (auto &face: faces)
    {
        if ((face[0] == orderedBoundary[0].start && face[1] == orderedBoundary[1].start &&
             face[2] == orderedBoundary[2].start) ||
            (face[0] == orderedBoundary[1].start && face[1] == orderedBoundary[2].start &&
             face[2] == orderedBoundary[0].start) ||
            (face[0] == orderedBoundary[2].start && face[1] == orderedBoundary[0].start &&
             face[2] == orderedBoundary[1].start))
        {
            throw std::runtime_error("Last triangle is a duplicate");
        }
    }
    // add the last triangle
    faces.push_back({orderedBoundary[0].start, orderedBoundary[0].end, orderedBoundary[1].end});

}


int Simplification::findSmallestAngle(const std::vector<Edge> &boundary)
{
    int smallestAngleIndex = 0;
    float smallestAngle = 200.f;
    // loop through the boundary and test the edge and the next edge
    for (int i = 0; i < boundary.size(); i++)
    {
        // get the two edges
        Edge edge1 = boundary[i];
        Edge edge2 = boundary[(i + 1) % boundary.size()];
        // reverse the edges to test the angle in the other direction
        Edge edge1Reverse = {edge1.end, edge1.start};
        Edge edge2Reverse = {edge2.end, edge2.start};
        float angleBetween = getAngleBetweenEdges(edge1Reverse, edge2Reverse);
        // if the angle between the two edges is smaller than the smallest angle, then update the smallest angle
        if (angleBetween < smallestAngle)
        {
            smallestAngle = angleBetween;
            smallestAngleIndex = i;
        }
    }
    return smallestAngleIndex;
}

Cartesian3 Simplification::computeLaplaceBeltrami(int vertexIndex)
{
    // get the centre of gravity of the vertex and its one ring
    Vertex centreOfGravity = getCentreOfGravity(getOneRingVertices(vertexIndex));
    // get the vertex
    Vertex vertex = vertices[vertexIndex];
    // subtract the vertex from the centre of gravity
    Cartesian3 result = centreOfGravity - vertex;
    return result;
}

float Simplification::getAngleBetweenVectors(Cartesian3 vector1, Cartesian3 vector2)
{
    // the angle between two vectors is the acos of the dot product of the two vectors divided by the product of their magnitudes
    return std::acos(dotProduct(vector1, vector2) / (vector1.lengthSqrt() * vector2.lengthSqrt()));
}

bool Simplification::isEulerian()
{
    // calculate the genus of the mesh
    // the built in function won't work since we may have multiple components
    // genus = 1 - (V - E + F) / 2
    // F = faces.size()
    // E can be calculated as the number of unique undirected edges in the faces
    // V can be calculated as the number of unique vertices in the faces
    int F = (int) faces.size();
    std::unordered_set<int> uniqueVertices;
    std::unordered_set<Edge, EdgeHash> uniqueEdges;
    for (auto &face: faces)
    {
        for (int i = 0; i < 3; i++)
        {
            uniqueVertices.insert(face[i]);
            uniqueEdges.insert({face[i], face[(i + 1) % 3]});
        }
    }
    int E = (int) uniqueEdges.size();
    int V = (int) uniqueVertices.size();
    int genus = 1 - (V - E + F) / 2;
    // if the genus isn't the same as the starting genus then the mesh is not eulerian
    if (genus != startingGenus)
    {
        //std::cerr << "Genus is not the same as the starting genus" << std::endl;
        //exit(-2);
        return false;
    }

    return true;
}

std::vector<int> Simplification::getOneRingVertices(int vertexIndex)
{
    // get the one ring edges
    std::vector<Edge> oneRingEdges = getOneRing(vertexIndex);
    // get the one ring vertices
    std::vector<int> oneRingVertices;
    oneRingVertices.reserve(oneRingEdges.size());
for (auto &edge: oneRingEdges)
    {
        oneRingVertices.push_back(edge.end);
    }
    return oneRingVertices;
//
//    std::vector<int> oneRing;
//    // the first one ring vertex is the other side of the first directed edge from the vertex
//    int firstDirectedEdge = directedEdges[vertexIndex];
//    int faceIndex = firstDirectedEdge / 3;
//    oneRing.push_back(faces[faceIndex][(firstDirectedEdge + 1) % 3]);
//    // get the other half of the directed edge
//    int otherHalfIndex = otherHalf[firstDirectedEdge];
//    // get the next edge in the face
//    int otherHalfFaceIndex = otherHalfIndex / 3;
//    int otherHalfVertexIndex = ((otherHalfIndex % 3) + 1) % 3;
//    otherHalfIndex = otherHalfFaceIndex * 3 + otherHalfVertexIndex;
//    while (otherHalfIndex != firstDirectedEdge)
//    {
//        // get the vertex index
//        int vertex = faces[otherHalfFaceIndex][(otherHalfVertexIndex + 1) % 3];
//        // -debug-
//        // if the vertex without + 1 isn't the vertex index then we have a problem
//        if (faces[otherHalfFaceIndex][otherHalfVertexIndex] != vertexIndex)
//        {
//            throw std::runtime_error("Vertex index doesn't match");
//        }
//        // add the vertex to the one ring
//        oneRing.push_back(vertex);
//        // get the next directed edge
//        otherHalfIndex = otherHalf[otherHalfIndex];
//        // get the next edge in the face
//        otherHalfFaceIndex = otherHalfIndex / 3;
//        otherHalfVertexIndex = ((otherHalfIndex % 3) + 1) % 3;
//        otherHalfIndex = otherHalfFaceIndex * 3 + otherHalfVertexIndex;
//    }
//    return oneRing;
}

void Simplification::generateCurvatures()
{
    // loop through each vertex and get the curvature
    for (int i = 0; i < vertices.size(); i++)
    {
        curvatures[i] = findCurvature(i);
    }
}

void Simplification::updateCurvatures()
{
    // remove the curvature for the removed vertex
    curvatures.erase(removedVertexIndex);
    // update the curvatures for the vertices that have been affected by the removal
    // these are the vertices in the one ring of the removed vertex - stored in holeEdges - precisely the start of each edge
    for (auto &edge: holeEdges)
    {
        curvatures[edge.start] = findCurvature(edge.start);
    }
}

Cartesian3 Simplification::getHoleNormal(std::vector<Edge> boundary)
{
    Cartesian3 normal = {0, 0, 0};
    // for each edge in the boundary
    for (int i = 0; i < boundary.size(); i++)
    {
        // get the two vectors
        Cartesian3 edge1 = vertices[boundary[i].start] - vertices[boundary[i].end];
        Cartesian3 edge2 = vertices[boundary[(i + 1) % boundary.size()].end] - vertices[boundary[(i + 1) % boundary.size()].start];
        // get the normal of the two vectors
        Cartesian3 faceNormal = edge1.cross(edge2);
        // add the normal to the total normal
        normal = normal + faceNormal.normalise();
    }
    // don't need to normalise since we only need the direction
    return normal.normalise();
}

int Simplification::getConcaveEdge(std::vector<Edge> boundary, Cartesian3 boundaryNormal)
{
    // for each edge in the boundary
    for (int i = 0; i < boundary.size(); i++)
    {
        // get the two edges
        Edge edge1 = boundary[i];
        Edge edge2 = boundary[(i + 1) % boundary.size()];
        // get the two vectors
        Cartesian3 vector1 = vertices[edge1.start] - vertices[edge1.end];
        Cartesian3 vector2 = vertices[edge2.end] - vertices[edge2.start];
        // get the cross product of the two vectors
        Cartesian3 crossProduct = vector1.cross(vector2).normalise();
        // get the dot product of the cross product and the boundary normal
        float dp = dotProduct(crossProduct, boundaryNormal);
        // if the dot product is positive then the angle between the two edges is convex as the normal they create is pointing in the same direction as the boundary normal
        if (dp > 0)
        {
            return i;
        }
    }
    return 0;
    // it isn't possible for no edge to be convex
    std::cerr << "No concave edge found on boundary" << std::endl;
    // print out the boundary
    for (auto &edge: boundary)
    {
        std::cerr << "Edge: " << vertices[edge.start] << " " << vertices[edge.end] << std::endl;
    }
    exit(-4);
}

void Simplification::cleanUpNonManifoldEdges()
{
    std::cout << "cleaning" << std::endl;
//    // delete any vertices that are not in any face
//    std::vector<int> verticesToDelete;
//    for (int i = 0; i < vertices.size(); i++)
//    {
//        bool found = false;
//        for (auto &face: faces)
//        {
//            for (int j = 0; j < 3; j++)
//            {
//                if (face[j] == i)
//                {
//                    found = true;
//                    break;
//                }
//            }
//            if (found)
//            {
//                break;
//            }
//        }
//        if (!found)
//        {
//            verticesToDelete.push_back(i);
//        }
//    }
//    // for each vertex in each face
//    for(auto &face:faces)
//    {
//        for (int i = 0; i < 3; i++)
//        {
//            // for each vertex to delete
//            for (auto &vertexToDelete: verticesToDelete)
//            {
//                if (face[i] > vertexToDelete)
//                {
//                    face[i]--;
//                }
//            }
//        }
//    }
//    std::cout << "Deleting " << verticesToDelete.size() << " loose vertices" << std::endl;
//    for (int i = verticesToDelete.size() - 1; i >= 0; i--)
//    {
//        vertices.erase(vertices.begin() + verticesToDelete[i]);
//    }
    // search through the mesh for faces that have a back face and delete both faces
    // so the loop stays fixed size we will store the faces to delete in a vector
    std::vector<int> facesToDelete;
    // for each face thats newly added
    for (int i = faces.size() - facesAdded - 1; i < faces.size(); i++)
    {
        // get the back face
        Face backFace = {faces[i][1], faces[i][0], faces[i][2]};
        // for each face once again
        for (int j = 0; j < faces.size() - facesAdded; j++)
        {
            // if the back face is the same as the face
            if (faces[j] == backFace)
            {
                // add the faces to delete
                facesToDelete.push_back(i);
                facesToDelete.push_back(j);
            }
        }
    }
    std::cout << "There are " << facesToDelete.size() / 2 << " faces to delete" << std::endl;
    // delete the faces
    for (int i = 0; i < facesToDelete.size(); i++)
    {
        // can't do it nicely with a swap in case one we want to delete is at the end
        faces.erase(faces.begin() + facesToDelete[i]);
        for (int j = i + 1; j < facesToDelete.size(); j++)
        {
            if (facesToDelete[j] > facesToDelete[i])
            {
                facesToDelete[j]--;
            }
        }
    }
}
bool Simplification::testNonManifoldEdges()
{
    // only the edges that have recently been effected could be non-manifold
    // these are the edges that are in the newly added faces
    // we should see these edges exactly twice in the faces vector (order independent)
    std::vector<Edge> affectedEdges;
    for(int i = faces.size() - facesAdded - 1; i < faces.size(); i++)
    {
        for (int j = 0; j < 3; j++)
        {
            Edge edge = {faces[i][j], faces[i][(j + 1) % 3]};
            affectedEdges.push_back(edge);
        }
    }
    std::vector<int> edgeCounts(affectedEdges.size(), 0);
    // for each face
    for (auto &face: faces)
    {
        // for each edge in the face
        for (int i = 0; i < 3; i++)
        {
            Edge edgeInFace = {face[i], face[(i + 1) % 3]};
            // for each edge in the hole
            for (int j = 0; j < affectedEdges.size(); j++)
            {
                if (edgeInFace == affectedEdges[j])
                {
                    edgeCounts[j]++;
                }
            }
        }
    }
    // if the counts aren't all 2 then we have a non-manifold edge
    if(std::ranges::all_of(edgeCounts, [](int count){return count == 2;}))
    {
        return true;
    }
    // if we have a non-manifold edge then we need to backtrack
    return false;
}




