//
// Created by thomas on 18/11/24.
//

#include <algorithm>
#include <cmath>
#include "Simplification.h"

void Simplification::simplifyMesh()
{
    // TODO: work out exit conditions
    //  for now, if there is no vertex that can be removed without breaking the eulerian condition, then we exit
    bool removed = true;
    for (int i = 0; i < 1000; i++)
    {
        // -debug - write the mesh to a file with the iteration number
        // if its a multiple of 10
        if (i % 10 == 0)
        {
            //writeRepairedMeshTri("iteration" + std::to_string(i));
        }
        //writeRepairedMeshTri("iteration" + std::to_string(i));
        removed = false;
        int currentVertex = 0;
        int smallestCurvature = findSmallestCurvature(0);
        // print out the vertex with the smallest curvature
        while (!removed && currentVertex < vertices.size())
        {
            std::cout << "Smallest curvature: " << vertices[smallestCurvature] << std::endl;
            if (!removeVertex(smallestCurvature))
            {
                // debug - write the mesh to a file
                writeRepairedMeshTri("removed");
                std::cout << "backtracking" << std::endl;
                // if we can't remove the vertex then we need to backtrack
                backtrack();
                currentVertex++;
                smallestCurvature = findSmallestCurvature(currentVertex);
                // debug - write the mesh to a file
                writeRepairedMeshTri("de-removed");
                exit(0);
            } else
            {
                removed = true;
            }
        }
    }
}

int Simplification::findSmallestCurvature(int n)
{
    // find the nth smallest curvature
    std::vector<float> curvatures;
    std::vector<int> vertexIndices;
    for (int i = 0; i < vertices.size(); i++)
    {
        curvatures.push_back(findCurvature(i));
        vertexIndices.push_back(i);
    }
    // sort the vertices by curvature
    std::sort(vertexIndices.begin(), vertexIndices.end(), [&curvatures](int a, int b)
    {
        return curvatures[a] < curvatures[b];
    });
    return vertexIndices[n];
}

void Simplification::backtrack()
{
    // remove the last added faces
    faces.resize(faces.size() - facesAdded);
    // add the removed vertex back
    vertices.push_back(removedVertex);
    int vertexIndex = (int) vertices.size() - 1;
    // add the removed faces back
    for (auto &edge: holeEdges)
    {
        faces.push_back({edge.start, edge.end, vertexIndex});
    }

}

bool Simplification::removeVertex(int vertexIndex)
{
    // store the vertex that is being removed
    removedVertex = vertices[vertexIndex];
    holeEdges.clear();
    // remove the vertex from the mesh
    // remove the vertex from the vertices vector
    vertices.erase(vertices.begin() + vertexIndex);
    // remove any faces that contain the vertex, and add the edges of the hole to the holeEdges vector
    for (int i = 0; i < faces.size(); i++)
    {
        Face face = faces[i];
        if (face[0] == vertexIndex)
        {
            holeEdges.push_back({face[1], face[2]});
            // remove the face
            faces[i] = faces.back();
            faces.pop_back();
            i--;
        } else if (face[1] == vertexIndex)
        {
            holeEdges.push_back({face[2], face[0]});
            // remove the face
            faces[i] = faces.back();
            faces.pop_back();
            i--;
        } else if (face[2] == vertexIndex)
        {
            holeEdges.push_back({face[0], face[1]});
            // remove the face
            faces[i] = faces.back();
            faces.pop_back();
            i--;
        }
        else
        {
            // subtract the vertex index in the faces by 1 if it is greater than the removed vertex
            if (face[0] > vertexIndex)
            {
                faces[i][0]--;
            }
            if (face[1] > vertexIndex)
            {
                faces[i][1]--;
            }
            if (face[2] > vertexIndex)
            {
                faces[i][2]--;
            }
        }
    }
    // subtract the vertex indices in the holeEdges by 1 if they are greater than the removed vertex
    for (auto &edge: holeEdges)
    {
        if (edge.start > vertexIndex)
        {
            edge.start--;
        }
        if (edge.end > vertexIndex)
        {
            edge.end--;
        }

    }
    int previousFaces = (int) faces.size();
    std::cout << "filling hole with " << holeEdges.size() << " edges" << std::endl;
    // triangulate the hole
    triangulateHole(holeEdges);
    facesAdded = (int) faces.size() - previousFaces;

    // check that the mesh is still eulerian
    return isEulerian();
}


float Simplification::findCurvature(int vertexIndex)
{
    return findMeanCurvature(vertexIndex);
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
    return laplaceBeltrami.length() / 2;
}

float Simplification::findGaussianCurvature(int vertexIndex)
{
    // this is the "angular deficit"
    // its equal to 1 / Ai (what is Ai?) * (2pi - sum of angles around the vertex)
    // get the one ring vertices
    std::vector<int> oneRingVertices = getOneRingVertices(vertexIndex);
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

    // return the gaussian curvature
    return (float) (2.f * M_PI - sumOfAngles);
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
    // while the boundary has more than 2 edges
    while (orderedBoundary.size() > 3)
    {
        // get the two edges that form the smallest angle
        int smallestAngleIndex = findSmallestAngle(orderedBoundary);
        Edge smallestAngleEdge1 = orderedBoundary[smallestAngleIndex];
        Edge smallestAngleEdge2 = orderedBoundary[(smallestAngleIndex + 1) % orderedBoundary.size()];
        // add the triangle to the faces
        faces.push_back({smallestAngleEdge1.start, smallestAngleEdge1.end, smallestAngleEdge2.end});
        // remove the edges from the boundary and add the new edge
        // we need to maintain the order of the boundary, so we remove one edge and swap the other with the new edge
        // swap the new edge in with the second edge
        orderedBoundary[(smallestAngleIndex + 1) % orderedBoundary.size()] = {smallestAngleEdge1.start, smallestAngleEdge2.end};
        // remove the first edge
        orderedBoundary.erase(orderedBoundary.begin() + smallestAngleIndex);
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
    return std::acos(dotProduct(vector1, vector2) / (vector1.length() * vector2.length()));
}

bool Simplification::isEulerian()
{
    return true;
}

std::vector<int> Simplification::getOneRingVertices(int vertexIndex)
{
    // get the one ring edges
    std::vector<Edge> oneRing = getOneRing(vertexIndex);
    // get the one ring vertices (the start of each edge)
    std::vector<int> oneRingVertices;
    for (auto &edge: oneRing)
    {
        oneRingVertices.push_back(edge.start);
    }
    return oneRingVertices;
}




