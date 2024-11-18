//
// Created by thomas on 18/11/24.
//

#include <algorithm>
#include "Simplification.h"

void Simplification::simplifyMesh()
{
    // TODO: work out exit conditions
    //  for now, if there is no vertex that can be removed without breaking the eulerian condition, then we exit
    bool removed = true;
    while (removed)
    {
        removed = false;
        int currentVertex = 0;
        int smallestCurvature = findSmallestCurvature(0);
        while (!removed && currentVertex < vertices.size())
        {
            if (!removeVertex(smallestCurvature))
            {
                // if we can't remove the vertex then we need to backtrack
                backtrack();
                currentVertex++;
                smallestCurvature = findSmallestCurvature(currentVertex);
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
    holeEdges.clear();
    // remove the vertex from the mesh
    // remove the vertex from the vertices vector
    vertices.erase(vertices.begin() + vertexIndex);
    // remove any faces that contain the vertex, and add the edges of the hole to the holeEdges vector
    for (auto &face: faces)
    {
        if (face[0] == vertexIndex)
        {
            holeEdges.push_back({face[1], face[2]});
            facesAdded++;
            // remove the face
            face = faces.back();
            faces.pop_back();
        } else if (face[1] == vertexIndex)
        {
            holeEdges.push_back({face[0], face[2]});
            facesAdded++;
            // remove the face
            face = faces.back();
            faces.pop_back();
        } else if (face[2] == vertexIndex)
        {
            holeEdges.push_back({face[0], face[1]});
            facesAdded++;
            // remove the face
            face = faces.back();
            faces.pop_back();
        } else
        {
            // subtract the vertex index in the faces by 1 if it is greater than the removed vertex
            if (face[0] > vertexIndex)
            {
                face[0]--;
            }
            if (face[1] > vertexIndex)
            {
                face[1]--;
            }
            if (face[2] > vertexIndex)
            {
                face[2]--;
            }
        }
    }
    // triangulate the hole
    fillHole(holeEdges);

    // check that the mesh is still eulerian
    return isEulerian();
}

void Simplification::fillHole(std::vector<Edge>& boundary)
{
    // order the boundary edges so that they form a single cycle
    std::vector<Edge> orderedBoundary;
    orderedBoundary.push_back(boundary[0]);
    while(orderedBoundary.size() < boundary.size())
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
    // while the boundary has more than 3 edges
    while(boundary.size() > 3)
    {
        // get the two edges that form the smallest angle
        int smallestAngleIndex = findSmallestAngle(orderedBoundary);
        Edge smallestAngleEdge1 = orderedBoundary[smallestAngleIndex];
        Edge smallestAngleEdge2 = orderedBoundary[(smallestAngleIndex + 1) % orderedBoundary.size()];
        // add the triangle to the faces
        faces.push_back({smallestAngleEdge1.start, smallestAngleEdge1.end, smallestAngleEdge2.end});
    }
    // add the last triangle
    faces.push_back({orderedBoundary[0].start, orderedBoundary[0].end, orderedBoundary[1].end});

}

int Simplification::findSmallestAngle(const std::vector<Edge> &boundary)
{
    int smallestAngleIndex = 0;
    float smallestAngle = 200.f;
    for (int i = 0; i < boundary.size(); i++)
    {
        // get the two edges
        Edge edge1 = boundary[i];
        Edge edge2 = boundary[(i + 1) % boundary.size()];
        float angleBetween = getAngleBetweenEdges(edge1, edge2);
        if (angleBetween < smallestAngle)
        {
            smallestAngle = angleBetween;
            smallestAngleIndex = i;
        }
    }
    return smallestAngleIndex;
}

float Simplification::findCurvature(int vertexIndex)
{
    return 0;
}
