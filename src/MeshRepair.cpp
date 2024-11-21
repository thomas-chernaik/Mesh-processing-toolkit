//
// Created by thomas on 15/11/24.
//

#include "MeshRepair.h"

#include <cmath>

// find the first unpaired edge (pair is -1)
int MeshRepair::findLooseEdge()
{
    for (int i = 0; i < otherHalf.size(); i++)
    {
        if (otherHalf[i] == -1)
        {
            return i;
        }
    }
    return -1;
}

void MeshRepair::repairMesh()
{
    // calculate the different components
    testMultipleComponents();
    // if there is more than one component then remove all but the largest
    if (components.size() > 1)
    {
        removeAllButLargestComponent();
    }
    int unpairedEdge;
    // while there are still unpaired edges, patch them up
    while ((unpairedEdge = findLooseEdge()) != -1)
    {
        // walk around the edge
        auto boundary = walkAroundEdge(unpairedEdge);
        // see if there are any reflex angles in the boundary
        // fill the hole
        fillHole(boundary);

        // print out the number of unpaired edges
        int numUnpairedEdges = 0;
        for (int i: otherHalf)
        {
            if (i == -1)
            {
                numUnpairedEdges++;
            }
        }
        //std::cout << "Number of unpaired edges: " << numUnpairedEdges << std::endl;
        if (numUnpairedEdges == 1)
            break;
    }
    // compute the directed edges
    constructDirectedEdges();
}

void MeshRepair::removeAllButLargestComponent()
{
    // create new vertices and faces vectors
    std::vector<Vertex> newVertices;
    std::vector<Face> newFaces;
    // create a conversion dictionary from the oldVertexIndex to the newVertexIndex
    std::vector<int> conversion(vertices.size(), -1);
    // find the largest component
    int largestComponent = 0;
    for (int i = 1; i < components.size(); i++)
    {
        if (components[i].size() > components[largestComponent].size())
        {
            largestComponent = i;
        }
    }
    // for each face in the largest component
    for (int faceIndex: components[largestComponent])
    {
        // for each vertex in the face
        for (int i = 0; i < 3; i++)
        {
            // if we haven't already added the vertex to the newVertices vector
            if (conversion[faces[faceIndex][i]] == -1)
            {
                // add the vertex to the newVertices vector
                conversion[faces[faceIndex][i]] = (int) newVertices.size();
                newVertices.push_back(vertices[faces[faceIndex][i]]);
            }
        }
        // add the face to the newFaces vector
        newFaces.push_back({conversion[faces[faceIndex][0]], conversion[faces[faceIndex][1]], conversion[faces[faceIndex][2]]});
    }
    // set the new vertices and faces
    vertices = newVertices;
    faces = newFaces;
    int count = 0;
    for (auto i: otherHalf)
    { if (i == -1) count++; }
    //std::cout << "Number of unpaired edges: " << count << std::endl;
    // generate the directed edges and other half
    constructDirectedEdges();
    // calculate the different components
    testMultipleComponents();
    count = 0;
    for (auto i: otherHalf)
    { if (i == -1) count++; }
    //std::cout << "Number of unpaired edges: " << count << std::endl;
    // if there is still more than one component I did smth wrong with the code
    if (components.size() > 1)
    {
        std::cerr << "Failed to remove all but the largest component !!!! should not happen !!!!!" << std::endl;
        exit(-5);
    }
}

void MeshRepair::writeRepairedMeshTri(const std::string &filename)
{
    // open the file to write
    std::ofstream fileWriter(filename + ".tri");
    if (!fileWriter.is_open())
    {
        std::cerr << "Failed to write to " << filename << ".tri" << std::endl;
        exit(-6);
    }

    std::cout << "Writing the repaired mesh to a .tri file" << std::endl;
    // add the number of faces to the header
    fileWriter << faces.size() << std::endl;
    // add each face to the file
    for (auto face: faces)
    {
        fileWriter << vertices[face[0]].x << " " << vertices[face[0]].y << " " << vertices[face[0]].z << " ";
        fileWriter << vertices[face[1]].x << " " << vertices[face[1]].y << " " << vertices[face[1]].z << " ";
        fileWriter << vertices[face[2]].x << " " << vertices[face[2]].y << " " << vertices[face[2]].z << std::endl;
    }
    std::cout << "Successfully wrote the repaired mesh to " << filename << ".tri" << std::endl;
    // close the file
    fileWriter.close();
}


void MeshRepair::writeRepairedMeshFace(const std::string &filename)
{
    // open the file to write
    std::ofstream fileWriter(filename + ".face");
    if (!fileWriter.is_open())
    {
        std::cerr << "Failed to write to " << filename << ".face" << std::endl;
        exit(-6);
    }

    std::cout << "Writing the repaired mesh to a .face file" << std::endl;
    // add the header comments in according to the specification
    fileWriter << "# University of Leeds 2024-2025" << "\n" << "# COMP 5893M Assignment 1" << "\n"
               << "# Thomas Chernaik" << "\n" << "# sc21trc" << "\n" << "#" << "\n" << "# Object Name: " << filename
               << "\n" << "# Vertices=" << vertices.size() << " Faces=" << faces.size() << "\n" << "#" << "\n";
    // write the vertices in
    for (int vertexIndex = 0; vertexIndex < vertices.size(); vertexIndex++)
    {
        fileWriter << "Vertex " << vertexIndex << "  " << vertices[vertexIndex].x << "  " << vertices[vertexIndex].y
                   << "  " << vertices[vertexIndex].z << "\n";
    }
    // write the faces in
    for (int faceIndex = 0; faceIndex < faces.size(); faceIndex++)
    {
        fileWriter << "Face " << faceIndex << "  " << faces[faceIndex][0] << "  " << faces[faceIndex][1] << "  "
                   << faces[faceIndex][2] << "\n";
    }
    // close the file
    fileWriter.close();
    std::cout << "Successfully wrote the repaired mesh to " << filename << ".face" << std::endl;
}

void MeshRepair::fillHole(const std::vector<Edge> &boundary)
{
    int previousFaces = (int) faces.size();
    std::cout << "filling hole with " << boundary.size() << " edges" << std::endl;
//    // if the boundary is a triangle then we can just add the face
//    if (boundary.size() == 3)
//    {
//        faces.push_back({boundary[0].end, boundary[0].start, boundary[1].end});
//        // create the 3 new edges
//        Edge edge1 = {boundary[0].start, boundary[0].end};
//        Edge edge2 = {boundary[1].end, boundary[0].start};
//        Edge edge3 = {boundary[0].end, boundary[1].end};
//        otherHalf.resize(otherHalf.size() + 3, 1);
//        edge currentEdge = 0;
//        for (auto face: faces)
//        {
//            for (int i = 0; i < 3; i++)
//            {
//                if (face[i] == edge1.start && face[(i + 1) % 3] == edge1.end)
//                {
//                    otherHalf[currentEdge] = otherHalf.size() + 1;
//                    otherHalf[otherHalf.size() - 3 + 1] = currentEdge;
//                }
//                if (face[i] == edge2.start && face[(i + 1) % 3] == edge2.end)
//                {
//                    otherHalf[currentEdge] = otherHalf.size() + 2;
//                    otherHalf[otherHalf.size() - 3 + 2] = currentEdge;
//                }
//                if (face[i] == edge3.start && face[(i + 1) % 3] == edge3.end)
//                {
//                    otherHalf[currentEdge] = otherHalf.size() + 3;
//                    otherHalf[otherHalf.size() - 3 + 3] = currentEdge;
//                }
//                currentEdge++;
//            }
//        }
//
//        return;
//    }
    // get the set of vertices in the boundary
    std::vector<int> boundaryVertices;
    boundaryVertices.reserve(boundary.size());
    for (auto edge: boundary)
    {
        boundaryVertices.push_back(edge.start);
    }
    // get the centre of gravity of the boundary
    Vertex centreOfGravity = getCentreOfGravity(boundaryVertices);
    // add the centre of gravity to the vertices
    vertices.push_back(centreOfGravity);
    // get the index of the centre of gravity
    int centreOfGravityIndex = (int) vertices.size() - 1;
    // add a new face for each edge in the boundary
    for (auto edge: boundary)
    {
        // the edge that pairs to one on the boundary will be winding the other way
        faces.push_back({edge.end, edge.start, centreOfGravityIndex});
    }
    // add the other half edges
    int facesAdded = (int) faces.size() - previousFaces;
    addNewDirectedEdges(facesAdded);
}

std::vector<Edge> MeshRepair::walkAroundEdge(int edgeIndex)
{
    std::vector<Edge> boundary;
    // get the first edge
    Edge firstEdge = {faces[edgeIndex / 3][edgeIndex % 3], faces[edgeIndex / 3][(edgeIndex + 1) % 3]};
    boundary.push_back(firstEdge);
    // get the options for the next edge
    std::vector<Edge> options = getUnpairedEdges(firstEdge.end);
    if (options.empty())
    {
        std::cerr << "No edges to choose from" << std::endl;
        exit(-4);
    }
    // choose the option with the smallest angle
    Edge nextEdge = chooseSmallestAngledEdge(options, firstEdge);
    // while we haven't returned to the start
    while (!(nextEdge == firstEdge))
    {
        // if we have already visited this vertex then we have a cycle and want to return from the previously visited vertex to the end
        for (int i = 1; i < boundary.size(); i++)
        {
            if (boundary[i].start == nextEdge.end)
            {
                // remove the edges from the boundary
                boundary.erase(boundary.begin(), boundary.begin() + i - 1);
                boundary.push_back(nextEdge);
                return boundary;
            }
        }
        // add the edge to the boundary
        boundary.push_back(nextEdge);
        // get the options for the next edge
        options = getUnpairedEdges(nextEdge.end);
        // remove the edges already in the boundary
        removeItems(options, boundary);
        if (options.empty())
        {
            std::cerr << "No edges to choose from" << std::endl;
            exit(-4);
        }
        // choose the option with the smallest angle
        nextEdge = chooseSmallestAngledEdge(options, nextEdge);
    }
    return boundary;
}

std::vector<Edge> MeshRepair::getUnpairedEdges(int vertexIndex)
{
    std::vector<Edge> unpairedEdgesAtVertex;
    int edgeIndex = 0;
    // go through each face, and then edge in the face and see if the edge is unpaired
    // where the edge is unpaired, add it to the list of unpaired edges
    for (auto face: faces)
    {
        for (int i = 0; i < 3; i++)
        {
            if (face[i] == vertexIndex)
            {
                Edge edge = {face[i], face[(i + 1) % 3]};
                if (otherHalf[edgeIndex] == -1)
                {
                    unpairedEdgesAtVertex.push_back(edge);
                }
            }
            edgeIndex++;
        }
    }
    return unpairedEdgesAtVertex;
}

Edge MeshRepair::chooseSmallestAngledEdge(std::vector<Edge> edges, Edge currentEdge)
{
    // we should have at least one edge
    if (edges.empty())
    {
        std::cerr << "No edges to choose from" << std::endl;
        exit(-4);
    }
    // if we only have one edge then return it
    if (edges.size() == 1)
    {
        return edges[0];
    }
    //std::cout << "Choosing smallest angle from " << edges.size() << " edges" << std::endl;
    Edge smallestAngleEdge = edges[0];
    float smallestAngle = getAngleBetweenEdges(currentEdge, edges[0]);
    // loop through edges 2 - n and find the smallest angle
    for (int i = 1; i < edges.size(); i++)
    {
        float angle = getAngleBetweenEdges(currentEdge, edges[i]);
        if (angle < smallestAngle)
        {
            smallestAngle = angle;
            smallestAngleEdge = edges[i];
        }
    }
    return smallestAngleEdge;
}

float MeshRepair::getAngleBetweenEdges(Edge edge1, Edge edge2)
{
    // use the dot product to find the angle between two vectors without acos
    Cartesian3 edge1Vector = vertices[edge1.end] - vertices[edge1.start];
    Cartesian3 edge2Vector = vertices[edge2.start] - vertices[edge2.end];
    auto result = (float) -(dotProduct(edge1Vector, edge2Vector) / (edge1Vector.length() * edge2Vector.length()));
    // test if the angle is reflex, if it is we want to return the -result + 2
    if (isReflexAngle(edge1, edge2))
    {
        return -result + 2;
    }
    return result;
}

Vertex MeshRepair::getCentreOfGravity(const std::vector<int> &vertices)
{
    Vertex centreOfGravity{};
    // sum up the vertices locations
    for (int i: vertices)
    {
        centreOfGravity = centreOfGravity + this->vertices[i];
    }
    // divide by the number of vertices
    centreOfGravity = centreOfGravity / (float) vertices.size();
    return centreOfGravity;
}

void MeshRepair::removeItems(std::vector<Edge> &vector1, const std::vector<Edge> &vector2)
{
    // remove items in vector2 from vector1
    for (int j = 1; j < vector2.size(); j++)
    {
        for (int i = 0; i < vector1.size(); i++)
        {
            if (vector1[i] == vector2[j])
            {
                vector1.erase(vector1.begin() + i);
                break;
            }
        }
    }

}


bool MeshRepair::isReflexAngle(Edge edge1, Edge edge2)
{
//    std::cout << "reflexing ";
    // early exit if the edges belong to the same face
    // find the third point of the face that edge1 belongs to
    int thirdPoint = -1;
    for (auto & face : faces)
    {
        if ((face[0] == edge1.start && face[1] == edge1.end))
        {
            thirdPoint = face[2];
            break;
        } else if (face[1] == edge1.start && face[2] == edge1.end)
        {
            thirdPoint = face[0];
            break;
        } else if (face[2] == edge1.start && face[0] == edge1.end)
        {
            thirdPoint = face[1];
            break;
        }
    }
    if (thirdPoint == -1)
    {
        std::cerr << "Failed to find the third point of the face" << std::endl;
        exit(-4);
    }
    // if edge2 belongs to the same face then it is a reflex angle
    if (edge2.end == thirdPoint)
    {
        return true;
    }
    //
    // if edge2.end is to the right of edge1 then it is a reflex angle
    // what is left:
    // we have a forward vector - edge1.end - edge1.start
    // we have an up vector - the normal to the plane of the face of edge1
    // we can therefore compute a left vector - the cross product of the forward and up vectors
    // we can then compute the dot product of the left vector and the forward vector of edge2
    // if this is positive then edge2 is to the left of edge1 (convex angle)
    // if this is negative then edge2 is to the right of edge1 (reflex angle)
    // if this is zero then edge2 is colinear with edge1 (handle as reflex angle)
    Cartesian3 forward = vertices[edge1.end] - vertices[edge1.start];
    // normalise the forward vector
    forward = forward.normalise();
    Cartesian3 up = (vertices[thirdPoint] - vertices[edge1.start]).cross(forward);
    up = up.normalise();
    Cartesian3 left = forward.cross(up);
    left = left.normalise();
    Cartesian3 forward2 = vertices[edge2.end] - vertices[edge2.start];
    forward2 = forward2.normalise();
    float dp = dotProduct(left, forward2);
//    if (dp > 1e-5)
//    {
//        std::cout << "Reflex angle detected" << std::endl;
//        std::cout << "Dot product: " << dp << std::endl;
//        std::cout << "Triangle 1: " << vertices[edge1.start] << " " << vertices[edge1.end] << " "
//                  << vertices[thirdPoint]
//                  << std::endl;
//        std::cout << "Edge 2: " << vertices[edge2.start] << " " << vertices[edge2.end] << std::endl;
//        std::cout << "Vectors: forward: " << forward << " up: " << up << " left: " << left << " forward2: "
//                  << forward2 << std::endl;
//    }
    return dp > 1e-5;
}

void MeshRepair::addNewDirectedEdges(int numFacesAdded)
{
    // expand the other half vector
    otherHalf.resize(faces.size() * 3, 1);
    // add the new first directed edge - this will be the third edge of the first new face
    int firstDirectedEdge = (int) faces.size() - numFacesAdded;
    directedEdges.push_back(firstDirectedEdge * 3 + 2);
    // the first edge of each new face will be paired to a previously unpaired edge
    // add these first
    for (int i = (int)faces.size() - numFacesAdded; i < faces.size(); i++)
    {
        // get the first edge of the face
        Edge firstEdge = {faces[i][0], faces[i][1]};
        // find the edge that pairs to this edge - there is guaranteed to be one
        int otherHalfIndex;
        for (int j = 0; j < faces.size() - numFacesAdded; j++)
        {
            // get the three edges of the face
            Edge edge1 = {faces[j][0], faces[j][1]};
            Edge edge2 = {faces[j][1], faces[j][2]};
            Edge edge3 = {faces[j][2], faces[j][0]};
            if (firstEdge == edge1)
            {
                otherHalfIndex = j * 3;
                otherHalf[otherHalfIndex] = i * 3;
                break;
            } else if (firstEdge == edge2)
            {
                otherHalfIndex = j * 3 + 1;
                otherHalf[otherHalfIndex] = i * 3;
                break;
            } else if (firstEdge == edge3)
            {
                otherHalfIndex = j * 3 + 2;
                otherHalf[otherHalfIndex] = i * 3;
                break;
            }
        }

    }

}
