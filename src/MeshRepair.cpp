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
        std::cout << "Number of unpaired edges: " << numUnpairedEdges << std::endl;
    }
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
    // for each vertex in the largest component
    for (int i: components[largestComponent])
    {
        // add the vertex to the new vertices
        newVertices.push_back(vertices[i]);
        // set the conversion from the old index to the new index
        conversion[i] = (int)newVertices.size() - 1;
    }
    // for each vertex in the largest component
    for (int i: components[largestComponent])
    {
        // for each face
        for (auto &face: faces)
        {
            // if the face contains this vertex, and this is the lowest index vertex in the face
            if ((face[0] == i || face[1] == i || face[2] == i) && (face[0] >= i && face[1] >= i && face[2] >= i))
            {
                // add the face to the new faces
                newFaces.push_back({conversion[face[0]], conversion[face[1]], conversion[face[2]]});
            }
        }
    }
    // set the new vertices and faces
    vertices = newVertices;
    faces = newFaces;
    int count = 0;for(auto i:otherHalf){if(i == -1) count++;}std::cout << "Number of unpaired edges: " << count << std::endl;
    // generate the directed edges and other half
    constructDirectedEdges();
    // calculate the different components
    testMultipleComponents();
    count = 0;for(auto i:otherHalf){if(i == -1) count++;}std::cout << "Number of unpaired edges: " << count << std::endl;
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
    for (auto face : faces)
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

void MeshRepair::fillHole(std::vector<Edge> boundary)
{
    std::cout << "filling hole with " << boundary.size() << " edges" << std::endl;
    // get the set of vertices in the boundary
    std::vector<int> boundaryVertices;
    for (auto edge: boundary)
    {
        boundaryVertices.push_back(edge.start);
    }
    // get the centre of gravity of the boundary
    Vertex centreOfGravity = getCentreOfGravity(boundaryVertices);
    // add the centre of gravity to the vertices
    vertices.push_back(centreOfGravity);
    // get the index of the centre of gravity
    int centreOfGravityIndex = (int)vertices.size() - 1;
    // add a new face for each edge in the boundary
    for (auto edge: boundary)
    {
        // the edge that pairs to one on the boundary will be winding the other way
        faces.push_back({edge.end, edge.start, centreOfGravityIndex});
    }
    // generate the directed edges and other half
    constructDirectedEdges();
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
    while(!(nextEdge == firstEdge))
    {
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
    for(auto face:faces)
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
    if (edges.empty())
    {
        std::cerr << "No edges to choose from" << std::endl;
        exit(-4);
    }
    Edge smallestAngleEdge = edges[0];
    float smallestAngle = getAngleBetweenEdges(currentEdge, edges[0]);
    // loop through edges 2 - n
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
    Cartesian3 edge1Vector = vertices[edge1.end] - vertices[edge1.start];
    Cartesian3 edge2Vector = vertices[edge2.end] - vertices[edge2.start];
    return (float) (dotProduct(edge1Vector, edge2Vector) / (edge1Vector.length() * edge2Vector.length()));
}

Vertex MeshRepair::getCentreOfGravity(const std::vector<int> &vertices)
{
    Vertex centreOfGravity{};
    for (int i: vertices)
    {
        centreOfGravity = centreOfGravity + this->vertices[i];
    }
    centreOfGravity = centreOfGravity / (int)vertices.size();
    return centreOfGravity;
}

void MeshRepair::removeItems(std::vector<Edge> &vector1, const std::vector<Edge> &vector2)
{
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
