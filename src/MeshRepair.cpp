//
// Created by thomas on 15/11/24.
//

#include "MeshRepair.h"

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
            // if the face contains this vertex
            if (face[0] == i || face[1] == i || face[2] == i)
            {
                // add the face to the new faces
                newFaces.push_back({conversion[face[0]], conversion[face[1]], conversion[face[2]]});
            }
        }
    }
    // set the new vertices and faces
    vertices = newVertices;
    faces = newFaces;
    // generate the directed edges and other half
    constructDirectedEdges();
    // calculate the different components
    testMultipleComponents();
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