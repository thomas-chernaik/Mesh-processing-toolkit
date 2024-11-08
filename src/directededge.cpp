//
// Created by thomas on 07/11/24.
//

#include <sstream>
#include "directededge.h"

void DirectedEdge::readFile(std::string filename)
{
    // clear the faces and vertices vectors in case they already contain a mesh
    faces.clear();
    vertices.clear();
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
    Vertex currentVertex;
    int currentIndex;
    // we are reading vertices while the line's first word is Vertex
    while (line.substr(0, 6) == "Vertex")
    {
        // create a new vertex
        Vertex vertex;
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
    // we are reading faces while the line's first word is Face
    while (line.substr(0, 4) == "Face")
    {
        // create a new face
        Face face;
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
    // check we don't have anything else left but white space
    char whiteSpaceCheckerChar;
    while (fileReader.get(whiteSpaceCheckerChar))
    {
        if (!isspace(whiteSpaceCheckerChar))
        {
            std::cerr << "File contains more than just faces and vertices" << std::endl;
            exit(-3);
        }
    }
    // close the file
    fileReader.close();

    // construct the directed edges
    constructDirectedEdges();
}

void DirectedEdge::constructDirectedEdges()
{
    // clear the directed edges and other half vectors
    directedEdges.clear();
    otherHalf.clear();
    // we will have numVertices directed edges
    directedEdges.resize(vertices.size(), -1);
    // we will have numFaces * 3 other half edges
    otherHalf.resize(faces.size() * 3, -1);
    // TODO: optimise - this is currently O(v * e) where v is the number of vertices and e is the number of edges
    // for each directed edge
    for (int i = 0; i < directedEdges.size(); i++)
    {
        // find the first face that contains this vertex
        for (int f = 0; f < faces.size(); f++)
        {
            // for each vertex in the face
            for (int v = 0; v < 3; v++)
            {
                // if the vertex is the one we are looking for
                if (faces[f][v] == i)
                {
                    // set the directed edge to the edge index for faces[f][v] -> faces[f][(v + 1) % 3]
                    // this is faces * 3 + v
                    // only set the directed edge if it hasn't been set yet, we want the first one (for funsies, technically we could use any)
                    if (directedEdges[i] == -1)
                        directedEdges[i] = f * 3 + v;

                }
            }
        }
        // if we didn't find a directed edge then we don't have a manifold mesh
        if (directedEdges[i] == -1)
        {
            std::cerr << "Mesh is not manifold as vertex " << i << " is not in any face" << std::endl;
            exit(-4);
        }
    }

    // for each edge
    for (int e = 0; e < otherHalf.size(); e++)
    {
        // find the face that contains this edge
        int f = e / 3;
        // find the vertex that is the start of this edge
        int vStart = faces[f][e % 3];
        // find the vertex that is the end of this edge
        int vEnd = faces[f][(e + 1) % 3];
        // TODO: optimise so that it adds both pairs of edges at the same time
        // search through all the edges again to find the other half
        for (int e2 = 0; e2 < otherHalf.size(); e2++)
        {
            // if the edge is the other half
            if (e2 != e && faces[e2 / 3][e2 % 3] == vEnd && faces[e2 / 3][(e2 + 1) % 3] == vStart)
            {
                // if there is already an other half then we don't have a manifold mesh
                if (otherHalf[e] != -1)
                {
                    std::cerr << "Mesh is not manifold as edge " << e << " has multiple other halves" << std::endl;
                    exit(-4);
                }
                // set the other half to be the edge index
                otherHalf[e] = e2;
            }
        }
        // if we didn't find an other half then we don't have a manifold mesh
        if (otherHalf[e] == -1)
        {
            std::cerr << "Mesh is not manifold as edge " << e << " has no other half" << std::endl;
            exit(-4);
        }
    }
}

void DirectedEdge::writeFile(std::string filename)
{
    // get the object name to put in the header data, and give easier to read outputs
    std::string objectName = filename;
    // strip data before the last / or \ in the filename
    size_t lastSlash = objectName.find_last_of(std::string{'\\', '/'});
    // if there isn't a slash at all lastSlash will be npos (-1)
    if (lastSlash != std::string::npos)
    {
        objectName = objectName.substr(lastSlash + 1);
    }
    // open the file to write
    std::ofstream fileWriter(filename + ".diredge");
    if (!fileWriter.is_open())
    {
        std::cerr << "Failed to write to " << filename << ".face" << std::endl;
        exit(-6);
    }

    std::cout << "Writing the " << objectName << " to a .diredge file" << std::endl;
    // add the header comments in according to the specification
    fileWriter << "# University of Leeds 2024-2025" << "\n" << "# COMP 5893M Assignment 1" << "\n"
               << "# Thomas Chernaik" << "\n" << "# sc21trc" << "\n" << "#" << "\n" << "# Object Name: " << objectName
               << "\n" << "# Vertices=" << vertices.size() << " Faces=" << faces.size() << "\n" << "#" << "\n";

    // write the vertices in
    for (int vertexIndex = 0; vertexIndex < vertices.size(); vertexIndex++)
    {
        fileWriter << "Vertex " << vertexIndex << "  " << vertices[vertexIndex].x << "  " << vertices[vertexIndex].y
                   << "  " << vertices[vertexIndex].z << "\n";
    }
    // write the first directed edges in
    for (int vertexIndex = 0; vertexIndex < directedEdges.size(); vertexIndex++)
    {
        fileWriter << "FirstDirectedEdge " << vertexIndex << "  " << directedEdges[vertexIndex] << "\n";
    }
    // write the faces in
    for (int faceIndex = 0; faceIndex < faces.size(); faceIndex++)
    {
        fileWriter << "Face " << faceIndex << "  " << faces[faceIndex][0] << "  " << faces[faceIndex][1] << "  "
                   << faces[faceIndex][2] << "\n";
    }
    // write the other half edges in
    for (int edgeIndex = 0; edgeIndex < otherHalf.size(); edgeIndex++)
    {
        fileWriter << "OtherHalf " << edgeIndex << "  " << otherHalf[edgeIndex] << "\n";
    }

}