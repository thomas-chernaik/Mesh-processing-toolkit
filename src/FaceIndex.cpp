#include "FaceIndex.h"


int FaceIndex::getVertexIndex(Vertex vertex)
{
    int index = -1;
    // TODO: optimise this to use a hash map or similar, alternately parallelise this
    // do a linear search through the vertices and see if any match
    for (int i = 0; i < vertices.size(); i++)
    {
        if (vertices[i] == vertex)
        {
            index = i;
            break;
        }
    }
    // if the index is still -1 then this vertex hasn't been found yet
    // so we add it to the vertex list as a new vertex
    if (index == -1)
    {
        vertices.push_back(vertex);
        index = vertices.size() - 1;
    }
    return index;
}

void FaceIndex::readFile(std::string filename)
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
    int numFaces;
    fileReader >> numFaces;
    std::cout << "Reading in " << numFaces << " faces from " << filename << std::endl;
    faces.resize(numFaces);
    Vertex currentVertex;
    // read in each face
    for (int faceIndex = 0; faceIndex < numFaces; faceIndex++)
    {
        // for each vertex
        for (int vertexIndex = 0; vertexIndex < 3; vertexIndex++)
        {
            // read in the 3 points
            fileReader >> currentVertex.x;
            fileReader >> currentVertex.y;
            fileReader >> currentVertex.z;
            // get the vertex index for this vertex and assign it to the face
            faces[faceIndex][vertexIndex] = getVertexIndex(currentVertex);
        }
        // check that the face has 3 unique vertices
        if (faces[faceIndex][0] == faces[faceIndex][1] || faces[faceIndex][0] == faces[faceIndex][2] ||
            faces[faceIndex][1] == faces[faceIndex][2])
        {
            // TODO: decide whether this should be an error or just a warning
            std::cerr << "Face " << faceIndex << " has repeated vertices and is not a valid triangle" << std::endl;
            // exit(-3);
            // we can discard the face as it isn't a triangle
            numFaces--;
            faces.pop_back();
            faceIndex--;
        }
    }
    // check we don't have anything else left but white space
    char whiteSpaceCheckerChar;
    while (fileReader.get(whiteSpaceCheckerChar))
    {
        if (!std::isspace(whiteSpaceCheckerChar))
        {
            std::cerr << "Finished reading the file but there are still non-whitespace characters left" << std::endl;
            exit(-5);
        }
    }
    std::cout << "Successfully read in " << numFaces << " faces from " << filename << std::endl;
}

void FaceIndex::writeFile(std::string filename)
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
    std::ofstream fileWriter(filename + ".face");
    if (!fileWriter.is_open())
    {
        std::cerr << "Failed to write to " << filename << ".face" << std::endl;
        exit(-6);
    }

    std::cout << "Writing the " << objectName << " to a .face file" << std::endl;
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

    // write the faces in
    for (int faceIndex = 0; faceIndex < faces.size(); faceIndex++)
    {
        fileWriter << "Face " << faceIndex << "  " << faces[faceIndex][0] << "  " << faces[faceIndex][1] << "  "
                   << faces[faceIndex][2] << "\n";
    }

    fileWriter.close();
    std::cout << "succesfully written the " << objectName << " to a .face file" << std::endl;
}

