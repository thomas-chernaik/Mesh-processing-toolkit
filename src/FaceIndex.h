#ifndef FACE2FACEI
#define FACE2FACEI

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <limits>
#include <algorithm>
#include "vertex.h"
#include "face.h"

class FaceIndex
{
public:

    // called to read a .tri file into the instance
    void readFile(const std::string& filename);

    // called to write the mesh stored in the instance to a .face file
    void writeFile(const std::string& filename);


private:


    // gets the index of the vertex with the same value as the parameter from vertices, appending the vertex to the list if needed.
    int getVertexIndex(Vertex vertex);

    std::vector<Vertex> vertices;
    std::vector<Face> faces;

};


#endif