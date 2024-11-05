#ifndef FACE2FACEI
#define FACE2FACEI

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <limits>
#include <algorithm>

class Face2FaceIndex
{
    public:

    // called to read a .tri file into the instance
    void readFile(std::string filename);
    // called to write the mesh stored in the instance to a .face file
    void writeFile(std::string filename);

    

    private:
    // stores the data in a face (the indices of the three vertices that make it up)
    struct Face
    {
        int vertex1;
        int vertex2;
        int vertex3;
        int& operator[](size_t index);
    };
    // stores the data in a vertex (the 3D coordinate it exists at)
    struct Vertex
    {
        float x;
        float y;
        float z;
        bool operator==(const Vertex& other);
    };

    // gets the index of the vertex with the same value as the parameter from vertices, appending the vertex to the list if needed.
    int getVertexIndex(Vertex vertex);

    std::vector<Vertex> vertices;
    std::vector<Face> faces;

};



#endif