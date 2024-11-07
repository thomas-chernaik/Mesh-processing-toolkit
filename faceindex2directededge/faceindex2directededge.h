//
// Created by thomas on 07/11/24.
//

#ifndef MODELLINGCWK1_FACEINDEX2DIRECTEDEDGE_H
#define MODELLINGCWK1_FACEINDEX2DIRECTEDEDGE_H
using edge = int;

#include "../utilities/vertex.h"
#include "../utilities/face.h"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>



class FaceIndex2DirectedEdge
{
public:
    // called to read a .face file into the instance
    void readFile(std::string filename);

    // called to write the mesh stored in the instance to a .directededge file
    void writeFile(std::string filename);

private:

    // called to compute the directed edges from the indexed faces
    void constructDirectedEdges();

    std::vector<Vertex> vertices;
    std::vector<Face> faces;
    std::vector<edge> directedEdges;
    std::vector<edge> otherHalf;
};


#endif //MODELLINGCWK1_FACEINDEX2DIRECTEDEDGE_H
