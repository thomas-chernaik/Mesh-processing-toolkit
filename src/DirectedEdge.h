//
// Created by thomas on 07/11/24.
//

#ifndef MODELLINGCWK1_DIRECTEDEDGE_H
#define MODELLINGCWK1_DIRECTEDEDGE_H
using edge = int;

#include "vertex.h"
#include "face.h"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "progressBar.h"



class DirectedEdge
{
public:

    // constructor with an optional argument to set the readWithErrors flag
    explicit DirectedEdge(bool readWithErrors = false) : readWithErrors(readWithErrors) {}

    // called to read a .face file into the instance
    void readFile(const std::string& filename);

    // called to write the mesh stored in the instance to a .directededge file
    void writeFile(const std::string& filename);

protected:

    // function used as debug to write the file to an obj file
    void writeObjFile(const std::string& filename);

    // called to compute the directed edges from the indexed faces
    void constructDirectedEdges();

    bool readWithErrors = false;

    void error(const std::string& message, int exitCode = -1) const
    {
        if (!readWithErrors)
        {
            std::cerr << message << std::endl;
            exit(exitCode);
        }
    }

    std::vector<Vertex> vertices;
    std::vector<Face> faces;
    std::vector<edge> directedEdges;
    std::vector<edge> otherHalf;
};


#endif //MODELLINGCWK1_DIRECTEDEDGE_H
