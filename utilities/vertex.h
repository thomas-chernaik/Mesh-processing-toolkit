//
// Created by thomas on 07/11/24.
//

#ifndef FACEINDEX2DIRECTEDEDGE_VERTEX_H
#define FACEINDEX2DIRECTEDEDGE_VERTEX_H

#include <limits>
#include <cmath>
// stores the data in a vertex (the 3D coordinate it exists at)

struct Vertex
{
    float x;
    float y;
    float z;

    bool operator==(const Vertex &other) const
    {
        // checks that the 3 values in the vertex are equivalent
        float epsilon = std::numeric_limits<float>::epsilon();
        return std::abs(x - other.x) < epsilon && std::abs(y - other.y) < epsilon && std::abs(z - other.z) < epsilon;
    }
};
extern Vertex nullVertex;
#endif //FACEINDEX2DIRECTEDEDGE_VERTEX_H