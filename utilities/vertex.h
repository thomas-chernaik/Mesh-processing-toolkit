//
// Created by thomas on 07/11/24.
//

#ifndef FACEINDEX2DIRECTEDEDGE_VERTEX_H
#define FACEINDEX2DIRECTEDEDGE_VERTEX_H

#include <limits>
#include <cmath>
#include <ostream>
#include "../triangle_renderer/Cartesian3.h"
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

    Cartesian3 operator-(const Vertex &other)
    {
        return {x - other.x, y - other.y, z - other.z};
    }
    Vertex operator+(const Vertex &other) const
    {
        return {x + other.x, y + other.y, z + other.z};
    }

    // operator to print out the vertex
    friend std::ostream &operator<<(std::ostream &os, const Vertex &vertex)
    {
        os << "(" << vertex.x << ", " << vertex.y << ", " << vertex.z << ")";
        return os;
    }

    Vertex operator/(float factor) const
    {
        return {x / factor, y / factor, z / factor};
    }

};
extern Vertex nullVertex;
#endif //FACEINDEX2DIRECTEDEDGE_VERTEX_H