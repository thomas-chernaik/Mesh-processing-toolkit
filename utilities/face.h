//
// Created by thomas on 07/11/24.
//

#ifndef FACE2FACEINDEX_FACE_H
#define FACE2FACEINDEX_FACE_H

#include <iostream>
// stores the data in a face (the indices of the three vertices that make it up)
struct Face
{
    int vertex1;
    int vertex2;
    int vertex3;

    int& operator[](size_t index)
    {
        switch(index)
        {
            case 0:
                return vertex1;
            case 1:
                return vertex2;
            case 2:
                return vertex3;
            default:
                std::cerr << "You're indexing a face wrong." << std::endl;
                exit(-4);
        }
    }
    bool operator==(const Face &other) const
    {
        return vertex1 == other.vertex1 && vertex2 == other.vertex2 && vertex3 == other.vertex3;
    }
};

extern Face nullFace;
#endif //FACE2FACEINDEX_FACE_H
