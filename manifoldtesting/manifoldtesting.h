//
// Created by thomas on 08/11/24.
//

#ifndef MODELLINGCWK1_MANIFOLDTESTING_H
#define MODELLINGCWK1_MANIFOLDTESTING_H

#include "../utilities/vertex.h"
#include "../utilities/face.h"
#include "../src/directededge.h"
#include "../src/faceindex.h"

class ManifoldTesting
{
public:
    // called to load in a mesh from a .tri file
    void LoadMesh(std::string filename);
    void CheckForPinchPoints();



private:
    Directed
};


#endif //MODELLINGCWK1_MANIFOLDTESTING_H
