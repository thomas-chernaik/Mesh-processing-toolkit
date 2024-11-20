//
// Created by thomas on 20/11/24.
//

#ifndef MODELLINGCWK1_TRIANGLEINTERSECTIONRENDERING_H
#define MODELLINGCWK1_TRIANGLEINTERSECTIONRENDERING_H

enum Axis
{
    X,
    Y,
    Z
};

class TriangleIntersectionRendering
{
public:
    // function to initialise the openGL context and create the window
    void initOpenGL(int width, int height);

    // function to test if two triangles intersect
    bool trianglesIntersect(float *triangle1, float *triangle2);

private:
    // function to generate an orthographic projection matrix, taking the axis it looks down and the size of the view
    void generateOrthographicProjectionMatrix(Axis axis, float width, float height);

    // the screen width and height
    int screenWidth, screenHeight;

};


#endif //MODELLINGCWK1_TRIANGLEINTERSECTIONRENDERING_H
