//
// Created by thomas on 08/11/24.
//

#ifndef MODELLINGCWK1_EDGE_H
#define MODELLINGCWK1_EDGE_H

struct Edge
{
    int start;
    int end;

    bool operator==(const Edge &other) const
    {
        return (start == other.start && end == other.end) || (start == other.end && end == other.start);
    }
    bool exactlyEqual(const Edge &other) const
    {
        return start == other.start && end == other.end;
    }
};

#endif //MODELLINGCWK1_EDGE_H
