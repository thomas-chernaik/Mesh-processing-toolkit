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

// hashing struct for the edge
// I want edges to be equal if they have the same start and end, regardless of the order
// so for the hash we will add the two together
struct EdgeHash
{
    std::size_t operator()(const Edge &edge) const
    {
        return std::hash<int>()(edge.start + edge.end);
    }
};

#endif //MODELLINGCWK1_EDGE_H
