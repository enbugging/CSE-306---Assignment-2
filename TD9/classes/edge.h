#ifndef EDGE_H
#define EDGE_H

#include <algorithm>

class Edge {
public:
    Edge() {}

    Edge(int a, int b) : a(a), b(b), min_ab(std::min(a, b)), max_ab(std::max(a, b)) {}

    bool operator<(const Edge& other) const {
        return min_ab == other.min_ab ? max_ab < other.max_ab : min_ab < other.min_ab;
    }
    
    int a, b, min_ab, max_ab;
};

#endif