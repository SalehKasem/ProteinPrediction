#pragma once

#include "Node.h"
#include <vector>

class NodePCN : public Node {
private:
    std::vector<Node> neighbors;
public:
#pragma region Getters and setters

    std::vector<Node> getNeighbors();

    void setNeighbors(std::vector<Node> neighbors);

#pragma endregion Getters and setters

    NodePCN();

    NodePCN(long protein, int index);

    NodePCN(Node node);

    ~NodePCN();
};

