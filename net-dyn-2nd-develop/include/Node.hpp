#ifndef _NODE_H_
#define _NODE_H_

#include <vector>
#include <memory>

class Node
{
private:
public:
    int x, y; // Coordinates (just for plotting)
    int id;
    std::vector<std::shared_ptr<Node>> neighbors;
    double z;     // State
    double z_old; // Previous state (for discrete time computations)
    double z_dot{0};
    double z_dot_old{0};

    Node(int xx, int yy, int id);
    Node(int xx, int yy, int id, double zz);

    double getDist(const Node &n) const;
    double isNear(const Node &n) const;

    bool isNeighbor(std::shared_ptr<Node> n) const;

    void print(const std::string head) const;
};

#endif