#ifndef _GRAPH_H_
#define _GRAPH_H_

#include "Node.hpp"
#include <Eigen/Dense>

class Graph
{
private:
public:
    std::vector<std::shared_ptr<Node>> nodes;
    Eigen::MatrixXf adjacencyMatrix, degreeMatrix, laplacianMatrix;

    void constructSimpleGraph(int x, int y);

    void computeMatrices();

    void simulateDynamics(int tMax);
};

#endif