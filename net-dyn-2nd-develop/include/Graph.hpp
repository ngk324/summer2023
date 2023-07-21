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

    void constructSimpleGraph(const int x, const int y);

    void computeMatrices();

    void computeMatrices2();

    void simulateDynamics(const int tMax);

    void calc_grad_descent();
};

#endif