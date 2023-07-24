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
    Eigen::MatrixXi connectivityMatrix;
    Eigen::VectorXf eigenValues;
    Eigen::MatrixXf eigenVectors;

    void constructSimpleGraph(const int x, const int y);

    std::shared_ptr<Graph> applyGradient(const Eigen::MatrixXf &newAdjacencyMatrix) const;

    void computeMatrices();

    void simulateDynamics(const int tMax);
};

#endif
