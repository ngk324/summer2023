#include "Graph.hpp"
#include <iostream>
#include <memory>

void Graph::constructSimpleGraph(const int x, const int y)
{

    // Randomly assign some state value and generate nodes
    srand(time(NULL));
    int idx{0};
    for (int i{0}; i < x; i++)
    {
        for (int j{0}; j < y; j++)
        {
            double rand_state{((double)rand() / (double)RAND_MAX)};
            std::shared_ptr<Node> ptr{new Node(i, j, idx, rand_state)};
            nodes.push_back(ptr);
            idx++;
        }
    }
    // Establish connections between nodes based on distance
    for (int i{0}; i < nodes.size(); i++)
    {
        for (int j{i + 1}; j < nodes.size(); j++)
        {
            if (nodes[i]->isNear(*nodes[j]))
            {
                nodes[i]->neighbors.push_back(nodes[j]);
                nodes[j]->neighbors.push_back(nodes[i]);
            }
        }
    }
}

std::shared_ptr<Graph> applyGradient(const Eigen::MatrixXf &newAdjacencyMatrix) const{
    auto newGraph = std::make_shared<Graph>();
    newGraph->nodes = nodes;
    newGraph->adjacencyMatrix = newAdjacencyMatrix;
    newGraph->connectivityMatrix = connectivityMatrix;
    int nNodes = nodes.size();
    Eigen::MatrixXf D(nNodes, nNodes);
    for (int i{0}; i < nNodes; i++)
    {
        double sum = 0;
        for (int j{0}; j < nNodes; j++)
        {
            D(i, j) = 0;
            sum += (double)newAdjacencyMatrix(i, j);
        }
        D(i, i) = sum;
    }
    newGraph->degreeMatrix = D;
    Eigen::MatrixXf L(nNodes, nNodes);
    newGraph->laplacianMatrix = newGraph->degreeMatrix - newGraph->adjacencyMatrix;
}

void Graph::computeMatrices()
{
    int nNodes = nodes.size();
    Eigen::MatrixXf A(nNodes, nNodes);
    Eigen::MatrixXf C(nNodes, nNodes);
    for (int i{0}; i < nNodes; i++)
    {
        for (int j{0}; j < nNodes; j++)
        {
            A(i, j) = 0;
            if (nodes[i]->isNeighbor(nodes[j]))
            {
                A(i, j) = 1;
                A(j, i) = 1;
                C(i, j) = 1;
                C(j, i) = 1;
            }
        }
    }
    adjacencyMatrix = A;
    connectivityMatrix = C;
    Eigen::MatrixXf D(nNodes, nNodes);
    for (int i{0}; i < nNodes; i++)
    {
        double sum = 0;
        for (int j{0}; j < nNodes; j++)
        {
            D(i, j) = 0;
            sum += (double)A(i, j);
        }
        D(i, i) = sum;
    }
    degreeMatrix = D;
    Eigen::MatrixXf L(nNodes, nNodes);
    L = degreeMatrix - adjacencyMatrix;
    laplacianMatrix = L;
}

void eigenDecompose()
{
    Eigen::MatrixXf matrixToSolve = laplacianMatrix + eps * Eigen::MatrixXf::Identity(laplacianMatrix.size());
    Eigen::EigenSolver<Eigen::MatrixXf> solver(matrixToSolve);
    eigenValues = solver.eigenvalues().real();
    eigenVectors = solver.eigenvectors().real();

    // Combine eigenvalues and eigenvectors into a std::vector of pairs
    std::vector<std::pair<double, Eigen::VectorXd>> eigenPairs;
    for (int i = 0; i < eigenValues.size(); ++i) {
        eigenPairs.push_back(std::make_pair(eigenValues[i], eigenVectors.col(i)));
    }

    // Sort the vector of pairs based on eigenvalues in ascending order
    std::sort(eigenPairs.begin(), eigenPairs.end(), [](const std::pair<double, Eigen::VectorXd>& a,
                const std::pair<double, Eigen::VectorXd>& b) {
            return a.first < b.first;
            });

    // Extract the sorted eigenvalues and eigenvectors
    for (int i = 0; i < eigenValues.size(); ++i) {
        eigenValues[i] = eigenPairs[i].first;
        eigenVectors.col(i) = eigenPairs[i].second;
    }
}


