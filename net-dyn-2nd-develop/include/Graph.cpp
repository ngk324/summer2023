#include "Graph.hpp"
#include <iostream>

void Graph::constructSimpleGraph(int x, int y)
{

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

void Graph::computeMatrices()
{
    int nNodes = nodes.size();
    Eigen::MatrixXf A(nNodes, nNodes);
    for (int i{0}; i < nNodes; i++)
    {
        for (int j{0}; j < nNodes; j++)
        {
            A(i, j) = 0;
            if (nodes[i]->isNeighbor(nodes[j]))
            {
                A(i, j) = 1;
                A(j, i) = 1;
            }
        }
    }
    adjacencyMatrix = A;
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

void Graph::computeMatrices2()
{
    int nNodes = nodes.size();
    Eigen::MatrixXf A(nNodes, nNodes);
    for (int i{0}; i < nNodes; i++)
    {
        for (int j{0}; j < nNodes; j++)
        {
            A(i, j) = 0;
            if (nodes[i]->isNeighbor(nodes[j]))
            {
                std::cout << i << ", " << j << "\n";
                A(i, j) = 1;
                A(j, i) = 1;
            }
        }
    }


    adjacencyMatrix = A;
    Eigen::MatrixXf D(nNodes, nNodes);
    for (int i{0}; i < nNodes; i++)
    {
        int sum = 0;
        for (int j{0}; j < nNodes; j++)
        {
            D(i, j) = 0;
            sum += A(i, j);
        }
        D(i, i) = sum;
    }
    degreeMatrix = D;
    Eigen::MatrixXf L(nNodes, nNodes);
    L = degreeMatrix - adjacencyMatrix;
    laplacianMatrix = L;
}