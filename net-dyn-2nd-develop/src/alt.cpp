#include "include/Graph.hpp"
#include <algorithm>
#include <memory>
#include <vector>
#include <Eigen/Eigenvalues>

struct MatrixIdx{
    int j;
    int k;
    MatrixIdx(int jj, int kk): j{jj}, k{kk}{}
};

// =====================================================================================
MatrixIdx findObjMaximizingEigIndex(const std::shared_ptr<Graph> currentGraph, const double h);
std::shared_ptr<Graph> graphGradientDescent(const std::shared_ptr<Graph> currentGraph);
double evalObjectiveFunctionAtJK(const double h, const double lambda_j, const double lambda_k);
std::vector<MatrixIdx> getNegativeMatrixIndices(const Eigen::MatrixXf &adjMatrix);
double evaluateObjGradientAtJK(const double h, const double lambda_j, const double lambda_k);
Eigen::MatrixXf computeObjGradient(const std::shared_ptr<Graph> currentGraph, const std::vector<MatrixIdx> &weightsToAvoid);
// =====================================================================================
double evalObjectiveFunctionAtJK(const double h, const double lambda_j, const double lambda_k){
    return ((4*h*h)/pow(pow(sqrt(lambda_j)-sqrt(lambda_k),2)+h*h,2));
}

MatrixIdx findObjMaximizingEigIndex(const std::shared_ptr<Graph> currentGraph, const double h){

    int max_j{0};
    int max_k{0};
    double maxVal{0.0};

    for(int j{0}; j<currentGraph->eigenValues.size(); j++){
        for(int k{j+1}; k<currentGraph->eigenValues.size(); k++){
            double currentVal = evalObjectiveFunctionAtJK(h,currentGraph->eigenValues[j],currentGraph->eigenValues[k]);
            if(currentVal>maxVal){
                max_j = j;
                max_k = k;
                maxVal = currentVal;
            }
        }
    }
    return MatrixIdx(max_j,max_k);
}

std::vector<MatrixIdx> getNegativeMatrixIndices(const Eigen::MatrixXf &adjMatrix){

    std::vector<MatrixIdx> res;
    for(int j{0}; j<adjMatrix.rows(); j++){
        for(int k{j+1}; k<adjMatrix.cols(); k++){
            if(adjMatrix.coeff(j,k)<0)
                res.push_back(MatrixIdx(j,k));
        }
    }
    return res;
}

double evaluateObjGradientAtJK(const double h, const double lambda_j, const double lambda_k){
    
    return 8*h*h*pow(sqrt(lambda_j)-sqrt(lambda_k),2)/(sqrt(lambda_j*lambda_k)*pow(pow(sqrt(lambda_j)-sqrt(lambda_k),2)+h*h,3));
}

Eigen::MatrixXf computeObjGradient(const std::shared_ptr<Graph> currentGraph, const std::vector<MatrixIdx> &weightsToAvoid){

    // if the objective function contains a max function find the necessary indices for the eigenvalues that maximize
    double h{1};
    auto maxEigIdx = findObjMaximizingEigIndex(currentGraph,h);
    // Compute the gradient using the derivation
    Eigen::MatrixXf gradient = Eigen::MatrixXf::Zero(currentGraph->adjacencyMatrix.size());
    double lambda_j = currentGraph->eigenValues[maxEigIdx.j];
    double lambda_k = currentGraph->eigenValues[maxEigIdx.k];
    gradient.coeff(maxEigIdx.j,maxEigIdx.k) = evaluateObjGradientAtJK(h,lambda_j,lambda_k);
    gradient.coeff(maxEigIdx.k,maxEigIdx.j) = evaluateObjGradientAtJK(h,lambda_j,lambda_k);

    // Remove the components corresponding to non existing edges
    gradient = gradient.array() * currentGraph->connectivityMatrix.array(); 

    // Remove the components corresponding to weightsToAvoid
    for(auto &el: weightsToAvoid){
        gradient.coeff(el.j,el.k) = 0;
        gradient.coeff(el.k,el.j) = 0;
    }

    return gradient;
}

std::shared_ptr<Graph> graphGradientDescent(const  std::shared_ptr<Graph> currentGraph){

    currentGraph.eigenDecompose();
    EigenMatrixXf oldAdjacencyMatrix = currentGraph->adjacencyMatrix;

    // initialize the new adjacency matrix to something that does not satisfy the check
    Eigen::MatrixXf newAdjacencyMatrix = Eigen::MatrixXf::Constant(oldAdjacencyMatrix.size(),-1); 

    std::vector<MatrixIdx> weightsToAvoid;
    // compute a gradient that will satisfy the check
    while(weightsToAvoid.size()){
        // compute the gradient
        auto graphGradient = computeObjGradient(currentGraph, weightsToAvoid);
        // apply the gradient
        auto newAdjacencyMatrix = oldAdjacencyMatrix + eps * graphGradient;
        // scale the weights (projection step)
        newAdjacencyMatrix = newAdjacencyMatrix * oldAdjacencyMatrix.lpNorm<1>() / newAdjacencyMatrix.lpNorm<1>();
        weightsToAvoid = getNegativeMatrixIndices(newAdjacencyMatrix);
    }

    // Need to generate a graph from the current graph with modified adjacency matrix
    return currentGraph->applyGradient(newAdjacencyMatrix);
}


int main(){

    // Initialize graph
    int x{10};
    int y{10};
    std::shared_ptr<Graph> graphInit = std::make_shared<Graph>();
    graphInit->constructSimpleGraph(x,y);
    // For loop
    std::vector<std::shared_ptr<Graph>> graphHistory;
    graphHistory.push_back(graphInit);
    int nIter{5};
    for (int i{0}; i<nIter; i++){
        // Gradient descent
        graphHistory.push_back(graphGradientDescent(graphHistory.back()));    
        // (Keep a history of the graph)
        // generate necessary plots
        // generatePlots(graphHistory.back());
        // save plots
        // savePlots();
    }


    return 0;
}
