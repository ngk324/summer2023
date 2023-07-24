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
Eigen::MatrixXf computeObjGradient(const std::shared_ptr<Graph> currentGraph, const MatrixIdx maximizingIndices, const std::vector<MatrixIdx> &weightsToAvoid);

// =====================================================================================
double evalObjectiveFunctionAtJK(const double h, const double lambda_j, const double lambda_k){
    return ((4*h*h)/pow(pow(sqrt(lambda_j)-sqrt(lambda_k),2)+h*h,2));
}

MatrixIdx findObjMaximizingEigIndex(const std::shared_ptr<Graph> currentGraph, const double h){

    int max_j{0};
    int max_k{0};
    double maxVal{0.0};

    for(int j{0}; j<currentGraph->Eigenvalues; j++){
        for(int k{j+1}; k<currentGraph->Eigenvalues; k++){
            double currentVal = evalObjectiveFunctionAtJK(h,currentGraph->Eigenvalues[j],currentGraph->Eigenvalues[k]);
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

Eigen::MatrixXf computeObjGradient(const std::shared_ptr<Graph> currentGraph, const MatrixIdx maximizingIndices, const std::vector<MatrixIdx> &weightsToAvoid){
    // Compute the gradient using the derivation
    //
    // Remove the components corresponding to non existing edges
    //
    // Remove the components corresponding to weightsToAvoid
}

std::shared_ptr<Graph> graphGradientDescent(const  std::shared_ptr<Graph> currentGraph){

    double h{1};
    // if the objective function contains a max function find the necessary indices for the eigenvalues that maximize
    auto max_eig_idx = findObjMaximizingEigIndex(currentGraph,h);

    EigenMatrixXf oldAdjacencyMatrix = currentGraph->adjacencyMatrix;

    // initialize the new adjacency matrix to something that does not satisfy the check
    Eigen::MatrixXf newAdjacencyMatrix = Eigen::MatrixXf::Constant(oldAdjacencyMatrix.size(),-1); 

    std::vector<MatrixIdx> weightsToAvoid;
    // compute a gradient that will satisfy the check
    while(weightsToAvoid.size()){
        // compute the gradient
        auto graphGradient = computeObjGradient(currentGraph, max_eig_idx, weightsToAvoid);
        // apply the gradient
        auto newAdjacencyMatrix = oldAdjacencyMatrix + eps * graphGradient;
        // scale the weights (projection step)
        newAdjacencyMatrix = newAdjacencyMatrix * l1Norm(oldAdjacencyMatrix) / l1Norm(newAdjacencyMatrix);
        weightsToAvoid = getNegativeMatrixIndices(newAdjacencyMatrix);
    }

    auto descendedGraph = applyGradient(currentGraph, newAdjacencyMatrix);

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
        generatePlots(graphHistory.back());
        // save plots
        savePlots();
    }


    return 0;
}
