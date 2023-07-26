#include "../include/Graph.hpp"
#include <algorithm>
#include <memory>
#include <numeric>
#include <vector>
#include <Eigen/Eigenvalues>
#include "../include/Histogram.hpp"

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
// double l1NormOfMatrix(const Eigen::MatrixXf &mat);
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
            if(adjMatrix(j,k)<0.5 && adjMatrix(j,k)>0)
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
    // std::cout << "Chosen indices: j=" << maxEigIdx.j << ", k=" << maxEigIdx.k << std::endl;
    // Compute the gradient using the derivation
    Eigen::MatrixXf gradient = Eigen::MatrixXf::Zero(currentGraph->adjacencyMatrix.rows(),currentGraph->adjacencyMatrix.cols());
    double lambda_j = currentGraph->eigenValues[maxEigIdx.j];
    double lambda_k = currentGraph->eigenValues[maxEigIdx.k];
    // std::cout << "lambda_j=" << lambda_j << std::endl;
    // std::cout << "lambda_k=" << lambda_k << std::endl;
    // gradient(maxEigIdx.j,maxEigIdx.k) = evaluateObjGradientAtJK(h,lambda_j,lambda_k);
    // gradient(maxEigIdx.k,maxEigIdx.j) = evaluateObjGradientAtJK(h,lambda_j,lambda_k);
    Eigen::VectorXf u_j = currentGraph->eigenVectors.col(maxEigIdx.j);
    Eigen::VectorXf u_k = currentGraph->eigenVectors.col(maxEigIdx.k);
    // std::cout << u_j << std::endl;
    // std::cout << u_k << std::endl;
    for(int i{0}; i<gradient.rows(); i++){
        for(int l{i+1}; l<gradient.cols(); l++){
//             gradient(i,l) = pow(u_j(i)-u_j(l),2) - pow(u_k(i)-u_k(l),2);
//             gradient(l,i) = pow(u_j(i)-u_j(l),2) - pow(u_k(i)-u_k(l),2);
//             if(lambda_j<lambda_k){
//                 gradient(i,l) = -gradient(i,l);
//                 gradient(l,i) = -gradient(l,i);
//             }
            if(lambda_j>lambda_k){
                gradient(i,l) = pow(u_j(i)-u_j(l),2);
                gradient(l,i) = gradient(i,l);
            }
            else{
                gradient(i,l) = pow(u_k(i)-u_k(l),2);
                gradient(l,i) = gradient(i,l);
            }
        }
    }

    // Remove the components corresponding to non existing edges
    // std::cout << gradient << std::endl;
    // std::cout << "gradient sum(1): " << gradient.sum() << std::endl;
    // std::cout << "gradient=\n" << gradient << std::endl;
    // std::cout << "connectivity=\n" << currentGraph->connectivityMatrix << std::endl;
    Eigen::MatrixXf gradientOfExistingWeights = gradient.array() * currentGraph->connectivityMatrix.array(); 
    // std::cout << "gradientofexisting=\n" << gradientOfExistingWeights << std::endl;
    // Remove the components corresponding to weightsToAvoid
    for(auto &el: weightsToAvoid){
        gradientOfExistingWeights(el.j,el.k) = 0;
        gradientOfExistingWeights(el.k,el.j) = 0;
    }

    return gradientOfExistingWeights;
}

std::shared_ptr<Graph> graphGradientDescent(const std::shared_ptr<Graph> currentGraph){
    double eps{1};
    // std::cout << "eigVals= " << currentGraph->eigenValues << std::endl;
    Eigen::MatrixXf oldAdjacencyMatrix = currentGraph->adjacencyMatrix;
    // std::cout << "old:\n" << oldAdjacencyMatrix << std::endl;
    // initialize the new adjacency matrix to something that does not satisfy the check
    Eigen::MatrixXf newAdjacencyMatrix = Eigen::MatrixXf::Constant(oldAdjacencyMatrix.rows(), oldAdjacencyMatrix.cols(), -1); 
    Eigen::MatrixXf scaledAdjacencyMatrix = newAdjacencyMatrix;
    std::vector<MatrixIdx> weightsToAvoid;
    // compute a gradient that will satisfy the check
    do{
        // compute the gradient
        auto graphGradient = computeObjGradient(currentGraph, weightsToAvoid);
        // std::cout << "grad:\n" << graphGradient << std::endl;
        // std::cout << graphGradient.sum() << std::endl;
        // apply the gradient
        auto newAdjacencyMatrix = oldAdjacencyMatrix + eps * graphGradient;
        // std::cout << "new:\n" << newAdjacencyMatrix << std::endl;
        // scale the weights (projection step)
        // scaledAdjacencyMatrix = newAdjacencyMatrix * oldAdjacencyMatrix.sum()/ newAdjacencyMatrix.sum();
        scaledAdjacencyMatrix = newAdjacencyMatrix;
        // std::cout << "scaled:\n" << scaledAdjacencyMatrix << std::endl;
        weightsToAvoid = getNegativeMatrixIndices(scaledAdjacencyMatrix);
        //         for(auto &el: weightsToAvoid){
        //             std::cout << "avoid j,k: " << el.j << ", " << el.k << std::endl;
        //         }
    } while(!weightsToAvoid.empty());

    // Need to generate a graph from the current graph with modified adjacency matrix
    return currentGraph->applyGradient(scaledAdjacencyMatrix);
}


int main(){

    // Initialize graph
    int x{10};
    int y{10};
    std::shared_ptr<Graph> graphInit = std::make_shared<Graph>();
    graphInit->constructSimpleGraph(x,y);
    graphInit->computeMatrices();
    graphInit->eigenDecompose();
    // For loop
    std::vector<std::shared_ptr<Graph>> graphHistory;
    graphHistory.push_back(graphInit);
    int nIter{1000000};
    Gnuplot gp;
    for (int i{0}; i<nIter; i++){
        std::cout << "Iter: " << i << std::endl;
        // Gradient descent
        auto nextGraph = graphGradientDescent(graphHistory.back());
        graphHistory.push_back(nextGraph);
        // std::cout << "new adj=\n" << nextGraph->adjacencyMatrix << std::endl;
        // std::cout << "new conn=\n" << nextGraph->connectivityMatrix << std::endl;
        nextGraph->eigenDecompose();
        histogram::generateHistogram(gp,nextGraph->eigenValues);
    }

    return 0;
}
