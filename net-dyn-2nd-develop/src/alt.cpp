#include "../include/Graph.hpp"
#include <algorithm>
#include <memory>
#include <numeric>
#include <vector>
#include <Eigen/Eigenvalues>
#include "../include/Histogram.hpp"
#include "../include/Plot.hpp"

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
    double lambda_step{0.1};
    double maxVal{currentGraph->eigenValues[currentGraph->eigenValues.size()-1]};

    for(int j{0}; j<currentGraph->eigenValues.size(); j++){
        double lambda_j = currentGraph->eigenValues[j];
        if(lambda_j<lambda_step)
            continue;
        for(int k{j+1}; k<currentGraph->eigenValues.size(); k++){
            // double currentVal = evalObjectiveFunctionAtJK(h,currentGraph->eigenValues[j],currentGraph->eigenValues[k]);
            double lambda_k = currentGraph->eigenValues[k];
            double currentVal = lambda_k-lambda_j;
            if(currentVal<maxVal){
                max_j = j;
                max_k = k;
                maxVal = currentVal;
            }
        }
    }
    std::cout << "min separation: " << maxVal << std::endl;
    return MatrixIdx(max_j,max_k);
}

std::vector<MatrixIdx> getNegativeMatrixIndices(const Eigen::MatrixXf &adjMatrix){

    std::vector<MatrixIdx> res;
    for(int j{0}; j<adjMatrix.rows(); j++){
        for(int k{j+1}; k<adjMatrix.cols(); k++){
            if(adjMatrix(j,k)<0.2 && adjMatrix(j,k)>0)
                res.push_back(MatrixIdx(j,k));
        }
    }
    return res;
}

double evaluateObjGradientAtJK(const double h, const double lambda_j, const double lambda_k){
    
    return 8*h*h*pow(sqrt(lambda_j)-sqrt(lambda_k),2)/(sqrt(lambda_j*lambda_k)*pow(pow(sqrt(lambda_j)-sqrt(lambda_k),2)+h*h,3));
}

Eigen::MatrixXf computeObjGradient(const std::shared_ptr<Graph> currentGraph, const std::vector<MatrixIdx> &weightsToAvoid){

//     // if the objective function contains a max function find the necessary indices for the eigenvalues that maximize
//     double h{1};
//     auto maxEigIdx = findObjMaximizingEigIndex(currentGraph,h);
//     // std::cout << "Chosen indices: j=" << maxEigIdx.j << ", k=" << maxEigIdx.k << std::endl;
//     // Compute the gradient using the derivation
//     Eigen::MatrixXf gradient = Eigen::MatrixXf::Zero(currentGraph->adjacencyMatrix.rows(),currentGraph->adjacencyMatrix.cols());
//     double lambda_j = currentGraph->eigenValues[maxEigIdx.j];
//     double lambda_k = currentGraph->eigenValues[maxEigIdx.k];
//     // std::cout << "lambda_j=" << lambda_j << std::endl;
//     // std::cout << "lambda_k=" << lambda_k << std::endl;
//     // gradient(maxEigIdx.j,maxEigIdx.k) = evaluateObjGradientAtJK(h,lambda_j,lambda_k);
//     // gradient(maxEigIdx.k,maxEigIdx.j) = evaluateObjGradientAtJK(h,lambda_j,lambda_k);
//     Eigen::VectorXf u_j = currentGraph->eigenVectors.col(maxEigIdx.j);
//     Eigen::VectorXf u_k = currentGraph->eigenVectors.col(maxEigIdx.k);
//     // std::cout << u_j << std::endl;
//     // std::cout << u_k << std::endl;
//     for(int i{0}; i<gradient.rows(); i++){
//         for(int l{i+1}; l<gradient.cols(); l++){
//             gradient(i,l) = pow(u_k(i)-u_k(l),2) - pow(u_j(i)-u_j(l),2);
//             gradient(l,i) = gradient(i,l);
// 
// //             if(lambda_j>lambda_k){
// //                 gradient(i,l) = pow(u_j(i)-u_j(l),2);
// //                 gradient(l,i) = gradient(i,l);
// //             }
// //             else{
// //                 gradient(i,l) = pow(u_k(i)-u_k(l),2);
// //                 gradient(l,i) = gradient(i,l);
// //             }
//         }
//     }
// 
//     // Remove the components corresponding to non existing edges
//     // std::cout << gradient << std::endl;
//     // std::cout << "gradient sum(1): " << gradient.sum() << std::endl;
//     // std::cout << "gradient=\n" << gradient << std::endl;
//     // std::cout << "connectivity=\n" << currentGraph->connectivityMatrix << std::endl;
//
    double objSum{0};
    for(int j{0};j<currentGraph->eigenValues.size(); j++){
        for(int k{0};k<currentGraph->eigenValues.size(); k++){
            objSum += pow(currentGraph->eigenValues[j]-currentGraph->eigenValues[k],2);

        }
    }
    std::cout << "objSum: " << objSum << std::endl;

    Eigen::MatrixXf gradient = Eigen::MatrixXf::Zero(currentGraph->adjacencyMatrix.rows(), currentGraph->adjacencyMatrix.cols());
    double eigSum = std::accumulate(currentGraph->eigenValues.begin(), currentGraph->eigenValues.end(),0);
    int eigSize = currentGraph->eigenValues.size();
    for(int j{0}; j<gradient.rows(); j++){
        for(int k{j+1}; k<gradient.cols(); k++){
            double gradSum{0};
            for(int i{0}; i<currentGraph->eigenValues.size(); i++){
                gradSum += (4*eigSize*currentGraph->eigenValues[i] - 4*eigSum)*pow(currentGraph->eigenVectors(j,i)-currentGraph->eigenVectors(k,i),2);
            }
            gradient(j,k) = gradSum;
            gradient(k,j) = gradSum;
        }
    }
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
    double eps{0.0001};
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
        if(weightsToAvoid.size() != 0){
            // for(auto &el: weightsToAvoid)
                 // std::cout << "modified gradient element at " << el.j << ", " << el.k << " = " << graphGradient(el.j,el.k) << std::endl;
            // std::cout << scaledAdjacencyMatrix << std::endl;
        }
        // std::cout << "grad:\n" << graphGradient << std::endl;
        // std::cout << graphGradient.sum() << std::endl;
        // apply the gradient
        std::cout << "eps: " << eps << std::endl;
        auto newAdjacencyMatrix = oldAdjacencyMatrix + eps * graphGradient;
        // std::cout << "new:\n" << newAdjacencyMatrix << std::endl;
        // scale the weights (projection step)
        scaledAdjacencyMatrix = newAdjacencyMatrix * oldAdjacencyMatrix.sum()/ newAdjacencyMatrix.sum();
        // scaledAdjacencyMatrix = newAdjacencyMatrix;
        // std::cout << "scaled:\n" << scaledAdjacencyMatrix << std::endl;
        weightsToAvoid = getNegativeMatrixIndices(scaledAdjacencyMatrix);
        //         for(auto &el: weightsToAvoid){
        //             std::cout << "avoid j,k: " << el.j << ", " << el.k << std::endl;
        //         }
        if(weightsToAvoid.size() != 0){
            // for(auto &el: weightsToAvoid)
                // std::cout << "invalid A element at " << el.j << ", " << el.k << " = " << scaledAdjacencyMatrix(el.j,el.k) << std::endl;
            // std::cout << scaledAdjacencyMatrix << std::endl;
            // std::cout << oldAdjacencyMatrix << std::endl;
            eps /= 10;
        }
        if(eps < 1E-5){
            return currentGraph;
        }
    } while(!weightsToAvoid.empty());
    if(abs(scaledAdjacencyMatrix.sum() - oldAdjacencyMatrix.sum())>0.1){
        std::cout << "old sum: " << oldAdjacencyMatrix.sum() << std::endl;
        std::cout << "new sum: " << scaledAdjacencyMatrix.sum() << std::endl;
    }
    // Need to generate a graph from the current graph with modified adjacency matrix
    return currentGraph->applyGradient(scaledAdjacencyMatrix);
}


int main(int argc, char *argv[]){
    int size{10};
    if(argc>1){
        size = std::stoi(argv[1]);
    }

    // Initialize graph
    int x{size};
    int y{size};
    std::shared_ptr<Graph> graphInit = std::make_shared<Graph>();
    graphInit->constructSimpleGraph(x,y);
    graphInit->computeMatrices();
    graphInit->eigenDecompose();
    Gnuplot gp_init;
    histogram::generateHistogram(gp_init,graphInit->eigenValues);
    Plot my_plot_init("Graph Plot Init", 500/size, 2, 2, size, size);
    my_plot_init.plotGraph(*graphInit);
    my_plot_init.displayPlot(true);
    // For loop
    std::vector<std::shared_ptr<Graph>> graphHistory;
    graphHistory.push_back(graphInit);
    int nIter{1000000};
    Gnuplot gp;
    for (int i{0}; i<nIter; i++){
        std::cout << "Iter: " << i << std::endl;
        // Gradient descent
        auto nextGraph = graphGradientDescent(graphHistory.back());
        if(nextGraph == graphHistory.back())
            break;
        graphHistory.push_back(nextGraph);
        // std::cout << "new adj=\n" << nextGraph->adjacencyMatrix << std::endl;
        // std::cout << "new conn=\n" << nextGraph->connectivityMatrix << std::endl;
        nextGraph->eigenDecompose();
        std::cout << "EigenSum: " << std::accumulate(nextGraph->eigenValues.begin(),nextGraph->eigenValues.end(),0) << std::endl;
        histogram::generateHistogram(gp,nextGraph->eigenValues);
    }

    Plot my_plot("Graph Plot", 500/size, 2, 2, size, size);
    my_plot.plotGraph(*graphHistory.back());
    my_plot.displayPlot(true);

    gp << "exit\n";
    return 0;
}
