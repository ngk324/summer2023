#include "GradientDescent.hpp"
#include "Histogram.hpp"
#include "Plot.hpp"
#include <memory>
#include <numeric>
#include <vector>

GradientDescent::GradientDescent(std::shared_ptr<Graph> initGraph): graphGridSize{initGraph->gridSize}{
    graphHistory.push_back(initGraph);
}

void GradientDescent::runNStepDescent(const int nIter){
    
    plotHistogram();
    plotGraph();
    for(int i{0}; i<nIter; i++){
        printIterInfo(i);
        runOneStepDescent();
        if(graphHistory.size() < i+1){
            std::cout << "Stuck at iter: " << i+1 << std::endl;
            break;
        }
        plotHistogram();
    }
    plotGraph();
}

void GradientDescent::runOneStepDescent(){
    Eigen::MatrixXf oldAdjacencyMatrix = graphHistory.back()->adjacencyMatrix; 
    Eigen::MatrixXf newAdjacencyMatrix = oldAdjacencyMatrix; 
    Eigen::MatrixXf scaledAdjacencyMatrix = newAdjacencyMatrix; 
    std::vector<MatrixIdx> weightsToAvoid; 
    int nRecompute{0};
    do{
        if(nRecompute > maxRecompute)
            return;
        auto adjGradient = computeAdjGradient(graphHistory.back(), weightsToAvoid); 
        std::cout << adjGradient << std::endl;
        newAdjacencyMatrix = oldAdjacencyMatrix + gradientStep * adjGradient; 
        scaledAdjacencyMatrix = newAdjacencyMatrix * oldAdjacencyMatrix.sum() / newAdjacencyMatrix.sum(); 
        nRecompute++;
        weightsToAvoid = getInvalidWeightIdx(scaledAdjacencyMatrix);
    }while(!weightsToAvoid.empty());
    graphHistory.push_back(graphHistory.back()->applyGradient(scaledAdjacencyMatrix));
}

Eigen::MatrixXf GradientDescent::computeAdjGradient(const std::shared_ptr<Graph> graph, std::vector<MatrixIdx> &weightsToAvoid) const{

    // double sum gradient
    Eigen::MatrixXf gradientMat = Eigen::MatrixXf::Zero(graph->adjacencyMatrix.rows(),graph->adjacencyMatrix.cols());
    double sumEig = std::accumulate(graph->eigenValues.begin(), graph->eigenValues.end(), 0);
    int nEig = graph->eigenValues.size();
    for(int j{0}; j<gradientMat.rows(); j++){
        auto u_row_j = graph->eigenVectors.row(j);
        for(int k{j+1}; k<gradientMat.cols(); k++){
            auto u_row_k = graph->eigenVectors.row(k);
            if(graph->connectivityMatrix(j,k) == 0)
                continue;
            double gradAtJK{0};
            for(int i{0}; i<nEig; i++){
                double lambda_i = graph->eigenValues[i];
                gradAtJK += 4*pow(u_row_j[i]-u_row_k[i],2)*(nEig*lambda_i - sumEig); 
            }
            gradientMat(j,k) = gradAtJK;
            gradientMat(k,j) = gradAtJK;
        }
    }
    for(auto &el: weightsToAvoid){
        gradientMat(el.j,el.k) = 0;
        gradientMat(el.k,el.j) = 0;
    }
    return gradientMat;
}

std::vector<MatrixIdx> GradientDescent::getInvalidWeightIdx(const Eigen::MatrixXf &adjMat) const{
    std::vector<MatrixIdx> res;
    for(int j{0}; j<adjMat.rows(); j++){
        for(int k{j+1}; k<adjMat.cols(); k++){
            if(graphHistory.back()->connectivityMatrix(j,k) && adjMat(j,k) < minEdgeWeight)
                res.push_back(MatrixIdx(j,k));
        }
    }
    return res; 
}

void GradientDescent::plotHistogram(){

    if(!graphHistory.empty())
        histogram::generateHistogram(histogramStream,graphHistory.back()->eigenValues);

}

void GradientDescent::destroyHistogram(){
    histogramStream << "set term x11 close\n";
}

void GradientDescent::plotGraph(){

    Plot graphPlotter("Graph Plot", 500/graphGridSize, 2, 2, graphGridSize, graphGridSize);
    graphPlotter.plotGraph(*graphHistory.back());
    graphPlotter.displayPlot(true);

}

void GradientDescent::printIterInfo(const int iterNo) const{
    std::cout << "Iter #: " << iterNo+1 << std::endl;
    auto eigenvalues = graphHistory.back()->eigenValues;
    std::cout << "Eigenvalue sum: " << std::accumulate(eigenvalues.begin(), eigenvalues.end(), 0) << std::endl;
    double objValue{0};
    for(int j{0}; j<eigenvalues.size(); j++){
        for(int k{0}; k<eigenvalues.size(); k++){
            objValue += pow(eigenvalues[j] - eigenvalues[k], 2);
        }
    }
    std::cout << "Objective value: " << objValue << std::endl;
}
