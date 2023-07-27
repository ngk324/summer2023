#include "GradientDescent.hpp"
#include "Histogram.hpp"
#include "Plot.hpp"
#include <memory>

GradientDescent::GradientDescent(std::shared_ptr<Graph> initGraph){
    graphHistory.push_back(initGraph);
    graphGridSize = initGraph->gridSize;
}

void GradientDescent::runNStepDescent(const int nIter){
    
    plotHistogram();
    plotGraph();
    for(int i{0}; i<nIter; i++){
        runOneStepDescent();
        plotHistogram();
    }
    plotGraph();
}

void GradientDescent::runOneStepDescent(){
    // initialize old adj matrix
    //
    // initialize new adj matrix
    //
    // initialize scaled adj matrix
    //
    // initialize weights to avoid vector
    //
    // do-while loop (condition: there are weights to avoid)
    //
    // compute gradient
    //
    // add gradient to old adj matrix to obtain new adj matrix
    //
    // scale the adj matrix
    //
    // check if scaled adj matrix satisfies weight constraints (obtain weights to avoid if not)
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
