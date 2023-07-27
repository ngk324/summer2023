#ifndef _GRADIENTDESCENT_H_
#define _GRADIENTDESCENT_H_

#include <Eigen/Dense>
#include <memory>

#include "Graph.hpp"
#include "gnuplot-iostream.h"

class GradientDescent{

public:
    GradientDescent(std::shared_ptr<Graph> initGraph);
    void runNStepDescent(const int nIter);
    void destroyHistogram();
private:
    int graphGridSize;
    Gnuplot histogramStream;
    double gradientStep;    
    void decreaseGradientStep();
    std::vector<std::shared_ptr<Graph>> graphHistory;
    void runOneStepDescent();
    void plotHistogram();
    void plotGraph();
};



#endif
