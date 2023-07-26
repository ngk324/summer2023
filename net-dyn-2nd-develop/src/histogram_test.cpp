#include "../include/Graph.hpp"
#include <algorithm>
#include <memory>
#include <vector>
#include "../include/Histogram.hpp"

int main(){

    // Initialize graph
    int x{10};
    int y{10};
    std::shared_ptr<Graph> graphInit = std::make_shared<Graph>();
    graphInit->constructSimpleGraph(x,y);
    graphInit->computeMatrices();
    graphInit->eigenDecompose();
    std::cout << graphInit->eigenValues << std::endl;

    Gnuplot gp;
    histogram::generateHistogram(gp,graphInit->eigenValues);    
    return 0;
}
