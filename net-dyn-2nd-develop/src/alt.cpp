#include "include/Graph.hpp"
#include <algorithm>
#include <memory>
#include <vector>

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
    
    // (Keep a history of the graph)
    }


    return 0;
}
