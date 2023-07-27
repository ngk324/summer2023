#include "../include/GradientDescent.hpp"
#include "../include/Graph.hpp"

int main(int argc, char *argv[]){

    int gridSize{10};
    int nIterFirst{100};

    if(argc==2){
        gridSize = std::stoi(argv[1]);
    }
    else if(argc>2){
        gridSize = std::stoi(argv[1]);
        nIterFirst = std::stoi(argv[2]);
    }

    // Initialize graph
    std::shared_ptr<Graph> graphInit = std::make_shared<Graph>();
    graphInit->constructSimpleGraph(gridSize);

    GradientDescent gdObj(graphInit);
    gdObj.runNStepDescent(nIterFirst);
    gdObj.destroyHistogram();

    return 0;
}
