#include "../include/GradientDescent.hpp"
#include "../include/Graph.hpp"

// My headers
#include "include/Plot.hpp"
#include "include/Dynamics.hpp"
#include "include/Force.hpp"

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string.h>
#include <string>

int main(int argc, char *argv[]){
    int gridSize;
    double alpha, freq;
    bool resultsGiven;
    bool weightConstraint{true};

    if(argc==5){
        gridSize = std::stoi(argv[1]);
        resultsGiven = std::stoi(argv[2]);
        alpha = std::stof(argv[3]);
        freq = std::stof(argv[4]);
    }
    else{
        return 0;
    }
    std::cout << alpha << " " << freq;
    std::shared_ptr<Graph> graphInit = std::make_shared<Graph>();
    graphInit->constructSimpleGraph(gridSize);
    std::vector<double> ev;

    if(resultsGiven){
        std::string line;
        std::ifstream myFile("LapResults.txt");
        while(getline(myFile, line))
        {
            std::istringstream lineStream(line);
            double first;
            lineStream >> first;
            ev.push_back(first);
        }
        int counter = 0;
        for (int i{0}; i < gridSize*gridSize; i++)
        {
            for (int j = 0; j < gridSize*gridSize; j++)
            {
                if(i == j){
                    graphInit->laplacianMatrix(i,j) = ev[counter];
                    graphInit->degreeMatrix(i,j) = ev[counter];
                }
                else{
                    graphInit->laplacianMatrix(i,j) = ev[counter];
                    graphInit->adjacencyMatrix(i,j) = -ev[counter];
                }
                counter++;
            }
        }
    }

    graphInit->eigenDecompose();
    GradientDescent gdObj(graphInit,weightConstraint);
    gdObj.plotHistogram();

    int MAX_X = gridSize;
    int MAX_Y = gridSize;
    int PLOT_SCALE = 40;
    int vPad = 2;
    int hPad = 2;
    double damping{0.1}, stiffness{5}, epsilon{0.01};
    double amp{1};
    bool decentralizedAlg = false;
    
    freq = sqrt(freq);
 
    // Generate plot
    Plot my_plot("State Plot - Chosen EigVal: " + std::to_string(freq), PLOT_SCALE, vPad, hPad, MAX_X, MAX_Y);
    my_plot.plotGraph(*graphInit);
    my_plot.displayPlot(true);
         
    Force my_force(amp, freq, graphInit->nodes.size(), alpha);
    my_force.insertForceElement(1);
     
    // Simulate dynamics
    int simulationTime{5000};
    int simulationSteps{simulationTime * 100};
    Dynamics my_sim(simulationTime, simulationSteps, damping, stiffness);
    if (decentralizedAlg)
        my_sim.runDecentralizedDynamics(graphInit->nodes, my_force, my_plot);
    else
        my_sim.runCentralizedDynamics(*graphInit, my_force, my_plot);
    my_plot.displayPlot(true);

    return 0;
}
