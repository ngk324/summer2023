#include "Dynamics.hpp"
#include "EnergyPlot.hpp"
#include "XPlot.hpp"

#include <unistd.h>
#include <algorithm>
#include <memory>
#include <numeric>
#include <vector>
#include <string>
#include <string.h>

Dynamics::Dynamics(int sim_time, int sim_steps, double damping, double stiffness, double epsilon, int simNum, int seed)
    : simTime{sim_time}, simSteps{sim_steps}, dampingCoeff{damping}, stiffnessCoeff{stiffness}, epsilon{epsilon}, simNum{simNum}, seed{seed} {}

Eigen::VectorXf Dynamics::getStateVector(Graph &g) const
{
    int nNodes = g.nodes.size();
    Eigen::VectorXf vec(nNodes);
    for (int i{0}; i < nNodes; i++)
    {
        vec(i) = g.nodes[i]->z;
    }
    return vec;
}

void Dynamics::setNodeStates(Graph &g, Eigen::VectorXf &states) const
{
    for (int i{0}; i < states.size(); i++)
    {
        g.nodes[i]->z_old = g.nodes[i]->z;
        g.nodes[i]->z = states(i);
    }
}

void Dynamics::writeNodeAvgFile(std::vector<double> nodeValsMax, double avg){
    std::string fileName = "NodeAvgResults";
    fileName.append(std::to_string((simNum +1)%2));
    fileName = fileName + "-";
    fileName.append(std::to_string(seed));
    fileName = fileName + ".txt";
    std::string line;
    std::ofstream myFile(fileName);
    for(int i{0}; i <= nodeValsMax.size(); i++){
        if(i < nodeValsMax.size()){
            myFile << "Node " << i << ": " << nodeValsMax[i] << "\n";
        }
        else{
            myFile << "Node Avg: " << ": " << avg;
        }
    }
    myFile.close();
}

void Dynamics::writeNodeValuesFile(std::vector<std::vector<double>> XValueHistory, int nodeSize, int simSteps){
    std::string fileName = "NodeValueResults";
    fileName.append(std::to_string((simNum + 1)%2));
    fileName = fileName + "-";
    fileName.append(std::to_string(seed));
    fileName = fileName + ".txt";
    std::string line;
    std::ofstream myFile(fileName);

    for(int j{0}; j < nodeSize; j++){
        myFile << "Node " << j << ":\n";
        for (int i{0}; i < simSteps; i++)
        {
            if(i % 99 == 0){
                myFile << XValueHistory[j][i] << ", ";
            }
        }
        myFile << "\n\n\n";
    }
    myFile.close();
}

void Dynamics::runCentralizedDynamics(Graph &g, Force &force, Plot &plot)
{
    plot.displayMethod("Centralized");
    int nNodes = g.nodes.size();
    Eigen::MatrixXf A_Matrix = dampingCoeff * (Eigen::MatrixXf::Identity(nNodes, nNodes));
    Eigen::MatrixXf B_Matrix = (g.laplacianMatrix + epsilon * Eigen::MatrixXf::Identity(nNodes, nNodes));
    Eigen::VectorXf x = getStateVector(g);
    Eigen::VectorXf x_dot = Eigen::VectorXf::Zero(nNodes);
    Eigen::VectorXf x_ddot(nNodes);
    double timeStep = double(simTime) / simSteps;    

    //std::vector<double> XValueHistory[g.nodes.size()][simSteps];
    std::vector<std::vector<double>> XValueHistory(g.nodes.size(), std::vector<double> (simSteps, 0));

    for (int i{0}; i < simSteps + 1; i++)
    {
        //x_ddot = force.sinCauchyForce(i * timeStep) - A_Matrix * x_dot - B_Matrix * x;
        x_ddot = force.sinCauchyForce(i * timeStep) - A_Matrix * x_dot - B_Matrix * x;
        x_dot += (x_ddot * timeStep);
        x += (x_dot * timeStep);

        double energyVal = (.5 * x_dot.transpose()*Eigen::MatrixXf::Identity(nNodes, nNodes)*x_dot);
        energyVal += (.5 * x.transpose()*g.laplacianMatrix*x); 

        energyValueHistory.push_back(energyVal);

        for(int j{0}; j < g.nodes.size(); j++){
            XValueHistory[j][i] = (double)x[j];
        }

        XValueHistory1.push_back(x[1]);

        setNodeStates(g, x);

        double minEV = g.nodes[0]->z;
        double maxEV = g.nodes[0]->z;

        for(int i = 0; i < g.nodes.size(); i++){
            if(g.nodes[i]->z > maxEV){
                double maxEV = g.nodes[i]->z;
            }
            if(g.nodes[i]->z < minEV){
                double minEV = g.nodes[i]->z;
            }
        }
        for (int j{0}; j < nNodes; j++)
        {
            plot.plotNode(*g.nodes[j], maxEV, minEV);
            plot.displayState(*g.nodes[j]);
        }
        // std::cout << x << std::endl
        //           << std::endl;
        plot.displayTime(std::to_string(i * timeStep) + " s");
        plot.displayPlot();
        usleep(1E+3 * timeStep);
        /*if(i > 25000){
            std::cout << "\n" << 100 * ((energyValueHistory[energyValueHistory.size() - 25000] - energyValueHistory.back()) / energyValueHistory.back());
            //energyPlot::generateEnergyPlot(energyPlotStream, energyValueHistory);
        }*/
    }
    
    //int numOfWindows = 5;
    //int widthOfWindow = 20000;
    //int startTime = 400000;
    int numOfWindows = 5;
    int widthOfWindow = 20000;
    int startTime = 400000;
    std::vector<double> nodeValsMax(nNodes, 0);
    for(int i{0}; i < numOfWindows; i++){
        std::vector<double> nodeMaxVector = calculateNodeVals(XValueHistory, startTime + i * widthOfWindow, widthOfWindow);
        for(int j{0}; j < nodeValsMax.size(); j++){
            nodeValsMax[j] += (nodeMaxVector[j]/numOfWindows);
        }
    }
    double avg = 0;
    for(int i{0}; i <= nodeValsMax.size(); i++){
        if(i < nodeValsMax.size()){
            avg += nodeValsMax[i];
        }
        else{
            avg = avg / nodeValsMax.size();
        }
    }
    writeNodeAvgFile(nodeValsMax, avg);
    writeNodeValuesFile(XValueHistory, g.nodes.size(), simSteps+1);

    energyPlot::generateEnergyPlot(energyPlotStream, energyValueHistory, (simNum + 1)%2, seed);
    XPlot::generateXPlot(energyPlotStream, XValueHistory1, (simNum + 1)%2, seed);
}

std::vector<double> Dynamics::calculateNodeVals(std::vector<std::vector<double>> XValueHistory, int startTime, int windowSize){ // end time = numOfWindows*windowSize + startTime
    std::vector<double> nodeMax(XValueHistory.size(), 0);
    for(int i = startTime; i < startTime + windowSize; i++){
        for(int j{0}; j < XValueHistory.size(); j++){
            if(nodeMax[j] < abs(XValueHistory[j][i])){
                nodeMax[j] = abs(XValueHistory[j][i]);
            }
        }
    }
    return nodeMax;
}

bool Dynamics::determineSteadyState(std::vector<double> energyValueHistory, int iterationRange, double percentDifference){
    bool withinPercent = false;
    double changeFromRange = ((energyValueHistory[energyValueHistory.size() - iterationRange - 1] - energyValueHistory.back()) / energyValueHistory.back());
    double changeFromPrevious = ((energyValueHistory[energyValueHistory.size() - 2] - energyValueHistory.back()) / energyValueHistory.back());
    if((changeFromRange - changeFromPrevious) / changeFromPrevious < 0.1){
        withinPercent = true;
    }
    return withinPercent;
}

void Dynamics::runDecentralizedDynamics(std::vector<std::shared_ptr<Node>> &nodes, Force &force, Plot &plot) const
{
    plot.displayMethod("Decentralized");
    double timeStep = double(simTime) / simSteps;
    for (int i{0}; i < simSteps + 1; i++)
    {
        Eigen::VectorXf force_vec = force.sinusoidalForce(i * timeStep);
        for (int j{0}; j < nodes.size(); j++)
        {
            double neighbor_z_sum{0}, neighbor_zdot_sum{0};
            for (int k{0}; k < nodes[j]->neighbors.size(); k++)
            {
                neighbor_z_sum += nodes[j]->neighbors[k]->z_old;
                neighbor_zdot_sum += nodes[j]->neighbors[k]->z_dot_old;
            }
            double z_ddot = force_vec(j) - dampingCoeff * (nodes[j]->neighbors.size() * nodes[j]->z_dot - neighbor_zdot_sum + epsilon * nodes[j]->z_dot) - stiffnessCoeff * (nodes[j]->neighbors.size() * nodes[j]->z - neighbor_z_sum + epsilon * nodes[j]->z);
            nodes[j]->z += (nodes[j]->z_dot * timeStep);
            nodes[j]->z_dot += (z_ddot * timeStep);
        }
        for (int j{0}; j < nodes.size(); j++)
        {
            nodes[j]->z_old = nodes[j]->z;
            nodes[j]->z_dot_old = nodes[j]->z_dot;
            //plot.plotNode(*nodes[j]);
            plot.displayState(*nodes[j]);
        }

        plot.displayTime(std::to_string(i * timeStep) + " s");
        plot.displayPlot();
        usleep(1E+2 * timeStep);
    }
}
