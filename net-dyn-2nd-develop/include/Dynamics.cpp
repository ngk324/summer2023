#include "Dynamics.hpp"
#include "EnergyPlot.hpp"
#include "XPlot.hpp"

#include <unistd.h>
#include <algorithm>
#include <memory>
#include <numeric>
#include <vector>

Dynamics::Dynamics(int sim_time, int sim_steps, double damping, double stiffness)
    : simTime{sim_time}, simSteps{sim_steps}, dampingCoeff{damping}, stiffnessCoeff{stiffness} {}

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

void Dynamics::runCentralizedDynamics(Graph &g, Force &force, Plot &plot)
{
    plot.displayMethod("Centralized");
    int nNodes = g.nodes.size();
    Eigen::MatrixXf A_Matrix = dampingCoeff * (epsilon * Eigen::MatrixXf::Identity(nNodes, nNodes));
    Eigen::MatrixXf B_Matrix = (g.laplacianMatrix + epsilon * Eigen::MatrixXf::Identity(nNodes, nNodes));
    Eigen::VectorXf x = getStateVector(g);
    Eigen::VectorXf x_dot = Eigen::VectorXf::Zero(nNodes);
    Eigen::VectorXf x_ddot(nNodes);
    double timeStep = double(simTime) / simSteps;
    for (int i{0}; i < simSteps + 1; i++)
    {
        x_ddot = force.sinCauchyForce(i * timeStep) - A_Matrix * x_dot - B_Matrix * x;
        // x_ddot = getForcing(i * timeStep, 1, 6.13602, nNodes) - A_Matrix * x_dot - B_Matrix * x;
        x_dot += (x_ddot * timeStep);
        x += (x_dot * timeStep);

        /*double energyVal = (.5 * x_dot.transpose()*x_ddot);
        energyVal += + (.5 * x_dot.transpose()*g.laplacianMatrix*x); */

        double energyVal = (.5 * x_dot.transpose()*Eigen::MatrixXf::Identity(nNodes, nNodes)*x_dot);
        energyVal += (.5 * x.transpose()*g.laplacianMatrix*x); 
        //energyValueHistory.push_back(energyVal);
        energyValueHistory.push_back(energyVal);
        XValueHistory.push_back(x[1]);
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
        if(i > 25000){
            std::cout << "\n" << 100 * ((energyValueHistory[energyValueHistory.size() - 25000] - energyValueHistory.back()) / energyValueHistory.back());
            //energyPlot::generateEnergyPlot(energyPlotStream, energyValueHistory);
        }
    }
    //plotEnergy();
    energyPlot::generateEnergyPlot(energyPlotStream, energyValueHistory);
    XPlot::generateXPlot(energyPlotStream, XValueHistory);

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
/*
void Dynamics::plotEnergy(){
    energyPlot::generateEnergyPlot(energyPlotStream, energyValueHistory);
}*/