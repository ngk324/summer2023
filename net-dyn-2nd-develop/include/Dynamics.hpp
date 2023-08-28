#ifndef _DYNAMICS_H_
#define _DYNAMICS_H_

#include "Graph.hpp"
#include "Plot.hpp"
#include "Force.hpp"

#include <memory>

#include "EnergyPlot.hpp"
#include "gnuplot-iostream.h"

class Dynamics
{
public:
    int simTime;
    int simSteps;
    double dampingCoeff;
    double stiffnessCoeff;
    double epsilon{0.01};

    Dynamics(int sim_time, int sim_steps, double damping, double stiffness);

    Eigen::VectorXf getStateVector(Graph &g) const;
    void setNodeStates(Graph &g, Eigen::VectorXf &states) const;

    void runCentralizedDynamics(Graph &g, Force &force, Plot &plot);
    void runDecentralizedDynamics(std::vector<std::shared_ptr<Node>> &nodes, Force &force, Plot &plot) const;
    Gnuplot energyPlotStream;
    std::vector<double> energyValueHistory;
    std::vector<double> XValueHistory;

};

#endif
