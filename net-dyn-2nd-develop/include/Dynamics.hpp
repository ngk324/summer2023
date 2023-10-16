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
    double epsilon;

    Dynamics(int sim_time, int sim_steps, double damping, double stiffness, double epsilon);

    Eigen::VectorXf getStateVector(Graph &g) const;
    void setNodeStates(Graph &g, Eigen::VectorXf &states) const;

    void writeNodeAvgFile(std::vector<double> nodeValsMax, double avg);
    void writeTwoNormAvgFile(double avg);
    void writeNodeValuesFile(std::vector<std::vector<double>> XValueHistory, int nodeSize, int simSteps);
    double calculateTwoNormVals(std::vector<double> XValueHistory, int startTime, int windowSize);
    double inverse_of_normal_cdf(const double p, const double mu, const double sigma);
    void write_file_results(std::string print);

    void runCentralizedDynamics(Graph &g, Force &force, Plot &plot);
    void runDecentralizedDynamics(std::vector<std::shared_ptr<Node>> &nodes, Force &force, Plot &plot) const;
    std::vector<double> calculateNodeVals(std::vector<std::vector<double>> XHistory, int startTime, int windowSize);
    bool determineSteadyState(std::vector<double> energyValueHistory, int iterationRange, double percentDifference);
    Gnuplot energyPlotStream;
    std::vector<double> energyValueHistory;
    std::vector<double> XValueHistory1;

};

#endif
