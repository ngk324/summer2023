#include "EnergyPlot.hpp"
#include <numeric>
#include <vector>

namespace energyPlot{
    void generateEnergyPlot(Gnuplot &plotEnergy, std::vector<double> &energyValues){

        double minVal = *std::min_element(energyValues.begin(), energyValues.end());
        double maxVal = *std::max_element(energyValues.begin(), energyValues.end());

        plotEnergy << "reset\n";
        plotEnergy << "min=" << minVal << "\n";
        plotEnergy << "max=" << maxVal << "\n";

        plotEnergy << "set yrange [min:max]\n";
        plotEnergy << "set xrange [0:" << energyValues.size()+1 << "]\n";
        plotEnergy << "set xlabel 'Iteration'\n";
        plotEnergy << "set ylabel 'Energy Value'\n";

        plotEnergy << "plot '-' u 1:2 with lines lc rgb'red' notitle\n"; 

        for (size_t i = 0; i < energyValues.size(); ++i) {
            plotEnergy << i << " " << energyValues[i] << "\n";
        }

        plotEnergy << "e\n";

        plotEnergy << "set terminal png\n";
        plotEnergy << "set output 'EnergyPlot.png'\n";
        plotEnergy << "replot\n";
       // plotEnergy << "set output\n";

        plotEnergy.flush();
        usleep(10000);
    }
}