#include "EnergyPlot.hpp"
#include <numeric>
#include <vector>
#include <string>
#include <string.h>

namespace energyPlot{
    void generateEnergyPlot(Gnuplot &plotEnergy, std::vector<double> &energyValues, int simNum, int seed){

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

        std::string plotName= "set output ";

        if(simNum == 0){
            plotName += "'Before-'";
        } 
        else{
            plotName += "'After-'";
        }


        if(seed == 10){
            plotName += "'seed10-'";
        }
        else if(seed == 11){
            plotName += "'seed11-'";
        }
        else if(seed == 12){
            plotName += "'seed12-'";
        }
        else if(seed == 13){
            plotName += "'seed13-'";
        }
        else if(seed == 14){
            plotName += "'seed14-'";
        }
        else if(seed == 15){
            plotName += "'seed15-'";
        }
        else if(seed == 16){
            plotName += "'seed16-'";
        }
        else if(seed == 17){
            plotName += "'seed17-'";
        }
        else if(seed == 18){
            plotName += "'seed18-'";
        }
        else if(seed == 19){
            plotName += "'seed19-'";
        }

        plotName += "'EnergyPlot.png'\n";

        

        plotEnergy << "set terminal png\n";
        //plotEnergy << "set output 'EnergyPlot.png'\n";
        plotEnergy << plotName;
        plotEnergy << "replot\n";

        plotEnergy.flush();
        usleep(10000);
    }
}