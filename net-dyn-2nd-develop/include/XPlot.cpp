#include "XPlot.hpp"
#include <numeric>
#include <vector>

namespace XPlot{
    void generateXPlot(Gnuplot &plotX, std::vector<double> &XValues, int simNum, int seed){

        double minVal = *std::min_element(XValues.begin(), XValues.end());
        double maxVal = *std::max_element(XValues.begin(), XValues.end());

        plotX << "reset\n";
        plotX << "min=" << minVal << "\n";
        plotX << "max=" << maxVal << "\n";

        plotX << "set yrange [min:max]\n";
        plotX << "set xrange [0:" << XValues.size()+1 << "]\n";
        plotX << "set xlabel 'Iteration'\n";
        plotX << "set ylabel 'X'\n";

        plotX << "plot '-' u 1:2 with lines lc rgb'red' notitle\n"; 

        for (size_t i = 0; i < XValues.size(); ++i) {
            plotX << i << " " << XValues[i] << "\n";
        }

        plotX << "e\n";

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

        plotName += "'XPlot.png'\n";

        

        plotX << "set terminal png\n";
        //plotEnergy << "set output 'EnergyPlot.png'\n";
        plotX << plotName;
        plotX << "replot\n";


        plotX.flush();
        usleep(10000);

    }
}