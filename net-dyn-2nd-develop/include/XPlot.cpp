#include "XPlot.hpp"
#include <numeric>
#include <vector>

namespace XPlot{
    void generateXPlot(Gnuplot &plotX, std::vector<double> &XValues){

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

        plotX << "set terminal png\n";
        plotX << "set output 'XPlot.png'\n";
        plotX << "replot\n";
        //plotX << "set output\n";
        

        plotX.flush();
        usleep(10000);

    }
}