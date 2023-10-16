#include "TwoNormPlot.hpp"
#include <numeric>
#include <vector>
#include <string>
#include <string.h>

#include "gdMain.hpp"


namespace twoNormPlot{
    void generateTwoNormPlot(Gnuplot &plotTwoNorm, std::vector<double> &twoNormValues){

        double minVal = *std::min_element(twoNormValues.begin(), twoNormValues.end());
        double maxVal = *std::max_element(twoNormValues.begin(), twoNormValues.end());

        plotTwoNorm << "reset\n";
        plotTwoNorm << "min=" << minVal << "\n";
        plotTwoNorm << "max=" << maxVal << "\n";

        plotTwoNorm << "set yrange [min:max]\n";
        plotTwoNorm << "set xrange [0:" << twoNormValues.size()+1 << "]\n";
        plotTwoNorm << "set xlabel 'Iteration'\n";
        plotTwoNorm << "set ylabel 'Two Norm Value'\n";

        plotTwoNorm << "plot '-' u 1:2 with lines lc rgb'red' notitle\n"; 

        for (size_t i = 0; i < twoNormValues.size(); ++i) {
            plotTwoNorm << i << " " << twoNormValues[i] << "\n";
        }

        plotTwoNorm << "e\n";

        std::string plotName;

        if(beforeGrad){
            plotName="set output 'Before-" + std::to_string(simNum) +  "-TwoNorm-Plot.png'\n";
        } 
        else{
            plotName= "set output 'After-" + std::to_string(simNum) +  "-TwoNorm-Plot.png'\n";
        }

        plotTwoNorm << "set terminal png\n";
        //plotTwoNorm << "set output 'TwoNormPlot.png'\n";
        plotTwoNorm << plotName;
        plotTwoNorm << "replot\n";

        plotTwoNorm.flush();
        usleep(10000);
    }
}