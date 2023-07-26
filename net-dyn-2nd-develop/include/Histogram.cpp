#include "Histogram.hpp"
#include <numeric>

namespace histogram{
    void generateHistogram(Gnuplot &plotObj, Eigen::VectorXf &eigenvalues){
        // std::cout << eigenvalues << std::endl;
        double eigen_sum = std::accumulate(eigenvalues.begin(), eigenvalues.end(), 0);

        plotObj << "reset\n";
        plotObj << "n=50\n";
        // plotObj << "max=" << eigen_sum/20 << "\n";
        plotObj << "max=500\n"; 
        plotObj << "min=0\n";
        plotObj << "width=(max-min)/n\n";

        plotObj << "hist(x,width)=width*floor(x/width)+width/2.0\n";
        // plotObj << "set term png\n";
        // plotObj << "set output 'histogram.png'\n";
        plotObj << "set xrange [min:max]\n";
        plotObj << "set yrange [0:" << eigenvalues.size()+5 << "]\n";

        plotObj << "set offset graph 0.05,0.05,0.05,0.0\n";
        plotObj << "set xtics min,(max-min)/5,max\n";
        plotObj << "set boxwidth width*0.9\n";
        plotObj << "set style fill solid 0.5\n";
        plotObj << "set tics out nomirror\n";
        plotObj << "set xlabel 'x'\n";
        plotObj << "set ylabel 'Frequency'\n";


        plotObj << "plot '-' u (hist($1,width)):(1.0) smooth freq w boxes lc rgb'green' notitle\n"; 

        for (size_t i = 0; i < eigenvalues.size(); ++i) {
            plotObj << eigenvalues[i] << "\n";
        }

        plotObj << "e\n";

        plotObj.flush();
        usleep(10000);
    }
}
