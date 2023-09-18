#include "../include/GradientDescent.hpp"
#include "../include/Graph.hpp"

// My headers
#include "include/Plot.hpp"
#include "include/Dynamics.hpp"
#include "include/Force.hpp"

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string.h>
#include <string>

double plotHistogramFirst(const std::vector<double>& values)
{
    double bins = 10;
    // Find the minimum and maximum values in the vector
    double minValue = *std::min_element(values.begin(), values.end());
    double maxValue = *std::max_element(values.begin(), values.end());

    // Calculate the histogram
    std::vector<int> histogram(bins, 0);
    float binSize = static_cast<float>(maxValue - minValue + .1) / bins;

    for (double value : values) {
        int binIndex = static_cast<int>((value - minValue) / binSize);
        histogram[binIndex]++;
    }

    // Find the maximum frequency in the histogram
    int maxFrequency = *std::max_element(histogram.begin(), histogram.end());

    int maxVal = histogram[0];
    int maxIndex = 0;
    for(int i = 0; i < bins; i++){
        if(histogram[i] > maxVal){
            maxIndex = i;
        }
    }
    double avg = 0;
    double count = 0;
    for (double value : values) {
        int binIndex = static_cast<int>((value - minValue) / binSize);
        if(binIndex == maxIndex){
            avg += value;
            count++;
        }
    }
    avg = avg / count;


    // Create a histogram visualization image
    int histWidth = 512*2, histHeight = 400*2;
    int histWidth2 = histWidth*1.25,histHeight2 = histHeight *1.25;
    cv::Mat histImage(histHeight2, histWidth2, CV_8UC3, cv::Scalar(255, 255, 255));

    // Draw the histogram bars
    int binWidth = cvRound(static_cast<double>(histWidth) / bins);
    for (int i = 0; i < bins; i++) {
        int barHeight = cvRound(static_cast<double>(histogram[i]) / maxFrequency * histHeight);
        cv::rectangle(histImage, cv::Point(i * binWidth + ((histWidth2-histWidth)/2), histHeight + ((histHeight2-histHeight)/2)), cv::Point((i + 1) * binWidth + ((histWidth2-histWidth)/2), histHeight + ((histHeight2-histHeight)/2) - barHeight), cv::Scalar(0, 0, 0), -1);

        // Add bin values to the x-axis of the graph
        std::ostringstream binValue;
        binValue << std::fixed << std::setprecision(0) << minValue + binSize * i;
        cv::putText(histImage, binValue.str(), cv::Point(i * binWidth + ((histWidth2-histWidth)/2), histHeight2 - 10), cv::FONT_HERSHEY_SIMPLEX, 0.2, cv::Scalar(0, 0, 0), 1, cv::LINE_AA);
    }

    // Add a scale to the y-axis to show bin height
    for (int i = 0; i <= 10; i += 2) {
        int y = histHeight2 - cvRound(static_cast<double>(i) / 10 * histHeight) - ((histHeight2-histHeight)/2);
        cv::putText(histImage, std::to_string(i * maxFrequency / 10), cv::Point(5, y), cv::FONT_HERSHEY_SIMPLEX, 0.4, cv::Scalar(0, 0, 0), 1, cv::LINE_AA);
        // Draw the dotted lines
        cv::Point pt1(30, y);
        cv::Point pt2(histWidth2 - 30, y);
        cv::line(histImage, pt1, pt2, cv::Scalar(200, 200, 200), 1, cv::LINE_4);
    }
/*
    // Create a window and display the histogram
    cv::namedWindow("Histogram", cv::WINDOW_NORMAL);
    cv::imshow("Histogram", histImage);
    std::__cxx11::basic_string<char> name = "tests/Iteration_scaled2"+std::to_string(itNum) + ".jpg";
    cv::imwrite(name, histImage);  
    cv::waitKey(0);
    cv::destroyAllWindows();*/
    //return maxFrequency;

    return avg;
}


std::shared_ptr<Graph> createGraphFromTxt(std::shared_ptr<Graph> graphInit, int gridSize){

    std::vector<double> ev;
    std::string line;
    std::ifstream myFile("LapResults.txt");
    while(getline(myFile, line))
    {
        std::istringstream lineStream(line);
        double first;
        lineStream >> first;
        ev.push_back(first);
    }
    int counter = 0;
    for (int i{0}; i < gridSize*gridSize; i++)
    {
        for (int j = 0; j < gridSize*gridSize; j++)
        {
            if(i == j){
                graphInit->laplacianMatrix(i,j) = ev[counter];
                graphInit->degreeMatrix(i,j) = ev[counter];
            }
            else{
                graphInit->laplacianMatrix(i,j) = ev[counter];
                graphInit->adjacencyMatrix(i,j) = -ev[counter];
            }
            counter++;
        }
    }
    return graphInit;
}


void oneRun(bool resultsGiven, int gridSize, double alpha, bool weightConstraint, bool runManyTimes, int seedNum, int simNum){
    std::shared_ptr<Graph> graphInit = std::make_shared<Graph>();
    graphInit->constructSimpleGraph(gridSize);
    std::vector<double> ev;

    if(resultsGiven){
       graphInit = createGraphFromTxt(graphInit, gridSize);
    }

    graphInit->eigenDecompose();
    GradientDescent gdObj(graphInit,weightConstraint);
    //gdObj.plotHistogram();

    double sum = 0;
    for(int i{0}; i < graphInit->eigenValues.size(); i++){
        ev.push_back(graphInit->eigenValues[i]);
        sum += graphInit->eigenValues[i];
    }
    double mean = sum / graphInit->eigenValues.size();

    int MAX_X = gridSize;
    int MAX_Y = gridSize;
    int PLOT_SCALE = 40;
    int vPad = 2;
    int hPad = 2;
    double damping{0.001}, stiffness{1}, epsilon{0.1};
    double amp{1};
    
    //double freq = sqrt(mean);
    double freq = sqrt(plotHistogramFirst(ev));
    std::cout << freq;
 
    // Generate plot
    Plot my_plot("State Plot - Chosen EigVal: " + std::to_string(freq), PLOT_SCALE, vPad, hPad, MAX_X, MAX_Y);
    my_plot.plotGraph(*graphInit);
    if(runManyTimes){
        my_plot.displayPlot(false);
    }
    else{
        my_plot.displayPlot(true);
    }
         
    Force my_force(amp, freq, graphInit->nodes.size(), alpha, seedNum);
    my_force.insertForceElement(1);
     
    // Simulate dynamics
    int simulationTime{5000};
    int simulationSteps{simulationTime * 100};
    Dynamics my_sim(simulationTime, simulationSteps, damping, stiffness, epsilon, simNum, seedNum);
    my_sim.runCentralizedDynamics(*graphInit, my_force, my_plot);
    if(runManyTimes){
        my_plot.displayPlot(false);
    }
    else{
        my_plot.displayPlot(true);
    }
}

void runNTimes(int runs, int gridSize, double alpha, bool weightConstraint){
    int counter = 0;
    for(int i = 1; i < runs*2 + 1; i++){
        oneRun((i-1)%2, gridSize, alpha, weightConstraint, true, 10+counter, i);
        if(i % 2 == 0){
            counter++;
        }
    }
}

int main(int argc, char *argv[]){
    int gridSize;
    double alpha;
    bool resultsGiven;
    bool weightConstraint{true};

    if(argc==3){
        gridSize = std::stoi(argv[1]);
        alpha = std::stof(argv[2]);
        runNTimes(1, gridSize, alpha, weightConstraint);
    }
    else if(argc==4){
        gridSize = std::stoi(argv[1]);
        resultsGiven = std::stoi(argv[2]);
        alpha = std::stof(argv[3]);
        oneRun(resultsGiven, gridSize, alpha, weightConstraint, false, 15, 0);
    }
    else{
        return 0;
    }
}
