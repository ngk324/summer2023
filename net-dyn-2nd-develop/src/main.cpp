// OpenCV libraries
#include <opencv2/opencv.hpp>
#include <opencv2/highgui.hpp>

#include <Eigen/Eigenvalues>

// My headers
#include "include/Plot.hpp"
#include "include/Dynamics.hpp"
#include "include/Force.hpp"

#include <opencv2/opencv.hpp>

#include <cmath>

#include <math.h>

#include <numeric>

/*double compute_var_grad(const std::vector<double>& values, double mean, Eigen::Matrix<float, -1, 1> vki, Eigen::Matrix<float, -1, 1> vkj){
    double d_var = 0; 
    for(int i = 0; i < values.size(); i++){
        d_var = d_var + (2.0/values.size()) * (values[i]-mean)*(vki[i] * vki[i] + vkj[i]*vkj[i] - 2.0*vki[i]*vkj[i]);
    }
    return d_var;
}*/

/*double compute_var_grad(const std::vector<double>& values, double mean, Eigen::Matrix<float, -1, 1> vki, Eigen::Matrix<float, -1, 1> vkj){
    double d_var = 0; 
    for(int i = 0; i < values.size(); i++){
        d_var = d_var + (2.0/values.size()) * (values[i]-mean)*(vki[i] * vki[i] + vkj[i]*vkj[i] - 2.0*vki[i]*vkj[i]);
    }
    return d_var;
}*/

double compute_var_grad(const std::vector<double>& values, Eigen::Matrix<float, -1, 1> vki, Eigen::Matrix<float, -1, 1> vkj, double stdev){
    double d_var = 0; 

    //double h = 1.09;

    //double h = 0.54975 * stdev;

    double h = stdev;

    //(arctan((stdev)/x) / M_PI) - (arctan((-stdev)/x) / M_PI) = .68;

    

    for(int i = 1; i < values.size(); i++){
        int index = i;
        double maxVal;
        int maxIndex;
        if(index != 1){
            //maxVal = ((4*h*h)/values[1]) * (1/(pow(pow(sqrt(values[1])-sqrt(values[index]),2) + h * h, 2)));   
            maxVal = ((4*h*h)/values[index]) * (1/(pow(pow(sqrt(values[index])-sqrt(values[1]),2) + h * h, 2)));
            maxIndex = 1;
            //std::cout << " \nTest i,j " << i << " , " << 0 << " " << (a / (2 * values[1])) * exp(-b * (pow(values[1] - values[index],2)));
            
        }
        else{
            //maxVal = ((4*h*h)/values[2]) * (1/(pow(pow(sqrt(values[2])-sqrt(values[index]),2) + h * h, 2)));   
            maxVal = ((4*h*h)/values[index]) * (1/(pow(pow(sqrt(values[index])-sqrt(values[2]),2) + h * h, 2)));
            maxIndex = 2;
            //std::cout << " \nTest i,j " << i << " , " << 0 << " " << (a / (2 * values[1])) * exp(-b * (pow(values[1] - values[index],2)));
            
        }
        for(int j = 0; j < values.size(); j++){
            if(j != index && j != maxIndex){
                //double testVal = ((4*h*h)/values[j]) * (1/(pow(pow(sqrt(values[j])-sqrt(values[index]),2) + h * h, 2))); 
                double testVal = ((4*h*h)/values[index]) * (1/(pow(pow(sqrt(values[index])-sqrt(values[j]),2) + h * h, 2)));
                if(testVal > maxVal){
                    maxVal = testVal;
                    maxIndex = j;
                }
            }
            //std::cout << " \nTest i,j " << i << " , " << j << " " << (a / (2 * values[j])) * exp(-b * (pow(values[j] - values[index],2)));
        }

        //std::cout << "Val " << maxVal << " Index " << maxIndex << "\n ";

        //double min_func_val =  -((8 * h * h) / (values[maxIndex] * sqrt(values[index]))) * ((sqrt(values[index])-sqrt(values[maxIndex])) / pow(pow(sqrt(values[index]-sqrt(values[maxIndex])),2) + h * h,3));
        //double min_func_val =  ((8 * h * h) / (values[maxIndex] * sqrt(values[index]))) * ((sqrt(values[maxIndex])-sqrt(values[index])) / pow(pow(sqrt(values[maxIndex])-sqrt(values[index]),2) + h * h,3));

        //double min_func_val = ((4 * h*h)/(sqrt(values[index])*pow(values[index],1.5))) * ((-3 * values[index] + 4 * sqrt(values[maxIndex])*sqrt(values[index])-values[maxIndex]-h*h)/(pow(h*h+pow(-sqrt(values[index])+sqrt(values[maxIndex]),2),3)));
        double min_func_val = -((4 * h*h)/(sqrt(values[index])*pow(values[index],1.5))) * ((3 * values[index] - 4 * sqrt(values[maxIndex])*sqrt(values[index])+values[maxIndex]+h*h)/(pow(h*h+pow(sqrt(values[index])-sqrt(values[maxIndex]),2),3)));
        
        d_var = d_var + min_func_val * (vki[i] * vki[i] + vkj[i]*vkj[i] - 2.0*vki[i]*vkj[i]);
    }


    /*for(int i = 0; i < values.size(); i++){
        d_var = d_var + (2.0/values.size()) * (values[i]-mean)*(vki[i] * vki[i] + vkj[i]*vkj[i] - 2.0*vki[i]*vkj[i]);
    }*/

    return d_var;
}

double round(double var, int place)
{
    double value = (int)(var * 100 + .5);
    return (double)value / 100;
}

void calc_steady_state(const std::vector<double>& values, double freq){
    double sum = 0;
    int counter = 0;
    for(int i =0; i < values.size(); i++){
        if(values[i] == freq){
            //sum = sum + 1000000;
            counter++;
        }
        else{
            sum = sum + 1.0/(values[i] - freq);
        }
    }
    std::cout << "\nSteady state solution: " << sum << " counter: " << counter << "\n";
}


void plotHeatMap(Eigen::Matrix<float, -1, 1> arr, int size)
{
    // Define dimensions of the heatmap
    int width = 400;
    int height = 400;

    // Create a 2D array of values
    double heatmap[width][height];

    // Populate the heatmap array with values (example values shown)
    int index = 0;
    for(int k = 0; k < size * size; k++){
        for (int i = 0; i < width / size; i++){
            for (int j = 0; j < height / size; j++){
                // Calculate value based on position (example calculation shown)
                //heatmap[i][j] = static_cast<double>(arr.row(count));
                heatmap[(k / size) * (width / size) + i][(k%size) * (width / size) + j] = arr[index];
            }
        }
        index++;
    }
    // Find the minimum and maximum values in the heatmap array
    double minValue = heatmap[0][0];
    double maxValue = heatmap[0][0];
    for (int i = 0; i < width; ++i)
    {
        for (int j = 0; j < height; ++j)
        {
            if (heatmap[i][j] < minValue)
                minValue = heatmap[i][j];
            if (heatmap[i][j] > maxValue)
                maxValue = heatmap[i][j];
        }
    }

    // Create a Mat object to store the heatmap image
    cv::Mat image(height, width, CV_8UC3);

    // Iterate through each pixel and set the color based on the corresponding value in the array
    for (int i = 0; i < width; ++i)
    {
        for (int j = 0; j < height; ++j)
        {
            // Normalize the value to the range [0, 255]
            int intensity = static_cast<int>((heatmap[i][j] - minValue) / (maxValue - minValue) * 255);

            // Set the color of the current pixel based on the intensity
            cv::Vec3b color(intensity, intensity, intensity);

            // Set the color of the current pixel
            image.at<cv::Vec3b>(j, i) = color;
        }
    }

    // Display the heatmap
    cv::imshow("Heatmap", image);

    // Wait for a key press
    cv::waitKey(0);

    // Close the window
    cv::destroyAllWindows();
}


void plotHistogram(const std::vector<double>& values, double bins = 8.0)
{
    // Find the minimum and maximum values in the vector
    int minValue = *std::min_element(values.begin(), values.end());
    int maxValue = *std::max_element(values.begin(), values.end());

    // Calculate the histogram
    std::vector<int> histogram(bins, 0);
    float binSize = static_cast<float>(maxValue - minValue + 1) / bins;

    for (int value : values) {
        int binIndex = static_cast<int>((value - minValue) / binSize);
        histogram[binIndex]++;
    }

    // Find the maximum frequency in the histogram
    int maxFrequency = *std::max_element(histogram.begin(), histogram.end());

    // Create a histogram visualization image
    int histWidth = 512, histHeight = 400;
    cv::Mat histImage(histHeight, histWidth, CV_8UC3, cv::Scalar(255, 255, 255));

    // Draw the histogram bars
    int binWidth = cvRound(static_cast<double>(histWidth) / bins);
    for (int i = 0; i < bins; i++) {
        int barHeight = cvRound(static_cast<double>(histogram[i]) / maxFrequency * histHeight);
        cv::rectangle(histImage, cv::Point(i * binWidth, histHeight), cv::Point((i + 1) * binWidth, histHeight - barHeight), cv::Scalar(0, 0, 0), -1);
    }

    // Create a window and display the histogram
    cv::namedWindow("Histogram", cv::WINDOW_NORMAL);
    cv::imshow("Histogram", histImage);
    cv::waitKey(0);
    cv::destroyAllWindows();
}


std::vector<std::pair<float, Eigen::VectorXf>> get_eigen_pairs(Eigen::MatrixXf &matrix)
{
    std::vector<std::pair<float, Eigen::VectorXf>> eigen_pairs;
    //Eigen::SelfAdjointEigenSolver<Eigen::MatrixXf> solver(matrix);
    Eigen::EigenSolver<Eigen::MatrixXf> solver(matrix);
    Eigen::VectorXcf eigVals = solver.eigenvalues();
    Eigen::MatrixXcf eigVecs = solver.eigenvectors();

    for (int i{0}; i < eigVals.size(); i++)
    {
        /*if (eigVals(i).imag() != 0)
        {
            std::cout << "Complex eigenvalue!" << std::endl;
            exit(EXIT_FAILURE);
        }*/
        eigen_pairs.push_back(std::make_pair(eigVals(i).real(), eigVecs.col(i).real()));
    }
    std::sort(eigen_pairs.begin(), eigen_pairs.end(), [&](const std::pair<float, Eigen::VectorXf> &a, const std::pair<float, Eigen::VectorXf> &b) -> bool
              { return a.first < b.first; });
    return eigen_pairs;
}

std::vector<std::pair<float, Eigen::VectorXf>> get_eigen_pairs2(Eigen::MatrixXf &matrix)
{
    std::vector<std::pair<float, Eigen::VectorXf>> eigen_pairs;
    //Eigen::SelfAdjointEigenSolver<Eigen::MatrixXf> solver(matrix);
    Eigen::EigenSolver<Eigen::MatrixXf> solver(matrix);
    Eigen::VectorXcf eigVals = solver.eigenvalues();
    Eigen::MatrixXcf eigVecs = solver.eigenvectors();

    //eigVecs = eigVecs.transpose();

    for (int i{0}; i < eigVals.size(); i++)
    {
        /*if (eigVals(i).imag() != 0)
        {
            std::cout << "Complex eigenvalue!" << std::endl;
            exit(EXIT_FAILURE);
        }*/
        eigen_pairs.push_back(std::make_pair(eigVals(i).real(), eigVecs.row(i).real()));
    }
    std::sort(eigen_pairs.begin(), eigen_pairs.end(), [&](const std::pair<float, Eigen::VectorXf> &a, const std::pair<float, Eigen::VectorXf> &b) -> bool
              { return a.first < b.first; });
    return eigen_pairs;
}



int findModeIndex(const std::vector<double>& nums) {
    std::unordered_map<double, double> freq;
    int maxFreq = 0;
    int modeIndex = -1;

    // Count the frequency of each element
    for (int i = 0; i < nums.size(); i++) {
        freq[nums[i]]++;
        if (freq[nums[i]] > maxFreq) {
            maxFreq = freq[nums[i]];
            modeIndex = i;
        }
    }
    std::cout << "\nFreqCount: " << maxFreq << "\n";
    return modeIndex;
}

int main(int argc, char *argv[])
{
    bool decentralizedAlg = false;
    if (argc == 2)
    {
        decentralizedAlg = std::stoi(argv[1]);
    }
    // Initialize and construct simple graph
    int xGrid{20}, yGrid{20};
    Graph my_graph;
    my_graph.constructSimpleGraph(xGrid, yGrid);
    my_graph.computeMatrices();

    int MAX_X = 20;
    int MAX_Y = 20;
    int PLOT_SCALE = 40;
    int vPad = 2;
    int hPad = 2;

    /*
    for(int i = 0; i < MAX_X*MAX_X; i++){
        for(int j = 0; j < MAX_Y*MAX_Y; j++){
            std::cout << i << ", " << j << ": " << my_graph.laplacianMatrix(i,j) << "\n";
        }
    }
    
    std::cout << "\n\n";*/
    double damping{0.1}, stiffness{5}, epsilon{0.01};

    double amp{1};
    std::vector<int> modes;
    //Eigen::MatrixXf B_Matrix = stiffness * (my_graph.laplacianMatrix + epsilon * Eigen::MatrixXf::Identity(my_graph.nodes.size(), my_graph.nodes.size()));
    Eigen::MatrixXf B_Matrix = my_graph.laplacianMatrix;
    auto eigen_pairs = get_eigen_pairs(B_Matrix);
    float chosenEigVal{0};
    Eigen::VectorXf chosenEigVec;

    //std::cout << "end eign comp." << std::endl;

    std::vector<double> ev;
    for (int i{0}; i < eigen_pairs.size(); i++)
    {
        ev.push_back(round(eigen_pairs[i].first,2));
        std::cout << ev[i] << ", ";
    }

    int modeIndex = findModeIndex(ev);

    //calc_steady_state(ev, ev[modeIndex]);

    
    if (modeIndex != -1) {
        std::cout << "\nMode Index: " << modeIndex << " Mode Value: " << ev[modeIndex] << std::endl;
    } else {
        std::cout << "\nNo mode found." << std::endl;
    }

    chosenEigVal = eigen_pairs[modeIndex].first;

    plotHistogram(ev);
    plotHeatMap(eigen_pairs[modeIndex].second, MAX_X);
    
    
    // Generate plot
    Plot my_plot("State Plot - Chosen EigVal: " + std::to_string(chosenEigVal), PLOT_SCALE, vPad, hPad, MAX_X, MAX_Y);
    my_plot.plotGraph(my_graph);
    my_plot.displayPlot(true);

    //for(int f = 0; f < 150.0 / (.1 / sqrt(2*(MAX_X*MAX_X) - 2 * MAX_X)); f++){
    for(int f = 0; f < 225; f++){
        Eigen::MatrixXf var_grad(MAX_X, MAX_Y);
        std::vector<double> var_gradient;

        double mean = std::accumulate(ev.begin(), ev.end(), 0.0) / ev.size();

        std::vector<double> diff(ev.size());
        std::transform(ev.begin(), ev.end(), diff.begin(),
        std::bind2nd(std::minus<double>(), mean));
        double sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
        double stdev = std::sqrt(sq_sum / ev.size());

        //std::cout << "\n" << (((double)f)  / (150.0 / (.1 / sqrt(2*(MAX_X*MAX_X) - 2 * MAX_X)))) * 100 <<  "% Variance : " << stdev * stdev << "\n";

        std::cout << "\n" << "Variance : " << stdev * stdev << "\n";


        auto eigen_pairs3 = get_eigen_pairs2(B_Matrix);


        for(int i = 0; i < my_graph.nodes.size(); i++){
            for(int j = 0; j < (my_graph.nodes[i])->neighbors.size(); j++){
                /*if(my_graph.nodes[i]->isNeighbor(my_graph.nodes[j])){
                    var_grad(i,j) = compute_var_grad(ev, mean, eigen_pairs[(*my_graph.nodes[0]).id].second, eigen_pairs[(*my_graph.nodes[0]->neighbors[1]).id].second);
                */
                if((*my_graph.nodes[i]).id < (*my_graph.nodes[i]->neighbors[j]).id){
                    var_gradient.push_back(compute_var_grad(ev, eigen_pairs3[(*my_graph.nodes[i]).id].second, eigen_pairs3[(*my_graph.nodes[i]->neighbors[j]).id].second, stdev));
                    //std::cout << (*my_graph.nodes[i]).id << ", " << (*my_graph.nodes[i]->neighbors[j]).id << ": " << (compute_var_grad(ev, mean, eigen_pairs3[(*my_graph.nodes[i]).id].second, eigen_pairs3[(*my_graph.nodes[i]->neighbors[j]).id].second)) << "\n";
                    //var_gradient.push_back(compute_var_grad(ev, mean, eigen_pairs[].second, eigen_pairs[].second));
                }
            }
        }

        //std::cout << "\n\n";


        double vec_dot_unit_norm = 0;
        std::vector<double> dot_vec;

        for(int i = 0; i < var_gradient.size(); i++){
            vec_dot_unit_norm = vec_dot_unit_norm + (var_gradient[i] * (1.0/sqrt(var_gradient.size())));
        }

        for(int i = 0; i < var_gradient.size(); i++){
            dot_vec.push_back(var_gradient[i] - (vec_dot_unit_norm * (1.0/sqrt(var_gradient.size()))));
        }

        std::cout << "\nv' vector \n";

        for(int i = 0; i < var_gradient.size(); i++){
           // std::cout << dot_vec[i] << "\n";
        }


        double ep = .001 / sqrt(2*(MAX_X*MAX_X) - 2 * MAX_X);
        //double ep = 1;
        int counter = 0;
        for(int i = 0; i < my_graph.nodes.size(); i++){
            for(int j = 0; j < (my_graph.nodes[i])->neighbors.size(); j++){
                /*if(my_graph.nodes[i]->isNeighbor(my_graph.nodes[j])){
                    var_grad(i,j) = compute_var_grad(ev, mean, eigen_pairs[(*my_graph.nodes[0]).id].second, eigen_pairs[(*my_graph.nodes[0]->neighbors[1]).id].second);
                }*/
                if((*my_graph.nodes[i]).id < (*my_graph.nodes[i]->neighbors[j]).id){
                    my_graph.adjacencyMatrix((*my_graph.nodes[i]).id, (*my_graph.nodes[i]->neighbors[j]).id) = my_graph.adjacencyMatrix((*my_graph.nodes[i]).id, (*my_graph.nodes[i]->neighbors[j]).id) + ep * dot_vec[counter];
                    my_graph.adjacencyMatrix((*my_graph.nodes[i]->neighbors[j]).id, (*my_graph.nodes[i]).id) = my_graph.adjacencyMatrix((*my_graph.nodes[i]->neighbors[j]).id, (*my_graph.nodes[i]).id) + ep * dot_vec[counter];
                    //std::cout << (*my_graph.nodes[i]).id << ", " << (*my_graph.nodes[i]->neighbors[j]).id << ": " << (compute_var_grad(ev, mean, eigen_pairs[(*my_graph.nodes[i]).id].second, eigen_pairs[(*my_graph.nodes[i]->neighbors[j]).id].second)) << "\n";
                    //std::cout << ((*my_graph.nodes[i]).id) << ", " << (*my_graph.nodes[i]->neighbors[j]).id << " - " << dot_vec[counter] << "\n";
                    counter++;
                }
            }
        }

        my_graph.computeMatrices2();
        Eigen::MatrixXf B_Matrix2 = my_graph.laplacianMatrix;
        auto eigen_pairs2 = get_eigen_pairs(B_Matrix2);

        std::vector<double> ev2;
        for (int i{0}; i < eigen_pairs2.size(); i++)
        {
            ev2.push_back(round(eigen_pairs2[i].first,2));
            //std::cout << i << ": " << ev[i] << "\n";
        }
    /*
        std::cout << "Trace vales\n";
        double trace = 0;
        for(int i = 0; i < MAX_X*MAX_X; i++){
            trace = trace + my_graph.laplacianMatrix(i,i);
        }

        std::cout << "Trace" << trace << "\n";
    
    */
        modeIndex = findModeIndex(ev2);

        //plotHistogram(ev2);

        //Plot my_plot2("State Plot - Chosen EigVal: " + std::to_string(chosenEigVal), PLOT_SCALE, vPad, hPad, MAX_X, MAX_Y);
        //my_plot2.plotGraph(my_graph);
        //my_plot2.displayPlot(true);

        //plotHeatMap(eigen_pairs2[modeIndex].second, MAX_X);

        ev = ev2;
        B_Matrix = B_Matrix2;
    }

    double mean = std::accumulate(ev.begin(), ev.end(), 0.0) / ev.size();

        std::vector<double> diff(ev.size());
        std::transform(ev.begin(), ev.end(), diff.begin(),
        std::bind2nd(std::minus<double>(), mean));
        double sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
        double stdev = std::sqrt(sq_sum / ev.size());

    std::cout << "Var: " << stdev * stdev << "\n";
/*
    for(int i = 0; i < MAX_X*MAX_X; i++){
        for(int j = 0; j < MAX_X*MAX_X; j++){
            std::cout << "( " << i << ", " << j << ") :" << my_graph.laplacianMatrix(i,j) << "\n";
        }
    }*/

    plotHistogram(ev);

        Plot my_plot2("State Plot - Chosen EigVal: " + std::to_string(chosenEigVal), PLOT_SCALE, vPad, hPad, MAX_X, MAX_Y);
        my_plot2.plotGraph(my_graph);
        my_plot2.displayPlot(true);

    for(int i = 0; i < ev.size(); i++){
        std::cout << ev[i] << ", ";
    }

    if (modeIndex != -1) {
        std::cout << "\nMode Index: " << modeIndex << " Mode Value: " << ev[modeIndex] << std::endl;
    } else {
        std::cout << "\nNo mode found." << std::endl;
    }

    

        //plotHeatMap(eigen_pairs2[modeIndex].second, MAX_X);
    /*
    double freq{sqrt(chosenEigVal)};
    Force my_force(amp, freq, my_graph.nodes.size());
    my_force.insertForceElement(1);

    // Generate plot
    for (int i{0}; i < modes.size(); i++)
    {
        Plot eig_plot("Eigen Plot " + std::to_string(i + 1) + " - Eigval: " + std::to_string(eigen_pairs[modes[0]].first), PLOT_SCALE, vPad, hPad, MAX_X, MAX_Y);
        eig_plot.plotGraph(my_graph);
        eig_plot.displayPlot(false);
        eig_plot.plotEigen(my_graph, eigen_pairs[modes[0]].second);
        eig_plot.displayPlot(true);
    }
    */
    /*
    // Simulate dynamics
    int simulationTime{400};
    int simulationSteps{simulationTime * 100};
    Dynamics my_sim(simulationTime, simulationSteps, damping, stiffness);
    if (decentralizedAlg)
        my_sim.runDecentralizedDynamics(my_graph.nodes, my_force, my_plot);
    else
        my_sim.runCentralizedDynamics(my_graph, my_force, my_plot);
    my_plot.displayPlot(true);*/

    return 0;
}