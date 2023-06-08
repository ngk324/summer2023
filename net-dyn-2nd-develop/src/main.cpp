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


void plotHistogram(const std::vector<float>& values, int bins = 10)
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
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXf> solver(matrix);
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

int findModeIndex(const std::vector<float>& nums) {
    std::unordered_map<int, int> freq;
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
    std::cout << "FreqCount: " << maxFreq << "\n";
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
    int xGrid{50}, yGrid{50};
    Graph my_graph;
    my_graph.constructSimpleGraph(xGrid, yGrid);
    my_graph.computeMatrices();

    int MAX_X = 50;
    int MAX_Y = 50;
    int PLOT_SCALE = 70;
    int vPad = 2;
    int hPad = 2;

    double damping{0.1}, stiffness{5}, epsilon{0.01};

    double amp{1};
    std::vector<int> modes;
    Eigen::MatrixXf B_Matrix = stiffness * (my_graph.laplacianMatrix + epsilon * Eigen::MatrixXf::Identity(my_graph.nodes.size(), my_graph.nodes.size()));
    auto eigen_pairs = get_eigen_pairs(B_Matrix);
    float chosenEigVal{0};
    Eigen::VectorXf chosenEigVec;

    std::cout << "end eign comp." << std::endl;


    std::vector<float> ev;
    for (int i{0}; i < eigen_pairs.size(); i++)
    {
        ev.push_back(eigen_pairs[i].first);
        //std::cout << i << ": " << eigen_pairs[i].first << "\n";
    }

    int modeIndex = findModeIndex(ev);

    if (modeIndex != -1) {
        std::cout << "Mode Index: " << modeIndex << " Mode Value: " << eigen_pairs[modeIndex].first << std::endl;
    } else {
        std::cout << "No mode found." << std::endl;
    }

    chosenEigVal = eigen_pairs[modes.front()].first;

    plotHistogram(ev);
    plotHeatMap(eigen_pairs[modeIndex].second, MAX_X);
    //std::cout << eigen_pairs[0].second << "\n";
    
    /*
    // Generate plot
    Plot my_plot("State Plot - Chosen EigVal: " + std::to_string(chosenEigVal), PLOT_SCALE, vPad, hPad, MAX_X, MAX_Y);
    my_plot.plotGraph(my_graph);
    my_plot.displayPlot(true);
    

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