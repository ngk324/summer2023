#include "Plot.hpp"
#include <algorithm>

CvScalar Plot::defaultNodeColor = cvScalar(150, 150, 150);
CvScalar Plot::defaultEdgeColor = cvScalar(150, 150, 150);

int Plot::defaultNodeSize = 1;
double Plot::defaultEdgeThickness = 1;

Plot::Plot(std::string name, int scale, int vPadding, int hPadding, int xMax, int yMax)
    : windowName{name}, plotScale{scale}, verticalPadding{vPadding}, horizontalPadding{hPadding}, maxX{xMax}, maxY{yMax}
{
    defaultNodeSize = plotScale * 0.2;
    defaultEdgeThickness = plotScale * 0.05;
    //blankImg = cv::Mat(cv::Size(maxX + 2 * verticalPadding - 1, maxY + 2 * horizontalPadding - 1), CV_8UC3, cv::Scalar(255, 255, 255));
    blankImg = cv::Mat(cv::Size(maxX + 2 * verticalPadding - 1, maxY + 2 * horizontalPadding - 1), CV_8UC3, cv::Scalar(0, 165, 255));
    cv::resize(blankImg, blankImg, cv::Size(), plotScale, plotScale);
    currentImg = blankImg.clone();
    initWindow();
    displayPlot();
}

cv::Point Plot::transformGraphToPlot(const Node &n) const
{
    return (cvPoint(plotScale * (n.x + verticalPadding), plotScale * (n.y + horizontalPadding)));
}

void Plot::displayTime(const std::string &time_str)
{
    cv::Size textSize = cv::getTextSize(time_str, cv::FONT_HERSHEY_SIMPLEX, 1, 1, 0);
    cv::Point rightTop(currentImg.size().width - textSize.width - 10, 30);
    rectangle(currentImg, rightTop, cv::Point(textSize.width, -textSize.height) + rightTop,
              cv::Scalar(0, 0, 255), cv::FILLED, cv::LINE_AA);
    putText(currentImg, time_str,
            rightTop, cv::FONT_HERSHEY_SIMPLEX, 1, CV_RGB(255, 255, 255), 1, cv::LINE_AA);
}
void Plot::displayMethod(const std::string &method_str)
{
    cv::Size textSize = cv::getTextSize(method_str, cv::FONT_HERSHEY_SIMPLEX, 1, 1, 0);
    cv::Point leftTop(10, 30);
    rectangle(currentImg, leftTop, cv::Point(textSize.width, -textSize.height) + leftTop,
              cv::Scalar(0, 0, 255), cv::FILLED, cv::LINE_AA);
    putText(currentImg, method_str,
            leftTop, cv::FONT_HERSHEY_SIMPLEX, 1, CV_RGB(255, 255, 255), 1, cv::LINE_AA);
}

void Plot::displayState(const Node &n)
{
    cv::Point plot_pt = transformGraphToPlot(n) + cv::Point(15, -20);
    std::string state_str = std::to_string(n.z).substr(0, 5);
    cv::Size textSize = cv::getTextSize(state_str, cv::FONT_HERSHEY_SIMPLEX, 0.3, 1, 0);
    rectangle(currentImg, plot_pt, cv::Point(textSize.width, -textSize.height) + plot_pt,
              cv::Scalar(255, 255, 255), cv::FILLED, cv::LINE_AA);
    putText(currentImg, state_str,
            plot_pt, cv::FONT_HERSHEY_SIMPLEX, 0.3, CV_RGB(0, 0, 0), 0.5, cv::LINE_AA);
}

void Plot::displayEigenvectorVal(const Node &n, float val)
{
    cv::Point plot_pt = transformGraphToPlot(n) + cv::Point(15, -20);
    std::string state_str = std::to_string(val).substr(0, 5);
    cv::Size textSize = cv::getTextSize(state_str, cv::FONT_HERSHEY_SIMPLEX, 0.3, 1, 0);
    rectangle(currentImg, plot_pt, cv::Point(textSize.width, -textSize.height) + plot_pt,
              cv::Scalar(255, 255, 255), cv::FILLED, cv::LINE_AA);
    putText(currentImg, state_str,
            plot_pt, cv::FONT_HERSHEY_SIMPLEX, 0.3, CV_RGB(0, 0, 0), 0.5, cv::LINE_AA);
}

void Plot::plotEigen(const Graph &g, Eigen::VectorXf &vec, int size)
{
    CvScalar color;
    for (int i{0}; i < vec.size(); i++)
    {
        double val = abs(vec(i)) * 500;
        if (vec(i) > 0)
            color = cvScalar(val, 0, 0);
        else
            color = cvScalar(0, 0, val);
        cv::circle(currentImg, transformGraphToPlot(*g.nodes[i]), size, color, -1);
        displayEigenvectorVal(*g.nodes[i], vec(i));
    }
}

void Plot::plotNode(Node &n, double max, double min, int size)
{
    CvScalar color;
    double val;
    if (n.z > 0){
        val = abs(n.z) * (255.0/abs(max));
        color = cvScalar(val, 0, 0);
    }
    else{
        val = abs(n.z) * (255.0/abs(min));
        color = cvScalar(0, 0, val);
    }
    //cv::circle(currentImg, transformGraphToPlot(n), size, cvScalar(0, 0, 0), -1);
    cv::circle(currentImg, transformGraphToPlot(n), size, color, -1);
}

//void Plot::plotEdge(Node &n1, Node &n2, CvScalar color, double thickness)
void Plot::plotEdge(Node &n1, Node &n2,  double thickness, CvScalar color)
{
    cv::line(currentImg, transformGraphToPlot(n1), transformGraphToPlot(n2), color, thickness);
}


double Plot::round(double var, int place)
{
    double value = (int)(var * 10000 + .5);
    return (double)value / 10000;
}

void Plot::plotGraph(Graph &g)
{   
    double max = std::max<double>(g.adjacencyMatrix.maxCoeff(),0.0);
    double min = std::min<double>(g.adjacencyMatrix.minCoeff(),0.0);

    double maxEV = std::max<double>(g.eigenValues.maxCoeff(),0.0);
    double minEV = std::min<double>(g.eigenValues.minCoeff(),0.0);
    //std::cout << "edges plotted \n";
    for (int i{0}; i < g.nodes.size(); i++)
    {
        for (int j{0}; j < g.nodes[i]->neighbors.size(); j++)
        {
            double edgeWeight = g.adjacencyMatrix(g.nodes[i]->id,g.nodes[i]->neighbors[j]->id);
            CvScalar edgeColor = cvScalar(0,0,0);
            edgeColor.val[0] = 255-(((int)(g.adjacencyMatrix.coeff((*g.nodes[i]).id,(*g.nodes[i]->neighbors[j]).id)*255))/(max-min));
            edgeColor.val[1] = edgeColor.val[0];
            edgeColor.val[2] = edgeColor.val[0];
            if((*g.nodes[i]).id < (*g.nodes[i]->neighbors[j]).id){
                plotEdge(*g.nodes[i], *g.nodes[i]->neighbors[j], 75/sqrt(g.nodes.size()), edgeColor);
            }
        }
    }

    for (int i{0}; i < g.nodes.size(); i++)
    {
        plotNode(*g.nodes[i], maxEV, minEV);
    }
}

void Plot::plotEdgeCircle(double value1, double value2,  double thickness, CvScalar color)
{
    // alpha masking 
    cv::Mat overlay = currentImg.clone();
    cv::Mat original = currentImg.clone();
    cv::line(overlay, transformGraphToPlotCircle(value1), transformGraphToPlotCircle(value2), color, thickness);
    double alpha = 0.5;
    cv::addWeighted(overlay, alpha, original, 1 - alpha, 0, currentImg);

}

cv::Point Plot::transformGraphToPlotCircle(double value) const
{
    double x = (cos(value) + 1.1) * (double) plotScale / 2.0;
    double y = (sin(value) + 1.025) * (double) plotScale / 2.0;

    return (cvPoint(x*12,y*12));
}

void Plot::plotNodeCircle(Node &n, double max, double min, int num, int size)
{
    double value = ((double)num)*((3.14159265358979323846264338327 * 2.0) / 100.0);
    CvScalar color;
    double val;
    if (n.z > 0){
        val = abs(n.z) * (255.0/abs(max));
        color = cvScalar(val, 0, 0);
    }
    else{
        val = abs(n.z) * (255.0/abs(min));
        color = cvScalar(0, 0, val);
    }
    //cv::circle(currentImg, transformGraphToPlot(n), size, cvScalar(0, 0, 0), -1);
    //cv::circle(currentImg, cvPoint(plotScale * (cos(value) * 0.95), plotScale * (sin(value))*0.95) , size, color, -1);
    cv::circle(currentImg, transformGraphToPlotCircle(value), size, color, -1);

    //cv::circle(currentImg, cvPoint(plotScale * (cos(value) + verticalPadding), plotScale * (sin(value) + horizontalPadding)) , size, color, -1);
}


void Plot::plotGraphCircle(Graph &g)
{   
    double max = std::max<double>(g.adjacencyMatrix.maxCoeff(),0.0);
    double min = std::min<double>(g.adjacencyMatrix.minCoeff(),0.0);

    double maxEV = std::max<double>(g.eigenValues.maxCoeff(),0.0);
    double minEV = std::min<double>(g.eigenValues.minCoeff(),0.0);
    //std::cout << "edges plotted \n";
    
    for (int i{0}; i < g.nodes.size()-1; i++)
    {
        for (int j{i+1}; j < g.nodes.size(); j++)
        {
            double edgeWeight = g.adjacencyMatrix(g.nodes[i]->id,g.nodes[j]->id);
            CvScalar edgeColor = cvScalar(0,0,0);
            edgeColor.val[0] = 255-(((int)(g.adjacencyMatrix.coeff((*g.nodes[i]).id,g.nodes[j]->id)*255))/(max-min));
            edgeColor.val[1] = edgeColor.val[0];
            edgeColor.val[2] = edgeColor.val[0];
            //if((*g.nodes[i]).id < (*g.nodes[i]->neighbors[j]).id){
                plotEdgeCircle(g.nodes[i]->id, g.nodes[j]->id, 75/sqrt(g.nodes.size()), edgeColor);
            //}
        }
    }

    for (int i{0}; i < g.nodes.size(); i++)
    {
        plotNodeCircle(*g.nodes[i], maxEV, minEV,i);
    }
}


void Plot::initWindow()
{
    cv::namedWindow(windowName, cv::WINDOW_AUTOSIZE);
    cv::moveWindow(windowName, 50, 50);
}

void Plot::displayPlot(bool waitForKey)
{
    cv::imshow(windowName, currentImg);
    if (waitForKey){
        char key = (char)cv::waitKey(); // explicit cast
        while (key != 27)
        {
            key = (char)cv::waitKey();
        }
    }
    else
        cv::waitKey(1);
}
