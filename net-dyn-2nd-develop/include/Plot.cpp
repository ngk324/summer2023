#include "Plot.hpp"

CvScalar Plot::defaultNodeColor = cvScalar(150, 150, 150);
CvScalar Plot::defaultEdgeColor = cvScalar(150, 150, 150);

int Plot::defaultNodeSize = 1;
double Plot::defaultEdgeThickness = 1;

Plot::Plot(std::string name, int scale, int vPadding, int hPadding, int xMax, int yMax)
    : windowName{name}, plotScale{scale}, verticalPadding{vPadding}, horizontalPadding{hPadding}, maxX{xMax}, maxY{yMax}
{
    defaultNodeSize = plotScale * 0.2;
    defaultEdgeThickness = plotScale * 0.05;
    blankImg = cv::Mat(cv::Size(maxX + 2 * verticalPadding - 1, maxY + 2 * horizontalPadding - 1), CV_8UC3, cv::Scalar(255, 255, 255));
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

void Plot::plotNode(Node &n, int size)
{
    CvScalar color;
    double val = abs(n.z) * 500;
    if (n.z > 0)
        color = cvScalar(val, 0, 0);
    else
        color = cvScalar(0, 0, val);
    cv::circle(currentImg, transformGraphToPlot(n), size, color, -1);
}

void Plot::plotEdge(Node &n1, Node &n2, CvScalar color, double thickness)
{
    cv::line(currentImg, transformGraphToPlot(n1), transformGraphToPlot(n2), color, thickness);
}

void Plot::plotGraph(Graph &g)
{
    for (int i{0}; i < g.nodes.size(); i++)
    {
        for (int j{0}; j < g.nodes[i]->neighbors.size(); j++)
        {
         //   plotEdge(*g.nodes[i], *g.nodes[i]->neighbors[j]);
        }
    }

    for (int i{0}; i < g.nodes.size(); i++)
    {
        plotNode(*g.nodes[i]);
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
    if (waitForKey)
        cv::waitKey(0);
    else
        cv::waitKey(1);
}