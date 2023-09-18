#ifndef _FORCE_H_
#define _FORCE_H_

#include <vector>
#include <Eigen/Dense>

class Force
{
public:
    double amplitude;
    double frequency;
    int size;
    double alpha;
    std::vector<int> nonZeroElements;
    int seedNum;

    Force(double amp, double freq, int size, double alpha, int seedNum);

    void insertForceElement(int idx);

    Eigen::VectorXf sinusoidalForce(double t);
    Eigen::VectorXf sinCauchyForce(double t);
    Eigen::VectorXf sinCauchyForceNew(double t);
    double inverse_of_normal_cdf(const double p, const double mu, const double sigma);

};

#endif