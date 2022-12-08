#pragma once
#include <cmath>
#include <iostream>
#include <Eigen/Dense>

using std::cout;
using std::endl;
using std::fmax;

namespace Equations {
namespace Scalar1D {
class Burgers
{
public:
	double f(const double& u);
	double df(const double& u);
	double Rusanov(const double& ul, const double& ur);
	double dRusanov(const double& ul, const double& ur,
			const double& dul, const double& dur);
    double LxF(const double& ul, const double& ur, const double& cmax);
    double dLxF(const double& ul, const double& ur,
                const double& dul, const double& dur, const double& cmax);
};
}
}

