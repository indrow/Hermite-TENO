#pragma once
#include <cmath>
#include <iostream>
#include <Eigen/Dense>

using std::cout;
using std::endl;
using std::fmax;
using Eigen::Matrix;
typedef Matrix<double, 9, 1> Vector9d;

namespace Equations {
namespace Scalar1D {
class LinearAdvection
{
public:
	double sine(const double& x, const double& dx, const Vector9d& weights, const Vector9d& absis);
	double sinedx(const double& x, const double& dx, const Vector9d& weights, const Vector9d& absis);

private:
	const double pi = 4.0 * atan(1.0);
	inline double fSine(const double& x) {
		return 0.5 + sin(pi * x);
//		return 0.5 + pow(sin(pi * x), 4);
	}

	inline double dfSine(const double& x) {
		return cos(pi * x) * pi;
//		return 4.0 * pi * pow(sin(pi * x), 3) * cos(pi * x);
	}
};
}
}
