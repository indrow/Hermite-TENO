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
namespace Scalar2D {
class LinearAdvection {
public:
	double sine(const double& x, const double& y, const double& dx, const double& dy,
			const Vector9d& weights, const Vector9d& absis);
	double sinedx(const double& x, const double& y, const double& dx, const double& dy,
			const Vector9d& weights, const Vector9d& absis);
	double sinedy(const double& x, const double& y, const double& dx, const double& dy,
			const Vector9d& weights, const Vector9d& absis);

	double f(const double& u);
	double df(const double& u);
	double Rusanov(const double& ul, const double& ur);
	double dRusanov(const double& ul, const double& ur,
			const double& dul, const double& dur);

    double LxF(const double& ul, const double& ur, const double& cmax);
    double dLxF(const double& ul, const double& ur,
                    const double& dul, const double& dur, const double& cmax);


private:
	const double pi = 4.0 * atan(1.0);
	inline double fSine(const double& x, const double& y) {
		return 0.5 + sin(pi * (x + y) / 2.0);
	}

	inline double dfSine(const double& x, const double& y) {
		return cos(pi * (x + y) / 2.0) * pi / 2.0;
	}

};
}
}
