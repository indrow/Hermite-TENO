#pragma once
#include <Eigen/Dense>
#include <iostream>

using Eigen::Matrix;
using Eigen::Dynamic;
using std::cout;

typedef Matrix<double, 5, 1> Vector5d;

namespace Schemes {
	class Alternative
	{
	public:
		Alternative();
		double reCalcConservative(const Vector5d& q);
		double reCalcFlux(const Vector5d& flux);

	private:
		Vector5d cmcu, cmcf;
	};

}