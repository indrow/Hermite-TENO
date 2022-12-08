#include "Equations_Scalar1D_LinearAdvection.h"

double Equations::Scalar1D::LinearAdvection::sine(const double& x,
		const double& dx, const Vector9d& weights, const Vector9d& absis) {

	double uAverage = 0.0;

	for (int i = 0; i < 9; i++) {
		uAverage += weights(i) *
				this->fSine(x + absis(i) * dx / 2.0) / 2.0;
	}

	return uAverage;
}
double Equations::Scalar1D::LinearAdvection::sinedx(const double& x,
		const double& dx, const Vector9d& weights, const Vector9d& absis) {

	double vAverage = 0.0;
	for (int i = 0; i < 9; i++) {
		vAverage += weights(i) *
				this->dfSine(x + absis(i) * dx / 2.0) / 2.0;
	}

	return vAverage;
}
