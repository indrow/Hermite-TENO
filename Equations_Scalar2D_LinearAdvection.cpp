#include "Equations_Scalar2D_LinearAdvection.h"

double Equations::Scalar2D::LinearAdvection::sine(const double& x, const double& y,
		const double& dx, const double& dy,
		const Vector9d& weights, const Vector9d& absis) {

	int i, j;
	double uAverage = 0.0;
	for (i = 0; i < 9; i++) {
		for (j = 0; j < 9; j++) {
			uAverage += weights(i) * weights(j) *
					this->fSine(x + absis(i) * dx / 2.0, y + absis(j) * dy / 2.0) / 4.0;
		}
	}
	return uAverage;
}
double Equations::Scalar2D::LinearAdvection::sinedx(const double& x, const double& y,
		const double& dx, const double& dy,
		const Vector9d& weights, const Vector9d& absis) {

	int i, j;
	double vAverage = 0.0;
	for (i = 0; i < 9; i++) {
		for (j = 0; j < 9; j++) {
			vAverage += weights(i) * weights(j) *
					this->dfSine(x + absis(i) * dx / 2.0, y + absis(j) * dy / 2.0) / 4.0;
		}
	}

	return vAverage;
}
double Equations::Scalar2D::LinearAdvection::sinedy(const double& x, const double& y,
		const double& dx, const double& dy,
		const Vector9d& weights, const Vector9d& absis) {

	int i, j;
	double wAverage = 0.0;
	for (i = 0; i < 9; i++) {
		for (j = 0; j < 9; j++) {
			wAverage += weights(i) * weights(j) *
					this->dfSine(x + absis(i) * dx / 2.0, y + absis(j) * dy / 2.0) / 4.0;
		}
	}
	return wAverage;
}

double Equations::Scalar2D::LinearAdvection::f(const double& u) {
	return u;
}

double Equations::Scalar2D::LinearAdvection::df(const double& u) {
	return 1.0;
}

double Equations::Scalar2D::LinearAdvection::Rusanov(const double& ul, const double& ur) {
	return 0.5 * (f(ur) + f(ul) - fmax(fabs(df(ul)), fabs(df(ur))) * (ur - ul));
}

double Equations::Scalar2D::LinearAdvection::dRusanov(const double& ul, const double& ur,
		const double& dul, const double& dur) {

	return 0.5 * (df(ur) * dur + df(ul) * dul - fmax(fabs(df(ul)), fabs(df(ur))) * (dur - dul));
}

double Equations::Scalar2D::LinearAdvection::LxF(const double& ul, const double& ur, const double& cmax) {
    return 0.5 * (f(ur) + f(ul) - cmax * (ur - ul));
}

double Equations::Scalar2D::LinearAdvection::dLxF(const double& ul, const double& ur,
                                                      const double& dul, const double& dur, const double& cmax) {

    return 0.5 * (df(ur) * dur + df(ul) * dul - cmax * (dur - dul));
}