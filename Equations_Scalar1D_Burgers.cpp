#include "Equations_Scalar1D_Burgers.h"

double Equations::Scalar1D::Burgers::f(const double& u) {
	return 0.5 * u * u;
}

double Equations::Scalar1D::Burgers::df(const double& u) {
	return u;
}

double Equations::Scalar1D::Burgers::Rusanov(const double& ul, const double& ur) {
	return 0.5 * (f(ur) + f(ul) - fmax(fabs(df(ul)), fabs(df(ur))) * (ur - ul));
}

double Equations::Scalar1D::Burgers::dRusanov(const double& ul, const double& ur,
		const double& dul, const double& dur) {

	return 0.5 * (df(ur) * dur + df(ul) * dul - fmax(fabs(df(ul)), fabs(df(ur))) * (dur - dul));
}

double Equations::Scalar1D::Burgers::LxF(const double& ul, const double& ur, const double& cmax) {
    return 0.5 * (f(ur) + f(ul) - cmax * (ur - ul));
}

double Equations::Scalar1D::Burgers::dLxF(const double& ul, const double& ur,
                                              const double& dul, const double& dur, const double& cmax) {

    return 0.5 * (df(ur) * dur + df(ul) * dul - cmax * (dur - dul));
}