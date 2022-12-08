#include "Schemes_SSPRK.h"

template <typename T> T Schemes::SSPRK::RungeKuttaIII(const T& q0, const T& qn, const T& rhs, \
	const double& dt, const int& stage) {

	if (stage == 1) {
		return q0 + dt * rhs;
	}
	else if (stage == 2) {
		return 0.75 * q0 + 0.25 * (qn + rhs * dt);
	}
	else {
		return (q0 + 2.0 * (qn + dt * rhs)) / 3.0;
	}
}

template double Schemes::SSPRK::RungeKuttaIII<double>(const double&, const double&, const double&, \
	const double&, const int&);
template Vector4d Schemes::SSPRK::RungeKuttaIII<Vector4d>(const Vector4d&, const Vector4d&, const Vector4d&, \
	const double&, const int&);
template Vector3d Schemes::SSPRK::RungeKuttaIII<Vector3d>(const Vector3d&,
		const Vector3d&, const Vector3d&, const double&, const int&);
