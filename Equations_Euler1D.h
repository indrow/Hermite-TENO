#pragma once
#include <cmath>
#include <iostream>
#include <Eigen/Dense>

using std::cout;
using std::endl;
using std::fmax;
using Eigen::Matrix;
using Eigen::Vector3d;
using Eigen::Matrix3d;
typedef Matrix<double, 9, 1> Vector9d;

namespace Equations {
class Euler1D
{
public:
	Euler1D(const double& specificHeatRatio);
	Vector3d primitive(const Vector3d& q);
	Vector3d conservative(const Vector3d& w);
	Vector3d flux(const Vector3d& q, const Vector3d& w);
	Matrix3d dflux(const Vector3d& q, const Vector3d& w);
	Matrix3d rightEigenVector(const double& vxm, const double& hm,
			const double& qm, const double& cm, const double& t0);
	Matrix3d leftEigenVector(const double& rcm, const double& b1,
			const double& b2, const double& t0, const double& t1, const double& t2);
	double maxSpeed(const Vector3d& w);

	Vector3d maxLocalWaveSpeed(const Vector3d& wl, const Vector3d& wr);
	Vector3d localWaveSpeed(const Vector3d& w);
	Vector3d Rusanov(const Vector3d& ql, const Vector3d& qr, const Vector3d& wl, const Vector3d& wr, const Matrix3d& leftEigenVec,
			const Matrix3d& rightEigenVec, const Vector3d& maxSpeed);
	Vector3d CUFlux(const Vector3d& ql, const Vector3d& qr, const Matrix3d& leftEigenVec, const Matrix3d& rightEigenVec);
	Vector3d HLLC(const Vector3d& ql, const Vector3d& qr);
	Vector3d dRusanov(const Vector3d& dql, const Vector3d& dqr, const Vector3d& wl,
			const Vector3d& wr, const Matrix3d& leftEigenVec, const Matrix3d& rightEigenVec,
			const Matrix3d& dfluxl, const Matrix3d& dfluxr, const Vector3d& maxSpeed);
	Vector3d dCUFlux(const Vector3d& ql, const Vector3d& qr, const Matrix3d& leftEigenVec, const Matrix3d& rightEigenVec,
			const Vector3d& dql, const Vector3d& dqr, const Matrix3d& dfluxl, const Matrix3d& dfluxr);

	Vector3d entropyWaves(const double& rhoInf, const double& uInf,
			const double& pInf, const double& amplitude, const double& x,
			const double& dx, const Matrix<double, 9, 1>& weights, const Matrix<double, 9, 1>& absis);

	Vector3d entropyWavesdx(const double& rhoInf, const double& uInf,
			const double& pInf, const double& amplitude, const double& x,
			const double& dx, const Matrix<double, 9, 1>& weights, const Matrix<double, 9, 1>& absis);

	Vector3d interactingBlast(const double& x);
	Vector3d sodShock(const double& x);
	Vector3d shuOsher(const double& x);

private:
	Vector3d conservativeOfEntropyWave(const double& rhoInf, const double& uInf,
			const double& pInf, const double& amplitude, const double& x);
	Vector3d dconservativeOfEntropyWave(const double& rhoInf, const double& uInf,
			const double& pInf, const double& amplitude, const double& x);

	const double gamma, gamma_1, gamma_3, pi;
	double squareVelocityx, dynVelocityx, enthalpyx;

	double soundSpeedLx, soundSpeedRx, velocityLx, velocityRx;

	double soundSpeed;

	double soundSpeedLloc, soundSpeedRloc, velocityLloc, velocityRloc;
	double soundSpeedLocal, velocity;

	Vector3d initEntWave;
	Vector3d outhx;
	Vector3d wLoc, qLoc, fxLoc, maxWave, rusanovFlux, rusanovFluxH, waveSpeed;
	Matrix3d dfxLoc, eigRx, eigLx;
};
}
