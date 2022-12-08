#pragma once
#include <cmath>
#include <iostream>
#include <Eigen/Dense>

using std::cout;
using std::endl;
using std::fmax;
using std::fmin;
using std::isnan;
using Eigen::Matrix;
using Eigen::Vector2d;
using Eigen::Vector4d;
using Eigen::Matrix4d;
typedef Matrix<double, 9, 1> Vector9d;

namespace Equations {
class Euler2D
{
public:
	Euler2D(const double& specificHeatRatio);
	Vector4d primitive(const Vector4d& q);
	Vector4d conservative(const Vector4d& w);
	Vector4d fluxx(const Vector4d& q, const Vector4d& w);
	Vector4d fluxy(const Vector4d& q, const Vector4d& w);
	Matrix4d dfluxx(const Vector4d& q, const Vector4d& w);
	Matrix4d dfluxy(const Vector4d& q, const Vector4d& w);
	Matrix4d rightEigenVectorx(const double& vxm, const double& vym, const double& hm,
			const double& qm, const double& cm, const double& t0);
	Matrix4d leftEigenVectorx(const double& rcm, const double& vym, const double& b1,
			const double& b2, const double& t0, const double& t1, const double& t2);
	Matrix4d rightEigenVectory(const double& vxm, const double& vym, const double& hm,
			const double& qm, const double& cm, const double& t0);
	Matrix4d leftEigenVectory(const double& rcm, const double& vxm, const double& b1,
			const double& b2, const double& t0, const double& t1, const double& t2);
	double maxSpeed(const Vector4d& w, int dir);

	Matrix4d maxLocalWaveSpeed(const Vector4d& wl, const Vector4d& wr, const Matrix4d& leftEigenVec, const Matrix4d& rightEigenVec, int dir);
    Vector4d maxLocalWaveSpeed(const Vector4d& wl, const Vector4d& wr, int dir);
	Vector4d localWaveSpeed(const Vector4d& w, int dir);

	Vector4d Rusanov(const Vector4d& ql, const Vector4d& qr, const Vector4d& wl, const Vector4d& wr, const Matrix4d& leftEigenVec,
			const Matrix4d& rightEigenVec, const Vector4d& maxSpeed, int dir);
    Vector4d Rusanov(const Vector4d& ql, const Vector4d& qr, const Vector4d& wl, const Vector4d& wr,
                     const Matrix4d& maxSpeed, int dir);

	Vector4d dRusanov(const Vector4d& dql, const Vector4d& dqr, const Vector4d& wl,
			const Vector4d& wr, const Matrix4d& leftEigenVec, const Matrix4d& rightEigenVec,
			const Matrix4d& dfluxl, const Matrix4d& dfluxr, const Vector4d& maxSpeed);
    Vector4d dRusanov(const Vector4d& dql, const Vector4d& dqr, const Vector4d& wl,
                      const Vector4d& wr, const Matrix4d& dfluxl, const Matrix4d& dfluxr,
                      const Matrix4d& maxSpeed);

	Vector4d HLLC(const Vector4d& ql, const Vector4d& qr, int dir);

	Vector4d entropyWaves(const double& rhoInf, const double& uInf, const double& vInf,
			const double& pInf, const double& amplitude, const double& x, const double& y,
			const double& dx, const double& dy, const Matrix<double, 9, 1>& weights,
			const Matrix<double, 9, 1>& absis);
	Vector4d entropyWavesdx(const double& rhoInf, const double& uInf, const double& vInf,
			const double& pInf, const double& amplitude, const double& x, const double& y,
			const double& dx, const double& dy, const Matrix<double, 9, 1>& weights,
			const Matrix<double, 9, 1>& absis);
	Vector4d entropyWavesdy(const double& rhoInf, const double& uInf, const double& vInf,
			const double& pInf, const double& amplitude, const double& x, const double& y,
			const double& dx, const double& dy, const Matrix<double, 9, 1>& weights,
			const Matrix<double, 9, 1>& absis);

	Vector4d isentropicVortex(const double& x, const double& y, const double& dx, const double& dy,
                           const Matrix<double, 9, 1>& weights, const Matrix<double, 9, 1>& absis);
    Vector4d dxIsentropicVortex(const double& x, const double& y, const double& dx, const double& dy,
                                const Matrix<double, 9, 1>& weights, const Matrix<double, 9, 1>& absis);
    Vector4d dyIsentropicVortex(const double& x, const double& y, const double& dx, const double& dy,
                                const Matrix<double, 9, 1>& weights, const Matrix<double, 9, 1>& absis);


	Vector4d conservativeOfRiemann2D(const double& x, const double& y);
	Vector4d conservativeOfTaylorIns2D(const double& x, const double& y);
	Vector4d conservativeOfDoubleMach2D(const double& x, const double& y);

private:
	const double gamma, gamma_1, gamma_3, pi;

	Vector4d conservativeOfEntropyWave(const double& rhoInf, const double& uInf, const double& vInf,
			const double& pInf, const double& amplitude, const double& x, const double& y);
	Vector4d dconservativeOfEntropyWave(const double& rhoInf, const double& uInf, const double& vInf,
			const double& pInf, const double& amplitude, const double& x, const double& y);

    Vector4d conservativeOfIsentropicVortex( const double& x, const double& y);
    Vector4d dxconservativeOfIsentropicVortex(const double& x, const double& y);
    Vector4d dyconservativeOfIsentropicVortex(const double& x, const double& y);

	double squareVelocityx, dynVelocityx, enthalpyx;
	double squareVelocityy, dynVelocityy, enthalpyy;

	double soundSpeedLx, soundSpeedRx, velocityLx, velocityRx;
	double soundSpeedLy, soundSpeedRy, velocityLy, velocityRy;

	double soundSpeed;

	double soundSpeedLloc, soundSpeedRloc, velocityLloc, velocityRloc;
	double soundSpeedLocal, velocity;

	Vector4d initEntWave;
	Vector4d outhx, outhy;
	Vector4d wLoc, qLoc, fxLoc, fyLoc, maxWave, rusanovFlux, rusanovFluxH, waveSpeed;
	Matrix4d dfxLoc, dfyLoc, eigRx, eigLx, eigRy, eigLy, dfdu;
};
}

