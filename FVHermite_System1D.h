#pragma once
#include "Equations_Euler1D.h"
#include <fstream>
#include <iomanip>
#include "Schemes_Hermite.h"
#include "Schemes_Lagrange.h"
#include "Schemes_SSPRK.h"

namespace FVHermite {

class System1D
{
public:
	System1D(const double& xl, const double& xr, const int& nx,
			const double& cflnumber, const double& tmax, const double& tInit,
			const double& specificHeatRatio);
	void entropyWave();
	void blastWave();
	void timeStepping();
	void printResult();

	void calcPrimitiveVariable();
	void periodicBc();
	void zeroGradBc();
	void reflectiveBc();
	void eigen1DEuler();
	void rhs();
	void RungeKutta3();

protected:
	struct SystemData
	{
		double x{ 0 };

		Vector3d pvar;

		Eigen::Matrix<bool, 3, 1> discontinue;

		Vector3d u, u0, u1, u2;
		Vector3d uxl, uxr;
		Vector3d fx;
		Vector3d rhsu;

		Vector3d v, v0, v1, v2;
		Vector3d vxl, vxr;
		Vector3d gx;
		Vector3d rhsv;

		bool c;

		double sqrtRho{ 0 }, enthalpy{ 0 };

		Matrix3d leftEigenVectorX, rightEigenVectorX;

		Vector3d exact;

		void transformToConservative() {
			uxl = rightEigenVectorX * uxl;
			vxl = rightEigenVectorX * vxl;

			uxr = rightEigenVectorX * uxr;
			vxr = rightEigenVectorX * vxr;
		}


		void setZeroFlux() {
			fx.setZero(); gx.setZero();
		}

		void save() { u0 = u; v0 = v; }
		void update1() { u = u1; v = v1; }
		void update2() { u = u2; v = v2; }
	}; // end of structure

private:
	const double gamma, gamma_1;

	void write(const int& filenum);

	int nx, k, s, nxx, n;
	float comp_time;
	double h1, h2, Ii, circleRadMod;
	double dx, area;
	double cfl, dt, t, tFinal;

	Matrix<FVHermite::System1D::SystemData, Dynamic, 1> cell;

	Equations::Euler1D eulObj;
	Schemes::Hermite hermiteObj;
	Schemes::Lagrange lagrangeObj;
	Schemes::SSPRK ssprkObj;

	double t0, t1, vxm, hm, qm, cm;
	double rcm, b1, b2, t2;

	Matrix<double, 3, 4> interpolateA, interpolateB;
	Vector4d valueA;

	Vector3d pvarL, pvarR, aM;
	Matrix3d dFluxL, dFluxR;

	void limitDerivatives();
};
}
