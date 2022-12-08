#pragma once
#include "Equations_Euler2D.h"
#include <fstream>
#include <iomanip>
#include "Schemes_Hermite.h"
#include "Schemes_Lagrange.h"
#include "Schemes_SSPRK.h"
#include <omp.h>

using std::isnan;

namespace FVHermite {
class System2D
{
public:
	System2D(const double& xl, const double& xr, const int& nx,
			const double& yb, const double& yt, const int& ny,
			const double& cflnumber, const double& tmax, const double& tInit,
			const double& specificHeatRatio);
	void entropyWave();
    void isentropicVortex();
	void Riemann2D();
	void timeStepping();
	void printResult();

protected:
	struct SystemData
	{
		double x{ 0 }, y{ 0 };

		Vector4d pvar;

		Vector4d u, u0, u1, u2;
		Vector4d uxl, uxr, uyl, uyr;
		Vector4d fx, fy;
		Vector4d rhsu;
		Matrix<double, 4, 3> uxgl, uxgr, uygl, uygr;

		Vector4d v, v0, v1, v2;
		Vector4d vxl, vxr, vyl, vyr;
		Vector4d gx, gy;
		Vector4d rhsv;
		Matrix<double, 4, 3> vxgl, vxgr, vygl, vygr;

		Vector4d w, w0, w1, w2;
		Vector4d wxl, wxr, wyl, wyr;
		Vector4d hx, hy;
		Vector4d rhsw;
		Matrix<double, 4, 3> wxgl, wxgr, wygl, wygr;

		double sqrtRho{ 0 }, enthalpy{ 0 };

		Matrix4d leftEigenVectorX, rightEigenVectorX,
		leftEigenVectorY, rightEigenVectorY;

		Vector4d exact;

		void transformToConservative() {
			uxl = rightEigenVectorX * uxl;
			vxl = rightEigenVectorX * vxl;
			wxl = rightEigenVectorX * wxl;

			uxr = rightEigenVectorX * uxr;
			vxr = rightEigenVectorX * vxr;
			wxr = rightEigenVectorX * wxr;

			uyl = rightEigenVectorY * uyl;
			vyl = rightEigenVectorY * vyl;
			wyl = rightEigenVectorY * wyl;

			uyr = rightEigenVectorY * uyr;
			vyr = rightEigenVectorY * vyr;
			wyr = rightEigenVectorY * wyr;
		}

		void transformToConservativeAtGQ() {
			for (int i = 0; i < 3; i++) {
				uxgl.col(i) = rightEigenVectorX * uxgl.col(i);
				vxgl.col(i) = rightEigenVectorX * vxgl.col(i);
				wxgl.col(i) = rightEigenVectorX * wxgl.col(i);

				uxgr.col(i) = rightEigenVectorX * uxgr.col(i);
				vxgr.col(i) = rightEigenVectorX * vxgr.col(i);
				wxgr.col(i) = rightEigenVectorX * wxgr.col(i);

				uygl.col(i) = rightEigenVectorY * uygl.col(i);
				vygl.col(i) = rightEigenVectorY * vygl.col(i);
				wygl.col(i) = rightEigenVectorY * wygl.col(i);

				uygr.col(i) = rightEigenVectorY * uygr.col(i);
				vygr.col(i) = rightEigenVectorY * vygr.col(i);
				wygr.col(i) = rightEigenVectorY * wygr.col(i);
			}
		}

		void setZeroFlux() {
			fx.setZero(); fy.setZero(); gx.setZero();
			gy.setZero(); hx.setZero(); hy.setZero();
		}

		void save() { u0 = u; v0 = v; w0 = w; }
		void update1() { u = u1; v = v1; w = w1; }
		void update2() { u = u2; v = v2; w = w2; }
	}; // end of structure

private:
	const double gamma, gamma_1;

	void calcPrimitiveVariable();
	void periodicBc();
	void zeroGradBc();
	void TaylorInstBc();
	void doubleMachBc();
	void eigen2DEuler();
	void limitDerivatives();
	void rhs();
	void RungeKutta3();
	void write(const int& filenum);

	int nx, ny, k, s, sh, nxx, nyy, n, nt;
	double dx, dy, area;
	double cfl, dt, t, tFinal;

	Matrix<FVHermite::System2D::SystemData, Dynamic, Dynamic> cell;

	Equations::Euler2D eulObj;
	Schemes::Hermite hermiteObj;
	Schemes::Lagrange lagrangeObj;
	Schemes::SSPRK ssprkObj;

	Vector2d twoPointsGQWeights;
	Vector3d threePointsGQWeights;

	double t0, t1, vxm, vym, hm, qm, cm;
	double rcm, b1, b2, t2;

	Matrix<double, 4, 4> interpolateA, interpolateB;
	Vector4d valueA, fourPointsGQWeights;

	Matrix<double, 4, 6> interpolateC;
	Vector2d valueB;

	Matrix<double, 4, 5> aa, bb, cc, dd;
	Matrix<double, 4, 5> ee, ff;
	Matrix<double, 8, 1> valueC;

	Vector4d source, pvarL, pvarR, aM;
	Matrix4d dFluxL, dFluxR, dfdu;
    Matrix<double, 4, 3> temporaryFluxUx, temporaryFluxUy;
    Matrix<double, 4, 3> temporaryFluxVx, temporaryFluxVy;
    Matrix<double, 4, 3> temporaryFluxWx, temporaryFluxWy;
};
}
