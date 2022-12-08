#pragma once
#include "Equations_Euler2D.h"
#include <fstream>
#include <iomanip>
#include "Schemes_Lagrange.h"
#include "Schemes_SSPRK.h"

namespace FV {
	class System2D
	{
	public:
		System2D(const double& xl, const double& xr, const int& nx, \
			const double& yb, const double& yt, const int& ny, \
			const double& cflnumber, const double& tmax, const double& tInit, \
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

			double sqrtRho{ 0 }, enthalpy{ 0 };

			Matrix4d leftEigenVectorX, rightEigenVectorX,
				leftEigenVectorY, rightEigenVectorY;

			Vector4d exact;

			void transformToConservative() {
				uxl = rightEigenVectorX * uxl;
				uxr = rightEigenVectorX * uxr;

				uyl = rightEigenVectorY * uyl;
				uyr = rightEigenVectorY * uyr;
			}

			void transformToConservativeAtGQ() {
				for (int i = 0; i < 3; i++) {
					uxgl.col(i) = rightEigenVectorX * uxgl.col(i);
					uxgr.col(i) = rightEigenVectorX * uxgr.col(i);

					uygl.col(i) = rightEigenVectorY * uygl.col(i);
					uygr.col(i) = rightEigenVectorY * uygr.col(i);
				}
			}

			void setZeroFlux() {
				fx.setZero(); fy.setZero();
			}

			void save() { u0 = u; }
			void update1() { u = u1; }
			void update2() { u = u2; }

		};

	private:
		const double gamma, gamma_1;

		void calcPrimitiveVariable();
		void periodicBc();
		void zeroGradBc();
		void TaylorInstBc();
		void doubleMachBc();
		void eigen2DEuler();
		void rhs();
		void RungeKutta3();
		void write(const int& filenum);

		int nx, ny, k, s, sh, nxx, nyy, n, nt;
		double dx, dy, area;
		double cfl, dt, t, tFinal;

		Matrix<FV::System2D::SystemData, Dynamic, Dynamic> cell;

		Equations::Euler2D eulObj;
		Schemes::Lagrange lagrangeObj;
		Schemes::SSPRK ssprkObj;

		Vector2d twoPointsGQWeights;
		Vector3d threePointsGQWeights;
		Vector4d fourPointsGQWeights;

		double t0, t1, vxm, vym, hm, qm, cm;
		double rcm, b1, b2, t2;

		Matrix<double, 4, 6> interpolateC;
		Vector2d valueB;

		Matrix<double, 4, 5> ee, ff;
		Matrix<double, 4, 1> valueC;
        Matrix<double, 4, 3> temporaryFluxUx, temporaryFluxUy;

		Vector4d source, pvarL, pvarR, aM, eigenValueX, eigenValueY;
	};
}