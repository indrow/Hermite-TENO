#pragma once
#include "Equations_Scalar2D_LinearAdvection.h"
#include "Schemes_Hermite.h"
#include "Schemes_Lagrange.h"
#include "Schemes_SSPRK.h"

using Eigen::Vector2d;
using Eigen::Dynamic;

namespace FVHermite {
	class Scalar2D {
	public:
		Scalar2D(const double& xl, const double& xr, const int& nx, \
			const double& yb, const double& yt, const int& ny, \
			const double& cflnumber, const double& tmax, const double& tInit);
		void initialize();
		void periodicBc();
		void rhs();
		void RungeKutta3();
		void timeStepping();
		void Exact();
		void printResult();

	protected:
		struct ScalarData
		{
			double x{ 0 }, y{ 0 };

			double u{ 0 }, u0{ 0 }, u1{ 0 }, u2{ 0 };
			double uxl{ 0 }, uxr{ 0 }, uyl{ 0 }, uyr{ 0 };
			double fx{ 0 }, fy{ 0 };
			double rhsu{ 0 };
			Vector2d uxgl, uxgr, uygl, uygr;

			double v{ 0 }, v0{ 0 }, v1{ 0 }, v2{ 0 };
			double vxl{ 0 }, vxr{ 0 }, vyl{ 0 }, vyr{ 0 };
			double gx{ 0 }, gy{ 0 };
			double rhsv{ 0 };
			Vector2d vxgl, vxgr, vygl, vygr;

			double w{ 0 }, w0{ 0 }, w1{ 0 }, w2{ 0 };
			double wxl{ 0 }, wxr{ 0 }, wyl{ 0 }, wyr{ 0 };
			double hx{ 0 }, hy{ 0 };
			double rhsw{ 0 };
			Vector2d wxgl, wxgr, wygl, wygr;

			double exact{ 0 };
			void save() { u0 = u; v0 = v; w0 = w; }
			void update1() { u = u1; v = v1; w = w1; }
			void update2() { u = u2; v = v2; w = w2; }
			void setZeroFlux() { fx = 0; fy = 0; gx = 0; gy = 0; hx = 0; hy = 0; }
		};

		private:
			int nx, ny, k, s, sh, nxx, nyy, n, nt;
			double dx, dy, area;
			double cfl, dt, t, tFinal;
			Matrix<FVHermite::Scalar2D::ScalarData, Dynamic, Dynamic> cell;
			Equations::Scalar2D::LinearAdvection linAdvObj;
			Schemes::Hermite hermiteObj;
			Schemes::Lagrange lagrangeObj;
			Schemes::SSPRK ssprkObj;

			Vector2d twoPointsGQWeights;

			Vector4d interpolateA, interpolateB;
			Vector4d valueA;
			Vector6d interpolateC;
			Vector2d valueB;

			Vector3d aa, bb, cc, dd;
			Vector5d ee, ff;
		};
}