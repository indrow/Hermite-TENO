#pragma once
#include "Equations_Scalar2D_LinearAdvection.h"
#include "Schemes_Hermite.h"
#include "Schemes_Lagrange.h"
#include "Schemes_SSPRK.h"

namespace FV {
class Scalar2D
{
public:
	Scalar2D(const double& xl, const double& xr, const int& nx,
			const double& yb, const double& yt, const int& ny,
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
		double fx{ 0 }, fy{ 0 };
		double rhs{ 0 };
		double exact{ 0 };

		Vector2d ux, uy;
		Vector2d uxg2l, uxg2r, uyg2l, uyg2r;
		Vector4d uxg4l, uxg4r, uyg4l, uyg4r;

		void save() { u0 = u; }
		void update1() { u = u1; }
		void update2() { u = u2; }
	}; // End of struct of ScalarData

private:
	Matrix <FV::Scalar2D::ScalarData, Dynamic, Dynamic> cell;
	int nx, ny, k, s, nxx, nyy, n, nt;
	double dx, dy, area;
	double cfl, dt, t, tFinal;
	Equations::Scalar2D::LinearAdvection linAdvObj;
	Schemes::Hermite hermiteObj;
	Schemes::Lagrange lagrangeObj;
	Schemes::SSPRK ssprkObj;

	Vector2d twoPointsGQWeights;
	Vector4d fourPointsGLegendreQWeights;

	Vector6d interpolateC;
	Vector2d valueB;
	Vector5d ee, ff;

	void GaussLegendre2Points();
	void GaussLobatto4Points();
};
}

