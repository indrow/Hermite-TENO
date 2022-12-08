#pragma once
#include "Equations_Scalar1D_Burgers.h"
#include "Equations_Scalar1D_LinearAdvection.h"
#include "Equations_Scalar2D_LinearAdvection.h"
#include <fstream>
#include <iomanip>
#include "Schemes_Hermite.h"
#include "Schemes_Lagrange.h"
#include "Schemes_SSPRK.h"

namespace FV {
class Scalar1D
{
public:
	Scalar1D(const double& xl, const double& xr, const int& nx,
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
		double x{ 0 };

		double u{ 0 }, u0{ 0 }, u1{ 0 }, u2{ 0 };
		double uxl{ 0 }, uxr{ 0 };
		double fx{ 0 };
		double rhsu{ 0 };

		double exact{ 0 };
		void save() { u0 = u; }
		void update1() { u = u1; }
		void update2() { u = u2; }
		void setZeroFlux() { fx = 0; }
	};

private:
	int nx, k, s, nxx, n;
	double dx, cmax, area;
	double cfl, dt, t, tFinal;

	Matrix<FV::Scalar1D::ScalarData, Dynamic, 1> cell;
	Equations::Scalar1D::LinearAdvection linAdv1DObj;
	Equations::Scalar2D::LinearAdvection linAdvObj;
	Equations::Scalar1D::Burgers burgersObj;
	Schemes::Lagrange lagrangeObj;
	Schemes::SSPRK ssprkObj;

	Vector6d interpolateA;
	Vector2d valueA;
};
}

