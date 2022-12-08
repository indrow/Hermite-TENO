#include "FVHermite_Scalar2D.h"

FVHermite::Scalar2D::Scalar2D(const double& xl, const double& xr, const int& nx, \
	const double& yb, const double& yt, const int& ny, \
	const double& cflnumber, const double& tmax, const double& tInit) {

	int i, j;

	this->cfl = cflnumber;this->tFinal = tmax; this->t = tInit; dt = 0.0;
	this->nx = nx; this->ny = ny;
	dx = (xr - xl) / nx;
	dy = (yt - yb) / ny;
	area = dx * dy;
	k = 5; s = (k - 1) / 2; sh = 1;
	k += 3;
	nxx = k + k + nx;
	nyy = k + k + ny;
	n = nx * ny;
	nt = nxx * nyy;
	hermiteObj.setGridSize(dx);
	twoPointsGQWeights = lagrangeObj.getQuadrature2Weights() / 2.0;

	cout << "Creating domain... ";
	cell.resize(nxx, nyy);
	for (j = 0; j < nyy; j++) {
		for (i = 0; i < nxx; i++) {
			cell(i, j).x = (double(i) - double(k) + 0.5) * dx + xl;
			cell(i, j).y = (double(j) - double(k) + 0.5) * dy + yb;
		}
	}

	cout << "done.\n";
	cout << "Number of cells: " << n << "\n";
}

void FVHermite::Scalar2D::initialize() {
	int i, j;
	cout << "Initializing variable...";
	for (i = k; i < nxx - k; i++) {
		for (j = k; j < nyy - k; j++) {
			cell(i, j).u = linAdvObj.sine(cell(i, j).x, cell(i, j).y, dx, dy, \
				lagrangeObj.getQuadrature9Weights(), lagrangeObj.getQuadrature9Points());
			cell(i, j).v = linAdvObj.sinedx(cell(i, j).x, cell(i, j).y, dx, dy, \
				lagrangeObj.getQuadrature9Weights(), lagrangeObj.getQuadrature9Points());
			cell(i, j).w = linAdvObj.sinedy(cell(i, j).x, cell(i, j).y, dx, dy, \
				lagrangeObj.getQuadrature9Weights(), lagrangeObj.getQuadrature9Points());
		}
	}
	cout << "done.\n";
}

void FVHermite::Scalar2D::periodicBc() {
	int i, j, m;

	for (j = k; j < ny + k; j++) {
		for (m = 0; m < k; m++) {
			cell(k - 1 - m, j).u = cell(nx + k - 1 - m, j).u;
			cell(nxx - 1 - m, j).u = cell(k + k - 1 - m, j).u;

			cell(k - 1 - m, j).v = cell(nx + k - 1 - m, j).v;
			cell(nxx - 1 - m, j).v = cell(k + k - 1 - m, j).v;

			cell(k - 1 - m, j).w = cell(nx + k - 1 - m, j).w;
			cell(nxx - 1 - m, j).w = cell(k + k - 1 - m, j).w;
		}
	}

	for (i = 0; i < nxx; i++) {
		for (m = 0; m < k; m++) {
			cell(i, k - 1 - m).u = cell(i, ny + k - 1 - m).u;
			cell(i, nyy - 1 - m).u = cell(i, k + k - 1 - m).u;

			cell(i, k - 1 - m).v = cell(i, ny + k - 1 - m).v;
			cell(i, nyy - 1 - m).v = cell(i, k + k - 1 - m).v;

			cell(i, k - 1 - m).w = cell(i, ny + k - 1 - m).w;
			cell(i, nyy - 1 - m).w = cell(i, k + k - 1 - m).w;
		}
	}
}

void FVHermite::Scalar2D::rhs() {
	int i, j, m;
	periodicBc();

	for (i = s; i < nxx - s - 1; ++i) {
		for (j = 0; j < nyy; ++j) {
			for (m = 0; m < 4; ++m) {
				interpolateA(m) = cell(i - 1 + m, j).u;
				interpolateB(m) = cell(i - 1 + m, j).v;
			}
			valueA = hermiteObj.hweno(interpolateA, interpolateB);
			cell(i, j).uxl = valueA(0); cell(i, j).uxr = valueA(1);
			cell(i, j).vxl = valueA(2); cell(i, j).vxr = valueA(3);
			for (m = 0; m < 6; ++m) {
				interpolateC(m) = cell(i - 2 + m, j).w;
			}
			valueB = lagrangeObj.weno5JS(interpolateC);
			cell(i, j).wxl = valueB(0); cell(i, j).wxr = valueB(1);
		}
	}

	for (j = s; j < nyy - s - 1; ++j) {
		for (i = 0; i < nxx; ++i) {
			for (m = 0; m < 4; ++m) {
				interpolateA(m) = cell(i, j - 1 + m).u;
				interpolateB(m) = cell(i, j - 1 + m).w;
			}
			valueA = hermiteObj.hweno(interpolateA, interpolateB);
			cell(i, j).uyl = valueA(0); cell(i, j).uyr = valueA(1);
			cell(i, j).wyl = valueA(2); cell(i, j).wyr = valueA(3);

			for (m = 0; m < 6; ++m) {
				interpolateC(m) = cell(i, j - 2 + m).v;
			}
			valueB = lagrangeObj.weno5JS(interpolateC);
			cell(i, j).vyl = valueB(0); cell(i, j).vyr = valueB(1);
		}
	}

	for (i = s + 2; i < nxx - s - 3; ++i) {
		for (j = s + 2; j < nyy - s - 3; ++j) {
			for (m = 0; m < 3; ++m) {
				aa(m) = cell(i, j - 1 + m).uxl;
				bb(m) = cell(i, j - 1 + m).vxl;
				cc(m) = cell(i, j - 1 + m).uxr;
				dd(m) = cell(i, j - 1 + m).vxr;
			}
			valueA = hermiteObj.hweno2PointsGLQ(aa, bb);
			cell(i, j).uxgl = valueA.segment(0, 2);
			cell(i, j).vxgl = valueA.segment(2, 2);
			valueA = hermiteObj.hweno2PointsGLQ(cc, dd);
			cell(i, j).uxgr = valueA.segment(0, 2);
			cell(i, j).vxgr = valueA.segment(2, 2);
			for (m = 0; m < 5; ++m) {
				ee(m) = cell(i, j - 2 + m).wxl;
				ff(m) = cell(i, j - 2 + m).wxr;
			}
			cell(i, j).wxgl = lagrangeObj.weno5JS2PointsQuadrature(ee);
			cell(i, j).wxgr = lagrangeObj.weno5JS2PointsQuadrature(ff);

			for (m = 0; m < 3; ++m) {
				aa(m) = cell(i - 1 + m, j).uyl;
				bb(m) = cell(i - 1 + m, j).wyl;
				cc(m) = cell(i - 1 + m, j).uyr;
				dd(m) = cell(i - 1 + m, j).wyr;
			}
			valueA = hermiteObj.hweno2PointsGLQ(aa, bb);
			cell(i, j).uygl = valueA.segment(0, 2);
			cell(i, j).wygl = valueA.segment(2, 2);
			valueA = hermiteObj.hweno2PointsGLQ(cc, dd);
			cell(i, j).uygr = valueA.segment(0, 2);
			cell(i, j).wygr = valueA.segment(2, 2);
			for (m = 0; m < 5; ++m) {
				ee(m) = cell(i - 2 + m, j).vyl;
				ff(m) = cell(i - 2 + m, j).vyr;
			}
			cell(i, j).vygl = lagrangeObj.weno5JS2PointsQuadrature(ee);
			cell(i, j).vygr = lagrangeObj.weno5JS2PointsQuadrature(ff);
		}
	}


	for (i = s + 2; i < nxx - s - 3; ++i) {
		for (j = s + 2; j < nyy - s - 3; ++j) {
			cell(i, j).setZeroFlux();
			for (m = 0; m < 2; m++) {
				cell(i, j).fx += twoPointsGQWeights(m) * \
					linAdvObj.Rusanov(cell(i, j).uxgl(m), cell(i, j).uxgr(m));

				cell(i, j).fy += twoPointsGQWeights(m) * \
					linAdvObj.Rusanov(cell(i, j).uygl(m), cell(i, j).uygr(m));

				cell(i, j).gx += twoPointsGQWeights(m) * \
					linAdvObj.dRusanov(cell(i, j).uxgl(m), cell(i, j).uxgr(m), \
					cell(i, j).vxgl(m), cell(i, j).vxgr(m));

				cell(i, j).gy += twoPointsGQWeights(m) * \
					linAdvObj.dRusanov(cell(i, j).uygl(m), cell(i, j).uygr(m), \
					cell(i, j).vygl(m), cell(i, j).vygr(m));

				cell(i, j).hx += twoPointsGQWeights(m) * \
					linAdvObj.dRusanov(cell(i, j).uxgl(m), cell(i, j).uxgr(m), \
					cell(i, j).wxgl(m), cell(i, j).wxgr(m));

				cell(i, j).hy += twoPointsGQWeights(m) * \
					linAdvObj.dRusanov(cell(i, j).uygl(m), cell(i, j).uygr(m), \
					cell(i, j).wygl(m), cell(i, j).wygr(m));
			}
		}
	}

	for (i = k; i < nxx - k; i++) {
		for (j = k; j < nyy - k; j++) {
			cell(i, j).rhsu = ((cell(i - 1, j).fx - cell(i, j).fx) * dy + \
				(cell(i, j - 1).fy - cell(i, j).fy) * dx) / area;

			cell(i, j).rhsv = ((cell(i - 1, j).gx - cell(i, j).gx) * dy + \
				(cell(i, j - 1).gy - cell(i, j).gy) * dx) / area;

			cell(i, j).rhsw = ((cell(i - 1, j).hx - cell(i, j).hx) * dy + \
				(cell(i, j - 1).hy - cell(i, j).hy) * dx) / area;
		}
	}
}

void FVHermite::Scalar2D::RungeKutta3() {
	int i, j;
	int stage;

	stage = 1;
	rhs();
	for (i = k; i < nxx - k; i++) {
		for (j = k; j < nyy - k; j++) {
			cell(i, j).save();
			cell(i, j).u1 = ssprkObj.RungeKuttaIII<double>(cell(i, j).u0, \
				cell(i, j).u0, cell(i, j).rhsu, dt, stage);
			cell(i, j).v1 = ssprkObj.RungeKuttaIII<double>(cell(i, j).v0, \
				cell(i, j).v0, cell(i, j).rhsv, dt, stage);
			cell(i, j).w1 = ssprkObj.RungeKuttaIII<double>(cell(i, j).w0, \
				cell(i, j).w0, cell(i, j).rhsw, dt, stage);
			cell(i, j).update1();
		}
	}

	stage += 1;
	rhs();
	for (i = k; i < nxx - k; i++) {
		for (j = k; j < nyy - k; j++) {
			cell(i, j).u2 = ssprkObj.RungeKuttaIII<double>(cell(i, j).u0, \
				cell(i, j).u1, cell(i, j).rhsu, dt, stage);
			cell(i, j).v2 = ssprkObj.RungeKuttaIII<double>(cell(i, j).v0, \
				cell(i, j).v1, cell(i, j).rhsv, dt, stage);
			cell(i, j).w2 = ssprkObj.RungeKuttaIII<double>(cell(i, j).w0, \
				cell(i, j).w1, cell(i, j).rhsw, dt, stage);
			cell(i, j).update2();
		}
	}

	stage += 1;
	rhs();
	for (i = k; i < nxx - k; i++) {
		for (j = k; j < nyy - k; j++) {
			cell(i, j).u = ssprkObj.RungeKuttaIII<double>(cell(i, j).u0, \
				cell(i, j).u2, cell(i, j).rhsu, dt, stage);
			cell(i, j).v = ssprkObj.RungeKuttaIII<double>(cell(i, j).v0, \
				cell(i, j).v2, cell(i, j).rhsv, dt, stage);
			cell(i, j).w = ssprkObj.RungeKuttaIII<double>(cell(i, j).w0, \
				cell(i, j).w2, cell(i, j).rhsw, dt, stage);
		}
	}
}

void FVHermite::Scalar2D::timeStepping() {
	int i, j;
	double maxSpeedx, maxSpeedy;

	while (t < tFinal) {
		maxSpeedx = 0.0; maxSpeedy = 0.0;
		for (i = k; i < nxx - k; i++) {
			for (j = k; j < nyy - k; j++) {
				maxSpeedx = fmax(maxSpeedx, fabs(linAdvObj.df(cell(i, j).u)));
				maxSpeedy = fmax(maxSpeedy, fabs(linAdvObj.df(cell(i, j).u)));
			}
		}

		dt = cfl * pow(fmin(dx, dy), 5.0 / 3.0) / sqrt(maxSpeedx * maxSpeedx + maxSpeedy * maxSpeedy);
		dt = t + dt > tFinal ? tFinal - t : dt;

		RungeKutta3();
		t += dt;
		cout << "tnum: " << t << endl;
	}
}

void FVHermite::Scalar2D::Exact() {
	int i, j;

	for (i = k; i < nxx - k; i++) {
		for (j = k; j < nyy - k; j++) {
			cell(i, j).exact = linAdvObj.sine(cell(i, j).x - tFinal, \
				cell(i, j).y - tFinal, dx, dy, lagrangeObj.getQuadrature9Weights(), \
				lagrangeObj.getQuadrature9Points());
		}
	}
}

void FVHermite::Scalar2D::printResult() {
	int i, j;
	double err;

	err = 0.0;
	for (j = k; j < nyy - k; j++) {
		for (i = k; i < nxx - k; i++) {
			err += fabs(cell(i, j).u - cell(i, j).exact);
		}
	}
	cout << "Numerical error:" << err / n << endl;
}
