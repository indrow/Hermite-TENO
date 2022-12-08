#include "FV_Scalar2D.h"

FV::Scalar2D::Scalar2D(const double& xl, const double& xr, const int& nx,
		const double& yb, const double& yt, const int& ny,
		const double& cflnumber, const double& tmax, const double& tInit) {
	int i, j;

	this->cfl = cflnumber;this->tFinal = tmax; this->t = tInit; dt = 0.0;
	this->nx = nx; this->ny = ny;
	dx = (xr - xl) / nx;
	dy = (yt - yb) / ny;
	area = dx * dy;
	k = 5; s = (k - 1) / 2;
	nxx = k + k + nx;
	nyy = k + k + ny;
	n = nx * ny;
	nt = nxx * nyy;
	hermiteObj.setGridSize(dx);
	twoPointsGQWeights = lagrangeObj.getQuadrature2Weights() / 2.0;
	fourPointsGLegendreQWeights = lagrangeObj.getGaussLobatto4Weights() / 2.0;

	cout << "Creating domain... ";
	cell.resize(nxx, nyy);
	for (j = 0; j < nyy; j++) {
		for (i = 0; i < nxx; i++) {
			cell(i, j).x = ((double)i - (double)k + 0.5) * dx + xl;
			cell(i, j).y = ((double)j - (double)k + 0.5) * dy + yb;
		}
	}

	cout << "done.\n";
	cout << "Number of cells: " << n << "n";


}

void FV::Scalar2D::initialize() {
	int i, j;
	cout << "Initializing variable...";
	for (i = k; i < nxx - k; i++) {
		for (j = k; j < nyy - k; j++) {
			cell(i, j).u = linAdvObj.sine(cell(i, j).x, cell(i, j).y, dx, dy,
					lagrangeObj.getQuadrature9Weights(), lagrangeObj.getQuadrature9Points());
		}
	}
	cout << "done.\n";
}

void FV::Scalar2D::periodicBc() {
	int i, j, m;

	for (j = k; j < ny + k; j++) {
		for (m = 0; m < k; m++) {
			cell((int)(k - 1 - m), j).u = cell(nx + k - 1 - m, j).u;
			cell(nxx - 1 - m, j).u = cell(k + k - 1 - m, j).u;
		}
	}

	for (i = 0; i < nxx; i++) {
		for (m = 0; m < k; m++) {
			cell(i, k - 1 - m).u = cell(i, ny + k - 1 - m).u;
			cell(i, nyy - 1 - m).u = cell(i, k + k - 1 - m).u;
		}
	}
}

void FV::Scalar2D::rhs() {
	int i, j, m;
	periodicBc();

	// Find the interpolated points in x dir
	for (i = s; i < nxx - s - 1; i++) {
		for (j = 0; j < nyy; j++) {
			for (m = 0; m < 6; m++) {
				interpolateC(m) = cell(i - s + m, j).u;
			}
			cell(i, j).ux = lagrangeObj.teno5(interpolateC);
		}
	}

	// Find the interpolated points in y dir
	for (i = 0; i < nxx; i++) {
		for (j = s; j < nyy - s - 1; j++) {
			for (m = 0; m < 6; m++) {
				interpolateC(m) = cell(i, j - s + m).u;
			}
			cell(i, j).uy = lagrangeObj.teno5(interpolateC);
		}
	}

	GaussLegendre2Points();

	for (i = k; i < nxx - k; i++) {
		for (j = k; j < nyy - k; j++) {
			cell(i, j).rhs = ((cell(i - 1, j).fx - cell(i, j).fx) * dy +
					(cell(i, j - 1).fy - cell(i, j).fy) * dx) / area;
		}
	}
}

void FV::Scalar2D::GaussLegendre2Points() {
	int i, j, m;

	for (i = k - 1; i < nxx - k; i++) {
		for (j = k - 1; j < nyy - k; j++) {
			for (m = 0; m < 5; m++) {
				ee(m) = cell(i, j - s + m).ux(0);
				ff(m) = cell(i, j - s + m).ux(1);
			}
			cell(i, j).uxg2l = lagrangeObj.teno5th2PointsQuadrature(ee);
			cell(i, j).uxg2r = lagrangeObj.teno5th2PointsQuadrature(ff);

			for (m = 0; m < 5; m++) {
				ee(m) = cell(i - s + m, j).uy(0);
				ff(m) = cell(i - s + m, j).uy(1);
			}
			cell(i, j).uyg2l = lagrangeObj.teno5th2PointsQuadrature(ee);
			cell(i, j).uyg2r = lagrangeObj.teno5th2PointsQuadrature(ff);
		}
	}

	for (i = k - 1; i < nxx - k; i++) {
		for (j = k - 1; j < nyy - k; j++) {
			cell(i, j).fx = 0.0;
			cell(i, j).fy = 0.0;
			for (m = 0; m < 2; m++) {
				cell(i, j).fx += twoPointsGQWeights(m) *
						linAdvObj.Rusanov(cell(i, j).uxg2l(m),
								cell(i, j).uxg2r(m));
				cell(i, j).fy += twoPointsGQWeights(m) *
						linAdvObj.Rusanov(cell(i, j).uyg2l(m),
								cell(i, j).uyg2r(m));
			}
		}
	}
}

void FV::Scalar2D::GaussLobatto4Points() {
	int i, j, m;

	for (i = k - 1; i < nxx - k; i++) {
		for (j = k - 1; j < nyy - k; j++) {
			for (m = 0; m < 5; m++) {
				ee(m) = cell(i, j - s + m).ux(0);
				ff(m) = cell(i, j - s + m).ux(1);
			}
			cell(i, j).uxg4l = lagrangeObj.weno5JS4PointsGaussLobatto(ee);
			cell(i, j).uxg4r = lagrangeObj.weno5JS4PointsGaussLobatto(ff);

			for (m = 0; m < 5; m++) {
				ee(m) = cell(i - s + m, j).uy(0);
				ff(m) = cell(i - s + m, j).uy(1);
			}
			cell(i, j).uyg4l = lagrangeObj.weno5JS4PointsGaussLobatto(ee);
			cell(i, j).uyg4r = lagrangeObj.weno5JS4PointsGaussLobatto(ff);
		}
	}

	for (i = k - 1; i < nxx - k; i++) {
		for (j = k - 1; j < nyy - k; j++) {
			cell(i, j).fx = 0.0;
			cell(i, j).fy = 0.0;
			for (m = 0; m < 4; m++) {
				cell(i, j).fx += fourPointsGLegendreQWeights(m) *
						linAdvObj.Rusanov(cell(i, j).uxg4l(m),
								cell(i, j).uxg4r(m));
				cell(i, j).fy += fourPointsGLegendreQWeights(m) *
						linAdvObj.Rusanov(cell(i, j).uyg4l(m),
								cell(i, j).uyg4r(m));
			}
		}
	}
}

void FV::Scalar2D::RungeKutta3() {
	int i, j;
	int stage;

	stage = 1;
	rhs();
	for (i = k; i < nxx - k; i++) {
		for (j = k; j < nyy - k; j++) {
			cell(i, j).save();
			cell(i, j).u1 = ssprkObj.RungeKuttaIII<double>(cell(i, j).u0,
					cell(i, j).u0, cell(i, j).rhs, dt, stage);
			cell(i, j).update1();
		}
	}

	stage += 1;
	rhs();
	for (i = k; i < nxx - k; i++) {
		for (j = k; j < nyy - k; j++) {
			cell(i, j).u2 = ssprkObj.RungeKuttaIII<double>(cell(i, j).u0,
					cell(i, j).u1, cell(i, j).rhs, dt, stage);
			cell(i, j).update2();
		}
	}

	stage += 1;
	rhs();
	for (i = k; i < nxx - k; i++) {
		for (j = k; j < nyy - k; j++) {
			cell(i, j).u = ssprkObj.RungeKuttaIII<double>(cell(i, j).u0,
					cell(i, j).u2, cell(i, j).rhs, dt, stage);
		}
	}
}

void FV::Scalar2D::timeStepping() {
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

void FV::Scalar2D::Exact() {
	int i, j;

	for (i = k; i < nxx - k; i++) {
		for (j = k; j < nyy - k; j++) {
			cell(i, j).exact = linAdvObj.sine(cell(i, j).x - tFinal,
					cell(i, j).y - tFinal, dx, dy, lagrangeObj.getQuadrature9Weights(),
					lagrangeObj.getQuadrature9Points());
		}
	}
}

void FV::Scalar2D::printResult() {
	int i, j;
	double err;

	err = 0.0;
	for (j = k; j < nyy - k; j++) {
		for (i = k; i < nxx - k; i++) {
			err += fabs(cell(i, j).u - cell(i, j).exact);
		}
	}
	cout << "L1 error:" << err / n << endl;
}
