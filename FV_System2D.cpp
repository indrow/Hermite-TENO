#include "FV_System2D.h"

FV::System2D::System2D(const double &xl, const double &xr, const int &nx,
		const double &yb, const double &yt, const int &ny,
		const double &cflnumber, const double &tmax, const double &tInit,
		const double &specificHeatRatio) :
		gamma(specificHeatRatio), gamma_1(specificHeatRatio - 1.0), eulObj(specificHeatRatio) {

	int i, j;

	this->cfl = cflnumber;
	this->tFinal = tmax;
	this->t = tInit;
	dt = 0.0;
	this->nx = nx;
	this->ny = ny;
	dx = (xr - xl) / nx;
	dy = (yt - yb) / ny;
	area = dx * dy;
	k = 6;
	s = (k - 1) / 2;
	//k += 3;
	nxx = k + k + nx;
	nyy = k + k + ny;
	n = nx * ny;
	nt = nxx * nyy;
	twoPointsGQWeights = lagrangeObj.getQuadrature2Weights() / 2.0;
    threePointsGQWeights = lagrangeObj.getGaussLegendre3Weights() / 2.0;
	fourPointsGQWeights = lagrangeObj.getGaussLegendre4Weights() / 2.0;
	sh = 1;
	cout << "Creating domain... ";
	cell.resize(nxx, nyy);

#pragma omp parallel for private(i, j)
	for (j = 0; j < nyy; j++) {
		for (i = 0; i < nxx; i++) {
			cell(i, j).x = ((double) i - (double) k + 0.5) * dx + xl;
			cell(i, j).y = ((double) j - (double) k + 0.5) * dy + yb;
		}
	}

	cout << "done.n";
	cout << "Number of cells: " << n << "n";

	t0 = 0.0, t1 = 0.0, vxm = 0.0, vym = 0.0, hm = 0.0, qm = 0.0, cm = 0.0;
	rcm = 0.0, b1 = 0.0, b2 = 0.0, t2 = 0.0;
}

void FV::System2D::entropyWave() {
	int i, j;
	double rho = 1.0, u = 1.0, v = 1.0, p = 1.0;
	double A = 0.5;

	cout << "Initializing variable... case: Entropy Waves...";
#pragma omp parallel for private(i, j) firstprivate(lagrangeObj, eulObj)
	for (i = 0; i < nxx; i++) {
		for (j = 0; j < nyy; j++) {
			cell(i, j).u = eulObj.entropyWaves(rho, u, v, p, A, cell(i, j).x, cell(i, j).y,
					dx, dy, lagrangeObj.getQuadrature9Weights(), lagrangeObj.getQuadrature9Points());
		}
	}

	cout << "done\n";
#pragma omp parallel for private(i, j) firstprivate(lagrangeObj, eulObj)
	for (i = k; i < nxx - k; i++) {
		for (j = k; j < nyy - k; j++) {
			cell(i, j).exact = eulObj.entropyWaves(rho, u, v, p, A,
					cell(i, j).x - u * tFinal, cell(i, j).y - v * tFinal,
					dx, dy, lagrangeObj.getQuadrature9Weights(), lagrangeObj.getQuadrature9Points());
			cell(i, j).exact = eulObj.primitive(cell(i, j).exact);
		}
	}
	cout << "Exact solutions has been computed...n";
}

void FV::System2D::isentropicVortex() {
    int i, j;

    cout << "Initializing variable... case: Entropy Waves...";
#pragma omp parallel for private(i, j) firstprivate(lagrangeObj, eulObj)
    for (i = 0; i < nxx; i++) {
        for (j = 0; j < nyy; j++) {
            cell(i, j).u = eulObj.isentropicVortex(cell(i, j).x, cell(i, j).y,
                                               dx, dy, lagrangeObj.getQuadrature9Weights(), lagrangeObj.getQuadrature9Points());
        }
    }
    periodicBc();

    cout << " done\n";
#pragma omp parallel for private(i, j) firstprivate(lagrangeObj, eulObj)
    for (i = k; i < nxx - k; i++) {
        for (j = k; j < nyy - k; j++) {
            cell(i, j).exact = eulObj.isentropicVortex(cell(i, j).x, cell(i, j).y,
                                                       dx, dy, lagrangeObj.getQuadrature9Weights(), lagrangeObj.getQuadrature9Points());
            cell(i,j).exact = eulObj.primitive(cell(i,j).exact);
        }
    }
    cout << "Exact solutions has been computed...\n";
}

void FV::System2D::Riemann2D() {
	int i, j;
	cout << "Case: Riemann 2D config. 5 problem.n";
	for (i = 0; i < nxx; i++) {
		for (j = 0; j < nyy; j++) {
			cell(i, j).u = eulObj.conservativeOfTaylorIns2D(cell(i, j).x, cell(i, j).y);
		}
	}
	zeroGradBc();
//	TaylorInstBc();
//	doubleMachBc();
}

void FV::System2D::periodicBc() {
	int i, j, m;

#pragma omp parallel for private(j, m)
	for (j = k; j < ny + k; j++) {
		for (m = 0; m < k; m++) {
			cell(k - 1 - m, j).u = cell(nx + k - 1 - m, j).u;
			cell(nxx - 1 - m, j).u = cell(k + k - 1 - m, j).u;
		}
	}

#pragma omp parallel for private(i, m)
	for (i = 0; i < nxx; i++) {
		for (m = 0; m < k; m++) {
			cell(i, k - 1 - m).u = cell(i, ny + k - 1 - m).u;
			cell(i, nyy - 1 - m).u = cell(i, k + k - 1 - m).u;
		}
	}
}

void FV::System2D::zeroGradBc() {
	int i, j, m;

#pragma omp parallel for private (j, m)
	for (j = k; j < ny + k; ++j) {
		for (m = 0; m < k; ++m) {
			cell(m, j).u = cell(k, j).u;
			cell(nx + k + m, j).u = cell(nx + k - 1, j).u;
		}
	}

#pragma omp parallel for private (i, m)
	for (i = 0; i < nxx; ++i) {
		for (m = 0; m < k; ++m) {
			cell(i, m).u = cell(i, k).u;
			cell(i, ny + k + m).u = cell(i, ny + k - 1).u;
		}
	}
}

void FV::System2D::TaylorInstBc() {
	int i, j, m;

#pragma omp parallel for private (j, m)
	for (j = k; j < nyy - k; ++j) {
		for (m = 0; m < k; ++m) {
			cell(k - 1 - m, j).u = cell(k + m, j).u;
			cell(k - 1 - m, j).u(1) = -cell(k + m, j).u(1);

			cell(nx + k + m, j).u = cell(nx + k - 1 - m, j).u;
			cell(nx + k + m, j).u(1) = -cell(nx + k - 1 - m, j).u(1);
		}
	}

	Vector4d wu;
	wu << 1.0, 0.0, 0.0, 2.5;
	Vector4d wb;
	wb << 2.0, 0.0, 0.0, 1.0;
	Vector4d qu = eulObj.conservative(wu);
	Vector4d qb = eulObj.conservative(wb);

#pragma omp parallel for private (i, m)
	for (i = 0; i < nxx; ++i) {
		for (m = 0; m < k; ++m) {
			cell(i, m).u = qb;
			cell(i, ny + k + m).u = qu;
		}
	}
}

void FV::System2D::doubleMachBc() {
	int i, j, m;
	Vector4d qa;
	qa = eulObj.conservativeOfDoubleMach2D(4.0, 0.0);
	double pi = 4.0 * atan(1.0);
	double ctan = tan(pi / 6.0);
	double hx = 0.0;

#pragma omp parallel for private (i, m, hx) firstprivate(qa)
	for (i = k; i < nxx - k; ++i) {
		if (cell(i, k).x <= 1.0 / 6.0) {
			for (m = 0; m < k; ++m) {
				cell(i, m).u = eulObj.conservativeOfDoubleMach2D(cell(i, m).x, cell(i, m).y);
			}
		} else {
			for (m = 0; m < k; ++m) {
				cell(i, k - 1 - m).u = cell(i, k + m).u;
				cell(i, k - 1 - m).u(2) = -cell(i, k - 1 - m).u(2);
			}
		}

		for (m = 0; m < k; ++m) {
			hx = ctan * (cell(i, nyy - k - 1).y + (1.0 + (double) m) * dy + 20.0 * t) + 1.0 / 6.0;

			if (cell(i, ny - 1).x < hx) {
				cell(i, ny + k + m).u = cell(k, 2 * k).u;
			} else {
				cell(i, ny + k + m).u = cell(nx - 1, 2 * k).u;
			}

		}
	}

#pragma omp parallel for private (j, m) firstprivate(qa)
	for (j = 0; j < nyy; ++j) {
		for (m = 0; m < k; ++m) {
			cell(m, j).u = cell(k, j).u;
			cell(nx + k + m, j).u = cell(nx + k - 1, j).u;
		}
	}
}

void FV::System2D::calcPrimitiveVariable() {
	int i, j;
	eigenValueX.setZero();
	eigenValueY.setZero();

#pragma omp parallel for private(i, j) firstprivate(eulObj)
	for (i = 0; i < nxx; i++) {
		for (j = 0; j < nyy; j++) {
			cell(i, j).pvar = eulObj.primitive(cell(i, j).u);

			/*try {
				if (cell(i, j).pvar(0) <= 0.0)
					throw (-1);
			}
			catch (int x) {
				cout << "Exception occurred, exception number: " << x << endl;
			}

			try {
				if (cell(i, j).pvar(3) <= 0.0)
					throw (-2);
			}
			catch (int x) {
				cout << "Exception occurred, exception number: " << x << endl;
			}*/

			cell(i, j).sqrtRho = sqrt(fabs(cell(i, j).pvar(0)));
			cell(i, j).enthalpy = (cell(i, j).u(3) + cell(i, j).pvar(3)) / cell(i, j).pvar(0);

			//#pragma omp critical
			//			{
			//				eigenValueX = eigenValueX.cwiseMax(eulObj.localWaveSpeed(cell(i, j).u, 1));
			//				eigenValueY = eigenValueY.cwiseMax(eulObj.localWaveSpeed(cell(i, j).u, 2));
			//			}

		}
	}
}

void FV::System2D::eigen2DEuler() {
	int i, j;

#pragma omp parallel for private(i, j, t0, t1, t2, b1, b2, vxm, vym, hm, qm, cm, rcm) firstprivate(eulObj)
	for (i = 0; i < nxx - 1; i++) {
		for (j = 0; j < nyy; j++) {
			t0 = cell(i, j).sqrtRho / (cell(i, j).sqrtRho + cell(i + 1, j).sqrtRho);
			t1 = 1.0 - t0;
			vxm = t0 * cell(i, j).pvar(1) + t1 * cell(i + 1, j).pvar(1);
			vym = t0 * cell(i, j).pvar(2) + t1 * cell(i + 1, j).pvar(2);
			hm = t0 * cell(i, j).enthalpy + t1 * cell(i + 1, j).enthalpy;
			qm = 0.5 * (vxm * vxm + vym * vym);

			/*try {
				if ((hm - qm) < 0.0)
					throw (-3);
			}
			catch (int x) {
				cout << "Exception occurred, exception number: " << x << endl;
				exit(-1);
			}*/

			cm = sqrt(fabs(gamma_1 * (hm - qm)));
			t0 = vxm * cm;
			cell(i, j).rightEigenVectorX = eulObj.rightEigenVectorx(vxm, vym, hm, qm, cm, t0);

			rcm = 1.0 / cm;
			b1 = gamma_1 * rcm * rcm;
			b2 = qm * b1;
			t0 = vxm * rcm;
			t1 = b1 * vxm;
			t2 = 0.5 * b1;

			cell(i, j).leftEigenVectorX = eulObj.leftEigenVectorx(rcm, vym, b1, b2, t0, t1, t2);
		}
	}

	eigenValueX = -aM.setOnes();

#pragma omp parallel for private(i,j)
    for (i = 0; i < nxx - 1; i++) {
        for (j = 0; j < nyy; j++) {
#pragma omp critical
            {
                eigenValueX = eigenValueX.cwiseMax(eulObj.localWaveSpeed(cell(i, j).pvar, 1));
            }
        }
    }

    eigenValueX /= 1.1;
    eigenValueY = -aM.setOnes();

#pragma omp parallel for private(i,j)
    for (i = 0; i < nxx; i++) {
        for (j = 0; j < nyy - 1; j++) {
#pragma omp critical
            {
                eigenValueY = eigenValueY.cwiseMax(eulObj.localWaveSpeed(cell(i, j).pvar, 2));
            }
        }
    }

    eigenValueY /= 1.1;

#pragma omp parallel for private(i, j, t0, t1, t2, b1, b2, vxm, vym, hm, qm, cm, rcm) firstprivate(eulObj)
	for (i = 0; i < nxx; i++) {
		for (j = 0; j < nyy - 1; j++) {
			t0 = cell(i, j).sqrtRho / (cell(i, j).sqrtRho + cell(i, j + 1).sqrtRho);
			t1 = 1.0 - t0;
			vxm = t0 * cell(i, j).pvar(1) + t1 * cell(i, j + 1).pvar(1);
			vym = t0 * cell(i, j).pvar(2) + t1 * cell(i, j + 1).pvar(2);
			hm = t0 * cell(i, j).enthalpy + t1 * cell(i, j + 1).enthalpy;
			qm = 0.5 * (vxm * vxm + vym * vym);

			/*try {
				if ((hm - qm) < 0.0)
					throw (-4);
			}
			catch (int x) {
				cout << "Exception occurred, exception number: " << x << endl;
			}*/

			cm = sqrt(fabs(gamma_1 * (hm - qm)));
			t0 = vym * cm;
			cell(i, j).rightEigenVectorY = eulObj.rightEigenVectory(vxm, vym, hm, qm, cm, t0);

			rcm = 1.0 / cm;
			b1 = gamma_1 * rcm * rcm;
			b2 = qm * b1;
			t0 = vym * rcm;
			t1 = b1 * vym;
			t2 = 0.5 * b1;

			cell(i, j).leftEigenVectorY = eulObj.leftEigenVectory(rcm, vxm, b1, b2, t0, t1, t2);
		}
	}
}

void FV::System2D::rhs() {
	int i, j, m, eq;
//	periodicBc();
	zeroGradBc();
//	TaylorInstBc();
//	doubleMachBc();
	calcPrimitiveVariable();
	eigen2DEuler();

#pragma omp parallel for private(i, j, m, eq, interpolateC,valueB) firstprivate(lagrangeObj)
	for (i = s; i < nxx - s - 1; ++i) {
		for (j = 0; j < nyy; ++j) {
			for (m = 0; m < 6; ++m) {
				interpolateC.col(m) = cell(i, j).leftEigenVectorX * cell(i - 2 + m, j).u;
			}

			for (eq = 0; eq < 4; eq++) {
				valueB = lagrangeObj.teno5(interpolateC.row(eq));
				cell(i, j).uxl(eq) = valueB(0);
				cell(i, j).uxr(eq) = valueB(1);
			}
		}
	}

#pragma omp parallel for private(i, j, m, eq, interpolateC,valueB) firstprivate(lagrangeObj)
	for (j = s; j < nyy - s - 1; ++j) {
		for (i = 0; i < nxx; ++i) {
			for (m = 0; m < 6; ++m) {
				interpolateC.col(m) = cell(i, j).leftEigenVectorY * cell(i, j - 2 + m).u;
			}

			for (eq = 0; eq < 4; eq++) {
				valueB = lagrangeObj.teno5(interpolateC.row(eq));
				cell(i, j).uyl(eq) = valueB(0);
				cell(i, j).uyr(eq) = valueB(1);
			}
		}
	}

#pragma omp parallel for private(i,j,m,eq,ee,ff) firstprivate(lagrangeObj)
	for (i = s + 2; i < nxx - s - 3; ++i) {
		for (j = s + 2; j < nyy - s - 3; ++j) {
			for (m = 0; m < 5; ++m) {
				ee.col(m) = cell(i, j - 2 + m).uxl;
				ff.col(m) = cell(i, j - 2 + m).uxr;
			}
			for (eq = 0; eq < 4; ++eq) {
				cell(i, j).uxgl.row(eq) = lagrangeObj.weno5JS3PointsQuadrature(ee.row(eq));
				cell(i, j).uxgr.row(eq) = lagrangeObj.weno5JS3PointsQuadrature(ff.row(eq));
			}

			for (m = 0; m < 5; ++m) {
				ee.col(m) = cell(i - 2 + m, j).uyl;
				ff.col(m) = cell(i - 2 + m, j).uyr;
			}
			for (eq = 0; eq < 4; ++eq) {
				cell(i, j).uygl.row(eq) = lagrangeObj.weno5JS3PointsQuadrature(ee.row(eq));
				cell(i, j).uygr.row(eq) = lagrangeObj.weno5JS3PointsQuadrature(ff.row(eq));
			}
		}
	}

#pragma omp parallel for private(i, j)
	for (i = k - 1; i < nxx - k; ++i) {
		for (j = k - 1; j < nyy - k; ++j) {
			cell(i, j).transformToConservativeAtGQ();
		}
	}

#pragma omp parallel for private(i,j,m,temporaryFluxUx,temporaryFluxUy,pvarL,pvarR,aM) firstprivate(threePointsGQWeights,eulObj,eigenValueX,eigenValueY)
	for (i = k - 1; i < nxx - k; ++i) {
		for (j = k - 1; j < nyy - k; ++j) {
			cell(i, j).setZeroFlux();

			for (m = 0; m < 3; ++m) {
				pvarL = eulObj.primitive(cell(i, j).uxgl.col(m));
                pvarR = eulObj.primitive(cell(i, j).uxgr.col(m));

				/*try
                {
                    if (pvarL(0) <= 0.0 || pvarL(3) <= 0.0) {
                        cout << pvarL.transpose() << endl;
                        throw (-10);
                    }

                    if (pvarR(0) <= 0.0 || pvarR(3) <= 0.0) {
                        cout << "Cell(" << i << "," << j << ") : " << pvarR.transpose() << endl;
                        throw (-11);
                    }
                }
                catch (const int& x)
                {
                    cout << "Exception occurred, exception number: " << x << endl;
                    exit(-1);
                }*/

//                aM = eulObj.maxLocalWaveSpeed(pvarL.cwiseMax(pvarR), pvarR.cwiseMin(pvarL), 1) / 1.1;
                aM = eulObj.maxLocalWaveSpeed(pvarL, pvarR, 1) / 1.0;
//
                temporaryFluxUx.col(m) =
                        threePointsGQWeights(m) * eulObj.Rusanov(cell(i, j).uxgl.col(m), cell(i, j).uxgr.col(m), pvarL, pvarR,
                                                                cell(i, j).leftEigenVectorX, cell(i, j).rightEigenVectorX, aM, 1);
//                  temporaryFluxUx.col(m) =
//                        threePointsGQWeights(m) * eulObj.HLLC(cell(i, j).uxgl.col(m), cell(i, j).uxgr.col(m), 1);

				pvarL = eulObj.primitive(cell(i, j).uygl.col(m));
                pvarR = eulObj.primitive(cell(i, j).uygr.col(m));

				/*try
                {
                    if (pvarL(0) <= 0.0 || pvarL(3) <= 0.0)
                        throw (-12);
                    if (pvarR(0) <= 0.0 || pvarR(3) <= 0.0)
                        throw (-13);
                }
                catch (const int& x)
                {
                    cout << "Exception occurred, exception number: " << x << endl;
                    exit(-1);
                }*/

//                aM = eulObj.maxLocalWaveSpeed(pvarL.cwiseMax(pvarR), pvarR.cwiseMin(pvarL), 2) / 1.1;
                aM = eulObj.maxLocalWaveSpeed(pvarL, pvarR, 2) / 1.0;
//
                temporaryFluxUy.col(m) =
                        threePointsGQWeights(m) * eulObj.Rusanov(cell(i, j).uygl.col(m), cell(i, j).uygr.col(m), pvarL, pvarR,
                    cell(i, j).leftEigenVectorY, cell(i, j).rightEigenVectorY, aM, 2);
//                temporaryFluxUy.col(m) =
//                        threePointsGQWeights(m) * eulObj.HLLC(cell(i, j).uygl.col(m), cell(i, j).uygr.col(m), 2);
			}
            cell(i, j).fx = temporaryFluxUx.rowwise().sum();
            cell(i, j).fy = temporaryFluxUy.rowwise().sum();

		}
	}

//#pragma omp parallel for private(i, j)
//    for (i = k - 1; i < nxx - k; ++i) {
//        for (j = k - 1; j < nyy - k; ++j) {
//            cell(i, j).transformToConservative();
//        }
//    }
//#pragma omp parallel for private(i,j,m,temporaryFluxUx,temporaryFluxUy) firstprivate(pvarL,pvarR,aM,fourPointsGQWeights,eulObj)
//    for (i = k - 1; i < nxx - k; ++i) {
//        for (j = k - 1; j < nyy - k; ++j) {
//
//            pvarL = eulObj.primitive(cell(i, j).uxl);
//            pvarR = eulObj.primitive(cell(i, j).uxr);
//
//                /*try
//                {
//                    if (pvarL(0) <= 0.0 || pvarL(3) <= 0.0) {
//                        cout << pvarL.transpose() << endl;
//                        throw (-10);
//                    }
//
//                    if (pvarR(0) <= 0.0 || pvarR(3) <= 0.0) {
//                        cout << "Cell(" << i << "," << j << ") : " << pvarR.transpose() << endl;
//                        throw (-11);
//                    }
//                }
//                catch (const int& x)
//                {
//                    cout << "Exception occurred, exception number: " << x << endl;
//                    exit(-1);
//                }*/
//
//            aM = eulObj.maxLocalWaveSpeed(pvarL, pvarR, 1) / 1.1;
//
//            cell(i, j).fx =
//                    eulObj.Rusanov(cell(i, j).uxl, cell(i, j).uxr, pvarL, pvarR,
//                                   cell(i, j).leftEigenVectorX, cell(i, j).rightEigenVectorX, aM, 1);
//                pvarL = eulObj.primitive(cell(i, j).uyl);
//                pvarR = eulObj.primitive(cell(i, j).uyr);
//
//                /*try
//                {
//                    if (pvarL(0) <= 0.0 || pvarL(3) <= 0.0)
//                        throw (-12);
//                    if (pvarR(0) <= 0.0 || pvarR(3) <= 0.0)
//                        throw (-13);
//                }
//                catch (const int& x)
//                {
//                    cout << "Exception occurred, exception number: " << x << endl;
//                    exit(-1);
//                }*/
//
//                aM = eulObj.maxLocalWaveSpeed(pvarL, pvarR, 2) / 1.1;
//
//            cell(i, j).fy =
//                    eulObj.Rusanov(cell(i, j).uyl, cell(i, j).uyr, pvarL, pvarR,
//                                                                cell(i, j).leftEigenVectorY, cell(i, j).rightEigenVectorY, aM, 2);
////                    eulObj.HLLC(cell(i, j).uygl.col(m), cell(i, j).uygr.col(m), 2);
//
//        }
//    }


#pragma omp parallel for private(i, j, source)
	for (i = k; i < nxx - k; i++) {
		for (j = k; j < nyy - k; j++) {
			source.setZero();
//			source(2) = cell(i, j).u(0);
//			source(3) = cell(i, j).u(2);

			cell(i, j).rhsu = ((cell(i - 1, j).fx - cell(i, j).fx) * dy +
					(cell(i, j - 1).fy - cell(i, j).fy) * dx) / area + source;
		}
	}
}

void FV::System2D::RungeKutta3() {
	int i, j;
	int stage;

	stage = 1;
	rhs();
#pragma omp parallel for private(i, j)
	for (i = k; i < nxx - k; i++) {
		for (j = k; j < nyy - k; j++) {
			cell(i, j).save();
			cell(i, j).u1 = ssprkObj.RungeKuttaIII<Vector4d>(cell(i, j).u0,
					cell(i, j).u0, cell(i, j).rhsu, dt, stage);
			cell(i, j).update1();
		}
	}

	stage += 1;
	rhs();
#pragma omp parallel for private(i, j)
	for (i = k; i < nxx - k; i++) {
		for (j = k; j < nyy - k; j++) {
			cell(i, j).u2 = ssprkObj.RungeKuttaIII<Vector4d>(cell(i, j).u0,
					cell(i, j).u1, cell(i, j).rhsu, dt, stage);
			cell(i, j).update2();
		}
	}

	stage += 1;
	rhs();
#pragma omp parallel for private(i, j)
	for (i = k; i < nxx - k; i++) {
		for (j = k; j < nyy - k; j++) {
			cell(i, j).u = ssprkObj.RungeKuttaIII<Vector4d>(cell(i, j).u0,
					cell(i, j).u2, cell(i, j).rhsu, dt, stage);
		}
	}
	calcPrimitiveVariable();
}

void FV::System2D::timeStepping() {
	int i, j;
	double maxSpeedx, maxSpeedy;
	periodicBc();
	calcPrimitiveVariable();

//	bool a1 = true, a2 = true;

	while (t < tFinal) {
		maxSpeedx = 0.0;
		maxSpeedy = 0.0;
#pragma omp parallel for private(i, j) reduction(max:maxSpeedx, maxSpeedy)
		for (i = k; i < nxx - k; i++) {
			for (j = k; j < nyy - k; j++) {
				maxSpeedx = fmax(maxSpeedx, eulObj.maxSpeed(cell(i, j).pvar, 1));
				maxSpeedy = fmax(maxSpeedy, eulObj.maxSpeed(cell(i, j).pvar, 2));
			}
		}

		dt = cfl * pow(fmin(dx, dy), 3.0 / 3.0) / fmax(maxSpeedx, maxSpeedy);

//		if (t + dt > 100.0 && a1) {
//		    dt = 100.0 - t;
//		    write(100);
//		    a1 = false;
//		} else if (t + dt > tFinal && a2) {
//		    dt = tFinal - t;
//		    write(1000);
//		    a2 = false;
//		}

		dt = t + dt > tFinal ? tFinal - t : dt;

		RungeKutta3();
		t += dt;
		cout << "tnum: " << t << endl;
	}
}

void FV::System2D::printResult() {
	int i, j;
	double err;

	err = 0.0;
//#pragma omp parallel for private(i, j) reduction(+:err)
	for (i = k; i < nxx - k; i++) {
		for (j = k; j < nyy - k; j++) {
			err += fabs(cell(i, j).u(0) - cell(i, j).exact(0));
		}
	}
	cout << "Numerical error:" << err / n << endl;
	write(0);
}

void FV::System2D::write(const int &filenum) {

	std::string filename =
			"../hasil/2D_RTI_U0_" + std::to_string(nx) + "_" + std::to_string(ny) +
			"_" + "Final.dat";

	std::ofstream outfile;
	int i, j;
	const char separator = ' ';
	const int strwidth = 20;
	const int numwidth = 20;
	const int pn = 8;

	outfile.open(filename);
	if (outfile.good()) {
		outfile << "t: " << t << "\n";
		outfile << std::right << std::setw(strwidth) << std::setfill(separator) << "x";
		outfile << std::right << std::setw(strwidth) << std::setfill(separator) << "y";
		outfile << std::right << std::setw(strwidth) << std::setfill(separator) << "Density";
		outfile << std::right << std::setw(strwidth) << std::setfill(separator) << "Ux";
		outfile << std::right << std::setw(strwidth) << std::setfill(separator) << "Uy";
		outfile << std::right << std::setw(strwidth) << std::setfill(separator) << "Pressure";
		outfile << "\n";
	} else {
		exit(1);
	}

	for (j = k; j < nyy - k; ++j) {
		for (i = k; i < nxx - k; ++i) {
			outfile << std::right << std::setw(numwidth) << std::setfill(separator) << std::setprecision(pn)
					<< cell(i, j).x;
			outfile << std::right << std::setw(numwidth) << std::setfill(separator) << std::setprecision(pn)
					<< cell(i, j).y;
			outfile << std::right << std::setw(numwidth) << std::setfill(separator) << std::setprecision(pn)
					<< cell(i, j).u(0);
			outfile << std::right << std::setw(numwidth) << std::setfill(separator) << std::setprecision(pn)
					<< cell(i, j).u(1);
			outfile << std::right << std::setw(numwidth) << std::setfill(separator) << std::setprecision(pn)
					<< cell(i, j).u(2);
			outfile << std::right << std::setw(numwidth) << std::setfill(separator) << std::setprecision(pn)
					<< cell(i, j).u(3);
			outfile << "\n";
		}
	}

	std::cout << "n Solution is written in " << filename << "\n";
	outfile.close();
}
