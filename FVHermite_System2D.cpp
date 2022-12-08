#include "FVHermite_System2D.h"

FVHermite::System2D::System2D(const double &xl, const double &xr, const int &nx,
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
    k = 5;
    s = (k - 1) / 2;
    sh = 1;
    k += 3;
    nxx = k + k + nx;
    nyy = k + k + ny;
    n = nx * ny;
    nt = nxx * nyy;
    hermiteObj.setGridSize(dx);
    twoPointsGQWeights = lagrangeObj.getQuadrature2Weights() / 2.0;
    threePointsGQWeights = lagrangeObj.getGaussLegendre3Weights() / 2.0;
    fourPointsGQWeights = lagrangeObj.getGaussLegendre4Weights() / 2.0;

    cout << "Creating domain... ";
    cell.resize(nxx, nyy);

#pragma omp parallel for private(i, j)
    for (j = 0; j < nyy; j++) {
        for (i = 0; i < nxx; i++) {
            cell(i, j).x = ((double) i - (double) k + 0.5) * dx + xl;
            cell(i, j).y = ((double) j - (double) k + 0.5) * dy + yb;
        }
    }

    cout << "done.\n";
    cout << "Number of cells: " << n << "\n";

    t0 = 0.0, t1 = 0.0, vxm = 0.0, vym = 0.0, hm = 0.0, qm = 0.0, cm = 0.0;
    rcm = 0.0, b1 = 0.0, b2 = 0.0, t2 = 0.0;
    source.setZero();
}

void FVHermite::System2D::entropyWave() {
    int i, j;
    double rho = 1.0, u = 1, v = 1, p = 1.0;
    double A = 0.5;

    cout << "Initializing variable... case: Entropy Waves...";
#pragma omp parallel for private(i, j) firstprivate(lagrangeObj, eulObj)
    for (i = 0; i < nxx; i++) {
        for (j = 0; j < nyy; j++) {
            cell(i, j).u = eulObj.entropyWaves(rho, u, v, p, A, cell(i, j).x, cell(i, j).y,
                                               dx, dy, lagrangeObj.getQuadrature9Weights(),
                                               lagrangeObj.getQuadrature9Points());

            cell(i, j).v = eulObj.entropyWavesdx(rho, u, v, p, A, cell(i, j).x, cell(i, j).y,
                                                 dx, dy, lagrangeObj.getQuadrature9Weights(),
                                                 lagrangeObj.getQuadrature9Points());

            cell(i, j).w = eulObj.entropyWavesdy(rho, u, v, p, A, cell(i, j).x, cell(i, j).y,
                                                 dx, dy, lagrangeObj.getQuadrature9Weights(),
                                                 lagrangeObj.getQuadrature9Points());
        }
    }

    cout << "done\n";
#pragma omp parallel for private(i, j) firstprivate(lagrangeObj, eulObj)
    for (i = k; i < nxx - k; i++) {
        for (j = k; j < nyy - k; j++) {
            cell(i, j).exact = eulObj.entropyWaves(rho, u, v, p, A,
                                                   cell(i, j).x - u * tFinal, cell(i, j).y - v * tFinal,
                                                   dx, dy, lagrangeObj.getQuadrature9Weights(),
                                                   lagrangeObj.getQuadrature9Points());
            cell(i, j).exact = eulObj.primitive(cell(i, j).exact);
        }
    }
    cout << "Exact solutions has been computed...\n";
}

void FVHermite::System2D::Riemann2D() {
    int i, j;
    cout << "Case: Riemann 2D config. 5 problem.\n";
    for (i = 0; i < nxx; i++) {
        for (j = 0; j < nyy; j++) {
            cell(i, j).u = eulObj.conservativeOfRiemann2D(cell(i, j).x, cell(i, j).y);

            cell(i, j).v = (eulObj.conservativeOfRiemann2D(cell(i, j).x + dx / 2.0, cell(i, j).y) -
                            eulObj.conservativeOfRiemann2D(cell(i, j).x - dx / 2.0, cell(i, j).y)) / dx;

            cell(i, j).w = (eulObj.conservativeOfRiemann2D(cell(i, j).x, cell(i, j).y + dy / 2.0) -
                            eulObj.conservativeOfRiemann2D(cell(i, j).x, cell(i, j).y - dy / 2.0)) / dy;

        }
    }
    zeroGradBc();
//    TaylorInstBc();
//	doubleMachBc();
}

void FVHermite::System2D::periodicBc() {
    int i, j, m;

#pragma omp parallel for private(j, m)
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

#pragma omp parallel for private(i, m)
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

void FVHermite::System2D::zeroGradBc() {
    int i, j, m;

#pragma omp parallel for private (j, m)
    for (j = k; j < ny + k; ++j) {
        for (m = 0; m < k; ++m) {
            cell(m, j).u = cell(k, j).u;
            cell(m, j).v = cell(k, j).v;
            cell(m, j).w = cell(k, j).w;

            cell(nx + k + m, j).u = cell(nx + k - 1, j).u;
            cell(nx + k + m, j).v = cell(nx + k - 1, j).v;
            cell(nx + k + m, j).w = cell(nx + k - 1, j).w;
        }
    }

#pragma omp parallel for private (i, m)
    for (i = 0; i < nxx; ++i) {
        for (m = 0; m < k; ++m) {
            cell(i, m).u = cell(i, k).u;
            cell(i, m).v = cell(i, k).v;
            cell(i, m).w = cell(i, k).w;

            cell(i, ny + k + m).u = cell(i, ny + k - 1).u;
            cell(i, ny + k + m).v = cell(i, ny + k - 1).v;
            cell(i, ny + k + m).w = cell(i, ny + k - 1).w;
        }
    }
}

void FVHermite::System2D::TaylorInstBc() {
    int i, j, m;

#pragma omp parallel for private (j, m)
    for (j = k; j < nyy - k; ++j) {
        for (m = 0; m < k; ++m) {
            cell(k - 1 - m, j).u = cell(k + m, j).u;
            cell(k - 1 - m, j).v = -cell(k + m, j).v;
            cell(k - 1 - m, j).w = cell(k + m, j).w;

            cell(k - 1 - m, j).u(1) = -cell(k + m, j).u(1);
            cell(k - 1 - m, j).v(1) = -cell(k + m, j).v(1);// du/dx
            cell(k - 1 - m, j).w(1) = -cell(k + m, j).w(1);// du/dy


//            cell(k - 1 - m, j).w(2) = -cell(k + m, j).w(2); // dv/dy

            cell(nx + k + m, j).u = cell(nx + k - 1 - m, j).u;
            cell(nx + k + m, j).v = -cell(nx + k - 1 - m, j).v;
            cell(nx + k + m, j).w = cell(nx + k - 1 - m, j).w;

            cell(nx + k + m, j).u(1) = -cell(nx + k - 1 - m, j).u(1);
            cell(nx + k + m, j).v(1) = -cell(nx + k - 1 - m, j).v(1);
            cell(nx + k + m, j).w(1) = -cell(nx + k - 1 - m, j).w(1);

//            cell(nx + k + m, j).w(2) = -cell(nx + k - 1 - m, j).w(2);
        }
    }

    Vector4d wu;
    wu << 1.0, 0.0, 0.0, 2.5;
    Vector4d wb;
    wb << 2.0, 0.0, 0.0, 1.0;
    Vector4d qu = eulObj.conservative(wu);
    Vector4d qb = eulObj.conservative(wb);
    Vector4d v, w;
    v.setZero();
    w.setZero();

#pragma omp parallel for private (i, m) firstprivate(qu, qb, v, w)
    for (i = 0; i < nxx; ++i) {
        for (m = 0; m < k; ++m) {
            cell(i, m).u = qb;
            cell(i, k - 1 - m).v = v;
            cell(i, k - 1 - m).w = w;

            cell(i, ny + k + m).u = qu;
            cell(i, ny + k + m).v = v;
            cell(i, ny + k + m).w = w;
        }
    }
}

void FVHermite::System2D::doubleMachBc() {
    int i, j, m;
    Vector4d qa(4), va(4), wa(4);
    qa = eulObj.conservativeOfDoubleMach2D(4.0, 0.0);
    va.setZero();
    wa.setZero();
    double pi = 4.0 * atan(1.0);
    double ctan = tan(pi / 6.0);
    double hx = 0.0;

#pragma omp parallel for private (i, m, hx) firstprivate(qa, va, wa)
    for (i = k; i < nxx - k; ++i) {
        if (cell(i, k).x <= 1.0 / 6.0) {
            for (m = 0; m < k; ++m) {
                cell(i, m).u = qa;
                cell(i, m).v = va;
                cell(i, m).w = wa;
            }
        } else {
            for (m = 0; m < k; ++m) {
                cell(i, k - 1 - m).u = cell(i, k + m).u;
                cell(i, k - 1 - m).v = cell(i, k + m).v;
                cell(i, k - 1 - m).w = -cell(i, k + m).w;
                cell(i, k - 1 - m).u(2) = -cell(i, k + m).u(2);
                cell(i, k - 1 - m).v(2) = -cell(i, k + m).v(2);
                cell(i, k - 1 - m).w(2) = -cell(i, k + m).w(2);
            }
        }
        for (m = 0; m < k; ++m) {
            hx = ctan * (cell(i, nyy - k - 1).y + (1.0 + (double) m) * dy + 20.0 * t) + 1.0 / 6.0;

            if (cell(i, ny - 1).x < hx) {
                cell(i, ny + k + m).u = cell(k, 2 * k).u;
                cell(i, ny + k + m).v = cell(k, 2 * k).v;
                cell(i, ny + k + m).w = cell(k, 2 * k).w;
            } else {
                cell(i, ny + k + m).u = cell(nx - 1, 2 * k).u;
                cell(i, ny + k + m).v = cell(nx - 1, 2 * k).v;
                cell(i, ny + k + m).w = cell(nx - 1, 2 * k).w;
            }
        }
    }
    qa = eulObj.conservativeOfDoubleMach2D(0.0, 0.0);

#pragma omp parallel for private (j, m) firstprivate(qa, va, wa)
    for (j = 0; j < nyy; ++j) {
        for (m = 0; m < k; ++m) {
            cell(m, j).u = qa;
            cell(m, j).v = va;
            cell(m, j).w = wa;
            cell(nx + k + m, j).u = cell(nx + k - 1, j).u;
            cell(nx + k + m, j).v = cell(nx + k - 1, j).v;
            cell(nx + k + m, j).w = cell(nx + k - 1, j).w;
        }
    }
}

void FVHermite::System2D::calcPrimitiveVariable() {
    int i, j;

#pragma omp parallel for private(i, j) firstprivate(eulObj)
    for (i = 0; i < nxx; i++) {
        for (j = 0; j < nyy; j++) {
            cell(i, j).pvar = eulObj.primitive(cell(i, j).u);

            /*try {
                if (cell(i, j).pvar(0) <= 0.0)
                    throw (-1);
                else if (isnan(cell(i, j).pvar(0))) {
                    throw (-101);
                }

            }
            catch (int x) {
                cout << "Exception occurred, exception number: " << x << endl;
                exit(-1);
            }

            try {
                if (cell(i, j).pvar(3) <= 0.0)
                    throw (-2);
                else if (isnan(cell(i, j).pvar(3))) {
                    throw (-102);
                }
            }
            catch (int x) {
                cout << "Exception occurred, exception number: " << x << endl;
                exit(-1);
            }*/

            cell(i, j).sqrtRho = sqrt(fabs(cell(i, j).pvar(0)));
            cell(i, j).enthalpy = (cell(i, j).u(3) + cell(i, j).pvar(3)) / cell(i, j).pvar(0);
        }
    }
}

void FVHermite::System2D::eigen2DEuler() {
    int i, j;

#pragma omp parallel for private(i, j, t0, t1, t2, b1, b2, vxm, vym, hm, qm, cm, rcm) \
        firstprivate(eulObj)
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

#pragma omp parallel for private(i, j, t0, t1, t2, b1, b2, vxm, vym, hm, qm, cm, rcm) \
        firstprivate(eulObj)
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

void FVHermite::System2D::limitDerivatives() {
    int i, j, eq, m;
    Eigen::Matrix<double, 4, 3> uu, vv, ww;

#pragma omp parallel for private(i, j, m, eq,uu, vv, ww) firstprivate(hermiteObj)
    for (i = k; i < nxx - k + 1; ++i) {
        for (j = 0; j < nyy; ++j) {
            for (m = 0; m < 3; ++m) {
                uu.col(m) = cell(i - 1 + m, j).u;
                vv.col(m) = cell(i - 1 + m, j).v;
                //ww.col(m) = cell(i - 1 + m, j).w;

            }
            for (eq = 0; eq < 4; eq++) {
                cell(i, j).v(eq) = hermiteObj.hermiteLimiter(uu.row(eq), vv.row(eq));
                //cell(i,j).w(eq) = hermiteObj.hermiteLimiter(uu.row(eq), ww.row(eq));
            }
        }
    }

#pragma omp parallel for private(i, j, m, eq,uu, vv, ww) firstprivate(hermiteObj)
    for (i = 0; i < nxx; ++i) {
        for (j = k; j < nyy - k + 1; ++j) {
            for (m = 0; m < 3; ++m) {
                uu.col(m) = cell(i, j - 1 + m).u;
                //vv.col(m) = cell(i, j - 1 + m).v;
                ww.col(m) = cell(i, j - 1 + m).w;
            }
            for (eq = 0; eq < 4; eq++) {
                //cell(i,j).v(eq) = hermiteObj.hermiteLimiter(uu.row(eq), vv.row(eq));
                cell(i, j).w(eq) = hermiteObj.hermiteLimiter(uu.row(eq), ww.row(eq));
            }
        }
    }
}

void FVHermite::System2D::rhs() {
    int i, j, m, eq;
//    TaylorInstBc();
    zeroGradBc();
//    periodicBc();
//    doubleMachBc();
    limitDerivatives();
//    periodicBc();
    zeroGradBc();
//    TaylorInstBc();
//	doubleMachBc();
    calcPrimitiveVariable();
    eigen2DEuler();

#pragma omp parallel for private(i, j, m, eq,interpolateA, interpolateB, interpolateC, \
        valueA, valueB) firstprivate(lagrangeObj, hermiteObj)
    for (i = s; i < nxx - s - 1; ++i) {
        for (j = 0; j < nyy; ++j) {
            for (m = 0; m < 4; ++m) {
                interpolateA.col(m) = cell(i, j).leftEigenVectorX * cell(i - 1 + m, j).u;
                interpolateB.col(m) = cell(i, j).leftEigenVectorX * cell(i - 1 + m, j).v;
            }

            for (m = 0; m < 6; ++m) {
                interpolateC.col(m) = cell(i, j).leftEigenVectorX * cell(i - 2 + m, j).w;
            }

            for (eq = 0; eq < 4; eq++) {
                valueA = hermiteObj.hteno(interpolateA.row(eq), interpolateB.row(eq));
                cell(i, j).uxl(eq) = valueA(0);
                cell(i, j).uxr(eq) = valueA(1);
                cell(i, j).vxl(eq) = valueA(2);
                cell(i, j).vxr(eq) = valueA(3);

                valueB = lagrangeObj.teno5(interpolateC.row(eq));
                cell(i, j).wxl(eq) = valueB(0);
                cell(i, j).wxr(eq) = valueB(1);
            }
        }
    }

#pragma omp parallel for private(i, j, m, eq, interpolateA, interpolateB, interpolateC, \
        valueA, valueB) firstprivate(lagrangeObj, hermiteObj)
    for (j = s; j < nyy - s - 1; ++j) {
        for (i = 0; i < nxx; ++i) {
            for (m = 0; m < 4; ++m) {
                interpolateA.col(m) = cell(i, j).leftEigenVectorY * cell(i, j - 1 + m).u;
                interpolateB.col(m) = cell(i, j).leftEigenVectorY * cell(i, j - 1 + m).w;
            }

            for (m = 0; m < 6; ++m) {
                interpolateC.col(m) = cell(i, j).leftEigenVectorY * cell(i, j - 2 + m).v;
            }

            for (eq = 0; eq < 4; eq++) {
                valueA = hermiteObj.hteno(interpolateA.row(eq), interpolateB.row(eq));
                cell(i, j).uyl(eq) = valueA(0);
                cell(i, j).uyr(eq) = valueA(1);
                cell(i, j).wyl(eq) = valueA(2);
                cell(i, j).wyr(eq) = valueA(3);

                valueB = lagrangeObj.teno5(interpolateC.row(eq));
                cell(i, j).vyl(eq) = valueB(0);
                cell(i, j).vyr(eq) = valueB(1);
            }
        }
    }

#pragma omp parallel for private(i, j, m, eq,aa, bb, cc, dd, ee, ff, valueC) firstprivate(lagrangeObj, hermiteObj)
    for (i = s + 2; i < nxx - s - 3; ++i) {
        for (j = s + 2; j < nyy - s - 3; ++j) {
            for (m = 0; m < 5; ++m) {
                aa.col(m) = cell(i, j - 2 + m).uxl;
                bb.col(m) = cell(i, j - 2 + m).vxl;
                cc.col(m) = cell(i, j - 2 + m).uxr;
                dd.col(m) = cell(i, j - 2 + m).vxr;
            }
            for (eq = 0; eq < 4; ++eq) {
                cell(i, j).uxgl.row(eq) = lagrangeObj.weno5JS3PointsQuadrature(aa.row(eq));
                cell(i, j).vxgl.row(eq) = lagrangeObj.weno5JS3PointsQuadrature(bb.row(eq));
                cell(i, j).uxgr.row(eq) = lagrangeObj.weno5JS3PointsQuadrature(cc.row(eq));
                cell(i, j).vxgr.row(eq) = lagrangeObj.weno5JS3PointsQuadrature(dd.row(eq));
            }
            for (m = 0; m < 5; ++m) {
                ee.col(m) = cell(i, j - 2 + m).wxl;
                ff.col(m) = cell(i, j - 2 + m).wxr;
            }
            for (eq = 0; eq < 4; ++eq) {
                cell(i, j).wxgl.row(eq) = lagrangeObj.weno5JS3PointsQuadrature(ee.row(eq));
                cell(i, j).wxgr.row(eq) = lagrangeObj.weno5JS3PointsQuadrature(ff.row(eq));
            }

            for (m = 0; m < 5; ++m) {
                aa.col(m) = cell(i - 2 + m, j).uyl;
                bb.col(m) = cell(i - 2 + m, j).wyl;
                cc.col(m) = cell(i - 2 + m, j).uyr;
                dd.col(m) = cell(i - 2 + m, j).wyr;
            }
            for (eq = 0; eq < 4; ++eq) {
                cell(i, j).uygl.row(eq) = lagrangeObj.weno5JS3PointsQuadrature(aa.row(eq));
                cell(i, j).wygl.row(eq) = lagrangeObj.weno5JS3PointsQuadrature(bb.row(eq));
                cell(i, j).uygr.row(eq) = lagrangeObj.weno5JS3PointsQuadrature(cc.row(eq));
                cell(i, j).wygr.row(eq) = lagrangeObj.weno5JS3PointsQuadrature(dd.row(eq));
            }
            for (m = 0; m < 5; ++m) {
                ee.col(m) = cell(i - 2 + m, j).vyl;
                ff.col(m) = cell(i - 2 + m, j).vyr;
            }
            for (eq = 0; eq < 4; ++eq) {
                cell(i, j).vygl.row(eq) = lagrangeObj.weno5JS3PointsQuadrature(ee.row(eq));
                cell(i, j).vygr.row(eq) = lagrangeObj.weno5JS3PointsQuadrature(ff.row(eq));
            }
        }
    }

#pragma omp parallel for private(i, j)
    for (i = k - 1; i < nxx - k; ++i) {
        for (j = k - 1; j < nyy - k; ++j) {
            cell(i, j).transformToConservativeAtGQ();
        }
    }

#pragma omp parallel for private(i, j, m, temporaryFluxUx, temporaryFluxUy, temporaryFluxVx, \
temporaryFluxVy, temporaryFluxWx, temporaryFluxWy,pvarL, pvarR, dFluxL, dFluxR, aM) firstprivate(eulObj, threePointsGQWeights)
    for (i = k - 1; i < nxx - k; ++i) {
        for (j = k - 1; j < nyy - k; ++j) {
            cell(i, j).setZeroFlux();

            for (m = 0; m < 3; ++m) {
                pvarL = eulObj.primitive(cell(i, j).uxgl.col(m));
                dFluxL = eulObj.dfluxx(cell(i, j).uxgl.col(m), pvarL);
                pvarR = eulObj.primitive(cell(i, j).uxgr.col(m));
                dFluxR = eulObj.dfluxx(cell(i, j).uxgr.col(m), pvarR);

                /*try
                {
                    if (pvarL(0) <= 0.0 || pvarL(3) < 0.0)
                        throw (-10);
                    if (pvarR(0) <= 0.0 || pvarR(3) < 0.0)
                        throw (-11);
                }
                catch (const int& x)
                {
                    cout << "Exception occurred, exception number: " << x << endl;
                    exit(-1);
                }*/
//                aM = eulObj.maxLocalWaveSpeed(pvarL.cwiseMin(pvarR), pvarR.cwiseMax(pvarL), 1);
                aM = eulObj.maxLocalWaveSpeed(pvarL, pvarR, 1) / 1.0;

                temporaryFluxUx.col(m) = threePointsGQWeights(m) *
                                        eulObj.Rusanov(cell(i, j).uxgl.col(m), cell(i, j).uxgr.col(m), pvarL, pvarR,
                                                       cell(i, j).leftEigenVectorX, cell(i, j).rightEigenVectorX,
                                                       aM, 1);

//                aM = eulObj.maxLocalWaveSpeed(pvarL.cwiseMin(pvarR), pvarR.cwiseMax(pvarL), 1);

                temporaryFluxVx.col(m) = threePointsGQWeights(m) *
                                 eulObj.dRusanov(cell(i, j).vxgl.col(m), cell(i, j).vxgr.col(m),
                                                 pvarL, pvarR, cell(i, j).leftEigenVectorX,
                                                 cell(i, j).rightEigenVectorX, dFluxL, dFluxR, aM);

                temporaryFluxWx.col(m) = threePointsGQWeights(m) *
                                 eulObj.dRusanov(cell(i, j).wxgl.col(m), cell(i, j).wxgr.col(m),
                                                 pvarL, pvarR, cell(i, j).leftEigenVectorX,
                                                 cell(i, j).rightEigenVectorX, dFluxL, dFluxR, aM);

                pvarL = eulObj.primitive(cell(i, j).uygl.col(m));
                dFluxL = eulObj.dfluxy(cell(i, j).uygl.col(m), pvarL);
                pvarR = eulObj.primitive(cell(i, j).uygr.col(m));
                dFluxR = eulObj.dfluxy(cell(i, j).uygr.col(m), pvarR);

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
//                aM = eulObj.maxLocalWaveSpeed(pvarL.cwiseMin(pvarR), pvarR.cwiseMax(pvarL), 2);
                aM = eulObj.maxLocalWaveSpeed(pvarL, pvarR, 2) / 1.0;

                temporaryFluxUy.col(m) = threePointsGQWeights(m) *
                                 eulObj.Rusanov(cell(i, j).uygl.col(m), cell(i, j).uygr.col(m), pvarL, pvarR,
                                                cell(i, j).leftEigenVectorY, cell(i, j).rightEigenVectorY, aM, 2);

//                aM = eulObj.maxLocalWaveSpeed(pvarL.cwiseMin(pvarR), pvarR.cwiseMax(pvarL), 2);

                temporaryFluxVy.col(m) = threePointsGQWeights(m) *
                                 eulObj.dRusanov(cell(i, j).vygl.col(m), cell(i, j).vygr.col(m),
                                                 pvarL, pvarR, cell(i, j).leftEigenVectorY,
                                                 cell(i, j).rightEigenVectorY, dFluxL, dFluxR, aM);

                temporaryFluxWy.col(m) = threePointsGQWeights(m) *
                                 eulObj.dRusanov(cell(i, j).wygl.col(m), cell(i, j).wygr.col(m),
                                                 pvarL, pvarR, cell(i, j).leftEigenVectorY,
                                                 cell(i, j).rightEigenVectorY, dFluxL, dFluxR, aM);

            }
            cell(i, j).fx = temporaryFluxUx.rowwise().sum();
            cell(i, j).gx = temporaryFluxVx.rowwise().sum();
            cell(i, j).hx = temporaryFluxWx.rowwise().sum();

            cell(i, j).fy = temporaryFluxUy.rowwise().sum();
            cell(i, j).gy = temporaryFluxVy.rowwise().sum();
            cell(i, j).hy = temporaryFluxWy.rowwise().sum();
        }
    }

//#pragma omp parallel for private(i, j)
//    for (i = k - 1; i < nxx - k; ++i) {
//        for (j = k - 1; j < nyy - k; ++j) {
//            cell(i, j).transformToConservative();
//        }
//    }
//
//#pragma omp parallel for private(i, j, m) firstprivate(pvarL, pvarR, dFluxL, dFluxR, aM, eulObj)
//    for (i = k - 1; i < nxx - k; ++i) {
//        for (j = k - 1; j < nyy - k; ++j) {
//            cell(i, j).setZeroFlux();
//
//            pvarL = eulObj.primitive(cell(i, j).uxl);
//            dFluxL = eulObj.dfluxx(cell(i, j).uxl, pvarL);
//            pvarR = eulObj.primitive(cell(i, j).uxr);
//            dFluxR = eulObj.dfluxx(cell(i, j).uxr, pvarR);
//
//            aM = eulObj.maxLocalWaveSpeed(pvarL, pvarR, 1);
////            dfdu = eulObj.maxLocalWaveSpeed(pvarL, pvarR, cell(i,j).leftEigenVectorX, cell(i,j).rightEigenVectorX, 1);
//
//            cell(i, j).fx = eulObj.Rusanov(cell(i, j).uxl, cell(i, j).uxr, pvarL, pvarR,
//                                           cell(i, j).leftEigenVectorX, cell(i, j).rightEigenVectorX, aM, 1);
//
//            cell(i, j).gx = eulObj.dRusanov(cell(i, j).vxl, cell(i, j).vxr,
//                                            pvarL, pvarR, cell(i, j).leftEigenVectorX,
//                                            cell(i, j).rightEigenVectorX, dFluxL, dFluxR, aM);
//
//            cell(i, j).hx = eulObj.dRusanov(cell(i, j).wxl, cell(i, j).wxr,
//                                            pvarL, pvarR, cell(i, j).leftEigenVectorX,
//                                            cell(i, j).rightEigenVectorX, dFluxL, dFluxR, aM);
////            cell(i, j).fx = eulObj.Rusanov(cell(i, j).uxl, cell(i, j).uxr, pvarL, pvarR, dfdu, 1);
////
////            cell(i, j).gx = eulObj.dRusanov(cell(i, j).vxl, cell(i, j).vxr,
////                                            pvarL, pvarR, dFluxL, dFluxR, dfdu);
////
////            cell(i, j).hx = eulObj.dRusanov(cell(i, j).wxl, cell(i, j).wxr,
////                                            pvarL, pvarR, dFluxL, dFluxR, dfdu);
//
//
//            pvarL = eulObj.primitive(cell(i, j).uyl);
//            dFluxL = eulObj.dfluxy(cell(i, j).uyl, pvarL);
//            pvarR = eulObj.primitive(cell(i, j).uyr);
//            dFluxR = eulObj.dfluxy(cell(i, j).uyr, pvarR);
//
//            aM = eulObj.maxLocalWaveSpeed(pvarL, pvarR, 2);
////            dfdu = eulObj.maxLocalWaveSpeed(pvarL, pvarR, cell(i,j).leftEigenVectorY, cell(i,j).rightEigenVectorY, 2);
//
//            cell(i, j).fy = eulObj.Rusanov(cell(i, j).uyl, cell(i, j).uyr, pvarL, pvarR,
//                                           cell(i, j).leftEigenVectorY, cell(i, j).rightEigenVectorY, aM, 2);
//
//            cell(i, j).gy = eulObj.dRusanov(cell(i, j).vyl, cell(i, j).vyr,
//                                            pvarL, pvarR, cell(i, j).leftEigenVectorY,
//                                            cell(i, j).rightEigenVectorY, dFluxL, dFluxR, aM);
//
//            cell(i, j).hy = eulObj.dRusanov(cell(i, j).wyl, cell(i, j).wyr,
//                                            pvarL, pvarR, cell(i, j).leftEigenVectorY,
//                                            cell(i, j).rightEigenVectorY, dFluxL, dFluxR, aM);
////            cell(i, j).fy = eulObj.Rusanov(cell(i, j).uyl, cell(i, j).uyr, pvarL, pvarR, dfdu, 2);
////
////            cell(i, j).gy = eulObj.dRusanov(cell(i, j).vyl, cell(i, j).vyr,
////                                            pvarL, pvarR, dFluxL, dFluxR, dfdu);
////
////            cell(i, j).hy = eulObj.dRusanov(cell(i, j).wyl, cell(i, j).wyr,
////                                            pvarL, pvarR, dFluxL, dFluxR, dfdu);
//
//        }
//    }

#pragma omp parallel for private(i, j, source)
    for (i = k; i < nxx - k; i++) {
        for (j = k; j < nyy - k; j++) {
            source.setZero();
//            source(2) = cell(i, j).u(0);
//            source(3) = cell(i, j).u(2);

            cell(i, j).rhsu = ((cell(i - 1, j).fx - cell(i, j).fx) * dy +
                               (cell(i, j - 1).fy - cell(i, j).fy) * dx) / area + source;
//
//            source(2) = cell(i, j).v(0);
//            source(3) = cell(i, j).v(2);
            cell(i, j).rhsv = ((cell(i - 1, j).gx - cell(i, j).gx) * dy +
                               (cell(i, j - 1).gy - cell(i, j).gy) * dx) / area + source;
//
//            source(2) = cell(i, j).w(0);
//            source(3) = cell(i, j).w(2);
            cell(i, j).rhsw = ((cell(i - 1, j).hx - cell(i, j).hx) * dy +
                               (cell(i, j - 1).hy - cell(i, j).hy) * dx) / area + source;
        }
    }
}

void FVHermite::System2D::RungeKutta3() {
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
            cell(i, j).v1 = ssprkObj.RungeKuttaIII<Vector4d>(cell(i, j).v0,
                                                             cell(i, j).v0, cell(i, j).rhsv, dt, stage);
            cell(i, j).w1 = ssprkObj.RungeKuttaIII<Vector4d>(cell(i, j).w0,
                                                             cell(i, j).w0, cell(i, j).rhsw, dt, stage);
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
            cell(i, j).v2 = ssprkObj.RungeKuttaIII<Vector4d>(cell(i, j).v0,
                                                             cell(i, j).v1, cell(i, j).rhsv, dt, stage);
            cell(i, j).w2 = ssprkObj.RungeKuttaIII<Vector4d>(cell(i, j).w0,
                                                             cell(i, j).w1, cell(i, j).rhsw, dt, stage);
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
            cell(i, j).v = ssprkObj.RungeKuttaIII<Vector4d>(cell(i, j).v0,
                                                            cell(i, j).v2, cell(i, j).rhsv, dt, stage);
            cell(i, j).w = ssprkObj.RungeKuttaIII<Vector4d>(cell(i, j).w0,
                                                            cell(i, j).w2, cell(i, j).rhsw, dt, stage);
        }
    }
    calcPrimitiveVariable();
}

void FVHermite::System2D::timeStepping() {
    int i, j;
    double maxSpeedx, maxSpeedy;
    calcPrimitiveVariable();

//    bool a1 = true, a2 = true;

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

        dt = cfl * pow(fmin(dx, dy), 3.0 / 3.0) / fmax(maxSpeedx, maxSpeedy) / 1.0;
        dt = t + dt > tFinal ? tFinal - t : dt;

//        if (t + dt > 100.0 && a1) {
//            dt = 100.0 - t;
//            write(100);
//            a1 = false;
//        } else if (t + dt > tFinal && a2) {
//            dt = tFinal - t;
//            write(1000);
//            a2 = false;
//        }

        RungeKutta3();
        t += dt;
        cout << "tnum: " << t << endl;
    }
}

void FVHermite::System2D::printResult() {
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
    write(99988);
}

void FVHermite::System2D::write(const int &filenum) {

    std::string filename = "../hasil/2D_RIEMANN_HTENO(Test)_" + std::to_string(nx) + "_" + std::to_string(ny) +
                           "_"  + "Final" + ".dat";

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

    std::cout << "\n Solution is written in " << filename << "\n";
    outfile.close();
}

void FVHermite::System2D::isentropicVortex() {
    int i, j;

    cout << "Initializing variable... case: Entropy Waves...";
#pragma omp parallel for private(i, j) firstprivate(lagrangeObj, eulObj)
    for (i = 0; i < nxx; i++) {
        for (j = 0; j < nyy; j++) {
            cell(i, j).u = eulObj.isentropicVortex(cell(i, j).x, cell(i, j).y,
                                                   dx, dy, lagrangeObj.getQuadrature9Weights(), lagrangeObj.getQuadrature9Points());
            cell(i, j).v = eulObj.dxIsentropicVortex(cell(i, j).x, cell(i, j).y,
                                                   dx, dy, lagrangeObj.getQuadrature9Weights(), lagrangeObj.getQuadrature9Points());
            cell(i, j).w = eulObj.dyIsentropicVortex(cell(i, j).x, cell(i, j).y,
                                                     dx, dy, lagrangeObj.getQuadrature9Weights(), lagrangeObj.getQuadrature9Points());
        }
    }
    periodicBc();

    cout << " done\n";
#pragma omp parallel for private(i, j) firstprivate(lagrangeObj, eulObj)
    for (i = k; i < nxx - k; i++) {
        for (j = k; j < nyy - k; j++) {
            cell(i, j).exact = eulObj.isentropicVortex(cell(i, j).x - tFinal, cell(i, j).y - tFinal,
                                                       dx, dy, lagrangeObj.getQuadrature9Weights(), lagrangeObj.getQuadrature9Points());
            cell(i,j).exact = eulObj.primitive(cell(i,j).exact);
        }
    }
    cout << "Exact solutions has been computed...\n";
}
