#include "FVHermite_System1D.h"

FVHermite::System1D::System1D(const double &xl, const double &xr, const int &nx,
                              const double &cflnumber, const double &tmax, const double &tInit,
                              const double &specificHeatRatio) :
        gamma(specificHeatRatio), gamma_1(specificHeatRatio - 1.0), eulObj(specificHeatRatio) {

    int i;

    this->cfl = cflnumber;
    this->tFinal = tmax;
    this->t = tInit;
    dt = 0.0;
    this->nx = nx;
    dx = (xr - xl) / nx;

    circleRadMod = pow(dx / 2.0, 1.5);

    area = dx;
    k = 4;
    s = 1;
    nxx = k + k + nx;
    n = nx;

    hermiteObj.setGridSize(dx);

    cout << "Creating domain... ";
    cell.resize(nxx);

    for (i = 0; i < nxx; i++) {
        cell(i).x = ((double) i - (double) k + 0.5) * dx + xl;
    }

    cout << "done.\n";
    cout << "Number of cells: " << n << "\n";

    t0 = 0.0, t1 = 0.0, vxm = 0.0, hm = 0.0, qm = 0.0, cm = 0.0;
    rcm = 0.0, b1 = 0.0, b2 = 0.0, t2 = 0.0;
}

void FVHermite::System1D::entropyWave() {
    int i;
    double rho = 1.0, u = 1.0, p = 1.0;
    double A = 0.5;

    cout << "Initializing variable... case: Entropy Waves...";

    for (i = 0; i < nxx; i++) {
        cell(i).u = eulObj.entropyWaves(rho, u, p, A, cell(i).x,
                                        dx, lagrangeObj.getQuadrature9Weights(), lagrangeObj.getQuadrature9Points());
        cell(i).v = eulObj.entropyWavesdx(rho, u, p, A, cell(i).x,
                                          dx, lagrangeObj.getQuadrature9Weights(), lagrangeObj.getQuadrature9Points());
    }

    cout << "done\n";

    for (i = k; i < nxx - k; i++) {
        cell(i).exact = eulObj.entropyWaves(rho, u, p, A,
                                            cell(i).x - u * tFinal, dx, lagrangeObj.getQuadrature9Weights(),
                                            lagrangeObj.getQuadrature9Points());

        cell(i).exact = eulObj.primitive(cell(i).exact);
    }
    cout << "Exact solutions has been computed...\n";
}

void FVHermite::System1D::blastWave() {
    int i;

    for (i = 0; i < nxx; i++) {
        cell(i).u = eulObj.conservative(eulObj.interactingBlast(cell(i).x));
        cell(i).v = eulObj.conservative(eulObj.interactingBlast(cell(i).x + dx / 2.0) -
                                        eulObj.interactingBlast(cell(i).x - dx / 2.0)) / dx;
    }
}

void FVHermite::System1D::periodicBc() {

    for (int m = 0; m < k; m++) {
        cell(k - 1 - m).u = cell(nx + k - 1 - m).u;
        cell(nxx - 1 - m).u = cell(k + k - 1 - m).u;

        cell(k - 1 - m).v = cell(nx + k - 1 - m).v;
        cell(nxx - 1 - m).v = cell(k + k - 1 - m).v;
    }
}

void FVHermite::System1D::zeroGradBc() {
    for (int m = 0; m < k; m++) {
        cell(k - 1 - m).u = cell(k).u;
        cell(nxx - 1 - m).u = cell(nxx - k - 1).u;

        cell(k - 1 - m).v = cell(k).v;
        cell(nxx - 1 - m).v = cell(nxx - k - 1).v;
    }
}

void FVHermite::System1D::reflectiveBc() {

    for (int m = 0; m < k; m++) {
        cell(k - 1 - m).u(0) = cell(k + m).u(0);
        cell(k - 1 - m).u(1) = -cell(k + m).u(1);
        cell(k - 1 - m).u(2) = cell(k + m).u(2);

        cell(nxx - 1 - m).u(0) = cell(nx - 1 + m).u(0);
        cell(nxx - 1 - m).u(1) = -cell(nx - 1 + m).u(1);
        cell(nxx - 1 - m).u(2) = cell(nx - 1 + m).u(2);

        cell(k - 1 - m).v(0) = -cell(k + m).v(0);
        cell(k - 1 - m).v(1) = cell(k + m).v(1);
        cell(k - 1 - m).v(2) = -cell(k + m).v(2);

        cell(nxx - 1 - m).v(0) = -cell(nx - 1 + m).v(0);
        cell(nxx - 1 - m).v(1) = cell(nx - 1 + m).v(1);
        cell(nxx - 1 - m).v(2) = -cell(nx - 1 + m).v(2);
    }
}

void FVHermite::System1D::calcPrimitiveVariable() {

    for (int i = 0; i < nxx; i++) {
        cell(i).pvar = eulObj.primitive(cell(i).u);

        /*try {
            if (cell(i).pvar(0) <= 0.0)
                throw (-1);
        }
        catch (int x) {
            cout << "Exception occurred, exception number: " << x << endl;
        }

        try {
            if (cell(i).pvar(2) <= 0.0)
                throw (-2);
        }
        catch (int x) {
            cout << "Exception occurred, exception number: " << x << endl;
        }*/

        cell(i).sqrtRho = sqrt(fabs(cell(i).pvar(0)));
        cell(i).enthalpy = (cell(i).u(2) + cell(i).pvar(2)) / cell(i).pvar(0);
    }
}

void FVHermite::System1D::eigen1DEuler() {
    int i;
//    aM = -aM.setOnes();

    for (i = 0; i < nxx - 1; i++) {
        t0 = cell(i).sqrtRho / (cell(i).sqrtRho + cell(i + 1).sqrtRho);
        t1 = 1.0 - t0;
        vxm = t0 * cell(i).pvar(1) + t1 * cell(i + 1).pvar(1);
        hm = t0 * cell(i).enthalpy + t1 * cell(i + 1).enthalpy;
        qm = 0.5 * vxm * vxm;

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
        cell(i).rightEigenVectorX = eulObj.rightEigenVector(vxm, hm, qm, cm, t0);

        rcm = 1.0 / cm;
        b1 = gamma_1 * rcm * rcm;
        b2 = qm * b1;
        t0 = vxm * rcm;
        t1 = b1 * vxm;
        t2 = 0.5 * b1;

        cell(i).leftEigenVectorX = eulObj.leftEigenVector(rcm, b1, b2, t0, t1, t2);

//        aM = aM.cwiseMax(eulObj.localWaveSpeed(cell(i).pvar));

    }
}

void FVHermite::System1D::limitDerivatives() {
    int i, j, eq;
    Matrix3d uu, vv;

    for (i = k - 1; i < nxx - k + 1; i++) {

        for (j = 0; j < 3; ++j) {
            uu.col(j) = cell(i - 1 + j).u;
            vv.col(j) = cell(i - 1 + j).v;
        }
        for (eq = 0; eq < 3; eq++) {
            cell(i).v(eq) = hermiteObj.hermiteLimiter(uu.row(eq), vv.row(eq));
        }
    }
}

void FVHermite::System1D::rhs() {
    int i, j, m, eq;
    Matrix3d uu, vv;
    periodicBc();
//    zeroGradBc();
//    reflectiveBc();
//    limitDerivatives();
//    periodicBc();
//    zeroGradBc();
//    reflectiveBc();

    calcPrimitiveVariable();
    eigen1DEuler();

    // Reset
//    for (i = s + 1; i < nxx - s - 1; ++i) {
//        cell(i).discontinue.setZero();
//    }

    // KXRCF
//    for (i = s + 1; i < nxx - s - 1; ++i) {
//        if (!cell(i).discontinue) {
//            if (cell(i).u(1) < 0.0) {
//                h1 = (2.0 * cell(i-1).u(0) + 5.0 * cell(i).u(0) - cell(i+1).u(0)) / 6.0;
//                h2 = (-cell(i-2).u(0) + 5.0 * cell(i-1).u(0) + 2.0 * cell(i).u(0)) / 6.0;
//            } else {
//                h1 = (-cell(i-1).u(0) + 5.0 * cell(i).u(0) + 2.0 * cell(i+1).u(0)) / 6.0;
//                h2 = (2.0 * cell(i).u(0) + 5.0 * cell(i+1).u(0) - cell(i+2).u(0)) / 6.0;
//            }
//
//            cell(i).discontinue = fabs(h1 - h2) / (circleRadMod * fabs(cell(i).u(0))) > 1.0;
//
//            if (cell(i).discontinue){
//                cell(i-1).discontinue = true;
//                cell(i+1).discontinue = true;
//            }
//        } else
//            continue;
//    }

//    LMPE
    double A, xc;

    for (eq = 0; eq < 3; eq++) {
        for (i = s + 1; i < nxx - s - 1; ++i) {
            if (!cell(i).discontinue(eq)) {
                A = fabs(-cell(i - 1).u(eq) + 2.0 * cell(i).u(eq) - cell(i + 1).u(eq)) / 2.0;
                if (A > 0.5 * dx) {
                    xc = (-2. * cell(i - 1).u(eq) + 3. * cell(i).u(eq) - cell(i + 1).u(eq)) /
                         (-cell(i - 1).u(eq) + 2.0 * cell(i).u(eq) - cell(i + 1).u(eq)) * dx;
                    if (0.0 <= xc && xc <= 3. * dx) {
                        cell(i).discontinue(eq) = true;
                        cell(i - 1).discontinue(eq) = true;
                        cell(i + 1).discontinue(eq) = true;
                    }
                }
            }
        }
    }

    for (i = s + 1; i < nxx - s - 1; ++i) {
        if (cell(i).discontinue(0) || cell(i).discontinue(1) || cell(i).discontinue(2)) {
            for (i = k - 1; i < nxx - k + 1; i++) {

                for (j = 0; j < 3; ++j) {
                    uu.col(j) = cell(i - 1 + j).u;
                    vv.col(j) = cell(i - 1 + j).v;
                }
                for (eq = 0; eq < 3; eq++) {

                    cell(i).v(eq) = hermiteObj.hermiteLimiter(uu.row(eq), vv.row(eq));
                }
            }
        }
    }

    periodicBc();
//    zeroGradBc();
//    reflectiveBc();

    for (i = s; i < nxx - s - 1; ++i) {
        for (m = 0; m < 4; ++m) {
            interpolateA.col(m) = cell(i).leftEigenVectorX * cell(i - 1 + m).u;
            interpolateB.col(m) = cell(i).leftEigenVectorX * cell(i - 1 + m).v;
        }

        for (eq = 0; eq < 3; eq++) {
            if (cell(i).discontinue(eq)) {
                valueA = hermiteObj.htenoi(interpolateA.row(eq), interpolateB.row(eq));
            } else {
                valueA = hermiteObj.linear(interpolateA.row(eq), interpolateB.row(eq));
            }

//            valueA = hermiteObj.htenoi(interpolateA.row(eq), interpolateB.row(eq));

            cell(i).uxl(eq) = valueA(0);
            cell(i).uxr(eq) = valueA(1);
            cell(i).vxl(eq) = valueA(2);
            cell(i).vxr(eq) = valueA(3);

        }
    }

    for (i = k - 1; i < nxx - k; ++i) {
        cell(i).transformToConservative();
    }

    for (i = k - 1; i < nxx - k; ++i) {
        cell(i).setZeroFlux();
        pvarL = eulObj.primitive(cell(i).uxl);
        dFluxL = eulObj.dflux(cell(i).uxl, pvarL);
        pvarR = eulObj.primitive(cell(i).uxr);
        dFluxR = eulObj.dflux(cell(i).uxr, pvarR);

        /*try {
            if (pvarL(0) <= 0.0 || pvarL(2) < 0.0)
                throw (-10);
            if (pvarR(0) <= 0.0 || pvarR(2) < 0.0)
                throw (-11);
        }
        catch (const int &x) {
            cout << "exception occurred, exception number: " << x << endl;
            exit(-1);
        }*/

//        aM = eulObj.maxLocalWaveSpeed(pvarL.cwiseMax(pvarR), pvarR.cwiseMin(pvarL));
        aM = eulObj.maxLocalWaveSpeed(pvarL, pvarR);

        cell(i).fx = eulObj.Rusanov(cell(i).uxl, cell(i).uxr, pvarL, pvarR,
                                    cell(i).leftEigenVectorX, cell(i).rightEigenVectorX, aM);
//		cell(i).fx = eulObj.HLLC(cell(i).uxl, cell(i).uxr);
//		cell(i).fx = eulObj.CUFlux(cell(i).uxl, cell(i).uxr, cell(i).leftEigenVectorX, cell(i).rightEigenVectorX);

        cell(i).gx = eulObj.dRusanov(cell(i).vxl, cell(i).vxr,
                                     pvarL, pvarR, cell(i).leftEigenVectorX,
                                     cell(i).rightEigenVectorX, dFluxL, dFluxR, aM);
//		cell(i).gx = eulObj.dCUFlux(cell(i).uxl, cell(i).uxr, cell(i).leftEigenVectorX, cell(i).rightEigenVectorX,
//				cell(i).vxl, cell(i).vxr, dFluxL, dFluxR);
    }

    for (i = k; i < nxx - k; i++) {
        cell(i).rhsu = (cell(i - 1).fx - cell(i).fx) / area;
        cell(i).rhsv = (cell(i - 1).gx - cell(i).gx) / area;
    }
}

void FVHermite::System1D::RungeKutta3() {
    int i;
    int stage;

    stage = 1;
    rhs();
    for (i = k; i < nxx - k; i++) {
        cell(i).save();
        cell(i).u1 = ssprkObj.RungeKuttaIII<Vector3d>(cell(i).u0,
                                                      cell(i).u0, cell(i).rhsu, dt, stage);
        cell(i).v1 = ssprkObj.RungeKuttaIII<Vector3d>(cell(i).v0,
                                                      cell(i).v0, cell(i).rhsv, dt, stage);

        cell(i).update1();
    }

    stage += 1;
    rhs();
    for (i = k; i < nxx - k; i++) {
        cell(i).u2 = ssprkObj.RungeKuttaIII<Vector3d>(cell(i).u0,
                                                      cell(i).u1, cell(i).rhsu, dt, stage);
        cell(i).v2 = ssprkObj.RungeKuttaIII<Vector3d>(cell(i).v0,
                                                      cell(i).v1, cell(i).rhsv, dt, stage);
        cell(i).update2();
    }

    stage += 1;
    rhs();
    for (i = k; i < nxx - k; i++) {
        cell(i).u = ssprkObj.RungeKuttaIII<Vector3d>(cell(i).u0,
                                                     cell(i).u2, cell(i).rhsu, dt, stage);
        cell(i).v = ssprkObj.RungeKuttaIII<Vector3d>(cell(i).v0,
                                                     cell(i).v2, cell(i).rhsv, dt, stage);
    }

    calcPrimitiveVariable();
}

void FVHermite::System1D::timeStepping() {
    int i;
    double maxSpeedx;
    calcPrimitiveVariable();

    const clock_t begin_time = clock();

    while (t < tFinal) {
        maxSpeedx = 0.0;
        for (i = k; i < nxx - k; i++) {
            maxSpeedx = fmax(maxSpeedx, eulObj.maxSpeed(cell(i).pvar));
        }

        dt = cfl * pow(dx, 5.0 / 3.0) / maxSpeedx;
        dt = t + dt > tFinal ? tFinal - t : dt;

        RungeKutta3();
        t += dt;
//        cout << "tnum: " << t << endl;
    }

    comp_time = float( clock () - begin_time ) /  CLOCKS_PER_SEC;
}

void FVHermite::System1D::printResult() {
    int i;
    double err;

    err = 0.0;
    for (i = k; i < nxx - k; i++) {
        err += fabs(cell(i).u(0) - cell(i).exact(0));
    }

    cout << "N: " << n << endl;
    cout << "L1 error: " << err / n << endl;
    cout << "CPU Time: " << comp_time << std::endl;
//    write(1);
}

void FVHermite::System1D::write(const int &filenum) {

    std::string filename = "/home/indro/Documents/1D_SHU_KXRCF_HTENO_" + std::to_string(nx) +
                           "_" + std::to_string(filenum) + ".dat";

    std::ofstream outfile;
    int i;
    const char separator = ' ';
    const int strwidth = 20;
    const int numwidth = 20;
    const int pn = 8;

    outfile.open(filename);
    if (outfile.good()) {
        outfile << "final time: " << t << "\n";
        outfile << "computational time: " << comp_time << "\n\n";
        outfile << std::right << std::setw(strwidth) << std::setfill(separator) << "x";
        outfile << std::right << std::setw(strwidth) << std::setfill(separator) << "Density";
        outfile << std::right << std::setw(strwidth) << std::setfill(separator) << "Ux";
        outfile << std::right << std::setw(strwidth) << std::setfill(separator) << "Pressure";
        outfile << "\n";
    } else {
        exit(1);
    }

    for (i = 0; i < nxx; ++i) {
        outfile << std::right << std::setw(numwidth) << std::setfill(separator) << std::setprecision(pn) << cell(i).x;
        outfile << std::right << std::setw(numwidth) << std::setfill(separator) << std::setprecision(pn)
                << cell(i).pvar(0);
        outfile << std::right << std::setw(numwidth) << std::setfill(separator) << std::setprecision(pn)
                << cell(i).pvar(1);
        outfile << std::right << std::setw(numwidth) << std::setfill(separator) << std::setprecision(pn)
                << cell(i).pvar(2);
        outfile << std::right << std::setw(numwidth) << std::setfill(separator) << std::setprecision(pn)
                << cell(i).discontinue.transpose();
        outfile << "\n";
    }

    std::cout << "\n Solution is written in " << filename << "\n";
    outfile.close();
}
