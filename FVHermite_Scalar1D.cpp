#include "FVHermite_Scalar1D.h"

FVHermite::Scalar1D::Scalar1D(const double &xl, const double &xr, const int &nx, \
    const double &cflnumber, const double &tmax, const double &tInit) {

    int i;

    this->cfl = cflnumber;
    this->tFinal = tmax;
    this->t = tInit;
    dt = 0.0;
    this->nx = nx;

    dx = (xr - xl) / nx;
    area = dx;

    k = 5;
    s = 1;
    nxx = k + k + nx;
    n = nx;

    hermiteObj.setGridSize(dx);

    cout << "Creating domain... ";
    cell.resize(nxx);

    for (i = 0; i < nxx; i++) {
        cell(i).x = (double(i) - double(k) + 0.5) * dx + xl;
    }

    cout << "done.\n";
    cout << "Number of cells: " << n << "\n";
}

void FVHermite::Scalar1D::initialize() {
    int i;
    cout << "Initializing variable...";
    for (i = k; i < nxx - k; i++) {
        cell(i).u = linAdv1DObj.sine(cell(i).x, dx, \
            lagrangeObj.getQuadrature9Weights(), lagrangeObj.getQuadrature9Points());
        cell(i).v = linAdv1DObj.sinedx(cell(i).x, dx, \
            lagrangeObj.getQuadrature9Weights(), lagrangeObj.getQuadrature9Points());
    }
    cout << "done.\n";
}

void FVHermite::Scalar1D::periodicBc() {
    for (int m = 0; m < k; m++) {
        cell(k - 1 - m).u = cell(nx + k - 1 - m).u;
        cell(nxx - 1 - m).u = cell(k + k - 1 - m).u;

        cell(k - 1 - m).v = cell(nx + k - 1 - m).v;
        cell(nxx - 1 - m).v = cell(k + k - 1 - m).v;
    }
}

void FVHermite::Scalar1D::limitDerivative() {
    int i, j;
    Vector3d uu, vv;

    for (i = k; i < nxx - k; i++) {
        for (j = 0; j < 3; ++j) {
            uu(j) = cell(i - 1 + j).u;
            vv(j) = cell(i - 1 + j).v;
        }
        cell(i).v = hermiteObj.hermiteLimiter(uu, vv);
    }
}

void FVHermite::Scalar1D::rhs() {
    int i, m;
//    periodicBc();
//    limitDerivative();
    periodicBc();

    for (i = s; i < nxx - s - 1; ++i) {
        for (m = 0; m < 4; ++m) {
            interpolateA(m) = cell(i - 1 + m).u;
            interpolateB(m) = cell(i - 1 + m).v;
        }
        valueA = hermiteObj.hweno(interpolateA, interpolateB);
        cell(i).uxl = valueA(0);
        cell(i).uxr = valueA(1);
        cell(i).vxl = valueA(2);
        cell(i).vxr = valueA(3);
    }

    for (i = s + 1; i < nxx - s - 2; ++i) {
        cell(i).fx = burgersObj.Rusanov(cell(i).uxl, cell(i).uxr);
        cell(i).gx = burgersObj.dRusanov(cell(i).uxl, cell(i).uxr, cell(i).vxl, cell(i).vxr);
    }

    for (i = k; i < nxx - k; i++) {
        cell(i).rhsu = (cell(i - 1).fx - cell(i).fx) / area;
        cell(i).rhsv = (cell(i - 1).gx - cell(i).gx) / area;
    }
}

void FVHermite::Scalar1D::RungeKutta3() {
    int i;
    int stage;

    stage = 1;
    rhs();
    for (i = k; i < nxx - k; i++) {
        cell(i).save();
        cell(i).u1 = ssprkObj.RungeKuttaIII<double>(cell(i).u0, \
            cell(i).u0, cell(i).rhsu, dt, stage);
        cell(i).v1 = ssprkObj.RungeKuttaIII<double>(cell(i).v0, \
            cell(i).v0, cell(i).rhsv, dt, stage);

        cell(i).update1();
    }

    stage += 1;
    rhs();
    for (i = k; i < nxx - k; i++) {
        cell(i).u2 = ssprkObj.RungeKuttaIII<double>(cell(i).u0, \
            cell(i).u1, cell(i).rhsu, dt, stage);
        cell(i).v2 = ssprkObj.RungeKuttaIII<double>(cell(i).v0, \
            cell(i).v1, cell(i).rhsv, dt, stage);
        cell(i).update2();
    }

    stage += 1;
    rhs();
    for (i = k; i < nxx - k; i++) {
        cell(i).u = ssprkObj.RungeKuttaIII<double>(cell(i).u0, \
            cell(i).u2, cell(i).rhsu, dt, stage);
        cell(i).v = ssprkObj.RungeKuttaIII<double>(cell(i).v0, \
            cell(i).v2, cell(i).rhsv, dt, stage);
    }
}

void FVHermite::Scalar1D::timeStepping() {
    int i;
    double maxSpeedx;

    while (t < tFinal) {
        maxSpeedx = 0.0;
        for (i = k; i < nxx - k; i++) {
            maxSpeedx = fmax(maxSpeedx, fabs(linAdvObj.df(cell(i).u)));
        }

        dt = cfl * pow(dx, 3.0 / 3.0) / maxSpeedx;
        dt = t + dt > tFinal ? tFinal - t : dt;

        RungeKutta3();
        t += dt;
        cout << "tnum: " << t << endl;
    }
}

void FVHermite::Scalar1D::Exact() {
    for (int i = k; i < nxx - k; i++) {
        cell(i).exact = linAdv1DObj.sine(cell(i).x - tFinal, dx, lagrangeObj.getQuadrature9Weights(), \
            lagrangeObj.getQuadrature9Points());
    }
}

void FVHermite::Scalar1D::printResult() {
    int i;
    double err;

    err = 0.0;
    for (i = k; i < nxx - k; i++) {
        err += fabs(cell(i).u - cell(i).exact);
    }
    cout << "Numerical error:" << err / n << endl;

    std::string filename = "../hasil/1D_BURGERS_HTENO(R3)_" + std::to_string(nx) + "_" + ".dat";
    std::ofstream outfile;

    const char separator = ' ';
    const int strwidth = 20;
    const int numwidth = 20;
    const int pn = 8;


    outfile.open(filename);
    if (outfile.good()) {
        outfile << "t: " << t << "\n";
        outfile << std::right << std::setw(strwidth) << std::setfill(separator) << "x";
        outfile << std::right << std::setw(strwidth) << std::setfill(separator) << "u";
        outfile << "\n";
    } else {
        exit(1);
    }

    for (i = k; i < nxx - k; ++i) {
        outfile << std::right << std::setw(numwidth) << std::setfill(separator) << std::setprecision(pn) << cell(i).x;
        outfile << std::right << std::setw(numwidth) << std::setfill(separator) << std::setprecision(pn) << cell(i).u;
        outfile << "\n";
    }

    std::cout << "\n Solution is written in " << filename << "\n";
    outfile.close();
}
