#include "FV_System1D.h"

FV::System1D::System1D(const double& xl, const double& xr, const int& nx,
		const double& cflnumber, const double& tmax, const double& tInit,
		const double& specificHeatRatio) :
		gamma(specificHeatRatio), gamma_1(specificHeatRatio - 1.0), eulObj(specificHeatRatio) {

	int i;

	this->cfl = cflnumber; this->tFinal = tmax; this->t = tInit; dt = 0.0;
	this->nx = nx;
	dx = (xr - xl) / nx;

	area = dx;
	k = 5; s = 2;
	nxx = k + k + nx;
	n = nx;

	circleRadMod = pow(dx / 2.0, 1.5);

	cout << "Creating domain... ";
	cell.resize(nxx);

	for (i = 0; i < nxx; i++) {
		cell(i).x = ((double)i - (double)k + 0.5) * dx + xl;
	}

	cout << "done.\n";
	cout << "Number of cells: " << n << "\n";

	t0 = 0.0, t1 = 0.0, vxm = 0.0, hm = 0.0, qm = 0.0, cm = 0.0;
	rcm = 0.0, b1 = 0.0, b2 = 0.0, t2 = 0.0;
}

void FV::System1D::entropyWave() {
	int i;
	double rho = 1.0, u = 1.0, p = 1.0;
	double A = 0.2;

	cout << "Initializing variable... case: Entropy Waves...\n";
	for (i = 0; i < nxx; i++) {
		cell(i).u = eulObj.entropyWaves(rho, u, p, A, cell(i).x,
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

void FV::System1D::blastWave() {

	for (int i = 0; i < nxx; i++) {
		cell(i).u = eulObj.conservative(eulObj.shuOsher(cell(i).x));
	}
	zeroGradBc();
	rhs();
}

void FV::System1D::periodicBc() {

	for (int m = 0; m < k; m++) {
		cell(k - 1 - m).u = cell(nx + k - 1 - m).u;
		cell(nxx - 1 - m).u = cell(k + k - 1 - m).u;
	}
}

void FV::System1D::zeroGradBc() {
	for (int m = 0; m < k; m++) {
		cell(k - 1 - m).u = cell(k).u;
		cell(nxx - 1 - m).u = cell(nxx - k - 1).u;
//        cell(nx + k + m).u = eulObj.conservative(eulObj.shuOsher(cell(nx + k + m).x));
	}
}

void FV::System1D::reflectiveBc() {

	for (int m = 0; m < k; m++) {
		cell(k - 1 - m).u(0) = cell(k + m).u(0);
		cell(k - 1 - m).u(1) = -cell(k + m).u(1);
		cell(k - 1 - m).u(2) = cell(k + m).u(2);

		cell(nxx - 1 - m).u(0) = cell(nx - 1 + m).u(0);
		cell(nxx - 1 - m).u(1) = -cell(nx - 1 + m).u(1);
		cell(nxx - 1 - m).u(2) = cell(nx - 1 + m).u(2);
	}
}

void FV::System1D::calcPrimitiveVariable() {

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

void FV::System1D::eigen1DEuler() {
	int i;
//	aM = -aM.setOnes();

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

//		aM = aM.cwiseMax(eulObj.localWaveSpeed(cell(i).pvar));
	}
}

void FV::System1D::rhs() {
	int i, m, eq;
//	periodicBc();
	zeroGradBc();
//	reflectiveBc();
	calcPrimitiveVariable();
	eigen1DEuler();

	// Reset
	for (i = s + 1; i < nxx - s - 1; ++i) {
	    cell(i).discontinue = false;
	}

	// KXRCF
	for (i = s + 1; i < nxx - s - 1; ++i) {
	    if (!cell(i).discontinue) {
	        if (cell(i).u(1) < 0.0) {
	            h1 = (2.0 * cell(i-1).u(0) + 5.0 * cell(i).u(0) - cell(i+1).u(0)) / 6.0;
	            h2 = (-cell(i-2).u(0) + 5.0 * cell(i-1).u(0) + 2.0 * cell(i).u(0)) / 6.0;
	        } else {
	            h2 = (2.0 * cell(i).u(0) + 5.0 * cell(i+1).u(0) - cell(i+2).u(0)) / 6.0;
	            h1 = (-cell(i-1).u(0) + 5.0 * cell(i).u(0) + 2.0 * cell(i+1).u(0)) / 6.0;
	        }

	        cell(i).discontinue = fabs(h1 - h2) / (circleRadMod * fabs(cell(i).u(0))) > 1.0;

	        if (cell(i).discontinue){
//	            cell(i-2).discontinue = true;
//	            cell(i-1).discontinue = true;
//	            cell(i+1).discontinue = true;
//	            cell(i+2).discontinue = true;
	        }
	    } else
	        continue;
	}

	for (i = s; i < nxx - s - 1; ++i) {
		for (m = 0; m < 6; ++m) {
			interpolateA.col(m) = cell(i).leftEigenVectorX * cell(i - 2 + m).u;
		}

		for (eq = 0; eq < 3; eq++) {
            valueA = cell(i).discontinue ? lagrangeObj.teno5i(interpolateA.row(eq)) : lagrangeObj.teno5i(
                    interpolateA.row(eq));

			cell(i).uxl(eq) = valueA(0); cell(i).uxr(eq) = valueA(1);
		}
	}

	for (i = k - 1; i < nxx - k; ++i) {
		cell(i).transformToConservative();
	}

	for (i = k - 1; i < nxx - k; ++i) {
		cell(i).setZeroFlux();
		pvarL = eulObj.primitive(cell(i).uxl);
		pvarR = eulObj.primitive(cell(i).uxr);
//        std::cout << i << " : " << cell(i).uxr.transpose() << "\n";
		/*try
		{
			if (pvarL(0) <= 0.0 || pvarL(2) < 0.0)
				throw (-10);
			if (pvarR(0) <= 0.0 || pvarR(2) < 0.0)
				throw (-11);
		}
		catch (const int& x)
		{
			cout << "Exception occurred, exception number: " << x << endl;
			exit(-1);
		}*/
//        aM = eulObj.maxLocalWaveSpeed(pvarL.cwiseMax(pvarR), pvarR.cwiseMin(pvarL));
		aM = eulObj.maxLocalWaveSpeed(pvarL, pvarR);
		cell(i).fx = eulObj.Rusanov(cell(i).uxl, cell(i).uxr, pvarL, pvarR,
			cell(i).leftEigenVectorX, cell(i).rightEigenVectorX, aM);

//		cell(i).fx = eulObj.HLLC(cell(i).uxl, cell(i).uxr);
//        cell(i).fx = eulObj.CUFlux(cell(i).uxl, cell(i).uxr, cell(i).leftEigenVectorX, cell(i).rightEigenVectorX);
	}

	for (i = k; i < nxx - k; i++) {
		cell(i).rhsu = (cell(i - 1).fx - cell(i).fx) / area;
//		std::cout << cell(i).fx.transpose() << "\n";
	}
}

void FV::System1D::RungeKutta3() {
	int i;
	int stage;

	stage = 1;
	rhs();
	for (i = k; i < nxx - k; i++) {
		cell(i).save();
		cell(i).u1 = ssprkObj.RungeKuttaIII<Vector3d>(cell(i).u0,
				cell(i).u0, cell(i).rhsu, dt, stage);
		cell(i).update1();
	}

	stage += 1;
	rhs();
	for (i = k; i < nxx - k; i++) {
		cell(i).u2 = ssprkObj.RungeKuttaIII<Vector3d>(cell(i).u0,
				cell(i).u1, cell(i).rhsu, dt, stage);
		cell(i).update2();
	}

	stage += 1;
	rhs();
	for (i = k; i < nxx - k; i++) {
		cell(i).u = ssprkObj.RungeKuttaIII<Vector3d>(cell(i).u0,
				cell(i).u2, cell(i).rhsu, dt, stage);
	}

	calcPrimitiveVariable();
}

void FV::System1D::timeStepping() {
	int i;
	double maxSpeedx;
	calcPrimitiveVariable();

	while (t < tFinal) {
		maxSpeedx = 0.0;
		for (i = k; i < nxx - k; i++) {
			maxSpeedx = fmax(maxSpeedx, eulObj.maxSpeed(cell(i).pvar));
		}

		dt = cfl * pow(dx, 3.0 / 3.0) / maxSpeedx;
		dt = t + dt > tFinal ? tFinal - t : dt;

		RungeKutta3();
		t += dt;
		cout << "tnum: " << t << endl;
	}
}

void FV::System1D::printResult() {
	int i;
	double err, merr;

	err = 0.0;
	merr = -1.0;

	for (i = k; i < nxx - k; i++) {
		err += fabs(cell(i).u(0) - cell(i).exact(0));
		merr = fmax(merr, fabs(cell(i).u(0) - cell(i).exact(0)));

	}

	cout << "Numerical error:" << err / n << endl;
	cout <<"L inf: " << merr << endl;
	write(0);
}

void FV::System1D::write(const int& filenum) {
    std::string filename = "/home/indro/Documents/1D_SOD_KXRCF_TENO_" + std::to_string(nx) +
			"_" + std::to_string(filenum) + ".dat";
	std::ofstream outfile;
	int i;
	const char separator = ' ';
	const int strwidth = 30;
	const int numwidth = 30;
	const int pn = 15;

	outfile.open(filename);
	if (outfile.good()) {
		outfile << "t: " << t << "n";
		outfile << std::right << std::setw(strwidth) << std::setfill(separator) << "x";
		outfile << std::right << std::setw(strwidth) << std::setfill(separator) << "Density";
		outfile << std::right << std::setw(strwidth) << std::setfill(separator) << "Ux";
		outfile << std::right << std::setw(strwidth) << std::setfill(separator) << "Pressure";
		outfile << "\n";
	}
	else {
		exit(1);
	}

	for (i = 0; i < nxx; ++i) {
		outfile << std::right << std::setw(numwidth) << std::setfill(separator) << std::setprecision(pn) << cell(i).x;
		outfile << std::right << std::setw(numwidth) << std::setfill(separator) << std::setprecision(pn) << cell(i).pvar(0);
		outfile << std::right << std::setw(numwidth) << std::setfill(separator) << std::setprecision(pn) << cell(i).pvar(1);
		outfile << std::right << std::setw(numwidth) << std::setfill(separator) << std::setprecision(pn) << cell(i).pvar(2);
		outfile << std::right << std::setw(numwidth) << std::setfill(separator) << std::setprecision(pn)
		<< cell(i).discontinue;
		outfile << "\n";
	}

	std::cout << "\nSolution is written in " << filename << "\n";
	outfile.close();
}
