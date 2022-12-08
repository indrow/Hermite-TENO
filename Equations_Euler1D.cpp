#include "Equations_Euler1D.h"

Equations::Euler1D::Euler1D(const double& specificHeatRatio) :
gamma(specificHeatRatio), gamma_1(specificHeatRatio - 1.0),
gamma_3(specificHeatRatio - 3.0), pi(4.0 * atan(1.0)) {

	squareVelocityx = 0.0, dynVelocityx = 0.0, enthalpyx = 0.0;

	soundSpeedLx = 0.0, soundSpeedRx = 0.0, velocityLx = 0.0, velocityRx = 0.0;

	soundSpeed = 0.0;

	soundSpeedLloc = 0.0, soundSpeedRloc = 0.0, velocityLloc = 0.0, velocityRloc = 0.0;
	soundSpeedLocal = 0.0, velocity = 0.0;

}

Vector3d Equations::Euler1D::primitive(const Vector3d& q) {
	wLoc(0) = q(0);
	wLoc(1) = q(1) / q(0);
	wLoc(2) = (q(2) - 0.5 * q(1) * wLoc(1)) * gamma_1;
	return wLoc;
}

Vector3d Equations::Euler1D::conservative(const Vector3d& w) {
	qLoc(0) = w(0);
	qLoc(1) = w(0) * w(1);
	qLoc(2) = w(2) / gamma_1 + 0.5 * qLoc(1) * w(1);
	return qLoc;
}

Vector3d Equations::Euler1D::flux(const Vector3d& q, const Vector3d& w) {
	fxLoc(0) = q(1);
	fxLoc(1) = q(1) * w(1) + w(2);
	fxLoc(2) = w(1) * (q(2) + w(2));
	return fxLoc;
}

Matrix3d Equations::Euler1D::dflux(const Vector3d& q, const Vector3d& w) {
	dfxLoc(0, 0) = 0.0;
	dfxLoc(0, 1) = 1.0;
	dfxLoc(0, 2) = 0.0;

	dfxLoc(1, 0) = 0.5 * gamma_3 * w(1) * w(1);
	dfxLoc(1, 1) = -gamma_3 * w(1);
	dfxLoc(1, 2) = gamma_1;

	dfxLoc(2, 0) = (gamma_1 * pow(w(1), 2) - gamma * q(2) / w(0)) * w(1);
	dfxLoc(2, 1) = gamma * q(2) / w(0) - 1.5 * gamma_1 * w(1) * w(1);
	dfxLoc(2, 2) = gamma * w(1);

	return dfxLoc;
}

Matrix3d Equations::Euler1D::rightEigenVector(const double& vxm, const double& hm,
		const double& qm, const double& cm, const double& t0) {

	eigRx(0, 0) = 1.0;
	eigRx(0, 1) = 1.0;
	eigRx(0, 2) = 1.0;

	eigRx(1, 0) = vxm - cm;
	eigRx(1, 1) = vxm;
	eigRx(1, 2) = vxm + cm;

	eigRx(2, 0) = hm - t0;
	eigRx(2, 1) = qm;
	eigRx(2, 2) = hm + t0;

	return eigRx;
}

Matrix3d Equations::Euler1D::leftEigenVector(const double& rcm, const double& b1, const double& b2,
		const double& t0, const double& t1, const double& t2) {

	eigLx(0, 0) = 0.5 * (b2 + t0);
	eigLx(0, 1) = -0.5 * (t1 + rcm);
	eigLx(0, 2) = t2;

	eigLx(1, 0) = 1.0 - b2;
	eigLx(1, 1) = t1;
	eigLx(1, 2) = -b1;

	eigLx(2, 0) = 0.5 * (b2 - t0);
	eigLx(2, 1) = -0.5 * (t1 - rcm);
	eigLx(2, 2) = t2;

	return eigLx;
}

double Equations::Euler1D::maxSpeed(const Vector3d& w) {
	soundSpeed = sqrt(fabs(gamma * w(2) / w(0)));
	return fmax(fabs(w(1) - soundSpeed), fabs(w(1) + soundSpeed));
}

Vector3d Equations::Euler1D::maxLocalWaveSpeed(const Vector3d& wl, const Vector3d& wr) {
	soundSpeedLloc = sqrt(fabs(gamma * wl(2) / wl(0)));
	soundSpeedRloc = sqrt(fabs(gamma * wr(2) / wr(0)));
	velocityLloc = wl(1);
	velocityRloc = wr(1);

	maxWave(0) = fmax(fabs(velocityLloc - soundSpeedLloc), fabs(velocityRloc - soundSpeedRloc));
	maxWave(1) = fmax(fabs(velocityLloc), fabs(velocityRloc));
	maxWave(2) = fmax(fabs(velocityLloc + soundSpeedLloc), fabs(velocityRloc + soundSpeedRloc));

	return maxWave;
}

Vector3d Equations::Euler1D::localWaveSpeed(const Vector3d& w) {
	/*try {
		if (w(0) <= 0.0 || w(2) < 0.0)
			throw (-5);
	}
	catch (int x) {
		cout << "Exception occurred, exception number: " << x << endl;
		exit(-1);
	}*/

	soundSpeedLocal = sqrt(fabs(gamma * w(2) / w(0)));
	velocity = w(1);

	waveSpeed(0) = fabs(velocity - soundSpeedLocal);
	waveSpeed(1) = fabs(velocity);
	waveSpeed(2) = fabs(velocity + soundSpeedLocal);

	return waveSpeed;
}

Vector3d Equations::Euler1D::Rusanov(const Vector3d& ql, const Vector3d& qr, const Vector3d& wl, const Vector3d& wr,
		const Matrix3d& leftEigenVec, const Matrix3d& rightEigenVec, const Vector3d& maxSpeed) {

	rusanovFlux = 0.5 * (flux(ql, wl) + flux(qr, wr) -
			rightEigenVec * maxSpeed.cwiseProduct(leftEigenVec * (qr - ql)));
	return rusanovFlux;
}

Vector3d Equations::Euler1D::CUFlux(const Vector3d& ql, const Vector3d& qr, const Matrix3d& leftEigenVec, const Matrix3d& rightEigenVec) {
	Vector3d qmax, qmin, spdmax, spdmin, aplus, aminus, cuflux;
	qmax = ql.cwiseMax(qr);
	qmin = ql.cwiseMin(qr);

	spdmax = this->localWaveSpeed(this->primitive(qmax));
	spdmin = this->localWaveSpeed(this->primitive(qmin));

	aplus = (spdmax.cwiseMax(spdmin)).cwiseMax(0.0);
	aminus = (spdmax.cwiseMin(spdmin)).cwiseMin(0.0);

	cuflux = (aplus.cwiseProduct(flux(ql,primitive(ql))) -
			aminus.cwiseProduct(flux(qr,primitive(qr)))).cwiseProduct((aplus - aminus).cwiseInverse());

	cuflux  += rightEigenVec * (aplus.cwiseProduct(aminus)).cwiseProduct((aplus - aminus).cwiseInverse()).cwiseProduct(leftEigenVec * (qr - ql));

	return cuflux;
}

Vector3d Equations::Euler1D::HLLC(const Vector3d& ql, const Vector3d& qr) {
	Vector3d fl, fr, wl, wr, qs, hllcflux;

	double rhoRatio, rrhoRatio, hl, hr, vm, hm, cm, qm;

	wl = primitive(ql);
	wr = primitive(qr);
	hl = (ql(2) + wl(2)) / wl(0);
	hr = (qr(2) + wr(2)) / wr(0);
	rhoRatio = sqrt(fabs(wr(0) / wl(0)));
	rrhoRatio = 1.0 / (1.0 + rhoRatio);

	vm = (wl(1) + wr(1) * rhoRatio) * rrhoRatio;
	qm = 0.5 * vm * vm;
	hm = (hl + hr * rhoRatio) * rrhoRatio;
	cm = sqrt(fabs(gamma_1 * (hm - qm)));

	double cl, cr, sl, sr, sm, ps;
	cl = sqrt(fabs(wl(2) / wl(0)));
	cr = sqrt(fabs(wr(2) / wr(0)));

	sl = fmin(wl(1) - cl, vm - cm);
	sr = fmax(wr(1) + cr, vm + cm);
	sm = (qr(1) * (sr - wr(1)) - ql(1) * (sl - wl(1)) + wl(2) - wr(2)) /
			(wr(0) * (sr - wr(1)) - wl(0) * (sl - wl(1)));
	ps = wl(2) + wl(0) * (wl(1) - sl) * (wl(1) - sm);

	fl = flux(ql, wl);
	fr = flux(qr, wr);

	double nom, denom;

	if (sl > 0.0) {
		return fl;
	}
	else if (sl <= 0.0 && sm > 0.0) {
		nom = sl - wl(1);
		denom = sl - sm;
		qs(0) = wl(0) * nom / denom;
		qs(1) = (nom * ql(1) + ps - wl(2)) / denom;
		qs(2) = (nom * ql(2) - wl(2) * wl(1) + ps * sm) / denom;
		return fl + sl * (qs - ql);
	}
	else if (sm <= 0.0 && sr >= 0.0) {
		nom = sr - wr(1);
		denom = sr - sm;
		qs(0) = wr(0) * nom / denom;
		qs(1) = (nom * qr(1) + ps - wr(2)) / denom;
		qs(2) = (nom * qr(2) - wl(2) * wr(1) + ps * sm) / denom;
		return fr + sr * (qs - qr);
	}
	else {
		return fr;
	}
}

Vector3d Equations::Euler1D::dRusanov(const Vector3d& dql, const Vector3d& dqr, const Vector3d& wl,
		const Vector3d& wr, const Matrix3d& leftEigenVec, const Matrix3d& rightEigenVec,
		const Matrix3d& dfluxl, const Matrix3d& dfluxr, const Vector3d& maxSpeed) {

	rusanovFluxH = 0.5 * (dfluxl * dql + dfluxr * dqr -
			rightEigenVec * maxSpeed.cwiseProduct(leftEigenVec * (dqr - dql)));

	return rusanovFluxH;
}

Vector3d Equations::Euler1D::dCUFlux(const Vector3d& ql, const Vector3d& qr, const Matrix3d& leftEigenVec, const Matrix3d& rightEigenVec,
		const Vector3d& dql, const Vector3d& dqr, const Matrix3d& dfluxl, const Matrix3d& dfluxr) {
	Vector3d qmax, qmin, spdmax, spdmin, aplus, aminus, cuflux;
	qmax = ql.cwiseMax(qr);
	qmin = ql.cwiseMin(qr);

	spdmax = this->localWaveSpeed(this->primitive(qmax));
	spdmin = this->localWaveSpeed(this->primitive(qmin));

	aplus = (spdmax.cwiseMax(spdmin)).cwiseMax(0.0);
	aminus = (spdmax.cwiseMin(spdmin)).cwiseMin(0.0);

	cuflux = aplus.cwiseProduct(dfluxl * dql) -
			aminus.cwiseProduct(dfluxr * dqr) +
			rightEigenVec * aplus.cwiseProduct(aminus).cwiseProduct((dqr - dql));

	cuflux = cuflux.cwiseProduct((aplus - aminus).cwiseInverse());

	return cuflux;
}

Vector3d Equations::Euler1D::conservativeOfEntropyWave(const double& rhoInf, const double& uInf,
		const double& pInf, const double& amplitude, const double& x) {

	Vector3d w;
	w(0) = rhoInf + amplitude * sin(pi * (x));
	w(1) = uInf;
	w(2) = pInf;
	return conservative(w);
}

Vector3d Equations::Euler1D::dconservativeOfEntropyWave(const double& rhoInf, const double& uInf,
		const double& pInf, const double& amplitude, const double& x) {

	Vector3d q;
	q(0) = amplitude * cos(pi * (x)) * pi;
	q(1) = uInf * q(0);
	q(2) = (0.5 * (uInf * uInf)) * q(0);
	return q;
}

Vector3d Equations::Euler1D::entropyWaves(const double& rhoInf, const double& uInf,
		const double& pInf, const double& amplitude, const double& x,
		const double& dx, const Matrix<double, 9, 1>& weights, const Matrix<double, 9, 1>& absis) {

	int i;
	initEntWave.setZero();

	for (i = 0; i < 9; i++) {
		initEntWave += weights(i) *
				conservativeOfEntropyWave(rhoInf, uInf, pInf, amplitude,
						x + dx * absis(i) / 2.0) / 2.0;
	}

	return initEntWave;
}

Vector3d Equations::Euler1D::entropyWavesdx(const double& rhoInf, const double& uInf,
		const double& pInf, const double& amplitude, const double& x,
		const double& dx, const Matrix<double, 9, 1>& weights, const Matrix<double, 9, 1>& absis) {

	int i;
	outhx.setZero();

	for (i = 0; i < 9; i++) {
		outhx += weights(i) *
				this->dconservativeOfEntropyWave(rhoInf, uInf, pInf, amplitude,
						x + dx * absis(i) / 2.0) / 2.0;
	}

	return outhx;
}

Vector3d Equations::Euler1D::interactingBlast(const double& x) {
	Vector3d w0;

	if (x < 0.1) {
		w0(0) = 1.0;
		w0(1) = 0.0;
		w0(2) = 1000.0;
	}
	else if (x < 0.9) {
		w0(0) = 1.0;
		w0(1) = 0.0;
		w0(2) = 0.01;
	}
	else {
		w0(0) = 1.0;
		w0(1) = 0.0;
		w0(2) = 100.0;
	}

	return w0;
}

Vector3d Equations::Euler1D::sodShock(const double& x) {
	Vector3d w0;

	if (x < 0.5) {
		w0(0) = 1.0;
		w0(1) = 0.0;
		w0(2) = 1.0;
	}
	else {
		w0(0) = 0.125;
		w0(1) = 0.0;
		w0(2) = 0.1;
	}

	return w0;
}

Vector3d Equations::Euler1D::shuOsher(const double& x) {
	Vector3d w0;

	if (x < -4.0) {
		w0 << 3.857143, 2.629369, 10.333333;
	}
	else {
		w0 << 1.0 + 0.2 * sin(5.0 * x), 0.0, 1.0;
	}

	return w0;
}
