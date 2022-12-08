#include "Equations_Euler2D.h"

Equations::Euler2D::Euler2D(const double &specificHeatRatio) :
        gamma(specificHeatRatio), gamma_1(specificHeatRatio - 1.0),
        gamma_3(specificHeatRatio - 3.0), pi(4.0 * atan(1.0)) {

    squareVelocityx = 0.0, dynVelocityx = 0.0, enthalpyx = 0.0;
    squareVelocityy = 0.0, dynVelocityy = 0.0, enthalpyy = 0.0;

    soundSpeedLx = 0.0, soundSpeedRx = 0.0, velocityLx = 0.0, velocityRx = 0.0;
    soundSpeedLy = 0.0, soundSpeedRy = 0.0, velocityLy = 0.0, velocityRy = 0.0;

    soundSpeed = 0.0;

    soundSpeedLloc = 0.0, soundSpeedRloc = 0.0, velocityLloc = 0.0, velocityRloc = 0.0;
    soundSpeedLocal = 0.0, velocity = 0.0;

}

Vector4d Equations::Euler2D::primitive(const Vector4d &q) {
    wLoc(0) = q(0);
    wLoc(1) = q(1) / q(0);
    wLoc(2) = q(2) / q(0);
    wLoc(3) = (q(3) - 0.5 * (q(1) * wLoc(1) + q(2) * wLoc(2))) * gamma_1;
    return wLoc;
}

Vector4d Equations::Euler2D::conservative(const Vector4d &w) {
    qLoc(0) = w(0);
    qLoc(1) = w(0) * w(1);
    qLoc(2) = w(0) * w(2);
    qLoc(3) = w(3) / gamma_1 + 0.5 * (qLoc(1) * w(1) + qLoc(2) * w(2));
    return qLoc;
}

Vector4d Equations::Euler2D::fluxx(const Vector4d &q, const Vector4d &w) {
    fxLoc(0) = q(1);
    fxLoc(1) = q(1) * w(1) + w(3);
    fxLoc(2) = q(1) * w(2);
    fxLoc(3) = w(1) * (q(3) + w(3));
    return fxLoc;
}

Vector4d Equations::Euler2D::fluxy(const Vector4d &q, const Vector4d &w) {
    fyLoc(0) = q(2);
    fyLoc(1) = q(2) * w(1);
    fyLoc(2) = q(2) * w(2) + w(3);
    fyLoc(3) = w(2) * (q(3) + w(3));
    return fyLoc;
}

Matrix4d Equations::Euler2D::dfluxx(const Vector4d &q, const Vector4d &w) {
    squareVelocityx = w(1) * w(1);
    dynVelocityx = 0.5 * gamma_1 * (squareVelocityx + w(2) * w(2));
    enthalpyx = (q(3) + w(3)) / w(0);

    dfxLoc(0, 0) = 0.0;
    dfxLoc(0, 1) = 1.0;
    dfxLoc(0, 2) = 0.0;
    dfxLoc(0, 3) = 0.0;
    dfxLoc(1, 0) = dynVelocityx - squareVelocityx;
    dfxLoc(1, 1) = -gamma_3 * w(1);
    dfxLoc(1, 2) = -gamma_1 * w(2);
    dfxLoc(1, 3) = gamma_1;
    dfxLoc(2, 0) = -w(1) * w(2);
    dfxLoc(2, 1) = w(2);
    dfxLoc(2, 2) = w(1);
    dfxLoc(2, 3) = 0.0;
    dfxLoc(3, 0) = (dynVelocityx - enthalpyx) * w(1);
    dfxLoc(3, 1) = enthalpyx - gamma_1 * squareVelocityx;
    dfxLoc(3, 2) = -gamma_1 * w(1) * w(2);
    dfxLoc(3, 3) = gamma * w(1);
    return dfxLoc;
}

Matrix4d Equations::Euler2D::dfluxy(const Vector4d &q, const Vector4d &w) {
    squareVelocityy = w(2) * w(2);
    dynVelocityy = 0.5 * gamma_1 * (squareVelocityy + w(1) * w(1));
    enthalpyy = (q(3) + w(3)) / w(0);

    dfyLoc(0, 0) = 0.0;
    dfyLoc(0, 1) = 0.0;
    dfyLoc(0, 2) = 1.0;
    dfyLoc(0, 3) = 0.0;
    dfyLoc(1, 0) = -w(1) * w(2);
    dfyLoc(1, 1) = w(2);
    dfyLoc(1, 2) = w(1);
    dfyLoc(1, 3) = 0.0;
    dfyLoc(2, 0) = dynVelocityy - squareVelocityy;
    dfyLoc(2, 1) = -gamma_1 * w(1);
    dfyLoc(2, 2) = -gamma_3 * w(2);
    dfyLoc(2, 3) = gamma_1;
    dfyLoc(3, 0) = (dynVelocityy - enthalpyy) * w(2);
    dfyLoc(3, 1) = -gamma_1 * w(1) * w(2);
    dfyLoc(3, 2) = enthalpyy - gamma_1 * squareVelocityy;
    dfyLoc(3, 3) = gamma * w(2);
    return dfyLoc;
}

Matrix4d Equations::Euler2D::rightEigenVectorx(const double &vxm,
                                               const double &vym, const double &hm,
                                               const double &qm, const double &cm, const double &t0) {

    eigRx(0, 0) = 1.0;
    eigRx(0, 1) = 0.0;
    eigRx(0, 2) = 1.0;
    eigRx(0, 3) = 1.0;

    eigRx(1, 0) = vxm - cm;
    eigRx(1, 1) = 0.0;
    eigRx(1, 2) = vxm;
    eigRx(1, 3) = vxm + cm;

    eigRx(2, 0) = vym;
    eigRx(2, 1) = 1.0;
    eigRx(2, 2) = vym;
    eigRx(2, 3) = vym;

    eigRx(3, 0) = hm - t0;
    eigRx(3, 1) = vym;
    eigRx(3, 2) = qm;
    eigRx(3, 3) = hm + t0;

    return eigRx;
}

Matrix4d Equations::Euler2D::leftEigenVectorx(const double &rcm,
                                              const double &vym, const double &b1, const double &b2,
                                              const double &t0, const double &t1, const double &t2) {
    eigLx(0, 0) = 0.5 * (b2 + t0);
    eigLx(0, 1) = -0.5 * (t1 + rcm);
    eigLx(0, 2) = -t2 * vym;
    eigLx(0, 3) = t2;

    eigLx(1, 0) = -vym;
    eigLx(1, 1) = 0.0;
    eigLx(1, 2) = 1.0;
    eigLx(1, 3) = 0.0;

    eigLx(2, 0) = 1.0 - b2;
    eigLx(2, 1) = t1;
    eigLx(2, 2) = b1 * vym;
    eigLx(2, 3) = -b1;

    eigLx(3, 0) = 0.5 * (b2 - t0);
    eigLx(3, 1) = -0.5 * (t1 - rcm);
    eigLx(3, 2) = -t2 * vym;
    eigLx(3, 3) = t2;

    return eigLx;
}

Matrix4d Equations::Euler2D::rightEigenVectory(const double &vxm,
                                               const double &vym, const double &hm,
                                               const double &qm, const double &cm, const double &t0) {

    eigRy(0, 0) = 1.0;
    eigRy(0, 1) = 0.0;
    eigRy(0, 2) = 1.0;
    eigRy(0, 3) = 1.0;

    eigRy(1, 0) = vxm;
    eigRy(1, 1) = 1.0;
    eigRy(1, 2) = vxm;
    eigRy(1, 3) = vxm;

    eigRy(2, 0) = vym - cm;
    eigRy(2, 1) = 0.0;
    eigRy(2, 2) = vym;
    eigRy(2, 3) = vym + cm;

    eigRy(3, 0) = hm - t0;
    eigRy(3, 1) = vxm;
    eigRy(3, 2) = qm;
    eigRy(3, 3) = hm + t0;

    return eigRy;
}

Matrix4d Equations::Euler2D::leftEigenVectory(const double &rcm,
                                              const double &vxm, const double &b1, const double &b2, const double &t0,
                                              const double &t1, const double &t2) {

    eigLy(0, 0) = 0.5 * (b2 + t0);
    eigLy(0, 1) = -t2 * vxm;
    eigLy(0, 2) = -0.5 * (t1 + rcm);
    eigLy(0, 3) = t2;

    eigLy(1, 0) = -vxm;
    eigLy(1, 1) = 1.0;
    eigLy(1, 2) = 0.0;
    eigLy(1, 3) = 0.0;

    eigLy(2, 0) = 1.0 - b2;
    eigLy(2, 1) = b1 * vxm;
    eigLy(2, 2) = t1;
    eigLy(2, 3) = -b1;

    eigLy(3, 0) = 0.5 * (b2 - t0);
    eigLy(3, 1) = -t2 * vxm;
    eigLy(3, 2) = -0.5 * (t1 - rcm);
    eigLy(3, 3) = t2;

    return eigLy;
}

double Equations::Euler2D::maxSpeed(const Vector4d &w, int dir) {
    soundSpeed = sqrt(fabs(gamma * w(3) / w(0)));
    try {
        if (dir != 1 && dir != 2) {
            throw 2;
        }
    }
    catch (int x) {
        cout << "Exception occurred, exception number: " << x << endl;
    }

    return fmax((w(dir) - soundSpeed), fabs(w(dir) + soundSpeed));
}

Vector4d Equations::Euler2D::maxLocalWaveSpeed(const Vector4d &wl, const Vector4d &wr, int dir) {
    /*try {
        if (wl(0) <= 0.0 || wr(0) <= 0.0 || wl(3) < 0.0 || wr(3) < 0.0) {
            throw (-5);
        }
        else if (isnan(wl(0)) || isnan(wr(0))) {
            throw (-100);
        }
    }
    catch (int x) {
        cout << "Exception occurred, exception number: " << x << endl;
        exit(-1);
    }*/

    soundSpeedLloc = sqrt(fabs(gamma * wl(3) / wl(0)));
    soundSpeedRloc = sqrt(fabs(gamma * wr(3) / wr(0)));
    velocityLloc = wl(dir);
    velocityRloc = wr(dir);

    maxWave(0) = fmax(fabs(velocityLloc - soundSpeedLloc), fabs(velocityRloc - soundSpeedRloc));
    maxWave(1) = fmax(fabs(velocityLloc), fabs(velocityRloc));
    maxWave(2) = maxWave(1);
    maxWave(3) = fmax(fabs(velocityLloc + soundSpeedLloc), fabs(velocityRloc + soundSpeedRloc));

    return maxWave;
}

Matrix4d Equations::Euler2D::maxLocalWaveSpeed(const Vector4d& wl, const Vector4d& wr, const Matrix4d& leftEigenVec, const Matrix4d& rightEigenVec, int dir) {
    /*try {
        if (wl(0) <= 0.0 || wr(0) <= 0.0 || wl(3) < 0.0 || wr(3) < 0.0) {
            throw (-5);
        }
        else if (isnan(wl(0)) || isnan(wr(0))) {
            throw (-100);
        }
    }
    catch (int x) {
        cout << "Exception occurred, exception number: " << x << endl;
        exit(-1);
    }*/

    soundSpeedLloc = sqrt(fabs(gamma * wl(3) / wl(0)));
    soundSpeedRloc = sqrt(fabs(gamma * wr(3) / wr(0)));
    velocityLloc = wl(dir);
    velocityRloc = wr(dir);

    waveSpeed(0) = fabs(velocityLloc - soundSpeedLloc);
    waveSpeed(1) = fabs(velocityLloc);
    waveSpeed(2) = waveSpeed(1);
    waveSpeed(3) = fabs(velocityLloc + soundSpeedLloc);

    dfdu = rightEigenVec * waveSpeed.asDiagonal() * leftEigenVec;

    waveSpeed(0) = fabs(velocityRloc - soundSpeedRloc);
    waveSpeed(1) = fabs(velocityRloc);
    waveSpeed(2) = waveSpeed(1);
    waveSpeed(3) = fabs(velocityRloc + soundSpeedRloc);

    dfdu = dfdu.cwiseMax(rightEigenVec * waveSpeed.asDiagonal() * leftEigenVec);

    return dfdu;
}

Vector4d Equations::Euler2D::localWaveSpeed(const Vector4d &w, int dir) {
//    try {
//        if (w(0) <= 0.0 || w(3) < 0.0)
//            throw (-5);
//        else if (dir != 1 && dir != 2)
//            throw 2;
//        else if (isnan(w(0))) {
//            throw (-100);
//        }
//    }
//    catch (int x) {
//        cout << "Exception occurred, exception number: " << x << endl;
//        exit(-1);
//    }

    soundSpeedLocal = sqrt(fabs(gamma * w(3) / w(0)));
    velocity = w(dir);

    waveSpeed(0) = fabs(velocity - soundSpeedLocal);
    waveSpeed(1) = fabs(velocity);
    waveSpeed(2) = waveSpeed(1);
    waveSpeed(3) = fabs(velocity + soundSpeedLocal);

    return waveSpeed;
}

Vector4d Equations::Euler2D::Rusanov(const Vector4d &ql, const Vector4d &qr, const Vector4d &wl, const Vector4d &wr,
                                     const Matrix4d &leftEigenVec, const Matrix4d &rightEigenVec,
                                     const Vector4d &maxSpeed, int dir) {
    try {
        if (dir != 1 && dir != 2)
            throw 1;
    }
    catch (int x) {
        cout << "Exception occurred, exception number: " << x << endl;
    }

    if (dir == 1) {
        rusanovFlux = 0.5 * (fluxx(ql, wl) + fluxx(qr, wr) -
                             rightEigenVec * maxSpeed.cwiseProduct(leftEigenVec * (qr - ql)));
    } else if (dir == 2) {
        rusanovFlux = 0.5 * (fluxy(ql, wl) + fluxy(qr, wr) -
                             rightEigenVec * maxSpeed.cwiseProduct(leftEigenVec * (qr - ql)));
    }

    return rusanovFlux;
}

Vector4d Equations::Euler2D::Rusanov(const Vector4d &ql, const Vector4d &qr, const Vector4d &wl, const Vector4d &wr,
                                     const Matrix4d &maxSpeed, int dir) {
    try {
        if (dir != 1 && dir != 2)
            throw 1;
    }
    catch (int x) {
        cout << "Exception occurred, exception number: " << x << endl;
    }

    if (dir == 1) {
        rusanovFlux = 0.5 * (fluxx(ql, wl) + fluxx(qr, wr) -
                             maxSpeed*(qr - ql));
    } else if (dir == 2) {
        rusanovFlux = 0.5 * (fluxy(ql, wl) + fluxy(qr, wr) -
                             maxSpeed*(qr - ql));
    }
    return rusanovFlux;
}

Vector4d Equations::Euler2D::dRusanov(const Vector4d &dql, const Vector4d &dqr, const Vector4d &wl,
                                      const Vector4d &wr, const Matrix4d &leftEigenVec, const Matrix4d &rightEigenVec,
                                      const Matrix4d &dfluxl, const Matrix4d &dfluxr, const Vector4d &maxSpeed) {

    rusanovFluxH = 0.5 * (dfluxl * dql + dfluxr * dqr -
                          rightEigenVec * maxSpeed.cwiseProduct(leftEigenVec * (dqr - dql)));

    return rusanovFluxH;
}

Vector4d Equations::Euler2D::dRusanov(const Vector4d &dql, const Vector4d &dqr, const Vector4d &wl,
                                      const Vector4d &wr, const Matrix4d &dfluxl, const Matrix4d &dfluxr,
                                      const Matrix4d &maxSpeed) {

    rusanovFluxH = 0.5 * (dfluxl * dql + dfluxr * dqr -
                          maxSpeed*(dqr - dql));

    return rusanovFluxH;
}

Vector4d Equations::Euler2D::HLLC(const Vector4d &ql, const Vector4d &qr, int dir) {
    int dir_1 = dir - 1;
    Vector4d fl, fr, wl, wr, qs, hllcflux;

    double rhoRatio, rrhoRatio, hl, hr, hm, cm, qm;
    Vector2d velm;

    wl = primitive(ql);
    wr = primitive(qr);
    hl = (ql(3) + wl(3)) / wl(0);
    hr = (qr(3) + wr(3)) / wr(0);
    rhoRatio = sqrt(fabs(wr(0) / wl(0)));
    rrhoRatio = 1.0 / (1.0 + rhoRatio);

    velm(0) = (wl(1) + wr(1) * rhoRatio) * rrhoRatio;
    velm(1) = (wl(2) + wr(2) * rhoRatio) * rrhoRatio;
    qm = 0.5 * (velm(0) * velm(0) + velm(1) * velm(1));
    hm = (hl + hr * rhoRatio) * rrhoRatio;
    cm = sqrt(fabs(gamma_1 * (hm - qm)));

    double cl, cr, sl, sr, sm, ps;
    cl = sqrt(fabs(wl(3) / wl(0)));
    cr = sqrt(fabs(wr(3) / wr(0)));

    sl = fmin(wl(dir) - cl, velm(dir_1) - cm);
    sr = fmax(wr(dir) + cr, velm(dir_1) + cm);
    sm = (qr(dir) * (sr - wr(dir)) - ql(dir) * (sl - wl(dir)) + wl(3) - wr(3)) /
         (wr(0) * (sr - wr(dir)) - wl(0) * (sl - wl(dir)));
    ps = wl(3) + wl(0) * (wl(dir) - sl) * (wl(dir) - sm);

    fl = dir == 1 ? fluxx(ql, wl) : fluxy(ql, wl);
    fr = dir == 1 ? fluxx(qr, wr) : fluxy(qr, wr);

    double nom, denom;

    if (sl > 0.0) {
        return fl;
    } else if (sl <= 0.0 && sm > 0.0) {
        nom = sl - wl(dir);
        denom = sl - sm;
        qs(0) = wl(0) * nom / denom;
        qs(1) = dir == 1 ? (nom * ql(1) + ps - wl(3)) / denom : nom * ql(1) / denom;
        qs(2) = dir == 2 ? (nom * ql(2) + ps - wl(3)) / denom : nom * ql(2) / denom;
        qs(3) = (nom * ql(3) - wl(3) * wl(dir) + ps * sm) / denom;
        return fl + sl * (qs - ql);
    } else if (sm <= 0.0 && sr >= 0.0) {
        nom = sr - wr(dir);
        denom = sr - sm;
        qs(0) = wr(0) * nom / denom;
        qs(1) = dir == 1 ? (nom * qr(1) + ps - wr(3)) / denom : nom * qr(1) / denom;
        qs(2) = dir == 2 ? (nom * qr(2) + ps - wr(3)) / denom : nom * qr(2) / denom;
        qs(3) = (nom * qr(3) - wl(3) * wr(dir) + ps * sm) / denom;
        return fr + sr * (qs - qr);
    } else {
        return fr;
    }
}

Vector4d Equations::Euler2D::conservativeOfEntropyWave(const double &rhoInf, const double &uInf, const double &vInf,
                                                       const double &pInf, const double &amplitude, const double &x,
                                                       const double &y) {

    Vector4d w;
    w(0) = rhoInf + amplitude * sin(pi * (x + y));
    w(1) = uInf;
    w(2) = vInf;
    w(3) = pInf;
    return conservative(w);
}

Vector4d Equations::Euler2D::dconservativeOfEntropyWave(const double &rhoInf, const double &uInf, const double &vInf,
                                                        const double &pInf, const double &amplitude, const double &x,
                                                        const double &y) {

    Vector4d q;
    q(0) = amplitude * cos(pi * (x + y)) * pi;
    q(1) = uInf * q(0);
    q(2) = vInf * q(0);
    q(3) = (0.5 * (uInf * uInf + vInf * vInf)) * q(0);
    return q;
}

Vector4d Equations::Euler2D::entropyWaves(const double &rhoInf, const double &uInf, const double &vInf,
                                          const double &pInf, const double &amplitude, const double &x, const double &y,
                                          const double &dx, const double &dy, const Matrix<double, 9, 1> &weights,
                                          const Matrix<double, 9, 1> &absis) {

    int i, j;
    initEntWave.setZero();

    for (i = 0; i < 9; i++) {
        for (j = 0; j < 9; j++) {
            initEntWave += weights(i) * weights(j) *
                           conservativeOfEntropyWave(rhoInf, uInf, vInf, pInf, amplitude,
                                                     x + dx * absis(i) / 2.0, y + dy * absis(j) / 2.0) / 4.0;
        }
    }

    return initEntWave;
}

Vector4d Equations::Euler2D::entropyWavesdx(const double &rhoInf, const double &uInf, const double &vInf,
                                            const double &pInf, const double &amplitude, const double &x,
                                            const double &y,
                                            const double &dx, const double &dy, const Matrix<double, 9, 1> &weights,
                                            const Matrix<double, 9, 1> &absis) {

    int i, j;
    outhx.setZero();

    for (i = 0; i < 9; i++) {
        for (j = 0; j < 9; j++) {
            outhx += weights(i) * weights(j) *
                     this->dconservativeOfEntropyWave(rhoInf, uInf, vInf, pInf, amplitude,
                                                      x + dx * absis(i) / 2.0, y + dy * absis(j) / 2.0) / 4.0;
        }
    }

    return outhx;
}

Vector4d Equations::Euler2D::entropyWavesdy(const double &rhoInf, const double &uInf, const double &vInf,
                                            const double &pInf, const double &amplitude, const double &x,
                                            const double &y,
                                            const double &dx, const double &dy, const Matrix<double, 9, 1> &weights,
                                            const Matrix<double, 9, 1> &absis) {

    int i, j;
    outhy.setZero();

    for (i = 0; i < 9; i++) {
        for (j = 0; j < 9; j++) {
            outhy += weights(i) * weights(j) *
                     this->dconservativeOfEntropyWave(rhoInf, uInf, vInf, pInf, amplitude,
                                                      x + dx * absis(i) / 2.0, y + dy * absis(j) / 2.0) / 4.0;
        }
    }

    return outhy;
}

Vector4d Equations::Euler2D::isentropicVortex(const double &x, const double &y, const double &dx, const double &dy,
                                              const Matrix<double, 9, 1> &weights, const Matrix<double, 9, 1> &absis)
{
    int i, j;
    initEntWave.setZero();

    for (i = 0; i < 9; i++) {
        for (j = 0; j < 9; j++) {
            initEntWave += weights(i) * weights(j) *
                    this->conservativeOfIsentropicVortex(x + dx * absis(i) / 2.0, y + dy * absis(j) / 2.0) / 4.0;
        }
    }
//    initEntWave = conservativeOfIsentropicVortex(x, y);

    return initEntWave;
}

Vector4d Equations::Euler2D::conservativeOfRiemann2D(const double &x, const double &y) {
    Vector4d w0;

    if ((x > 0.5) && (y > 0.5)) {
        w0(0) = 1.0;
        w0(1) = 0.75;
        w0(2) = -0.5;
        w0(3) = 1.0;

        //w0(0) = 1.5;
        //w0(1) = 0.0;
        //w0(2) = 0.0;
        //w0(3) = 1.5;

        //w0(0) = 1.0;
        //w0(1) = -0.75;
        //w0(2) = -0.5;
        //w0(3) = 1.0;

//		w0(0) = 1.0;
//		w0(1) = 0.0;
//		w0(2) = -0.4;
//		w0(3) = 1.0;
    } else if ((x < 0.5) && (y > 0.5)) {
        w0(0) = 2.0;
        w0(1) = 0.75;
        w0(2) = 0.5;
        w0(3) = 1.0;

        //w0(0) = 0.5323;
        //w0(1) = 1.206;
        //w0(2) = 0.0;
        //w0(3) = 0.3;

        //w0(0) = 2.0;
        //w0(1) = -0.75;
        //w0(2) = 0.5;
        //w0(3) = 1.0;

//		w0(0) = 2.0;
//		w0(1) = 0.0;
//		w0(2) = -0.3;
//		w0(3) = 1.0;

    } else if ((x < 0.5) && (y < 0.5)) {
        w0(0) = 1.0;
        w0(1) = -0.75;
        w0(2) = 0.5;
        w0(3) = 1.0;

        //w0(0) = 0.138;
        //w0(1) = 1.206;
        //w0(2) = 1.206;
        //w0(3) = 0.029;

        //w0(0) = 1.0;
        //w0(1) = 0.75;
        //w0(2) = 0.5;
        //w0(3) = 1.0;

//		w0(0) = 1.0625;
//		w0(1) = 0.0;
//		w0(2) = 0.2145;
//		w0(3) = 0.4;
    } else {
        w0(0) = 3.0;
        w0(1) = -0.75;
        w0(2) = -0.5;
        w0(3) = 1.0;

        //w0(0) = 0.5323;
        //w0(1) = 0.0;
        //w0(2) = 1.206;
        //w0(3) = 0.3;

        //w0(0) = 3.0;
        //w0(1) = 0.75;
        //w0(2) = -0.5;
        //w0(3) = 1.0;

//		w0(0) = 0.5197;
//		w0(1) = 0.0;
//		w0(2) = -1.1259;
//		w0(3) = 0.4;
    }

    return conservative(w0);
}

Vector4d Equations::Euler2D::conservativeOfTaylorIns2D(const double &x, const double &y) {
    Vector4d w0;

    if (y < 0.5) {
        w0(0) = 2.0;
        w0(1) = 0.0;
        w0(3) = 1.0 + 2.0 * y;
        w0(2) = -0.025 * sqrt(gamma * w0(3) / w0(0)) * cos(8.0 * pi * x);
    } else {
        w0(0) = 1.0;
        w0(1) = 0.0;
        w0(3) = y + 1.5;
        w0(2) = -0.025 * sqrt(gamma * w0(3) / w0(0)) * cos(8.0 * pi * x);
    }

    return conservative(w0);
}

Vector4d Equations::Euler2D::conservativeOfDoubleMach2D(const double &x, const double &y) {
    Vector4d w0;
    double lim = 1.0 / 6.0 + y / sqrt(3.0);

    if (x < lim) {
        w0(0) = 8.0;
        w0(1) = 8.25 * 0.5 * sqrt(3.0);
        w0(2) = -8.25 * 0.5;
        w0(3) = 116.5;
    } else {
        w0(0) = 1.4;
        w0(1) = 0.0;
        w0(2) = 0.0;
        w0(3) = 1.0;
    }

    return conservative(w0);
}

Vector4d Equations::Euler2D::conservativeOfIsentropicVortex(const double &x, const double &y) {

    Vector4d w;
    double r = pow(x, 2) + pow(y, 2), b=5.0;
    w(0) = 1.0 - gamma_1 * pow(b,2) / 8.0 / gamma / pow(pi,2) * exp(1 - r);
    w(0) = pow(w(0), 1.0 / gamma_1);
    w(1) = 1.0 - b * y / 2.0 / pi * exp(0.5 * (1 - r));
    w(2) = 1.0 + b * x / 2.0 / pi * exp(0.5 * (1 - r));;
    w(3) = pow(w(0), gamma);
    return conservative(w);
}

Vector4d Equations::Euler2D::dxconservativeOfIsentropicVortex(const double &x, const double &y) {
    Vector4d w, q;
    double r = pow(x, 2) + pow(y, 2);

    w(0) = 1.0 - gamma_1 * 25.0 / 8.0 / gamma / pow(pi,2) * exp(1 - r);
    w(0) = pow(w(0), 1.0 / gamma_1);
    w(1) = 1.0 - 5.0 * y / 2.0 / pi * exp(0.5 * (1 - r));
    w(2) = 1.0 + 5.0 * x / 2.0 / pi * exp(0.5 * (1 - r));;
    w(3) = pow(w(0), gamma);

    q(0) = 1.0 - gamma_1 * 25.0 / 8.0 / gamma / pi * exp(1 - r);
    q(0) = 25.0 * x * exp(1 - r) / 4.0 / pow(pi, 2) /
           gamma * pow(q(0), (2.0 - gamma) / gamma_1);
    q(1) = 5.0 * x * y / 2.0 / pi * exp(0.5 * (1 - r));
    q(2) = 5.0 * (1 - pow(x, 2)) / 2.0 / pi * exp(0.5 * (1 - r));
    q(3) = gamma * pow(w(0), gamma_1) * q(0) / gamma_1 +
            0.5 * q(0) * (pow(w(1), 2) + pow(w(2), 2)) +
           0.5 * w(0) * (2.0 * w(1) * q(1) + 2.0 * w(2) * q(2));

    q(1) = w(0) * q(1) + q(0) * w(1);
    q(2) = w(0) * q(2) + q(0) * w(2);
    return q;
}

Vector4d Equations::Euler2D::dyconservativeOfIsentropicVortex(const double &x, const double &y) {
    Vector4d w, q;
    double r = pow(x, 2) + pow(y, 2);

    w(0) = 1.0 - gamma_1 * 25.0 / 8.0 / gamma / pow(pi,2) * exp(1 - r);
    w(0) = pow(w(0), 1.0 / gamma_1);
    w(1) = 1.0 - 5.0 * y / 2.0 / pi * exp(0.5 * (1 - r));
    w(2) = 1.0 + 5.0 * x / 2.0 / pi * exp(0.5 * (1 - r));;
    w(3) = pow(w(0), gamma);

    q(0) = 1.0 - gamma_1 * 25.0 / 8.0 / gamma / pow(pi,2) * exp(1 - r);
    q(0) = 25.0 * y * exp(1 - r) / 4.0 / pow(pi, 2) /
           gamma * pow(q(0), (2.0 - gamma) / gamma_1);
    q(1) = 5.0 * (pow(y,2) - 1.0) / 2.0 / pi * exp(0.5 * (1 - r));
    q(2) = -5.0 * x * y / 2.0 / pi * exp(0.5 * (1 - r));
    q(3) = gamma * pow(w(0), gamma_1) * q(0) / gamma_1 + 0.5 * q(0) * (pow(w(1), 2) + pow(w(2), 2)) +
           0.5 * w(0) * (2.0 * w(1) * q(1) + 2.0 * w(2) * q(2));

    q(1) = w(0) * q(1) + q(0) * w(1);
    q(2) = w(0) * q(2) + q(0) * w(2);
    return q;
}

Vector4d Equations::Euler2D::dxIsentropicVortex(const double &x, const double &y, const double &dx, const double &dy,
                                                const Matrix<double, 9, 1> &weights,
                                                const Matrix<double, 9, 1> &absis) {
    int i, j;
    outhx.setZero();

    for (i = 0; i < 9; i++) {
        for (j = 0; j < 9; j++) {
            outhx += weights(i) * weights(j) *
                     this->dxconservativeOfIsentropicVortex(x + dx * absis(i) / 2.0, y + dy * absis(j) / 2.0) / 4.0;
        }
    }

    return outhx;
}

Vector4d Equations::Euler2D::dyIsentropicVortex(const double &x, const double &y, const double &dx, const double &dy,
                                                const Matrix<double, 9, 1> &weights,
                                                const Matrix<double, 9, 1> &absis) {
    int i, j;
    outhy.setZero();

    for (i = 0; i < 9; i++) {
        for (j = 0; j < 9; j++) {
            outhy += weights(i) * weights(j) *
                     this->dyconservativeOfIsentropicVortex(x + dx * absis(i) / 2.0, y + dy * absis(j) / 2.0) / 4.0;
        }
    }

    return outhy;
}
