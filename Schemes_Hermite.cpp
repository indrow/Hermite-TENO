#include "Schemes_Hermite.h"

using std::isnan;

Schemes::Hermite::Hermite() :
        eps(1E-6), machine_eps(1E-15), epsz(1E-40), epszd(1E-40), ct(1E-3), ctd(1E-3),
        dif(0.0), c1(13.0 / 12.0), c2(sqrt(3.0) / 3.0), c3(37.0 * sqrt(3.0) / 108.0),
        c4(sqrt(3.0) / 3.0), c5(37.0 * sqrt(3.0) / 108.0),
        c6(sqrt(5.0) / 5.0), c7(61.0 * sqrt(5.0) / 300.0),
        sq3(sqrt(3.0)), sq5(sqrt(5.0)),
        brk(0.0), brkd(0.0), tau(0.0) {
    cout << "Hermite says: please set the grid size using function object.setGridSize(dif)... \n";
}

void Schemes::Hermite::setGridSize(const double &gridSize) {
    try {
        if (gridSize == 0.0)
            throw 21;

        dif = gridSize;
    }
    catch (int x) {
        cout << "Exception occurred, exception number: " << x << "\n";
        exit(1);
    }

    cout << "Hermite says: done.\n";

    double denom;

    crj(0, 0) = -7.0 / 6.0;
    crj(1, 0) = 13.0 / 6.0;
    crj(2, 0) = -2.0 / 3.0 * dif;
    crj(0, 1) = 1.0 / 6.0;
    crj(1, 1) = 5.0 / 6.0;
    crj(2, 1) = -1.0 / 3.0 * dif;
    crj(0, 2) = -1.0 / 6.0;
    crj(1, 2) = 5.0 / 6.0;
    crj(2, 2) = 1.0 / 3.0;

    crjd(0, 0) = 4.0 / dif;
    crjd(1, 0) = -4.0 / dif;
    crjd(2, 0) = 1.5;
    crjd(3, 0) = 3.5;
    crjd(0, 1) = -2.0 / dif;
    crjd(1, 1) = 2.0 / dif;
    crjd(2, 1) = -0.5;
    crjd(3, 1) = -0.5;
    crjd(0, 2) = 0.25 / dif;
    crjd(1, 2) = -1.0 / dif;
    crjd(2, 2) = 0.75 / dif;
    crjd(3, 2) = 0.5;

    dkm(0) = 9.0 / 80.0;
    dkm(1) = 42.0 / 80.0;
    dkm(2) = 29.0 / 80.0;

    dkp = dkm;
    dkp.row(0).swap(dkp.row(1));

    dkmd(0) = 1.0 / 18.0;
    dkmd(1) = 5.0 / 6.0;
    dkmd(2) = 1.0 / 9.0;

    dkpd = dkmd;
    dkpd.row(0).swap(dkpd.row(1));

    crjlin(0, 0) = -23.0 / 120.0;
    crjlin(1, 0) = 19.0 / 30.0;
    crjlin(2, 0) = 67.0 / 120.0;
    crjlin(3, 0) = -3.0 / 40.0 * dif;
    crjlin(4, 0) = -7.0 / 40.0 * dif;

    crjlin(0,1) = 67.0 / 120.0;
    crjlin(1,1) = 19.0 / 30.0;
    crjlin(2,1) = -23.0 / 120.0;
    crjlin(3,1) = 7.0 / 40.0 * dif;
    crjlin(4,1) = 3.0 / 40.0 * dif;

    crjlind(0,0) = 0.25 / dif;
    crjlind(1,0) = -2.0 / dif;
    crjlind(2,0) = 1.75 / dif;
    crjlind(3,0) = 1.0 / 12.0;
    crjlind(4,0) = -1.0 / 6.0;
    crjlind(5,0) = -5.0 / 12.0;

    crjlind(0,1) = -7.0 / 4.0 / dif;
    crjlind(1,1) = 2.0 / dif;
    crjlind(2,1) = -0.25 / dif;
    crjlind(3,1) = -5.0 / 12.0;
    crjlind(4,1) = -1.0 / 6.0;
    crjlind(5,1) = 1.0 / 12.0;

    /* Two point quadrature */
    crjgp(0, 0) = -sq3 / 3.0;
    crjgp(1, 0) = 1.0 + sq3 / 3.0;
    crjgp(2, 0) = -dif * sq3 / 6.0;
    crjgp(0, 1) = 1.0 - sq3 / 3.0;
    crjgp(1, 1) = sq3 / 3.0;
    crjgp(2, 1) = -dif * sq3 / 6.0;
    crjgp(0, 2) = -sq3 / 12.0;
    crjgp(1, 2) = 1.0;
    crjgp(2, 2) = sq3 / 12.0;

    crjgm(0, 0) = sq3 / 3.0;
    crjgm(1, 0) = 1.0 - sq3 / 3.0;
    crjgm(2, 0) = dif * sq3 / 6.0;

    crjgm(0, 1) = 1.0 + sq3 / 3.0;
    crjgm(1, 1) = -sq3 / 3.0;
    crjgm(2, 1) = dif * sq3 / 6.0;

    crjgm(0, 2) = sq3 / 12.0;
    crjgm(1, 2) = 1.0;
    crjgm(2, 2) = -sq3 / 12.0;

    denom = 360.0;
    dkpg(0) = -sq3 + 105.0;
    dkpg(1) = sq3 + 105.0;
    dkpg = dkpg / denom;
    dkpg(2) = 5.0 / 12.0;

    dkmg(0) = sq3 + 105.0;
    dkmg(1) = -sq3 + 105.0;
    dkmg = dkmg / denom;
    dkmg(2) = 5.0 / 12.0;

    crjdgp(0, 0) = sq3 / dif;
    crjdgp(1, 0) = -sq3 / dif;
    crjdgp(2, 0) = sq3 / 3.0;
    crjdgp(3, 0) = (2.0 * sq3 + 3.0) / 3.0;
    crjdgp(0, 1) = -sq3 / dif;
    crjdgp(1, 1) = sq3 / dif;
    crjdgp(2, 1) = (-2.0 * sq3 + 3.0) / 3.0;
    crjdgp(3, 1) = -sq3 / 3.0;
    crjdgp(0, 2) = sq3 / 6.0 / dif;
    crjdgp(1, 2) = -sq3 / 3.0 / dif;
    crjdgp(2, 2) = sq3 / 6.0 / dif;
    crjdgp(3, 2) = 1.0;

    crjdgm(0, 0) = -sq3 / dif;
    crjdgm(1, 0) = sq3 / dif;
    crjdgm(2, 0) = -sq3 / 3.0;
    crjdgm(3, 0) = (-2.0 * sq3 + 3.0) / 3.0;
    crjdgm(0, 1) = sq3 / dif;
    crjdgm(1, 1) = -sq3 / dif;
    crjdgm(2, 1) = (2.0 * sq3 + 3.0) / 3.0;
    crjdgm(3, 1) = sq3 / 3.0;
    crjdgm(0, 2) = -sq3 / 6.0 / dif;
    crjdgm(1, 2) = sq3 / 3.0 / dif;
    crjdgm(2, 2) = -sq3 / 6.0 / dif;
    crjdgm(3, 2) = 1.0;

    dkpdg(0) = -sq3 / 144.0 + 1.0 / 3.0;
    dkpdg(1) = sq3 / 144.0 + 1.0 / 3.0;
    dkpdg(2) = 1.0 / 3.0;

    dkmdg(0) = sq3 / 144.0 + 1.0 / 3.0;
    dkmdg(1) = -sq3 / 144.0 + 1.0 / 3.0;
    dkmdg(2) = 1.0 / 3.0;

    /* Four point quadrature */
    crjg4p(0, 0) = -sq5 / 5.0 + 1.0 / 30.0;
    crjg4p(1, 0) = sq5 / 5.0 + 29.0 / 30.0;
    crjg4p(2, 0) = (-sq5 / 10.0 + 1.0 / 30.0) * dif;
    crjg4p(0, 1) = -sq5 / 5.0 + 29.0 / 30.0;
    crjg4p(1, 1) = sq5 / 5.0 + 1.0 / 30.0;
    crjg4p(2, 1) = -(sq5 / 10.0 + 1.0 / 30.0) * dif;
    crjg4p(0, 2) = -(sq5 / 20.0 + 1.0 / 60.0);
    crjg4p(1, 2) = 31.0 / 30.0;
    crjg4p(2, 2) = sq5 / 20.0 - 1.0 / 60.0;

    crjg4m(0, 0) = sq5 / 5.0 + 1.0 / 30.0;
    crjg4m(1, 0) = -sq5 / 5.0 + 29.0 / 30.0;
    crjg4m(2, 0) = (sq5 / 10.0 + 1.0 / 30.0) * dif;
    crjg4m(0, 1) = sq5 / 5.0 + 29.0 / 30.0;
    crjg4m(1, 1) = -sq5 / 5.0 + 1.0 / 30.0;
    crjg4m(2, 1) = (sq5 / 10.0 - 1.0 / 30.0) * dif;
    crjg4m(0, 2) = sq5 / 20.0 - 1.0 / 60.0;
    crjg4m(1, 2) = 31.0 / 30.0;
    crjg4m(2, 2) = -(sq5 / 20.0 + 1.0 / 60.0);

    denom = 330.0 * sq5 - 110.0;
    dkpg4(0) = 99.0 * sq5 + 33.0;
    dkpg4(1) = 108.0 * sq5 - 102.0;
    dkpg4(2) = 123.0 * sq5 - 41.0;
    dkpg4 /= denom;

    denom = 330.0 * sq5 + 110.0;
    dkmg4(0) = 99.0 * sq5 - 33.0;
    dkmg4(1) = 108.0 * sq5 + 102.0;
    dkmg4(2) = 123.0 * sq5 + 41.0;
    dkmg4 /= denom;

    crjdg4p(0, 0) = (6.0 * sq5 - 2.0) / 10.0 / dif;
    crjdg4p(1, 0) = (-6.0 * sq5 + 2.0) / 10.0 / dif;
    crjdg4p(2, 0) = 0.2 * sq5 - 0.1;
    crjdg4p(3, 0) = 0.4 * sq5 + 0.9;
    crjdg4p(0, 1) = -(6.0 * sq5 + 2.0) / 10.0 / dif;
    crjdg4p(1, 1) = (6.0 * sq5 + 2.0) / 10.0 / dif;
    crjdg4p(2, 1) = -0.4 * sq5 + 0.9;
    crjdg4p(3, 1) = -0.2 * sq5 - 0.1;
    crjdg4p(0, 2) = (2.0 * sq5 + 1.0) / 20.0 / dif;
    crjdg4p(1, 2) = -sq5 / 5.0 / dif;
    crjdg4p(2, 2) = (2.0 * sq5 - 1.0) / 20.0 / dif;
    crjdg4p(3, 2) = 1.1;

    crjdg4m(0, 0) = -(6.0 * sq5 + 2.0) / 10.0 / dif;
    crjdg4m(1, 0) = (6.0 * sq5 + 2.0) / 10.0 / dif;
    crjdg4m(2, 0) = -0.2 * sq5 - 0.1;
    crjdg4m(3, 0) = -0.4 * sq5 + 0.9;
    crjdg4m(0, 1) = (6.0 * sq5 - 2.0) / 10.0 / dif;
    crjdg4m(1, 1) = (-6.0 * sq5 + 2.0) / 10.0 / dif;
    crjdg4m(2, 1) = 0.4 * sq5 + 0.9;
    crjdg4m(3, 1) = 0.2 * sq5 - 0.1;
    crjdg4m(0, 2) = (-2.0 * sq5 + 1.0) / 20.0 / dif;
    crjdg4m(1, 2) = sq5 / 5.0 / dif;
    crjdg4m(2, 2) = (-2.0 * sq5 - 1.0) / 20.0 / dif;
    crjdg4m(3, 2) = 1.1;

    denom = 1140.0 * sq5 - 570.0;
    dkpdg4(0) = 399.0 * sq5 + 190.0;
    dkpdg4(1) = 481.0 * sq5 - 630.0;
    dkpdg4(2) = 260.0 * sq5 - 130.0;
    dkpdg4 /= denom;

    denom = 1140.0 * sq5 + 570.0;
    dkmdg4(0) = 399.0 * sq5 - 190.0;
    dkmdg4(1) = 481.0 * sq5 + 630.0;
    dkmdg4(2) = 260.0 * sq5 + 130.0;
    dkmdg4 /= denom;
}

Vector4d Schemes::Hermite::hweno(const Vector4d &u, const Vector4d &derivatives) {
    ud = derivatives * dif;

//     Reconstruction u^-_{i+1/2}
    p(0) = crj(0, 0) * u(0) + crj(1, 0) * u(1) +
           crj(2, 0) * derivatives(0);
    p(1) = crj(0, 1) * u(1) + crj(1, 1) * u(2) +
           crj(2, 1) * derivatives(2);
    p(2) = crj(0, 2) * u(0) + crj(1, 2) * u(1) +
           crj(2, 2) * u(2);

    smoothnessMeasure(u.segment(0, 3), ud.segment(0, 3));
    weights = dkm.array() / (SI.array() + eps).pow(2);
    weights = weights.array() / weights.sum();
    out(0) = p.dot(weights);
    out(0) = catchnan(out(0));

//     Reconstruction u^+_{i+1/2}
    p(0) = crj(1, 1) * u(1) + crj(0, 1) * u(2) -
           crj(2, 1) * derivatives(1);
    p(1) = crj(1, 0) * u(2) + crj(0, 0) * u(3) -
           crj(2, 0) * derivatives(3);
    p(2) = crj(2, 2) * u(1) + crj(1, 2) * u(2) + crj(0, 2) * u(3);

    smoothnessMeasure(u.segment(1, 3), ud.segment(1, 3));
    weights = dkp.array() / (SI.array() + eps).pow(2);
    weights = weights.array() / weights.sum();
    out(1) = p.dot(weights);
    out(1) = catchnan(out(1));

//     Reconstruction derivatives^-_{i+1/2}
    p(0) = crjd(0, 0) * u(0) + crjd(1, 0) * u(1) +
           crjd(2, 0) * derivatives(0) + crjd(3, 0) * derivatives(1);
    p(1) = crjd(0, 1) * u(1) + crjd(1, 1) * u(2) +
           crjd(2, 1) * derivatives(1) + crjd(3, 1) * derivatives(2);
    p(2) = crjd(0, 2) * u(0) + crjd(1, 2) * u(1) +
           crjd(2, 2) * u(2) + crjd(3, 2) * derivatives(1);

    smoothnessMeasureDerivative(u.segment(0, 3), ud.segment(0, 3));
    weights = dkmd.array() / (SI.array() + eps).pow(2);
    weights = weights.array() / weights.sum();
    out(2) = p.dot(weights);
    out(2) = catchnan(out(2));

//     Reconstruction derivatives^+_{i+1/2}
    p(0) = crjd(0, 1) * u(1) + crjd(1, 1) * u(2) +
           crjd(3, 1) * derivatives(1) + crjd(2, 1) * derivatives(2);
    p(1) = crjd(0, 0) * u(2) + crjd(1, 0) * u(3) +
           crjd(3, 0) * derivatives(2) + crjd(2, 0) * derivatives(3);
    p(2) = -(crjd(2, 2) * u(1) + crjd(1, 2) * u(2) + crjd(0, 2) * u(3)) +
           crjd(3, 2) * derivatives(2);

    smoothnessMeasureDerivative(u.segment(1, 3), ud.segment(1, 3));
    weights = dkpd.array() / (SI.array() + eps).pow(2);
    weights = weights.array() / weights.sum();
    out(3) = p.dot(weights);
    out(3) = catchnan(out(3));

    return out;
}

Vector4d Schemes::Hermite::hteno(const Vector4d &u, const Vector4d &derivatives) {
    ud = derivatives * dif;
    double C = 1.0;
    double Cd = 1.0;
    int power = 6;
    int powerD = 6;
    Vector3d undevidedDiff, undevidedDiffDer;//, powerq;

    /* Reconstruction u^-_{i+1/2} */
    p(0) = crj(0, 0) * u(0) + crj(1, 0) * u(1) +
           crj(2, 0) * derivatives(0);
    p(1) = crj(0, 1) * u(1) + crj(1, 1) * u(2) +
           crj(2, 1) * derivatives(2);
    p(2) = crj(0, 2) * u(0) + crj(1, 2) * u(1) +
           crj(2, 2) * u(2);

    smoothnessMeasure(u.segment(0, 3), ud.segment(0, 3));
    fullPointsSI(u.segment(0, 3), ud.segment(0, 3));
    tau = fabs(brk - 0.25 * (SI(0) + SI(1) + 2.0 * SI(2)));
    undevidedDiff = tau / (SI.array() + epsz);
    weights = (C + undevidedDiff.array()).pow(power);
    cutOff();
    weights = sigm.cwiseProduct(dkm);
    weights = weights.array() / weights.sum();
    out(0) = p.dot(weights);
    out(0) = catchnan(out(0));

    /* Reconstruction derivatives^-_{i+1/2} */
    p(0) = crjd(0, 0) * u(0) + crjd(1, 0) * u(1) +
           crjd(2, 0) * derivatives(0) + crjd(3, 0) * derivatives(1);
    p(1) = crjd(0, 1) * u(1) + crjd(1, 1) * u(2) +
           crjd(2, 1) * derivatives(1) + crjd(3, 1) * derivatives(2);
    p(2) = crjd(0, 2) * u(0) + crjd(1, 2) * u(1) +
           crjd(2, 2) * u(2) + crjd(3, 2) * derivatives(1);

    smoothnessMeasureDerivative(u.segment(0, 3), ud.segment(0, 3));
    fullPointsSId(u.segment(0, 3), ud.segment(0, 3));
    tau = fabs(brkd - 0.25 * (SI(0) + SI(1) + 2.0 * SI(2)));
    undevidedDiffDer = tau / (SI.array() + epszd);
    weights = (Cd + undevidedDiffDer.array()).pow(powerD);
    cutOffd();
    weights = sigm.cwiseProduct(dkmd);
    weights = weights.array() / weights.sum();
    out(2) = p.dot(weights);
    out(2) = catchnan(out(2));

    /* Reconstruction u^+_{i+1/2} */
    p(0) = crj(1, 1) * u(1) + crj(0, 1) * u(2) -
           crj(2, 1) * derivatives(1);
    p(1) = crj(1, 0) * u(2) + crj(0, 0) * u(3) -
           crj(2, 0) * derivatives(3);
    p(2) = crj(2, 2) * u(1) + crj(1, 2) * u(2) + crj(0, 2) * u(3);

    smoothnessMeasure(u.segment(1, 3), ud.segment(1, 3));
    fullPointsSI(u.segment(1, 3), ud.segment(1, 3));
    tau = fabs(brk - 0.25 * (SI(0) + SI(1) + 2.0 * SI(2)));
    undevidedDiff = tau / (SI.array() + epsz);
    weights = (C + undevidedDiff.array()).pow(power);
    cutOff();
    weights = sigm.cwiseProduct(dkp);
    weights = weights.array() / weights.sum();
    out(1) = p.dot(weights);
    out(1) = catchnan(out(1));

    /* Reconstruction derivatives^+_{i+1/2} */
    p(0) = crjd(0, 1) * u(1) + crjd(1, 1) * u(2) +
           crjd(3, 1) * derivatives(1) + crjd(2, 1) * derivatives(2);
    p(1) = crjd(0, 0) * u(2) + crjd(1, 0) * u(3) +
           crjd(3, 0) * derivatives(2) + crjd(2, 0) * derivatives(3);
    p(2) = -(crjd(2, 2) * u(1) + crjd(1, 2) * u(2) + crjd(0, 2) * u(3)) +
           crjd(3, 2) * derivatives(2);

    smoothnessMeasureDerivative(u.segment(1, 3), ud.segment(1, 3));
    fullPointsSId(u.segment(1, 3), ud.segment(1, 3));
    tau = fabs(brkd - 0.25 * (SI(0) + SI(1) + 2.0 * SI(2)));
    undevidedDiffDer = tau / (SI.array() + epszd);
    weights = (Cd + undevidedDiffDer.array()).pow(powerD);
    cutOffd();
    weights = sigm.cwiseProduct(dkpd);
    weights = weights.array() / weights.sum();
    out(3) = p.dot(weights);
    out(3) = catchnan(out(3));

    return out;
}

Vector4d Schemes::Hermite::hweno2PointsGLQ(const Vector3d &u, const Vector3d &derivatives) {
    udglq = derivatives * dif;
    smoothnessMeasure2PointsGLQ(u, udglq);

    p(0) = crjgp(0, 0) * u(0) + crjgp(1, 0) * u(1) +
           crjgp(2, 0) * derivatives(0);
    p(1) = crjgp(0, 1) * u(1) + crjgp(1, 1) * u(2) +
           crjgp(2, 1) * derivatives(2);
    p(2) = crjgp(0, 2) * u(0) + crjgp(1, 2) * u(1) + crjgp(2, 2) * u(2);
    weights = dkpg.array() / (SIGLQ.array() + eps).pow(2);
    weights = weights.array() / weights.sum();
    out(0) = p.dot(weights);

    p(0) = crjgm(0, 0) * u(0) + crjgm(1, 0) * u(1) +
           crjgm(2, 0) * derivatives(0);
    p(1) = crjgm(0, 1) * u(1) + crjgm(1, 1) * u(2) +
           crjgm(2, 1) * derivatives(2);
    p(2) = crjgm(0, 2) * u(0) + crjgm(1, 2) * u(1) + crjgm(2, 2) * u(2);
    weights = dkmg.array() / (SIGLQ.array() + eps).pow(2);
    weights = weights.array() / weights.sum();
    out(1) = p.dot(weights);

    smoothnessMeasureDerivatives2PointsGLQ(u, udglq);
    p(0) = crjdgp(0, 0) * u(0) + crjdgp(1, 0) * u(1) +
           crjdgp(2, 0) * derivatives(0) + crjdgp(3, 0) * derivatives(1);
    p(1) = crjdgp(0, 1) * u(1) + crjdgp(1, 1) * u(2) +
           crjdgp(2, 1) * derivatives(1) + crjdgp(3, 1) * derivatives(2);
    p(2) = crjdgp(0, 2) * u(0) + crjdgp(1, 2) * u(1) +
           crjdgp(2, 2) * u(2) + crjdgp(3, 2) * derivatives(1);
    weights = dkpdg.array() / (SIGLQ.array() + eps).pow(2);
    weights = weights.array() / weights.sum();
    out(2) = p.dot(weights);

    p(0) = crjdgm(0, 0) * u(0) + crjdgm(1, 0) * u(1) +
           crjdgm(2, 0) * derivatives(0) + crjdgm(3, 0) * derivatives(1);
    p(1) = crjdgm(0, 1) * u(1) + crjdgm(1, 1) * u(2) +
           crjdgm(2, 1) * derivatives(1) + crjdgm(3, 1) * derivatives(2);
    p(2) = crjdgm(0, 2) * u(0) + crjdgm(1, 2) * u(1) +
           crjdgm(2, 2) * u(2) + crjdgm(3, 2) * derivatives(1);
    weights = dkmdg.array() / (SIGLQ.array() + eps).pow(2);
    weights = weights.array() / weights.sum();
    out(3) = p.dot(weights);

    return out;
}

Matrix<double, 8, 1> Schemes::Hermite::hweno4PointsGLQ(const Vector3d &u, const Vector3d &derivatives) {
    udglq = derivatives * dif;

    // i+-1/2
    smoothnessMeasure4PointsGLQ(u, udglq);
    p(0) = crj(0, 0) * u(0) + crj(1, 0) * u(1) +
           crj(2, 0) * derivatives(0);
    p(1) = crj(0, 1) * u(1) + crj(1, 1) * u(2) +
           crj(2, 1) * derivatives(2);
    p(2) = crj(0, 2) * u(0) + crj(1, 2) * u(1) + crj(2, 2) * u(2);
    weights = dkm.array() / (SIGLQ.array() + eps).pow(2);
    weights = weights.array() / weights.sum();
    outg4(0) = p.dot(weights);

    p(0) = crj(1, 1) * u(0) + crj(0, 1) * u(1) -
           crj(2, 1) * derivatives(0);
    p(1) = crj(1, 0) * u(1) + crj(0, 0) * u(2) -
           crj(2, 0) * derivatives(2);
    p(2) = crj(2, 2) * u(0) + crj(1, 2) * u(1) + crj(0, 2) * u(2);
    weights = dkp.array() / (SIGLQ.array() + eps).pow(2);
    weights = weights.array() / weights.sum();
    outg4(1) = p.dot(weights);

    // i+-1/sqrt(5)/2
    p(0) = crjg4p(0, 0) * u(0) + crjg4p(1, 0) * u(1) +
           crjg4p(2, 0) * derivatives(0);
    p(1) = crjg4p(0, 1) * u(1) + crjg4p(1, 1) * u(2) +
           crjg4p(2, 1) * derivatives(2);
    p(2) = crjg4p(0, 2) * u(0) + crjg4p(1, 2) * u(1) + crjg4p(2, 2) * u(2);
    weights = dkpg4.array() / (SIGLQ4.array() + eps).pow(2);
    weights = weights.array() / weights.sum();
    outg4(2) = p.dot(weights);

    p(0) = crjg4m(0, 0) * u(0) + crjg4m(1, 0) * u(1) +
           crjg4m(2, 0) * derivatives(0);
    p(1) = crjg4m(0, 1) * u(1) + crjg4m(1, 1) * u(2) +
           crjg4m(2, 1) * derivatives(2);
    p(2) = crjg4m(0, 2) * u(0) + crjg4m(1, 2) * u(1) + crjg4m(2, 2) * u(2);
    weights = dkmg4.array() / (SIGLQ4.array() + eps).pow(2);
    weights = weights.array() / weights.sum();
    outg4(3) = p.dot(weights);

    /* Reconstruction derivatives^-_{i+-1/2} */
    smoothnessMeasureDerivatives4PointsGLQ(u, udglq);
    p(0) = crjd(0, 0) * u(0) + crjd(1, 0) * u(1) +
           crjd(2, 0) * derivatives(0) + crjd(3, 0) * derivatives(1);
    p(1) = crjd(0, 1) * u(1) + crjd(1, 1) * u(2) +
           crjd(2, 1) * derivatives(1) + crjd(3, 1) * derivatives(2);
    p(2) = crjd(0, 2) * u(0) + crjd(1, 2) * u(1) +
           crjd(2, 2) * u(2) + crjd(3, 2) * derivatives(1);
    weights = dkmd.array() / (SIGLQ.array() + eps).pow(2);
    weights = weights.array() / weights.sum();
    outg4(4) = p.dot(weights);

    p(0) = crjd(0, 1) * u(0) + crjd(1, 1) * u(1) +
           crjd(3, 1) * derivatives(0) + crjd(2, 1) * derivatives(1);
    p(1) = crjd(0, 0) * u(1) + crjd(1, 0) * u(2) +
           crjd(3, 0) * derivatives(1) + crjd(2, 0) * derivatives(2);
    p(2) = -(crjd(2, 2) * u(0) + crjd(1, 2) * u(1) + crjd(0, 2) * u(2)) +
           crjd(3, 2) * derivatives(1);
    weights = dkpd.array() / (SIGLQ.array() + eps).pow(2);
    weights = weights.array() / weights.sum();
    outg4(5) = p.dot(weights);

    // i+-1/sqrt(5)/2
    p(0) = crjdg4p(0, 0) * u(0) + crjdg4p(1, 0) * u(1) +
           crjdg4p(2, 0) * derivatives(0) + crjdg4p(3, 0) * derivatives(1);
    p(1) = crjdg4p(0, 1) * u(1) + crjdg4p(1, 1) * u(2) +
           crjdg4p(2, 1) * derivatives(1) + crjdg4p(3, 1) * derivatives(2);
    p(2) = crjdg4p(0, 2) * u(0) + crjdg4p(1, 2) * u(1) +
           crjdg4p(2, 2) * u(2) + crjdg4p(3, 2) * derivatives(1);
    weights = dkpdg4.array() / (SIGLQ4.array() + eps).pow(2);
    weights = weights.array() / weights.sum();
    outg4(6) = p.dot(weights);

    p(0) = crjdg4m(0, 0) * u(0) + crjdg4m(1, 0) * u(1) +
           crjdg4m(2, 0) * derivatives(0) + crjdg4m(3, 0) * derivatives(1);
    p(1) = crjdg4m(0, 1) * u(1) + crjdg4m(1, 1) * u(2) +
           crjdg4m(2, 1) * derivatives(1) + crjdg4m(3, 1) * derivatives(2);
    p(2) = crjdg4m(0, 2) * u(0) + crjdg4m(1, 2) * u(1) +
           crjdg4m(2, 2) * u(2) + crjdg4m(3, 2) * derivatives(1);
    weights = dkmdg4.array() / (SIGLQ4.array() + eps).pow(2);
    weights = weights.array() / weights.sum();
    outg4(7) = p.dot(weights);

    return outg4;
}

void Schemes::Hermite::smoothnessMeasure(const Vector3d &u, const Vector3d &v) {
    b(0) = -2.0 * (u(0) - u(1)) - v(0);
    c(0) = -2.0 * (u(0) - u(1) + v(0));

    b(1) = -2.0 * (u(1) - u(2)) - v(2);
    c(1) = 2.0 * (u(1) - u(2)) + 2.0 * v(2);

    b(2) = (-u(0) + u(2)) / 2.0;
    c(2) = u(0) - 2.0 * u(1) + u(2);

    SI = c1 * c.array().pow(2) + b.array().pow(2);
}

void Schemes::Hermite::smoothnessMeasureDerivative(const Vector3d &u, const Vector3d &udf) {
    b(0) = 6.0 * (u(0) - u(1)) + (2.0 * udf(0) + 4.0 * udf(1));
    c(0) = 12.0 * (u(0) - u(1)) + 6.0 * (udf(0) + udf(1));

    b(1) = 6.0 * (u(2) - u(1)) - (4.0 * udf(1) + 2.0 * udf(2));
    c(1) = 12.0 * (u(1) - u(2)) + 6.0 * (udf(1) + udf(2));

    b(2) = u(0) - 2.0 * u(1) + u(2);
    c(2) = 3.0 * (u(2) - u(0)) - 6.0 * udf(1);

//    c=c/2.0;

    SI = c1 * c.array().pow(2) + b.array().pow(2);
}

void Schemes::Hermite::smoothnessMeasure2PointsGLQ(const Vector3d &u, const Vector3d &udf) {
    b(0) = -2.0 * (u(0) - u(1)) - udf(0);
    c(0) = -2.0 * (u(0) - u(1)) - 2.0 * udf(0);

    b(1) = -2.0 * (u(1) - u(2)) - udf(2);
    c(1) = 2.0 * (u(1) - u(2)) + 2.0 * udf(2);

    b(2) = -(u(0) - u(2)) / 2.0;
    c(2) = u(0) - 2.0 * u(1) + u(2);

    SIGLQ = c2 * b.array().pow(2) + c3 * c.array().pow(2);
}

void Schemes::Hermite::smoothnessMeasureDerivatives2PointsGLQ(const Vector3d &u, const Vector3d &udf) {
    b(0) = 6.0 * (u(0) - u(1)) + (2.0 * udf(0) + 4.0 * udf(1));
    c(0) = 12.0 * (u(0) - u(1)) + 6.0 * (udf(0) + udf(1));

    b(1) = 6.0 * (u(2) - u(1)) - (4.0 * udf(1) + 2.0 * udf(2));
    c(1) = 12.0 * (u(1) - u(2)) + 6.0 * (udf(1) + udf(2));

    b(2) = u(0) - 2.0 * u(1) + u(2);
    c(2) = 3.0 * (u(2) - u(0)) - 6.0 * udf(1);

    SIGLQ = c4 * b.array().pow(2) + c5 * c.array().pow(2);
}

void Schemes::Hermite::smoothnessMeasure4PointsGLQ(const Vector3d &u, const Vector3d &udf) {
    b(0) = -2.0 * (u(0) - u(1)) - udf(0);
    c(0) = -2.0 * (u(0) - u(1)) - 2.0 * udf(0);

    b(1) = -2.0 * (u(1) - u(2)) - udf(2);
    c(1) = 2.0 * (u(1) - u(2)) + 2.0 * udf(2);

    b(2) = -(u(0) - u(2)) / 2.0;
    c(2) = u(0) - 2.0 * u(1) + u(2);

    SIGLQ = c1 * c.array().pow(2) + b.array().pow(2);
    SIGLQ4 = c6 * b.array().pow(2) + c7 * c.array().pow(2);
}

void Schemes::Hermite::smoothnessMeasureDerivatives4PointsGLQ(const Vector3d &u, const Vector3d &udf) {
    b(0) = 6.0 * (u(0) - u(1)) + (2.0 * udf(0) + 4.0 * udf(1));
    c(0) = 12.0 * (u(0) - u(1)) + 6.0 * (udf(0) + udf(1));

    b(1) = 6.0 * (u(2) - u(1)) - (4.0 * udf(1) + 2.0 * udf(2));
    c(1) = 12.0 * (u(1) - u(2)) + 6.0 * (udf(1) + udf(2));

    b(2) = u(0) - 2.0 * u(1) + u(2);
    c(2) = 3.0 * (u(2) - u(0)) - 6.0 * udf(1);

    SIGLQ = c1 * c.array().pow(2) + b.array().pow(2);
    SIGLQ4 = c6 * b.array().pow(2) + c7 * c.array().pow(2);
}

void Schemes::Hermite::fullPointsSI(const Vector3d &u, const Vector3d &udf) {
    brk = fabs(
            (1029270.0 * u(2) - 2129184.0 * u(1) + 1099914.0 * u(0) + 281121.0 * udf(0) - 494958.0 * udf(2)) * udf(0) +
            (-1029270.0 * u(0) + 2129184.0 * u(1) - 1099914.0 * u(2) + 281121.0 * udf(2)) * udf(2) +
            (1099445.0 * u(2) - 4317056.0 * u(1) + 2118166.0 * u(0)) * u(2) +
            4317056.0 * (u(1) - u(0)) * u(1) + 1099445.0 * u(0) * u(0)) / 6720.0;
}

void Schemes::Hermite::fullPointsSId(const Vector3d &u, const Vector3d &udf) {
    brkd = fabs((3156862.0 * udf(2) + 13189240.0 * udf(1) + 1710479.0 * udf(0) - 9617754.0 * u(2) - 531552.0 * u(1) +
                 10149306.0 * u(0)) * udf(0) +
                (13189240.0 * udf(2) + 26511200.0 * udf(1) - 39700440.0 * u(2) + 39700440.0 * u(0)) * udf(1) +
                (1710479.0 * udf(2) - 10149306.0 * u(2) + 531552.0 * u(1) + 9617754.0 * u(0)) * udf(2) +
                (15136011.0 * u(2) - 1076544.0 * u(1) - 29195478.0 * u(0)) * u(2) +
                1076544.0 * (u(1) - u(0)) * u(1) + 15136011.0 * u(0) * u(0)) / 1680.0;
}

void Schemes::Hermite::cutOff() {
    int i;
    weights /= weights.sum();
    for (i = 0; i < 3; i++) {
        sigm(i) = weights(i) < ct ? 0.0 : 1.0;
    }
}

void Schemes::Hermite::cutOffd() {
    int i;
    weights /= weights.sum();
    for (i = 0; i < 3; i++) {
        sigm(i) = weights(i) < ctd ? 0.0 : 1.0;
    }
}

double Schemes::Hermite::catchnan(const double &input) {
    return isnan(input) ? machine_eps : input;
}

double Schemes::Hermite::hermiteLimiter(const Vector3d &u, const Vector3d &derivatives) {
    p(0) = 0.75 * (u(2) - u(0)) / dif -
           0.25 * (derivatives(0) + derivatives(2));
    p(1) = (u(1) - u(0)) / dif;
    p(2) = (u(2) - u(1)) / dif;

    Vector3d linWeights, smoothness;
//    linWeights.setRandom();
//    linWeights = linWeights.cwiseAbs();
//    linWeights = linWeights / linWeights.sum();
//      Default hteno
    linWeights(0) = 1.0 - 2E-5, linWeights(1) = 1E-5, linWeights(2) = 1E-5;

//  Hybrid htenoi
//    linWeights(0) = 1.0 - 2E-5, linWeights(1) = 1E-5, linWeights(2) = 1E-5;

    fullPointsSI(u, derivatives * dif);
    smoothness(0) = brk;
    smoothness(1) = pow(u(0) - u(1), 2);
    smoothness(2) = pow(u(1) - u(2), 2);

    tau = pow(fabs(smoothness(0) - smoothness(1)) +
              fabs(smoothness(0) - smoothness(2)), 2) / 4.0;
    weights = linWeights.array() * (1.0 + tau / (smoothness.array() + eps));
    weights = weights / weights.sum();

    return weights(0) * (p(0) / linWeights(0) - (linWeights.segment(1, 2) / linWeights(0)).dot(p.segment(1, 2))) +
           weights.segment(1, 2).dot(p.segment(1, 2));
}

bool Schemes::Hermite::KXRCF(const double &uj, const double &unbj, const double &unorm) const {
    double Ikxrcf;
    Ikxrcf = fabs(uj - unbj) / (fabs(unorm) * pow(dif, 1.5));
    return Ikxrcf > 1.0;
}

Vector4d Schemes::Hermite::htenoi(const Vector4d &u, const Vector4d &derivatives) {
    ud = derivatives * dif;
    double C = 1.0;
    double Cd = 1.0;
    int power = 6;
    int powerD = 6;
    Vector3d undevidedDiff, undevidedDiffDer;//, powerq;

    /* Reconstruction u^-_{i+1/2} */
    p(0) = crj(0, 0) * u(0) + crj(1, 0) * u(1) +
           crj(2, 0) * derivatives(0);
    p(1) = crj(0, 1) * u(1) + crj(1, 1) * u(2) +
           crj(2, 1) * derivatives(2);
    p(2) = crj(0, 2) * u(0) + crj(1, 2) * u(1) +
           crj(2, 2) * u(2);

    smoothnessMeasure(u.segment(0, 3), ud.segment(0, 3));
    //fullPointsSI(u.segment(0, 3), ud.segment(0, 3));
    //tau = fabs(brk - 0.25 * (SI(0) + SI(1) + 2.0 * SI(2)));
    tau = SI.maxCoeff();
    undevidedDiff = tau / (SI.array() + epsz);
    weights = (C + undevidedDiff.array()).pow(power);
    cutOff();
    weights = sigm.cwiseProduct(dkm);
    weights = weights.array() / weights.sum();
    out(0) = p.dot(weights);
    out(0) = catchnan(out(0));

    /* Reconstruction derivatives^-_{i+1/2} */
    p(0) = crjd(0, 0) * u(0) + crjd(1, 0) * u(1) +
           crjd(2, 0) * derivatives(0) + crjd(3, 0) * derivatives(1);
    p(1) = crjd(0, 1) * u(1) + crjd(1, 1) * u(2) +
           crjd(2, 1) * derivatives(1) + crjd(3, 1) * derivatives(2);
    p(2) = crjd(0, 2) * u(0) + crjd(1, 2) * u(1) +
           crjd(2, 2) * u(2) + crjd(3, 2) * derivatives(1);

    smoothnessMeasureDerivative(u.segment(0, 3), ud.segment(0, 3));
    //fullPointsSId(u.segment(0, 3), ud.segment(0, 3));
    //tau = fabs(brkd - 0.25 * (SI(0) + SI(1) + 2.0 * SI(2)));
    tau = SI.maxCoeff();
    undevidedDiffDer = tau / (SI.array() + epszd);
    weights = (Cd + undevidedDiffDer.array()).pow(powerD);
    cutOffd();
    weights = sigm.cwiseProduct(dkmd);
    weights = weights.array() / weights.sum();
    out(2) = p.dot(weights);
    out(2) = catchnan(out(2));

    /* Reconstruction u^+_{i+1/2} */
    p(0) = crj(1, 1) * u(1) + crj(0, 1) * u(2) -
           crj(2, 1) * derivatives(1);
    p(1) = crj(1, 0) * u(2) + crj(0, 0) * u(3) -
           crj(2, 0) * derivatives(3);
    p(2) = crj(2, 2) * u(1) + crj(1, 2) * u(2) + crj(0, 2) * u(3);

    smoothnessMeasure(u.segment(1, 3), ud.segment(1, 3));
//    fullPointsSI(u.segment(1, 3), ud.segment(1, 3));
//    tau = fabs(brk - 0.25 * (SI(0) + SI(1) + 2.0 * SI(2)));
    tau = SI.maxCoeff();
    undevidedDiff = tau / (SI.array() + epsz);
    weights = (C + undevidedDiff.array()).pow(power);
    cutOff();
    weights = sigm.cwiseProduct(dkp);
    weights = weights.array() / weights.sum();
    out(1) = p.dot(weights);
    out(1) = catchnan(out(1));

    /* Reconstruction derivatives^+_{i+1/2} */
    p(0) = crjd(0, 1) * u(1) + crjd(1, 1) * u(2) +
           crjd(3, 1) * derivatives(1) + crjd(2, 1) * derivatives(2);
    p(1) = crjd(0, 0) * u(2) + crjd(1, 0) * u(3) +
           crjd(3, 0) * derivatives(2) + crjd(2, 0) * derivatives(3);
    p(2) = -(crjd(2, 2) * u(1) + crjd(1, 2) * u(2) + crjd(0, 2) * u(3)) +
           crjd(3, 2) * derivatives(2);

    smoothnessMeasureDerivative(u.segment(1, 3), ud.segment(1, 3));
//    fullPointsSId(u.segment(1, 3), ud.segment(1, 3));
//    tau = fabs(brkd - 0.25 * (SI(0) + SI(1) + 2.0 * SI(2)));
    tau = SI.maxCoeff();
    undevidedDiffDer = tau / (SI.array() + epszd);
    weights = (Cd + undevidedDiffDer.array()).pow(powerD);
    cutOffd();
    weights = sigm.cwiseProduct(dkpd);
    weights = weights.array() / weights.sum();
    out(3) = p.dot(weights);
    out(3) = catchnan(out(3));

    return out;
}

Vector4d Schemes::Hermite::linear(const Vector4d &u, const Vector4d &derivatives) {
    out(0) = crjlin.col(0).segment(0,3).dot(u.segment(0,3)) +
            crjlin(3,0) * derivatives(0) + crjlin(4,0) * derivatives(2);
    out(1) = crjlin.col(1).segment(0,3).dot(u.segment(1,3)) +
            crjlin(3,1) * derivatives(1) + crjlin(4,1) * derivatives(3);
    out(2) = crjlind.col(0).segment(0,3).dot(u.segment(0,3)) +
            crjlind.col(0).segment(3,3).dot(derivatives.segment(0,3));
    out(3) = crjlind.col(1).segment(0,3).dot(u.segment(1,3)) +
            crjlind.col(1).segment(3,3).dot(derivatives.segment(1,3));

//    ud = derivatives * dif;
//
//    /* Reconstruction u^-_{i+1/2} */
//    p(0) = crj(0, 0) * u(0) + crj(1, 0) * u(1) +
//            crj(2, 0) * derivatives(0);
//    p(1) = crj(0, 1) * u(1) + crj(1, 1) * u(2) +
//            crj(2, 1) * derivatives(2);
//    p(2) = crj(0, 2) * u(0) + crj(1, 2) * u(1) +
//            crj(2, 2) * u(2);
//
//    out(0) = p.dot(dkm);
//
//    /* Reconstruction derivatives^-_{i+1/2} */
//    p(0) = crjd(0, 0) * u(0) + crjd(1, 0) * u(1) +
//            crjd(2, 0) * derivatives(0) + crjd(3, 0) * derivatives(1);
//    p(1) = crjd(0, 1) * u(1) + crjd(1, 1) * u(2) +
//            crjd(2, 1) * derivatives(1) + crjd(3, 1) * derivatives(2);
//    p(2) = crjd(0, 2) * u(0) + crjd(1, 2) * u(1) +
//            crjd(2, 2) * u(2) + crjd(3, 2) * derivatives(1);
//
//    out(2) = p.dot(dkmd);
//
//    /* Reconstruction u^+_{i+1/2} */
//    p(0) = crj(1, 1) * u(1) + crj(0, 1) * u(2) -
//            crj(2, 1) * derivatives(1);
//    p(1) = crj(1, 0) * u(2) + crj(0, 0) * u(3) -
//            crj(2, 0) * derivatives(3);
//    p(2) = crj(2, 2) * u(1) + crj(1, 2) * u(2) + crj(0, 2) * u(3);
//
//    out(1) = p.dot(dkp);
//
//    /* Reconstruction derivatives^+_{i+1/2} */
//    p(0) = crjd(0, 1) * u(1) + crjd(1, 1) * u(2) +
//            crjd(3, 1) * derivatives(1) + crjd(2, 1) * derivatives(2);
//    p(1) = crjd(0, 0) * u(2) + crjd(1, 0) * u(3) +
//            crjd(3, 0) * derivatives(2) + crjd(2, 0) * derivatives(3);
//    p(2) = -(crjd(2, 2) * u(1) + crjd(1, 2) * u(2) + crjd(0, 2) * u(3)) +
//            crjd(3, 2) * derivatives(2);
//
//    out(3) = p.dot(dkpd);

    return out;
}
