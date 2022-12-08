#pragma once
#include <cmath>
#include <Eigen/Dense>
#include <iostream>

using Eigen::Matrix;
using Eigen::Dynamic;
using std::cout;

using Eigen::Vector2d;
using Eigen::Vector3d;
using Eigen::Vector4d;
typedef Matrix<double, 5, 1> Vector5d;
typedef Matrix<double, 6, 1> Vector6d;
typedef Matrix<double, 9, 1> Vector9d;
typedef Matrix<double, 8, 1> Vector8d;
typedef Matrix<double, 12, 1> Vector12d;


namespace Schemes {
	class Hermite {
	public:
		Hermite();
		void setGridSize(const double& gridSize);
		Vector4d hweno(const Vector4d& u, const Vector4d& derivatives);
		Vector4d hteno(const Vector4d& u, const Vector4d& derivatives);
		Vector4d linear(const Vector4d& u, const Vector4d& derivatives);

        Vector4d htenoi(const Vector4d& u, const Vector4d& derivatives);

		double hermiteLimiter(const Vector3d& u, const Vector3d& derivatives);
		Vector4d hweno2PointsGLQ(const Vector3d& u, const Vector3d& derivatives);
		Matrix<double, 8, 1> hweno4PointsGLQ(const Vector3d& u, const Vector3d& derivatives);
		bool KXRCF(const double& u0, const double& ul, const double& unorm) const;

	private:
		void smoothnessMeasure(const Vector3d& u, const Vector3d& v);
		void smoothnessMeasureDerivative(const Vector3d& u, const Vector3d& ud);

		void smoothnessMeasure2PointsGLQ(const Vector3d& u, const Vector3d& ud);
		void smoothnessMeasureDerivatives2PointsGLQ(const Vector3d& u, const Vector3d& ud);

		void smoothnessMeasure4PointsGLQ(const Vector3d& u, const Vector3d& ud);
		void smoothnessMeasureDerivatives4PointsGLQ(const Vector3d& u, const Vector3d& ud);

		void fullPointsSI(const Vector3d& u, const Vector3d& ud);
		void fullPointsSId(const Vector3d& u, const Vector3d& ud);

		void cutOff();
		void cutOffd();

		double catchnan(const double& input);

		Matrix<double, 3, 3> crj, crjgp, crjgm, crjg4p, crjg4m;
		Matrix<double, 4, 3> crjd, crjdgp, crjdgm, crjdg4p, crjdg4m;

		Matrix<double, 5, 2> crjlin;
		Matrix<double, 6, 2> crjlind;

		Vector3d dkm, dkp, dkmd, dkpd, dkmg, dkpg, dkmdg, dkpdg, dkmg4, dkpg4, dkmdg4, dkpdg4;
		Vector3d udglq, b, c, p, weights;
		Vector3d SI, SIGLQ, SIGLQ4, sigm;
		Vector4d ud, out;
		Vector8d outg4;

        const double eps, machine_eps;
        const double epsz, epszd, ct, ctd;
        double dif;
		double c1, c2, c3, c4, c5, c6, c7;
		double sq3, sq5;
		double brk, brkd, tau;
	};
}