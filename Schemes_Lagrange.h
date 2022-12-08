#pragma once
#include <cmath>
#include <Eigen/Dense>
#include <iostream>
#include <limits>

using Eigen::Matrix;
using Eigen::Dynamic;

using Eigen::Vector2d;
using Eigen::Vector3d;
using Eigen::Vector4d;
typedef Matrix<double, 5, 1> Vector5d;
typedef Matrix<double, 6, 1> Vector6d;
typedef Matrix<double, 9, 1> Vector9d;

using Eigen::Matrix2d;
using Eigen::Matrix4d;

namespace Schemes {
	class Lagrange
	{
	public:
		Lagrange();
		Vector2d weno3JS(const Vector4d& u);
		Vector2d weno3JS2PointsQuadrature(const Vector3d& u);

		Vector2d weno5JS(const Vector6d& u);
		Vector2d linear(const Vector6d& u);

		Vector2d weno5Z(const Vector6d& u);
		Vector2d weno5Z(const Vector6d& u, const double& diff);

		Vector2d weno5N(const Vector6d& u);

		Vector2d teno5(const Vector6d& u);
		Vector2d teno5a(const Vector6d& u);
		Vector2d teno5i(const Vector6d& u);
		Vector2d teno5(const Vector6d& u, const double& diff);

		Vector2d teno5N(const Vector6d& u);

		Vector2d weno5IS(const Vector6d& u);
		Vector2d weno5ISmod(const Vector6d& u, const double& diff);
		Vector2d weno5IS_Extended(const Vector6d& u);

		Vector2d weno5JS2PointsQuadrature(const Vector5d& u);
		Vector2d weno5Z2PointsQuadrature(const Vector5d& u);
		Vector2d teno5th2PointsQuadrature(const Vector5d& u);

        Vector3d weno5JS3PointsQuadrature(const Vector5d& u);

		Vector4d weno5JS4PointsGaussLobatto(const Vector5d& u);

		Vector4d weno5Z4PointsGaussLobatto(const Vector5d& u);
		Vector4d teno5th4PointsGaussLobatto(const Vector5d& u);
		
		Vector4d weno5JS4PointsGaussLegendre(const Vector5d& u);
		Vector4d weno5Z4PointsGaussLegendre(const Vector5d& u);

		Vector9d getQuadrature9Weights(), getQuadrature9Points();

		inline Vector2d getQuadrature2Weights() {
			return gaussianQuadrature2Weights;
		}

        inline Vector3d getGaussLegendre3Weights() {
            return gaussLegendreQuadrature3Weights;
        }

		inline Vector4d getGaussLobatto4Weights() {
			return gaussLobattoQuadrature4Weights;
		}

		inline Vector4d getGaussLegendre4Weights() {
			return gaussLegendreQuadrature4Weights;
		}

		~Lagrange();

	protected:
		Vector2d gaussianQuadrature2Weights;
		Vector3d gaussLegendreQuadrature3Weights;
		Vector4d gaussLobattoQuadrature4Weights;
		Vector4d gaussLegendreQuadrature4Weights;
		Vector9d gaussianQuadrature9Weights, gaussianQuadrature9Points;

	private:
		void smoothnessMeasure3rd(const Vector3d& u);
		void smoothnessMeasure2PointsGaussLegendre3th(const Vector3d& u);

		void smoothnessMeasure(const Vector5d& u);
		void smoothnessMeasureIS(const Vector5d& u);
		void smoothnessMeasureISmod(const Vector5d& u, const double& diff);
		void smoothnessMeasureISExtended(const Vector5d& u);
		void smoothnessMeasure2PointsGaussLegendre(const Vector5d& u);
		void smoothnessMeasure4PointsGaussLobatto(const Vector5d& u);
		void smoothnessMeasure4PointsGaussLegendre(const Vector5d& u);

		void fullpointsSmootness(const Vector5d& u);
		
		void cutOff();

		double catchnan(const double& input);

		/* WENO 5th coeff */
		Matrix<double, 3, 4> crj;
		Matrix<double, 3, 3> crjgm, crjgp, crjgmLobatto, crjgpLobatto, crjgmLegendre_12, crjgpLegendre_12,
			crjgmLegendre_34, crjgpLegendre_34, crjgp3L, crjgm3L;
		Vector3d dkm, dkp, dkmg, dkpg, dkmLobatto, dkpLobatto, dkmLegendre_12, dkpLegendre_12,
			dkmLegendre_34, dkpLegendre_34;

		/* WENO 3rd coeff */
		Matrix2d crjm3rd, crjp3rd, crjm3rdg, crjp3rdg;

		/* for WENO 5th */
		Vector3d SI, b, c, weights, p;
		Vector4d pis, dkmis, dkpis;
		Vector3d SIGLQ, bg, cg, SZG;
		Vector3d SILobattoA, SILobattoB, SILegendreA, SILegendreB;
		Vector6d pise, dkmise, dkpise, dkmise_pos, dkmise_neg, dkpise_pos, dkpise_neg, weightsise;
		Vector6d wise_pos, wise_neg;
		double dkmise_pos_sum, dkmise_neg_sum, dkpise_neg_sum, dkpise_pos_sum;

		/* TENO 5th */
		Vector3d sigm;

		/* for WENO 3rd */
		Vector2d SI3rd, p3rd, weights3rd, dkm3rd, dkp3rd, dk3rdg;

		/* Output vector */
		Vector2d out, outglq;

		Vector4d outGaussLobatto, SI4is, weightsis, outGaussLegendre;

		// three point gauss legendre
		double c3L1, c3L2, sumw0p_3L, sumw0m_3L;
        Vector3d outGausLegendre3n, linw0p_3L, linw0m_3L, linwp_3L, linwm_3L, SI3L;
        void smoothnessMeasure3L(const Vector5d& u);

		/* constants */
		const double eps, machine_eps;
		double sq3, sq5, sq15, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10;
		const double epsz, ct, ctg, theta, paramJ;
		double tau, b12, b35;
	};
}
