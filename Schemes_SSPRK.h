#pragma once
#include "Schemes_Hermite.h"

namespace Schemes {
	class SSPRK
	{
	public:
		template <typename T> T RungeKuttaIII(const T& q0, const T& qn, const T& rhs, \
			const double& dt, const int& stage);
	};
}
