#include "Schemes_Alternative.h"

Schemes::Alternative::Alternative() {
	cmcu(0) = -3.0 / 640.0;
	cmcu(1) = 29.0 / 480.0;
	cmcu(2) = -107.0 / 960.0;
	cmcu(3) = cmcu(1);
	cmcu(4) = cmcu(0);

	cmcf(0) = -17.0 / 5760.0;
	cmcf(1) = 77.0 / 1440.0;
	cmcf(2) = -97.0 / 960.0;
	cmcf(3) = cmcf(1);
	cmcf(4) = cmcf(0);
}

double Schemes::Alternative::reCalcConservative(const Vector5d& q) {
	return q(2) - (q.dot(cmcu));
}

double Schemes::Alternative::reCalcFlux(const Vector5d& flux) {
	return flux(2) + flux.dot(cmcf);
}