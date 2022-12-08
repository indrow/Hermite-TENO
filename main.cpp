#include "FV_Scalar1D.h"
#include "FV_Scalar2D.h"
#include "FV_System1D.h"
#include "FV_System2D.h"
#include "FVHermite_Scalar1D.h"
#include "FVHermite_Scalar2D.h"
#include "FVHermite_System1D.h"
#include "FVHermite_System2D.h"
#include <omp.h>
#include "Schemes_Lagrange.h"

int main() {
	omp_set_dynamic(false);
	omp_set_num_threads(omp_get_max_threads());

//	const clock_t begin_time = clock();

    /* Linear Adv - WENO/TENO */
	/*FV::Scalar1D scalarObj(-1.0, 1.0, 320, 0.4, 2.0, 0.0);
	scalarObj.initialize();
	scalarObj.timeStepping();
	scalarObj.Exact();
	scalarObj.printResult();*/

    /* Linear Adv - HWENO/HTENO */
	/*FVHermite::Scalar1D scalarObj(-1.0, 1.0, 10, 0.4, 2.0, 0.0);
	scalarObj.initialize();
	scalarObj.timeStepping();
	scalarObj.Exact();
	scalarObj.printResult();*/

    /* Linear Burger's - WENO/TENO */
//    double pi = 4.0 * atan(1.0);
//    FV::Scalar1D scalarObj(0.0, 2.0, 80, 0.4, 1.5/pi, 0.0);
//    scalarObj.initialize();
//    scalarObj.timeStepping();
//    scalarObj.Exact();
//    scalarObj.printResult();

    /* Linear Burger's - HWENO/HTENO */
//    double pi = 4.0 * atan(1.0);
//    FVHermite::Scalar1D scalarObj(0.0, 2.0, 80, 0.4, 1.5/pi, 0.0);
//    scalarObj.initialize();
//    scalarObj.timeStepping();
//    scalarObj.Exact();
//    scalarObj.printResult();

    /* entropy wave - HWENO/HTENO */
    FVHermite::System1D systemObj(-1.0, 1.0, 960, 1.0, 0.5, 0.0, 1.4);
    systemObj.entropyWave();
    systemObj.timeStepping();
    systemObj.printResult();

    /* Euler Eqs - SOD shock tube - WENO/TENO */
//	FV::System1D systemObj(0.0, 1.0, 100, 0.4, 0.2, 0.0, 1.4);
//	systemObj.blastWave();
//	systemObj.timeStepping();
//	systemObj.printResult();

    /* Euler Eqs - SOD shock tube - HWENO/HTENO */
//    FVHermite::System1D systemObj(0.0, 1.0, 100, 0.4, 0.0, 0.0, 1.4);
//    systemObj.blastWave();
//    systemObj.timeStepping();
//    systemObj.printResult();

    /* Euler Eqs - Shu-Osher problem - WENO/TENO */
//	FV::System1D systemObj(-5.0, 5.0, 300, 0.4, 1.8, 0.0, 1.4);
//	systemObj.blastWave();
//	systemObj.timeStepping();
//	systemObj.printResult();

    /* Euler Eqs - Shu-Osher problem - HWENO/HTENO */
//    FVHermite::System1D systemObj(-5.0, 5.0, 300, 0.4, 1.8, 0.0, 1.4);
//    systemObj.blastWave();
//    systemObj.timeStepping();
//    systemObj.printResult();

    /* Euler Eqs - Interacting Blast - WENO/TENO */
//    FV::System1D systemObj(0.0, 1.0, 400, 0.4, 0.038, 0.0, 1.4);
//    systemObj.blastWave();
//    systemObj.timeStepping();
//    systemObj.printResult();

    /* Euler Eqs - Interacting Blast - HWENO/HTENO */
//    FVHermite::System1D systemObj(0.0, 1.0, 800, 0.4, 0.038, 0.0, 1.4);
//    systemObj.blastWave();
//    systemObj.timeStepping();
//    systemObj.printResult();

    /*FV::Scalar2D scalarObj(-2.0, 2.0, 320, -2.0, 2.0, 320, 0.8, 1.0, 0.0);
    scalarObj.initialize();
    scalarObj.timeStepping();
    scalarObj.Exact();
    scalarObj.printResult();*/

    /*FVHermite::Scalar2D scalarObj(-2.0, 2.0, 10, -2.0, 2.0, 10, 0.8, 1.0, 0.0);
	scalarObj.initialize();
	scalarObj.timeStepping();
	scalarObj.Exact();
	scalarObj.printResult();*/

	//FVHermite::System2D systemObj(-1.0, 1.0, 160, -1.0, 1.0, 160, 0.4, 1.0, 0.0, 1.4);
	//systemObj.entropyWave();
	//systemObj.timeStepping();
	//systemObj.printResult();

//	FVHermite::System2D systemObj1(0.0, 4.0, 1024, 0.0, 1.0, 1024/4, 0.4, 0.2, 0.0, 1.4);
//	systemObj1.Riemann2D();
//	systemObj1.timeStepping();
//	systemObj1.printResult();

//    FVHermite::System2D systemObj2(0.0, 4.0, 1024, 0.0, 1.0, 1024/4, 0.4, 0.2, 0.0, 1.4);
//    systemObj2.Riemann2D();
//    systemObj2.timeStepping();
//    systemObj2.printResult();

//	FVHermite::System2D systemObj(0.0, 0.25, 128, 0.0, 1.0, 512, 0.4, 1.95, 0.0, 5.0 / 3.0);
//	systemObj.Riemann2D();
//	systemObj.timeStepping();
//	systemObj.printResult();

//    FVHermite::System2D systemObj(0.0, 1.0, 512, 0.0, 1.0, 512, 0.4, 0.3, 0.0, 1.4);
//    systemObj.Riemann2D();
//    systemObj.timeStepping();
//    systemObj.printResult();

//    FV::System2D systemObj2(0.0, 1.0, 512, 0.0, 1.0, 512, 0.4, 0.3, 0.0, 1.4);
//    systemObj2.Riemann2D();
//    systemObj2.timeStepping();
//    systemObj2.printResult();

	/*FV::System2D systemObj(0.0, 1.0, 720, 0.0, 1.0, 720, 0.4, 0.3, 0.0, 1.4);
	systemObj.Riemann2D();
	systemObj.timeStepping();
	systemObj.printResult();*/

	/*FV::System2D systemObj(-1.0, 1.0, 40, -1.0, 1.0, 40, 0.6, 2.0, 0.0, 1.4);
	systemObj.entropyWave();
	systemObj.timeStepping();
	systemObj.printResult();*/

//	FVHermite::System2D systemObj(-5.0, 5.0, 81, -5.0, 5.0, 81, 0.4, 1000.0, 0.0, 1.4);
//    systemObj.isentropicVortex();
//    systemObj.timeStepping();
//    systemObj.printResult();


//    std::cout << float( clock () - begin_time ) /  CLOCKS_PER_SEC;

	/*Schemes::Lagrange lgrObj;
	Schemes::Hermite hmtObj;
	hmtObj.setGridSize(1.0);
	Vector6d a;
	Vector4d b, c;
	a << 1, 1, 1, 5, 5, 5;
	b << 1, 1, 10, 10;
	c << 0, 0, 0, 0;
	Vector5d d;
	d << 1, 2, 3, 4, 5;
	Vector3d e, f;
	e << 1, 2, 3;
	f << 1, 1, 1;

	std::cout << hmtObj.hteno(b, c) << endl;
	std::cout << lgrObj.weno5JS3PointsQuadrature(d);*/

//	FV::System2D systemObj(0.0, 0.25, 128, 0.0, 1.0, 512, 0.4, 0.0, 0.0, 5.0 / 3.0);
//	systemObj.Riemann2D();
////	systemObj.timeStepping();
//	systemObj.printResult();

//    FV::System1D systemObj(0.0, 1.0, 100, 0.4, 0.0, 0.0, 1.4);
//    systemObj.blastWave();
//    systemObj.printResult();

}
