#include "Schemes_Lagrange.h"

using std::isnan;

// Constructor
Schemes::Lagrange::Lagrange() :
	eps(1E-6), machine_eps(1E-15), sq3(sqrt(3.0)), sq5(sqrt(5.0)), c1(13.0 / 12.0),
	epsz(1E-40), ct(5E-4), ctg(1E-3), theta(0.1), paramJ(2.0)
{
	double denom;

	c2 = sq3 / 12.0;
	c3 = 37.0 / 108.0 * sq3;
	c4 = sq5 / 20.0;
	c5 = 61.0 * sq5 / 300.0;
	c6 = sq3 / 3.0;
	c7 = sqrt(0.525e3 - 0.70e2 * sqrt(0.30e2)) / 0.140e3;
	c8 = pow(0.525e3 - 0.70e2 * sqrt(0.30e2), 0.3e1 / 0.2e1) / 0.514500e6 + sqrt(0.525e3 - 0.70e2 * sqrt(0.30e2)) / 0.35e2;
	c9 = sqrt(0.525e3 + 0.70e2 * sqrt(0.30e2)) / 0.140e3;
	c10 = pow(0.525e3 + 0.70e2 * sqrt(0.30e2), 0.3e1 / 0.2e1) / 0.514500e6 + sqrt(0.525e3 + 0.70e2 * sqrt(0.30e2)) / 0.35e2;

	gaussianQuadrature9Points <<
		0,
		-0.8360311073266357942994297880697348765441067181246759961043719796394550068815901188939461970258575402563758103910561868767921700399852813493611963795348388298072683628655858714286307690921827503279179493378017903390282931287792638170061442346288416366768259295268522725491437592698775616386150514997832,
		0.8360311073266357942994297880697348765441067181246759961043719796394550068815901188939461970258575402563758103910561868767921700399852813493611963795348388298072683628655858714286307690921827503279179493378017903390282931287792638170061442346288416366768259295268522725491437592698775616386150514997832,
		-0.9681602395076260898355762029036728700494048004919253295500233118490803743966007530618737492268941116024875911233178159906522811969602509341080006111457157352577320594030742939105200742221799581448832412180479160165668557217628253178605064255816845030589843605433053781978726946425719821479150514997832,
		0.9681602395076260898355762029036728700494048004919253295500233118490803743966007530618737492268941116024875911233178159906522811969602509341080006111457157352577320594030742939105200742221799581448832412180479160165668557217628253178605064255816845030589843605433053781978726946425719821479150514997832,
		-0.3242534234038089290385380146433366085719562607369730888270474768421865795351242491930986016984975672077778257173507373911718045575238432394572865005705333805025491599132630235053630398924931286361909328940173345187813296193687231694926973637651870715469270935223550274475117654585286698075150514997832,
		0.3242534234038089290385380146433366085719562607369730888270474768421865795351242491930986016984975672077778257173507373911718045575238432394572865005705333805025491599132630235053630398924931286361909328940173345187813296193687231694926973637651870715469270935223550274475117654585286698075150514997832,
		-0.6133714327005903973087020393414741847857206049405646928728129422812673464910011985832400139035685845782334895968597685619397117528519746872458346040371559996202334828312987463516926466812888532978280620182027590531371274017229787367921934803381534015176954113597402763904697814697273286917150514997832,
		0.6133714327005903973087020393414741847857206049405646928728129422812673464910011985832400139035685845782334895968597685619397117528519746872458346040371559996202334828312987463516926466812888532978280620182027590531371274017229787367921934803381534015176954113597402763904697814697273286917150514997832;

	gaussianQuadrature9Weights <<
		0.3302393550012597631645250692869740488788107835726883345930964978584026203073822121441169060216679264298311917359536407155454774502393550012597631645250692869740488788107835726883345930964978584026203073822121441169060216679264298311917359536407155454774502393550012597631645250692869740489,
		0.1806481606948574040584720312429128095143378217320404844983359064713572905449462697645949773031997041476074679602577937226796268460630127231790100804745577374812973964868278705556370432288860477148539230329025541102198218481213990057413494800065234875808239968200871271576666111786816983312150514997832,
		0.1806481606948574040584720312429128095143378217320404844983359064713572905449462697645949773031997041476074679602577937226796268460630127231790100804745577374812973964868278705556370432288860477148539230329025541102198218481213990057413494800065234875808239968200871271576666111786816983312150514997832,
		0.081274388361574411971892158110523650675661720782410750711107676880686686308452062945578554702942576957794073317963038094590048795093955759528141378844750853767333972349567507324558127938133868301667395157245896802611234739695631672003334674766636592975299135275084311484994311087346192507215051499783203,
		0.081274388361574411971892158110523650675661720782410750711107676880686686308452062945578554702942576957794073317963038094590048795093955759528141378844750853767333972349567507324558127938133868301667395157245896802611234739695631672003334674766636592975299135275084311484994311087346192507215051499783203,
		0.3123470770400028400686304065844436655987548612619046455540111655991438973240193165701219218880063538522954773181646973116391818098875271459600370901478405885572589090757645984059641355722376816546561522245422024969266380802745127735793790292136245228820749357799614002097074181144513901973150514997832,
		0.3123470770400028400686304065844436655987548612619046455540111655991438973240193165701219218880063538522954773181646973116391818098875271459600370901478405885572589090757645984059641355722376816546561522245422024969266380802745127735793790292136245228820749357799614002097074181144513901973150514997832,
		0.2606106964029354623187428694186328497718402044372999519399970021196108156688912446476460930950174018273873855356376505133184038238358268707029298682703161767070852826824482373696733967124934731275123758942032745317892944979452416330800688391928576238230768124473665313152599422632809847998150514997832,
		0.2606106964029354623187428694186328497718402044372999519399970021196108156688912446476460930950174018273873855356376505133184038238358268707029298682703161767070852826824482373696733967124934731275123758942032745317892944979452416330800688391928576238230768124473665313152599422632809847998150514997832;

	gaussianQuadrature2Weights <<
		1.0,
		1.0;

	gaussLobattoQuadrature4Weights <<
		1.0 / 6.0,
		1.0 / 6.0,
		5.0 / 6.0,
		5.0 / 6.0;

	gaussLegendreQuadrature4Weights <<
		(18.0 + sqrt(30.0)) / 36.0,
		(18.0 + sqrt(30.0)) / 36.0,
		(18.0 - sqrt(30.0)) / 36.0,
		(18.0 - sqrt(30.0)) / 36.0;

	//+++++++++++++++++++++++++++++++++
	gaussLegendreQuadrature3Weights <<
	    5.0/9.0, 8.0/9.0, 5.0/9.0;

	sq15 = sqrt(15.0);

    crjgp3L(0,0) = 1.0/30.0 + sq15/20.0;
    crjgp3L(1,0) = -1.0/15.0 - sq15/5.0;
    crjgp3L(2,0) = 31.0/30.0 + 3.0*sq15/20.0;
    crjgp3L(0,1) = 1.0/30.0 - sq15/20.0;
    crjgp3L(1,1) = 14.0/15.0;
    crjgp3L(2,1) = 1.0/30.0 + sq15/20.0;
    crjgp3L(0,2) = 31.0/30.0 - 3.0*sq15/20.0;
    crjgp3L(1,2) = -1.0/15.0 + sq15/5.0;
    crjgp3L(2,2) = 1.0/30.0 - sq15/20.0;

    crjgm3L(0,0) = 1.0/30.0 - sq15/20.0;
    crjgm3L(1,0) = -1.0/15.0 + sq15/5.0;
    crjgm3L(2,0) = 31.0/30.0 - 3.0*sq15/20.0;
    crjgm3L(0,1) = 1.0/30.0 + sq15/20.0;
    crjgm3L(1,1) = 14.0/15.0;
    crjgm3L(2,1) = 1.0/30.0 - sq15/20.0;
    crjgm3L(0,2) = 31.0/30.0 + 3.0*sq15/20.0;
    crjgm3L(1,2) = -1.0/15.0 - sq15/5.0;
    crjgm3L(2,2) = 1.0/30.0 + sq15/20.0;

    linwp_3L << -9.0/80.0, 49.0/40.0, -9.0/80.0;
    for (int i = 0; i < 3; ++i) {
        linw0p_3L(i) = 0.5 * (linwp_3L(i) + 3.0 * fabs(linwp_3L(i)));
        linw0m_3L(i) = linw0p_3L(i) - linwp_3L(i);
    }
    sumw0p_3L = linw0p_3L.sum(), sumw0m_3L = linw0m_3L.sum();
    linw0p_3L = linw0p_3L / sumw0p_3L;
    linw0m_3L = linw0m_3L / sumw0m_3L;

    denom = 15720.0*sq15 + 10480.0;
    linwp_3L(0) = 2882.0*sq15 - 1179.0;
    linwp_3L(1) = 9672.0*sq15 + 6448.0;
    linwp_3L(2) = 3166.0*sq15 + 5211.0;
    linwp_3L = linwp_3L / denom;

    denom = 15720.0 * sq15 - 10480;
    linwm_3L(0) = 2882.0*sq15 + 1179.0;
    linwm_3L(1) = 9672.0*sq15 - 6448.0;
    linwm_3L(2) = 3166.0*sq15 - 5211.0;
    linwm_3L = linwm_3L / denom;

    c3L1 = sq15/5.0; c3L2 = 21.0/100.0*sq15;
	//

	// coeff for 3rd weno
	crjm3rd(0, 0) = -0.5;
	crjm3rd(1, 0) = 1.5;
	crjm3rd(0, 1) = 0.5;
	crjm3rd(1, 1) = 0.5;
	crjp3rd = crjm3rd;
	crjp3rd.row(0).swap(crjp3rd.row(1));
	crjp3rd.col(0).swap(crjp3rd.col(1));

	dkm3rd(0) = 1.0 / 3.0;
	dkm3rd(1) = 2.0 / 3.0;
	dkp3rd(0) = dkm3rd(1);
	dkp3rd(1) = dkm3rd(0);

	// coeff for 5th weno
	crj(0, 0) = 11.0 / 6.0;
	crj(1, 0) = -7.0 / 6.0;
	crj(2, 0) = 1.0 / 3.0;
	crj(0, 1) = 1.0 / 3.0;
	crj(1, 1) = 5.0 / 6.0;
	crj(2, 1) = -1.0 / 6.0;
	crj(0, 2) = -1.0 / 6.0;
	crj(1, 2) = 5.0 / 6.0;
	crj(2, 2) = 1.0 / 3.0;
	crj(0, 3) = 1.0 / 3.0;
	crj(1, 3) = -7.0 / 6.0;
	crj(2, 3) = 11.0 / 6.0;

	// notation m -> x^{-}_{i+1/2}
	// notation p -> x^{+}_{i-1/2}
	dkm(0) = 0.3; dkm(1) = 0.6; dkm(2) = 0.1;
	dkp = dkm;
	dkp.row(0).swap(dkp.row(2));

	// WENO5-IS optimal weights
	dkmis(0) = 0.4;
	dkmis(1) = 0.2;
	dkmis(2) = 0.3;
	dkmis(3) = 0.1;

	dkpis(0) = 0.2;
	dkpis(1) = 0.4;
	dkpis(2) = 0.1;
	dkpis(3) = 0.3;

	// notation p -> x_{i+xG}
	crjgp(0, 0) = sq3 / 12.0;
	crjgp(1, 0) = -sq3 / 3.0;
	crjgp(2, 0) = 1.0 + sq3 / 4.0;
	crjgp(0, 1) = -sq3 / 12.0;
	crjgp(1, 1) = 1.0;
	crjgp(2, 1) = sq3 / 12.0;
	crjgp(0, 2) = -sq3 / 4.0 + 1.0;
	crjgp(1, 2) = sq3 / 3.0;
	crjgp(2, 2) = -sq3 / 12.0;

	// notation m -> x_{i-xG}
	crjgm = crjgp;
	crjgm.col(0).swap(crjgm.col(2));
	crjgm.row(0).swap(crjgm.row(2));

	// Third order WENO two-points quadrature
	crjp3rdg(0, 0) = -sq3 / 6.0;
	crjp3rdg(1, 0) = sq3 / 6.0 + 1.0;
	crjp3rdg(0, 1) = -sq3 / 6.0 + 1.0;
	crjp3rdg(1, 1) = sq3 / 6.0;
	crjm3rdg = crjp3rdg;
	crjm3rdg.row(0).swap(crjm3rdg.row(1));
	crjm3rdg.col(0).swap(crjm3rdg.col(1));

	dk3rdg(0) = 0.5;
	dk3rdg(1) = 0.5;

	// Gauss-Lobatto Quadrature for x_{i+-sqrt(5.0)/10.0}
	crjgpLobatto(0, 0) = sq5 / 20.0 - 1.0 / 60.0;
	crjgpLobatto(1, 0) = -sq5 / 5.0 + 1.0 / 30.0;
	crjgpLobatto(2, 0) = 3.0 * sq5 / 20.0 + 59 / 60.0;
	crjgpLobatto(0, 1) = -sq5 / 20.0 - 1.0 / 60.0;
	crjgpLobatto(1, 1) = 31.0 / 30.0;
	crjgpLobatto(2, 1) = sq5 / 20.0 - 1.0 / 60.0;
	crjgpLobatto(0, 2) = -3.0 * sq5 / 20.0 + 59.0 / 60.0;
	crjgpLobatto(1, 2) = sq5 / 5.0 + 1.0 / 30.0;
	crjgpLobatto(2, 2) = -sq5 / 20.0 - 1.0 / 60.0;

	crjgmLobatto = crjgpLobatto;
	crjgmLobatto.row(0).swap(crjgmLobatto.row(2));
	crjgmLobatto.col(0).swap(crjgmLobatto.col(2));

	denom = 660.0 * sq5 - 220.0;
	dkpLobatto(0) = 132.0 * sq5 + 22.0;
	dkpLobatto(1) = 387.0 * sq5 - 129.0;
	dkpLobatto(2) = 141.0 * sq5 - 113.0;
	dkpLobatto = dkpLobatto / denom;

	denom = 660.0 * sq5 + 220.0;
	dkmLobatto(0) = 132.0 * sq5 - 22.0;
	dkmLobatto(1) = 387.0 * sq5 + 129.0;
	dkmLobatto(2) = 141.0 * sq5 + 113.0;
	dkmLobatto = dkmLobatto / denom;

	denom = 1080;
	dkpg(0) = (210.0 - sq3) / denom;
	dkpg(1) = 11.0 / 18.0;
	dkpg(2) = (210.0 + sq3) / denom;

	dkmg = dkpg;
	dkmg.row(0).swap(dkmg.row(2));

	// Four-points Gauss-Legendre quadrature
	crjgpLegendre_12(0, 0) = sqrt(0.525e3 - 0.70e2 * sqrt(0.30e2)) / 0.140e3 - sqrt(0.30e2) / 0.140e3 + 0.1e1 / 0.84e2;
	crjgpLegendre_12(1, 0) = -sqrt(0.525e3 - 0.70e2 * sqrt(0.30e2)) / 0.35e2 + sqrt(0.30e2) / 0.70e2 - 0.1e1 / 0.42e2;
	crjgpLegendre_12(2, 0) = 0.3e1 / 0.140e3 * sqrt(0.525e3 - 0.70e2 * sqrt(0.30e2)) - sqrt(0.30e2) / 0.140e3 + 0.85e2 / 0.84e2;
	crjgpLegendre_12(0, 1) = -sqrt(0.525e3 - 0.70e2 * sqrt(0.30e2)) / 0.140e3 - sqrt(0.30e2) / 0.140e3 + 0.1e1 / 0.84e2;
	crjgpLegendre_12(1, 1) = sqrt(0.30e2) / 0.70e2 + 0.41e2 / 0.42e2;
	crjgpLegendre_12(2, 1) = sqrt(0.525e3 - 0.70e2 * sqrt(0.30e2)) / 0.140e3 - sqrt(0.30e2) / 0.140e3 + 0.1e1 / 0.84e2;
	crjgpLegendre_12(0, 2) = -0.3e1 / 0.140e3 * sqrt(0.525e3 - 0.70e2 * sqrt(0.30e2)) - sqrt(0.30e2) / 0.140e3 + 0.85e2 / 0.84e2;
	crjgpLegendre_12(1, 2) = sqrt(0.525e3 - 0.70e2 * sqrt(0.30e2)) / 0.35e2 + sqrt(0.30e2) / 0.70e2 - 0.1e1 / 0.42e2;
	crjgpLegendre_12(2, 2) = -sqrt(0.525e3 - 0.70e2 * sqrt(0.30e2)) / 0.140e3 - sqrt(0.30e2) / 0.140e3 + 0.1e1 / 0.84e2;

	crjgmLegendre_12 = crjgpLegendre_12;
	crjgmLegendre_12.row(0).swap(crjgmLegendre_12.row(2));
	crjgmLegendre_12.col(0).swap(crjgmLegendre_12.col(2));

	dkpLegendre_12(2) = ((0.257250000e9 * sqrt(0.30e2) - 0.1899362500e10) * sqrt(0.525e3 - 0.70e2 * sqrt(0.30e2)) +
		(-0.42000e5 * sqrt(0.30e2) + 0.310100e6) * pow(0.525e3 - 0.70e2 * sqrt(0.30e2), 0.3e1 / 0.2e1) - 0.1197315000e10 *
		sqrt(0.30e2) + 0.5597025000e10) / ((0.1234800000e10 * sqrt(0.30e2) - 0.9116940000e10) * sqrt(0.525e3 - 0.70e2 *
			sqrt(0.30e2)) + 0.11174940000e11 * sqrt(0.30e2) - 0.52238900000e11);
	dkpLegendre_12(1) = ((0.1761795000e10 * sqrt(0.30e2) - 0.13038532500e11) * sqrt(0.525e3 - 0.70e2 * sqrt(0.30e2)) +
		(-0.88200e5 * sqrt(0.30e2) + 0.16096500e8) * pow(0.525e3 - 0.70e2 * sqrt(0.30e2), 0.3e1 / 0.2e1) - 0.2520e4 *
		pow(0.525e3 - 0.70e2 * sqrt(0.30e2), 0.5e1 / 0.2e1) + 0.6757590000e10 * sqrt(0.30e2) - 0.32331425000e11) /
		((0.1234800000e10 * sqrt(0.30e2) - 0.9116940000e10) * sqrt(0.525e3 - 0.70e2 * sqrt(0.30e2)) + 0.11174940000e11 *
			sqrt(0.30e2) - 0.52238900000e11);
	dkpLegendre_12(0) = ((0.26092500e8 * sqrt(0.30e2) - 0.256576250e9) * sqrt(0.525e3 - 0.70e2 * sqrt(0.30e2)) +
		(-0.2100e4 * sqrt(0.30e2) - 0.3838100e7) * pow(0.525e3 - 0.70e2 * sqrt(0.30e2), 0.3e1 / 0.2e1) + 0.630e3 *
		pow(0.525e3 - 0.70e2 * sqrt(0.30e2), 0.5e1 / 0.2e1) + 0.5614665000e10 * sqrt(0.30e2) - 0.25504500000e11) /
		((0.1234800000e10 * sqrt(0.30e2) - 0.9116940000e10) * sqrt(0.525e3 - 0.70e2 * sqrt(0.30e2)) + 0.11174940000e11 *
			sqrt(0.30e2) - 0.52238900000e11);

	dkmLegendre_12 = dkpLegendre_12;
	dkmLegendre_12.row(0).swap(dkmLegendre_12.row(2));

	crjgpLegendre_34(0, 0) = sqrt(0.525e3 + 0.70e2 * sqrt(0.30e2)) / 0.140e3 + sqrt(0.30e2) / 0.140e3 + 0.1e1 / 0.84e2;
	crjgpLegendre_34(1, 0) = -sqrt(0.525e3 + 0.70e2 * sqrt(0.30e2)) / 0.35e2 - sqrt(0.30e2) / 0.70e2 - 0.1e1 / 0.42e2;
	crjgpLegendre_34(2, 0) = 0.3e1 / 0.140e3 * sqrt(0.525e3 + 0.70e2 * sqrt(0.30e2)) + sqrt(0.30e2) / 0.140e3 + 0.85e2 / 0.84e2;
	crjgpLegendre_34(0, 1) = -sqrt(0.525e3 + 0.70e2 * sqrt(0.30e2)) / 0.140e3 + sqrt(0.30e2) / 0.140e3 + 0.1e1 / 0.84e2;
	crjgpLegendre_34(1, 1) = -sqrt(0.30e2) / 0.70e2 + 0.41e2 / 0.42e2;
	crjgpLegendre_34(2, 1) = sqrt(0.525e3 + 0.70e2 * sqrt(0.30e2)) / 0.140e3 + sqrt(0.30e2) / 0.140e3 + 0.1e1 / 0.84e2;
	crjgpLegendre_34(0, 2) = -0.3e1 / 0.140e3 * sqrt(0.525e3 + 0.70e2 * sqrt(0.30e2)) + sqrt(0.30e2) / 0.140e3 + 0.85e2 / 0.84e2;
	crjgpLegendre_34(1, 2) = sqrt(0.525e3 + 0.70e2 * sqrt(0.30e2)) / 0.35e2 - sqrt(0.30e2) / 0.70e2 - 0.1e1 / 0.42e2;
	crjgpLegendre_34(2, 2) = -sqrt(0.525e3 + 0.70e2 * sqrt(0.30e2)) / 0.140e3 + sqrt(0.30e2) / 0.140e3 + 0.1e1 / 0.84e2;

	crjgmLegendre_34 = crjgpLegendre_34;
	crjgmLegendre_34.row(0).swap(crjgmLegendre_34.row(2));
	crjgmLegendre_34.col(0).swap(crjgmLegendre_34.col(2));

	dkpLegendre_34(2) = ((0.26092500e8 * sqrt(0.30e2) + 0.256576250e9) * sqrt(0.525e3 + 0.70e2 * sqrt(0.30e2)) +
		(-0.2100e4 * sqrt(0.30e2) + 0.3838100e7) * pow(0.525e3 + 0.70e2 * sqrt(0.30e2), 0.3e1 / 0.2e1) -
		0.630e3 * pow(0.525e3 + 0.70e2 * sqrt(0.30e2), 0.5e1 / 0.2e1) - 0.5614665000e10 * sqrt(0.30e2) -
		0.25504500000e11) / ((0.1234800000e10 * sqrt(0.30e2) + 0.9116940000e10) * sqrt(0.525e3 + 0.70e2 * sqrt(0.30e2))
			- 0.11174940000e11 * sqrt(0.30e2) - 0.52238900000e11);
	dkpLegendre_34(1) = ((0.1761795000e10 * sqrt(0.30e2) + 0.13038532500e11) * sqrt(0.525e3 + 0.70e2 * sqrt(0.30e2)) +
		(-0.88200e5 * sqrt(0.30e2) - 0.16096500e8) * pow(0.525e3 + 0.70e2 * sqrt(0.30e2), 0.3e1 / 0.2e1) +
		0.2520e4 * pow(0.525e3 + 0.70e2 * sqrt(0.30e2), 0.5e1 / 0.2e1) - 0.6757590000e10 * sqrt(0.30e2) -
		0.32331425000e11) / ((0.1234800000e10 * sqrt(0.30e2) + 0.9116940000e10) * sqrt(0.525e3 + 0.70e2 * sqrt(0.30e2)) -
			0.11174940000e11 * sqrt(0.30e2) - 0.52238900000e11);
	dkpLegendre_34(0) = ((0.257250000e9 * sqrt(0.30e2) + 0.1899362500e10) * sqrt(0.525e3 + 0.70e2 * sqrt(0.30e2)) +
		(-0.42000e5 * sqrt(0.30e2) - 0.310100e6) * pow(0.525e3 + 0.70e2 * sqrt(0.30e2), 0.3e1 / 0.2e1) +
		0.1197315000e10 * sqrt(0.30e2) + 0.5597025000e10) / ((0.1234800000e10 * sqrt(0.30e2) +
			0.9116940000e10) * sqrt(0.525e3 + 0.70e2 * sqrt(0.30e2)) - 0.11174940000e11 * sqrt(0.30e2) - 0.52238900000e11);

	dkmLegendre_34 = dkpLegendre_34;
	dkmLegendre_34.row(0).swap(dkmLegendre_34.row(2));

	// WENO-IS
	dkmise << -1.0 / 45.0, 11.0 / 90.0, 1.0 / 5.0, 2.0 / 5.0, 1.0 / 5.0, 1.0 / 10.0;
	dkpise << 1.0 / 10.0, 1.0 / 5.0, 4.0 / 10.0, 2.0 / 10.0, 11.0 / 90.0, -1.0 / 45.0;

	for (int i = 0; i < 6; ++i) {
		dkmise_pos(i) = 0.5 * (dkmise(i) + 3.0 * fabs(dkmise(i)));
		dkmise_neg(i) = dkmise_pos(i) - dkmise(i);

		dkpise_pos(i) = 0.5 * (dkpise(i) + 3.0 * fabs(dkpise(i)));
		dkpise_neg(i) = dkpise_pos(i) - dkpise(i);
	}
	dkmise_pos_sum = dkmise_pos.sum();
	dkmise_pos = dkmise_pos / dkmise_pos_sum;
	dkmise_neg_sum = dkmise_neg.sum();
	dkmise_neg = dkmise_neg / dkmise_neg_sum;
	dkpise_pos_sum = dkpise_pos.sum();
	dkpise_pos = dkpise_pos / dkpise_pos_sum;
	dkpise_neg_sum = dkpise_neg.sum();
	dkpise_neg = dkpise_neg / dkpise_neg_sum;

	tau = 0.0; b12 = 0.0; b35 = 0.0;
}

Vector9d Schemes::Lagrange::getQuadrature9Weights() {
	return gaussianQuadrature9Weights;
}

Vector9d Schemes::Lagrange::getQuadrature9Points() {
	return gaussianQuadrature9Points;
}

// Third order WENO-JS scheme
Vector2d Schemes::Lagrange::weno3JS(const Vector4d& u) {
	p3rd(0) = crjm3rd.col(0).dot(u.segment(0, 2));
	p3rd(1) = crjm3rd.col(1).dot(u.segment(1, 2));
	smoothnessMeasure3rd(u.segment(0, 3));

	weights3rd = dkm3rd.array() / (eps + SI3rd.array()).pow(2);
	weights3rd = weights3rd.array() / weights3rd.sum();
	out(0) = p3rd.dot(weights3rd);
	out(0) = catchnan(out(0));

	p3rd(0) = crjp3rd.col(0).dot(u.segment(1, 2));
	p3rd(1) = crjp3rd.col(1).dot(u.segment(2, 2));
	smoothnessMeasure3rd(u.segment(1, 3));

	weights3rd = dkp3rd.array() / (eps + SI3rd.array()).pow(2);
	weights3rd = weights3rd.array() / weights3rd.sum();
	out(1) = p3rd.dot(weights3rd);
	out(1) = catchnan(out(0));

	return out;
}

Vector2d Schemes::Lagrange::weno3JS2PointsQuadrature(const Vector3d& u) {
	smoothnessMeasure2PointsGaussLegendre3th(u);

	p3rd(0) = crjp3rdg.col(0).dot(u.segment(0, 2));
	p3rd(1) = crjp3rdg.col(1).dot(u.segment(1, 2));

	weights3rd = dk3rdg.array() / (eps + SI3rd.array()).pow(2);
	weights3rd = weights3rd.array() / weights3rd.sum();
	out(0) = p3rd.dot(weights3rd);

	p3rd(0) = crjm3rdg.col(0).dot(u.segment(0, 2));
	p3rd(1) = crjm3rdg.col(1).dot(u.segment(1, 2));

	weights3rd = dk3rdg.array() / (eps + SI3rd.array()).pow(2);
	weights3rd = weights3rd.array() / weights3rd.sum();
	out(1) = p3rd.dot(weights3rd);

	return out;
}

// Fifth order WENO-JS scheme
Vector2d Schemes::Lagrange::weno5JS(const Vector6d& u) {
	p(0) = crj.col(1).dot(u.segment(2, 3));
	p(1) = crj.col(2).dot(u.segment(1, 3));
	p(2) = crj.col(3).dot(u.segment(0, 3));

	smoothnessMeasure(u.segment(0, 5));
	weights = dkm.array() / (eps + SI.array()).pow(2);
	weights = weights.array() / weights.sum();
	out(0) = p.dot(weights);
	out(0) = catchnan(out(0));

	p(0) = crj.col(0).dot(u.segment(3, 3));
	p(1) = crj.col(1).dot(u.segment(2, 3));
	p(2) = crj.col(2).dot(u.segment(1, 3));

	smoothnessMeasure(u.segment(1, 5));
	weights = dkp.array() / (eps + SI.array()).pow(2);
	weights = weights.array() / weights.sum();
	out(1) = p.dot(weights);
	out(1) = catchnan(out(1));

	return out;
}

// Fifth order WENO-Z scheme
Vector2d Schemes::Lagrange::weno5Z(const Vector6d& u) {
	p(0) = crj.col(1).dot(u.segment(2, 3));
	p(1) = crj.col(2).dot(u.segment(1, 3));
	p(2) = crj.col(3).dot(u.segment(0, 3));

	smoothnessMeasure(u.segment(0, 5));
	tau = fabs(SI(2) - SI(0));
	weights = dkm.array() * (1.0 + (tau / (SI.array() + epsz)));
	weights = weights.array() / weights.sum();
	out(0) = p.dot(weights);
	out(0) = catchnan(out(0));

	p(0) = crj.col(0).dot(u.segment(3, 3));
	p(1) = crj.col(1).dot(u.segment(2, 3));
	p(2) = crj.col(2).dot(u.segment(1, 3));

	smoothnessMeasure(u.segment(1, 5));
	tau = fabs(SI(2) - SI(0));
	weights = dkp.array() * (1.0 + (tau / (SI.array() + epsz)));
	weights = weights.array() / weights.sum();
	out(1) = p.dot(weights);
	out(1) = catchnan(out(1));

	return out;
}

// Fifth order WENO-Z scheme
Vector2d Schemes::Lagrange::weno5Z(const Vector6d& u, const double& diff) {
	p(0) = crj.col(1).dot(u.segment(2, 3));
	p(1) = crj.col(2).dot(u.segment(1, 3));
	p(2) = crj.col(3).dot(u.segment(0, 3));

	smoothnessMeasure(u.segment(0, 5));
	fullpointsSmootness(u.segment(0, 5));
	tau = fabs(b35 - (SI(0) + SI(2) + 4.0 * SI(1)) / 6.0);
	weights = (1.0 + tau / pow(SI.array() + epsz, 0.75));
	weights = weights.cwiseProduct(dkm);
	weights /= weights.sum();
	out(0) = p.dot(weights);
	out(0) = catchnan(out(0));

	p(0) = crj.col(0).dot(u.segment(3, 3));
	p(1) = crj.col(1).dot(u.segment(2, 3));
	p(2) = crj.col(2).dot(u.segment(1, 3));

	smoothnessMeasure(u.segment(1, 5));
	fullpointsSmootness(u.segment(0, 5));
	tau = fabs(b35 - (SI(0) + SI(2) + 4.0 * SI(1)) / 6.0);
	weights = (1.0 + tau / pow(SI.array() + epsz, 0.75));
	weights = weights.cwiseProduct(dkp);
	weights /= weights.sum();
	out(1) = p.dot(weights);
	out(1) = catchnan(out(1));

	return out;
}

Vector2d Schemes::Lagrange::weno5N(const Vector6d& u) {
	p(0) = crj.col(1).dot(u.segment(2, 3));
	p(1) = crj.col(2).dot(u.segment(1, 3));
	p(2) = crj.col(3).dot(u.segment(0, 3));

	SI(2) = theta * fabs(u(0) - 3.0 * u(1) + 2.0 * u(2)) +
			fabs(u(0) - 2.0 * u(1) + u(2));
	SI(1) = theta * fabs(u(3) - u(2)) +
			fabs(u(1) - 2.0 * u(2) + u(3));
	SI(0) = theta * fabs(u(3) - u(2)) +
			fabs(u(2) - 2.0 * u(3) + u(4));

	tau = pow(fabs(u(0) - 4.0 * u(1) + 6.0 * u(2) - 4.0 * u(3) + u(4)), 2.5);

	weights = dkm.array() * (1.0 + (tau / pow(SI.array() + eps, 2)));
	weights = weights.array() / weights.sum();
	out(0) = p.dot(weights);
	out(0) = catchnan(out(0));

	p(0) = crj.col(0).dot(u.segment(3, 3));
	p(1) = crj.col(1).dot(u.segment(2, 3));
	p(2) = crj.col(2).dot(u.segment(1, 3));

	SI(2) = theta * fabs(u(1) - 3.0 * u(2) + 2.0 * u(3)) +
			fabs(u(1) - 2.0 * u(2) + u(3));
	SI(1) = theta * fabs(u(4) - u(3)) +
			fabs(u(2) - 2.0 * u(3) + u(4));
	SI(0) = theta * fabs(u(4) - u(3)) +
			fabs(u(3) - 2.0 * u(4) + u(5));

	tau = pow(fabs(u(1) - 4.0 * u(2) + 6.0 * u(3) - 4.0 * u(4) + u(5)), 2.5);

	weights = dkp.array() * (1.0 + (tau / pow(SI.array() + eps, 2)));
	weights = weights.array() / weights.sum();
	out(1) = p.dot(weights);
	out(1) = catchnan(out(1));

	return out;
}

// Fifth order TENO scheme
Vector2d Schemes::Lagrange::teno5(const Vector6d& u) {
	p(0) = crj.col(1).dot(u.segment(2, 3));
	p(1) = crj.col(2).dot(u.segment(1, 3));
	p(2) = crj.col(3).dot(u.segment(0, 3));

	smoothnessMeasure(u.segment(0, 5));
	fullpointsSmootness(u.segment(0, 5));
	tau = fabs(b35 - (SI(0) + SI(2) + 4.0 * SI(1)) / 6.0);
	weights = (1.0 + tau / (SI.array() + epsz)).pow(6);
	cutOff();
	weights = sigm.cwiseProduct(dkm);
	weights /= weights.sum();
	out(0) = p.dot(weights);
	out(0) = catchnan(out(0));

	p(0) = crj.col(0).dot(u.segment(3, 3));
	p(1) = crj.col(1).dot(u.segment(2, 3));
	p(2) = crj.col(2).dot(u.segment(1, 3));

	smoothnessMeasure(u.segment(1, 5));
	fullpointsSmootness(u.segment(1, 5));
	tau = fabs(b35 - (SI(0) + SI(2) + 4.0 * SI(1)) / 6.0);
	weights = (1.0 + tau / (SI.array() + epsz)).pow(6);
	cutOff();
	weights = sigm.cwiseProduct(dkp);
	weights /= weights.sum();
	out(1) = p.dot(weights);
	out(1) = catchnan(out(1));

	return out;
}

// Fifth order TENO scheme
Vector2d Schemes::Lagrange::teno5(const Vector6d& u, const double& diff) {
//	double df = epsz;
	double df = pow(diff, 4);
	p(0) = crj.col(1).dot(u.segment(2, 3));
	p(1) = crj.col(2).dot(u.segment(1, 3));
	p(2) = crj.col(3).dot(u.segment(0, 3));

	smoothnessMeasure(u.segment(0, 5));
	fullpointsSmootness(u.segment(0, 5));
	tau = fabs(b35 - (SI(0) + SI(2) + 4.0 * SI(1)) / 6.0);
	//	tau = fabs(SI(0) - SI(2));

	weights = (1.0 + pow(tau / (SI.array() + df), 2)).pow(3);
	cutOff();
	weights = sigm.cwiseProduct(dkm);
	weights /= weights.sum();
	out(0) = p.dot(weights);
	out(0) = catchnan(out(0));

	p(0) = crj.col(0).dot(u.segment(3, 3));
	p(1) = crj.col(1).dot(u.segment(2, 3));
	p(2) = crj.col(2).dot(u.segment(1, 3));

	smoothnessMeasure(u.segment(1, 5));
	fullpointsSmootness(u.segment(1, 5));
	tau = fabs(b35 - (SI(0) + SI(2) + 4.0 * SI(1)) / 6.0);
	//	tau = fabs(SI(0) - SI(2));

	weights = (1.0 + pow(tau / (SI.array() + df), 2)).pow(3);
	cutOff();
	weights = sigm.cwiseProduct(dkp);
	weights /= weights.sum();
	out(1) = p.dot(weights);
	out(1) = catchnan(out(1));

	return out;
}

// Fifth order TENO-N scheme
Vector2d Schemes::Lagrange::teno5N(const Vector6d& u) {
    p(0) = crj.col(1).dot(u.segment(2, 3));
    p(1) = crj.col(2).dot(u.segment(1, 3));
    p(2) = crj.col(3).dot(u.segment(0, 3));

    SI(2) = theta * fabs(u(0) - 3.0 * u(1) + 2.0 * u(2)) +
            fabs(u(0) - 2.0 * u(1) + u(2));
    SI(1) = theta * fabs(u(3) - u(2)) +
            fabs(u(1) - 2.0 * u(2) + u(3));
    SI(0) = theta * fabs(u(3) - u(2)) +
            fabs(u(2) - 2.0 * u(3) + u(4));

    tau = pow(fabs(u(0) - 4.0 * u(1) + 6.0 * u(2) - 4.0 * u(3) + u(4)), 2);
    weights = (1.0 + tau / pow(SI.array() + eps, 2)).pow(6);
    cutOff();
    weights = sigm.cwiseProduct(dkm);
    weights /= weights.sum();
    out(0) = p.dot(weights);
    out(0) = catchnan(out(0));

    p(0) = crj.col(0).dot(u.segment(3, 3));
    p(1) = crj.col(1).dot(u.segment(2, 3));
    p(2) = crj.col(2).dot(u.segment(1, 3));

    SI(2) = theta * fabs(u(1) - 3.0 * u(2) + 2.0 * u(3)) +
            fabs(u(1) - 2.0 * u(2) + u(3));
    SI(1) = theta * fabs(u(4) - u(3)) +
            fabs(u(2) - 2.0 * u(3) + u(4));
    SI(0) = theta * fabs(u(4) - u(3)) +
            fabs(u(3) - 2.0 * u(4) + u(5));

    tau = pow(fabs(u(1) - 4.0 * u(2) + 6.0 * u(3) - 4.0 * u(4) + u(5)), 2);
    weights = (1.0 + tau / pow(SI.array() + eps, 2)).pow(6);
    cutOff();
    weights = sigm.cwiseProduct(dkp);
    weights /= weights.sum();
    out(1) = p.dot(weights);
    out(1) = catchnan(out(1));

    return out;
}

// Fifth order WENO-IS scheme
Vector2d Schemes::Lagrange::weno5IS(const Vector6d& u) {
	pis(0) = 0.5 * (u(2) + u(3));
	pis(1) = -0.5 * u(1) + 1.5 * u(2);
	pis(2) = crj.col(1).dot(u.segment(2, 3));
	pis(3) = crj.col(3).dot(u.segment(0, 3));

	smoothnessMeasureIS(u.segment(0, 5));
	weightsis = dkmis.cwiseProduct(weightsis);
	weightsis /= weightsis.sum();
	out(0) = pis.dot(weightsis);
	out(0) = catchnan(out(0));

	pis(0) = 1.5 * u(3) - 0.5 * u(4);
	pis(1) = 0.5 * (u(2) + u(3));
	pis(2) = crj.col(0).dot(u.segment(3, 3));
	pis(3) = crj.col(2).dot(u.segment(1, 3));

	smoothnessMeasureIS(u.segment(1, 5));
	weightsis = dkpis.cwiseProduct(weightsis);
	weightsis /= weightsis.sum();
	out(1) = pis.dot(weightsis);
	out(1) = catchnan(out(1));

	return out;
}

Vector2d Schemes::Lagrange::weno5ISmod(const Vector6d& u, const double& diff) {
	double df = pow(diff, 3);
	pis(0) = 0.5 * (u(2) + u(3));
	pis(1) = -0.5 * u(1) + 1.5 * u(2);
	pis(2) = crj.col(1).dot(u.segment(2, 3));
	pis(3) = crj.col(3).dot(u.segment(0, 3));

	smoothnessMeasureISmod(u.segment(0, 5), df);
	weightsis = dkmis.cwiseProduct(weightsis);
	weightsis /= weightsis.sum();
	out(0) = pis.dot(weightsis);
	out(0) = catchnan(out(0));

	pis(0) = 1.5 * u(3) - 0.5 * u(4);
	pis(1) = 0.5 * (u(2) + u(3));
	pis(2) = crj.col(0).dot(u.segment(3, 3));
	pis(3) = crj.col(2).dot(u.segment(1, 3));

	smoothnessMeasureISmod(u.segment(1, 5), df);
	weightsis = dkpis.cwiseProduct(weightsis);
	weightsis /= weightsis.sum();
	out(1) = pis.dot(weightsis);
	out(1) = catchnan(out(1));

	return out;
}

Vector2d Schemes::Lagrange::weno5IS_Extended(const Vector6d& u) {
	pise(0) = -1.5 * u(0) + 2.5 * u(1);
	pise(1) = -0.5 * u(1) + 1.5 * u(2);
	pise(2) = pise(1);
	pise(3) = 0.5 * (u(2) + u(3));
	pise(4) = pise(3);
	pise(5) = 1.5 * u(3) - 0.5 * u(4);

	smoothnessMeasureISExtended(u.segment(0, 5));
	tau *= tau;

	weightsise(0) = 1.0 + tau / (SI(0) + eps) / (SI4is(0) + eps);
	weightsise(1) = 1.0 + tau / (SI(0) + eps) / (SI4is(1) + eps);
	weightsise(2) = 1.0 + tau / (SI(1) + eps) / (SI4is(1) + eps);
	weightsise(3) = 1.0 + tau / (SI(1) + eps) / (SI4is(2) + eps);
	weightsise(4) = 1.0 + tau / (SI(2) + eps) / (SI4is(2) + eps);
	weightsise(5) = 1.0 + tau / (SI(2) + eps) / (SI4is(3) + eps);
	weightsise = weightsise.array().pow(4);
	weightsise /= weightsise.sum();
	for (int i = 0; i < 6; i++) {
		weightsise(i) = weightsise(i) < 1E-5 ? 0.0 : 1.0;
	}

	wise_pos = dkmise_pos.cwiseProduct(weightsise);
	wise_pos = wise_pos / wise_pos.sum();

	wise_neg = dkmise_neg.cwiseProduct(weightsise);
	wise_neg = wise_neg / wise_neg.sum();

	out(0) = dkmise_pos_sum * pise.dot(wise_pos) -
		dkmise_neg_sum * pise.dot(wise_neg);
	out(0) = catchnan(out(0));

	pise(0) = -0.5 * u(1) + 1.5 * u(2);
	pise(1) = 0.5 * (u(2) + u(3));
	pise(2) = pise(1);
	pise(3) = 1.5 * u(3) - 0.5 * u(4);
	pise(4) = pise(3);
	pise(5) = 2.5 * u(4) - 1.5 * u(5);

	smoothnessMeasureISExtended(u.segment(1, 5));
	tau *= tau;

	weightsise(0) = 1.0 + tau / (SI(0) + eps) / (SI4is(0) + eps);
	weightsise(1) = 1.0 + tau / (SI(0) + eps) / (SI4is(1) + eps);
	weightsise(2) = 1.0 + tau / (SI(1) + eps) / (SI4is(1) + eps);
	weightsise(3) = 1.0 + tau / (SI(1) + eps) / (SI4is(2) + eps);
	weightsise(4) = 1.0 + tau / (SI(2) + eps) / (SI4is(2) + eps);
	weightsise(5) = 1.0 + tau / (SI(2) + eps) / (SI4is(3) + eps);
	weightsise = weightsise.array().pow(4);
	weightsise /= weightsise.sum();

	for (int i = 0; i < 6; i++) {
		weightsise(i) = weightsise(i) < 1E-5 ? 0.0 : 1.0;
	}

	wise_pos = dkpise_pos.cwiseProduct(weightsise);
	wise_pos = wise_pos / wise_pos.sum();

	wise_neg = dkpise_neg.cwiseProduct(weightsise);
	wise_neg = wise_neg / wise_neg.sum();

	out(1) = dkpise_pos_sum * pise.dot(wise_pos) -
		dkpise_neg_sum * pise.dot(wise_neg);
	out(1) = catchnan(out(1));

	return out;
}

// WENO 5th JS for interpolation at integration points using 2 points
Matrix<double, 2, 1> Schemes::Lagrange::weno5JS2PointsQuadrature(const Vector5d& u) {
	smoothnessMeasure2PointsGaussLegendre(u);

	p(0) = crjgp.col(0).dot(u.segment(0, 3));
	p(1) = crjgp.col(1).dot(u.segment(1, 3));
	p(2) = crjgp.col(2).dot(u.segment(2, 3));

	weights = dkpg.array() / (eps + SIGLQ.array()).pow(2);
	weights = weights.array() / weights.sum();
	outglq(0) = p.dot(weights);

	p(0) = crjgm.col(0).dot(u.segment(0, 3));
	p(1) = crjgm.col(1).dot(u.segment(1, 3));
	p(2) = crjgm.col(2).dot(u.segment(2, 3));

	weights = dkmg.array() / (eps + SIGLQ.array()).pow(2);
	weights = weights.array() / weights.sum();
	outglq(1) = p.dot(weights);

	return outglq;
}

// WENO 5th Z for interpolation at integration points using 2 points
Matrix<double, 2, 1> Schemes::Lagrange::weno5Z2PointsQuadrature(const Vector5d& u) {
	smoothnessMeasure2PointsGaussLegendre(u);
	tau = fabs(SIGLQ(0) - SIGLQ(2));
	SZG = (1.0 + tau / (SIGLQ.array() + epsz));

	p(0) = crjgp.col(0).dot(u.segment(0, 3));
	p(1) = crjgp.col(1).dot(u.segment(1, 3));
	p(2) = crjgp.col(2).dot(u.segment(2, 3));

	weights = dkpg.cwiseProduct(SZG);
	weights = weights.array() / weights.sum();
	outglq(0) = p.dot(weights);

	p(0) = crjgm.col(0).dot(u.segment(0, 3));
	p(1) = crjgm.col(1).dot(u.segment(1, 3));
	p(2) = crjgm.col(2).dot(u.segment(2, 3));

	weights = dkmg.cwiseProduct(SZG);
	weights = weights.array() / weights.sum();
	outglq(1) = p.dot(weights);

	return outglq;
}

// WENO 5th JS for interpolation at integration points using 2 points
Matrix<double, 2, 1> Schemes::Lagrange::teno5th2PointsQuadrature(const Vector5d& u) {
	smoothnessMeasure2PointsGaussLegendre(u);

	p(0) = crjgp.col(0).dot(u.segment(0, 3));
	p(1) = crjgp.col(1).dot(u.segment(1, 3));
	p(2) = crjgp.col(2).dot(u.segment(2, 3));
	tau = fabs(SIGLQ(0) - SIGLQ(2));
	weights = (1.0 + tau / (SIGLQ.array() + epsz)).pow(6);
	cutOff();
	weights = sigm.cwiseProduct(dkpg);
	weights /= weights.sum();
	outglq(0) = p.dot(weights);

	p(0) = crjgm.col(0).dot(u.segment(0, 3));
	p(1) = crjgm.col(1).dot(u.segment(1, 3));
	p(2) = crjgm.col(2).dot(u.segment(2, 3));

	weights = sigm.cwiseProduct(dkmg);
	weights /= weights.sum();
	outglq(1) = p.dot(weights);

	return outglq;
}

// Quadrature points of 4 Gauss Lobatto absicas
// Index 0->i+1/2, 1->i-1/2, 2->i+sqrt(5)/10, 3->i-sqrt(5)/10
Matrix<double, 4, 1> Schemes::Lagrange::weno5JS4PointsGaussLobatto(const Vector5d& u) {
	smoothnessMeasure4PointsGaussLobatto(u);

	p(2) = crj.col(1).dot(u.segment(2, 3));
	p(1) = crj.col(2).dot(u.segment(1, 3));
	p(0) = crj.col(3).dot(u.segment(0, 3));
	weights = dkm.array() / (eps + SILobattoA.array()).pow(2);
	weights = weights / weights.sum();
	outGaussLobatto(0) = p.dot(weights);
	outGaussLobatto(0) = catchnan(outGaussLobatto(0));

	p(2) = crj.col(0).dot(u.segment(2, 3));
	p(1) = crj.col(1).dot(u.segment(1, 3));
	p(0) = crj.col(2).dot(u.segment(0, 3));
	weights = dkp.array() / (eps + SILobattoA.array()).pow(2);
	weights = weights / weights.sum();
	outGaussLobatto(1) = p.dot(weights);
	outGaussLobatto(1) = catchnan(outGaussLobatto(1));

	p(0) = crjgpLobatto.col(0).dot(u.segment(0, 3));
	p(1) = crjgpLobatto.col(1).dot(u.segment(1, 3));
	p(2) = crjgpLobatto.col(2).dot(u.segment(2, 3));

	weights = dkpLobatto.array() / (eps + SILobattoB.array()).pow(2);
	weights = weights.array() / weights.sum();
	outGaussLobatto(2) = p.dot(weights);
	outGaussLobatto(2) = catchnan(outGaussLobatto(2));

	p(0) = crjgmLobatto.col(0).dot(u.segment(0, 3));
	p(1) = crjgmLobatto.col(1).dot(u.segment(1, 3));
	p(2) = crjgmLobatto.col(2).dot(u.segment(2, 3));

	weights = dkmLobatto.array() / (eps + SILobattoB.array()).pow(2);
	weights = weights.array() / weights.sum();
	outGaussLobatto(3) = p.dot(weights);
	outGaussLobatto(3) = catchnan(outGaussLobatto(3));

	return outGaussLobatto;
}

Matrix<double, 4, 1> Schemes::Lagrange::weno5JS4PointsGaussLegendre(const Vector5d& u) {
	smoothnessMeasure4PointsGaussLegendre(u);

	p(0) = crjgpLegendre_12.col(0).dot(u.segment(0, 3));
	p(1) = crjgpLegendre_12.col(1).dot(u.segment(1, 3));
	p(2) = crjgpLegendre_12.col(2).dot(u.segment(2, 3));
	weights = dkpLegendre_12.array() / (eps + SILegendreA.array()).pow(2);
	weights = weights / weights.sum();
	outGaussLegendre(0) = p.dot(weights);
	outGaussLegendre(0) = catchnan(outGaussLegendre(0));

	p(0) = crjgmLegendre_12.col(0).dot(u.segment(0, 3));
	p(1) = crjgmLegendre_12.col(1).dot(u.segment(1, 3));
	p(2) = crjgmLegendre_12.col(2).dot(u.segment(2, 3));
	weights = dkmLegendre_12.array() / (eps + SILegendreA.array()).pow(2);
	weights = weights / weights.sum();
	outGaussLegendre(1) = p.dot(weights);
	outGaussLegendre(1) = catchnan(outGaussLegendre(1));

	p(0) = crjgpLegendre_34.col(0).dot(u.segment(0, 3));
	p(1) = crjgpLegendre_34.col(1).dot(u.segment(1, 3));
	p(2) = crjgpLegendre_34.col(2).dot(u.segment(2, 3));
	weights = dkpLegendre_34.array() / (eps + SILegendreB.array()).pow(2);
	weights = weights.array() / weights.sum();
	outGaussLegendre(2) = p.dot(weights);
	outGaussLegendre(2) = catchnan(outGaussLegendre(2));

	p(0) = crjgmLegendre_34.col(0).dot(u.segment(0, 3));
	p(1) = crjgmLegendre_34.col(1).dot(u.segment(1, 3));
	p(2) = crjgmLegendre_34.col(2).dot(u.segment(2, 3));
	weights = dkmLegendre_34.array() / (eps + SILegendreB.array()).pow(2);
	weights = weights.array() / weights.sum();
	outGaussLegendre(3) = p.dot(weights);
	outGaussLegendre(3) = catchnan(outGaussLegendre(3));

	return outGaussLegendre;
}

Matrix<double, 4, 1> Schemes::Lagrange::weno5Z4PointsGaussLegendre(const Vector5d& u) {
	smoothnessMeasure4PointsGaussLegendre(u);
	tau = fabs(SILegendreA(0) - SILegendreA(2));

	p(0) = crjgpLegendre_12.col(0).dot(u.segment(0, 3));
	p(1) = crjgpLegendre_12.col(1).dot(u.segment(1, 3));
	p(2) = crjgpLegendre_12.col(2).dot(u.segment(2, 3));
	weights = dkpLegendre_12.array() * (1.0 + tau / (SILegendreA.array() + epsz));
	weights = weights / weights.sum();
	outGaussLegendre(0) = p.dot(weights);
	outGaussLegendre(0) = catchnan(outGaussLegendre(0));

	p(0) = crjgmLegendre_12.col(0).dot(u.segment(0, 3));
	p(1) = crjgmLegendre_12.col(1).dot(u.segment(1, 3));
	p(2) = crjgmLegendre_12.col(2).dot(u.segment(2, 3));
	weights = dkmLegendre_12.array() * (1.0 + tau / (SILegendreA.array() + epsz));
	weights = weights / weights.sum();
	outGaussLegendre(1) = p.dot(weights);
	outGaussLegendre(1) = catchnan(outGaussLegendre(1));

	tau = fabs(SILegendreB(0) - SILegendreB(2));

	p(0) = crjgpLegendre_34.col(0).dot(u.segment(0, 3));
	p(1) = crjgpLegendre_34.col(1).dot(u.segment(1, 3));
	p(2) = crjgpLegendre_34.col(2).dot(u.segment(2, 3));
	weights = dkpLegendre_34.array() * (1.0 + tau / (SILegendreB.array() + epsz));
	weights = weights.array() / weights.sum();
	outGaussLegendre(2) = p.dot(weights);
	outGaussLegendre(2) = catchnan(outGaussLegendre(2));

	p(0) = crjgmLegendre_34.col(0).dot(u.segment(0, 3));
	p(1) = crjgmLegendre_34.col(1).dot(u.segment(1, 3));
	p(2) = crjgmLegendre_34.col(2).dot(u.segment(2, 3));
	weights = dkmLegendre_34.array() * (1.0 + tau / (SILegendreB.array() + epsz));
	weights = weights.array() / weights.sum();
	outGaussLegendre(3) = p.dot(weights);
	outGaussLegendre(3) = catchnan(outGaussLegendre(3));

	return outGaussLegendre;
}

// Quadrature points of 4 Gauss Lobatto absicas
// Index 0->i+1/2, 1->i-1/2, 2->i+sqrt(5)/10, 3->i-sqrt(5)/10
Matrix<double, 4, 1> Schemes::Lagrange::weno5Z4PointsGaussLobatto(const Vector5d& u) {
	smoothnessMeasure4PointsGaussLobatto(u);
	tau = fabs(SILobattoA(0) - SILobattoA(2));
	SZG = (1.0 + tau / (SILobattoA.array() + epsz));

	p(2) = crj.col(1).dot(u.segment(2, 3));
	p(1) = crj.col(2).dot(u.segment(1, 3));
	p(0) = crj.col(3).dot(u.segment(0, 3));
	weights = dkm.cwiseProduct(SZG);
	weights = weights / weights.sum();
	outGaussLobatto(0) = p.dot(weights);

	p(2) = crj.col(0).dot(u.segment(2, 3));
	p(1) = crj.col(1).dot(u.segment(1, 3));
	p(0) = crj.col(2).dot(u.segment(0, 3));
	weights = dkp.cwiseProduct(SZG);
	weights = weights / weights.sum();
	outGaussLobatto(1) = p.dot(weights);

	tau = fabs(SILobattoB(0) - SILobattoB(2));
	SZG = (1.0 + tau / (SILobattoB.array() + epsz));

	p(0) = crjgpLobatto.col(0).dot(u.segment(0, 3));
	p(1) = crjgpLobatto.col(1).dot(u.segment(1, 3));
	p(2) = crjgpLobatto.col(2).dot(u.segment(2, 3));

	weights = dkpLobatto.cwiseProduct(SZG);
	weights = weights.array() / weights.sum();
	outGaussLobatto(2) = p.dot(weights);

	p(0) = crjgmLobatto.col(0).dot(u.segment(0, 3));
	p(1) = crjgmLobatto.col(1).dot(u.segment(1, 3));
	p(2) = crjgmLobatto.col(2).dot(u.segment(2, 3));

	weights = dkmLobatto.cwiseProduct(SZG);
	weights = weights.array() / weights.sum();
	outGaussLobatto(3) = p.dot(weights);

	return outGaussLobatto;
}

// Quadrature points of 4 Gauss Lobatto absicas
// Index 0->i+1/2, 1->i-1/2, 2->i+sqrt(5)/10, 3->i-sqrt(5)/10
Matrix<double, 4, 1> Schemes::Lagrange::teno5th4PointsGaussLobatto(const Vector5d& u) {
	smoothnessMeasure4PointsGaussLobatto(u);
	tau = fabs(SILobattoA(0) - SILobattoA(2));
	weights = (1.0 + tau / (SILobattoA.array() + epsz)).pow(6);
	weights = weights / weights.sum();
	cutOff();

	p(2) = crj.col(1).dot(u.segment(2, 3));
	p(1) = crj.col(2).dot(u.segment(1, 3));
	p(0) = crj.col(3).dot(u.segment(0, 3));
	weights = dkm.cwiseProduct(sigm);
	weights = weights / weights.sum();
	outGaussLobatto(0) = p.dot(weights);
	outGaussLobatto(0) = catchnan(outGaussLobatto(0));

	p(2) = crj.col(0).dot(u.segment(2, 3));
	p(1) = crj.col(1).dot(u.segment(1, 3));
	p(0) = crj.col(2).dot(u.segment(0, 3));
	weights = dkp.cwiseProduct(sigm);
	weights = weights / weights.sum();
	outGaussLobatto(1) = p.dot(weights);
	outGaussLobatto(1) = catchnan(outGaussLobatto(1));

	tau = fabs(SILobattoB(0) - SILobattoB(2));
	weights = (1.0 + tau / (SILobattoB.array() + epsz)).pow(6);
	weights = weights / weights.sum();
	cutOff();

	p(0) = crjgpLobatto.col(0).dot(u.segment(0, 3));
	p(1) = crjgpLobatto.col(1).dot(u.segment(1, 3));
	p(2) = crjgpLobatto.col(2).dot(u.segment(2, 3));

	weights = dkpLobatto.cwiseProduct(sigm);
	weights = weights.array() / weights.sum();
	outGaussLobatto(2) = p.dot(weights);
	outGaussLobatto(2) = catchnan(outGaussLobatto(2));

	p(0) = crjgmLobatto.col(0).dot(u.segment(0, 3));
	p(1) = crjgmLobatto.col(1).dot(u.segment(1, 3));
	p(2) = crjgmLobatto.col(2).dot(u.segment(2, 3));

	weights = dkmLobatto.cwiseProduct(sigm);
	weights = weights.array() / weights.sum();
	outGaussLobatto(3) = p.dot(weights);
	outGaussLobatto(3) = catchnan(outGaussLobatto(3));

	return outGaussLobatto;
}

void Schemes::Lagrange::smoothnessMeasure3rd(const Vector3d& u) {
	SI3rd(0) = u(0) - u(1);
	SI3rd(1) = u(1) - u(2);
	SI3rd = SI3rd.array().pow(2);
}

// Smoothness measurements of fifth order WENO-JS scheme
void Schemes::Lagrange::smoothnessMeasure(const Vector5d& u) {
	c(0) = u(2) - 2.0 * u(3) + u(4);
	b(0) = 3.0 * u(2) - 4.0 * u(3) + u(4);

	c(1) = u(1) - 2.0 * u(2) + u(3);
	b(1) = u(1) - u(3);

	c(2) = u(0) - 2.0 * u(1) + u(2);
	b(2) = u(0) - 4.0 * u(1) + 3.0 * u(2);

	b = b / 2.0;

	SI = b.array() * b.array() + c1 * c.array() * c.array();
}

// full 5-points smoothness
void Schemes::Lagrange::fullpointsSmootness(const Vector5d& u) {
	b35 = u(4) * (6908.0 * u(4) - 51001.0 * u(3) + 67923.0 * u(2) - 38947.0 * u(1) + 8209.0 * u(0)) +
		u(3) * (104963.0 * u(3) - 299076.0 * u(2) + 179098.0 * u(1) - 38947.0 * u(0)) + u(2) * (231153.0 * u(2) -
			299076.0 * u(1) + 67923.0 * u(0)) + u(1) * (104963.0 * u(1) - 51001.0 * u(0)) + 6908.0 * u(0) * u(0);
	b35 /= 5040.0;
}

void Schemes::Lagrange::smoothnessMeasureIS(const Vector5d& u) {
	SI4is(0) = pow(u(3) - u(2), 2);
	SI4is(1) = pow(u(2) - u(1), 2);
	SI4is(2) = 0.25 * pow(3.0 * u(2) - 4.0 * u(3) + u(4), 2) + c1 * pow(u(2) - 2.0 * u(3) + u(4), 2);
	SI4is(3) = 0.25 * pow(u(0) - 4.0 * u(1) + 3.0 * u(2), 2) + c1 * pow(u(0) - 2.0 * u(1) + u(2), 2);

	b12 = 0.25 * pow(u(1) - u(3), 2) + c1 * pow(u(1) - 2.0 * u(2) + u(3), 2);
	tau = 0.25 * pow(u(4) - 2.0 * (u(3) - u(1)) - u(0), 2) + \
		c1 * pow(u(4) - 4.0 * u(3) + 6.0 * u(2) - 4.0 * u(1) + u(0), 2);

	weightsis(0) = 1.0 + pow(tau, 2) / (SI4is(0) + eps) / (b12 + eps);
	weightsis(1) = 1.0 + pow(tau, 2) / (SI4is(1) + eps) / (b12 + eps);
	weightsis(2) = 1.0 + tau / (SI4is(2) + eps);
	weightsis(3) = 1.0 + tau / (SI4is(3) + eps);
}

void Schemes::Lagrange::smoothnessMeasureISmod(const Vector5d& u, const double& diff) {
	SI4is(0) = pow(u(3) - u(2), 2);
	SI4is(1) = pow(u(2) - u(1), 2);
	SI4is(2) = 0.25 * pow(3.0 * u(2) - 4.0 * u(3) + u(4), 2) + c1 * pow(u(2) - 2.0 * u(3) + u(4), 2);
	SI4is(3) = 0.25 * pow(u(0) - 4.0 * u(1) + 3.0 * u(2), 2) + c1 * pow(u(0) - 2.0 * u(1) + u(2), 2);

	b12 = 0.25 * pow(u(1) - u(3), 2) + c1 * pow(u(1) - 2.0 * u(2) + u(3), 2);
	tau = 0.25 * pow(u(4) - 2.0 * (u(3) - u(1)) - u(0), 2) +
		c1 * pow(u(4) - 4.0 * u(3) + 6.0 * u(2) - 4.0 * u(1) + u(0), 2);

	weightsis(0) = 1.0 + tau * tau / (SI4is(0) * b12 + diff);
	weightsis(1) = 1.0 + tau * tau / (SI4is(1) * b12 + diff);
	weightsis(2) = 1.0 + tau / (SI4is(2) + diff);
	weightsis(3) = 1.0 + tau / (SI4is(3) + diff);
}

void Schemes::Lagrange::smoothnessMeasureISExtended(const Vector5d& u) {
	c(2) = u(2) - 2.0 * u(3) + u(4);
	b(2) = 3.0 * u(2) - 4.0 * u(3) + u(4);
	c(1) = u(1) - 2.0 * u(2) + u(3);
	b(1) = u(1) - u(3);
	c(0) = u(0) - 2.0 * u(1) + u(2);
	b(0) = u(0) - 4.0 * u(1) + 3.0 * u(2);
	b = b / 2.0;
	SI = b.array() * b.array() + c1 * c.array() * c.array();

	for (int i = 0; i < 4; ++i) {
		SI4is(i) = pow(u(i + 1) - u(i), 2);
	}
	tau = 0.25 * pow(u(4) - 2.0 * (u(3) - u(1)) - u(0), 2) + \
		c1 * pow(u(4) - 4.0 * u(3) + 6.0 * u(2) - 4.0 * u(1) + u(0), 2);
}

// Smoothness measurement of 2 points gaussian integration
void Schemes::Lagrange::smoothnessMeasure2PointsGaussLegendre(const Vector5d& u) {
	bg(0) = u(0) - 4.0 * u(1) + 3.0 * u(2);
	cg(0) = u(0) - 2.0 * u(1) + u(2);

	bg(1) = u(1) - u(3);
	cg(1) = u(1) - 2.0 * u(2) + u(3);

	bg(2) = 3.0 * u(2) - 4.0 * u(3) + u(4);
	cg(2) = u(2) - 2.0 * u(3) + u(4);

	SIGLQ = c2 * bg.array().pow(2) + c3 * cg.array().pow(2);
}

// Smoothness measurements of four points Gauss-Lobatto quadrature
// A for x_{i+1/2} x_{i-1/2} and B for x_{i+sqrt(5)/10} x{i-sqrt(5)/10}
void Schemes::Lagrange::smoothnessMeasure4PointsGaussLobatto(const Vector5d& u) {
	bg(0) = u(0) - 4.0 * u(1) + 3.0 * u(2);
	cg(0) = u(0) - 2.0 * u(1) + u(2);

	bg(1) = u(1) - u(3);
	cg(1) = u(1) - 2.0 * u(2) + u(3);

	bg(2) = 3.0 * u(2) - 4.0 * u(3) + u(4);
	cg(2) = u(2) - 2.0 * u(3) + u(4);

	bg = bg.array().pow(2);
	cg = cg.array().pow(2);

	SILobattoA = 0.25 * bg + c3 * cg;
	SILobattoB = c4 * bg + c5 * cg;
}

void Schemes::Lagrange::smoothnessMeasure4PointsGaussLegendre(const Vector5d& u) {
	bg(0) = u(0) - 4.0 * u(1) + 3.0 * u(2);
	cg(0) = u(0) - 2.0 * u(1) + u(2);

	bg(1) = u(1) - u(3);
	cg(1) = u(1) - 2.0 * u(2) + u(3);

	bg(2) = 3.0 * u(2) - 4.0 * u(3) + u(4);
	cg(2) = u(2) - 2.0 * u(3) + u(4);

	bg = bg.array().pow(2);
	cg = cg.array().pow(2);

	SILegendreA = c7 * bg + c8 * cg;
	SILegendreB = c9 * bg + c10 * cg;
}

void Schemes::Lagrange::smoothnessMeasure2PointsGaussLegendre3th(const Vector3d& u) {
	SI3rd(0) = u(0) - u(1);
	SI3rd(1) = u(1) - u(2);
	SI3rd = c6 * SI3rd.array().pow(2);
}

void Schemes::Lagrange::cutOff() {
	int i;

    weights /= weights.sum();

	for (i = 0; i < 3; i++) {
		sigm(i) = weights(i) < ct ? 0.0 : 1.0;
	}
}

double Schemes::Lagrange::catchnan(const double& input) {
	return isnan(input) ? machine_eps : input;
}

void Schemes::Lagrange::smoothnessMeasure3L(const Vector5d &u) {
    bg(0) = u(0) - 4.0 * u(1) + 3.0 * u(2);
    cg(0) = u(0) - 2.0 * u(1) + u(2);

    bg(1) = u(1) - u(3);
    cg(1) = u(1) - 2.0 * u(2) + u(3);

    bg(2) = 3.0 * u(2) - 4.0 * u(3) + u(4);
    cg(2) = u(2) - 2.0 * u(3) + u(4);

    bg = 0.25 * bg.array().pow(2);
    cg = c1 * cg.array().pow(2);

    SI3L = bg + cg;
}

Vector3d Schemes::Lagrange::weno5JS3PointsQuadrature(const Vector5d &u) {
    smoothnessMeasure3L(u);
//    tau = fabs(SI3L(0) - SI3L(2));
//    Vector3d wg;
//    weights = (1.0 + (tau / (SI3L.array() + epsz))).pow(6);
//    cutOff();

    p(0) = crjgp3L.col(0).dot(u.segment(0, 3));
    p(1) = crjgp3L.col(1).dot(u.segment(1, 3));
    p(2) = crjgp3L.col(2).dot(u.segment(2, 3));

    weights = linwp_3L.array() / (eps + SI3L.array()).pow(2);
//    weights = linwp_3L.array() * wg.array();
//    weights = linwp_3L.cwiseProduct(sigm);
    weights = weights.array() / weights.sum();
    outGausLegendre3n(0) = p.dot(weights);

    p(0) = (-u(0) + 2*u(1) + 23.0*u(2))/24.0;
    p(1) = (-u(1) + 26.0*u(2) - u(3))/24.0;
    p(2) = (23.0*u(2) + 2.0*u(3) - u(4))/24.0;
    weights = linw0p_3L.array() / (eps + SI3L.array()).pow(2);
//    weights = linw0p_3L.array() * wg.array();
//    weights = linw0p_3L.cwiseProduct(sigm);
    weights = weights.array() / weights.sum();
    outGausLegendre3n(1) = sumw0p_3L * p.dot(weights);
    weights = linw0m_3L.array() / (eps + SI3L.array()).pow(2);
//    weights = linw0m_3L.array() * wg.array();
//    weights = linw0m_3L.cwiseProduct(sigm);
    weights = weights.array() / weights.sum();
    outGausLegendre3n(1) -= sumw0m_3L * p.dot(weights);

    p(0) = crjgm3L.col(0).dot(u.segment(0, 3));
    p(1) = crjgm3L.col(1).dot(u.segment(1, 3));
    p(2) = crjgm3L.col(2).dot(u.segment(2, 3));

    weights = linwm_3L.array() / (eps + SI3L.array()).pow(2);
//    weights = linwm_3L.array() * wg.array();
//    weights = linwm_3L.cwiseProduct(sigm);
    weights = weights.array() / weights.sum();
    outGausLegendre3n(2) = p.dot(weights);

    for(int i = 0; i < 3; ++i) {
        outGausLegendre3n(i) = catchnan(outGausLegendre3n(i));
    }

    return outGausLegendre3n;
}

Vector2d Schemes::Lagrange::linear(const Vector6d &u) {
    p(0) = crj.col(1).dot(u.segment(2, 3));
    p(1) = crj.col(2).dot(u.segment(1, 3));
    p(2) = crj.col(3).dot(u.segment(0, 3));

    out(0) = p.dot(dkm);

    p(0) = crj.col(0).dot(u.segment(3, 3));
    p(1) = crj.col(1).dot(u.segment(2, 3));
    p(2) = crj.col(2).dot(u.segment(1, 3));

    out(1) = p.dot(dkp);

    return out;
}

Vector2d Schemes::Lagrange::teno5i(const Vector6d &u) {
    p(0) = crj.col(1).dot(u.segment(2, 3));
    p(1) = crj.col(2).dot(u.segment(1, 3));
    p(2) = crj.col(3).dot(u.segment(0, 3));

    smoothnessMeasure(u.segment(0, 5));
    tau = SI.maxCoeff();
    weights = (1.0 + tau / (SI.array() + epsz)).pow(6);
    cutOff();
    weights = sigm.cwiseProduct(dkm);
    weights /= weights.sum();
    out(0) = p.dot(weights);
    out(0) = catchnan(out(0));

    p(0) = crj.col(0).dot(u.segment(3, 3));
    p(1) = crj.col(1).dot(u.segment(2, 3));
    p(2) = crj.col(2).dot(u.segment(1, 3));

    smoothnessMeasure(u.segment(1, 5));
    tau = SI.maxCoeff();
    weights = (1.0 + tau / (SI.array() + epsz)).pow(6);
    cutOff();
    weights = sigm.cwiseProduct(dkp);
    weights /= weights.sum();
    out(1) = p.dot(weights);
    out(1) = catchnan(out(1));

    return out;
}

Vector2d Schemes::Lagrange::teno5a(const Vector6d &u) {
    p(0) = crj.col(1).dot(u.segment(2, 3));
    p(1) = crj.col(2).dot(u.segment(1, 3));
    p(2) = crj.col(3).dot(u.segment(0, 3));

    smoothnessMeasure(u.segment(0, 5));
    tau = 1.0;
    weights = (tau / (SI.array() + epsz)).pow(6);
    cutOff();
    weights = sigm.cwiseProduct(dkm);
    weights /= weights.sum();
    out(0) = p.dot(weights);
    out(0) = catchnan(out(0));

    p(0) = crj.col(0).dot(u.segment(3, 3));
    p(1) = crj.col(1).dot(u.segment(2, 3));
    p(2) = crj.col(2).dot(u.segment(1, 3));

    smoothnessMeasure(u.segment(1, 5));
    tau = 1.0;
    weights = (tau / (SI.array() + epsz)).pow(6);
    cutOff();
    weights = sigm.cwiseProduct(dkp);
    weights /= weights.sum();
    out(1) = p.dot(weights);
    out(1) = catchnan(out(1));

    return out;
}

// Deconstructor
Schemes::Lagrange::~Lagrange() = default;
