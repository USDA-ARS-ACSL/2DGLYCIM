#pragma once
#include "initInfo.h"
double QuadSolnUpper (double a, double b, double c );
double QuadSolnLower (double a, double b, double c );

class CGas_exchange_new
{
public:
	CGas_exchange_new(void);
	~CGas_exchange_new(void);

private:
	TInitInfo initInfo;
	double PhotoFluxDensity, R_abs, Tair, CO2, RH, wind, age, SLA, width, Press, N;
	double psileaf; //leaf water potential, MPa
	double psileaf_stress; //0 to 1 factor for stomatal closure effect of leaf water potential
	double leaf_age;
	//void GasEx_psil(double psileaf, double etsupply, const TInitInfo info);
	void GasEx_psil(double psileaf, const TInitInfo info);
	void Photosynthesis(double Ci, double psileaf, const TInitInfo);
	//void EnergyBalance(double pressure);
	void EnergyBalance();
	//void getParms(double Tubmod, double Feedback);
	void getParms();
	double SearchCi(double CO2i, double psileaf, const TInitInfo info);
	double EvalCi(double Ci, double psileaf, const TInitInfo info);

	double gsw(double pressure, const TInitInfo info);
	double gbw();
	double Es(double T);
	double Slope(double T);
	double Rd();
	double set_PSIleafeffect(double pressure);
	double Ci_Ca;  //!< ratio of internal to external CO2, unitless
	double errTolerance; /*!< error tolerance for iterations */
	double eqlTolerance; /*!< equality tolerance */
	int iter_total;      //!< holds total number of iterations */
	int iter1, iter2;         //!< holds iteration counters
	int  iter_Ci;   /*!< iteration value for Ci */
	bool isCiConverged; /*!< true if Ci iterations have converged */

public:
	//void SetVal(double PhotoFluxDensity, const TInitInfo info, double Tair, double CO2, double RH, double wind, double Press, double width, double Tubmod, double Feedback, double leafage, double psil, double etsupply);
	void SetVal(double PhotoFluxDensity, const TInitInfo info, double Tair, double CO2, double RH, double wind, double Press, double width, double psil);
	double get_VPD() {return VPD;}
	double get_gs() {return gs;}
	double get_ci() {return Ci;}
	double get_psileaf_Pn_stress() {return psileaf_stress;}
	
	struct tparms
	{
		  double Eav, TPU25, Eap, Eaj, Jm25, Vcm25, Rd25, Sj, Hj, Sv, Hv, g0, g1, g2, d1, d2, beta_ABA, delta, a_ABA, lambda_r, lambda_l, K_max;
	} Parameters;
	enum gsModel {BBW, L, H};
	enum CalMethod {Stepwise, Simultaneous};
	double A_gross, A_net, ET, Tleaf, Ci, gs, gb, Rdc, VPD,temp;
  
};


	
	

 