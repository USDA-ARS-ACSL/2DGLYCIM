
#pragma once
#include "pch.h"
#include "Controller.h"
#include "initinfo.h"
#include "radtrans.h"
#include "solar.h"
#include "math.h"
#include <vector>
using namespace std;



CController::CController()
{

}

CController::CController(const TInitInfo info)
{

	initInfo = info;
	sunlitLAI = shadedLAI = 0.0;
	sunlitPFD = shadedPFD = 0.0;
	photosynthesis_net = photosynthesis_gross = transpiration = 0.0;
	VPD=conductance = 0.0;
	//initInfo.longitude = 20.0;
	//initInfo.altitude = 50.0;
	temperature = 20; //no meaning for now
	transpiration_sunlitleaf = transpiration_shadedleaf = 0.0;
	photosynthesis_netsunlitleaf = photosynthesis_netshadedleaf = 0.0;

	temp1 = 0.0;
	Ags = 0.0; //
}
CController::~CController()
{

}

int CController::run(const TWeather & weather, const TInitInfo info)

{
	//int redo = 1; 

	const double tau = 0.50; // atmospheric transmittance, to be implemented as a variable => done
	const double LAF = 0.81; // leaf angle factor for soybean leaves, Campbell and Norman (1998)(p.253)
	//Make leaf width a function of growing leaves as average width will increase as the plant grows
	// add it as a 
	const double leafwidth = 5.0; //to be calculated when implemented for individal leaves
	const double atmPressure = 101.3; //kPa, to be predicted using altitude
	//double poolmod = 1; //value to reduce photosynthetic parameters based on if demand < photosynthesisx0.8
	//double PhotosyntheticFeedback = 1;
	//double leafageEffect = 1;
	double LAI = weather.LAI;
	
	
	CGas_exchange_new * sunlit = new CGas_exchange_new ();
	CGas_exchange_new * shaded = new CGas_exchange_new ();

    CSolar* sun = new CSolar();
	CRadTrans* light = new CRadTrans();
   
	sun->SetVal(weather.DayOfYear, weather.time, initInfo.latitude, initInfo.longitude, initInfo.altitude, weather.solRad);
	light->SetVal(*sun, LAI, LAF);

	sunlitPFD = light->Qsl();
	shadedPFD = light->Qsh();
	sunlitLAI = light->LAIsl();
	shadedLAI = light->LAIsh();
	
	//sunlit->SetVal(sunlitPFD, info, weather.airT, weather.CO2, weather.RH, weather.wind, atmPressure, leafwidth, poolmod, PhotosyntheticFeedback, leafageEffect, weather.psil_, weather.ET_supply);
	//shaded->SetVal(shadedPFD, info, weather.airT, weather.CO2, weather.RH, weather.wind, atmPressure, leafwidth, poolmod, PhotosyntheticFeedback, leafageEffect, weather.psil_, weather.ET_supply);
	sunlit->SetVal(sunlitPFD, info, weather.airT, weather.CO2, weather.RH, weather.wind, atmPressure, leafwidth, weather.psil_);
	shaded->SetVal(shadedPFD, info, weather.airT, weather.CO2, weather.RH, weather.wind, atmPressure, leafwidth, weather.psil_);




	photosynthesis_gross = (sunlit->A_gross * sunlitLAI + shaded->A_gross * shadedLAI); 
		// photosynthesis_gross unit : umol co2 m-2 s-1

	photosynthesis_net = (sunlit->A_net * sunlitLAI + shaded->A_net * shadedLAI);
	transpiration = (sunlit->ET * sunlitLAI + shaded->ET * shadedLAI);
	//transpiration = (sunlit->ET * sunlitLAI + shaded->ET * shadedLAI)*LAI;
	//temperature = (sunlit->Tleaf * sunlitLAI + shaded->Tleaf * shadedLAI) / LAI;
	
	if (LAI != 0)
	{
		temperature = (sunlit->Tleaf * sunlitLAI + shaded->Tleaf * shadedLAI) / LAI;
		this->conductance = __max(0, ((sunlit->get_gs() * sunlitLAI + shaded->get_gs() * shadedLAI) / LAI));//average stomatal conductance
		//this->conductance = (sunlit->get_gs() * sunlitLAI() + shaded->get_gs() * shadedLAI())/LAI);//average stomatal conductance
		this->internal_CO2 = __max(0, ((sunlit->get_ci() * sunlitLAI + shaded->get_ci() * shadedLAI) / LAI)); //average internal CO2 concentrtion
	}
	else 
	{
		temperature = sunlit->Tleaf;
		this->conductance = 0;
	}
	
	set_sunPg(sunlit->A_gross);//check balance
	set_shadePg(shaded->A_gross);

	photosynthesis_netsunlitleaf = sunlit->A_net;
	photosynthesis_netshadedleaf = shaded->A_net;
	photosynthesis_grosssunlitleaf = sunlit->A_gross;
	photosynthesis_grossshadedleaf = shaded->A_gross;
	//psistress_gs_factor = sunlit->get_psileaf_Pn_stress();

	transpiration_sunlitleaf = sunlit->ET;
	transpiration_shadedleaf = shaded->ET;
	TLAI = sunlitLAI + shadedLAI;
	temp1 = sunlit->get_psileaf_Pn_stress();
	Ags = this->conductance;
	//Ags = sunlit->get_gs();
	//Ags = shaded->get_gs();
	delete sunlit;
	delete shaded;
	delete sun;
	delete light;


	return 0;
}