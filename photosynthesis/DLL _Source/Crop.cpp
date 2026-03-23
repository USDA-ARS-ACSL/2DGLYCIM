// 5/06/2020 WSun designed a gas exchange module using dll which connect the original GLYCIM
// 5/08/2020 WSun modified the program by using transpiration of photosynthesis module to calculate EO
// 5/09/2020 WSun deleted the ET_supply which not used in the photosynthesis
// 5/10/2020 WSun deleted the parameters (Tubmod, Feedback, leafageeffect) which was used in SPUDSIM (might added in future)
// 5/11/2020 WSun give a reasonable longitude and reset conductance
// 5/17/2020 WSun applied Ball-Berry model stomatal conductance values (g0 and g1)
// 5/18/2020 WSun revised the code in PNET rountines which calculate the USEDCS 

// crop.cpp : Defines the entry point for the DLL application
//#define MYPLANT_EXPORTS

#include "pch.h"
#include "crop.h"
#include "controller.h"
#include "weather.h"

#include <iostream>
//using namespace std;
#include <cmath>
#define endl "\n"
#define comma ","

// note that we have to dereference the variable in order to
// assign a value that can be passed back to GLYCIM. This is 
// because the FORTRAN program expects a pointer rather than
// a value. I don't think this applies to structures as it does 
// variables that may be in the arguments list.
// note use of lower case names. Upper and lower case conversions between
// fortran and C++ don't matter here because these are arguments and not
// a function name. GASEXCHANGER must be upper case because it is a function name

int compare(const void* arg1, const void* arg2)
{
	/* Compare all of both strings: */
	if (*(double*)arg1 > *(double*)arg2) return 1;
	else if (*(double*)arg1 < *(double*)arg2) return -1;
	return 0;
};
#ifdef _WIN32
void _stdcall GASEXCHANGER(struct
#else
void gasexchanger_(struct
#endif
	WeatherCommon* Weather,   // Define global variables block which connect GLYCIM variables
	PlantCommon* Plant

)

{


	static CController* pSC; //declare as static to ensure only one copy is instantiated during GLYCIM execution

//	TInitInfo initInfo;
	{

		initInfo.latitude = Weather->CLATUDE; 
		//initInfo.longitude = -13.0; // hard coded. Need to read in GLYCIM
		initInfo.longitude = -93.75;
		initInfo.altitude = 480;      // hard coded. Need to read in GLYCIM

		initInfo.NRATIO = Plant->NRATIO;
		//Weather->PSIL[Weather->ITIME-1] = -2.0;
		if (getInitialized() == 0)
		{
          pSC = new CController(initInfo); // RESET variables at the beginning of each run (1 hour)
		  setInitialized(1);
        }
		
	}
	

	 double Es;

	TWeather wthr;
	{
		wthr.DayOfYear = Weather->CDayOfYear;
		wthr.time = ((Weather->CITIME) - 1) * (1.0 / (Weather->CIPERD));  //not sure this will be used
		//wthr.time = ((Weather->ITIME) - 1) * 0.0416667;
		wthr.CO2 = Weather->CCO2;
		//wthr.airT = Weather->TAIR[Weather->ITIME];
		wthr.airT = Weather->CTAIR[Weather->CITIME-1];
		wthr.PFD = Weather->Cpar[Weather->CITIME-1] * 4.55; // conversion from PAR in W m-2 to umol s-1 m-2
		wthr.solRad = Weather->CWATTSM[Weather->CITIME-1]; //conversion from W m-2 total radiation to J m-2 in one hour 				
		double Es = (0.611 * exp(17.502 * wthr.airT / (240.97 + wthr.airT))); // saturated vapor pressure at airT
		wthr.RH = (1 - (Weather->CVPD[Weather->CITIME-1] / Es)) * 100.0; // relative humidity in percent
		wthr.wind = Weather->CWIND * (1000.0 / 3600.0); // conversion from km hr-1 to m s-1
		wthr.psil_ = Weather->CPSIL_/ 10;// convert current leaf water potential from bars to MPa
	    // wthr.psil_ = -0.2;
		wthr.LAI = Weather->CLAI;
		//wthr.LAI = Weather->LAREAT * 38.75196 * 1.0E-4; // LAI = LAREAT * POPARE * 1.0E-4
		//wthr.ET_supply = 0.0;

		/******************************/ 

	}


       // A gas exchange module object is initialized (calls function) here
       //  ***************************************************************************
 		

		pSC->run(wthr, initInfo);

		// Plant->photosynthesis_gross = pSC->get_photosynthesis_gross();//umol CO2 m-2 GROUND s-1 (already multiply LAI) 
		//I need to transfer this to mg CO2 m-2 GROUND s-1 
		
		Plant->photosynthesis_gross = pSC->get_photosynthesis_gross() * (44.0/1000.0);
		//std::cout << "Time= ";//<< wthr.time << "photosynthesis_gross = " << Plant->photosynthesis_gross << endl;
		
		//Plant->photosynthesis_gross = pSC->get_photosynthesis_gross()/(wthr.LAI)* (44.0 / 1000.0);

		//Plant->photosynthesis_net = pSC->get_photosynthesis_net()/(wthr.LAI)*(44.0/1000.0);
		
		Plant->photosynthesis_net = pSC->get_photosynthesis_net() * (44.0 / 1000.0);
		
		
		Plant->transpiration = pSC->get_transpiration();
		Plant->temperature = pSC->get_temperature(); // return back leaf temperature
		Plant->TLAI = pSC->get_TLAI(); 
		
		Plant->photosynthesis_netsunlitleaf = pSC->get_photosynthesis_netsunlitleaf()*Plant->sunlitLAI;
		Plant->photosynthesis_netshadedleaf = pSC->get_photosynthesis_netshadedleaf()*Plant->shadedLAI;

		Plant->sunlitLAI = pSC->get_sunlitLAI();
		Plant->shadedLAI = pSC->get_shadedLAI();

		Plant->transpiration_sunlitleaf = pSC->get_transpiration_sunlitleaf();
		Plant->transpiration_shadedleaf = pSC->get_transpiration_shadedleaf();
		Plant->temp1 = pSC->get_temp1();
		Plant->Ags = pSC->get_Ags();
		Plant->ARH = wthr.RH;
		return;
		
}