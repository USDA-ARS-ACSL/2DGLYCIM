// The following ifdef block is the standard way of creating macros which make exporting 
// from a DLL simpler. All files within this DLL are compiled with the PLANT_EXPORTS
// symbol defined on the command line. this symbol should not be defined on any project
// that uses this DLL. This way any other project whose source files include this file see 
// PLANT_API functions as being imported from a DLL, whereas this DLL sees symbols
// defined with this macro as being exported.

// See creating a pocket PC dll using C++ article for information
// on how to do this
#ifdef _WIN32
#ifdef MYPLANT_EXPORTS
#define PLANT_API __declspec(dllexport)
#else
#define PLANT_API __declspec(dllexport)
#endif
#else
#define PLANT_API
#endif
#include "initInfo.h"




TInitInfo	initInfo;

//Common Structures defined here
	
#pragma pack(2)
struct WeatherCommon {

	int   CDayOfYear, CITIME, CIPERD;
	float CWATTSM[24], Cpar[24], CTAIR[24], CCO2, CVPD[24], CWIND, CPSIL_, CLATUDE, CLAREAT,CLAI;

};

struct PlantCommon {

	float NRATIO, photosynthesis_gross, photosynthesis_net, transpiration, temperature, TLAI, sunlitLAI, shadedLAI, LightIC, 
		transpiration_sunlitleaf, transpiration_shadedleaf,temp1,Ags,ARH, photosynthesis_netsunlitleaf, photosynthesis_netshadedleaf;

};



int Initialized=0;



void setInitialized(int x) { Initialized = x; }
int getInitialized() { return Initialized; }
TInitInfo getInitInfo() { return initInfo; }



#pragma pack()


 
#ifdef __cplusplus
extern "C" {
#endif

// Your exported function headers go here 
// GASEXCHANGER must be upper and lower case because it is a function name
#ifdef _WIN32
	PLANT_API void _stdcall GASEXCHANGER(struct WeatherCommon    *, PlantCommon    *);
#else
	PLANT_API void gasexchanger_(struct WeatherCommon    *, PlantCommon    *);,
#endif
#ifdef __cplusplus
}
#endif


 

