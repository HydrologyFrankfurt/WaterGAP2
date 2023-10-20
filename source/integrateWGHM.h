//OE-Jan2018
#if !defined (_integrateWGHM_h_)
#define _integrateWGHM_h_  
#include "configFile.h"
#include "wghmStateFile.h"
#include "calib_param.h"
#include "additionalOutputInputFile.h"
#include "snowInElevationFile.h"
extern"C"
{
	// wird vom fortran aufgerufen
  void integrate_wghm_ (const char* s, ConfigFile *& configFile,WghmStateFile *& wghmState, calibParamClass *& calParam, AdditionalOutputInputFile *& additionalOutIn, SnowInElevationFile *& snow_in_elevation,long *step, long *total_steps,long *year, long *month,const char* s2);
}

// alias/wrapper zur funk_
// wird von cpp aufgerufen
void integrate_wghm(const char* s, ConfigFile *& configFile,WghmStateFile *& wghmState, calibParamClass *& calParam, AdditionalOutputInputFile *& additionalOutIn, SnowInElevationFile *& snow_in_elevation,long *step, long *total_steps,long *year, long *month,const char* s2);

#endif
