// OE 01-2018
//OE 02-2020 (added wghm temporal meanfield) -> path_wghmMean
#if !defined (_initializeWGHM_h_)
#define _initializeWGHM_h_
#include "configFile.h"
#include "wghmStateFile.h"
#include "calib_param.h"
#include "additionalOutputInputFile.h"
#include "snowInElevationFile.h"

extern"C"
{
	// wird vom fortran aufgerufen
  void initialize_wghm_ (const char* s, WghmStateFile *& initstate, calibParamClass *& initcal, AdditionalOutputInputFile *& initaddio, SnowInElevationFile *& initsnow, long *year, long *month,const char* s2, const char* s3, WghmStateFile *& wghmMean); //int nchar
}
// alias/wrapper zur funk_
// wird von cpp aufgerufen
void initialize_wghm(const char* s, WghmStateFile *& initstate, calibParamClass *& initcal, AdditionalOutputInputFile *& initaddio, SnowInElevationFile *& initsnow, long *year, long *month,const char* s2, const char* s3, WghmStateFile *& wghmMean);
//{
//
//
//	initialize_wghm_(s, initstate,initcal, initaddio, initsnow, year, month); //keine funkition wird hier definiert, sondern nur aufgerufen
//}

#endif
