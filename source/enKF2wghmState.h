#if !defined (_enKF2wghmState_h_)
#define _enKF2wghmState_h_
//OE-Jan2018
#include "configFile.h"
#include "wghmStateFile.h"
#include "additionalOutputInputFile.h"
#include "snowInElevationFile.h"
#include "parameterJsonFile.h"
extern"C"
{
  void enkf_wghmstate_(const char* s,double *field,double *prediction, ConfigFile *& configFile, WghmStateFile *& wghmState, AdditionalOutputInputFile *& additionalOutIn, SnowInElevationFile *& snow_in_elevation,long *step, long *total_steps,long *year, long *month,long *ny,double *& output,WghmStateFile *& wghmStateMean, const char* s2,calibParamClass *& calParam, long *calpar_size,const char* calparsample, const char* path_calpar_IDs,const char* calPar_filename, WghmStateFile *& wghmMean,double *calpar_range,long *oy, long *total_nr_calPar,int *calPar_index,int *groupmatrixindex);
}

// 2 remove
  void enkf_wghmstate(const char* s,double *field,double *prediction, ConfigFile *& configFile, WghmStateFile *& wghmState, AdditionalOutputInputFile *& additionalOutIn, SnowInElevationFile *& snow_in_elevation,long *step, long *total_steps,long *year, long *month,long *ny,double *& output,WghmStateFile *& wghmStateMean, const char* s2,calibParamClass *& calParam, long *calpar_size,const char* calparsample, const char* path_calpar_IDs,const char* calPar_filename, WghmStateFile *& wghmMean,double *calpar_range,long *oy, long *total_nr_calPar,int *calPar_index,int *groupmatrixindex);
#endif
