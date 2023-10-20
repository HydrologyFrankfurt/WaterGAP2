#if !defined (_extractsub_h_)
#define _extractsub_h_
//OE-Jan2018
#include "wghmStateFile.h"
extern"C"
{
  void extract_sub_(const char* s, WghmStateFile *& wghmState, double *& output,long *oy,calibParamClass *& calParam, long *total_nr_calPar, long *calpar_size,int *calPar_index,const char* calPar_filename, WghmStateFile *& wghmMean,int *groupmatrixindex);
}

// 2 remove
void extract_sub(const char* s, WghmStateFile *& wghmState, double *& output,long *oy,calibParamClass *& calParam, long *total_nr_calPar, long *calpar_size,int *calPar_index,const char* calPar_filename, WghmStateFile *& wghmMean,int *groupmatrixindex);

#endif
