// OE 01-2018
#if !defined (_initExtractsub_h_)
#define _initExtractsub_h_
#include "configFile.h"
#include "wghmStateFile.h"

extern"C"
{
  void init_extractsub_ (const char* string_ids,const char* string_config, WghmStateFile *& wghmState, double *& output, const char* str_calPar_file,long *oy,calibParamClass *& calParam, long *total_nr_calPar,long *calpar_size,int *calPar_index,const char* calPar_filename,int *groupmatrixindex,const char* pathTowghmMean,WghmStateFile *& wghmMean);
}

#endif
