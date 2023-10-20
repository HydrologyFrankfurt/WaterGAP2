// OE 01-2018 (based on AEMS lines 157-190 in watergap.cpp)
//OE 02-2020 (added wghm temporal meanfield) -> path_wghmMean
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include "def.h"

// AEMS: New include files
#include "common.h"
#include "initializeWGHM.h"
#include "configFile.h"
#include "wghmStateFile.h"
#include "calib_param.h"
#include "additionalOutputInputFile.h"
#include "snowInElevationFile.h"
#include "globals.h"


void initialize_wghm_ (const char* s,WghmStateFile *& initstate, calibParamClass *& initcal, AdditionalOutputInputFile *& initaddio, SnowInElevationFile *& initsnow, long *year, long *month,const char* s2,const char* s3, WghmStateFile *& wghmMean) // OE
{
    try{
        std::string fileName(s);std::string progName(s2);std::string pathTowghmMean(s3);
        // for measuring  execution time of watergap
        std::cout << "Running watergap with config file <" << fileName << ">" << std::endl;
        std::cout << "Using watergap temporal mean field <" << pathTowghmMean << ">" << std::endl;
        int yr;int mo; // OE
        yr=*year;mo=*month; // OE
        
        ConfigFile configFile(fileName,yr,mo,progName);  // OE: added yr,mo
        
        initstate=new WghmStateFile(ng,1);
        
        if(!configFile.startvaluefile.empty())
    {
        std::cout << "load start values from file <"<<configFile.startvaluefile<<">"<<std::endl;
        
        initstate->load(configFile.startvaluefile);
    }
    
    std::cout <<"initstate "<< initstate->cell(0).tws(0) << std::endl;
    
    
    wghmMean=new WghmStateFile(ng,1);
    std::cout << "load start values from file <"<<pathTowghmMean<<">"<<std::endl;
    wghmMean->load(pathTowghmMean);
    std::cout << "wghmMean ist eingelesen "<<std::endl;
    
    initcal=new calibParamClass;
    initcal->defineVNamesCalibParamJson();
    if(!configFile.parameterfile.empty())
    {
        std::cout << "load calibration parameter values from file <"<<configFile.parameterfile<<">"<<std::endl;
        initcal->readJson(configFile.parameterfile);
    }
   
    // AEMS: load additional file which contains G_days_since_start, G_growingStatus ,G_PrecSum
    // AdditionalOutputInputFile additionalOutIn;
    initaddio=new AdditionalOutputInputFile;
    if(!configFile.additionalfile.empty())
    {
      std::cout << "load G_days_since_start, G_growingStatus and G_PrecSum start values from file <"<<configFile.additionalfile<<">"<<std::endl;
      initaddio->load(configFile.additionalfile);
      initaddio->additionalfilestatus = 1;
      }
      else {
    initaddio->additionalfilestatus = 0;
      }
    
    // AEMS: load snow in elevation start values
    //SnowInElevationFile snow_in_elevation;
    initsnow=new SnowInElevationFile;
    if(!configFile.snowInElevationfile.empty())
    {
      std::cout << "load snow in elevation start values from file <"<<configFile.snowInElevationfile<<">"<<std::endl;
      initsnow->load(configFile.snowInElevationfile);
    }

 }
  catch(std::exception &e)
  {
    throw(Exception(std::string("In watergap::main():\n")+e.what()));
  }

    return;
}
void initialize_wghm(const char* s, WghmStateFile *& initstate, calibParamClass *& initcal, AdditionalOutputInputFile *& initaddio, SnowInElevationFile *& initsnow, long *year, long *month,const char* s2,const char* s3, WghmStateFile *& wghmMean){
	
	initialize_wghm_(s, initstate,initcal, initaddio, initsnow, year, month, s2, s3, wghmMean); //keine funktion wird hier definiert, sondern nur aufgerufen
}




