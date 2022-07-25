/***********************************************************************
*
*climateYear.cpp and .h introduced by Heike Hoffmann-Dobrev (HHD) to read in
*.365 climate input. A poster is showing the effect of daily precipitation
*input instead of monthly on modeled river discarge: Müller Schmied, H., H.
*Hoffmann-Dobrev & P. Döll (2011): Effect of temporal precipitation
*disaggregation for modeled discharge using a global hydrologial model. EGU
*general assembly, 4 - 8. April 2011, Vienna, Austria.
*http://presentations.copernicus.org/EGU2011-1949_presentation.pdf
*
***********************************************************************/
#include "gridio.h"
#include "option.h"
#include "climate.h"
#include "def.h"

#include "climateYear.h"

using namespace std;

extern climateClass climate;
extern gridioClass gridIO;
extern optionClass options;

climateYearClass::climateYearClass(void)
{
    //nothiung to do
}

void climateYearClass::init()
{
	// init year daily values
	G_temperature_d365 = new float [ng][365];
	G_precipitation_d365 = new float [ng][365];
	G_shortwave_d365 = new float[ng][365];
	G_longwave_d365 = new float[ng][365];
}

void climateYearClass::cleanup()
{
    // delete stored year daily values
    delete[] G_temperature_d365;
    G_temperature_d365 = 0;
	
    delete[] G_precipitation_d365;
    G_precipitation_d365 = 0;
    
    delete[] G_shortwave_d365;
    G_shortwave_d365 = 0;
	
    delete[] G_longwave_d365;
    G_longwave_d365     = 0;
}

// daily data (.365)
// read in daily data per year
void climateYearClass::read_climate_data_daily_per_year(short actual_year, short data_d)
{
  char filename[250];
	
	if (options.time_series == 3) {
        // read daily precipitation [mm]
	  sprintf(filename, "%s/G_GPCC_H08day_V20110128_%d.%d.UNF0", options.climate_dir, actual_year, data_d);
	  gridIO.readUnfFile(filename, ng * data_d, &G_precipitation_d365[0][0]);
	}
	else {
		// read temperature data
		// daily values [deg C]
		sprintf(filename, "%s/G_TEMP_H08_int_%d.%d.UNF0", options.climate_dir, actual_year, data_d);
		gridIO.readUnfFile(filename, ng * data_d, &G_temperature_d365[0][0]);
	
		// read daily precipitation [mm]
		sprintf(filename, "%s/G_GPCC_H08day_V20110128_%d.%d.UNF0", options.climate_dir, actual_year, data_d);
		gridIO.readUnfFile(filename, ng * data_d, &G_precipitation_d365[0][0]);
	
		// read daily shortwave radiation [W/m2]
		sprintf(filename, "%s/G_SSRD_H08_int_%d.%d.UNF0", options.climate_dir, actual_year, data_d);
		gridIO.readUnfFile(filename, ng * data_d, &G_shortwave_d365[0][0]);
	
		// read daily longwave radiation [W/m2]
		sprintf(filename, "%s/G_SLRD_H08_int_%d.%d.UNF0", options.climate_dir, actual_year, data_d);
		gridIO.readUnfFile(filename, ng * data_d, &G_longwave_d365[0][0]);
	}
	return;
}

// split yearly data to daily data
void climateYearClass::split_climate_data_from_year_to_daily(short dayOfTheYear, short numberDaysInMonth)
{
	short day = 0;
	short lastDayOfMonth = dayOfTheYear + numberDaysInMonth;
	for (short j = dayOfTheYear; j < lastDayOfMonth; j++) {
		for (long i = 0; i < ng; i++) {
	    if (options.time_series == 3) {
			  // read daily precipitation [mm]
			  climate.G_precipitation_d[i][day] = G_precipitation_d365[i][j];
			}
      else {
        // read temperature data
			  // daily values [oC]
			  climate.G_temperature_d[i][day] = G_temperature_d365[i][j];
			  
        // read daily precipitation [mm]
			  climate.G_precipitation_d[i][day] = G_precipitation_d365[i][j];
        
        // read daily shortwave radiation...
			  climate.G_shortwave_d[i][day] = G_shortwave_d365[i][j];
        
        // ... and longwave radiation      
			  climate.G_longwave_d[i][day] = G_longwave_d365[i][j];
			  // climate.G_longwave_d[i][day] = climateYear.G_longwave_d[i][j]; // as example.			  
		  }
		}
		day++; //[0..lastDayOfMonth-1]
	}  
  return;
}

climateYearClass::~climateYearClass(void)
{
	// nothing to do
}
