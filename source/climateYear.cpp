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
#include "common.h"
#include "def.h"
#include "globals.h"

#include "climateYear.h"

using namespace std;

extern climateClass climate;
extern optionClass options;

climateYearClass::climateYearClass(void)
{
    //nothing to do
}

void climateYearClass::init()
{
}

void climateYearClass::cleanup()
{
}

// daily data (.365)
// read in daily data per year
void climateYearClass::read_climate_data_daily_per_year(short actual_year)
{
	string climateDir = string(options.climate_dir);

	// read temperature data
	// daily values [deg C]
	G_temperature_d365.read(climateDir + "/G_TEMP_H08_int_" + to_string(actual_year) + ".365.UNF0");

	// read daily precipitation [mm]
	G_precipitation_d365.read(climateDir + "/G_GPCC_H08day_V20110128_" + to_string(actual_year) + ".365.UNF0");

	// read daily shortwave radiation [W/m2]
	G_shortwave_d365.read(climateDir + "/G_SSRD_H08_int_" + to_string(actual_year) + ".365.UNF0");

	// read daily longwave radiation [W/m2]
	G_longwave_d365.read(climateDir + "/G_SLRD_H08_int_" + to_string(actual_year) + ".365.UNF0");
}

// split yearly data to daily data
void climateYearClass::split_climate_data_from_year_to_daily(short dayOfTheYear, short numberDaysInMonth)
{
	short day = 0;
	short lastDayOfMonth = dayOfTheYear + numberDaysInMonth;
	for (short j = dayOfTheYear; j < lastDayOfMonth; j++) {
		for (long i = 0; i < ng; i++) {

			// read temperature data
			// daily values [oC]
			climate.G_temperature_d(i,day) = G_temperature_d365(i,j);

			// read daily precipitation [mm]
			climate.G_precipitation_d(i,day) = G_precipitation_d365(i,j);

			// read daily shortwave radiation...
			climate.G_shortwave_d(i,day) = G_shortwave_d365(i,j);

			// ... and longwave radiation
			climate.G_longwave_d(i,day) = G_longwave_d365(i,j);
		}
		day++; //[0..lastDayOfMonth-1]
	}
}

climateYearClass::~climateYearClass(void)
{
}
