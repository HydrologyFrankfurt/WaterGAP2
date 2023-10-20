/***********************************************************************
*
*precipitation correction after Adam and Lettenmeier instead of Leagates and
*Willmott, done by Kristina Fiedler in 2.1g. Here, correction factors are
*varying with time (temperature based). Method is described here: Döll, P., and K.
*Fiedler (2008), Global-scale modeling of groundwater recharge, Hydrology and
*Earth System Sciences, 12(3), 863-885, doi:10.5194/hess-12-863-2008. [online]
*Available from: http://www.hydrol-earth-syst-sci.net/12/863/2008/
*
*new time series options: read in of daily climate input by Kerstin
*Schulze/Verzano as .31; by Heike Hoffmann-Dobrev as .365, combination of daily
*precip and monthly climate input by Hannes Müller Schmied
*
*The time between 1971-2000 is now used for long term mean average climate
*input.
*
***********************************************************************/
//#include <iostream>
#include <cstdio>
#include <cmath>

#include "option.h"
#include "common.h"
#include "def.h"
#include "globals.h"
//#include "option.h"
//#include "def.h"
//#include "geo.h"

#include "climate.h"
using namespace std;

extern optionClass options;
extern geoClass geo;

template < class T > inline int round(const T value)
{
	return (int) floor(value + 0.5);
}


climateClass::climateClass(void)
{
	precCorrGridIsRead = 0;
}

void climateClass::init()
{
	// store daily values
	if ((options.time_series == 0) || (options.time_series == 1)) {

		G_shortwave_d.initialize();
		G_longwave_d.initialize();
		G_temperature_d.initialize();
		G_precipitation_d.initialize();
        if (options.clclOpt || options.permaOpt) {
			G_precipitation.initialize();
			G_temperature.initialize();
		}
	}
}

void climateClass::initPM()
{
	// initializations which are only necessary for Penman-Monteith
	G_windspeed_d.initialize();
	G_vappres_d.initialize();
	G_temprange_d.initialize();

}

void climateClass::initHG() //mw-H
{
	// initializations which are only necessary for Hargreaves
	G_temprange_d.initialize();
}

void climateClass::initMD()
{
    // initializations which are only necessary for Milly-Dunne
    G_temp_diff.initialize();
}

void climateClass::cleanup()
{
}



/***************** read climate data: start ***********************/
// Reading daily data from normal land mask
// daily data (.31)
void climateClass::read_climate_data_daily(short month, short actual_year)
{
	string climateDir = string(options.climate_dir);

	// new time series input option end - if time_series option == 0 (daily .31 input)
	// read temperature data
	// daily values [1 oC]
	G_temperature_d.read(climateDir + "/GTEMP_" + to_string(actual_year) + "_" + to_string(month) + ".31.UNF0");
	if ( (options.petOpt == 1 ) || (options.petOpt == 2) ){
			// read t_min data
		// daily values [1 oC]
		G_tmin_d.read(climateDir + "/GTMIN_" + to_string(actual_year) + "_" + to_string(month) + ".31.UNF0");
			// read t_max data
		// daily values [1 oC]
		G_tmax_d.read(climateDir + "/GTMAX_" + to_string(actual_year) + "_" + to_string(month) + ".31.UNF0");
	}
    if(options.calc_wtemp == 1 || options.calc_wtemp == 2)
        G_GW_Temperature_y.read(climateDir + "/GTEMP_MEAN_" + to_string(actual_year) + ".UNF0");

		// read daily precipitation [mm]
	G_precipitation_d.read(climateDir + "/GPREC_" + to_string(actual_year) + "_" + to_string(month) + ".31.UNF0");
		// read daily shortwave radiation [W/m2]...
	G_shortwave_d.read(climateDir + "/GSHORTWAVE_" + to_string(actual_year) + "_" + to_string(month) + ".31.UNF0");

	// ... and net longwave radiation [W/m2]
	if ( 0 == options.cloud ){
		G_longwave_d.read(climateDir + "/GLONGWAVE_" + to_string(actual_year) + "_" + to_string(month) + ".31.UNF0");
	}
		// ... or incoming longwave radiation for WFD [W/m2]
	if ( 1 == options.cloud ){
		G_longwave_d.read(climateDir + "/GLONGWAVE_DOWN_" + to_string(actual_year) + "_" + to_string(month) + ".31.UNF0");
	}
	// some PET-functions need vapour pressure ...
	if ((options.petOpt == 1) || (options.petOpt == 2) || (options.petOpt == 4) || (options.petOpt == 5) || (options.petOpt == 6)) {
		G_vappres_d.read(climateDir + "/GVAPPRES_" + to_string(actual_year) + "_" + to_string(month) + ".31.UNF0");
	}
	// ... some PET-functions need temperature range ...
	if ( (options.petOpt == 2) || (options.petOpt == 3)) {
		G_temprange_d.read(climateDir + "/GTEMPRANGE_" + to_string(actual_year) + "_" + to_string(month) + ".31.UNF0");
	}
	// and some PET-functions need wind speed
	if ( (options.petOpt == 1 ) || (options.petOpt == 2 ) || (options.petOpt == 4 ) || (options.petOpt == 5 )) {
			G_windspeed_d.read(climateDir + "/GWINDSPEED_" + to_string(actual_year) + "_" + to_string(month) + ".31.UNF0");
	}

}

// Reading daily data from climate land mask
void climateClass::read_climate_data_daily_withConvert(short month, short actual_year) {
    string climateDir = string(options.climate_dir);

    Daily31Grid<true, float, ng_climate> tmp_array;
    Daily31Grid<true, short, ng_climate> tmp_array_short;

    // read daily precipitation [mm]
    tmp_array.read(climateDir + "/GPREC_" + to_string(actual_year) + "_" + to_string(month) + ".31.UNF0");

    for (int cell = 0; cell < ng; cell++)
        for (int day = 0; day < 31; day++)
            G_precipitation_d(cell, day) = tmp_array(geo.G_cellnum_clm[cell] - 1, day);

    // new time series input option end - if time_series option == 0 (daily .31 input)
    // read temperature data
    // daily values [0.01 oC]
    tmp_array.read(climateDir + "/GTEMP_" + to_string(actual_year) + "_" + to_string(month) + ".31.UNF0");
    for (int cell = 0; cell < ng; cell++)
        for (int day = 0; day < 31; day++)
            G_temperature_d(cell, day) = tmp_array(geo.G_cellnum_clm[cell] - 1, day);

    if ((options.petOpt == 1) || (options.petOpt == 2)) {
        // read t_min data
        // daily values [0.01 oC]
        tmp_array.read(climateDir + "/GTMIN_" + to_string(actual_year) + "_" + to_string(month) + ".31.UNF0");
        for (int cell = 0; cell < ng; cell++)
            for (int day = 0; day < 31; day++)
                G_tmin_d(cell, day) = tmp_array(geo.G_cellnum_clm[cell] - 1, day);
        // read t_max data
        // daily values [0.01 oC]
        tmp_array.read(climateDir + "/GTMAX_" + to_string(actual_year) + "_" + to_string(month) + ".31.UNF0");
        for (int cell = 0; cell < ng; cell++)
            for (int day = 0; day < 31; day++)
                G_tmax_d(cell, day) = tmp_array(geo.G_cellnum_clm[cell] - 1, day);
    }

    // read daily shortwave radiation [W/m2]...
    tmp_array.read(climateDir + "/GSHORTWAVE_" + to_string(actual_year) + "_" + to_string(month) + ".31.UNF0");
    for (int cell = 0; cell < ng; cell++)
        for (int day = 0; day < 31; day++)
            G_shortwave_d(cell, day) = tmp_array(geo.G_cellnum_clm[cell] - 1, day);

    // ... and net longwave radiation [W/m2]
    if (0 == options.cloud) {
        tmp_array.read(climateDir + "/GLONGWAVE_" + to_string(actual_year) + "_" + to_string(month) + ".31.UNF0");
        for (int cell = 0; cell < ng; cell++)
            for (int day = 0; day < 31; day++)
                G_longwave_d(cell, day) = tmp_array(geo.G_cellnum_clm[cell] - 1, day);
    }
    // ... or incoming longwave radiation for WFD [W/m2]
    if (1 == options.cloud) {
        tmp_array.read(climateDir + "/GLONGWAVE_DOWN_" + to_string(actual_year) + "_" + to_string(month) + ".31.UNF0");
        for (int cell = 0; cell < ng; cell++)
            for (int day = 0; day < 31; day++)
                G_longwave_d(cell, day) = tmp_array(geo.G_cellnum_clm[cell] - 1, day);
    }
    // some PET-functions need vapour pressure ...
    if ((options.petOpt == 1) || (options.petOpt == 2) || (options.petOpt == 4) || (options.petOpt == 5) ||
        (options.petOpt == 6)) {
        tmp_array_short.read(climateDir + "/GVAPPRES_" + to_string(actual_year) + "_" + to_string(month) + ".31.UNF2");
        for (int cell = 0; cell < ng; cell++)
            for (int day = 0; day < 31; day++)
                G_vappres_d(cell, day) = tmp_array_short(geo.G_cellnum_clm[cell] - 1, day);
    }
    // ... some PET-functions need temperature range ...
    if ((options.petOpt == 2) || (options.petOpt == 3)) {
        tmp_array_short.read(
                       climateDir + "/GTEMPRANGE_" + to_string(actual_year) + "_" + to_string(month) + ".31.UNF2");
        for (int cell = 0; cell < ng; cell++)
            for (int day = 0; day < 31; day++)
                G_temprange_d(cell, day) = tmp_array_short(geo.G_cellnum_clm[cell] - 1, day);
    }
    // and some PET-functions need wind speed
    if (options.petOpt && (options.petOpt != 3) && (options.petOpt != 7)) {
        tmp_array_short.read(
                       climateDir + "/GWINDSPEED_" + to_string(actual_year) + "_" + to_string(month) + ".31.UNF1");
        for (int cell = 0; cell < ng; cell++)
            for (int day = 0; day < 31; day++)
                G_windspeed_d(cell, day) = tmp_array_short(geo.G_cellnum_clm[cell] - 1, day);
    }

}
// Reading long-term average data from normal land mask
void climateClass::read_climate_longtermAvg()
{
	string climateDir = string(options.climate_dir);

    if (options.clclOpt) {
		// time series input option

        // read temperature data
        // monthly values [0.01 oC]
		G_temperature.read(climateDir + "/GTEMP_1971_2000.12.UNF2");

        // read monthly precipitation [mm]
		G_precipitation.read(climateDir + "/GPREC_1971_2000.12.UNF2");

    }
}

// Reading long-term average data from climate land mask
void climateClass::read_climate_longtermAvg_withConvert()
{
	string climateDir = string(options.climate_dir);
    MonthlyGrid<true,short,ng_climate> tmp_array;

    if (options.clclOpt) {

		// time series input option
        // read temperature data
        // monthly values [0.01 oC]
		tmp_array.read(climateDir + "/GTEMP_1971_2000.12.UNF2");
        for (int cell=0; cell<ng; cell++)
            for (int mon=0;mon<12;mon++)
                G_temperature(cell,mon) = tmp_array(geo.G_cellnum_clm[cell] - 1,mon);

        // read monthly precipitation [mm]
		tmp_array.read(climateDir + "/GPREC_1971_2000.12.UNF2");
        for (int cell=0; cell<ng; cell++)
            for (int mon=0;mon<12;mon++)
                G_precipitation(cell,mon) = tmp_array(geo.G_cellnum_clm[cell] - 1,mon);

    }
}

// read temperature reduction  factor for PET-MD scheme (temp_diff), temp_diffi is the reduction factor for ith year.

// temp_diff is computed based on the following equations:
// temp_diffi= Ti(annual mean, (i-10)-(i+9) - T(annual mean, 1981-2000) for i = 2001-2090

//where:. temp_diffi is the temperature reduction factor for the ith year.
//Ti(annual mean, (i-10)-(i+9)  is annual mean temperature from the year i-10 to i+9 year (20 years)
//T(annual mean, 1981-2000) is annual mean temperature from baseline period (1981-2000)
// for the years i = 2091-2099, since the time window is smaller than 20 years, we used the annual mean temperature
// of the remaining years instead of 20 years.

void climateClass::read_annual_temp_reduction(short actual_year)
{
    string climateDir = string(options.climate_dir);
    if (options.petOpt == 7){
        // read temp reduction factors (temp_diff) for the running year
        // annual temp reduction factors for each grid cell [ oC]
        G_temp_diff.read(climateDir + "/G_TEMP_DIFF_" + to_string(actual_year) + ".UNF0");
    }
}

/***************** read climate data: end ***********************/

climateClass::~climateClass(void)
{
}
