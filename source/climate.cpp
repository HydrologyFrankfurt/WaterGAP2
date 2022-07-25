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
#include <iostream>
#include <cstdio>
#include <cmath>

#include "gridio.h"
#include "option.h"
#include "def.h"
#include "geo.h"

#include "climate.h"
using namespace std;

extern gridioClass gridIO;
// FP20161018N003 Enable reading input data from climate land mask
extern optionClass options;
extern geoClass geo;

template < class T > inline int round(const T value)
{
	return (int) floor(value + 0.5);
}


climateClass::climateClass(void)
{
	precCorrGridIsRead = 0;
	G_prec_correct = NULL;
	/*
	 * for (int n=0; n<=ng-1; n++) {
	 * for (short m=0;m<=11; m++) {
	 * G_temperature[n][m]   = 0;
	 * G_precipitation[n][m] = 0;
	 * G_raindays[n][m]      = 0;
	 * G_shortwave[n][m]     = 0;
	 * G_prec_correct[n][m]  = 100;
	 * }
	 * }
	 */
}

void climateClass::init()
{
	// store daily values
        if ((options.time_series == 0) || (options.time_series == 1) || (options.time_series == 2) || (options.time_series == 3)) {

		G_shortwave_d = new float[ng][31];
		G_longwave_d = new float[ng][31];
		G_temperature_d = new float [ng][31];
		G_precipitation_d = new float [ng][31];
        if (options.clclOpt || options.permaOpt) {
			G_precipitation = new short [ng][12];
			G_temperature = new short [ng][12];
		}
	}
  // store monthly values
    if (options.time_series > 1) {
		G_raindays = new unsigned char [ng][12];
		G_shortwave = new short[ng][12];
		G_longwave  = new short[ng][12];
		G_sunshine = new short[ng][12];
		G_cloud = new char[ng][12];
		G_temperature = new short [ng][12];
		G_tmin = new short [ng][12];
		G_tmax = new short [ng][12];
		G_precipitation = new short [ng][12];
	}
}

void climateClass::cleanup()
{
    // daily values stored
    if ((options.time_series == 0) || (options.time_series == 1) || (options.time_series == 2) || (options.time_series == 3)) {

		delete[] G_shortwave_d;      G_shortwave_d      = 0;
		delete[] G_longwave_d;       G_longwave_d     = 0;
		delete[] G_temperature_d;    G_temperature_d   = 0;
		delete[] G_precipitation_d;  G_precipitation_d = 0;
		if (options.clclOpt || options.permaOpt) {
			delete[] G_precipitation;  G_precipitation = 0;
			delete[] G_temperature;    G_temperature   = 0;
		}
	}
	// monthly values stored
    if (options.time_series > 1) {

		delete[] G_raindays;       G_raindays      = 0;
		delete[] G_shortwave;      G_shortwave     = 0;
		delete[] G_longwave;       G_longwave      = 0;
		delete[] G_sunshine;       G_sunshine      = 0;
		delete[] G_cloud;          G_cloud         = 0;
		delete[] G_temperature;    G_temperature   = 0;
		delete[] G_tmin;		   G_tmin 	       = 0;
		delete[] G_tmax;    	   G_tmax	  	   = 0;
		delete[] G_precipitation;  G_precipitation = 0;
	}
}

void climateClass::createTempAvg(const char *climate_dir, const short start_year, const short end_year)
{
	cout << "Averaging temperature data ...\n";

	short Grid[ng][12];

	int G_temperature_sum[ng][12];

	// has to be int instead of short
	// otherwise overflows might occur

	char filename[250];
	int n;

	for (n = 0; n <= ng - 1; n++) {
		for (short m = 0; m <= 11; m++) {
			G_temperature_sum[n][m] = 0;
		}
	}
	for (short year = start_year; year <= end_year; year++) {
		sprintf(filename, "%s/GTEMP_%d.12.UNF2", climate_dir, year);
		gridIO.readUnfFile(filename, ng * 12, &Grid[0][0]);
		for (n = 0; n <= ng - 1; n++) {
			for (short m = 0; m <= 11; m++) {
				G_temperature_sum[n][m] += Grid[n][m];
			}
		}
	}
	short numberOfYears = end_year - start_year + 1;

	for (n = 0; n <= ng - 1; n++) {
		for (short m = 0; m <= 11; m++) {
			G_temperature[n][m]
				= round(G_temperature_sum[n][m] / (float) numberOfYears);
		}
	}
}


void climateClass::createPrecAvg(const char *climate_dir, const short start_year, const short end_year)
{
	cout << "Averaging precipitation data ...\n";

	short Grid[ng][12];

	int G_precipitation_sum[ng][12];

	// has to be int instead of short
	// otherwise overflows might occur

	char filename[250];
	int n;

	for (n = 0; n <= ng - 1; n++) {
		for (short m = 0; m <= 11; m++) {
			G_precipitation_sum[n][m] = 0;
		}
	}
	for (short year = start_year; year <= end_year; year++) {
		sprintf(filename, "%s/GPREC_%d.12.UNF2", climate_dir, year);
		gridIO.readUnfFile(filename, ng * 12, &Grid[0][0]);
		for (n = 0; n <= ng - 1; n++) {
			for (short m = 0; m <= 11; m++) {
				G_precipitation_sum[n][m] += Grid[n][m];
			}
		}
	}
	short numberOfYears = end_year - start_year + 1;

	for (n = 0; n <= ng - 1; n++) {
		for (short m = 0; m <= 11; m++) {
			G_precipitation[n][m]
				= round(G_precipitation_sum[n][m] / (float) numberOfYears);
		}
	}
}


void climateClass::createRadiationAvg(const char *climate_dir, const short start_year, const short end_year)
{
	cout << "Averaging radiation data ...\n";

	short Grid[ng][12];

	int G_shortwave_sum[ng][12];

	// has to be int instead of short
	// otherwise overflows might occur

	char filename[250];
	int n;

	for (n = 0; n <= ng - 1; n++) {
		for (short m = 0; m <= 11; m++) {
			G_shortwave_sum[n][m] = 0;
		}
	}
	for (short year = start_year; year <= end_year; year++) {
		sprintf(filename, "%s/GRADIATION_%d.12.UNF2", climate_dir, year);
		gridIO.readUnfFile(filename, ng * 12, &Grid[0][0]);
		for (n = 0; n <= ng - 1; n++) {
			for (short m = 0; m <= 11; m++) {
				G_shortwave_sum[n][m] += Grid[n][m];
			}
		}
	}
	short numberOfYears = end_year - start_year + 1;

	for (n = 0; n <= ng - 1; n++) {
		for (short m = 0; m <= 11; m++) {
			G_shortwave[n][m] = round(G_shortwave_sum[n][m] / (float) numberOfYears);
		}
	}
}
void climateClass::createSunshineAvg(const char *climate_dir, const short start_year, const short end_year)
{
	cout << "Averaging sunshine data ...\n";

	short Grid[ng][12];

	int G_sunshine_sum[ng][12];

	// has to be int instead of short
	// otherwise overflows might occur

	char filename[250];
	int n;

	for (n = 0; n <= ng - 1; n++) {
		for (short m = 0; m <= 11; m++) {
			G_sunshine_sum[n][m] = 0;
		}
	}
	for (short year = start_year; year <= end_year; year++) {
		sprintf(filename, "%s/GSUNSHINE_%d.12.UNF2", climate_dir, year);
		gridIO.readUnfFile(filename, ng * 12, &Grid[0][0]);
		for (n = 0; n <= ng - 1; n++) {
			for (short m = 0; m <= 11; m++) {
				G_sunshine_sum[n][m] += Grid[n][m];
			}
		}
	}
	short numberOfYears = end_year - start_year + 1;

	for (n = 0; n <= ng - 1; n++) {
		for (short m = 0; m <= 11; m++) {
			G_sunshine[n][m] = round(G_sunshine_sum[n][m] / (float) numberOfYears);
		}
	}
}


void climateClass::createCloudAvg(const char *climate_dir, const short start_year, const short end_year)
{
	cout << "Averaging cloudiness data ...\n";

	char Grid[ng][12];

	short G_cloud_sum[ng][12];

	// has to be short instead of 'char'
	// otherwise overflows might occur

	char filename[250];
	int n;

	for (n = 0; n <= ng - 1; n++) {
		for (short m = 0; m <= 11; m++) {
			G_cloud_sum[n][m] = 0;
		}
	}
	for (short year = start_year; year <= end_year; year++) {
		sprintf(filename, "%s/GCLOUD_%d.12.UNF2", climate_dir, year);
		gridIO.readUnfFile(filename, ng * 12, &Grid[0][0]);
		for (int n = 0; n <= ng - 1; n++) {
			for (short m = 0; m <= 11; m++) {
				G_cloud_sum[n][m] += Grid[n][m];
			}
		}
	}
	short numberOfYears = end_year - start_year + 1;

	for (n = 0; n <= ng - 1; n++) {
		for (short m = 0; m <= 11; m++) {
			G_cloud[n][m] = round(G_cloud_sum[n][m] / (float) numberOfYears);
		}
	}
}


void climateClass::createRaindaysAvg(const char *climate_dir, const short start_year, const short end_year)
{
	cout << "Averaging 'number of raindays' ...\n";

	char Grid[ng][12];

	short G_raindays_sum[ng][12];

	// has to be short instead of 'char'
	// otherwise overflows might occur

	char filename[250];
	int n;

	for (n = 0; n <= ng - 1; n++) {
		for (short m = 0; m <= 11; m++) {
			G_raindays_sum[n][m] = 0;
		}
	}
	for (short year = start_year; year <= end_year; year++) {
		sprintf(filename, "%s/GNRD_%d.12.UNF1", climate_dir, year);
		gridIO.readUnfFile(filename, ng * 12, &Grid[0][0]);
		for (n = 0; n <= ng - 1; n++) {
			for (short m = 0; m <= 11; m++) {
				G_raindays_sum[n][m] += Grid[n][m];
			}
		}
	}
	short numberOfYears = end_year - start_year + 1;

	for (n = 0; n < ng; ++n) {
		for (short m = 0; m <= 11; ++m) {
			G_raindays[n][m] = round(G_raindays_sum[n][m] / (float) numberOfYears);
		}
	}
}


void climateClass::readPrecAdjustFactors(const char *input_dir, const char *climate_dir)  //Adam Lettenmaier precipcorr
{
	if (!precCorrGridIsRead) {
		char filename[250];

		sprintf(filename, "%s/G_PREC_CORR.12.UNF0", input_dir);

		G_prec_correct = new float[ng][12];

		gridIO.readUnfFile(filename, ng * 12, &G_prec_correct[0][0]);

		precCorrGridIsRead = true;
	}
        // reading long-term mean temperature data
        float Grid[ng][12];
	G_temp_mean = new float[ng][12];
	char filename[250];
	for (int n = 0; n < ng; n++) {
		for (short m = 0; m < 12; m++) {
			G_temp_mean[n][m] = 0.0;
		}
	}

        sprintf(filename, "%s/GTEMP_1971_2000.12.UNF0", climate_dir); //data needs to be with decimal places (not multiplied by 100)
	gridIO.readUnfFile(filename, ng * 12, &Grid[0][0]);
	for (int n = 0; n < ng; n++) {
		for (short m = 0; m < 12; m++) {
			G_temp_mean[n][m] += float (Grid[n][m]);
		}
	}
	
        // calculate snow percentage depending on
        // long-term mean temperature
	R_snow_mean = new float[ng][12];
	float R_exp;
	for (int n = 0; n < ng; ++n){
		for (short m = 0; m < 12; ++m){
			R_exp = pow(float(1.35), G_temp_mean[n][m]);
			R_snow_mean[n][m] = 1.0 / (1.0 + 1.61 * R_exp);
		}
	}
}

void climateClass::adjustPrecipitation()
{
        // new precipitation correction factors from Adam and Lettenmaier are used
        // monthly precipitation sum is adjusted according to monthly mean temperatures
	
        // calculate snow percentage depending on actual monthly mean temperature
	R_snow = new float[ng][12]; 
	G_temperature_float = new float[ng][12];
	float R_exp;

	for (int n = 0; n < ng; ++n){
		for (short m = 0; m < 12; ++m){
			G_temperature_float[n][m]= float(G_temperature[n][m] * 0.01);
			R_exp = pow(float(1.35), G_temperature_float[n][m]);
			R_snow[n][m] = 1.0 / (1.0 + 1.61 * R_exp);
		}
	}

        // correct precipitation
	for (int n = 0; n < ng; ++n){
		for (short i = 0; i < 12; ++i){
			Corr_factor_new[n][i] = (G_prec_correct[n][i]-100.0) * (R_snow[n][i] / R_snow_mean[n][i]) + 100.0;
			if (Corr_factor_new[n][i] > 300.0)
				Corr_factor_new[n][i] = 300.0; //max
			if (G_temp_mean[n][i] < 3.0) //correction factor will be adjusted
				G_precipitation[n][i] = (short) round (G_precipitation[n][i] * 0.01 * Corr_factor_new[n][i]);
			if (G_temp_mean[n][i] >= 3.0) //original correction factor will be used
				G_precipitation[n][i] = (short) round (G_precipitation[n][i] * 0.01 * G_prec_correct[n][i]);
		}
	}

	delete[]G_temperature_float;
	delete[]R_snow;
}

/***************** read climate data: start ***********************/
// Reading daily data from normal land mask
// daily data (.31)
void climateClass::read_climate_data_daily(short month, short actual_year, short data_m_d)
{
	char filename[250];
  // new time series input option start - if time_series option == 2 (only daily precip .31 input)
  if (options.time_series == 2) {
     	// read daily precipitation [mm]
	  sprintf(filename, "%s/GPREC_%d_%d.%d.UNF0", options.climate_dir, actual_year, month, data_m_d);
	  gridIO.readUnfFile(filename, ng * data_m_d, &G_precipitation_d[0][0]);
  } 
  else {  // new time series input option end - if time_series option == 0 (daily .31 input)
	  // read temperature data
          // daily values [1 oC]
	  sprintf(filename, "%s/GTEMP_%d_%d.%d.UNF0", options.climate_dir, actual_year, month, data_m_d);
	  gridIO.readUnfFile(filename, ng * data_m_d, &G_temperature_d[0][0]);
	  if ( (options.petOpt == 1 ) || (options.petOpt == 2) ){
		  // read t_min data
                  // daily values [1 oC]
		  sprintf(filename, "%s/GTMIN_%d_%d.%d.UNF0", options.climate_dir, actual_year, month, data_m_d);
		  gridIO.readUnfFile(filename, ng * data_m_d, &G_tmin_d[0][0]);
		  // read t_max data
                  // daily values [1 oC]
		  sprintf(filename, "%s/GTMAX_%d_%d.%d.UNF0", options.climate_dir, actual_year, month, data_m_d);
		  gridIO.readUnfFile(filename, ng * data_m_d, &G_tmax_d[0][0]);
	  }
	  // read daily precipitation [mm]
	  sprintf(filename, "%s/GPREC_%d_%d.%d.UNF0", options.climate_dir, actual_year, month, data_m_d);
	  gridIO.readUnfFile(filename, ng * data_m_d, &G_precipitation_d[0][0]);
	  // read daily shortwave radiation [W/m2]...
	  if ((2 == options.cloud)||( 3 == options.cloud)||( 4 == options.cloud)){
		  sprintf(filename, "%s/GSHORTWAVE_%d_%d.%d.UNF0", options.climate_dir, actual_year, month, data_m_d);
		  gridIO.readUnfFile(filename, ng * data_m_d, &G_shortwave_d[0][0]);
	  }
	  // ... and net longwave radiation [W/m2]
	  if ( 2 == options.cloud ){
	    sprintf(filename, "%s/GLONGWAVE_%d_%d.%d.UNF0", options.climate_dir, actual_year, month, data_m_d);
	    gridIO.readUnfFile(filename, ng * data_m_d, &G_longwave_d[0][0]);
	  }
          // ... or incoming longwave radiation for WFD [W/m2]
  	  if ( 3 == options.cloud ){
  	    sprintf(filename, "%s/GLONGWAVE_DOWN_%d_%d.%d.UNF0", options.climate_dir, actual_year, month, data_m_d);
  	    gridIO.readUnfFile(filename, ng * data_m_d, &G_longwave_d[0][0]);
  	  }
	  // some PET-functions need vapour pressure ...
	  if ((options.petOpt == 1) || (options.petOpt == 2) || (options.petOpt == 4) || (options.petOpt == 5) || (options.petOpt == 6)) {
		  sprintf(filename, "%s/GVAPPRES_%d_%d.%d.UNF2", options.climate_dir, actual_year, month, data_m_d);
		  gridIO.readUnfFile(filename, ng * data_m_d, &G_vappres_d[0][0]);
	  }
	  // ... some PET-functions need temperature range ...
	  if ( (options.petOpt == 2) || (options.petOpt == 3)) {
		  sprintf(filename, "%s/GTEMPRANGE_%d_%d.%d.UNF2", options.climate_dir, actual_year, month, data_m_d);
		  gridIO.readUnfFile(filename, ng * data_m_d, &G_temprange_d[0][0]);
	  }
	  // and some PET-functions need wind speed
	  if ( options.petOpt ) {
		  if ( options.petOpt != 3 ){
			  sprintf(filename, "%s/GWINDSPEED_%d_%d.%d.UNF1", options.climate_dir, actual_year, month, data_m_d);
			  gridIO.readUnfFile(filename, ng * data_m_d, &G_windspeed_d[0][0]);
		  }
	  }
  }
	return;
}
// end read_climate_data_daily


// FP20161018N003 Enable reading input data from climate land mask
// Reading daily data from climate land mask
void climateClass::read_climate_data_daily_withConvert(short month, short actual_year, short data_m_d)
{
    char filename[250];
    float (*tmp_array)[31] = new float[ng_climate][31];
    short (*tmp_array_short)[31] = new short[ng_climate][31];

    // read daily precipitation [mm]
    sprintf(filename, "%s/GPREC_%d_%d.%d.UNF0", options.climate_dir, actual_year, month, data_m_d);
    gridIO.readUnfFile(filename, ng_climate * data_m_d, &tmp_array[0][0]);
    for (int cell = 0; cell < ng; cell++)
        for (int day = 0; day < 31; day++)
            G_precipitation_d[cell][day] = tmp_array[geo.G_cellnum_clm[cell] - 1][day];

    if (options.time_series != 2){ // new time series input option end - if time_series option == 0 (daily .31 input)
        // read temperature data
        // daily values [0.01 oC]
        sprintf(filename, "%s/GTEMP_%d_%d.%d.UNF0", options.climate_dir, actual_year, month, data_m_d);
        gridIO.readUnfFile(filename, ng_climate * data_m_d, &tmp_array[0][0]);
        for (int cell = 0; cell < ng; cell++)
            for (int day = 0; day < 31; day++)
                G_temperature_d[cell][day] = tmp_array[geo.G_cellnum_clm[cell] - 1][day];

        if ((options.petOpt == 1) || (options.petOpt == 2)){
            // read t_min data
            // daily values [0.01 oC]
            sprintf(filename, "%s/GTMIN_%d_%d.%d.UNF0", options.climate_dir, actual_year, month, data_m_d);
            gridIO.readUnfFile(filename, ng_climate * data_m_d, &tmp_array[0][0]);
            for (int cell = 0; cell < ng; cell++)
                for (int day = 0; day < 31; day++)
                    G_tmin_d[cell][day] = tmp_array[geo.G_cellnum_clm[cell] - 1][day];
            // read t_max data
            // daily values [0.01 oC]
            sprintf(filename, "%s/GTMAX_%d_%d.%d.UNF0", options.climate_dir, actual_year, month, data_m_d);
            gridIO.readUnfFile(filename, ng_climate * data_m_d, &tmp_array[0][0]);
            for (int cell = 0; cell < ng; cell++)
                for (int day = 0; day < 31; day++)
                    G_tmax_d[cell][day] = tmp_array[geo.G_cellnum_clm[cell] - 1][day];
        }

        // read daily shortwave radiation [W/m2]...
        if ((2 == options.cloud) || (3 == options.cloud) || (4 == options.cloud)){
            sprintf(filename, "%s/GSHORTWAVE_%d_%d.%d.UNF0", options.climate_dir,	actual_year, month, data_m_d);
            gridIO.readUnfFile(filename, ng_climate * data_m_d, &tmp_array[0][0]);
            for (int cell = 0; cell < ng; cell++)
                for (int day = 0; day < 31; day++)
                    G_shortwave_d[cell][day] = tmp_array[geo.G_cellnum_clm[cell] - 1][day];
        }
        // ... and net longwave radiation [W/m2]
        if (2 == options.cloud){
            sprintf(filename, "%s/GLONGWAVE_%d_%d.%d.UNF0", options.climate_dir, actual_year, month, data_m_d);
            gridIO.readUnfFile(filename, ng_climate * data_m_d, &tmp_array[0][0]);
            for (int cell = 0; cell < ng; cell++)
                for (int day = 0; day < 31; day++)
                    G_longwave_d[cell][day] = tmp_array[geo.G_cellnum_clm[cell] - 1][day];
        }
        // ... or incoming longwave radiation for WFD [W/m2]
        if (3 == options.cloud) {
            sprintf(filename, "%s/GLONGWAVE_DOWN_%d_%d.%d.UNF0", options.climate_dir, actual_year, month, data_m_d);
            gridIO.readUnfFile(filename, ng_climate * data_m_d, &tmp_array[0][0]);
            for (int cell = 0; cell < ng; cell++)
                for (int day = 0; day < 31; day++)
                    G_longwave_d[cell][day] = tmp_array[geo.G_cellnum_clm[cell] - 1][day];
        }
        // some PET-functions need vapour pressure ...
        if ((options.petOpt == 1) || (options.petOpt == 2) || (options.petOpt == 4)	|| (options.petOpt == 5) || (options.petOpt == 6)){
            sprintf(filename, "%s/GVAPPRES_%d_%d.%d.UNF2", options.climate_dir,	actual_year, month, data_m_d);
            gridIO.readUnfFile(filename, ng_climate * data_m_d, &tmp_array_short[0][0]);
            for (int cell = 0; cell < ng; cell++)
                for (int day = 0; day < 31; day++)
                    G_vappres_d[cell][day] = tmp_array_short[geo.G_cellnum_clm[cell] - 1][day];
        }
        // ... some PET-functions need temperature range ...
        if ((options.petOpt == 2) || (options.petOpt == 3))	{
            sprintf(filename, "%s/GTEMPRANGE_%d_%d.%d.UNF2", options.climate_dir,	actual_year, month, data_m_d);
            gridIO.readUnfFile(filename, ng_climate * data_m_d, &tmp_array_short[0][0]);
            for (int cell = 0; cell < ng; cell++)
                for (int day = 0; day < 31; day++)
                    G_temprange_d[cell][day] = tmp_array_short[geo.G_cellnum_clm[cell] - 1][day];
        }
        // and some PET-functions need wind speed
        if (options.petOpt && (options.petOpt != 3)) {
            sprintf(filename, "%s/GWINDSPEED_%d_%d.%d.UNF1", options.climate_dir,	actual_year, month, data_m_d);
            gridIO.readUnfFile(filename, ng_climate * data_m_d, &tmp_array_short[0][0]);
            for (int cell = 0; cell < ng; cell++)
                for (int day = 0; day < 31; day++)
                    G_windspeed_d[cell][day] = tmp_array_short[geo.G_cellnum_clm[cell] - 1][day];

        }
    }

    delete[] tmp_array;					tmp_array = NULL;
    delete[] tmp_array_short;		tmp_array_short = NULL;

    return;
}
// end read_climate_data_daily_withConvert


// Reading monthly data from normal land mask
// monthly data
void climateClass::read_climate_data_monthly(short actual_year, short data_m_d)
{
	char filename[250];

	// read temperature data
	// daily values [0.01 oC]
	sprintf(filename, "%s/GTEMP_%d.%d.UNF2", options.climate_dir, actual_year, data_m_d);
	gridIO.readUnfFile(filename, ng * data_m_d, &G_temperature[0][0]);

	if ( (options.petOpt == 1 ) || (options.petOpt == 2) ){
		// read t_min data
		// daily values [0.01 oC]
		sprintf(filename, "%s/GTMIN_%d.%d.UNF2", options.climate_dir, actual_year, data_m_d);
		gridIO.readUnfFile(filename, ng * data_m_d, &G_tmin[0][0]);
		// read t_max data
		// daily values [0.01 oC]
		sprintf(filename, "%s/GTMAX_%d.%d.UNF2", options.climate_dir, actual_year, data_m_d);
		gridIO.readUnfFile(filename, ng * data_m_d, &G_tmax[0][0]);
	}
        if (4 == options.time_series) {    //new time series input option
	// read monthly precipitation [mm]
		sprintf(filename, "%s/GPREC_%d.%d.UNF2", options.climate_dir, actual_year, data_m_d);
		gridIO.readUnfFile(filename, ng * data_m_d, &G_precipitation[0][0]);
		if (1 == options.prec_correct)
			adjustPrecipitation();
		if (options.raindays == 1) {
			sprintf(filename, "%s/GNRD_%d.%d.UNF1", options.climate_dir, actual_year, data_m_d);
			gridIO.readUnfFile (filename, ng * data_m_d, &G_raindays[0][0]);
		}
	}
	// read daily sunshine...
	if (0 == options.cloud) {
		sprintf(filename, "%s/GSUNSHINE_%d.%d.UNF2", options.climate_dir, actual_year, data_m_d);
		gridIO.readUnfFile(filename, ng * data_m_d, &G_sunshine[0][0]);
	}
	// ... or cloudiness...
	if (1 == options.cloud) {
		sprintf(filename, "%s/GCLOUD_%d.%d.UNF1", options.climate_dir, actual_year, data_m_d);
		gridIO.readUnfFile(filename, ng * data_m_d, &G_cloud[0][0]);
		for (short i = 0; i < data_m_d; ++i){
			for (int n = 0; n < ng; ++n){
				G_sunshine[n][i] = (100 - G_cloud[n][i]) * 10;
				// make sure, that sunshine does not exceed the maximum of 1000 (=100%)
				if (G_sunshine[n][i] > 1000)
					G_sunshine[n][i] = 1000;
			}
		}
	}
	// ... or radiation
	if ((2 == options.cloud) || (3 == options.cloud) || (4 == options.cloud)) {
		sprintf(filename, "%s/GSHORTWAVE_%d.%d.UNF2", options.climate_dir, actual_year, data_m_d);
		gridIO.readUnfFile (filename, ng * data_m_d, &G_shortwave[0][0]);

		if (4 != options.cloud){
			sprintf(filename, "%s/GLONGWAVE_%d.%d.UNF2", options.climate_dir, actual_year, data_m_d);
			gridIO.readUnfFile (filename, ng * data_m_d, &G_longwave[0][0]);
		}
	}

	// some PET-functions need vapour pressure
	if ((options.petOpt == 1) || (options.petOpt == 2) || (options.petOpt == 4) || (options.petOpt == 5) || (options.petOpt == 6)) {
		sprintf(filename, "%s/GVAPPRES_%d.%d.UNF2", options.climate_dir, actual_year, data_m_d);
		gridIO.readUnfFile(filename, ng * data_m_d, &G_vappres[0][0]);
	}
	// and some PET-functions need temperature range 
	if ( (options.petOpt == 2) || (options.petOpt == 3)) {
		sprintf(filename, "%s/GTEMPRANGE_%d.%d.UNF2", options.climate_dir, actual_year, data_m_d);
		gridIO.readUnfFile(filename, ng * data_m_d, &G_temprange[0][0]);
	}

	return;
}
// end read_climate_data_monthly


// FP20161018N003 Enable reading input data from climate land mask
// Reading daily data from climate land mask
void climateClass::read_climate_data_monthly_withConvert(short actual_year, short data_m_d)
{
    char filename[250];
    short int (*tmp_array)[12] = new short int [ng_climate][12];

    // read temperature data
    // daily values [0.01 oC]
    sprintf(filename, "%s/GTEMP_%d.%d.UNF2", options.climate_dir, actual_year, data_m_d);
    gridIO.readUnfFile(filename, ng_climate * data_m_d, &tmp_array[0][0]);

        for (int cell=0; cell<ng; cell++)
          for (int mon=0;mon<12;mon++)
            G_temperature[cell][mon] = tmp_array[geo.G_cellnum_clm[cell] - 1][mon];

    if ( (options.petOpt == 1 ) || (options.petOpt == 2) ){
        // read t_min data
        // daily values [0.01 oC]
        sprintf(filename, "%s/GTMIN_%d.%d.UNF2", options.climate_dir, actual_year, data_m_d);
        gridIO.readUnfFile(filename, ng_climate * data_m_d, &tmp_array[0][0]);
        for (int cell=0; cell<ng; cell++)
          for (int mon=0;mon<12;mon++)
            G_tmin[cell][mon] = tmp_array[geo.G_cellnum_clm[cell] - 1][mon];

        // read t_max data
        // daily values [0.01 oC]
        sprintf(filename, "%s/GTMAX_%d.%d.UNF2", options.climate_dir, actual_year, data_m_d);
        gridIO.readUnfFile(filename, ng_climate * data_m_d, &tmp_array[0][0]);
        for (int cell=0; cell<ng; cell++)
          for (int mon=0;mon<12;mon++)
            G_tmax[cell][mon] = tmp_array[geo.G_cellnum_clm[cell] - 1][mon];
    }
    if (4 == options.time_series) {    //new time series input option
    // read monthly precipitation [mm]
        sprintf(filename, "%s/GPREC_%d.%d.UNF2", options.climate_dir, actual_year, data_m_d);
        gridIO.readUnfFile(filename, ng_climate * data_m_d, &tmp_array[0][0]);
        for (int cell=0; cell<ng; cell++)
          for (int mon=0;mon<12;mon++)
            G_precipitation[cell][mon] = tmp_array[geo.G_cellnum_clm[cell] - 1][mon];

        if (1 == options.prec_correct)
            adjustPrecipitation();

        if (options.raindays == 1) {
            sprintf(filename, "%s/GNRD_%d.%d.UNF1", options.climate_dir, actual_year, data_m_d);
            gridIO.readUnfFile (filename, ng_climate * data_m_d, &tmp_array[0][0]);
            for (int cell=0; cell<ng; cell++)
              for (int mon=0;mon<12;mon++)
                 G_raindays[cell][mon] = tmp_array[geo.G_cellnum_clm[cell] - 1][mon];
        }
    }
    // read daily sunshine...
    if (0 == options.cloud) {
        sprintf(filename, "%s/GSUNSHINE_%d.%d.UNF2", options.climate_dir, actual_year, data_m_d);
        gridIO.readUnfFile(filename, ng_climate * data_m_d, &tmp_array[0][0]);
        for (int cell=0; cell<ng; cell++)
            for (int mon=0;mon<12;mon++)
                G_sunshine[cell][mon] = tmp_array[geo.G_cellnum_clm[cell] - 1][mon];
    }
    // ... or cloudiness...
    if (1 == options.cloud) {
        sprintf(filename, "%s/GCLOUD_%d.%d.UNF1", options.climate_dir, actual_year, data_m_d);
        gridIO.readUnfFile(filename, ng_climate * data_m_d, &tmp_array[0][0]);
        for (short i = 0; i < data_m_d; ++i){
            for (int n = 0; n < ng; ++n){
                G_sunshine[n][i] = (100 - tmp_array[geo.G_cellnum_clm[n] - 1][i]) * 10;
                // make sure, that sunshine does not exceed the maximum of 1000 (=100%)
                if (G_sunshine[n][i] > 1000)
                    G_sunshine[n][i] = 1000;
            }
        }
    }
    // ... or radiation
    if ((2 == options.cloud) || (3 == options.cloud) || (4 == options.cloud)) {
        sprintf(filename, "%s/GSHORTWAVE_%d.%d.UNF2", options.climate_dir, actual_year, data_m_d);
        gridIO.readUnfFile (filename, ng_climate * data_m_d, &tmp_array[0][0]);
        for (int cell=0; cell<ng; cell++)
            for (int mon=0;mon<12;mon++)
                G_shortwave[cell][mon] = tmp_array[geo.G_cellnum_clm[cell] - 1][mon];

        if (4 != options.cloud){
            sprintf(filename, "%s/GLONGWAVE_%d.%d.UNF2", options.climate_dir, actual_year, data_m_d);
            gridIO.readUnfFile (filename, ng_climate * data_m_d, &tmp_array[0][0]);
            for (int cell=0; cell<ng; cell++)
                for (int mon=0;mon<12;mon++)
                    G_longwave[cell][mon] = tmp_array[geo.G_cellnum_clm[cell] - 1][mon];
        }
    }

    // some PET-functions need vapour pressure
    if ((options.petOpt == 1) || (options.petOpt == 2) || (options.petOpt == 4) || (options.petOpt == 5) || (options.petOpt == 6)) {
        sprintf(filename, "%s/GVAPPRES_%d.%d.UNF2", options.climate_dir, actual_year, data_m_d);
        gridIO.readUnfFile(filename, ng_climate * data_m_d, &tmp_array[0][0]);
        for (int cell=0; cell<ng; cell++)
            for (int mon=0;mon<12;mon++)
                G_vappres[cell][mon] = tmp_array[geo.G_cellnum_clm[cell] - 1][mon];
    }
    // and some PET-functions need temperature range
    if ( (options.petOpt == 2) || (options.petOpt == 3)) {
        sprintf(filename, "%s/GTEMPRANGE_%d.%d.UNF2", options.climate_dir, actual_year, data_m_d);
        gridIO.readUnfFile(filename, ng_climate * data_m_d, &tmp_array[0][0]);
        for (int cell=0; cell<ng; cell++)
            for (int mon=0;mon<12;mon++)
                G_temprange[cell][mon] = tmp_array[geo.G_cellnum_clm[cell] - 1][mon];
    }

    delete[] tmp_array;

    return;

}
// end read_climate_longtermAvg_withConvert


// FP20161018N003 Enable reading input data from climate land mask
// Reading long-term average data from normal land mask
void climateClass::read_climate_longtermAvg()
{
    char filename[250];

    if (5 == options.time_series || options.clclOpt) { //HMS new time series input option

        // read temperature data
        // monthly values [0.01 oC]
        sprintf(filename, "%s/GTEMP_1971_2000.12.UNF2", options.climate_dir);
        gridIO.readUnfFile(filename, ng * 12, &G_temperature[0][0]);
        // read monthly precipitation [mm]
        sprintf(filename, "%s/GPREC_1971_2000.12.UNF2", options.climate_dir);
        gridIO.readUnfFile(filename, ng * 12, &G_precipitation[0][0]);

        if (1 == options.prec_correct) {
            adjustPrecipitation();
        }

    }
    // end if (5 == options.time_series || options.clclOpt)

    if (5 == options.time_series) {

        // number of rain days per month
        if (1 == options.raindays) {
            sprintf(filename, "%s/GNRD_1971_2000.12.UNF1", options.climate_dir);
            gridIO.readUnfFile(filename, ng * 12, &G_raindays[0][0]);
        }
        // read sunshine hours per month
        if (0 == options.cloud) {
            sprintf(filename, "%s/GSUNSHINE_1971_2000.12.UNF2", options.climate_dir);
            gridIO.readUnfFile(filename, ng * 12, &G_sunshine[0][0]);
        }
        // read cloudiness per month
        if (1 == options.cloud) {
            sprintf(filename, "%s/GCLOUD_1971_2000.12.UNF1", options.climate_dir);
            gridIO.readUnfFile(filename, ng * 12, &G_cloud[0][0]);
            for (short m = 0; m < 12; m++)
                for (int n = 0; n < ng; n++)
                    G_sunshine[n][m]
                            = (100 - G_cloud[n][m]) * 10;
        }
        // read radiation data per month
        if (2 >= options.cloud) {

            sprintf(filename, "%s/GSHORTWAVE_1971_2000.12.UNF2", options.climate_dir);
            gridIO.readUnfFile(filename, ng * 12, &G_shortwave[0][0]);

            if (2 == options.cloud) {
                sprintf(filename, "%s/GLONGWAVE_1971_2000.12.UNF2", options.climate_dir);
                gridIO.readUnfFile(filename, ng * 12, &G_longwave[0][0]);
            }
            else if (3 == options.cloud) {
                sprintf(filename, "%s/GLONGWAVE_1971_2000.12.UNF2", options.climate_dir);
                gridIO.readUnfFile(filename, ng * 12, &G_longwave[0][0]);
            }

        }
        // end if (2 >= options.cloud)

    }
    // end if (5 == options.time_series)

    // several PET-functions need windspeed data
    // if monthly climate data is used, only data from climate normal period is available
    // daily data will be read later...
    if (options.time_series) {
        if (options.petOpt) {
            if (options.petOpt != 3) {
                sprintf(filename, "%s/GWINDSPEED_1971_2000.12.UNF1", options.climate_dir);
                gridIO.readUnfFile(filename, ng * 12, &G_windspeed[0][0]);
            }
            // end if (options.petOpt != 3)
        }
        // end if (options.petOpt)
    }
    // end if (options.time_series)

    // return command to position of call
    return;

}
// end read_climate_longtermAvg


// FP20161018N003 Enable reading input data from climate land mask
// Reading long-term average data from climate land mask
void climateClass::read_climate_longtermAvg_withConvert()
{
    char filename[250];
    short (*tmp_array)[12] = new short [ng_climate][12];

    if (5 == options.time_series || options.clclOpt) { //HMS new time series input option

        // read temperature data
        // monthly values [0.01 oC]
        sprintf(filename, "%s/GTEMP_1971_2000.12.UNF2", options.climate_dir);
        gridIO.readUnfFile(filename, ng_climate * 12, &tmp_array[0][0]);
        for (int cell=0; cell<ng; cell++)
            for (int mon=0;mon<12;mon++)
                G_temperature[cell][mon] = tmp_array[geo.G_cellnum_clm[cell] - 1][mon];

        // read monthly precipitation [mm]
        sprintf(filename, "%s/GPREC_1971_2000.12.UNF2", options.climate_dir);
        gridIO.readUnfFile(filename, ng_climate * 12, &tmp_array[0][0]);
        for (int cell=0; cell<ng; cell++)
            for (int mon=0;mon<12;mon++)
                G_precipitation[cell][mon] = tmp_array[geo.G_cellnum_clm[cell] - 1][mon];

        if (1 == options.prec_correct) {
            adjustPrecipitation();
        }

    }
    // end if (5 == options.time_series || options.clclOpt)

    if (5 == options.time_series){

        // number of rain days per month
        if (options.raindays == 1) {
            sprintf(filename, "%s/GNRD_1971_2000.12.UNF1", options.climate_dir);
            gridIO.readUnfFile(filename, ng_climate * 12, &tmp_array[0][0]);
            for (int cell=0; cell<ng; cell++)
                for (int mon=0;mon<12;mon++)
                    G_raindays[cell][mon] = tmp_array[geo.G_cellnum_clm[cell] - 1][mon];
        }
        // read sunshine hours per month
        if (0 == options.cloud) {
            sprintf(filename, "%s/GSUNSHINE_1971_2000.12.UNF2", options.climate_dir);
            gridIO.readUnfFile(filename, ng_climate * 12, &tmp_array[0][0]);
            for (int cell=0; cell<ng; cell++)
                for (int mon=0;mon<12;mon++)
                    G_sunshine[cell][mon] = tmp_array[geo.G_cellnum_clm[cell] - 1][mon];
        }
        // read cloudiness per month
        if (1 == options.cloud) {
            sprintf(filename, "%s/GCLOUD_1971_2000.12.UNF1", options.climate_dir);
            gridIO.readUnfFile(filename, ng_climate * 12, &tmp_array[0][0]);

            for (int cell=0; cell<ng; cell++)
                for (int mon=0;mon<12;mon++)
                    G_sunshine[cell][mon] = (100 - tmp_array[geo.G_cellnum_clm[cell] - 1][mon] * 10);
        }
        // read radiation data per month
        if (2 >= options.cloud) {

            sprintf(filename, "%s/GSHORTWAVE_1971_2000.12.UNF2", options.climate_dir);
            gridIO.readUnfFile(filename, ng_climate * 12, &tmp_array[0][0]);
            for (int cell=0; cell<ng; cell++)
                for (int mon=0;mon<12;mon++)
                    G_shortwave[cell][mon] = tmp_array[geo.G_cellnum_clm[cell] - 1][mon];

            if (2 == options.cloud) {
                sprintf(filename, "%s/GLONGWAVE_1971_2000.12.UNF2",	options.climate_dir);
                gridIO.readUnfFile(filename, ng_climate * 12, &tmp_array[0][0]);
                for (int cell=0; cell<ng; cell++)
                    for (int mon=0;mon<12;mon++)
                        G_longwave[cell][mon] = tmp_array[geo.G_cellnum_clm[cell] - 1][mon];
            }
            else if (3 == options.cloud) {
                sprintf(filename, "%s/GLONGWAVE_DOWN_1971_2000.12.UNF2",	options.climate_dir);
                gridIO.readUnfFile(filename, ng_climate * 12, &tmp_array[0][0]);
                for (int cell=0; cell<ng; cell++)
                    for (int mon=0;mon<12;mon++)
                        G_longwave[cell][mon] = tmp_array[geo.G_cellnum_clm[cell] - 1][mon];
            }

        }
        // end if (2 >= options.cloud)

    }
    // end if (5 == options.time_series)

    // several PET-functions need windspeed data
    // if monthly climate data is used, only data from climate normal period is available
    // daily data will be read later...
    if (options.time_series) {
        if (options.petOpt) {
            if (options.petOpt != 3) {
                sprintf(filename, "%s/GWINDSPEED_1971_2000.12.UNF1", options.climate_dir);
                gridIO.readUnfFile(filename, ng * 12, &tmp_array[0][0]);
                for (int cell=0; cell<ng; cell++)
                    for (int mon=0;mon<12;mon++)
                        G_windspeed[cell][mon] = tmp_array[geo.G_cellnum_clm[cell] - 1][mon];
            }
            // end if (options.petOpt != 3)
        }
        // end if (options.petOpt != 3)
    }
    // end if (options.time_series)

    // delete temporary array
    delete[] tmp_array;     tmp_array = NULL;

    // return command to position of call
    return;

}
// end read_climate_longtermAvg_withConvert


/***************** read climate data: end ***********************/

void climateClass::initPM()
{
	// initializations which are only necessary for Penman-Monteith
	// store daily values 
  if ((options.time_series == 0) || (options.time_series == 1) || (options.time_series == 2) || (options.time_series == 3)) { //new time series input option
		G_windspeed_d = new unsigned char[ng][31];
		G_vappres_d = new unsigned short[ng][31];
		G_temprange_d = new unsigned short[ng][31];
	}
  // store monthly values 
	else {
		G_windspeed = new unsigned char[ng][12];
		G_vappres = new unsigned short[ng][12];
		G_temprange = new unsigned short[ng][12];
	}
}

void climateClass::initHG() //mw-H
{
	// initializations which are only necessary for Hargreaves 
	// daily values stored
  if ((options.time_series == 0) || (options.time_series == 1) || (options.time_series == 2) || (options.time_series == 3))  //new time series input option
		G_temprange_d = new unsigned short[ng][31];
	// monthly values stored
	else 
		G_temprange = new unsigned short[ng][12];
}

climateClass::~climateClass(void)
{
	delete[]G_prec_correct;
        delete[]G_temp_mean;
        delete[]R_snow_mean;
}
