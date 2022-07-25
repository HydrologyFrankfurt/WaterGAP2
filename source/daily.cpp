/***********************************************************************
 *
 *adding new PET options, arid/humid distinction using Köppen climates, max LAI
 *defined for land cover types based on LCT_22.DAT, daily input data and new
 *snow dynamic routine by Kerstin Verzano (former Schulze), to read at her
 *dissertation: Verzano, K. (2009), Climate Change Impacts on Flood Related
 *Hydrological Processes: Further Development and Application of a Global Scale
 *Hydrological Model, 166 pp., University of Kassel.
 *
 *read in daily values as .365 input and water use from groundwater by Heike
 *Hoffmann-Dobrev; added some output (radiation, daily precip by Hannes Müller
 *Schmied)
 *
 *changed snow water equivalent where albedo switches to snow albedo from > 0
 *mm to > 3 mm; added land use dependend snow albedos by Frank Voss
 *
 *initial calculation of storages and writing as outputs by Kristina Fiedler,
 *refinement and restruction by Stephanie Eisner and Heike Hoffmann-Dobrev
 *
 *change arid-humid query from "if" to "switch" by Stephanie Eisner to avoid
 *mistakes when using different input files. Note: CESR uses "2" for arid
 *areas ("1" for humid), IPG uses "1" for arid areas and "0" for humid. This
 *is one of the sparse differences in the code of CESR and IPG
 *
 *introduce land cover type dependent emissivity values (integrated in
 *LCT_22.DAT) when dealing with incoming (downward) longwave radiation with
 *values after Wilber, A. C., D. P. Kratz, and S. K. Gupta (1999), Surface
 *Emissivity Maps for Use in Satellite Retrievals of Longwave Radiation.
 *
 *enhancing and resturcturing radiation options to all different climate inputs
 *by Stephanie Eisner and Hannes Müller Schmied
 *
 *enhancing daily output options to .365 by Hannes Müller Schmied
 *
 * see former changes at file daily.cpp.versioninfos.txt
 *
 * Changes Felix Portmann (FP) 2015
 * Temporally introducing functions for check of arrays and scalars
 *
 * Items to check:
 * case 4: // FP 2015-06 ATTENTION! CHECK NEEDED: Why is there no definition and no "break" for petOpt=4 "PET: combination of Turc and Doorenbos&Pruitt"?
 *
 * Felix Portmann, Claudia Riedel, Petra Doell 2015-07
 * Corrections in treatment of fractional routing (surface storage changes),
 * e.g. limitation to maximum soil water storage after scaling
 *
***********************************************************************/
#include <cmath>
#include <cstdio>

#include "dailyprec.h"
#include "interp_day.h"
#include "calib_basins.h"
#include "option.h"
#include "lai.h"
#include "gridio.h"
#include "climate.h"
#include "climateYear.h" // for .365 (H08) data input
#include "geo.h"
#include "land.h"
#include "gw_frac.h"
#include "routing.h"	// for DDM30 land area
#include "def.h"
#include "clcl.h"
#include "s_max.h"

#include "daily.h"
using namespace std;

#ifdef _WATERGAP_CHECKS_GLOBAL_H
    // Check Felix Portmann 2015 - for invalid data
    extern short chk_year;
    extern short chk_month;
    extern short chk_day_in_month;

    extern void check_value_double_LT0(double value);
    extern void check_value_double_NaN(double value);
#endif

const double pi = 3.141592653589793;
extern optionClass options;

extern geoClass geo; // needed for geo.G_contcell[]

// #######################################
// Cell-specific calibration parameters
// FP
#include "calib_param.h"
extern calibParamClass calibParam; // JSON object and methods from watergap.cpp
// #######################################


dailyWaterBalanceClass::dailyWaterBalanceClass() {
	LCTdataIsRead = false;

}

template<class T> inline int round(const T value) {
	return (int) floor(value + 0.5);
}

void dailyWaterBalanceClass::calcNewDay(const short day, const short month,
		const short day_in_month, const short year, const int n) {
	extern laiClass lai;
	extern dailyPrecNrdClass dailyPrecNrd;
	extern climateClass climate;
        extern climateYearClass climateYear; //for .365 (H08) data input
	extern geoClass geo;
	extern landClass land;
	extern groundwaterFactorClass GW;
	extern clclClass clcl;
	extern soilWatCapClass maxSoilWaterCap;
	extern routingClass routing;

	extern short G_toBeCalculated[ng];
	extern signed short G_sbasin[ng]; // superbasins

	//------- ARID GROUNDWATER-------
	extern short G_aindex[ng]; // arid cells
	// arid_gw = true  -> calculation for arid areas
	// arid_gw = false -> calculation for humid areas
	bool arid_gw = false;

	//////////////////////////////////// temp & prec /////////////////////////////////////////////////////////////////////////////

	double dailyTempC = 0., dailyShortWave = 0., dailyLongWave = 0., dailySunshine = 0., dailyPrec = 0.;
	double dailyTmin = 0.;//, dailyTmax=0.;
    double dailyWindSpeed = 0., dailyVapPres = 0., dailySnowFall = 0., dailyRainFall = 0., dailySnowCoverFrac = 0., dailyTempRange = 0.;  // FP20161011N001
	double monthlyPrecipitation = 0.;
	char numberOfRainDays;
	// definitions for total radiation
	double net_short_wave_rad = 0.;
	double net_long_wave_rad = 0.;
	double long_wave_rad_in = 0.; // new time series input option - H08
	double long_wave_rad_out = 0.; // new time series input option - H08
	double openWaterNetShortWaveRad = 0.;
	double openWaterNetRad = 0.;
	// latent heat of evaporation
	double lat_heat = 0.;
	// solar radiation
	double solar_rad = 0.;
	// reflected solar radiation from earth (dependend on albedo) // HMS radiation output adjustment
	double refl_rad = 0.;
	// total net radiation [mm/day]
	double net_rad = 0.; // HMS 2013-11-21 defined here instead of below for single cell output
	// different surfaces with different albedo
	double albedo = 0.;
	// daily LAI development
	double dailyLai = 0.;
	// daily Kc development
	double dailyKc = 0.;
	// land area fraction for calculation of waterbalance parameters per grid cell...
	double landAreaFrac = routing.getLandAreaFrac(n);

    // variables for fractional routing (soil water content)
	double soil_saturation = 0.;
	double soil_water_overflow = 0.;


	if (1 == G_toBeCalculated[n]) {
		// added for cell AET calculation (2.1f) 
		dailyCellLandAET[n] = 0.;
		dailyCellLandAET_uncorr[n] = 0.;

		// Local copies instead of often used method getValue // FP
		// INCREASE CPU time -> NOT IMPLEMENTED
		//		double P_PTC_ARI__at__n = calibParam.getValue(P_PTC_ARI,n);
		//		double P_PTC_HUM__at__n = calibParam.getValue(P_PTC_HUM,n);
		//		double P_PET_MXDY__at__n = calibParam.getValue(P_PET_MXDY,n);

		//
		// do the interpolation of climatic variables from monthly to daily values
		// or use daily values if available
		//

		// use daily input data (.31)
		if ((options.time_series == 0) || (options.time_series == 2)) {

			// daily precipitation data
			dailyPrec = climate.G_precipitation_d[n][day_in_month - 1]; // if only .31 (WFD) precip is read in

			if (options.time_series == 0) { // if full .31 (WFD) are read in
				// daily temperature [°C]
				dailyTempC = climate.G_temperature_d[n][day_in_month - 1];

				if ((options.cloud == 2) || (options.cloud == 3) || (options.cloud == 4)) {

					dailyShortWave = climate.G_shortwave_d[n][day_in_month - 1];

					if ((options.cloud == 2) || (options.cloud == 3)) {
						dailyLongWave = climate.G_longwave_d[n][day_in_month- 1];
					}
				}
                // endif options

                // additional information depending on PET-calculation - data are only available as .31 and therefore not transformed to .365 input
				switch (options.petOpt) {
				case 1:	dailyWindSpeed = climate.G_windspeed_d[n][day_in_month - 1]	/ 10.; // [m/s]
						dailyVapPres = climate.G_vappres_d[n][day_in_month - 1]	/ 10.; // [kPa]
						dailyTmin = climate.G_tmin_d[n][day_in_month - 1]; // [oC]
						break;
				case 2:	dailyWindSpeed = climate.G_windspeed_d[n][day_in_month - 1]	/ 10.; // [m/s]
						dailyVapPres = climate.G_vappres_d[n][day_in_month - 1]	/ 10.; // [kPa]
						dailyTmin = climate.G_tmin_d[n][day_in_month - 1]; // [oC]
						break;
				case 3:	dailyTempRange = climate.G_temprange_d[n][day_in_month - 1]; // [oC]
						dailyTempRange = abs(dailyTempRange);
						break;
				case 4:
				case 5:	dailyWindSpeed = climate.G_windspeed_d[n][day_in_month - 1] / 10.; // [m/s]
						dailyVapPres = climate.G_vappres_d[n][day_in_month - 1]	/ 10.; // [kPa]
						break;
				case 6:	dailyVapPres = climate.G_vappres_d[n][day_in_month - 1]	/ 10.; //[kPa]
						break;
				}

			}
            // endif options (full .31 data)

		}
        // endif options (daily input data .31)

		// use daily input data per year (.365)
		if ((options.time_series == 1) || (options.time_series == 3)) {

			// daily precipitation data
			dailyPrec = climateYear.G_precipitation_d365[n][day - 1]; // if only .365 (H08) precip is read in

			if (options.time_series == 1) { // if full .365 (H08) are read in
				// daily temperature [°C]
                dailyTempC = climateYear.G_temperature_d365[n][day - 1];

                if ((options.cloud == 2) || (options.cloud == 3) || (options.cloud == 4)) {

					dailyShortWave = climateYear.G_shortwave_d365[n][day - 1];

					if ((options.cloud == 2) || (options.cloud == 3)) {
						dailyLongWave = climateYear.G_longwave_d365[n][day - 1];
					}
                    // endif options (cloud 2,3)

				}
                // endif options (cloud 2,3,4)

			}
            // endif options full .365 data

		}
        // endif options (daily input data per year (.365)

		// daily data interpolation from monthly data
		if ((options.time_series != 0) && (options.time_series != 1)) { // new time series input option - if somehow monthly data should be used

			dailyTempC = interp_day_temp(n, day, climate.G_temperature) / 100.; // [oC]
                        // G_temperature is [oC] exept from CRU [100 oC] and have therefore be divided by 100
			// Tmax values are never used rigth now 
			// dailyTmax = interp_day_temp(n, day, climate.G_tmax) / 100.;	// [oC]
			// G_temperature is [0.01 oC]

			//===========daily temperature variation for snowpack ====================
			//=== not used in 2.1f because effects are negligible
			//float dailyVarTempC;
			//dailyTempC = dailyTempC + dailyTempVar[day]; //comment out if no temp variation
			//=========== end daily temperature variation============================

			switch (options.cloud) {
			case 0:	dailySunshine = interp_day_sunshine(n, day, climate.G_sunshine);
					break;
			case 1:	dailySunshine = interp_day_sunshine(n, day, climate.G_sunshine);
					break;
			case 2: dailyShortWave = interp_day_sunshine(n, day, climate.G_shortwave); // W/m2
					dailyLongWave = interp_day_sunshine(n, day, climate.G_longwave); // W/m2
                                        // here, dailyLongWave is net longwave.
					break;
			case 3:	dailyShortWave = interp_day_sunshine(n, day, climate.G_shortwave);    // W/m2
					dailyLongWave = interp_day_sunshine(n, day, climate.G_longwave);    // W/m2
                                        // here, dailyLongWave is incoming (downward) longwave radiation.
					break;
			case 4:	dailyShortWave = interp_day_sunshine(n, day, climate.G_shortwave); // W/m2
					break;
			}

			switch (options.petOpt) {
			case 1:	dailyWindSpeed = interp_day_prec(n, day_in_month, month, climate.G_windspeed) / 10.0; // [m/s]
					dailyVapPres = interp_day_prec(n, day_in_month, month, climate.G_vappres) / 10.0; // [kPa]
					dailyTmin = interp_day_temp(n, day, climate.G_tmin) / 100.0; // [oC]
					break;
			case 2:	dailyWindSpeed = interp_day_prec(n, day_in_month, month, climate.G_windspeed) / 10.0; // [m/s]
					dailyVapPres = interp_day_prec(n, day_in_month, month, climate.G_vappres) / 10.0; // [kPa]
					dailyTmin = interp_day_temp(n, day, climate.G_tmin) / 100.0; // [oC]
					break;
			case 3:	dailyTempRange = interp_day_temp(n, day, climate.G_temprange)/ 100.0; // [oC]
					dailyTempRange = abs(dailyTempRange);
					break;
			case 4: // FP 2015-06 ATTENTION! CHECK NEEDED: Why is there no definition and no "break" for petOpt=4 "PET: combination of Turc and Doorenbos&Pruitt"?
			case 5:	dailyWindSpeed = interp_day_prec(n, day_in_month, month, climate.G_windspeed) / 10.0; // [m/s]
					dailyVapPres = interp_day_prec(n, day_in_month, month, climate.G_vappres) / 10.0; // [kPa]
					break;
			case 6:	dailyVapPres = interp_day_prec(n, day_in_month, month, climate.G_vappres) / 10.0; //[kPa]
					break;
			}

			// generate synthetic daily precipitation
			// (or interpolate it, depending on option)
			if ((options.time_series == 4) || (options.time_series == 5) || (options.time_series == 6)) { // new time series input option - if somehow monthly precip data should be used
				const short int nDaysPerMonth[12] = { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };
				if (0 == options.raindays)
					dailyPrec = double(interp_day_prec(n, day_in_month, month, climate.G_precipitation));
				else {
					numberOfRainDays = static_cast<int> (climate.G_raindays[n][month]);
					monthlyPrecipitation = (double) climate.G_precipitation[n][month];

					dailyPrec = double(dailyPrecNrd.getDailyPrec(day_in_month,
							month, numberOfRainDays, (short) monthlyPrecipitation));
				}
			}
			// end if synthetic daily precipitation
		}
        // endif options (daily data interpolation from monthly data)
        // end of daily data interpolation

		// Cell-specific calibration parameters - Apply multiplier  // FP
		dailyPrec = calibParam.getValue(M_PREC,n) * dailyPrec;

		// calculate saturated vapour pressure
		double temp2 = dailyTempC + 237.3;
		double e_s = 0.6108 * exp(17.27 * dailyTempC / temp2);
		// [kPa] (eq. 4.2.2 Shuttleworth)
		//e_s should be calculated from tmin and tmax /2! NEW mw
		//e_s = (((0.6108 * exp(17.27 * dailyTmin / (dailyTmin+237.3)))+(0.6108 * exp(17.27 * dailyTmax / (dailyTmax+237.3))))/2.0);

		// cloudiness coefficients and alpha are different in
		// humid and arid areas. also define arid_gw for use in ARID-HUMID GROUNDWATER.
		// we assume that steppe, scrubland, savanna and hot desert
		// are arid areas (with relative humidity < 60%)
		// based on a comparison with a global map
		// of relative humidity
		double alpha;
		double a_c, b_c; // cloudiness coefficients
		double maxDailyPET;// max. daily pot. evap.

		// definition with use of relative humidity
		if (6 == options.petOpt) {
			if (1 == options.clclOpt)
				cerr<< "PET-method 6: Definitions of 'arid-/humid-'parameters can not be performed according to KOEPPEN's map!!!  \n";

			double RH;
			RH = (dailyVapPres / e_s * 100.0); //RH = relative humidity in [%]
			if (RH > 100) RH = 100;
			if (RH < 0)	RH = 0;

			if (RH < 60) { // arid area
				// Cell-specific calibration parameters - Use parameter  // FP
				alpha = calibParam.getValue(P_PTC_ARI,n);
				// Cell-specific calibration parameters - Use parameter  // FP
				// maxDailyPET: No distinction anymore arid vs. humid
				maxDailyPET = calibParam.getValue(P_PET_MXDY,n);
				a_c = a_c_arid;
				b_c = b_c_arid;
				arid_gw = true;
			} else {
				// humid area
				// Cell-specific calibration parameters - Use parameter  // FP
				alpha = calibParam.getValue(P_PTC_HUM,n);
				// Cell-specific calibration parameters - Use parameter  // FP
				// maxDailyPET: No distinction anymore arid vs. humid
				maxDailyPET = calibParam.getValue(P_PET_MXDY,n);
				a_c = a_c_humid;
				b_c = b_c_humid;
			}
		}
		// definition with use of Koeppen-climate-classification
		else if (1 == options.clclOpt) {
			if ((clcl.cls[n] == 13) || (clcl.cls[n] == 20) || (clcl.cls[n] == 25)) { // arid area
				alpha = clcl.clcl_alpha[n]; //based on Class A Pan evaporation measurements
				// Cell-specific calibration parameters - Use parameter  // FP
				// maxDailyPET: No distinction anymore arid vs. humid
				maxDailyPET = calibParam.getValue(P_PET_MXDY,n);
				a_c = a_c_arid;
				b_c = b_c_arid;
				arid_gw = true;
			} else { // humid area
				// clcl-classes 11+12 should be treated as arid concerning max. daily pot. evap. 
				if ((clcl.cls[n] == 11) || (clcl.cls[n] == 12))
					// Cell-specific calibration parameters - Use parameter  // FP
					// maxDailyPET: No distinction anymore arid vs. humid
					maxDailyPET = calibParam.getValue(P_PET_MXDY,n);
				else
					// Cell-specific calibration parameters - Use parameter  // FP
					// maxDailyPET: No distinction anymore arid vs. humid
					maxDailyPET = calibParam.getValue(P_PET_MXDY,n);

				alpha = clcl.clcl_alpha[n];  //based on Class A Pan evaporation measurements
				a_c = a_c_humid;
				b_c = b_c_humid;
			}
		} else {
		// pre-defined arid-humid areas after G_ARID_HUMID.UNF2
			switch(G_aindex[n]){
				case 1: //arid area
						// Cell-specific calibration parameters - Use parameter  // FP
						alpha = calibParam.getValue(P_PTC_ARI,n);
						// Cell-specific calibration parameters - Use parameter  // FP
						// maxDailyPET: No distinction anymore arid vs. humid
						maxDailyPET = calibParam.getValue(P_PET_MXDY,n);
						a_c = a_c_arid;
						b_c = b_c_arid;
						arid_gw = true;
						break;
				case 0: //humid area
						// Cell-specific calibration parameters - Use parameter  // FP
						alpha = calibParam.getValue(P_PTC_HUM,n);
						// Cell-specific calibration parameters - Use parameter  // FP
						// maxDailyPET: No distinction anymore arid vs. humid
						maxDailyPET = calibParam.getValue(P_PET_MXDY,n);
						a_c = a_c_humid;
						b_c = b_c_humid;
						break;
				default:
						cerr<<"Error: Invalid value for Arid/humid index: "<< G_aindex[n] <<endl;
						exit(1);
			}
		}
		// end final else (cloudiness coefficients and alpha ...)

		// LAI -------------------------------------
		// use dailyPrec instead of monthly averaged (daily) input data
		//dailyLai = lai.getDailyLai(n, land.G_landCover[n], arid_gw, dailyTempC, (double) climate.G_precipitation[n][month] / number_of_days_in_month[month], dailyPET);
		// use threshold of precipitation sum instead of dailyPET for decision of growing conditions
		//dailyLai = lai.getDailyLai(n, land.G_landCover[n], arid_gw, dailyTempC, dailyPrec, dailyPET);
		dailyLai = lai.getDailyLai(n, land.G_landCover[n], arid_gw, dailyTempC,	dailyPrec);
		dailyKc = lai.getDailyKc(n, land.G_landCover[n], dailyLai);
		// end LAI ---------------------------------

		// albedo ----------------------------------
		// net short wave radiation depends on albedo
		// a different calculation will be done for
		// different surfaces (with different albedo)
		//
		// added for innerannual variation in albedo
		// for crop/land-cover specific evapotranspiration with kc_pet-method
		// fully snow covered when swe > 3mm => use var. snow albedo
		if (G_snow[n] > 3.)
			albedo = albedoSnow_lct[land.G_landCover[n] - 1];
		else {
			if (options.calc_albedo) {
				// use given albedo in the first step
				albedo = G_AlbedoCycle[n][month];
			} else {
				if (options.use_kc == 1)
					albedo = 0.23; // use gras-reference albedo
				else
					albedo = albedo_lct[land.G_landCover[n] - 1]; // use land cover-specific albedo given in LCT.DAT
			}
		}

		// end albedo ------------------------------

		// latent heat of evaporation
		if (dailyTempC > 0)
			// latent heat of vaporization of water
			lat_heat = 2.501 - 0.002361 * dailyTempC; // [MJ/kg]
		else
			// latent heat of sublimation
			lat_heat = 2.835; // 2.501 + 0.334

                // calculate extra terrestrial solar radiation [mm/day]
		double ext_rad;
		ext_rad = calc_ext_rad(day, n); //Anhang A.2 Dis Kaspar

		double solar_rad_0; // max. solar radiation
		solar_rad_0 = (a_s + b_s) * ext_rad; // mm/day =S0

		// radiation values are in W/m2 mm/day
		// and have to be converted to mm/day :
		// W/m2 = J/(m2*s)
		// *86400/1000000     -> MJ/(m2*day)
		// /latent heat       -> kg/(m2*day)
		// /(1000*dens_water) -> m3/(m2*day)
		// *1000              -> mm/day
		// with dens_water = 1 -> "/1000000.0"
		// to rationalise computing:  (86400.0 / 1000000.0) = 0.0864
		double conv_Wm2_to_mmd = 0.0864 / lat_heat; // about 0.034 for T = 0C and 0.0345 for T = 10C
        double conv_mmd_to_Wm2 = lat_heat / 0.0864; // HMS 2013-11-21 to convert mm/d back into W/m2 for comparison of modeled radiation outputs
		// calculate solar radiation and afterwards net radiation

		//Observed radiation
		if ((options.cloud == 2) || (options.cloud == 3) || (options.cloud == 4)) {
			// solar radiation from data set
			solar_rad = conv_Wm2_to_mmd * dailyShortWave; // [mm/day]

			if (options.cloud == 2) {
					// net long wave [mm/day], net long wave input
					net_long_wave_rad = conv_Wm2_to_mmd * dailyLongWave; // [mm/day]
			}
			if (options.cloud == 3) {
					// .365 (H08) or .31 (WFD) inc. long wave input
					long_wave_rad_in = conv_Wm2_to_mmd * dailyLongWave; // unit: mm/d
                                        //const double emissivity = 0.98; //value for long wave radiation taken from Wilber et al. 1999, NASA TP-1999-209362 // formerly for H08
                                        double emissivity = emissivity_lct[land.G_landCover[n] - 1]; // land use class dependent emissivity
					double temp_K = dailyTempC + 273.2; // [K]
					const double stefan_boltz_const = 0.000000004903; // MJ /(m2 * K4 * day)
					long_wave_rad_out = emissivity * stefan_boltz_const * pow(temp_K, 4.) / lat_heat; // unit: mm/d
					// subtract outgoing from incoming longwave radiation (at earth surface)
					net_long_wave_rad = long_wave_rad_in - long_wave_rad_out; // unit: mm/d
			}
		}

		//Calculate global radiation from sunshine percentage or cloudiness
		if ((options.cloud == 0) || (options.cloud == 1)) {

			solar_rad = (a_s + b_s * dailySunshine / 1000.0) * ext_rad; //kaspar gl (2.11) =St
			// total incoming short-wave radiation: [mm/day]
			// the factor of 1000.0 is used because the unit
			// of sunshine is 10*percent
		}

		//Calculate net longwave radiation
		if ((options.cloud != 2) && (options.cloud != 3)) {
			// adjustment for cloud cover
			double solar_fraction;
			if (ext_rad == 0.)
				solar_fraction = 0.;
			else
				solar_fraction = solar_rad / solar_rad_0;
			// this case can happen if solar_rad is derived
			// from measured data and solar_rad_0 is calculated
			if (solar_fraction > 1.)
				solar_fraction = 1.;

			double ff = a_c * solar_fraction + b_c;
			// net emissivity between the atmosphere and the ground
			double net_emissivity = -0.02 + 0.261 * exp(-0.000777 * dailyTempC * dailyTempC);
			double temp_K = dailyTempC + 273.2; // [K]
			const double stefan_boltz_const = 0.000000004903; // MJ /(m2 * K4 * day)
			net_long_wave_rad = (-ff * net_emissivity * stefan_boltz_const * pow(temp_K, 4.)) / lat_heat;
		}

		// net short wave [mm/day]  (gl 2.10 dis kaspar)
		net_short_wave_rad = solar_rad * (1. - albedo);

        // reflected solar radiation [mm/day] HMS radiation output adjustment
        refl_rad = solar_rad - net_short_wave_rad;
    
		// total net radiation [mm/day]
		// Cell-specific calibration parameters - Apply multiplier  // FP
		net_rad = calibParam.getValue(M_NETRAD,n) * (net_short_wave_rad + net_long_wave_rad);

		// open water surfaces have different albedo (0.08)
		// than the cell-specific land cover
		openWaterNetShortWaveRad = solar_rad * (1. - openWaterAlbedo);
		openWaterNetRad = openWaterNetShortWaveRad + net_long_wave_rad; // [mm/day]


		/////////////////////////////////// PET /////////////////////////////////////////////////////////////////////////////////////
		// PET-calculations: option 0-6
		// 0: PT landbed -old
		// 1: PM reference crop -old
		// 2: PM altitude-dependent pressure,e_s with temp range -old
		// 3: Hargreaves -new
		// 4: Turc/doorenbos -new
		// 5: Kimberly-Penman -new
		// 6: PT RH -new
		// begin of calculations for PET
		double dailyPET = 0.; // [mm]
		double dailyOpenWaterPET = 0.; // [mm]


		// psychrometric constant [kPa/oC]
		double atmos_pres = 101.3; // [kPa] atmospheric pressure
		double gamma, gammaStar;
		double c3;
		// increase of saturated vapor pressure with temperature
		double inc_svp = 4098. * e_s / (temp2 * temp2); //sollte man dann hier tmin+tmax/2 nehmen??
		// [kPa/oC]  (eq. 4.2.3 Shuttleworth)


		if (2 == options.petOpt)
			atmos_pres = 101.3 * pow(((293.0 - 0.0065 * land.G_altitude[n]) / 293.0), 5.26);
		c3 = 0.0016286 * atmos_pres;

		gamma = c3 / lat_heat;
		// eq. 4.2.32 Shuttleworth
		// values modified according to FAO 56
		// original value in Shuttleworth 0.33
		gammaStar = gamma * (1 + 0.34 * dailyWindSpeed);// grass
		//gammaStar = gamma * (1 + 0.38 * dailyWindSpeed);// alfalfa

		//
		// calculate potential evapotranspiration
		//
		// according to procedure proposed in the FAO-Report 56
		// 1. first calculate gras-reference evapotranspiration;
		// 2. calcuate crop/land-cover specific evapotranspiration

		// REFERENCE EVAPOTRANSPIRATION
		//
		// Priestley-Taylor
		// according to Shuttleworth (1993)
		if ((0 == options.petOpt) || (6 == options.petOpt)) {
			if (net_rad <= 0.)
				dailyPET = 0.;
			else
				dailyPET = alpha * (inc_svp * net_rad) / (inc_svp + gamma); // [mm/day]

			// calculation of PET for open water surfaces (lakes/wetlands)
			if (openWaterNetRad <= 0.)
				dailyOpenWaterPET = 0.;
			else
				dailyOpenWaterPET = alpha * (inc_svp * openWaterNetRad)	/ (inc_svp + gamma);// [mm/day]
		}

		// Penman-Monteith
		if ((1 == options.petOpt) || (2 == options.petOpt)) {

			// vapour pressure deficit
			double d_vp;
			double e_a;
			double Tdew;

			/*
			 d_vp = e_s - dailyVapPres;
			 if (d_vp < 0)
			 d_vp = 0;
			 */
			//I believe CRU vapour pressure data not to be reliable -mw-
			//therefore calculation acc to "an update for the definition of reference evapotranpiration" Allen, 1994
			//distinguish arid/humid
			// Cell-specific calibration parameters - Use parameter  // FP
			if (alpha == calibParam.getValue(P_PTC_ARI,n)) {
				Tdew = (dailyTmin - 3.0); //arid
			} else
				Tdew = dailyTmin; //humid

			e_a = (0.6108 * exp(17.27 * Tdew / (Tdew + 237.3)));
			d_vp = e_s - e_a;

			if (net_rad <= 0.)
				net_rad = 0.;

			// eq. 4.2.31 Shuttleworth
			// values modified according to FAO 56
			// original value in Shuttleworth 275
			dailyPET = ((inc_svp * net_rad) + (900.0 * gamma * dailyWindSpeed
					* d_vp) //grass
					//			+ (1600.0 * gamma * dailyWindSpeed * d_vp) // alfalfa
					/ (dailyTempC + 273.0)) / (inc_svp + gammaStar); // [mm/day]

			// calculation of PET for open water surfaces (lakes/wetlands)
			// (eq. 4.2.30 Shuttleworth)
			if (openWaterNetRad <= 0.)
				openWaterNetRad = 0.;
			dailyOpenWaterPET = ( (inc_svp * openWaterNetRad) + (gamma * 6.43 * (1 + 0.536 * dailyWindSpeed) * d_vp) / lat_heat)
								  / (inc_svp + gamma);	// [mm/day]
		}
		// Hargreaves
		if (3 == options.petOpt) {
			//Shuttleworth eq.4.2.44 (hat wurzel vergessen )
			if (ext_rad <= 0.) {
				dailyPET = 0.;
			}
			else {
				dailyPET = (0.0023 * ext_rad * (pow(dailyTempRange,0.5)) * (dailyTempC + 17.811)); //[mm/day]
				dailyPET = abs(dailyPET);
                if (std::isnan(dailyPET)) { // Felix Portmann for C++11 because of ambiguous call, result in C: int, C++11: boolean
				printf (" ext_rad %f, temprange %f, tempC%f \n", ext_rad, dailyTempRange, dailyTempC);
				}
			}

			// calculation of PET for open water surfaces (lakes/wetlands)
			// nach Priestley Taylor
			if (openWaterNetRad <= 0.)
				dailyOpenWaterPET = 0.;
			else {
				dailyOpenWaterPET = alpha * (inc_svp * openWaterNetRad) / (inc_svp + gamma);	// [mm/day]
			}

		}
		//Combination of Turc ivanov for humid areas and Doorenbos & Pruitt for arid areas
		if (4 == options.petOpt) {
			//Shuttleworth 4.2.40-4.2.42
			double RH;
			RH = (dailyVapPres / e_s * 100.0);   //RH = relative humidity in [%]
			if (RH > 100.) RH = 100.;
			if (RH < 0.) RH = 0.;

			// Cell-specific calibration parameters - Use parameter  // FP
			// Humid conditions
			if (alpha == calibParam.getValue(P_PTC_HUM,n)) {

				if (dailyTempC <= 0.) { //turc nur fuer t>0, bei T< 0 ivanov
					/*if (dailyTempC <= -25.0){
					 dailyPET = 0;
					 //printf("ivanov <-25\n");
					 }
					 else{
					 dailyPET = 0.000036 * (dailyTempC + 25.0) * (dailyTempC + 25.0) * (100.1 - RH);
					 //printf("ivanov\n");
					 }
					 //if (std::isinf( dailyPET )) {
					 //	printf("ivanov: tempC %f, RH %f\n", dailyTempC, RH);
					 //}
					 */
					dailyPET = 0.;

				} else { //temp>0
					if (net_short_wave_rad <= 0.)
						dailyPET = 0.;
					else {
						if (RH <= 50.0) { //Fallunterscheidung fuer Turc
							dailyPET = (0.31 * (dailyTempC / (dailyTempC + 15.0)) * (net_short_wave_rad + 2.09) * (1.0 + ((50.001 - (RH)) / 70.0)));
							//printf("turc rh <50 \n");
						} else {
							dailyPET = (0.31 * (dailyTempC
									/ (dailyTempC + 15.0))
									* (net_short_wave_rad + 2.09));
							printf("turc rh> 50 \n");
						}

					}
				}
			} else { // alpha == alphaArid -> Doorenbos
				if (net_short_wave_rad <= 0.)
					net_short_wave_rad = 0.;

				double bdp;
				bdp = (1.066 - 0.0013 * RH + 0.045 * dailyWindSpeed - 0.0002
						* RH * dailyWindSpeed - 0.000315 * pow(RH, 2) - 0.0011
						* pow(dailyWindSpeed, 2));

				dailyPET = (-0.3 + bdp * (inc_svp  / (inc_svp + gamma)) * net_short_wave_rad );
				//printf("doorenb\n");
                if (std::isinf(dailyPET)) { // Felix Portmann for C++11 because of ambiguous call, result in C: int, C++11: boolean
					printf(" Doorenbos: bdp %f, RH %f\n", bdp, RH);
				}
			}

			// OPEN WATER//
			if (openWaterNetShortWaveRad <= 0.)
				dailyOpenWaterPET = 0.;

			// Cell-specific calibration parameters - Use parameter  // FP
			// Humid conditions
			if (alpha == calibParam.getValue(P_PTC_HUM,n)) {
				if (RH <= 50.0){  //Fallunterscheidung fuer Turc
					dailyOpenWaterPET = (0.31 * (dailyTempC / (dailyTempC + 15.0)) * (openWaterNetShortWaveRad + 2.09) * (1.0 + ((50.0 - (RH)) / 70.0)));
				} else {
					dailyOpenWaterPET = (0.31 * (dailyTempC / (dailyTempC + 15.0)) * (openWaterNetShortWaveRad + 2.09));
				}
			} else { // alpha == alphaArid
				double bdp;
				bdp = (1.066 - 0.0013 * RH + 0.045 * dailyWindSpeed - 0.0002 * RH * dailyWindSpeed - 0.000315 * pow(RH,2) - 0.0011 * pow(dailyWindSpeed,2));
				dailyOpenWaterPET = (-0.3 + bdp * (inc_svp  / (inc_svp + gamma)) * openWaterNetShortWaveRad );
				}

		}
		//Kimberly Penman
		if (5 == options.petOpt) {
			//Shuttleworth eq.4.2.33
			if (net_rad <= 0.)
				dailyPET = 0.;
			else {

				double d_vp;// vapour pressure deficit
				d_vp = e_s - dailyVapPres;
				if (d_vp < 0) d_vp = 0.;

				double Wf;
				double awf;
				double bwf;
				short j_day;
				double j_a;
				double j_b;
				// Fallunterscheidung fuer northern and southern hemisphere
				if (-(geo.G_row[n] / 2.0 - 90.25) > 0.0) { //Umrechnung von Reihe auf Latitude in Grad; northern hem positiv
					j_day = day;
				} else { //southern lat
					if (day >= 182) {
						j_day = (day - 182);
					} else {
						j_day = (day + 183);
					}
				}

				j_a = (-(((j_day - 173.) / 58.) * ((j_day - 173.) / 58.)));
				j_b = (-(((j_day - 243.) / 80.) * ((j_day - 243.) / 80.)));
				awf = (0.4 + 1.4 * exp(j_a));
				// bwf = (0.605 + 0.345 * exp(j_b));// as in shuttleworth and asce
				bwf = ((0.007 + 0.004 * exp(j_b)) * 86.4); // agronomy journal better values in non-growing season
				Wf = (awf + bwf * dailyWindSpeed);
	
				dailyPET = (((inc_svp * net_rad) + (gamma * 6.43 * Wf * d_vp ) / lat_heat) / (inc_svp + gamma));	// [mm/day]
				
			}

			// OPEN WATER//
			// calculation of PET for open water surfaces (lakes/wetlands)
			// nach Priestley Taylor
			if (openWaterNetRad <= 0.)
				dailyOpenWaterPET = 0.;
			else {
				dailyOpenWaterPET = alpha * (inc_svp * openWaterNetRad)
						/ (inc_svp + gamma); // [mm/day]
			}

		}

		// crop/land-cover specific evapotranspiration using crop coefficients
		if ((G_snow[n] <= 3.) && (options.use_kc == 1)) {
			dailyPET *= dailyKc;
			dailyOpenWaterPET *= kc_OpenWater;
		}

        G_lakeBalance[n] = (dailyPrec - dailyOpenWaterPET) * G_cellCorrFact[n]; // HMS 2014-06-04 reintroduced

		// added for open water evaporation reduction (2.1f)
		G_openWaterPrec[n] = dailyPrec;
		G_openWaterPET[n] = dailyOpenWaterPET;


		////////// the following calculations do not need to be done for lake cells,		///////////////////////////////////////////////////
		////////// but for monthly summation we must not miss the calculations at the end...///////////////////////////////////////////////////

		// variables added to sum up change in canopy, snow, soil and groundwater storage
		// used for cell AET calculation (2.1f)
		double landStorageChangeSum = 0.;
        double initialStorage = 0.;
		// calculate interception
		double max_canopy_storage  = 0.;
		double dailyCanopyEvapo		= 0.;
		double daily_prec_to_soil	= 0.;
		double dailySoilPET			= 0.;
		double canopy_water_content	= 0.;
		// calculate snowpack
		double snowmelt				= 0.;
		double dailySnowEvapo		= 0.;
		double dailyEffPrec			= 0.;
		// calculate immediate (urban) runoff
		double immediate_runoff 	= 0.;
		// calculate surface runoff
		double surface_runoff 		= 0.;
		// calculate actual evapotranspiration
		double dailyAET 			= 0.;
		double daily_runoff 		= 0.;
		double total_daily_runoff 	= 0.;
		// calculate groundwater recharge
		double daily_gw_recharge 	= 0.;
		double pot_gw_recharge 	= 0.;
		// calculate groundwater runoff
        // double groundwater_runoff 	= 0.; // now in routing.cpp

/////////////////////!!!!!!!!!!!!! calculations only for land cells !!!!!!!!!!!!!!!!!////////////////////////////////////
        //if (n < ng_land) {
        // ! caluclations now for all cells (to avoid water balance errors in lake cells which are not (yet) anymore lake cells)
        if (n < ng) {

			// Efficiency enhancement through local copies instead of often used method getValue // FP
			double P_T_SNOWFZ__at__n = calibParam.getValue(P_T_SNOWFZ,n);
			double P_T_SNOWMT__at__n = calibParam.getValue(P_T_SNOWMT,n);
			double P_T_GRADNT__at__n = calibParam.getValue(P_T_GRADNT,n);
			double M_DEGDAY_F__at__n = calibParam.getValue(M_DEGDAY_F,n);

            double G_landAreaFracPrevTimestep__at__n = routing.G_landAreaFracPrevTimestep[n]; // FP20161005N001
///////////////////////////////////////// interception //////////////////////////////////////////////////////////////////

			if (1 == options.intercept) {
				double canopy_deficiency;
                // adapt canopy water content [mm] to changed land area frac // HMS 2015-04-29
                // if landAreaFrac == 0; add canopy water content to runoff and set it ot 0
                if (landAreaFrac == 0.) {
					
					G_dailyStorageTransfer[n] = G_canopyWaterContent[n];
                    G_canopyWaterContent[n] = 0.;
					dailyCanopyEvapo = 0.;
                }
                else { // adapt for changed land area frac and calculate canopy water content
//                    G_canopyWaterContent[n] *= routing.G_landAreaFracPrevTimestep[n]/ landAreaFrac; // obsolete
                    G_canopyWaterContent[n] *= G_landAreaFracPrevTimestep__at__n / landAreaFrac; // FP20161005N001

                    initialStorage = G_canopyWaterContent[n]; // calculation of cell AET (2.1f)

                    // calculation of dailyLAI has been moved before calculations of albedo starts
                    if (dailyLai > 0.00001) {
					   // double canopy_deficiency; // CR: moved to top to enable "cout << canopy_deficiency << endl;"
						// Cell-specific calibration parameters - Use parameter  // FP
						max_canopy_storage = calibParam.getValue(P_MCWH,n) * dailyLai;	// [mm]
                        canopy_deficiency = max_canopy_storage - G_canopyWaterContent[n];
                        if (dailyPrec < canopy_deficiency) {
                            G_canopyWaterContent[n] += dailyPrec;
                            daily_prec_to_soil = 0.;
                        } else {
                            G_canopyWaterContent[n] = max_canopy_storage;
                            daily_prec_to_soil = dailyPrec - canopy_deficiency;
                        }
                        // end treatment of canopy deficiency

                        canopy_water_content = G_canopyWaterContent[n];
                        dailyCanopyEvapo = dailyPET
                            * pow((canopy_water_content / max_canopy_storage), canopyEvapoExp); // canopyEvapoExp = 2/3

                        if (dailyCanopyEvapo > canopy_water_content) {
						
                            // All the water in the canopy is evaporated.
                            // dailyCanopyEvapo has to be reduced, because
                            // part of the energy is left and can lead to
                            // additional evapotranspiration from soil
                            // later in the program.
                            dailyCanopyEvapo = canopy_water_content;
                            dailySoilPET = dailyPET - canopy_water_content;
                            G_canopyWaterContent[n] = 0.0;
                        } else {
                            G_canopyWaterContent[n] -= dailyCanopyEvapo;
                            dailySoilPET = dailyPET - dailyCanopyEvapo;
                        }
                        // end treatment of dailyCanopyEvapo


                    } else {
                        // LAI == 0
                        daily_prec_to_soil = dailyPrec;
                        dailySoilPET = dailyPET;
                        dailyCanopyEvapo = 0.0;
                    }
                    // endelse
                    // end treatment of LAI


                    // added for cell AET calculation (2.1f)
                    landStorageChangeSum += G_canopyWaterContent[n] - initialStorage;

                }// endelse landAreaFrac > 0

            }
            else {
				// no interception
				daily_prec_to_soil = dailyPrec;
				dailySoilPET = dailyPET;
				dailyCanopyEvapo = 0.0;
			}
            // endelse (no interception)
            // end treatment of interception

			if (dailySoilPET < 0.)
				dailySoilPET = 0.0;

				
/////////////////////////////////////////////// snow ///////////////////////////////////////////////////////////////////////////////
			//
			// snowpack
			//
			double dailyEffPrecBeforeSnowMelt_elev;
			double dailyEffPrec_elev;
			double snowmelt_elev;
			double temp_elev; // Temperature at each elevation step
			double daily_snow_to_soil = 0.;
			double daily_snow_to_soil_elev;
			double dailySnowFall_elev;
            double dailyRainFall_elev = 0.;  // FP20161011N001
			double TempElevMax = 0.;
			// added for cell AET calculation (2.1f)
			double snowStorageChange 	= 0.;
	
			dailySnowFall 	= 0.;
			thresh_elev[n] 	= 0;
			G_snow[n] 		= 0.;
	
			//loop for all subgrids in cells with elevation >0m. At the end of the loop, snow cover of
			//all subgrids are added to the 0.5° Grid (G_Snow) again.
	
			for (short elev = 1; elev < 101; elev++) {	//count the subgrids
	
				temp_elev           = 0.;
				dailySnowFall_elev  = 0.;
                dailyRainFall_elev  = 0;  // FP20161011N001
				daily_snow_to_soil_elev = 0.;
				dailyEffPrec_elev   = 0.;
				snowmelt_elev       = 0.;
				dailyEffPrecBeforeSnowMelt_elev = 0.;

				// elevation dependent temperature in one grid cell
				// "0.006" = 0.6°C/100m after Semadeni-Davies (1997) and Dunn, Colohan (1999);
				// The first row of G_ELEV_RANGE.UNF2 (G_Elevation[n][0]) contains mean elevation;		
				// Cell-specific calibration parameters - Use parameter  // FP

				temp_elev = dailyTempC - ((G_Elevation[n][elev]	- G_Elevation[n][0]) * P_T_GRADNT__at__n);

	
				
				// adapt snow depth to changed land area fraction
				// and set snow arrays to 0 in case of land area fraction = 0.// HMS 2015-04-29
				if (landAreaFrac == 0.) {
					G_dailyStorageTransfer[n] += G_SnowInElevation[n][elev] / 100.;
					G_SnowInElevation[n][elev] = 0.;
					G_snow[n] = 0.;
				} 
				else {	// landAreaFrac > 0.
//					G_SnowInElevation[n][elev] *= routing.G_landAreaFracPrevTimestep[n]/ landAreaFrac; // obsolete
                    G_SnowInElevation[n][elev] *= G_landAreaFracPrevTimestep__at__n / landAreaFrac; // FP20161005N001

					// added for cell AET calculation (2.1f)
					initialStorage = G_SnowInElevation[n][elev];
					
					//-------------------------------------------------------------
					// If snow in subgrid exceeds 1000mm SWE, temperature doesn't decrease in upper subgrids.
					// This should prevent uncontrolled snow accumulation.		

					if (G_SnowInElevation[n][elev] > 1000.) {

						// define threshold elevation in cell
						if (thresh_elev[n] == 0.) {
							// remember threshold elevation of cell
							thresh_elev[n] = G_Elevation[n][elev];
						}
						// cell above threshold elevation
						else if (thresh_elev[n] > 0.) {
							//temp_elev = 0.;	// set to 0, because temp_elev will be recalculated
							// all upper elevations get same temperature calculated with remebered elev
							// Cell-specific calibration parameters - Use parameter  // FP
							temp_elev = dailyTempC - ((thresh_elev[n] - G_Elevation[n][0]) * P_T_GRADNT__at__n);
						}

						else {
							cerr << " Negative threshold elevation (thresh_elev[n])  \n";
						}

					}


					// added for cell AET calculation (2.1f)
					//CR 2015-08: initialStorage now moved to the beginning of snow storage calculation 
					//initialStorage = G_SnowInElevation[n][elev];

					// accumulation of new snow

					// Cell-specific calibration parameters - Use parameter  // FP
					if (temp_elev <= P_T_SNOWFZ__at__n) {
						dailySnowFall_elev = dailyPrec; // value above canopy, not used, only for output
						dailyEffPrecBeforeSnowMelt_elev = 0.;

						daily_snow_to_soil_elev = daily_prec_to_soil;
						G_SnowInElevation[n][elev] += daily_snow_to_soil_elev; //value below canopy

						if (G_SnowInElevation[n][elev] > dailySoilPET) {
							G_SnowInElevation[n][elev] -= dailySoilPET;
							dailySnowEvapo += dailySoilPET;
						} else {
							dailySnowEvapo += G_SnowInElevation[n][elev];
							G_SnowInElevation[n][elev] = 0.;
						}
					} else {
						dailyEffPrecBeforeSnowMelt_elev = daily_prec_to_soil;
                        dailyRainFall_elev = dailyPrec; // value above canopy, not used, only for output  // FP20161011N001
						// dailySnowEvapo += 0;
					}

					// melting of snow
					// Cell-specific calibration parameters - Use parameter  // FP
					if (temp_elev > P_T_SNOWMT__at__n) {

						if (G_SnowInElevation[n][elev] < 0.)
							cerr << "G_SnowInElevation[n][elev] < 0 \n";
						else {	//snowmelt with degree-day factor
							// Cell-specific calibration parameters - Apply multiplier & Use parameter  // FP
							snowmelt_elev = M_DEGDAY_F__at__n * ddf_lct[land.G_landCover[n] - 1] * (temp_elev - P_T_SNOWMT__at__n);
		
							if (snowmelt_elev > G_SnowInElevation[n][elev]) {
								snowmelt_elev = G_SnowInElevation[n][elev];
								G_SnowInElevation[n][elev] = 0.;
							} else {
								G_SnowInElevation[n][elev] -= snowmelt_elev;
							}
						}

					}
					// end of snowmelt
					
					dailyEffPrec_elev = dailyEffPrecBeforeSnowMelt_elev	+ snowmelt_elev;

					// added for cell AET calculation (2.1f)
					snowStorageChange += G_SnowInElevation[n][elev]	- initialStorage;

					// the information about the maximum temperature at the gridcell is used in AET calculations below
					if (elev == 1) {
						TempElevMax = temp_elev;
					}

					if (G_SnowInElevation[n][elev] > 0.)
					
						dailySnowCoverFrac += 1. / 100.; // sum up all snow covered subgrids to get the total fraction in cell
					G_snow[n] += G_SnowInElevation[n][elev]; // sum up all snow in subgrids to G_Snow
					dailyEffPrec += dailyEffPrec_elev; // sum up all precipitation in subgrids to grid precip.
					dailySnowFall += dailySnowFall_elev; // sum up all snowfall in subgrids to grid snowfall
                    dailyRainFall += dailyRainFall_elev; // sum up all rainfall in subgrids to grid rainfall // FP20161011N001
					snowmelt += snowmelt_elev;
					daily_snow_to_soil += daily_snow_to_soil_elev;
					
				} // end if landAreaFrac > 0
			}	// end subgrids (for)
	
				if (landAreaFrac > 0.) {
					G_snow[n] 			/= 100.;

					dailyEffPrec 		/= 100.;	// sum of all subgrids has to be divided by the number of land subgrids
					dailySnowFall 		/= 100.;	// within the cell (if only land-subgrids, value is 100).
                    dailyRainFall       /= 100.;  // FP20161011N001
					daily_snow_to_soil 	/= 100.;
					snowmelt 			/= 100.;
					// added by hunger (9/2006): in former snow pack version dailySnowEvapo of subgrids was not summed up
					dailySnowEvapo 		/= 100.;
			
					// added for cell AET calculation (2.1f)
					snowStorageChange 	/= 100.;
					landStorageChangeSum += snowStorageChange;
				}
					// end of snow calculations				

					
////////////////////////////////////////////////// immediate runoff ////////////////////////////////////////////////////////////////////////////////

			// over built-up areas 50 percent of precipitation is immediate runoff

			if (land.G_built_up[n] > 0.) {			
				immediate_runoff = runoffFracBuiltUp * dailyEffPrec * land.G_built_up[n];
				dailyEffPrec -= immediate_runoff;
			}

/////////////////////////////////////////////// SOIL AND AET ////////////////////////////////////////////////////////////////////////////////////////
			//
			// calculate actual evapotranspiration
			//

			// added for cell AET calculation (2.1f)
			// adapt soil water content [mm] to changed land area frac
			if (landAreaFrac == 0.) {
			
				// At this point, G_dailyStorageTransfer[n] already contains G_canopyWaterContent[n] and G_SnowInElevation[n].
				// if landAreaFrac == 0: G_dailyStorageTransfer[n] is used as local runoff in routing.cpp
				G_dailyStorageTransfer[n] += G_soilWaterContent[n];
				G_dailyStorageTransfer[n] *= G_cellCorrFact[n];
				G_soilWaterContent[n] = 0.;
				daily_gw_recharge = 0.;
				total_daily_runoff = 0.;
				G_dailyGwRecharge[n] = 0.;
				dailyCellLandAET[n] = 0.;
				dailyCellLandAET_uncorr[n] = 0.;
				total_daily_runoff = 0.;

			}
			else { // landAreaFrac > 0.

//				G_soilWaterContent[n] *= routing.G_landAreaFracPrevTimestep[n]/ landAreaFrac;  // obsolete
                G_soilWaterContent[n] *= G_landAreaFracPrevTimestep__at__n / landAreaFrac; // FP20161005N001

				initialStorage = G_soilWaterContent[n];

				// If scaling exceeds maximum soil storage,
				// let water drain as runoff
				// and set G_soilWaterContent[n] to maximum soil water content
				soil_water_overflow = 0;
				if (G_soilWaterContent[n] > maxSoilWaterCap.G_Smax[n]) {
					// calculate surplus
					soil_water_overflow = G_soilWaterContent[n] - maxSoilWaterCap.G_Smax[n];
					// set storage to maximum
					G_soilWaterContent[n] = maxSoilWaterCap.G_Smax[n];
					
					// CR 2015-08:
					// daily_runoff is re-calculated later (R=Peff *(S/Smax)^gamma) and cannot be used here.
					// Therefore, soil_water_overflow is increased by the "final" soil_water_overflow in the end.
					// daily_runoff += soil_water_overflow;
				}
				// end treatment of soil_water_overflow

				//CR 2015-08: initialStorage now moved to the beginning of soil water balance 
				//initialStorage = G_soilWaterContent[n];
				
				// was "if (dailyTempC > snowFreezeTemp)"
				// Changed to TempElevMax > snowFreezeTemp, because even if dailyTempC < 0, there might be runoff.
				// Now the Temperatur of the cell with lowest elevation is used as threshold.
				// This is caused by the subgrids in the snow algorithm, where some have a higher Temp
				// than the mean temperature of the cell and might produce runoff (here "dailyEffPrec").

				// Cell-specific calibration parameters - Use parameter  // FP
				if (TempElevMax > P_T_SNOWFZ__at__n) {

					//if (maxSoilWaterCap.G_Smax[n] > -1.) 
					if (maxSoilWaterCap.G_Smax[n] > 0.) {
					
						// CR ??? soil_saturation is computed in the following line anyway?
						// re-initialization
						//soil_saturation = 0;
						//CR 2015-08: soil_water_overflow is increased later and cannot be set to 0 here
						//soil_water_overflow = 0;
						
						soil_saturation = G_soilWaterContent[n] / maxSoilWaterCap.G_Smax[n];					
						daily_runoff = dailyEffPrec * pow(soil_saturation, (double) G_gammaHBV[n]);

						//check wether max daily PET should be limited or not (see maxDailyPET_arid or maxDailyPET_humid)
						if ( dailySoilPET > (maxDailyPET - dailyCanopyEvapo) * soil_saturation )
							dailyAET = (maxDailyPET - dailyCanopyEvapo)	* soil_saturation;
						else
							dailyAET = dailySoilPET;							

						G_soilWaterContent[n] += dailyEffPrec - dailyAET - daily_runoff;

						dailyEffPrec = 0.;

						if (G_soilWaterContent[n] < 0.) {

							// too much water has been taken out of the soil water storage:

							// correct dailyAET:
							// reduce it by G_soilWaterContent[n]
							// G_soilWaterContent[n] is negative! THerefore +=
							dailyAET += G_soilWaterContent[n];
							//soil water content is at most as negative as dailyAET!
							//if (dailyAET < 0.)
							//	dailyAET = 0.;

							// correct soil water storage
							G_soilWaterContent[n] = 0.;

						}

                        daily_runoff *= G_cellCorrFact[n];
                        immediate_runoff *= G_cellCorrFact[n];

						// as sublimation might have occured in parts of the cell snow evaporation
						// has to be added to dailyAET (hunger 9/2006)
						// dailyAET += dailySnowEvapo;

						//===============ARID GROUNDWATER===========================
						// arid gw calculations should be calculated depending on soil texture
                        // (semi)arid cells with texture class < 21

                        if ((arid_gw) && (GW.G_texture[n] < 21)) {
                            pot_gw_recharge = 0.;
                            if ((GW.getRgmax(n)/100.) < (GW.getgwFactor(n) * daily_runoff))
                                daily_gw_recharge = GW.getRgmax(n)/100.;
                            else
                                daily_gw_recharge = GW.getgwFactor(n) * daily_runoff;

                            if (dailyPrec <= calibParam.getValue(P_PCRITGWA,n)) {
                                pot_gw_recharge = daily_gw_recharge;
                                daily_gw_recharge = 0.;
                            }
                        }
                        else {
                            pot_gw_recharge = 0.;
                            if ((GW.getRgmax(n)/100.) < (GW.getgwFactor(n) * daily_runoff))
                                daily_gw_recharge = GW.getRgmax(n)/100.;
                            else
                                daily_gw_recharge = GW.getgwFactor(n) * daily_runoff;
                        }
                        // Removes doubling CFA (G_cellCorrFact[n]) from pot_gw_recharge
                        daily_runoff -= pot_gw_recharge;
                        pot_gw_recharge /= G_cellCorrFact[n];
                        G_soilWaterContent[n] += pot_gw_recharge ;

                        // If scaling exceeds maximum soil storage, let water drain as runoff
                        if (G_soilWaterContent[n] > maxSoilWaterCap.G_Smax[n]) {

                            //CR 2015-08: increase the value calculated in the beginning
                            soil_water_overflow += G_soilWaterContent[n] - maxSoilWaterCap.G_Smax[n];
                            G_soilWaterContent[n] = maxSoilWaterCap.G_Smax[n];
                           }


                        soil_water_overflow *= G_cellCorrFact[n]; // multiplied CFA to ensure CFA in all partions of total daily runoff
                        total_daily_runoff = daily_runoff + immediate_runoff + soil_water_overflow;

                    }
                        // end if (maxSoilWaterCap.G_Smax[n] > 0.)
                    else {
                        // no data for G_Smax (Greenland)
                        total_daily_runoff = 0.;
                        daily_gw_recharge = 0.;
                        //dailyAET = dailySnowEvapo;
                    }
                }
                    // end if (TempElevMax > snowFreezeTemp)
                    // all precipitation occurs as snow
                else {

                    // CR 2015-08: Even if (TempElevMax <= snowFreezeTemp), soil_water_overflow occurs as local runoff:
                    soil_water_overflow *= G_cellCorrFact[n];
                    // CR 2015-08: dailyEffPrec == 0 at this point if (TempElevMax <= snowFreezeTemp). In this case, dailyEffPrec has already been added to soilWaterContent.

                    dailyEffPrec *= G_cellCorrFact[n];
                    total_daily_runoff 	+= soil_water_overflow + dailyEffPrec;
                    daily_gw_recharge 	= 0.;
                    dailyAET = 0.;
                }

                G_dailyGwRecharge[n] = daily_gw_recharge; // to have gw_recharge in routing
			
				// added for cell AET calculation (2.1f)
				landStorageChangeSum += G_soilWaterContent[n] - initialStorage;

				// dailyCellLandAET is a corrected actual total evaporation (canopy, snow and soil)
				// it is consistent with cell-corrected runoff
				dailyCellLandAET[n] = landStorageChangeSum * (G_cellCorrFact[n] - 1.0)
								- dailyPrec * (G_cellCorrFact[n] - 1.0)
								+ (dailyAET + dailyCanopyEvapo + dailySnowEvapo) * G_cellCorrFact[n];
				dailyCellLandAET_uncorr[n] = (dailyAET + dailyCanopyEvapo + dailySnowEvapo);
			}
			// end if landAreaFrac > 0.

			
///////////////////////////////////////////// surface runoff /////////////////////////////////////////////////////////////
			// If landAreaFrac == 0: total_daily_runoff and daily_gw_recharge == 0.
			surface_runoff = total_daily_runoff - daily_gw_recharge;
            G_dailyLocalSurfaceRunoff[n] = surface_runoff;			
			

#ifdef _WATERGAP_CHECKS_GLOBAL_H
                // Check Felix Portmann 2015 - for invalid data
                if ( (year >= chk_year) && (month >= chk_month) && (day_in_month >= chk_day_in_month) ) {
                    check_value_double_LT0(total_daily_runoff);
                    check_value_double_LT0(daily_gw_recharge);
                    check_value_double_LT0(surface_runoff);

                    check_value_double_NaN(total_daily_runoff);
                    check_value_double_NaN(daily_gw_recharge);
                    check_value_double_NaN(surface_runoff);
                }
#endif

///////////////////////////////////////////// groundwater storage//////////////////////////////////////////////////////////
                        // groundwater_runoff is set to zero if groundwater storage drops below zero due to water withdrawals
            /*	if (G_groundwater[n] >= 0.0) // now in routing.cpp
				groundwater_runoff = k_g * G_groundwater[n]; // [mm]
			else
				groundwater_runoff = 0.0;
			G_groundwater[n] -= groundwater_runoff;
			G_dailyLocalGWRunoff[n] = groundwater_runoff;
			G_dailyLocalSurfaceRunoff[n] = surface_runoff;*/
		}
		//
		// end of calculations and store data
		//

		// monthly/yearly results
		if ((2 == options.grid_store) || (4 == options.grid_store) || (6 == options.grid_store)) {
			//G_monthlyPrecipitation[n][month] += ( dailyPrec * landAreaFrac / (((float)geo.G_landfreq[n] + (float)geo.G_fwaterfreq[n])/ 100.) );
			// use the whole amount of prec for all parts off the cell (land/open water)
			G_monthlyPrecipitation[n][month] += dailyPrec;
      		G_monthlyTemperature[n][month] += dailyTempC;
     		G_monthlySunshine[n][month] += dailySunshine / 10;	// sunshine percentage output [in %]
      		G_monthlyExtRad[n][month] += ext_rad * conv_mmd_to_Wm2; // HMS radiation output adjustment
            G_monthlyShortUpRad[n][month] += refl_rad * conv_mmd_to_Wm2; // HMS radiation output adjustment
            G_monthlyLongDownRad[n][month] += long_wave_rad_in * conv_mmd_to_Wm2; // HMS radiation output adjustment
            G_monthlyLongUpRad[n][month] += long_wave_rad_out * conv_mmd_to_Wm2; // HMS radiation output adjustment
            G_monthlyShortDownRad[n][month] += solar_rad * conv_mmd_to_Wm2;	// convert from mm/day to W/m2 // HMS radiation output adjustment
			G_monthlyNetShortRad[n][month] += net_short_wave_rad * conv_mmd_to_Wm2; //HMS net_short_wave_rad output
      		G_monthlyNetLongRad[n][month] += net_long_wave_rad * conv_mmd_to_Wm2; //HMS net_long_wave_rad output
			G_monthlyNetRadiation[n][month] += net_rad * conv_mmd_to_Wm2;	// convert from mm/day to W/m2
      
			// open water pet, if cell consists solely of open water bodies
			G_monthlyOpenWaterPET[n][month] += dailyOpenWaterPET;
			// pet over land area (whole cell is considered to be land)
			G_monthlyPET[n][month] += dailyPET;

            // FP 2015-06 // Only continental cells
            //if G_fwaterfreq + G_landfreq != 0.
            // FP20161010N001 Replacement: geo.G_contfreq[n] instead (geo.G_landfreq[n] + geo.G_fwaterfreq[n])
            // FP20161018N002 Reservoir operation start years: Replacement: routing.G_landfreq[n] instead geo.G_landfreq[n]; routing.G_fwaterfreq[n] instead geo.G_fwaterfreq[n]
             /*           if (n == 32879) {
                            cout << "G_landfreq " << routing.G_landfreq[n] << endl;
                            cout << "G_fwaterfreq " << routing.G_fwaterfreq[n] << endl;
                            cout << "G_contfreq " << geo.G_contfreq[n] << endl;
                            cout << "G_glolak" << routing.G_glo_lake[n] << endl;
                            cout << "G_loclak" << routing.G_loc_lake[n] << endl;
                            cout << "G_gw" << routing.G_glo_wetland[n] << endl;
                            cout << "G_lw" << routing.G_loc_wetland[n] << endl;
                            cout << "G_glores" << routing.G_glo_res[n] << endl;
                            cout << "G_lake" << routing.G_lake_area[n] << endl;
                            cout << "G_res" << routing.G_reservoir_area_full[n] << endl;
                        }*/
            if (geo.G_contcell[n]) {
                // consider pet and open water pet for land area and water fraction respectively
                G_monthlyTotalPET[n][month] += ( (routing.G_landfreq[n] * dailyPET + routing.G_fwaterfreq[n] * dailyOpenWaterPET)
                                                                               / geo.G_contfreq[n] );
                // grid cell aet; without canopy and snow evap
                G_monthlyAET[n][month] += ( dailyAET * landAreaFrac / geo.G_contfreq[n] );
                // for land AET in [mm]; accounts for land area fraction only (same as for pet and aet)
                G_monthlyLandAET[n][month] +=  ( dailyCellLandAET[n] * landAreaFrac / geo.G_contfreq[n] );
                G_monthlyLandAET_uncorr[n][month] +=  ( dailyCellLandAET_uncorr[n] * landAreaFrac / geo.G_contfreq[n] );
						
                G_monthlyLAI[n][month] += dailyLai;
                G_monthlyAlbedo[n][month] += (albedo * (landAreaFrac / 100.) + openWaterAlbedo * (1. - (landAreaFrac / 100.)));
                G_monthlyInterception[n][month]	+= (dailyCanopyEvapo * landAreaFrac	/ geo.G_contfreq[n]);
                // store canopy water content at end of month (consistent with water balance calculations)
                // *SE 2012-06* changed to monthly averages
                G_monthlycanopyWaterContent[n][month] += ( G_canopyWaterContent[n] * landAreaFrac / geo.G_contfreq[n] );
                //G_monthlycanopyWaterContent[n][month] = ( G_canopyWaterContent[n] * landAreaFrac / geo.G_contfreq[n] );
                G_monthlymaxcanopyWaterContent[n][month] += ( max_canopy_storage * landAreaFrac / geo.G_contfreq[n] );

                // snow parameters
                G_monthlySnowFall[n][month]	+= (dailySnowFall * landAreaFrac / geo.G_contfreq[n]);
                G_monthlyRainFall[n][month] += (dailyRainFall * landAreaFrac / geo.G_contfreq[n]); // FP20161011N001
                G_monthlySnowMelt[n][month]	+= (snowmelt * landAreaFrac / geo.G_contfreq[n]);
                G_monthlySnowEvap[n][month]	+= (dailySnowEvapo * landAreaFrac / geo.G_contfreq[n]);
                // store snow water equivalent at end of month (consistent with water balance calculations)
                // *SE 2012-06* changed to monthly averages
                G_monthlySWE[n][month] += ( G_snow[n] * landAreaFrac / geo.G_contfreq[n] );
                G_monthlySnowCoverAvg[n][month] += (dailySnowCoverFrac * landAreaFrac / geo.G_contfreq[n]);

                // soil water content
                // store soil water content at end of month (consistent with water balance calculations)
                // *SE 2012-06* changed to monthly averages
                G_monthlySoilWaterAvg[n][month] += ( G_soilWaterContent[n] * landAreaFrac / geo.G_contfreq[n] );
                //G_monthlySoilWaterAvg[n][month] = (G_soilWaterContent[n] * landAreaFrac / geo.G_contfreq[n]);
                // potential runoff from landcells, including the amount of gw-recharge
                G_monthlyRunoff[n][month] += (total_daily_runoff * landAreaFrac / geo.G_contfreq[n]);
                // urban runoff from landcells (accounts only for built-up area)
                G_monthlyUrbanRunoff[n][month] += immediate_runoff;
                // surface runoff from landcells
                G_monthlySurfaceRunoff[n][month] += (surface_runoff * landAreaFrac / geo.G_contfreq[n]);
                // goundwater recharge and gw-runoff into river network
                G_monthlyGwRecharge[n][month] += (daily_gw_recharge * landAreaFrac / geo.G_contfreq[n]);
                //G_monthlyGwRunoff[n][month]	+= (groundwater_runoff * landAreaFrac / ((float) geo.G_landfreq[n] + (float) geo.G_fwaterfreq[n]));

                // calculate monthly storage per grid cell
                // here: from vertical water balance (first part)
                // in km³
                if (options.grid_store_TypeForStorages == 1) {
                    G_monthlyCanopyStorage[n][month] += G_canopyWaterContent[n] * geo.areaOfCellByArrayPos(n) / 1000000.0 * landAreaFrac / 100.0; // G_canopyWaterContent in mm, G_monthlyCanopyStorage in km³
                    G_monthlySnowStorage[n][month] += G_snow[n]	* geo.areaOfCellByArrayPos(n) / 1000000.0 * landAreaFrac / 100.0; // G_snow in mm, G_monthlySnowStorage in km³
                    G_monthlySoilStorage[n][month] += G_soilWaterContent[n] * geo.areaOfCellByArrayPos(n) / 1000000.0 * landAreaFrac / 100.0; // G_soilWaterContent in mm, G_monthlySoilStorage in km³
                } else if (options.grid_store_TypeForStorages == 0) {
                    G_monthlyCanopyStorage[n][month] = G_canopyWaterContent[n] * geo.areaOfCellByArrayPos(n) / 1000000.0 * landAreaFrac / 100.0; // G_canopyWaterContent in mm, G_monthlyCanopyStorage in km³
                    G_monthlySnowStorage[n][month] = G_snow[n]	* geo.areaOfCellByArrayPos(n) / 1000000.0 * landAreaFrac / 100.0; // G_snow in mm, G_monthlySnowStorage in km³
                    G_monthlySoilStorage[n][month] = G_soilWaterContent[n] * geo.areaOfCellByArrayPos(n) / 1000000.0 * landAreaFrac / 100.0; // G_soilWaterContent in mm, G_monthlySoilStorage in km³
                }

            }
            // 100% water cells (e.g. of Caspian Sea) have now G_fwaterfreq + G_landfreq = 0. The water balance is calculated in the outlet cell of large lakes; therefore a calculation of e.g. precipitation for this cells leads to inconsistent quantification of global fluxes. A division through 0 has to be avoided now (alternative: define it as not land cells in ng).
            else {
                G_monthlyTotalPET[n][month] = 0.;
                G_monthlyAET[n][month] = 0.;
                G_monthlyLandAET[n][month] = 0.;
                G_monthlyLandAET_uncorr[n][month] = 0.;
                G_monthlyInterception[n][month] = 0.;
                G_monthlycanopyWaterContent[n][month] = 0.;
                G_monthlymaxcanopyWaterContent[n][month] = 0.;
                G_monthlySnowFall[n][month] = 0.;
                G_monthlyRainFall[n][month] = 0.; // FP20161011N001
                G_monthlySnowMelt[n][month] = 0.;
                G_monthlySnowEvap[n][month] = 0.;
                G_monthlySWE[n][month] = 0.;
                G_monthlySnowCoverAvg[n][month] = 0.;
                G_monthlySoilWaterAvg[n][month] = 0.;
                G_monthlyRunoff[n][month] = 0.;
                G_monthlySurfaceRunoff[n][month] = 0.;
                G_monthlyGwRecharge[n][month] = 0.;
                //G_monthlyGwRunoff[n][month] = 0.;
            }
            // end else

        }
        // endif options (monthly/yearly results)


        // daily 31 storage calculation
        if ((3 == options.grid_store) || (4 == options.grid_store)) {
            if ((options.outSingleStoragesDaily || options.outTotalWaterInStoragesDaily_km3 || options.outTotalWaterInStoragesDaily_mm)) G_daily31CanopyStorage[n][day_in_month - 1] = G_canopyWaterContent[n] * geo.areaOfCellByArrayPos(n) / 1000000.0 * landAreaFrac / 100.0; // G_canopyWaterContent in mm, G_dailyCanopyStorage in km³
            if ((options.outSingleStoragesDaily || options.outTotalWaterInStoragesDaily_km3 || options.outTotalWaterInStoragesDaily_mm)) G_daily31SnowStorage[n][day_in_month - 1] = G_snow[n]	* geo.areaOfCellByArrayPos(n) / 1000000.0 * landAreaFrac / 100.0; // G_snow in mm, G_dailySnowStorage in km³
            if ((options.outSingleStoragesDaily || options.outTotalWaterInStoragesDaily_km3 || options.outTotalWaterInStoragesDaily_mm)) G_daily31SoilStorage[n][day_in_month - 1] = G_soilWaterContent[n] * geo.areaOfCellByArrayPos(n) / 1000000.0 * landAreaFrac / 100.0; // G_soilWaterContent in mm, G_dailySoilStorage in km³
        }
        // endif options		
		
        if ((5 == options.grid_store) || (6 == options.grid_store)) {
            if ((options.outSingleStoragesDaily || options.outTotalWaterInStoragesDaily_km3 || options.outTotalWaterInStoragesDaily_mm)||(((options.scoutCanopyStor)||(options.scoutTWSkm3)||(options.scoutTWSmm))&&(2==options.day_store))) G_daily365CanopyStorage[n][day - 1] = G_canopyWaterContent[n] * geo.areaOfCellByArrayPos(n) / 1000000.0 * landAreaFrac / 100.0; // G_canopyWaterContent in mm, G_dailyCanopyStorage in km³
            if ((options.outSingleStoragesDaily || options.outTotalWaterInStoragesDaily_km3 || options.outTotalWaterInStoragesDaily_mm)||(((options.scoutSnowStor)||(options.scoutTWSkm3)||(options.scoutTWSmm))&&(2==options.day_store))) G_daily365SnowStorage[n][day - 1] = G_snow[n]	* geo.areaOfCellByArrayPos(n) / 1000000.0 * landAreaFrac / 100.0; // G_snow in mm, G_dailySnowStorage in km³
            if ((options.outSingleStoragesDaily || options.outTotalWaterInStoragesDaily_km3 || options.outTotalWaterInStoragesDaily_mm)||(((options.scoutSoilStor)||(options.scoutTWSkm3)||(options.scoutTWSmm))&&(2==options.day_store))) G_daily365SoilStorage[n][day - 1] = G_soilWaterContent[n] * geo.areaOfCellByArrayPos(n) / 1000000.0 * landAreaFrac / 100.0; // G_soilWaterContent in mm, G_dailySoilStorage in km³
        }
        // endif options


        // daily results
        // daily output option 31 start
		if ((3 == options.grid_store) || (4 == options.grid_store)) {

            // FP 2015-06 // Only continental cells
            // not 100% lake area
            // FP20161010N001 Replacement: geo.G_contfreq[n] instead (geo.G_landfreq[n] + geo.G_fwaterfreq[n])
            if (geo.G_contcell[n]) {

             // precipitation
                if (options.outPrecDaily) G_daily31Precip[n][day_in_month - 1] = dailyPrec;

             // PET and interception
                if (options.outPETDaily)			G_daily31PET[n][day_in_month - 1] = dailyPET;
                // FP20161018N002 Reservoir operation start years: Replacement: routing.G_landfreq[n] instead geo.G_landfreq[n]; routing.G_fwaterfreq[n] instead geo.G_fwaterfreq[n]
                if (options.outTotalPETDaily)		G_daily31TotalPET[n][day_in_month - 1] =
                                                                                        (routing.G_landfreq[n] * dailyPET + routing.G_fwaterfreq[n] * dailyOpenWaterPET)
                                                                                        / geo.G_contfreq[n];
                if (options.outLAIDaily)			G_daily31LAI[n][day_in_month - 1] = dailyLai;
                if (options.outLAIDaily)			G_daily31Kc[n][day_in_month - 1] = dailyKc;
                if (options.outAlbedoDaily)			G_daily31Albedo[n][day_in_month - 1] = albedo * (landAreaFrac / 100.) + openWaterAlbedo * (1. - (landAreaFrac / 100.));
                if (options.outInterceptionDaily)	G_daily31Interception[n][day_in_month - 1] = dailyCanopyEvapo * landAreaFrac / geo.G_contfreq[n];
                if (options.outCanopyWaterDaily)	G_daily31canopyWaterContent[n][day_in_month - 1] = G_canopyWaterContent[n] * landAreaFrac / geo.G_contfreq[n];
                if (options.outmaxCanopyWaterDaily)	G_daily31maxcanopyWaterContent[n][day_in_month - 1] = max_canopy_storage * landAreaFrac / geo.G_contfreq[n];
            // snow parameters
                if (options.outSnowFallDaily)		G_daily31SnowFall[n][day_in_month - 1] = dailySnowFall * landAreaFrac / geo.G_contfreq[n];
                if (options.outSnowMeltDaily)		G_daily31SnowMelt[n][day_in_month - 1] = snowmelt * landAreaFrac / geo.G_contfreq[n];
                if (options.outSWEDaily)			G_daily31SWE[n][day_in_month - 1] = G_snow[n] * landAreaFrac / geo.G_contfreq[n];
                if (options.outSnowCoverDaily)		G_daily31SnowCoverAvg[n][day_in_month - 1] = dailySnowCoverFrac * landAreaFrac / geo.G_contfreq[n];
            // soil water content
                if (options.outSoilWaterDaily)		G_daily31SoilWaterAvg[n][day_in_month - 1] = G_soilWaterContent[n] * landAreaFrac / geo.G_contfreq[n];
            // surface runoff from landcells
                if (options.outSurfaceRunoffDaily)	G_daily31SurfaceRunoff[n][day_in_month - 1] = surface_runoff * landAreaFrac / geo.G_contfreq[n];
            // goundwater recharge and gw-runoff into river network
                if (options.outGWRechargeDaily)		G_daily31GwRecharge[n][day_in_month - 1] = daily_gw_recharge * landAreaFrac / geo.G_contfreq[n];
                //if (options.outGWRunoffDaily)		G_daily31GwRunoff[n][day_in_month - 1] = groundwater_runoff * landAreaFrac / geo.G_contfreq[n];
               // if (options.outGWStorageDaily || options.outSingleStoragesDaily)		G_daily31GwStorage[n][day_in_month - 1] = G_groundwater[n] * landAreaFrac / geo.G_contfreq[n];
                // HMS radiation output adjustment start - daily radiation output
                if(options.outExtRadiationDaily)		G_daily31ExtRad[n][day_in_month - 1] = ext_rad * conv_mmd_to_Wm2;
                if(options.outShortDownRadiationDaily)		G_daily31ShortDownRad[n][day_in_month - 1] = solar_rad * conv_mmd_to_Wm2;
                if(options.outShortUpRadiationDaily)		G_daily31ShortUpRad[n][day_in_month - 1] = refl_rad * conv_mmd_to_Wm2;
                if(options.outNetShortWaveRadiationDaily)		G_daily31NetShortRad[n][day_in_month - 1] = net_short_wave_rad * conv_mmd_to_Wm2;
                if(options.outLongDownRadiationDaily)		G_daily31LongDownRad[n][day_in_month - 1] = long_wave_rad_in * conv_mmd_to_Wm2;
                if(options.outLongUpRadiationDaily)		G_daily31LongUpRad[n][day_in_month - 1] = long_wave_rad_out * conv_mmd_to_Wm2;
                if(options.outNetLongWaveRadiationDaily)		G_daily31NetLongRad[n][day_in_month -  1] = net_long_wave_rad * conv_mmd_to_Wm2;
                if(options.outNetRadiationDaily)		G_daily31NetRadiation[n][day_in_month - 1] = net_rad * conv_mmd_to_Wm2;
                // HMS radiation output adjustment end

            }
            // 100% water cells (e.g. of Caspian Sea) have now G_fwaterfreq + G_landfreq = 0. The water balance is calculated in the outlet cell of large lakes; therefore a calculation of e.g. precipitation for this cells leads to inconsistent quantification of global fluxes. A division through 0 has to be avoided now (alternative: define it as not land cells in ng).
            else {
                if (options.outTotalPETDaily)		G_daily31TotalPET[n][day_in_month - 1] = 0.;
                if (options.outInterceptionDaily)	G_daily31Interception[n][day_in_month - 1] = 0.;
                if (options.outCanopyWaterDaily)	G_daily31canopyWaterContent[n][day_in_month - 1] = 0.;
                if (options.outmaxCanopyWaterDaily)	G_daily31maxcanopyWaterContent[n][day_in_month - 1] = 0.;
                if (options.outSnowFallDaily)		G_daily31SnowFall[n][day_in_month - 1] = 0.;
                if (options.outSnowMeltDaily)		G_daily31SnowMelt[n][day_in_month - 1] = 0.;
                if (options.outSWEDaily)		G_daily31SWE[n][day_in_month - 1] = 0.;
                if (options.outSoilWaterDaily)		G_daily31SoilWaterAvg[n][day_in_month - 1] = 0.;
                if (options.outSurfaceRunoffDaily)	G_daily31SurfaceRunoff[n][day_in_month - 1] = 0.;
                if (options.outGWRechargeDaily)		G_daily31GwRecharge[n][day_in_month - 1] = 0.;
                //if (options.outGWRunoffDaily)		G_daily31GwRunoff[n][day_in_month - 1] = 0.;
               // if (options.outGWStorageDaily || options.outSingleStoragesDaily)		G_daily31GwStorage[n][day_in_month - 1] = 0.;
                if (options.outPrecDaily) G_daily31Precip[n][day_in_month - 1] = 0.;
                if (options.outLAIDaily)			G_daily31LAI[n][day_in_month - 1] = 0.;
                if (options.outLAIDaily)			G_daily31Kc[n][day_in_month - 1] = 0.;
                if (options.outAlbedoDaily)			G_daily31Albedo[n][day_in_month - 1] = 0.;
            }

        }
        // endif options
        // daily output option 31 end

    // daily output option 365 start
    if ((5 == options.grid_store) || (6 == options.grid_store)) {
        // 100% water cells (e.g. of Caspian Sea) have now G_fwaterfreq + G_landfreq = 0. The water balance is calculated in the outlet cell of large lakes; therefore a calculation of e.g. precipitation for this cells leads to inconsistent quantification of global fluxes. A division through 0 has to be avoided now (alternative: define it as not land cells in ng).
        // not 100% lake area
        // FP20161010N001 Replacement: geo.G_contfreq[n] instead (geo.G_landfreq[n] + geo.G_fwaterfreq[n])
        if (geo.G_contcell[n]) {
            // precipitation
            if ((options.outPrecDaily)||((options.scoutPrecip)&&(2==options.day_store)))		G_daily365Precip[n][day - 1] = dailyPrec;

        // PET and interception
            if ((options.outPETDaily)||((options.scoutLandPET)&&(2==options.day_store)))		G_daily365PET[n][day - 1] = dailyPET;
            // FP20161018N002 Reservoir operation start years: Replacement: routing.G_landfreq[n] instead geo.G_landfreq[n]; routing.G_fwaterfreq[n] instead geo.G_fwaterfreq[n]
            if ((options.outTotalPETDaily)||((options.scoutCellPET)&&(2==options.day_store)))	G_daily365TotalPET[n][day - 1] = (routing.G_landfreq[n] * dailyPET + routing.G_fwaterfreq[n] * dailyOpenWaterPET) / geo.G_contfreq[n];
            if ((options.outLAIDaily)||((options.scoutLAI)&&(2==options.day_store)))		G_daily365LAI[n][day - 1] = dailyLai;
            if ((options.outLAIDaily)||((options.scoutKc)&&(2==options.day_store)))		G_daily365Kc[n][day - 1] = dailyKc;
            if ((options.outAlbedoDaily)||((options.scoutAlbedo)&&(2==options.day_store)))		G_daily365Albedo[n][day - 1] = albedo * (landAreaFrac / 100.) + openWaterAlbedo * (1. - (landAreaFrac / 100.));
            if ((options.outInterceptionDaily)||((options.scoutInterception)&&(2==options.day_store)))	G_daily365Interception[n][day - 1] = dailyCanopyEvapo * landAreaFrac / geo.G_contfreq[n];
            if ((options.outCanopyWaterDaily)||((options.scoutCanopyWater)&&(2==options.day_store)))	G_daily365canopyWaterContent[n][day - 1] = G_canopyWaterContent[n] * landAreaFrac / geo.G_contfreq[n]; // [mm]
            if ((options.outmaxCanopyWaterDaily)||((options.scoutmaxCanopyWater)&&(2==options.day_store)))	G_daily365maxcanopyWaterContent[n][day - 1]	= max_canopy_storage * landAreaFrac / geo.G_contfreq[n];
        // snow parameters
            if ((options.outSnowFallDaily)||((options.scoutSnowfall)&&(2==options.day_store)))	G_daily365SnowFall[n][day - 1] = dailySnowFall * landAreaFrac / geo.G_contfreq[n];
            if ((options.outSnowMeltDaily)||((options.scoutSnowmelt)&&(2==options.day_store)))	G_daily365SnowMelt[n][day - 1] = snowmelt * landAreaFrac / geo.G_contfreq[n];
            if ((options.outSWEDaily)||((options.scoutSnowWater)&&(2==options.day_store)))		G_daily365SWE[n][day - 1] = G_snow[n] * landAreaFrac / geo.G_contfreq[n]; // [mm]
            if ((options.outSnowCoverDaily)||((options.scoutSnowCov)&&(2==options.day_store)))	G_daily365SnowCoverAvg[n][day - 1] = dailySnowCoverFrac * landAreaFrac / geo.G_contfreq[n];
        // soil water content
            if ((options.outSoilWaterDaily)||((options.scoutSoilWater)&&(2==options.day_store)))	G_daily365SoilWaterAvg[n][day - 1] = G_soilWaterContent[n] * landAreaFrac / geo.G_contfreq[n]; // [mm]
        // surface runoff from landcells
            if ((options.outSurfaceRunoffDaily)||((options.scoutSurfaceRunoff)&&(2==options.day_store)))	G_daily365SurfaceRunoff[n][day - 1] = surface_runoff * landAreaFrac / geo.G_contfreq[n];
        // groundwater recharge and gw-runoff into river network
            if ((options.outGWRechargeDaily)||((options.scoutGwRecharge)&&(2==options.day_store)))	G_daily365GwRecharge[n][day - 1] = daily_gw_recharge * landAreaFrac / geo.G_contfreq[n];
            //   if ((options.outGWRunoffDaily)||(options.scoutGwRunoff))	G_daily365GwRunoff[n][day - 1] = groundwater_runoff * landAreaFrac / geo.G_contfreq[n];
            // if ((options.outGWStorageDaily || options.outSingleStoragesDaily)||(2 == options.day_store))	G_daily365GwStorage[n][day - 1] = G_groundwater[n] * landAreaFrac / geo.G_contfreq[n];
        // HMS radiation output adjustment start - daily radiation output
            if((options.outExtRadiationDaily)||((options.scoutExtRad)&&(2==options.day_store)))		G_daily365ExtRad[n][day-1] = ext_rad * conv_mmd_to_Wm2;
            if((options.outShortDownRadiationDaily)||((options.scoutShortDownRad)&&(2==options.day_store)))		G_daily365ShortDownRad[n][day-1] = solar_rad * conv_mmd_to_Wm2;
            if((options.outShortUpRadiationDaily)||((options.scoutShortUpRad)&&(2==options.day_store)))		G_daily365ShortUpRad[n][day-1] = refl_rad * conv_mmd_to_Wm2;
            if((options.outNetShortWaveRadiationDaily)||((options.scoutNetShortRad)&&(2==options.day_store)))		G_daily365NetShortRad[n][day-1] = net_short_wave_rad * conv_mmd_to_Wm2;
            if((options.outLongDownRadiationDaily)||((options.scoutLongDownRad)&&(2==options.day_store)))		G_daily365LongDownRad[n][day-1] = long_wave_rad_in * conv_mmd_to_Wm2;
            if((options.outLongUpRadiationDaily)||((options.scoutLongUpRad)&&(2==options.day_store)))		G_daily365LongUpRad[n][day-1] = long_wave_rad_out * conv_mmd_to_Wm2;
            if((options.outNetLongWaveRadiationDaily)||((options.scoutNetLongRad)&&(2==options.day_store)))		G_daily365NetLongRad[n][day-1] = net_long_wave_rad * conv_mmd_to_Wm2;
            if((options.outNetRadiationDaily)||((options.scoutNetRad)&&(2==options.day_store)))		G_daily365NetRadiation[n][day-1] = net_rad * conv_mmd_to_Wm2;
            if(options.scoutTemp&&(2==options.day_store)) 			G_daily365TempC[n][day-1] = dailyTempC;
        // HMS radiation output adjustment end
        }
        else {
            if (options.outTotalPETDaily || (options.scoutCellPET&&(2==options.day_store)))	G_daily365TotalPET[n][day - 1] = 0.;
            if (options.outLAIDaily || (options.scoutLAI&&(2==options.day_store)))		G_daily365LAI[n][day - 1] = 0.;
            if (options.outLAIDaily || (options.scoutKc&&(2==options.day_store)))		G_daily365Kc[n][day - 1] = 0.;
            if (options.outAlbedoDaily || (options.scoutAlbedo&&(2==options.day_store)))		G_daily365Albedo[n][day - 1] = 0.;
            if (options.outInterceptionDaily || (options.scoutInterception&&(2==options.day_store)))	G_daily365Interception[n][day - 1] = 0.;
            if (options.outCanopyWaterDaily || (options.scoutCanopyWater&&(2==options.day_store)))	G_daily365canopyWaterContent[n][day - 1] = 0.;
            if (options.outmaxCanopyWaterDaily || (options.scoutmaxCanopyWater&&(2==options.day_store)))	G_daily365maxcanopyWaterContent[n][day - 1] = 0.;
            if (options.outSnowFallDaily || (options.scoutSnowfall&&(2==options.day_store)))	G_daily365SnowFall[n][day - 1] = 0.;
            if (options.outSnowMeltDaily || (options.scoutSnowmelt&&(2==options.day_store)))	G_daily365SnowMelt[n][day - 1] = 0.;
            if (options.outSWEDaily || (options.scoutSnowWater&&(2==options.day_store)))		G_daily365SWE[n][day - 1] = 0.;
            if (options.outSnowCoverDaily || (options.scoutSnowCov&&(2==options.day_store)))	G_daily365SnowCoverAvg[n][day - 1] = 0.;
            if (options.outSoilWaterDaily || (options.scoutSoilWater&&(2==options.day_store)))	G_daily365SoilWaterAvg[n][day - 1] = 0.;
            if (options.outSurfaceRunoffDaily || (options.scoutSurfaceRunoff&&(2==options.day_store)))	G_daily365SurfaceRunoff[n][day - 1] = 0.;
            if (options.outGWRechargeDaily || (options.scoutGwRecharge&&(2==options.day_store)))	G_daily365GwRecharge[n][day - 1] = 0.;
            //if (options.outGWRunoffDaily || options.scoutGwRunoff)	G_daily365GwRunoff[n][day - 1] = 0.;
         //   if (options.outGWStorageDaily || options.outSingleStoragesDaily)	G_daily365GwStorage[n][day - 1] = 0.;
            if (options.outPETDaily || (options.scoutLandPET&&(2==options.day_store)))		G_daily365PET[n][day - 1] = 0.;
            if (options.outPrecDaily || (options.scoutPrecip&&(2==options.day_store)))		G_daily365Precip[n][day - 1] = 0.;
        }
            // endelse

    }
    // daily output option 365 end

		// Superbasin calculation currently obsolete, probably not working correctly in WaterGAP2.2b, do in postprocessing  // FP 2015-10
		// store values for superbasins only if they are needed
//		if (options.outDailyValues || options.outDailyInterception || options.outSuperbasinClimate) {
//			if (G_sbasin[n] > 0) {
//				signed short superbasin;
//				const double landArea = landAreaFrac/100. * geo.areaOfCellByArrayPos(n);
	
//				superbasin = G_sbasin[n] - 1;

//				precAsSnowSb[superbasin] += dailySnowFall * landArea;
//				dailyTempSb[superbasin][day - 1] += dailyTempC * landArea;
//				dailyPrecSb[superbasin][day - 1] += dailyPrec * landArea;
//				dailyEffPrecSb[superbasin][day - 1] += dailyEffPrec * landArea;
//				dailySnowSb[superbasin][day - 1] += G_snow[n] * landArea;
//				dailyAET_Sb[superbasin][day - 1] += dailyAET * landArea;
//				dailyPET_Sb[superbasin][day - 1] += dailyPET * landArea;
//				dailyRunoffSb[superbasin][day - 1] += total_daily_runoff * landArea;
//				dailySaturationSb[superbasin][day - 1] += (G_soilWaterContent[n] / maxSoilWaterCap.G_Smax[n]) * landArea;
//				dailySunshineSb[superbasin][day - 1] += dailySunshine * landArea;
//				dailyLaiSb[superbasin][day - 1] += dailyLai * landArea;
//				dailyCanopyWaterSb[superbasin][day - 1] += canopy_water_content * landArea;
//				dailyCanopyEvapoSb[superbasin][day - 1] += dailyCanopyEvapo * landArea;
//			}
//		}

	}
	// end of loop for all cells (1 == G_toBeCalculated[n])

}
// end calcNewDay


void dailyWaterBalanceClass::init(const char *inputDir, const short nBasins) {
	extern gridioClass gridIO;
	char filename[250];

	// Superbasin calculation currently obsolete, probably not working correctly in WaterGAP2.2b, do in postprocessing  // FP 2015-10
	// number of superbasins varies.
	// therefore size of the arrays is allocated dynamically.
//	nSpecBasins = nBasins; // number of specified basins

//	if (options.outDailyValues || options.outDailyInterception || options.outSuperbasinClimate) {
//                snowDaysSb = new double[nSpecBasins];
//                precAsSnowSb = new double[nSpecBasins];

//                dailyTempSb = new double[nSpecBasins][365];
//                dailyPrecSb = new double[nSpecBasins][365];
//                dailyEffPrecSb = new double[nSpecBasins][365];
//                dailySnowSb = new double[nSpecBasins][365];
//                dailyAET_Sb = new double[nSpecBasins][365];
//                dailyPET_Sb = new double[nSpecBasins][365];
//                dailyRunoffSb = new double[nSpecBasins][365];
//                dailySaturationSb = new double[nSpecBasins][365];
//                dailySunshineSb = new double[nSpecBasins][365];
//                dailyLaiSb = new double[nSpecBasins][365];
//                dailyCanopyWaterSb = new double[nSpecBasins][365];
//                dailyCanopyEvapoSb = new double[nSpecBasins][365];
//        }

	// store monthly/yearly grids
	if ((2 == options.grid_store) || (4 == options.grid_store) || (6 == options.grid_store)) {
                G_monthlyPrecipitation = new double[ng][12];
                G_monthlyTemperature = new double[ng][12];
        G_monthlySunshine = new double[ng][12]; // sunshine percentage output
                G_monthlyNetRadiation = new double[ng][12];
                G_monthlyNetShortRad = new double[ng][12]; //HMS net_short_wave_rad output
                G_monthlyNetLongRad = new double[ng][12]; //HMS net_long_wave_rad output
                G_monthlyShortDownRad = new double[ng][12]; // HMS radiation output adjustment
                G_monthlyExtRad = new double[ng][12]; // HMS radiation output adjustment
                G_monthlyShortUpRad = new double[ng][12]; // HMS radiation output adjustment
                G_monthlyLongDownRad = new double[ng][12]; // HMS radiation output adjustment
                G_monthlyLongUpRad = new double[ng][12]; // HMS radiation output adjustment
                G_monthlyLAI = new double[ng][12];
                G_monthlyAlbedo = new double[ng][12];
                G_monthlyInterception = new double[ng][12];
                G_monthlycanopyWaterContent = new double[ng][12];
                G_monthlymaxcanopyWaterContent = new double[ng][12];
                G_monthlyAET = new double[ng][12];
                G_monthlyLandAET = new double[ng][12]; // land evap only
                G_monthlyLandAET_uncorr = new double[ng][12]; // land evap only
                G_monthlyPET = new double[ng][12];
                G_monthlyOpenWaterPET = new double[ng][12];
                G_monthlyTotalPET = new double[ng][12];
                G_monthlyRunoff = new double[ng][12];
                G_monthlyUrbanRunoff = new double[ng][12];
                G_monthlySurfaceRunoff = new double[ng][12];
	//G_monthlyGwRunoff = new float[ng][12];
                G_monthlyGwRecharge = new double[ng][12];
                G_monthlySnowCoverAvg = new double[ng][12];
                G_monthlySWE = new double[ng][12];
                G_monthlySnowFall = new double[ng][12];
                G_monthlyRainFall = new double[ng][12]; // FP20161011N001
                G_monthlySnowMelt = new double[ng][12];
                G_monthlySnowEvap = new double[ng][12];
                G_monthlySoilWaterAvg = new double[ng][12];
		// G_SnowInElevation = new float [ng_land][100]; //SNOW
		//G_Elevation = new short [ng_land][103]; //SNOW
	}

        //monthly storage components
        //G_monthlyGwStorage now calculated in routing.cpp
	//G_monthlyGwStorage = new double[ng][12];
	G_monthlyCanopyStorage = new double[ng][12];
	G_monthlySnowStorage = new double[ng][12];
	G_monthlySoilStorage = new double[ng][12];
	G_monthly_totalWaterInStorages = new double[ng][12];
	G_monthly_totalWaterInStorages_mm = new double[ng][12];

	// store daily grids

        // daily output option 31 start
	if ((3 == options.grid_store) || (4 == options.grid_store)) {

                if (options.outSurfaceRunoffDaily)	G_daily31SurfaceRunoff		= new double[ng][31];
                //if (options.outGWRunoffDaily)		G_daily31GwRunoff			= new float[ng][31];
                if (options.outGWRechargeDaily)		G_daily31GwRecharge			= new double[ng][31];
       // if (options.outGWStorageDaily || options.outSingleStoragesDaily)		G_daily31GwStorage			= new float[ng][31];
                if (options.outSoilWaterDaily)		G_daily31SoilWaterAvg		= new double[ng][31];
                if (options.outLAIDaily)			G_daily31LAI				= new double[ng][31];
                if (options.outLAIDaily)			G_daily31Kc					= new double[ng][31];
                if (options.outAlbedoDaily)			G_daily31Albedo				= new double[ng][31];
                if (options.outInterceptionDaily)	G_daily31Interception		= new double[ng][31];
                if (options.outCanopyWaterDaily)	G_daily31canopyWaterContent	= new double[ng][31];
                if (options.outmaxCanopyWaterDaily)	G_daily31maxcanopyWaterContent = new double[ng][31];
                if (options.outPrecDaily)			G_daily31Precip				= new double[ng][31];

                if (options.outPETDaily)			G_daily31PET				= new double[ng][31];
                if (options.outTotalPETDaily)		G_daily31TotalPET			= new double[ng][31];
                if (options.outSnowCoverDaily)		G_daily31SnowCoverAvg		= new double[ng][31];
                if (options.outSWEDaily)			G_daily31SWE				= new double[ng][31];
                if (options.outSnowFallDaily)		G_daily31SnowFall			= new double[ng][31];
                if (options.outSnowMeltDaily)		G_daily31SnowMelt			= new double[ng][31];
        if (options.outSingleStoragesDaily || options.outTotalWaterInStoragesDaily_km3 || options.outTotalWaterInStoragesDaily_mm)                G_daily31CanopyStorage			= new double[ng][31];
        if (options.outSingleStoragesDaily || options.outTotalWaterInStoragesDaily_km3 || options.outTotalWaterInStoragesDaily_mm)                G_daily31SnowStorage			= new double[ng][31];
        if (options.outSingleStoragesDaily || options.outTotalWaterInStoragesDaily_km3 || options.outTotalWaterInStoragesDaily_mm)                G_daily31SoilStorage			= new double[ng][31];

	    // HMS radiation output adjustment start - daily radiation output
        if(options.outExtRadiationDaily)		G_daily31ExtRad				= new double[ng][31];
        if(options.outShortDownRadiationDaily)		G_daily31ShortDownRad				= new double[ng][31];
        if(options.outShortUpRadiationDaily)		G_daily31ShortUpRad				= new double[ng][31];
        if(options.outNetShortWaveRadiationDaily)		G_daily31NetShortRad				= new double[ng][31];
        if(options.outLongDownRadiationDaily)		G_daily31LongDownRad				= new double[ng][31];
        if(options.outLongUpRadiationDaily)		G_daily31LongUpRad				= new double[ng][31];
        if(options.outNetLongWaveRadiationDaily)		G_daily31NetLongRad				= new double[ng][31];
        if(options.outNetRadiationDaily)		G_daily31NetRadiation				= new double[ng][31];
    
    	// HMS radiation output adjustment end - daily radiation output
        }
        // daily output option 31 end

        // daily output option 365 start
	if ((5 == options.grid_store) || (6 == options.grid_store)) {
                if ((options.outSurfaceRunoffDaily)||((options.scoutSurfaceRunoff)&&(2==options.day_store)))	G_daily365SurfaceRunoff = new double[ng][365];

                if ((options.outGWRechargeDaily)||((options.scoutGwRecharge)&&(2==options.day_store)))		G_daily365GwRecharge = new double[ng][365];
       // if ((options.outGWStorageDaily || options.outSingleStoragesDaily)||(2 == options.day_store))		G_daily365GwStorage = new float[ng][365];
                if ((options.outSoilWaterDaily)||((options.scoutSoilWater)&&(2==options.day_store)))		G_daily365SoilWaterAvg = new double[ng][365];
                if ((options.outLAIDaily)||((options.scoutLAI)&&(2==options.day_store)))			G_daily365LAI = new double[ng][365];
                if ((options.outLAIDaily)||((options.scoutKc)&&(2==options.day_store)))			G_daily365Kc = new double[ng][365];
                if ((options.outAlbedoDaily)||((options.scoutAlbedo)&&(2==options.day_store)))			G_daily365Albedo = new double[ng][365];
                if ((options.outInterceptionDaily)||((options.scoutInterception)&&(2==options.day_store)))	G_daily365Interception = new double[ng][365];
                if ((options.outCanopyWaterDaily)||((options.scoutCanopyWater)&&(2==options.day_store)))	G_daily365canopyWaterContent = new double[ng][365];
                if ((options.outmaxCanopyWaterDaily)||((options.scoutmaxCanopyWater)&&(2==options.day_store)))	G_daily365maxcanopyWaterContent = new double[ng][365];
        if ((options.outPrecDaily)||((options.scoutPrecip)&&(2==options.day_store)))			G_daily365Precip = new double[ng][365];
                if ((options.outPETDaily)||((options.scoutLandPET)&&(2==options.day_store)))			G_daily365PET = new double[ng][365];
                if ((options.outTotalPETDaily)||((options.scoutCellPET)&&(2==options.day_store)))		G_daily365TotalPET = new double[ng][365];
                if ((options.outSnowCoverDaily)||((options.scoutSnowCov)&&(2==options.day_store)))		G_daily365SnowCoverAvg = new double[ng][365];
                if ((options.outSWEDaily)||((options.scoutSnowWater)&&(2==options.day_store)))			G_daily365SWE = new double[ng][365];
                if ((options.outSnowFallDaily)||((options.scoutSnowfall)&&(2==options.day_store)))		G_daily365SnowFall = new double[ng][365];
                if ((options.outSnowMeltDaily)||((options.scoutSnowmelt)&&(2==options.day_store)))		G_daily365SnowMelt = new double[ng][365];
        if ((options.outSingleStoragesDaily || options.outTotalWaterInStoragesDaily_km3 || options.outTotalWaterInStoragesDaily_mm)||(((options.scoutCanopyStor)||(options.scoutTWSkm3)||(options.scoutTWSmm))&&(2==options.day_store)))                G_daily365CanopyStorage			= new double[ng][365];
        if ((options.outSingleStoragesDaily || options.outTotalWaterInStoragesDaily_km3 || options.outTotalWaterInStoragesDaily_mm)||(((options.scoutSnowStor)||(options.scoutTWSkm3)||(options.scoutTWSmm))&&(2==options.day_store)))                G_daily365SnowStorage			= new double[ng][365];
        if ((options.outSingleStoragesDaily || options.outTotalWaterInStoragesDaily_km3 || options.outTotalWaterInStoragesDaily_mm)||(((options.scoutSoilStor)||(options.scoutTWSkm3)||(options.scoutTWSmm))&&(2==options.day_store)))                G_daily365SoilStorage			= new double[ng][365];
		// HMS radiation output adjustment start - daily radiation output
        if((options.outExtRadiationDaily)||((options.scoutExtRad)&&(2==options.day_store)))		G_daily365ExtRad				= new double[ng][365];
        if((options.outShortDownRadiationDaily)||((options.scoutShortDownRad)&&(2==options.day_store)))		G_daily365ShortDownRad				= new double[ng][365];
        if((options.outShortUpRadiationDaily)||((options.scoutShortUpRad)&&(2==options.day_store)))		G_daily365ShortUpRad				= new double[ng][365];
        if((options.outNetShortWaveRadiationDaily)||((options.scoutNetShortRad)&&(2==options.day_store)))		G_daily365NetShortRad				= new double[ng][365];
        if((options.outLongDownRadiationDaily)||((options.scoutLongDownRad)&&(2==options.day_store)))		G_daily365LongDownRad				= new double[ng][365];
        if((options.outLongUpRadiationDaily)||((options.scoutLongUpRad)&&(2==options.day_store)))		G_daily365LongUpRad				= new double[ng][365];
        if((options.outNetLongWaveRadiationDaily)||((options.scoutNetLongRad)&&(2==options.day_store)))		G_daily365NetLongRad				= new double[ng][365];
        if((options.outNetRadiationDaily)||((options.scoutNetRad)&&(2==options.day_store)))		G_daily365NetRadiation				= new double[ng][365];
    	// HMS radiation output adjustment end - daily radiation output
        if ((options.scoutTemp&&(2==options.day_store)))		G_daily365TempC				= new double[ng][365];
    	//HMS daily output option 365 end

	}
        // daily output option 365 end

	// read specific data for different lct
	if (!LCTdataIsRead)
		readLCTdata(inputDir);

	// set parameters (except cell-specific calibration parameters)
	setParameters();

	// these two variables have to be initialized:
	// in daily.cpp they are used for land cells only, but
	// loop in routing.cpp accounts for all cells;
	for (int n = 0; n < ng; n++) {
		G_dailyLocalGWRunoff[n] = 0.;
		G_dailyLocalSurfaceRunoff[n] = 0.;
		G_dailyStorageTransfer[n] = 0.;
	}

	// added for innerannual variation in albedo
	// albedo ----------------------------------
	if (options.calc_albedo) {
		G_AlbedoCycle = new float[ng][12];

		sprintf(filename, "%s/G_ALBEDO_CYCLE.12.UNF0", options.input_dir);
		gridIO.readUnfFile(filename, ng * 12, &G_AlbedoCycle[0][0]);
	}
	// end albedo ------------------------------

}

void dailyWaterBalanceClass::annualInit() {
//	normalizedFlag = false; // Superbasin calculation currently obsolete, probably not working correctly in WaterGAP2.2b, do in postprocessing  // FP 2015-10

	if ((2 == options.grid_store) || (4 == options.grid_store) || (6 == options.grid_store)) {
		for (int n = 0; n < ng; ++n) {
			for (short month = 0; month < 12; month++) {
				G_monthlyPrecipitation[n][month] = 0.;
				G_monthlyTemperature[n][month] = 0.;
				G_monthlySunshine[n][month] = 0.;   // sunshine percentage output
				G_monthlyNetRadiation[n][month] = 0.;
				G_monthlyNetShortRad[n][month] = 0.; // HMS net_short_wave_rad output
        		G_monthlyNetLongRad[n][month] = 0.; // HMS net_long_wave_rad output
                                G_monthlyShortDownRad[n][month] = 0.; // HMS radiation output adjustment
				G_monthlyExtRad[n][month] = 0.; // HMS radiation output adjustment
                                G_monthlyShortUpRad[n][month] = 0.; // HMS radiation output adjustment
                                G_monthlyLongDownRad[n][month] = 0.; // HMS radiation output adjustment
                                G_monthlyLongUpRad[n][month] = 0.; // HMS radiation output adjustment
				G_monthlyLAI[n][month] = 0.;
				G_monthlyAlbedo[n][month] = 0.;
				G_monthlyInterception[n][month] = 0.;
				G_monthlycanopyWaterContent[n][month] = 0.;
				G_monthlymaxcanopyWaterContent[n][month] = 0.;
				G_monthlyAET[n][month] = 0.;
				G_monthlyLandAET[n][month] = 0.;
				G_monthlyLandAET_uncorr[n][month] = 0.;
				G_monthlyPET[n][month] = 0.;
				G_monthlyOpenWaterPET[n][month] = 0.;
				G_monthlyTotalPET[n][month] = 0.;
				G_monthlyRunoff[n][month] = 0.;
				G_monthlyUrbanRunoff[n][month] = 0.;
				G_monthlySurfaceRunoff[n][month] = 0.;
			//	G_monthlyGwRunoff[n][month] = 0.;
				G_monthlyGwRecharge[n][month] = 0.;
				G_monthlySnowCoverAvg[n][month] = 0.;
				G_monthlySWE[n][month] = 0.;
				G_monthlySnowMelt[n][month] = 0.;
				G_monthlySnowEvap[n][month] = 0.;
				G_monthlySnowFall[n][month] = 0.;
                G_monthlyRainFall[n][month] = 0.; // FP20161011N001
				G_monthlySoilWaterAvg[n][month] = 0.;

				G_monthlyCanopyStorage[n][month] = 0.;
				G_monthlySnowStorage[n][month] = 0.;
				G_monthlySoilStorage[n][month] = 0.;
				G_monthly_totalWaterInStorages[n][month] = 0.;
				G_monthly_totalWaterInStorages_mm[n][month] = 0.;

			}
		}
	}

	// Superbasin calculation currently obsolete, probably not working correctly in WaterGAP2.2b, do in postprocessing  // FP 2015-10
//	if (options.outDailyValues || options.outDailyInterception	|| options.outSuperbasinClimate) {
//		for (short i = 0; i < nSpecBasins; ++i) {
//			snowDaysSb[i] = 0;
//			precAsSnowSb[i] = 0;
//			for (short day = 0; day < 365; day++) {
//				dailyTempSb[i][day] = 0.;
//				dailyPrecSb[i][day] = 0.;
//				dailyEffPrecSb[i][day] = 0.;
//				dailySnowSb[i][day] = 0.;
//				dailyAET_Sb[i][day] = 0.;
//				dailyPET_Sb[i][day] = 0.;
//				dailyRunoffSb[i][day] = 0.;
//				dailySaturationSb[i][day] = 0.;
//				dailySunshineSb[i][day] = 0.;
//				dailyLaiSb[i][day] = 0.;
//				dailyCanopyWaterSb[i][day] = 0.;
//				dailyCanopyEvapoSb[i][day] = 0.;
//			}
//		}
//	}
	// end if options.outSuperbasinClimate
}

void dailyWaterBalanceClass::daily31outInit() {
	short d;
	int n;

        // daily output option 31 start
    //if ((3 == options.grid_store) || (4 == options.grid_store)) {
		for (n = 0; n < ng; n++)
			for (d = 0; d < 31; d++) {
				if (options.outSurfaceRunoffDaily)	G_daily31SurfaceRunoff[n][d] = 0.;
                                //if (options.outGWRunoffDaily)		G_daily31GwRunoff[n][d] = 0.;
				if (options.outGWRechargeDaily)		G_daily31GwRecharge[n][d] = 0.;
           //     if (options.outGWStorageDaily || options.outSingleStoragesDaily)		G_daily31GwStorage[n][d] = 0.;
				if (options.outSoilWaterDaily)		G_daily31SoilWaterAvg[n][d] = 0.;
				if (options.outLAIDaily)			G_daily31LAI[n][d] = 0.;
				if (options.outLAIDaily)			G_daily31Kc[n][d] = 0.;
				if (options.outAlbedoDaily)			G_daily31Albedo[n][d] = 0.;
				if (options.outInterceptionDaily)	G_daily31Interception[n][d] = 0.;
				if (options.outCanopyWaterDaily)	G_daily31canopyWaterContent[n][d] = 0.;
				if (options.outmaxCanopyWaterDaily)	G_daily31maxcanopyWaterContent[n][d] = 0.;
                                if (options.outPrecDaily)			G_daily31Precip[n][d] = 0.;
				if (options.outPETDaily)			G_daily31PET[n][d] = 0.;
				if (options.outTotalPETDaily)		G_daily31TotalPET[n][d] = 0.;
				if (options.outSnowCoverDaily)		G_daily31SnowCoverAvg[n][d] = 0.;
				if (options.outSWEDaily)			G_daily31SWE[n][d] = 0.;
				if (options.outSnowFallDaily)		G_daily31SnowFall[n][d] = 0.;
				if (options.outSnowMeltDaily)		G_daily31SnowMelt[n][d] = 0.;
                if (options.outSingleStoragesDaily || options.outTotalWaterInStoragesDaily_km3 || options.outTotalWaterInStoragesDaily_mm)                G_daily31CanopyStorage[n][d] = 0.;
                if (options.outSingleStoragesDaily || options.outTotalWaterInStoragesDaily_km3 || options.outTotalWaterInStoragesDaily_mm)                G_daily31SnowStorage[n][d] = 0.;
                if (options.outSingleStoragesDaily || options.outTotalWaterInStoragesDaily_km3 || options.outTotalWaterInStoragesDaily_mm)                G_daily31SoilStorage[n][d] = 0.;
				// HMS radiation output adjustment start - daily radiation output
                                if(options.outExtRadiationDaily)		G_daily31ExtRad[n][d]		= 0.;
                                if(options.outShortDownRadiationDaily)		G_daily31ShortDownRad[n][d]		= 0.;
                                if(options.outShortUpRadiationDaily)		G_daily31ShortUpRad[n][d]		= 0.;
				if(options.outNetShortWaveRadiationDaily)		G_daily31NetShortRad[n][d]		= 0.;
                                if(options.outLongDownRadiationDaily)		G_daily31LongDownRad[n][d]		= 0.;
                                if(options.outLongUpRadiationDaily)		G_daily31LongUpRad[n][d]		= 0.;
				if(options.outNetLongWaveRadiationDaily)		G_daily31NetLongRad[n][d]		= 0.;
				if(options.outNetRadiationDaily)		G_daily31NetRadiation[n][d]		= 0.;
				// HMS radiation output adjustment start - daily radiation output
		}
	}
//}
// daily output option 31 end

// daily output option 365 start
void dailyWaterBalanceClass::daily365outInit() {
	short d;
	int n;

       // if ((5 == options.grid_store) || (6 == options.grid_store)) {
		for (n = 0; n < ng; n++)
			for (d = 0; d < 365; d++) {
                                if ((options.outSurfaceRunoffDaily)||((options.scoutSurfaceRunoff)&&(2==options.day_store))) 	G_daily365SurfaceRunoff[n][d] = 0.;

                                if ((options.outGWRechargeDaily)||((options.scoutGwRecharge)&&(2==options.day_store)))		G_daily365GwRecharge[n][d] = 0.;
          //      if ((options.outGWStorageDaily || options.outSingleStoragesDaily)||(2 == options.day_store)) 		G_daily365GwStorage[n][d] = 0.;
                                if ((options.outSoilWaterDaily)||((options.scoutSoilWater)&&(2==options.day_store))) 		G_daily365SoilWaterAvg[n][d] = 0.;
                                if ((options.outLAIDaily)||((options.scoutLAI)&&(2==options.day_store))) 			G_daily365LAI[n][d] = 0.;
                                if ((options.outLAIDaily)||((options.scoutKc)&&(2==options.day_store))) 			G_daily365Kc[n][d] = 0.;
                                if ((options.outAlbedoDaily)||((options.scoutAlbedo)&&(2==options.day_store)))			G_daily365Albedo[n][d] = 0.;
                                if ((options.outInterceptionDaily)||((options.scoutInterception)&&(2==options.day_store))) 	G_daily365Interception[n][d] = 0.;
                                if ((options.outCanopyWaterDaily)||((options.scoutCanopyWater)&&(2==options.day_store))) 	G_daily365canopyWaterContent[n][d] = 0.;
                                if ((options.outmaxCanopyWaterDaily)||((options.scoutmaxCanopyWater)&&(2==options.day_store))) G_daily365maxcanopyWaterContent[n][d] = 0.;
                if ((options.outPrecDaily)||((options.scoutPrecip)&&(2==options.day_store))) 			G_daily365Precip[n][d] = 0.;
                                if ((options.outPETDaily)||((options.scoutLandPET)&&(2==options.day_store))) 			G_daily365PET[n][d] = 0.;
                                if ((options.outTotalPETDaily)||((options.scoutCellPET)&&(2==options.day_store))) 		G_daily365TotalPET[n][d] = 0.;
                                if ((options.outSnowCoverDaily)||((options.scoutSnowCov)&&(2==options.day_store))) 		G_daily365SnowCoverAvg[n][d] = 0.;
                                if ((options.outSWEDaily)||((options.scoutSnowWater)&&(2==options.day_store))) 			G_daily365SWE[n][d] = 0.;
                                if ((options.outSnowFallDaily)||((options.scoutSnowfall)&&(2==options.day_store))) 		G_daily365SnowFall[n][d] = 0.;
                                if ((options.outSnowMeltDaily)||((options.scoutSnowmelt)&&(2==options.day_store))) 		G_daily365SnowMelt[n][d] = 0.;
                if ((options.outSingleStoragesDaily || options.outTotalWaterInStoragesDaily_km3 || options.outTotalWaterInStoragesDaily_mm)||(((options.scoutCanopyStor)||(options.scoutTWSkm3)||(options.scoutTWSmm))&&(2==options.day_store)))                G_daily365CanopyStorage[n][d] = 0.;
                if ((options.outSingleStoragesDaily || options.outTotalWaterInStoragesDaily_km3 || options.outTotalWaterInStoragesDaily_mm)||(((options.scoutSnowStor)||(options.scoutTWSkm3)||(options.scoutTWSmm))&&(2==options.day_store)))                G_daily365SnowStorage[n][d] = 0.;
                if ((options.outSingleStoragesDaily || options.outTotalWaterInStoragesDaily_km3 || options.outTotalWaterInStoragesDaily_mm)||(((options.scoutSoilStor)||(options.scoutTWSkm3)||(options.scoutTWSmm))&&(2==options.day_store)))                G_daily365SoilStorage[n][d] = 0.;
				// HMS radiation output adjustment start - daily radiation output
                                if((options.outExtRadiationDaily)||((options.scoutExtRad)&&(2==options.day_store)))		G_daily365ExtRad[n][d]		= 0.;
                                if((options.outShortDownRadiationDaily)||((options.scoutShortDownRad)&&(2==options.day_store)))		G_daily365ShortDownRad[n][d]		= 0.;
                                if((options.outShortUpRadiationDaily)||((options.scoutShortUpRad)&&(2==options.day_store)))		G_daily365ShortUpRad[n][d]		= 0.;
                                if((options.outNetShortWaveRadiationDaily)||((options.scoutNetShortRad)&&(2==options.day_store)))		G_daily365NetShortRad[n][d]		= 0.;
                                if((options.outLongDownRadiationDaily)||((options.scoutLongDownRad)&&(2==options.day_store)))		G_daily365LongDownRad[n][d]		= 0.;
                                if((options.outLongUpRadiationDaily)||((options.scoutLongUpRad)&&(2==options.day_store)))		G_daily365LongUpRad[n][d]		= 0.;
                                if((options.outNetLongWaveRadiationDaily)||((options.scoutNetLongRad)&&(2==options.day_store)))		G_daily365NetLongRad[n][d]		= 0.;
                                if((options.outNetRadiationDaily)||((options.scoutNetRad)&&(2==options.day_store)))		G_daily365NetRadiation[n][d]		= 0.;
				// HMS radiation output adjustment start - daily radiation output
                                if ((options.scoutTemp&&(2==options.day_store)))		G_daily365TempC[n][d]		= 0.;
			}
	}

//}
// daily output option 365 end


void dailyWaterBalanceClass::setParameters() {
//	alphaArid = 1.74;  // OBSOLETE now read from calibration parameter file // FP
//	alphaHumid = 1.26;  // OBSOLETE now read from calibration parameter file // FP

//	maxCanopyStoragePerLAI = 0.3; // mm  // OBSOLETE now read from calibration parameter file // FP
	// value was 0.6 mm in WaterGAP 2.1c

	canopyEvapoExp      = 0.66666666;
	runoffFracBuiltUp   = 0.5;
//	maxDailyPET_humid   = 10.0;	  // [mm/day]
//	maxDailyPET_arid    = 20.0;   // [mm/day] **mw**
//    maxDailyPET_humid   = 15.0;	  // [mm/day]   //CR 2015-09-24  // OBSOLETE now read from calibration parameter file // FP
//    maxDailyPET_arid    = 15.0;   // [mm/day]   //CR 2015-09-24  // OBSOLETE now read from calibration parameter file // FP

	openWaterAlbedo     = 0.08;
	kc_OpenWater        = 1.05;
//	snowFreezeTemp      = 0.0;	// [degree celsius]  // ADAPT to new value 2 degC for comparison of outputs! // OBSOLETE now read from calibration parameter file // FP
//	snowMeltTemp        = 0.0;	// [degree celsius]  // OBSOLETE now read from calibration parameter file // FP

	a_s = 0.25; // a_s: fraction of extraterrestrial radiation on overcast days
	b_s = 0.5; // a_s + b_s: fraction of extraterrestrial radiation on clear days

	a_c_arid = 1.35; // int-wave radiation coefficients for clear skies
	b_c_arid = -0.35; // sum has to be 1.0
	a_c_humid = 1.00;
	b_c_humid = 0.00;

//	if (options.time_series == 0 || options.time_series == 2) // only for Watch Forcing Data precip input (see BSc thesis Anja Tögl)
//		pcrit = 12.5;  // ADAPT to new value 10.0 mm/d for comparison of outputs! // OBSOLETE now read from calibration parameter file // FP
//	else
//		pcrit = 10.0;  // OBSOLETE now read from calibration parameter file // FP

}

void dailyWaterBalanceClass::setStoragesToZero() {
	for (int n = 0; n <= ng - 1; n++) {
		G_soilWaterContent[n]	= 0.;
		G_canopyWaterContent[n]	= 0.;
		G_snow[n]				= 0.;
		G_dailyStorageTransfer[n] = 0.;
        //G_groundwater[n]		= 0.;
	}
        //for (int n = 0; n < ng_land; n++) {
        for (int n = 0; n < ng; n++) {
		for (short elev = 0; elev < 101; elev++) {
			G_SnowInElevation[n][elev] = 0.;
		}
	}
}

void dailyWaterBalanceClass::readLCTdata(const char *inputDir) {

	char filename[250];

	// open LCT_22.DAT
        sprintf(filename, "%s/LCT_22.DAT", inputDir);
	ifstream lct_file(filename);

	if (!lct_file) {
		cerr << "Can not open file " << filename << " for reading.\n";
		exit(-1);
	}
	// read commentlines indicated by #
	char string[250];

	while (!lct_file.eof() && lct_file.peek() == '#') {
		lct_file.getline(string, sizeof(string));
	}

	// read values
	short dummy_int;

	for (short i = 0; i < nlct; i++) {
		lct_file >> dummy_int;
		if ((dummy_int - 1) != i) {
			cerr << "Problem with " << filename;
			exit(-1);
		}
		lct_file >> rootingDepth_lct[i]; // [-]
		lct_file >> albedo_lct[i];
		lct_file >> albedoSnow_lct[i];
		lct_file >> ddf_lct[i];
                lct_file >> emissivity_lct[i]; // land use dependent emissivity
	}
	LCTdataIsRead = true;
}

// Superbasin calculation currently obsolete, probably not working correctly in WaterGAP2.2b, do in postprocessing  // FP 2015-10
//void dailyWaterBalanceClass::normalizeDailyResults() {
//	extern routingClass routing;

//	float basinLandArea;

//	for (short i = 0; i < nSpecBasins; ++i) {
//		basinLandArea = routing.getLandAreaOfBasin(i);
//		for (short day = 0; day < 365; day++) {
//			dailyPrecSb[i][day] /= basinLandArea;
//			dailyEffPrecSb[i][day] /= basinLandArea;
//			dailySnowSb[i][day] /= basinLandArea;
//			dailyAET_Sb[i][day] /= basinLandArea;
//			dailyPET_Sb[i][day] /= basinLandArea;
//			dailyRunoffSb[i][day] /= basinLandArea;
//			dailySaturationSb[i][day] /= basinLandArea;
//			dailySunshineSb[i][day] /= basinLandArea;
//			dailyLaiSb[i][day] /= basinLandArea;
//			dailyCanopyWaterSb[i][day] /= basinLandArea;
//			dailyCanopyEvapoSb[i][day] /= basinLandArea;
//			dailyTempSb[i][day] /= basinLandArea;
//			if (dailyTempSb[i][day] < snowFreezeTemp)  // not working in WG2.2b with cell-specific calibration parameter P_T_SNOWFZ
//				snowDaysSb[i] += 1;
//		}
//		precAsSnowSb[i] /= basinLandArea;
//	}
//	normalizedFlag = true;
//}

// Superbasin calculation currently obsolete, probably not working correctly in WaterGAP2.2b, do in postprocessing  // FP 2015-10
//void dailyWaterBalanceClass::storeDailyValuesToFile(const short actualYear) {
//	extern cbasinClass cbasin;

//	if (options.outDailyValues || options.outDailyInterception || options.outSuperbasinClimate) {
//		if (!normalizedFlag)
//			normalizeDailyResults();
//	}

//	// store daily values of superbasins
//	FILE *file_ptr;
//	char filename[250];

//	if (options.outDailyValues) { // new output options 2.1f

//		for (short i = 0; i < nSpecBasins; ++i) {
//			sprintf(filename, "%s/DAILY_VALUES_%d.%s", options.output_dir, actualYear, cbasin.name[i]);
//			file_ptr = fopen(filename, "w");
//			fprintf(file_ptr, "# Daily values of superbasin '%s':\n", cbasin.name[i]);
//			fprintf(file_ptr, "# 1. column: Number of the day\n");
//			fprintf(file_ptr, "# 2. column: Precipitation [mm]\n");
//			fprintf(file_ptr, "# 3. column: Effective precipitation [mm]\n");
//			fprintf(file_ptr, "# 4. column: Amount of snow [mm]\n");
//			fprintf(file_ptr, "# 5. column: Actual evapotransiration [mm]\n");
//			fprintf(file_ptr, "# 6. column: Potential evapotranspiration [mm]\n");
//			fprintf(file_ptr, "# 7. column: Runoff [mm]\n");
//			fprintf(file_ptr, "# 8. column: Soil saturation (S/Smax) [-]\n");
//			fprintf(file_ptr, "# 9. column: Sunshine [percent]\n#\n");
//			fprintf(file_ptr, "# Day\tPrecip.\tEff. prec.\tSnow\tAET\t");
//			fprintf(file_ptr, "PET\tRunoff\tSoil Sat.\tSunshine\n");
//			for (short day = 0; day <= 364; day++)
//				fprintf(
//						file_ptr,
//						"%3d\t%10.6f\t%10.6f\t%10.6f\t%10.6f\t%10.6f\t%10.6f\t%10.6f\t%10.6f\n",
//						day + 1, dailyPrecSb[i][day], dailyEffPrecSb[i][day],
//						dailySnowSb[i][day], dailyAET_Sb[i][day],
//						dailyPET_Sb[i][day], dailyRunoffSb[i][day],
//						dailySaturationSb[i][day], dailySunshineSb[i][day] / 10.0);
//			fclose(file_ptr);
//		}
//	}
//	if (1 == options.intercept && options.outDailyInterception) { // new output options 2.1f
//		for (short i = 0; i < nSpecBasins; ++i) {
//			sprintf(filename, "%s/DAILY_INTERCEPTION_%d.%s",
//					options.output_dir, actualYear, cbasin.name[i]);
//			file_ptr = fopen(filename, "w");
//			fprintf(file_ptr, "# Daily values of superbasin '%s' [mm]:\n", cbasin.name[i]);
//			fprintf(file_ptr, "# 1. column: Number of the day\n");
//			fprintf(file_ptr, "# 2. column: Leaf area index [-]\n");
//			fprintf(file_ptr, "# 3. column: AET of canopy layer [mm]\n");
//			fprintf(file_ptr, "# 4. column: Water content of canopy layer [mm]\n");
//			fprintf(file_ptr, "# Day\tLAI\tAET\tCanopy water con.\n");
//			for (short day = 0; day <= 364; day++)
//				fprintf(file_ptr,
//						"%3d\t%10.6f\t%10.6f\t%10.6f\n",
//						day,
//						dailyLaiSb[i][day], dailyCanopyEvapoSb[i][day], dailyCanopyWaterSb[i][day]);
//			fclose(file_ptr);
//		}
//	}
//}

// Superbasin calculation currently obsolete, probably not working correctly in WaterGAP2.2b, do in postprocessing  // FP 2015-10
//void dailyWaterBalanceClass::storeAnnualValuesToFile(const short actualYear) {
//	extern cbasinClass cbasin;

//	if (!normalizedFlag)
//		normalizeDailyResults();
//	FILE *file_ptr;
//	char filename[250];

//	sprintf(filename, "%s/SUPERBASIN_CLIMATE_%d.OUT", options.output_dir, actualYear);
//	file_ptr = fopen(filename, "w");
//	fprintf(file_ptr, "# Annual values for superbasins: \n");
//	fprintf(file_ptr, "#  2. column: name of main river of superbasin\n");
//	fprintf(file_ptr, "#  3. column: potential evapotranspiration [mm]\n");
//	fprintf(file_ptr, "#  4. column: precipitation [mm]\n");
//	fprintf(file_ptr, "#  5. column: effective precipitation [mm]\n");
//	fprintf(file_ptr, "#  6. column: precipitation as snow [mm]\n");
//	fprintf(file_ptr, "#  7. column: temperature [oC]\n");
//	fprintf(file_ptr, "#  8. column: actual evapotranspiration [mm]\n");
//	fprintf(file_ptr, "#  9. column: sunshine percentage\n");
//	fprintf(file_ptr, "# 10. column: Days with average temp. below 0 (Snow days)\n#\n");
//	fprintf(file_ptr, "# No\tRiver name\tPET\tPrecipitation\t");
//	fprintf(file_ptr, "Eff. precip.\tSnow\tTemperature\t");
//	fprintf(file_ptr, "AET\tSunshine\tSnow days\n");

//	// sum up daily values to annual values
//	float total_sb_pet = -1.0;
//	float total_sb_prec = -1.0;
//	float total_sb_eff_prec = -1.0;
//	float total_sb_aet = -1.0;
//	float total_sb_sunshine = -1.0;
//	float mean_sb_temp = -1.0;

//	for (short i = 0; i < nSpecBasins; ++i) {
//		total_sb_prec = 0;
//		total_sb_eff_prec = 0;
//		total_sb_pet = 0;
//		total_sb_aet = 0;
//		total_sb_sunshine = 0;
//		for (short day = 0; day <= 364; ++day) {
//			total_sb_prec += dailyPrecSb[i][day];
//			total_sb_eff_prec += dailyEffPrecSb[i][day];
//			total_sb_pet += dailyPET_Sb[i][day];
//			total_sb_aet += dailyAET_Sb[i][day];
//			total_sb_sunshine += dailySunshineSb[i][day];
//			mean_sb_temp += dailyTempSb[i][day];
//		}
//		total_sb_sunshine /= 3650.0;
//		mean_sb_temp /= 365.0;
//		// devided by number of days and factor 10 to get percentage
//		fprintf(file_ptr,
//				"%3d\t%-15s\t%12.3f\t%12.3f\t%12.3f\t%12.3f\t%12.3f\t%12.3f\t%12.3f\t%12.3f\n",
//				(i + 1), cbasin.name[i], total_sb_pet, total_sb_prec,
//				total_sb_eff_prec, precAsSnowSb[i], mean_sb_temp,
//				total_sb_aet, total_sb_sunshine, snowDaysSb[i]);
//	}
//	fclose(file_ptr);
//}

void dailyWaterBalanceClass::writeAnnualGrids(const char *output_dir, const short year) {
	//extern gridioClass gridIO;
	//char filename[250];

	// everything is stored as monthly grids:
	// nothing to do....

	// The routine is kept. Perhaps we want to store some of the 'storages'.
}

void dailyWaterBalanceClass::writeMonthlyGrids(const char *output_dir, const short year) {
	if (2 != options.grid_store && 4 != options.grid_store && 6	!= options.grid_store)
		return;
	extern gridioClass gridIO;
	extern routingClass routing;
	extern geoClass geo;
	char filename[250];
	double cellArea;
        double storages_routing = 0.;
        int n, month;

	if (options.outPrec) { // new output options 2.2
		sprintf(filename, "%s/G_PRECIPITATION_%d.12.UNF0", output_dir, year);
		gridIO.writeUnfFile(filename, ng * 12, &G_monthlyPrecipitation[0][0]);
	}
	if (options.outInterception) { // new output options 2.1f
		sprintf(filename, "%s/G_INTERCEPTION_%d.12.UNF0", output_dir, year);
		gridIO.writeUnfFile(filename, ng * 12, &G_monthlyInterception[0][0]);
	}
	if (options.outPET) { // new output options 2.1f
		sprintf(filename, "%s/G_PET_%d.12.UNF0", output_dir, year);
		gridIO.writeUnfFile(filename, ng * 12, &G_monthlyPET[0][0]);
	}
	if (options.outTotalPET) { // new output options 2.2
		sprintf(filename, "%s/G_TOTAL_PET_%d.12.UNF0", output_dir, year);
		gridIO.writeUnfFile(filename, ng * 12, &G_monthlyTotalPET[0][0]);
	}
	if (options.outAET) { // new output options 2.1f
		sprintf(filename, "%s/G_AET_%d.12.UNF0", output_dir, year);
		gridIO.writeUnfFile(filename, ng * 12, &G_monthlyAET[0][0]);
	}
	if (options.outCellAET) { // new output options 2.1f corrected land evapo
		sprintf(filename, "%s/G_LAND_AET_%d.12.UNF0", output_dir, year);
		gridIO.writeUnfFile(filename, ng * 12, &G_monthlyLandAET[0][0]);
		sprintf(filename, "%s/G_LAND_AET_UNCORR_%d.12.UNF0", output_dir, year);
		gridIO.writeUnfFile(filename, ng * 12, &G_monthlyLandAET_uncorr[0][0]);
	}
	if (options.outRunoff) { // new output options 2.1f
		sprintf(filename, "%s/G_RUNOFF_%d.12.UNF0", output_dir, year);
		gridIO.writeUnfFile(filename, ng * 12, &G_monthlyRunoff[0][0]);
	}
	if (options.outUrbanRunoff) { // new output options 2.2
		sprintf(filename, "%s/G_URBAN_RUNOFF_%d.12.UNF0", output_dir, year);
		gridIO.writeUnfFile(filename, ng * 12, &G_monthlyUrbanRunoff[0][0]);
	}
	if (options.outSurfaceRunoff) { // new output options 2.2
		sprintf(filename, "%s/G_SURFACE_RUNOFF_%d.12.UNF0", output_dir, year);
		gridIO.writeUnfFile(filename, ng * 12, &G_monthlySurfaceRunoff[0][0]);
	}


    if (options.outGWRecharge) { // new output options 2.2 // related to continental area
        sprintf(filename, "%s/G_GW_RECHARGE_mm_%d.12.UNF0", output_dir, year);
        gridIO.writeUnfFile(filename, ng * 12, &G_monthlyGwRecharge[0][0]);

        for (n = 0; n < ng; n++)
           for (short m = 0; m < 12; m++)
              G_monthlyArray[n][m] = G_monthlyGwRecharge[n][m]
                 * (geo.areaOfCellByArrayPos(n) * (geo.G_contfreq[n] / 100.)) / 1000000.0;
        sprintf(filename, "%s/G_GW_RECHARGE_km3_%d.12.UNF0", output_dir, year);
        gridIO.writeUnfFile(filename, ng * 12, &G_monthlyArray[0][0]);
	}
	if (options.outOpenWaterPET) { // new output options 2.1f
		sprintf(filename, "%s/G_OPEN_WATER_PET_%d.12.UNF0", output_dir, year);
		gridIO.writeUnfFile(filename, ng * 12, &G_monthlyOpenWaterPET[0][0]);
	}


//	int n, month;
	unsigned char numberOfDaysInMonth[12] = { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };
	for (n = 0; n <= ng - 1; n++) {
		for (month = 0; month <= 11; month++) {
			G_monthlySWE[n][month] /= numberOfDaysInMonth[month];
			G_monthlySnowCoverAvg[n][month] /= numberOfDaysInMonth[month];
			G_monthlySoilWaterAvg[n][month] /= numberOfDaysInMonth[month];
			G_monthlyNetRadiation[n][month] /= numberOfDaysInMonth[month]; // HMS radiation output adjustment
                        G_monthlyShortDownRad[n][month] /= numberOfDaysInMonth[month];  // HMS radiation output adjustment
			G_monthlyExtRad[n][month] /= numberOfDaysInMonth[month];  // HMS radiation output adjustment
                        G_monthlyShortUpRad[n][month] /= numberOfDaysInMonth[month];  // HMS radiation output adjustment
			G_monthlyNetShortRad[n][month] /= numberOfDaysInMonth[month];  // HMS radiation output adjustment
                        G_monthlyLongDownRad[n][month] /= numberOfDaysInMonth[month];  // HMS radiation output adjustment
                        G_monthlyLongUpRad[n][month] /= numberOfDaysInMonth[month];  // HMS radiation output adjustment
			G_monthlyNetLongRad[n][month] /= numberOfDaysInMonth[month];  // HMS radiation output adjustment
			
                        G_monthlyTemperature[n][month] /= numberOfDaysInMonth[month];  // monthly temperature output
		
			G_monthlySunshine[n][month] /= numberOfDaysInMonth[month];  // monthly sunshine output
		
		
			G_monthlyLAI[n][month] /= numberOfDaysInMonth[month];
			G_monthlyAlbedo[n][month] /= numberOfDaysInMonth[month];
			G_monthlycanopyWaterContent[n][month] /= numberOfDaysInMonth[month];
			G_monthlymaxcanopyWaterContent[n][month] /= numberOfDaysInMonth[month];

                        // GRACE application - divide by days in month to get monthly mean values
                        // routing.G_monthlyGwStorage[n][month] /= numberOfDaysInMonth[month]; // all calculations of G_monthlyGwStorage are moved to routing.cpp
			if (options.grid_store_TypeForStorages == 1) {
                G_monthlyCanopyStorage[n][month] /= numberOfDaysInMonth[month];
                G_monthlySnowStorage[n][month] /= numberOfDaysInMonth[month];
                G_monthlySoilStorage[n][month] /= numberOfDaysInMonth[month];
            }

                        // Storage summation restructured
                        storages_routing = routing.G_monthlyLocLakeStorage[n][month] // all values are in km3
			         + routing.G_monthlyLocWetlStorage[n][month]
			         + routing.G_monthlyGloLakeStorage[n][month]
			         + routing.G_monthlyGloWetlStorage[n][month]
			         + routing.G_monthlyRiverStorage[n][month]
                                 + routing.G_monthlyGwStorage[n][month];

			if (options.resOpt)
				storages_routing += routing.G_monthlyResStorage[n][month];

			if (options.grid_store_TypeForStorages == 1)
				storages_routing /= numberOfDaysInMonth[month];

			//add all waterstorages and subtract wateruse         // [km³]
			G_monthly_totalWaterInStorages[n][month] =
					G_monthlyCanopyStorage[n][month]
					+ G_monthlySnowStorage[n][month]
					+ G_monthlySoilStorage[n][month]
					+ storages_routing;
                        // end: Storage summation restructured

			//convert total storage from km³ to mm
			G_monthly_totalWaterInStorages_mm[n][month] = 
					G_monthly_totalWaterInStorages[n][month]
                    / (geo.areaOfCellByArrayPos(n) * (geo.G_contfreq[n] / 100.)) * 1000000;
                }
	}

        if (options.outSingleStorages) { // new output options 2.2  all storages are monthly means [km3]
            // FP20161010N001 Replacement: geo.G_contfreq[n] instead (geo.G_landfreq[n] + geo.G_fwaterfreq[n])

            for (n = 0; n < ng; n++)
                for (short m = 0; m < 12; m++)
                    G_monthlyArray[n][m] = G_monthlyCanopyStorage[n][m]; // TOTAL_STORAGES - conversion double -> float
            sprintf(filename, "%s/G_CANOPY_WATER_STORAGE_km3_%d.12.UNF0", output_dir, year);
            gridIO.writeUnfFile(filename, ng * 12, &G_monthlyArray[0][0]);

            for (n = 0; n < ng; n++)
                for (short m = 0; m < 12; m++)
                    if (0 == geo.G_contcell[n]) // ocean grid cells within land mask (e.g. Caspian Sea)
                        G_monthlyArray[n][m] = 0.;
                    else
                    G_monthlyArray[n][m] = G_monthlyCanopyStorage[n][m] // TOTAL_STORAGES - conversion double -> float
                            / (geo.areaOfCellByArrayPos(n) * (geo.G_contfreq[n] / 100.)) * 1000000;
            sprintf(filename, "%s/G_CANOPY_WATER_STORAGE_mm_%d.12.UNF0", output_dir, year);
            gridIO.writeUnfFile(filename, ng * 12, &G_monthlyArray[0][0]);

            for (n = 0; n < ng; n++)
                for (short m = 0; m < 12; m++)
                    G_monthlyArray[n][m] = G_monthlySnowStorage[n][m]; // TOTAL_STORAGES - conversion double -> float
            sprintf(filename, "%s/G_SNOW_WATER_STORAGE_km3_%d.12.UNF0", output_dir, year);
            gridIO.writeUnfFile(filename, ng * 12, &G_monthlyArray[0][0]);

            for (n = 0; n < ng; n++)
                for (short m = 0; m < 12; m++)
                    if (0 == geo.G_contcell[n]) // ocean grid cells within land mask (e.g. Caspian Sea)
                        G_monthlyArray[n][m] = 0.;
                    else
                    G_monthlyArray[n][m] = G_monthlySnowStorage[n][m] // TOTAL_STORAGES - conversion double -> float
                            / (geo.areaOfCellByArrayPos(n) * (geo.G_contfreq[n] / 100.)) * 1000000;
            sprintf(filename, "%s/G_SNOW_WATER_STORAGE_mm_%d.12.UNF0", output_dir, year);
            gridIO.writeUnfFile(filename, ng * 12, &G_monthlyArray[0][0]);

            for (n = 0; n < ng; n++)
                for (short m = 0; m < 12; m++)
                    G_monthlyArray[n][m] = G_monthlySoilStorage[n][m]; // TOTAL_STORAGES - conversion double -> float
            sprintf(filename, "%s/G_SOIL_WATER_STORAGE_km3_%d.12.UNF0", output_dir, year);
            gridIO.writeUnfFile(filename, ng * 12, &G_monthlyArray[0][0]);

            for (n = 0; n < ng; n++)
                for (short m = 0; m < 12; m++)
                    if (0 == geo.G_contcell[n]) // ocean grid cells within land mask (e.g. Caspian Sea)
                        G_monthlyArray[n][m] = 0.;
                    else
                    G_monthlyArray[n][m] = G_monthlySoilStorage[n][m] // TOTAL_STORAGES - conversion double -> float
                            / (geo.areaOfCellByArrayPos(n) * (geo.G_contfreq[n] / 100.)) * 1000000;
            sprintf(filename, "%s/G_SOIL_WATER_STORAGE_mm_%d.12.UNF0", output_dir, year);
            gridIO.writeUnfFile(filename, ng * 12, &G_monthlyArray[0][0]);

        }

        // write total monthly storage to file [km3]
        if (options.outTotalWaterInStorages_km3) {
		for (n = 0; n < ng; n++)
			for (short m = 0; m < 12; m++)
                                G_monthlyArray[n][m] = G_monthly_totalWaterInStorages[n][m]; // TOTAL_STORAGES - conversion double -> float
		sprintf(filename, "%s/G_TOTAL_STORAGES_km3_%d.12.UNF0", output_dir,	year);
		gridIO.writeUnfFile(filename, ng * 12, &G_monthlyArray[0][0]);
        }
        // [mm]
	if (options.outTotalWaterInStorages_mm) {
		for (n = 0; n < ng; n++)
			for (short m = 0; m < 12; m++)
                if (0 == geo.G_contcell[n]) // ocean grid cells within land mask (e.g. Caspian Sea)
                    G_monthlyArray[n][m] = 0.;
                else
                    G_monthlyArray[n][m] = G_monthly_totalWaterInStorages_mm[n][m]; // TOTAL_STORAGES in mm - conversion double -> float
		sprintf(filename, "%s/G_TOTAL_STORAGES_mm_%d.12.UNF0", output_dir, year);
		gridIO.writeUnfFile(filename, ng * 12, &G_monthlyArray[0][0]);
	}

	if (options.outLAI) { // new output options 2.2
		sprintf(filename, "%s/G_LAI_%d.12.UNF0", output_dir, year);
		gridIO.writeUnfFile(filename, ng * 12, &G_monthlyLAI[0][0]);
	}
	if (options.outAlbedo) { // new output options 2.2
		sprintf(filename, "%s/G_ALBEDO_%d.12.UNF0", output_dir, year);
		gridIO.writeUnfFile(filename, ng * 12, &G_monthlyAlbedo[0][0]);
	}
	if (options.outCanopyWater) { // new output options 2.2
		sprintf(filename, "%s/G_CANOPY_WATER_mm_%d.12.UNF0", output_dir, year);
		gridIO.writeUnfFile(filename, ng * 12, &G_monthlycanopyWaterContent[0][0]);
	}
	if (options.outmaxCanopyWater) { // new output options 2.2
		sprintf(filename, "%s/G_MAX_CANOPY_WATER_%d.12.UNF0", output_dir, year);
		gridIO.writeUnfFile(filename, ng * 12, &G_monthlymaxcanopyWaterContent[0][0]);
	}
	if (options.outSoilWater) { // new output options 2.1f
		sprintf(filename, "%s/G_SOIL_WATER_%d.12.UNF0", output_dir, year);
		gridIO.writeUnfFile(filename, ng * 12, &G_monthlySoilWaterAvg[0][0]);
	}
	if (options.outSnowCover) { // new output options 2.2
		sprintf(filename, "%s/G_SNOW_COVER_FRAC_%d.12.UNF0", output_dir, year);
		gridIO.writeUnfFile(filename, ng * 12, &G_monthlySnowCoverAvg[0][0]);
	}
	if (options.outSWE) { // new output options 2.2
		sprintf(filename, "%s/G_SNOW_WATER_EQUIV_%d.12.UNF0", output_dir, year);
		gridIO.writeUnfFile(filename, ng * 12, &G_monthlySWE[0][0]);
	}
	if (options.outSnowFall) { // new output options 2.2
		sprintf(filename, "%s/G_SNOW_FALL_%d.12.UNF0", output_dir, year);
		gridIO.writeUnfFile(filename, ng * 12, &G_monthlySnowFall[0][0]);
	}
    if (options.outSnowFall) { // new output option 2.2b (from 2.2 ISIMIP) // FP2061011N001
        sprintf(filename, "%s/G_RAIN_FALL_%d.12.UNF0", output_dir, year);
        gridIO.writeUnfFile(filename, ng * 12, &G_monthlyRainFall[0][0]);
    }
    if (options.outSnowMelt) { // new output options 2.1f
		sprintf(filename, "%s/G_SNOW_MELT_%d.12.UNF0", output_dir, year);
		gridIO.writeUnfFile(filename, ng * 12, &G_monthlySnowMelt[0][0]);
	}
	if (options.outSnowEvap) { // new output options 2.2
		sprintf(filename, "%s/G_SNOW_EVAP_%d.12.UNF0", output_dir, year);
		gridIO.writeUnfFile(filename, ng * 12, &G_monthlySnowEvap[0][0]);
	}

	if (options.outNetRadiation) { // new output options 2.1f
		sprintf(filename, "%s/G_NET_RADIATION_%d.12.UNF0", output_dir, year);
		gridIO.writeUnfFile(filename, ng * 12, &G_monthlyNetRadiation[0][0]);
	}
        if (options.outNetShortWaveRadiation) { // net_short_wave_rad output
                sprintf(filename, "%s/G_NET_SHORTWAVE_RAD_%d.12.UNF0", output_dir, year);
		gridIO.writeUnfFile(filename, ng * 12, &G_monthlyNetShortRad[0][0]);
		
		sprintf(filename, "%s/G_EXT_RAD_%d.12.UNF0", output_dir, year);
		gridIO.writeUnfFile(filename, ng * 12, &G_monthlyExtRad[0][0]);
		
                sprintf(filename, "%s/G_SHORTWAVE_UP_RAD_%d.12.UNF0", output_dir, year);
                gridIO.writeUnfFile(filename, ng * 12, &G_monthlyShortUpRad[0][0]);

                sprintf(filename, "%s/G_SHORTWAVE_DOWN_RAD_%d.12.UNF0", output_dir, year);
                gridIO.writeUnfFile(filename, ng * 12, &G_monthlyShortDownRad[0][0]);
	}
        if (options.outNetLongWaveRadiation) { // net_long_wave_rad output
                sprintf(filename, "%s/G_NET_LONGWAVE_RAD_%d.12.UNF0", output_dir, year);
		gridIO.writeUnfFile(filename, ng * 12, &G_monthlyNetLongRad[0][0]);
  		
                sprintf(filename, "%s/G_LONGWAVE_DOWN_RAD_%d.12.UNF0", output_dir, year);
                gridIO.writeUnfFile(filename, ng * 12, &G_monthlyLongDownRad[0][0]);
		
                sprintf(filename, "%s/G_LONGWAVE_UP_RAD_%d.12.UNF0", output_dir, year);
                gridIO.writeUnfFile(filename, ng * 12, &G_monthlyLongUpRad[0][0]);
	
	}
        if (options.outSunshine) {
        sprintf(filename, "%s/G_SUNSHINE_%d.12.UNF0", output_dir, year);   // monthly sunshine output in % (0 if not driven with CRU data)
        gridIO.writeUnfFile(filename, ng * 12, &G_monthlySunshine[0][0]);
        }
        if (options.outTemperature) {
        sprintf(filename, "%s/G_TEMPERATURE_%d.12.UNF0", output_dir, year);  // monthly temperature output in ?C
        gridIO.writeUnfFile(filename, ng * 12, &G_monthlyTemperature[0][0]);
        }
}

// daily output option 31 start

void dailyWaterBalanceClass::writeDaily31Grids(char *output_dir, const int year, const int month) {
	if (3 != options.grid_store && 4 != options.grid_store)
		return;
	extern gridioClass gridIO;
	char filename[250];

	if (options.outSurfaceRunoffDaily) { // new output options 2.2
		sprintf(filename, "%s/G_SURFACE_RUNOFF_%d_%d.31.UNF0", output_dir, year, month + 1);
		gridIO.writeUnfFile(filename, ng * 31, &G_daily31SurfaceRunoff[0][0]);
	}

	if (options.outGWRechargeDaily) { // new output options 2.2
		sprintf(filename, "%s/G_GW_RECHARGE_%d_%d.31.UNF0", output_dir, year, month + 1);
		gridIO.writeUnfFile(filename, ng * 31, &G_daily31GwRecharge[0][0]);
	}

	if (options.outSoilWaterDaily) { // new output options 2.2
		sprintf(filename, "%s/G_SOIL_WATER_%d_%d.31.UNF0", output_dir, year, month + 1);
		gridIO.writeUnfFile(filename, ng * 31, &G_daily31SoilWaterAvg[0][0]);
	}
	if (options.outLAIDaily) { // new output options 2.2
		sprintf(filename, "%s/G_LAI_%d_%d.31.UNF0", output_dir, year, month + 1);
		gridIO.writeUnfFile(filename, ng * 31, &G_daily31LAI[0][0]);
	}
	if ((options.outLAIDaily) && (options.use_kc == 1)) { // new output options 2.2
		sprintf(filename, "%s/G_KC_%d_%d.31.UNF0", output_dir, year, month + 1);
		gridIO.writeUnfFile(filename, ng * 31, &G_daily31Kc[0][0]);
	}
	if (options.outAlbedoDaily) { // new output options 2.2
		sprintf(filename, "%s/G_ALBEDO_%d_%d.31.UNF0", output_dir, year, month	+ 1);
		gridIO.writeUnfFile(filename, ng * 31, &G_daily31Albedo[0][0]);
	}
	if (options.outInterceptionDaily) { // new output options 2.2
		sprintf(filename, "%s/G_INTERCEPTION_%d_%d.31.UNF0", output_dir, year,	month + 1);
		gridIO.writeUnfFile(filename, ng * 31, &G_daily31Interception[0][0]);
	}
	if (options.outCanopyWaterDaily) { // new output options 2.2
		sprintf(filename, "%s/G_CANOPY_WATER_mm_%d_%d.31.UNF0", output_dir,	year, month + 1);
		gridIO.writeUnfFile(filename, ng * 31, &G_daily31canopyWaterContent[0][0]);
	}
	if (options.outmaxCanopyWaterDaily) { // new output options 2.2
		sprintf(filename, "%s/G_MAX_CANOPY_WATER_%d_%d.31.UNF0", output_dir, year, month + 1);
		gridIO.writeUnfFile(filename, ng * 31, &G_daily31maxcanopyWaterContent[0][0]);
	}
        if (options.outPrecDaily) { // daily precip output
		sprintf(filename, "%s/G_PRECIPITATION_%d_%d.31.UNF0", output_dir, year,	month + 1);
		gridIO.writeUnfFile(filename, ng * 31, &G_daily31Precip[0][0]);
    }
	if (options.outPETDaily) { // new output options 2.2
		sprintf(filename, "%s/G_PET_%d_%d.31.UNF0", output_dir, year, month + 1);
		gridIO.writeUnfFile(filename, ng * 31, &G_daily31PET[0][0]);
	}
	if (options.outTotalPETDaily) { // new output options 2.2
		sprintf(filename, "%s/G_TOTAL_PET_%d_%d.31.UNF0", output_dir, year,	month + 1);
		gridIO.writeUnfFile(filename, ng * 31, &G_daily31TotalPET[0][0]);
	}
	if (options.outSnowCoverDaily) { // new output options 2.2
		sprintf(filename, "%s/G_SNOW_COVER_FRAC_%d_%d.31.UNF0", output_dir,	year, month + 1);
		gridIO.writeUnfFile(filename, ng * 31, &G_daily31SnowCoverAvg[0][0]);
	}
	if (options.outSWEDaily) { // new output options 2.2
		sprintf(filename, "%s/G_SNOW_WATER_EQUIV_%d_%d.31.UNF0", output_dir, year, month + 1);
		gridIO.writeUnfFile(filename, ng * 31, &G_daily31SWE[0][0]);
	}
	if (options.outSnowFallDaily) { // new output options 2.2
		sprintf(filename, "%s/G_SNOW_FALL_%d_%d.31.UNF0", output_dir, year,	month + 1);
		gridIO.writeUnfFile(filename, ng * 31, &G_daily31SnowFall[0][0]);
	}
	if (options.outSnowMeltDaily) { // new output options 2.2
		sprintf(filename, "%s/G_SNOW_MELT_%d_%d.31.UNF0", output_dir, year,	month + 1);
		gridIO.writeUnfFile(filename, ng * 31, &G_daily31SnowMelt[0][0]);
	}
        if (options.outSingleStoragesDaily) { // new output options 2.2
                sprintf(filename, "%s/G_CANOPY_WATER_STORAGE_km3_%d_%d.31.UNF0", output_dir, year,	month + 1);
                gridIO.writeUnfFile(filename, ng * 31, &G_daily31CanopyStorage[0][0]);
        }
        if (options.outSingleStoragesDaily) { // new output options 2.2
                sprintf(filename, "%s/G_SNOW_WATER_STORAGE_km3_%d_%d.31.UNF0", output_dir, year,	month + 1);
                gridIO.writeUnfFile(filename, ng * 31, &G_daily31SnowStorage[0][0]);
        }
        if (options.outSingleStoragesDaily) { // new output options 2.2
                sprintf(filename, "%s/G_SOIL_WATER_STORAGE_km3_%d_%d.31.UNF0", output_dir, year,	month + 1);
                gridIO.writeUnfFile(filename, ng * 31, &G_daily31SoilStorage[0][0]);
        }
	// HMS radiation output adjustment start - daily rad output anyway if daily outputs are selected.
        if (options.outExtRadiationDaily) {
		sprintf(filename, "%s/G_EXT_RAD_%d_%d.31.UNF0", output_dir, year, month+1);
		gridIO.writeUnfFile(filename, ng * 31, &G_daily31ExtRad[0][0]);
  }
        if (options.outShortDownRadiationDaily) {
                sprintf(filename, "%s/G_SHORTWAVE_DOWN_RAD_%d_%d.31.UNF0", output_dir, year, month+1);
                gridIO.writeUnfFile(filename, ng * 31, &G_daily31ShortDownRad[0][0]);
  }
        if (options.outShortUpRadiationDaily) {
                sprintf(filename, "%s/G_SHORTWAVE_UP_RAD_%d_%d.31.UNF0", output_dir, year, month+1);
                gridIO.writeUnfFile(filename, ng * 31, &G_daily31ShortUpRad[0][0]);
  }
	if (options.outNetShortWaveRadiationDaily) {		
                sprintf(filename, "%s/G_NET_SHORTWAVE_RAD_%d_%d.31.UNF0", output_dir, year, month+1);
		gridIO.writeUnfFile(filename, ng * 31, &G_daily31NetShortRad[0][0]);
  }
        if (options.outLongDownRadiationDaily) {
                sprintf(filename, "%s/G_LONGWAVE_DOWN_RAD_%d_%d.31.UNF0", output_dir, year, month+1);
                gridIO.writeUnfFile(filename, ng * 31, &G_daily31LongDownRad[0][0]);
  }
        if (options.outLongUpRadiationDaily) {
                sprintf(filename, "%s/G_LONGWAVE_UP_RAD_%d_%d.31.UNF0", output_dir, year, month+1);
                gridIO.writeUnfFile(filename, ng * 31, &G_daily31LongUpRad[0][0]);
  }
	if (options.outNetLongWaveRadiationDaily) {		
		sprintf(filename, "%s/G_NET_LONG_WAVE_RAD_%d_%d.31.UNF0", output_dir, year, month+1);
		gridIO.writeUnfFile(filename, ng * 31, &G_daily31NetLongRad[0][0]);
  }
	if (options.outNetRadiationDaily) {		
		sprintf(filename, "%s/G_NET_RADIATION_%d_%d.31.UNF0", output_dir, year, month+1);
		gridIO.writeUnfFile(filename, ng * 31, &G_daily31NetRadiation[0][0]);
	}	
    // HMS radiation output adjustment end
}
// daily output option 31 end

// daily output option 365 start
void dailyWaterBalanceClass::writeDaily365Grids(char *output_dir, const int year) {
	if (5 != options.grid_store && 6 != options.grid_store)
		return;
	extern gridioClass gridIO;
    extern routingClass routing;
    extern geoClass geo;
	char filename[250];

	if (options.outSurfaceRunoffDaily) { // new output options 2.2
		sprintf(filename, "%s/G_SURFACE_RUNOFF_%d.365.UNF0", output_dir, year);
		gridIO.writeUnfFile(filename, ng * 365, &G_daily365SurfaceRunoff[0][0]);
	}

	if (options.outGWRechargeDaily) { // new output options 2.2
		sprintf(filename, "%s/G_GW_RECHARGE_%d.365.UNF0", output_dir, year);
		gridIO.writeUnfFile(filename, ng * 365, &G_daily365GwRecharge[0][0]);
	}

	if (options.outSoilWaterDaily) { // new output options 2.2
		sprintf(filename, "%s/G_SOIL_WATER_%d.365.UNF0", output_dir, year);
		gridIO.writeUnfFile(filename, ng * 365, &G_daily365SoilWaterAvg[0][0]);
	}
	if (options.outLAIDaily) { // new output options 2.2
		sprintf(filename, "%s/G_LAI_%d.365.UNF0", output_dir, year);
		gridIO.writeUnfFile(filename, ng * 365, &G_daily365LAI[0][0]);
	}
	if ((options.outLAIDaily) && (options.use_kc == 1)) { // new output options 2.2
		sprintf(filename, "%s/G_KC_%d.365.UNF0", output_dir, year);
		gridIO.writeUnfFile(filename, ng * 365, &G_daily365Kc[0][0]);
	}
	if (options.outAlbedoDaily) { // new output options 2.2
		sprintf(filename, "%s/G_ALBEDO_%d.365.UNF0", output_dir, year);
		gridIO.writeUnfFile(filename, ng * 365, &G_daily365Albedo[0][0]);
	}
	if (options.outInterceptionDaily) { // new output options 2.2
		sprintf(filename, "%s/G_INTERCEPTION_%d.365.UNF0", output_dir, year);
		gridIO.writeUnfFile(filename, ng * 365, &G_daily365Interception[0][0]);
	}
	if (options.outCanopyWaterDaily) { // new output options 2.2
		sprintf(filename, "%s/G_CANOPY_WATER_mm_%d.365.UNF0", output_dir, year);
		gridIO.writeUnfFile(filename, ng * 365, &G_daily365canopyWaterContent[0][0]);
	}
	if (options.outmaxCanopyWaterDaily) { // new output options 2.2
		sprintf(filename, "%s/G_MAX_CANOPY_WATER_%d.365.UNF0", output_dir, year);
		gridIO.writeUnfFile(filename, ng * 365,	&G_daily365maxcanopyWaterContent[0][0]);
	}
        if (options.outPrecDaily) { // daily precip output
		sprintf(filename, "%s/G_PRECIPITATION_%d.365.UNF0", output_dir, year);
		gridIO.writeUnfFile(filename, ng * 365, &G_daily365Precip[0][0]);
	}
	if (options.outPETDaily) { // new output options 2.2
		sprintf(filename, "%s/G_PET_%d.365.UNF0", output_dir, year);
		gridIO.writeUnfFile(filename, ng * 365, &G_daily365PET[0][0]);
	}
	if (options.outTotalPETDaily) { // new output options 2.2
		sprintf(filename, "%s/G_TOTAL_PET_%d.365.UNF0", output_dir, year);
		gridIO.writeUnfFile(filename, ng * 365, &G_daily365TotalPET[0][0]);
	}
	if (options.outSnowCoverDaily) { // new output options 2.2
		sprintf(filename, "%s/G_SNOW_COVER_FRAC_%d.365.UNF0", output_dir, year);
		gridIO.writeUnfFile(filename, ng * 365, &G_daily365SnowCoverAvg[0][0]);
	}
	if (options.outSWEDaily) { // new output options 2.2
		sprintf(filename, "%s/G_SNOW_WATER_EQUIV_%d.365.UNF0", output_dir, year);
		gridIO.writeUnfFile(filename, ng * 365, &G_daily365SWE[0][0]);
	}
	if (options.outSnowFallDaily) { // new output options 2.2
		sprintf(filename, "%s/G_SNOW_FALL_%d.365.UNF0", output_dir, year);
		gridIO.writeUnfFile(filename, ng * 365, &G_daily365SnowFall[0][0]);
	}
	if (options.outSnowMeltDaily) { // new output options 2.2
		sprintf(filename, "%s/G_SNOW_MELT_%d.365.UNF0", output_dir, year);
		gridIO.writeUnfFile(filename, ng * 365, &G_daily365SnowMelt[0][0]);
    }
 
    if (options.outSingleStoragesDaily) { // new output options 2.2
        sprintf(filename, "%s/G_CANOPY_WATER_STORAGE_km3_%d.365.UNF0", output_dir, year);
        gridIO.writeUnfFile(filename, ng * 365, &G_daily365CanopyStorage[0][0]);
    }
    if (options.outSingleStoragesDaily) { // new output options 2.2
        sprintf(filename, "%s/G_SNOW_WATER_STORAGE_km3_%d.365.UNF0", output_dir, year);
        gridIO.writeUnfFile(filename, ng * 365, &G_daily365SnowStorage[0][0]);
    }
    if (options.outSingleStoragesDaily) { // new output options 2.2
        sprintf(filename, "%s/G_SOIL_WATER_STORAGE_km3_%d.365.UNF0", output_dir, year);
        gridIO.writeUnfFile(filename, ng * 365, &G_daily365SoilStorage[0][0]);
    }
	// HMS radiation output adjustment start - daily rad output anyway if daily outputs are selected.

        if (options.outExtRadiationDaily) {
		sprintf(filename, "%s/G_EXT_RAD_%d.365.UNF0", output_dir, year);
		gridIO.writeUnfFile(filename, ng * 365, &G_daily365ExtRad[0][0]);
  	}
        if (options.outShortDownRadiationDaily) {
                sprintf(filename, "%s/G_SHORTWAVE_DOWN_RAD_%d.365.UNF0", output_dir, year);
                gridIO.writeUnfFile(filename, ng * 365, &G_daily365ShortDownRad[0][0]);
  	}
        if (options.outShortUpRadiationDaily) {
                sprintf(filename, "%s/G_SHORTWAVE_UP_RAD_%d.365.UNF0", output_dir, year);
                gridIO.writeUnfFile(filename, ng * 365, &G_daily365ShortUpRad[0][0]);
  	}
	if (options.outNetShortWaveRadiationDaily) {		
                sprintf(filename, "%s/G_NET_SHORTWAVE_RAD_%d.365.UNF0", output_dir, year);
		gridIO.writeUnfFile(filename, ng * 365, &G_daily365NetShortRad[0][0]);
  	}
        if (options.outLongDownRadiationDaily) {
                sprintf(filename, "%s/G_LONGWAVE_DOWN_RAD_%d.365.UNF0", output_dir, year);
                gridIO.writeUnfFile(filename, ng * 365, &G_daily365LongDownRad[0][0]);
  	}
        if (options.outLongUpRadiationDaily) {
                sprintf(filename, "%s/G_LONGWAVE_UP_RAD_%d.365.UNF0", output_dir, year);
                gridIO.writeUnfFile(filename, ng * 365, &G_daily365LongUpRad[0][0]);
  	}
	if (options.outNetLongWaveRadiationDaily) {		
                sprintf(filename, "%s/G_NET_LONGWAVE_RAD_%d.365.UNF0", output_dir, year);
		gridIO.writeUnfFile(filename, ng * 365, &G_daily365NetLongRad[0][0]);
  	}
	if (options.outNetRadiationDaily) {		
		sprintf(filename, "%s/G_NET_RADIATION_%d.365.UNF0", output_dir, year);
		gridIO.writeUnfFile(filename, ng * 365, &G_daily365NetRadiation[0][0]);
	}	
    // HMS radiation output adjustment end       
}
// daily output option 365 end

// Superbasin calculation currently obsolete, probably not working correctly in WaterGAP2.2b, do in postprocessing  // FP 2015-10
//void dailyWaterBalanceClass::getPETVector(vector<float>&annualPETVector) {
//	if (!normalizedFlag)
//		normalizeDailyResults();

//	float annualPET;

//	for (short n = 0; n < nSpecBasins; ++n) {
//		annualPET = 0;
//		for (short i = 0; i <= 364; i++) {
//			annualPET += dailyPET_Sb[n][i];
//		}
//		annualPETVector.push_back(annualPET);
//	}
//}
// Superbasin calculation currently obsolete, probably not working correctly in WaterGAP2.2b, do in postprocessing  // FP 2015-10
//void dailyWaterBalanceClass::getAETVector(vector<float>&annualAETVector) {
//	if (!normalizedFlag)
//		normalizeDailyResults();

//	float annualAET;

//	for (short n = 0; n < nSpecBasins; ++n) {
//		annualAET = 0;
//		for (short i = 0; i <= 364; i++) {
//			annualAET += dailyAET_Sb[n][i];
//		}
//		annualAETVector.push_back(annualAET);
//	}
//}
// Superbasin calculation currently obsolete, probably not working correctly in WaterGAP2.2b, do in postprocessing  // FP 2015-10
//void dailyWaterBalanceClass::getInterceptVector(vector<float>&annualInterceptVector) {
//	if (!normalizedFlag)
//		normalizeDailyResults();

//	float annualIntercept;

//	for (short n = 0; n < nSpecBasins; ++n) {
//		annualIntercept = 0;
//		for (short i = 0; i <= 364; i++) {
//			annualIntercept += dailyCanopyEvapoSb[n][i];
//		}
//		annualInterceptVector.push_back(annualIntercept);
//	}
//}
// Superbasin calculation currently obsolete, probably not working correctly in WaterGAP2.2b, do in postprocessing  // FP 2015-10
//void dailyWaterBalanceClass::getRunoffVector(vector<float>&annualRunoffVector) {
//	if (!normalizedFlag)
//		normalizeDailyResults();

//	float annualRunoff;

//	for (short n = 0; n < nSpecBasins; ++n) {
//		annualRunoff = 0;
//		for (short i = 0; i <= 364; i++) {
//			annualRunoff += dailyRunoffSb[n][i];
//		}
//		annualRunoffVector.push_back(annualRunoff);
//	}
//}
// Superbasin calculation currently obsolete, probably not working correctly in WaterGAP2.2b, do in postprocessing  // FP 2015-10
//void dailyWaterBalanceClass::getSnowVector(vector<float>&annualSnowVector) {
//	if (!normalizedFlag)
//		normalizeDailyResults();

//	float annualSnow;

//	for (short n = 0; n < nSpecBasins; ++n) {
//		annualSnow = 0;
//		for (short i = 0; i <= 364; i++) {
//			annualSnow += dailySnowSb[n][i];
//		}
//		annualSnow /= 365.0;
//		annualSnowVector.push_back(annualSnow);
//	}
//}
// Superbasin calculation currently obsolete, probably not working correctly in WaterGAP2.2b, do in postprocessing  // FP 2015-10
//void dailyWaterBalanceClass::getSoilSatVector(vector<float>&annualSoilSatVector) {
//	if (!normalizedFlag)
//		normalizeDailyResults();

//	float annualSoilSat;

//	for (short n = 0; n < nSpecBasins; ++n) {
//		annualSoilSat = 0;
//		for (short i = 0; i <= 364; i++) {
//			annualSoilSat += dailySaturationSb[n][i];
//		}
//		annualSoilSat /= 365.0;
//		annualSoilSatVector.push_back(annualSoilSat);
//	}
//}

//=====daily temperature variation==============
//=== not used in 2.1f because effects are negligible
/*
 void dailyWaterBalanceClass::readDailyTempVarFile()
 {
 //float DailyTempVar;
 float tempvar;

 //int m;
 short day;
 char filename[250];
 FILE *file_ptr;


 sprintf(filename, "%s/TempVar_25.dat", options.input_dir);
 file_ptr = fopen(filename, "r");
 if (file_ptr != NULL) {

 //gridIO.readUnfFile(filename,ng*12,&Grid[0][0]);
 cout << " start to put daily temperature variation in DailyTempVar[day] " << endl;
 //for (n=0; n<=ng_land-1; n++) {
 for (day = 1; day <= 365; day++) {
 fscanf(file_ptr, "%f", &tempvar);
 dailyTempVar[day] = tempvar;
 }
 //}
 cout << " finished putting daily temperature variation in DailyTempVar[day] " << endl;
 fclose(file_ptr);
 } else {
 cerr << "Unable to open file:" << filename << endl;
 exit(-1);
 }

 }
 */

// added for cell AET calculation (2.1f)
double dailyWaterBalanceClass::getDailyCellLandAET(const int cellNumber) {
	return dailyCellLandAET[cellNumber];
}

double dailyWaterBalanceClass::calc_ext_rad(const short day, const int n) {
	// compute solar declination and sunset hour angle
	// according to equations (1)-(3) in Forsythe et al.
	// Ecological Modelling 80, 1995, pp. 87-95
	// omega_s = D*pi/24, where D is length of day from
	// equation (3), (compare equation (6))
	//
	// the formula for calculation of omega_s from
	// Shuttleworth (1993) does not work for high latitudes:
	// delta = 0.4093 * sin((2*pi*day/365)-1.405);
	// omega_s = acos(-tan(theta)*tan(delta));

	extern geoClass geo;

	// solar declination angle (in radians)
	double declination_angle;

	//gleichung (1) in gleichung (2) eingesetzt Forsythe
	declination_angle = asin(0.39795 * cos(0.2163108 + 2. * atan(0.9671396
			* tan(0.00860 * (day - 186)))));

	// latitude of the site in radians
	const double pi_180 = pi / 180.0;
	double theta;

	theta = -(geo.G_row[n] / 2.0 - 90.25) * pi_180;

	// sunset hour angle (in radians)
	double omega_s; //eigentlich omega_1, so bezeichnet in dis kaspar A.3

	omega_s = (sin(theta) * sin(declination_angle)) / (cos(theta) * cos(
			declination_angle)); //gl(A.8)

	// for daylength coefficients <> 0 this equation has to be
	// modified according to equation (3) from Forsythe et al.
	if ((omega_s) < -1.)
		omega_s = -1;
	if ((omega_s) > 1.)
		omega_s = 1;
	omega_s = pi - acos(omega_s);//omega_s (stundenwinkel) wird nach kaspars konvention aus omega_1 berechnet

	// relative distance earth - sun
	const double pi2_365 = 2. * pi / 365.0;
	double dist_es;

	dist_es = 1. + 0.033 * cos(pi2_365 * day);

	// extraterrestrial radiation [mm/day]
	return (15.392 * dist_es * (omega_s * sin(theta) * sin(declination_angle)
			+ cos(theta) * cos(declination_angle) * sin(omega_s))); //Anhang A.2 Dis Kaspar

}

dailyWaterBalanceClass::~dailyWaterBalanceClass() {
	// delete monthly /yearly grids
	if ((2 == options.grid_store) || (4 == options.grid_store) || (6 == options.grid_store)) {
		delete[]G_monthlyPrecipitation;				G_monthlyPrecipitation 	= NULL;
		delete[]G_monthlyTemperature;				G_monthlyTemperature 	= NULL;
		delete[]G_monthlySunshine;				G_monthlySunshine 	= NULL;    		
		delete[]G_monthlyNetRadiation;	G_monthlyNetRadiation 	= NULL;
    	delete[]G_monthlyNetShortRad; 	G_monthlyNetShortRad 	= NULL; //HMS net_short_wave_rad output
		delete[]G_monthlyNetLongRad; 	G_monthlyNetLongRad 	= NULL; //HMS net_long_wave_rad output		
                delete[]G_monthlyShortDownRad;	G_monthlyShortDownRad 	= NULL;   // HMS radiation output adjustment
		delete[]G_monthlyExtRad;	G_monthlyExtRad 	= NULL;   // HMS radiation output adjustment
                delete[]G_monthlyShortUpRad;	G_monthlyShortUpRad 	= NULL;   // HMS radiation output adjustment
                delete[]G_monthlyLongDownRad;	G_monthlyLongDownRad 	= NULL;   // HMS radiation output adjustment
                delete[]G_monthlyLongUpRad;	G_monthlyLongUpRad 	= NULL;   // HMS radiation output adjustment
		delete[]G_monthlyLAI;			G_monthlyLAI 			= NULL;
		delete[]G_monthlyAlbedo;		G_monthlyAlbedo			= NULL;
		delete[]G_monthlyInterception; G_monthlyInterception 	= NULL;
		delete[]G_monthlycanopyWaterContent;		G_monthlycanopyWaterContent 	= NULL;
		delete[]G_monthlymaxcanopyWaterContent; 	G_monthlymaxcanopyWaterContent 	= NULL;
		delete[]G_monthlyAET;			G_monthlyAET 			= NULL;
		delete[]G_monthlyLandAET;		G_monthlyLandAET        = NULL;
		delete[]G_monthlyLandAET_uncorr;		G_monthlyLandAET_uncorr        = NULL;
		delete[]G_monthlyPET;			G_monthlyPET 			= NULL;
		delete[]G_monthlyOpenWaterPET;	G_monthlyOpenWaterPET 	= NULL;
		delete[]G_monthlyTotalPET;		G_monthlyTotalPET		= NULL;
		delete[]G_monthlyRunoff;		G_monthlyRunoff			= NULL;
		delete[]G_monthlyUrbanRunoff;	G_monthlyUrbanRunoff	= NULL;
		delete[]G_monthlySurfaceRunoff;	 			G_monthlySurfaceRunoff			= NULL;
	//	delete[]G_monthlyGwRunoff;		G_monthlyGwRunoff		= NULL;
		delete[]G_monthlyGwRecharge;	G_monthlyGwRecharge		= NULL;
		//delete[]G_monthlyGwContent;		G_monthlyGwContent		= NULL;
		delete[]G_monthlySnowCoverAvg;	G_monthlySnowCoverAvg	= NULL;
		delete[]G_monthlySWE;			G_monthlySWE			= NULL; 
		delete[]G_monthlySnowFall;		G_monthlySnowFall		= NULL; 
        delete[]G_monthlyRainFall;  G_monthlyRainFall  = NULL; // FP20161011N001
		delete[]G_monthlySnowMelt;		G_monthlySnowMelt		= NULL; 
		delete[]G_monthlySnowEvap;		G_monthlySnowEvap		= NULL; 
		delete[]G_monthlySoilWaterAvg;	G_monthlySoilWaterAvg	= NULL;

		//monthly storage
		delete[] G_monthlyCanopyStorage;			G_monthlyCanopyStorage = NULL;
		delete[] G_monthlySnowStorage;				G_monthlySnowStorage = NULL;
		delete[] G_monthlySoilStorage;				G_monthlySoilStorage = NULL;
		delete[] G_monthly_totalWaterInStorages;	G_monthly_totalWaterInStorages = NULL;
		delete[] G_monthly_totalWaterInStorages_mm;	G_monthly_totalWaterInStorages_mm = NULL;
	}
    // delete monthly grids

        // daily output option 31 start
	if ((3 == options.grid_store) || (4 == options.grid_store)) {
		if (options.outSurfaceRunoffDaily)	{ delete[] G_daily31SurfaceRunoff;			G_daily31SurfaceRunoff = NULL; }
	//	if (options.outGWRunoffDaily)		{ delete[] G_daily31GwRunoff;				G_daily31GwRunoff = NULL; }
		if (options.outGWRechargeDaily)		{ delete[] G_daily31GwRecharge;				G_daily31GwRecharge = NULL; }
       // if (options.outGWStorageDaily || options.outSingleStoragesDaily)		{ delete[] G_daily31GwStorage;				G_daily31GwStorage = NULL; }
		if (options.outSoilWaterDaily)		{ delete[] G_daily31SoilWaterAvg;			G_daily31SoilWaterAvg = NULL; }
		if (options.outLAIDaily)			{ delete[] G_daily31LAI;					G_daily31LAI = NULL; }
		if (options.outLAIDaily)			{ delete[] G_daily31Kc;						G_daily31Kc = NULL; }
		if (options.outAlbedoDaily)			{ delete[] G_daily31Albedo;					G_daily31Albedo = NULL; }
		if (options.outInterceptionDaily)	{ delete[] G_daily31Interception;			G_daily31Interception = NULL; }
		if (options.outCanopyWaterDaily)	{ delete[] G_daily31canopyWaterContent;		G_daily31canopyWaterContent = NULL; }
		if (options.outmaxCanopyWaterDaily)	{ delete[] G_daily31maxcanopyWaterContent;	G_daily31maxcanopyWaterContent = NULL; }
		if (options.outPETDaily)			{ delete[] G_daily31PET;					G_daily31PET = NULL; }
		if (options.outTotalPETDaily)		{ delete[] G_daily31TotalPET;				G_daily31TotalPET = NULL; }
		if (options.outSnowCoverDaily)		{ delete[] G_daily31SnowCoverAvg;			G_daily31SnowCoverAvg = NULL; }
		if (options.outSWEDaily)			{ delete[] G_daily31SWE;					G_daily31SWE = NULL; }
		if (options.outSnowFallDaily)		{ delete[] G_daily31SnowFall;				G_daily31SnowFall = NULL; }
		if (options.outSnowMeltDaily)		{ delete[] G_daily31SnowMelt;				G_daily31SnowMelt = NULL; }
        if (options.outSingleStoragesDaily || options.outTotalWaterInStoragesDaily_km3 || options.outTotalWaterInStoragesDaily_mm)               delete[] G_daily31CanopyStorage;				G_daily31CanopyStorage = NULL;
        if (options.outSingleStoragesDaily || options.outTotalWaterInStoragesDaily_km3 || options.outTotalWaterInStoragesDaily_mm) delete[] G_daily31SnowStorage;				G_daily31SnowStorage = NULL;
        if (options.outSingleStoragesDaily || options.outTotalWaterInStoragesDaily_km3 || options.outTotalWaterInStoragesDaily_mm) delete[] G_daily31SoilStorage;				G_daily31SoilStorage = NULL;
        // HMS radiation output adjustment start
                if(options.outExtRadiationDaily)		delete[] G_daily31ExtRad;					 G_daily31ExtRad			          = NULL;
                if(options.outShortDownRadiationDaily)		delete[] G_daily31ShortDownRad;					 G_daily31ShortDownRad			          = NULL;
                if(options.outShortUpRadiationDaily)		delete[] G_daily31ShortUpRad;					 G_daily31ShortUpRad			          = NULL;
	  	if(options.outNetShortWaveRadiationDaily)		delete[] G_daily31NetShortRad;					 G_daily31NetShortRad			          = NULL;
                if(options.outLongDownRadiationDaily)		delete[] G_daily31LongDownRad;					 G_daily31LongDownRad			          = NULL;
                if(options.outLongUpRadiationDaily)		delete[] G_daily31LongUpRad;					 G_daily31LongUpRad			          = NULL;
	  	if(options.outNetLongWaveRadiationDaily)		delete[] G_daily31NetLongRad;					 G_daily31NetLongRad			          = NULL;
	  	if(options.outNetRadiationDaily)		delete[] G_daily31NetRadiation;					 G_daily31NetRadiation			          = NULL;
	  	// HMS radiation output adjustment end
	}
        // daily output option 31 end

        // daily output option 365 start
	if ((5 == options.grid_store) || (6 == options.grid_store)) {
            if ((options.outSurfaceRunoffDaily)||((options.scoutSurfaceRunoff)&&(2==options.day_store)))	{ delete[] G_daily365SurfaceRunoff;			G_daily365SurfaceRunoff = NULL; }
           //     if ((options.outGWRunoffDaily)||(options.scoutGwRunoff))		{ delete[] G_daily365GwRunoff;				G_daily365GwRunoff = NULL; }
                if ((options.outGWRechargeDaily)||((options.scoutGwRecharge)&&(2==options.day_store)))		{ delete[] G_daily365GwRecharge;			G_daily365GwRecharge = NULL; }
     //   if ((options.outGWStorageDaily || options.outSingleStoragesDaily)||(2 == options.day_store))		{ delete[] G_daily365GwStorage;				G_daily365GwStorage = NULL; }
                if ((options.outSoilWaterDaily)||((options.scoutSoilWater)&&(2==options.day_store)))		{ delete[] G_daily365SoilWaterAvg;			G_daily365SoilWaterAvg = NULL; }
                if ((options.outLAIDaily)||((options.scoutLAI)&&(2==options.day_store)))			{ delete[] G_daily365LAI;					G_daily365LAI = NULL; }
                if ((options.outLAIDaily)||((options.scoutKc)&&(2==options.day_store)))			{ delete[] G_daily365Kc;					G_daily365Kc = NULL; }
                if ((options.outAlbedoDaily)||((options.scoutAlbedo)&&(2==options.day_store)))			{ delete[] G_daily365Albedo;				G_daily365Albedo = NULL; }
                if ((options.outInterceptionDaily)||((options.scoutInterception)&&(2==options.day_store)))	{ delete[] G_daily365Interception;			G_daily365Interception = NULL; }
                if ((options.outCanopyWaterDaily)||((options.scoutCanopyWater)&&(2==options.day_store)))	{ delete[] G_daily365canopyWaterContent;	G_daily365canopyWaterContent = NULL; }
                if ((options.outmaxCanopyWaterDaily)||((options.scoutmaxCanopyWater)&&(2==options.day_store)))	{ delete[] G_daily365maxcanopyWaterContent;	G_daily365maxcanopyWaterContent = NULL; }
                if ((options.outPETDaily)||((options.scoutLandPET)&&(2==options.day_store)))			{ delete[] G_daily365PET;					G_daily365PET = NULL; }
                if ((options.outTotalPETDaily)||((options.scoutCellPET)&&(2==options.day_store)))		{ delete[] G_daily365TotalPET;				G_daily365TotalPET = NULL; }
                if ((options.outSnowCoverDaily)||((options.scoutSnowCov)&&(2==options.day_store)))		{ delete[] G_daily365SnowCoverAvg;			G_daily365SnowCoverAvg = NULL; }
                if ((options.outSWEDaily)||((options.scoutSnowWater)&&(2==options.day_store)))			{ delete[] G_daily365SWE;					G_daily365SWE = NULL; }
                if ((options.outSnowFallDaily)||((options.scoutSnowfall)&&(2==options.day_store)))		{ delete[] G_daily365SnowFall;				G_daily365SnowFall = NULL; }
                if ((options.outSnowMeltDaily)||((options.scoutSnowmelt)&&(2==options.day_store)))		{ delete[] G_daily365SnowMelt;				G_daily365SnowMelt = NULL; }
        if ((options.outSingleStoragesDaily || options.outTotalWaterInStoragesDaily_km3 || options.outTotalWaterInStoragesDaily_mm)||(((options.scoutCanopyStor)||(options.scoutTWSkm3)||(options.scoutTWSmm))&&(2==options.day_store)))                 delete[] G_daily365CanopyStorage;				G_daily365CanopyStorage = NULL;
        if ((options.outSingleStoragesDaily || options.outTotalWaterInStoragesDaily_km3 || options.outTotalWaterInStoragesDaily_mm)||(((options.scoutSnowStor)||(options.scoutTWSkm3)||(options.scoutTWSmm))&&(2==options.day_store)))                 delete[] G_daily365SnowStorage;				G_daily365SnowStorage = NULL;
        if ((options.outSingleStoragesDaily || options.outTotalWaterInStoragesDaily_km3 || options.outTotalWaterInStoragesDaily_mm)||(((options.scoutSoilStor)||(options.scoutTWSkm3)||(options.scoutTWSmm))&&(2==options.day_store)))                 delete[] G_daily365SoilStorage;				G_daily365SoilStorage = NULL;
    	// HMS radiation output adjustment start
                if((options.outExtRadiationDaily)||((options.scoutExtRad)&&(2==options.day_store)))		delete[] G_daily365ExtRad;					 G_daily365ExtRad			          = NULL;
                if((options.outShortDownRadiationDaily)||((options.scoutShortDownRad)&&(2==options.day_store)))		delete[] G_daily365ShortDownRad;					 G_daily365ShortDownRad			          = NULL;
                if((options.outShortUpRadiationDaily)||((options.scoutShortUpRad)&&(2==options.day_store)))		delete[] G_daily365ShortUpRad;					 G_daily365ShortUpRad			          = NULL;
                if((options.outNetShortWaveRadiationDaily)||((options.scoutNetShortRad)&&(2==options.day_store)))		delete[] G_daily365NetShortRad;					 G_daily365NetShortRad			          = NULL;
                if((options.outLongDownRadiationDaily)||((options.scoutLongDownRad)&&(2==options.day_store)))		delete[] G_daily365LongDownRad;					 G_daily365LongDownRad			          = NULL;
                if((options.outLongUpRadiationDaily)||((options.scoutLongUpRad)&&(2==options.day_store)))		delete[] G_daily365LongUpRad;					 G_daily365LongUpRad			          = NULL;
                if((options.outNetLongWaveRadiationDaily)||((options.scoutNetLongRad)&&(2==options.day_store)))		delete[] G_daily365NetLongRad;					 G_daily365NetLongRad			          = NULL;
                if((options.outNetRadiationDaily)||((options.scoutNetRad)&&(2==options.day_store)))		delete[] G_daily365NetRadiation;					 G_daily365NetRadiation			          = NULL;
                // HMS radiation output adjustment end
}
        // daily output option 365 end
	// Superbasin calculation currently obsolete, probably not working correctly in WaterGAP2.2b, do in postprocessing  // FP 2015-10
//	if (options.outDailyValues || options.outDailyInterception	|| options.outSuperbasinClimate) {
//		delete[] dailyPrecSb;		dailyPrecSb = NULL;
//		delete[] dailyEffPrecSb;	dailyEffPrecSb = NULL;
//		delete[] dailySnowSb;		dailySnowSb = NULL;
//		delete[] dailyAET_Sb;		dailyAET_Sb = NULL;
//		delete[] dailyPET_Sb;		dailyPET_Sb = NULL;
//		delete[] dailyRunoffSb;		dailyRunoffSb = NULL;
//		delete[] dailySaturationSb;	dailySaturationSb = NULL;
//		delete[] dailySunshineSb;	dailySunshineSb = NULL;
//		delete[] dailyLaiSb;		dailyLaiSb = NULL;
//		delete[] dailyCanopyWaterSb;dailyCanopyWaterSb = NULL;
//		delete[] dailyCanopyEvapoSb;dailyCanopyEvapoSb = NULL;
//		delete[] dailyTempSb;		dailyTempSb = NULL;
//		delete[] snowDaysSb;		snowDaysSb = NULL;
//		delete[] precAsSnowSb;		precAsSnowSb = NULL;
//	}

	if (options.calc_albedo) {
		//delete[]G_AlbedoSoil;		G_AlbedoSoil = NULL;
		//delete[]G_AlbedoCanopy;		G_AlbedoCanopy = NULL;
		delete[] G_AlbedoCycle;		G_AlbedoCycle = NULL;
	}
}

