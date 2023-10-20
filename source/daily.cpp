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


#include "def.h"
//#include "clcl.h"
//#include "s_max.h"
//#include "snowInElevationFile.h"
//#include "wghmStateFile.h"
#include "globals.h"
#include "routing.h"
//#include "daily.h"
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
//extern optionClass options;
//
//extern geoClass geo; // needed for geo.G_contcell[]

extern geoClass geo; // needed for geo.G_contcell[]

// Cell-specific calibration parameters
#include "calib_param.h"
extern calibParamClass calibParam;


dailyWaterBalanceClass::dailyWaterBalanceClass() {
	LCTdataIsRead = false;

}

template<class T> inline int round(const T value) {
	return (int) floor(value + 0.5);
}

void dailyWaterBalanceClass::calcNewDay(const short day, const short month,const short day_in_month,const short last_day_in_month, const short year, const int n, WghmStateFile &wghmState, AdditionalOutputInputFile &additionalOutIn, SnowInElevationFile &snow_in_elevation, const short readinstatus, calibParamClass &calParam) //OE: calParam added)
{

//    extern laiClass lai;
//    extern climateClass climate;
//    extern climateYearClass climateYear; //for .365 (H08) data input
//    extern geoClass geo;
//    extern landClass land;
//    extern groundwaterFactorClass GW;
//    extern clclClass clcl;
//    extern soilWatCapClass maxSoilWaterCap;
//    extern routingClass routing;
//
//    extern Grid<short> G_toBeCalculated;
//    extern Grid<signed short> G_sbasin; // superbasins
//
//    //------- ARID GROUNDWATER-------
//    extern Grid<short> G_aindex; // arid cells

	//------- ARID GROUNDWATER-------
	extern Grid<short> G_aindex; // arid cells
    extern Grid<char> G_LDD;
	// arid_gw = true  -> calculation for arid areas
	// arid_gw = false -> calculation for humid areas
	bool arid_gw = false;

	//////////////////////////////////// temp & prec /////////////////////////////////////////////////////////////////////////////

	double dailyTempC = 0., dailyShortWave = 0., dailyLongWave = 0., dailySunshine = 0., dailyPrec = 0.;
	double dailyTmin = 0.;
    double dailyWindSpeed = 0., dailyVapPres = 0., dailySnowFall = 0., dailyRainFall = 0., dailySnowCoverFrac = 0., dailyTempRange = 0.;
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
	double lat_heat_md = 0. ;

	// solar radiation
	double solar_rad = 0.;

	// reflected solar radiation from earth (dependend on albedo) // HMS radiation output adjustment
	double refl_rad = 0.;

	// total net radiation [mm/day]
	double net_rad = 0.;

	// different surfaces with different albedo
	double albedo = 0.;

	// daily LAI development
	double dailyLai = 0.;

	// daily Kc development
	double dailyKc = 0.;

    // land area fraction for calculation of waterbalance parameters per grid cell...
    double landAreaFrac = routing.getLandAreaFrac(n);
    double G_landAreaFracPrevTimestep__at__n = 0.;
    // if not first day of model -> prevTimestep
    if (1 == routing.statusStarted_landAreaFracNextTimestep[n])
        G_landAreaFracPrevTimestep__at__n = routing.G_landAreaFracPrevTimestep[n];
    else {  // at first model day there is no prevTimestep
        if (additionalOutIn.additionalfilestatus == 1)  //not first day but restart of model
            G_landAreaFracPrevTimestep__at__n = routing.G_landAreaFracPrevTimestep[n];
        else    //first day
            G_landAreaFracPrevTimestep__at__n = landAreaFrac;
    }

    // variables for fractional routing (soil water content)
	double soil_saturation = 0.;
	double soil_water_overflow = 0.;
	double neg_land_aet = 0.;
    double neg_runoff = 0.;

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
		if (options.time_series == 0) {

			// daily precipitation data
			dailyPrec = climate.G_precipitation_d(n,day_in_month - 1); // if full .31 (WFD) is read in

			// daily temperature [°C]
			dailyTempC = climate.G_temperature_d(n,day_in_month - 1);
			dailyShortWave = climate.G_shortwave_d(n,day_in_month - 1);

			if ((options.cloud == 0) || (options.cloud == 1)) {
				dailyLongWave = climate.G_longwave_d(n,day_in_month- 1);
			}

            // additional information depending on PET-calculation - data are only available as .31 and therefore not transformed to .365 input
    		switch (options.petOpt) {
	    		case 1:	dailyWindSpeed = climate.G_windspeed_d(n,day_in_month - 1)	/ 10.; // [m/s]
		    			dailyVapPres = climate.G_vappres_d(n,day_in_month - 1)	/ 10.; // [kPa]
						dailyTmin = climate.G_tmin_d(n,day_in_month - 1); // [oC]
						break;
				case 2:	dailyWindSpeed = climate.G_windspeed_d(n,day_in_month - 1)	/ 10.; // [m/s]
						dailyVapPres = climate.G_vappres_d(n,day_in_month - 1)	/ 10.; // [kPa]
						dailyTmin = climate.G_tmin_d(n,day_in_month - 1); // [oC]
						break;
				case 3:	dailyTempRange = climate.G_temprange_d(n,day_in_month - 1); // [oC]
						dailyTempRange = abs(dailyTempRange);
						break;
				case 4:
				case 5:	dailyWindSpeed = climate.G_windspeed_d(n,day_in_month - 1) / 10.; // [m/s]
						dailyVapPres = climate.G_vappres_d(n,day_in_month - 1)	/ 10.; // [kPa]
						break;
				case 6:	dailyVapPres = climate.G_vappres_d(n,day_in_month - 1)	/ 10.; //[kPa]
						break;
			    case 7: break;
			}
		}

		// use daily input data per year (.365)
		if (options.time_series == 1) {

			// daily precipitation data
			dailyPrec = climateYear.G_precipitation_d365(n,day - 1); // if full .365 (H08) is read in

			// daily temperature [°C]
            dailyTempC = climateYear.G_temperature_d365(n,day - 1);

			dailyShortWave = climateYear.G_shortwave_d365(n,day - 1);

			if ((options.cloud == 0) || (options.cloud == 1)) {
				dailyLongWave = climateYear.G_longwave_d365(n,day - 1);
			}

		}

		// Cell-specific calibration parameters - Apply multiplier  // FP
		dailyPrec = calParam.getValue(M_PREC,n) * dailyPrec;
		// calculate saturated vapour pressure
		double temp2 = dailyTempC + 237.3;
		double e_s = 0.6108 * exp(17.27 * dailyTempC / temp2);

        // reduced daily temperature for PET-MD
        double annualtempdiff = climate.G_temp_diff[n];
        double dailyTempC_redu = dailyTempC - annualtempdiff;

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

			if (RH < 60) {
				// arid area
				// Cell-specific calibration parameters - Use parameter  // FP
				alpha = calParam.getValue(P_PTC_ARI,n);
				// Cell-specific calibration parameters - Use parameter  // FP
				// maxDailyPET: No distinction anymore arid vs. humid
				maxDailyPET = calParam.getValue(P_PET_MXDY,n);
				a_c = a_c_arid;
				b_c = b_c_arid;
				arid_gw = true;
			} else {
				// humid area
				// Cell-specific calibration parameters - Use parameter  // FP
				alpha = calParam.getValue(P_PTC_HUM,n);
				// Cell-specific calibration parameters - Use parameter  // FP
				// maxDailyPET: No distinction anymore arid vs. humid
				maxDailyPET = calParam.getValue(P_PET_MXDY,n);
				a_c = a_c_humid;
				b_c = b_c_humid;
			}
		}
		// definition with use of Koeppen-climate-classification
		else if (1 == options.clclOpt) {
			if ((clcl.cls[n] == 13) || (clcl.cls[n] == 20) || (clcl.cls[n] == 25)) {
				// arid area
				alpha = clcl.clcl_alpha[n]; //based on Class A Pan evaporation measurements
				// Cell-specific calibration parameters - Use parameter
				// maxDailyPET: No distinction anymore arid vs. humid
				maxDailyPET = calParam.getValue(P_PET_MXDY,n);
				a_c = a_c_arid;
				b_c = b_c_arid;
				arid_gw = true;
			} else {
				// humid area
				// clcl-classes 11+12 should be treated as arid concerning max. daily pot. evap.
				if ((clcl.cls[n] == 11) || (clcl.cls[n] == 12))
					// Cell-specific calibration parameters - Use parameter
					// maxDailyPET: No distinction anymore arid vs. humid
					maxDailyPET = calParam.getValue(P_PET_MXDY,n);
				else
					// Cell-specific calibration parameters - Use parameter  // FP
					// maxDailyPET: No distinction anymore arid vs. humid
					maxDailyPET = calParam.getValue(P_PET_MXDY,n);

				alpha = clcl.clcl_alpha[n];  //based on Class A Pan evaporation measurements
				a_c = a_c_humid;
				b_c = b_c_humid;
			}
		} else {
			// pre-defined arid-humid areas after G_ARID_HUMID.UNF2
			switch(G_aindex[n]){
				case 1: //arid area
						alpha = calParam.getValue(P_PTC_ARI,n);
						maxDailyPET = calParam.getValue(P_PET_MXDY,n);
						a_c = a_c_arid;
						b_c = b_c_arid;
						arid_gw = true;
						break;
				case 0: //humid area
						alpha = calParam.getValue(P_PTC_HUM,n);
						maxDailyPET = calParam.getValue(P_PET_MXDY,n);
						a_c = a_c_humid;
						b_c = b_c_humid;
						break;
				default:
						cerr<<"Error: Invalid value for Arid/humid index: "<< G_aindex[n] <<endl;
						exit(1);
			}
		}
		// LAI -------------------------------------
		// use dailyPrec instead of monthly averaged (daily) input data
		//dailyLai = lai.getDailyLai(n, land.G_landCover[n], arid_gw, dailyTempC, (double) climate.G_precipitation(n,month) / number_of_days_in_month[month], dailyPET);
		// use threshold of precipitation sum instead of dailyPET for decision of growing conditions
		//dailyLai = lai.getDailyLai(n, land.G_landCover[n], arid_gw, dailyTempC, dailyPrec, dailyPET);
		dailyLai = lai.getDailyLai(n, land.G_landCover[n], arid_gw, dailyTempC,	dailyPrec, day, last_day_in_month, additionalOutIn);
		dailyKc = lai.getDailyKc(n, land.G_landCover[n], dailyLai);

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
				albedo = G_AlbedoCycle(n,month);
			} else {
				if (options.use_kc == 1)
					albedo = 0.23; // use gras-reference albedo
				else
					albedo = albedo_lct[land.G_landCover[n] - 1]; // use land cover-specific albedo given in LCT.DAT
			}
		}

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
        double conv_mmd_to_Wm2 = lat_heat / 0.0864; // to convert mm/d back into W/m2 for comparison of modeled radiation outputs
		// calculate solar radiation and afterwards net radiation

		//Observed radiation
		// solar radiation from data set
		solar_rad = conv_Wm2_to_mmd * dailyShortWave; // [mm/day]

		if (options.cloud == 0) {
            // net long wave [mm/day], net long wave input
			net_long_wave_rad = conv_Wm2_to_mmd * dailyLongWave; // [mm/day]
        }
		if (options.cloud == 1) {
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

		//Calculate net longwave radiation
		if (options.cloud == 2) {
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
		net_rad = calParam.getValue(M_NETRAD,n) * (net_short_wave_rad + net_long_wave_rad);

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
		// 7: PT-MD - only for climate change studies
		// begin of calculations for PET
		double dailyPET = 0.; // [mm]
		double dailyOpenWaterPET = 0.; // [mm]


		// psychrometric constant [kPa/oC]
		double atmos_pres = 101.3; // [kPa] atmospheric pressure
		double gamma, gammaStar, gamma_md;
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

			// I believe CRU vapour pressure data not to be reliable -mw-
			// therefore calculation acc to "an update for the definition of reference evapotranpiration" Allen, 1994
			// distinguish arid/humid
			// Cell-specific calibration parameters - Use parameter  // FP
			if (alpha == calParam.getValue(P_PTC_ARI,n)) {
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

			// Cell-specific calibration parameters - Use parameter
			// Humid conditions
			if (alpha == calParam.getValue(P_PTC_HUM,n)) {

				if (dailyTempC <= 0.) {
					//turc nur fuer t>0, bei T< 0 ivanov
					dailyPET = 0.;

				} else {
					if (net_short_wave_rad <= 0.)
						dailyPET = 0.;
					else {
						if (RH <= 50.0) {
							//Fallunterscheidung fuer Turc
							dailyPET = (0.31 * (dailyTempC / (dailyTempC + 15.0)) * (net_short_wave_rad + 2.09) * (1.0 + ((50.001 - (RH)) / 70.0)));
						} else {
							dailyPET = (0.31 * (dailyTempC
									/ (dailyTempC + 15.0))
									* (net_short_wave_rad + 2.09));
							printf("turc rh> 50 \n");
						}

					}
				}
			} else {
				// alpha == alphaArid -> Doorenbos
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

			// Cell-specific calibration parameters - Use parameter
			// Humid conditions
			if (alpha == calParam.getValue(P_PTC_HUM,n)) {
				if (RH <= 50.0){
					//Fallunterscheidung fuer Turc
					dailyOpenWaterPET = (0.31 * (dailyTempC / (dailyTempC + 15.0)) * (openWaterNetShortWaveRad + 2.09) * (1.0 + ((50.0 - (RH)) / 70.0)));
				} else {
					dailyOpenWaterPET = (0.31 * (dailyTempC / (dailyTempC + 15.0)) * (openWaterNetShortWaveRad + 2.09));
				}
			} else {
				// alpha == alphaArid
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

				double d_vp;	// vapour pressure deficit
				d_vp = e_s - dailyVapPres;
				if (d_vp < 0) d_vp = 0.;

				double Wf;
				double awf;
				double bwf;
				short j_day;
				double j_a;
				double j_b;

				// Fallunterscheidung fuer northern and southern hemisphere
				if (-(geo.G_row[n] / 2.0 - 90.25) > 0.0) {
					//Umrechnung von Reihe auf Latitude in Grad; northern hem positiv
					j_day = day;
				} else {
					//southern lat
					if (day >= 182) {
						j_day = (day - 182);
					} else {
						j_day = (day + 183);
					}
				}

				j_a = (-(((j_day - 173.) / 58.) * ((j_day - 173.) / 58.)));
				j_b = (-(((j_day - 243.) / 80.) * ((j_day - 243.) / 80.)));
				awf = (0.4 + 1.4 * exp(j_a));
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


       /**
        * PET-MD : Milly-Dunne scheme is to reflect the effect of atmospheric CO2 concentrations on PET,
        This scheme is only applicable for the climate change studies
        When PET-MD is switched on, until the running year (actual_year) = 2001, watergap will compute daily PET
        using default (Priestley-Taylor (PT)) method. There after the daily Temperature that goes
        into default PT equation is adjusted using a reduction factor (annualtempdiff) and this is called
        reduced temperature (dailyTempC_redu --> temp_redu). annualtempdiff is a reduction factor which we computed outside
        and saved in to G_TEMP_DIFF_[YEAR]*.UNF0 (each file have 67420 values).
        When this PET option is switched on and when actual_year is 2001 or above, Watergap starts reading G_TEMP_DIFF_[YEAR]*.UNF0.
        The reduction factor is the difference between average annual temperature of i-9 to i+10 years (20 year window) and
        baseline average annual temperature (1981-2000) (adapted from Milly and Dunne, 2016).
        To compute PET from open water surfaces we again use default PT method.
        */

        if (7 == options.petOpt) {

            double temp_redu = dailyTempC_redu + 237.3;
            double e_s_md = 0.6108 * exp(17.27 * dailyTempC_redu / temp_redu);
            double inc_svp_md = 4098. * e_s_md / (temp_redu * temp_redu) ;

            //latent heat of evaporation for PET-MD
            if (dailyTempC_redu > 0)
                // latent heat of vaporization of water
                lat_heat_md = 2.501 - 0.002361 * dailyTempC_redu; // [MJ/kg]
            else
                // latent heat of sublimation
                lat_heat_md = 2.835; // 2.501 + 0.334

            gamma_md = c3 / lat_heat_md ;

           if ( year <= 2000) {  // until year 2000 we use default PT method
                if (net_rad <= 0.)
                    dailyPET = 0.;
                else
                    dailyPET = alpha * (inc_svp * net_rad) / (inc_svp + gamma); // [mm/day]

                // calculation of PET for open water surfaces (lakes/wetlands)
                if (openWaterNetRad <= 0.)
                    dailyOpenWaterPET = 0.;
                else
                    dailyOpenWaterPET = alpha * (inc_svp * openWaterNetRad)	/ (inc_svp + gamma);// [mm/day]
           } else {   // after year 2000 we change to milly-dunne but only for land cells, for open water we follow the same method as PT
                if (net_rad <= 0.)
                    dailyPET = 0.;
                else
                    dailyPET = alpha * (inc_svp_md * net_rad) / (inc_svp_md + gamma_md); // [mm/day]
                // calculation of PET for open water surfaces (lakes/wetlands)
                if (openWaterNetRad <= 0.)
                    dailyOpenWaterPET = 0.;
                else
                    dailyOpenWaterPET = alpha * (inc_svp * openWaterNetRad)	/ (inc_svp + gamma);// [mm/day]

           }

        }

        // crop/land-cover specific evapotranspiration using crop coefficients
        if ((G_snow[n] <= 3.) && (options.use_kc == 1)) {
            dailyPET *= dailyKc;
            dailyOpenWaterPET *= kc_OpenWater;
        }

        G_lakeBalance[n] = (dailyPrec - dailyOpenWaterPET) * G_cellCorrFact[n];

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

        if (n < ng) {

			// Efficiency enhancement through local copies instead of often used method getValue
			double P_T_SNOWFZ__at__n = calParam.getValue(P_T_SNOWFZ,n);
			double P_T_SNOWMT__at__n = calParam.getValue(P_T_SNOWMT,n);
			double P_T_GRADNT__at__n = calParam.getValue(P_T_GRADNT,n);
			double M_DEGDAY_F__at__n = calParam.getValue(M_DEGDAY_F,n);

///////////////////////////////////////// interception //////////////////////////////////////////////////////////////////

			if (1 == options.intercept) {
				double canopy_deficiency;
                // adapt canopy water content [mm] to changed land area frac
                // if landAreaFrac == 0; add canopy water content to runoff and set it ot 0
                if (landAreaFrac <= 0.) {

					G_dailyStorageTransfer[n] = G_canopyWaterContent[n];
                    G_canopyWaterContent[n] = 0.;
					dailyCanopyEvapo = 0.;
                }
                else {
					// adapt for changed land area frac and calculate canopy water content
                    G_canopyWaterContent[n] *= G_landAreaFracPrevTimestep__at__n / landAreaFrac;

                    if(abs(G_canopyWaterContent[n])<=routing.minStorVol) //to counter numerical inaccuracies
                        G_canopyWaterContent[n]=0.;

                    initialStorage = G_canopyWaterContent[n]; // calculation of cell AET (2.1f)

                    // calculation of dailyLAI has been moved before calculations of albedo starts
                    if (dailyLai > 0.00001) {
						// Cell-specific calibration parameters - Use parameter
						max_canopy_storage = calParam.getValue(P_MCWH,n) * dailyLai;	// [mm]
                        canopy_deficiency = max_canopy_storage - G_canopyWaterContent[n];
                        if (dailyPrec < canopy_deficiency) {
                            G_canopyWaterContent[n] += dailyPrec;
                            daily_prec_to_soil = 0.;
                        } else {
                            G_canopyWaterContent[n] = max_canopy_storage;
                            daily_prec_to_soil = dailyPrec - canopy_deficiency;
                        }

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
                    // added for cell AET calculation (2.1f)
                    landStorageChangeSum += G_canopyWaterContent[n] - initialStorage;
                }
            }
            else {
				// no interception
				daily_prec_to_soil = dailyPrec;
				dailySoilPET = dailyPET;
				dailyCanopyEvapo = 0.0;
			}

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
            double dailyRainFall_elev = 0.;
			double TempElevMax = 0.;
			// for cell AET calculation
			double snowStorageChange 	= 0.;

			dailySnowFall 	= 0.;
			thresh_elev[n] 	= 0.;
			G_snow[n] 		= 0.;

			// loop for all subgrids in cells with elevation >0m. At the end of the loop, snow cover of
			// all subgrids are added to the 0.5° Grid (G_Snow) again.

			for (short elev = 1; elev < 101; elev++) {

				//count the subgrids
				temp_elev           = 0.;
				dailySnowFall_elev  = 0.;
                dailyRainFall_elev  = 0;
				daily_snow_to_soil_elev = 0.;
				dailyEffPrec_elev   = 0.;
				snowmelt_elev       = 0.;
				dailyEffPrecBeforeSnowMelt_elev = 0.;

				// elevation dependent temperature in one grid cell
				// "0.006" = 0.6°C/100m after Semadeni-Davies (1997) and Dunn, Colohan (1999);
				// The first row of G_ELEV_RANGE.UNF2 (G_Elevation(n,0)) contains mean elevation;
				// Cell-specific calibration parameters - Use parameter

				temp_elev = dailyTempC - ((G_Elevation(n,elev)	- G_Elevation(n,0)) * P_T_GRADNT__at__n);

				// adapt snow depth to changed land area fraction
				// and set snow arrays to 0 in case of land area fraction = 0.
				if (landAreaFrac <= 0.) {
					G_dailyStorageTransfer[n] += G_SnowInElevation(n,elev) / 100.;
					G_SnowInElevation(n,elev) = 0.;
					G_snow[n] = 0.;
				}
				else {
					// landAreaFrac > 0.
                    G_SnowInElevation(n,elev) = G_SnowInElevation(n,elev) * G_landAreaFracPrevTimestep__at__n / landAreaFrac;
					if(abs(G_SnowInElevation(n,elev))<=routing.minStorVol)   //to counter numerical inaccuracies
                        G_SnowInElevation(n,elev) = 0.;

                    // added for cell AET calculation (2.1f)
					initialStorage = G_SnowInElevation(n,elev);

					//-------------------------------------------------------------
					// If snow in subgrid exceeds 1000mm SWE, temperature doesn't decrease in upper subgrids.
					// This should prevent uncontrolled snow accumulation.

					if (G_SnowInElevation(n,elev) > 1000.) {

						// define threshold elevation in cell
						if (thresh_elev[n] == 0.) {
							// remember threshold elevation of cell
							thresh_elev[n] = G_Elevation(n,elev);
						}
						// cell above threshold elevation
						else if (thresh_elev[n] > 0.) {
							// all upper elevations get same temperature calculated with remebered elev
							// Cell-specific calibration parameters - Use parameter
							temp_elev = dailyTempC - ((thresh_elev[n] - G_Elevation(n,0)) * P_T_GRADNT__at__n);
						}

						else {
							cerr << " Negative threshold elevation (thresh_elev[n])  \n";
						}

					}


					// added for cell AET calculation (2.1f)
					// accumulation of new snow
					// Cell-specific calibration parameters - Use parameter
					if (temp_elev <= P_T_SNOWFZ__at__n) {
						dailySnowFall_elev = dailyPrec; // value above canopy, not used, only for output
						dailyEffPrecBeforeSnowMelt_elev = 0.;

						daily_snow_to_soil_elev = daily_prec_to_soil;
						G_SnowInElevation(n,elev) += daily_snow_to_soil_elev; //value below canopy

						if (G_SnowInElevation(n,elev) > dailySoilPET) {
							G_SnowInElevation(n,elev) -= dailySoilPET;
							dailySnowEvapo += dailySoilPET;
						} else {
							dailySnowEvapo += G_SnowInElevation(n,elev);
							G_SnowInElevation(n,elev) = 0.;
						}
					} else {
						dailyEffPrecBeforeSnowMelt_elev = daily_prec_to_soil;
                        dailyRainFall_elev = dailyPrec; // value above canopy, not used, only for output
					}

					// melting of snow
					// Cell-specific calibration parameters - Use parameter
					if (temp_elev > P_T_SNOWMT__at__n) {

						if (G_SnowInElevation(n,elev) < 0.)
							cerr << "G_SnowInElevation(n,elev) < 0 \n";
						else {
							//snowmelt with degree-day factor
							// Cell-specific calibration parameters - Apply multiplier & Use parameter
							snowmelt_elev = M_DEGDAY_F__at__n * ddf_lct[land.G_landCover[n] - 1] * (temp_elev - P_T_SNOWMT__at__n);

							if (snowmelt_elev > G_SnowInElevation(n,elev)) {
								snowmelt_elev = G_SnowInElevation(n,elev);
								G_SnowInElevation(n,elev) = 0.;
							} else {
								G_SnowInElevation(n,elev) -= snowmelt_elev;
							}
						}
					}

					dailyEffPrec_elev = dailyEffPrecBeforeSnowMelt_elev	+ snowmelt_elev;

					// added for cell AET calculation (2.1f)
					snowStorageChange += G_SnowInElevation(n,elev)	- initialStorage;

					// the information about the maximum temperature at the gridcell is used in AET calculations below
					if (elev == 1) {
						TempElevMax = temp_elev;
					}

					if (G_SnowInElevation(n,elev) > 0.)
    					dailySnowCoverFrac += 1. / 100.; // sum up all snow covered subgrids to get the total fraction in cell

                    if (day == last_day_in_month+1) { // AEMS: write out snow in elevation
                        snow_in_elevation.snowInElevation(n, elev) = G_SnowInElevation(n,elev);
                    }
					G_snow[n] += G_SnowInElevation(n,elev); // sum up all snow in subgrids to G_Snow
					dailyEffPrec += dailyEffPrec_elev; // sum up all precipitation in subgrids to grid precip.
					dailySnowFall += dailySnowFall_elev; // sum up all snowfall in subgrids to grid snowfall
                    dailyRainFall += dailyRainFall_elev; // sum up all rainfall in subgrids to grid rainfall
					snowmelt += snowmelt_elev;
					daily_snow_to_soil += daily_snow_to_soil_elev;

				}
			}

			if (landAreaFrac > 0.) {
				G_snow[n] 			/= 100.;

				dailyEffPrec 		/= 100.;	// sum of all subgrids has to be divided by the number of land subgrids
				dailySnowFall 		/= 100.;	// within the cell (if only land-subgrids, value is 100).
				dailyRainFall       /= 100.;
				daily_snow_to_soil 	/= 100.;
				snowmelt 			/= 100.;

				// added by hunger (9/2006): in former snow pack version dailySnowEvapo of subgrids was not summed up
				dailySnowEvapo 		/= 100.;

				// for cell AET calculation
				snowStorageChange 	/= 100.;
				landStorageChangeSum += snowStorageChange;
			}

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

			// for cell AET calculation
			// adapt soil water content [mm] to changed land area frac
			if (landAreaFrac <= 0.) {

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
			else {
                G_soilWaterContent[n] *= G_landAreaFracPrevTimestep__at__n / landAreaFrac;
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

					// daily_runoff is re-calculated later (R=Peff *(S/Smax)^gamma) and cannot be used here.
					// Therefore, soil_water_overflow is increased by the "final" soil_water_overflow in the end.
					// daily_runoff += soil_water_overflow;
				}

				// was "if (dailyTempC > snowFreezeTemp)"
				// Changed to TempElevMax > snowFreezeTemp, because even if dailyTempC < 0, there might be runoff.
				// Now the Temperatur of the cell with lowest elevation is used as threshold.
				// This is caused by the subgrids in the snow algorithm, where some have a higher Temp
				// than the mean temperature of the cell and might produce runoff (here "dailyEffPrec").

				// Cell-specific calibration parameters - Use parameter
				if (TempElevMax > P_T_SNOWFZ__at__n) {

					if (maxSoilWaterCap.G_Smax[n] > 0.) {

						soil_saturation = G_soilWaterContent[n] / maxSoilWaterCap.G_Smax[n];
						daily_runoff = dailyEffPrec * pow(soil_saturation, (double) G_gammaHBV[n]);
						//check wether max daily PET should be limited or not (see maxDailyPET_arid or maxDailyPET_humid)
						if ( dailySoilPET > (maxDailyPET - dailyCanopyEvapo) * soil_saturation )
							dailyAET = (maxDailyPET - dailyCanopyEvapo)	* soil_saturation;
						else
							dailyAET = dailySoilPET;

						G_soilWaterContent[n] += dailyEffPrec - dailyAET - daily_runoff;

						if(abs(G_soilWaterContent[n])<=routing.minStorVol)    //to counter numerical inaccuracies
                            G_soilWaterContent[n]=0.;

						dailyEffPrec = 0.;

						if (G_soilWaterContent[n] < 0.) {

							// too much water has been taken out of the soil water storage:

							// correct dailyAET:
							// reduce it by G_soilWaterContent[n]
							// G_soilWaterContent[n] is negative! THerefore +=
							dailyAET += G_soilWaterContent[n];
							//soil water content is at most as negative as dailyAET!

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
                        if (((arid_gw) && (GW.G_texture[n] < 21)) && (G_LDD[n] >= 0)) {
                            pot_gw_recharge = 0.;
                            if ((GW.getRgmax(n)/100.) < (GW.getgwFactor(n) * daily_runoff))
                                daily_gw_recharge = GW.getRgmax(n)/100.;
                            else
                                daily_gw_recharge = GW.getgwFactor(n) * daily_runoff;

                            if (dailyPrec <= calParam.getValue(P_PCRITGWA,n)) {
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
                    else {
                        // no data for G_Smax (Greenland)
                        total_daily_runoff = 0.;
                        daily_gw_recharge = 0.;
                    }
                }
                    // all precipitation occurs as snow
                else {

                    // Even if (TempElevMax <= snowFreezeTemp), soil_water_overflow occurs as local runoff:
                    soil_water_overflow *= G_cellCorrFact[n];
                    // dailyEffPrec == 0 at this point if (TempElevMax <= snowFreezeTemp). In this case, dailyEffPrec has already been added to soilWaterContent.

                    dailyEffPrec *= G_cellCorrFact[n];
                    total_daily_runoff 	+= soil_water_overflow + dailyEffPrec;
                    daily_gw_recharge 	= 0.;
                    dailyAET = 0.;
                }

                G_dailyGwRecharge[n] = daily_gw_recharge; // to have gw_recharge in routing

				// for cell AET calculation
				landStorageChangeSum += G_soilWaterContent[n] - initialStorage;

				// dailyCellLandAET is a corrected actual total evaporation (canopy, snow and soil)
				// it is consistent with cell-corrected runoff

				dailyCellLandAET[n] = landStorageChangeSum * (G_cellCorrFact[n] - 1.0)
								- dailyPrec * (G_cellCorrFact[n] - 1.0)
								+ (dailyAET + dailyCanopyEvapo + dailySnowEvapo) * G_cellCorrFact[n];

				// to avoid negative aet output
				if (dailyCellLandAET[n] < 0.){
                    neg_land_aet = dailyCellLandAET[n];
                    dailyCellLandAET[n] = 0.;
                }
				dailyCellLandAET_uncorr[n] = (dailyAET + dailyCanopyEvapo + dailySnowEvapo);
			}


///////////////////////////////////////////// surface runoff /////////////////////////////////////////////////////////////
			// If landAreaFrac == 0: total_daily_runoff and daily_gw_recharge == 0.
			if (neg_land_aet < 0.) { // reduce total runoff by negative aet
                total_daily_runoff = total_daily_runoff + neg_land_aet;
			}
            if (total_daily_runoff < 0. ) { // to avoid negative runoff that can occur in some cases due to inaccuracies (>10⁻7)
                total_daily_runoff = 0.;
            }
			if ((total_daily_runoff - daily_gw_recharge) < 0.) {// to avoid negative surface_runoff after reduction of total runoff, all runoff becomes groundwater recharge, thus surface_runoff is 0.
                neg_runoff = total_daily_runoff - daily_gw_recharge;
                daily_gw_recharge = total_daily_runoff;
                G_soilWaterContent[n] += neg_runoff; // in case gw recharge is larger than total runoff, reduce soil storage by the difference}
                neg_runoff = 0.;
            }
			surface_runoff = total_daily_runoff - daily_gw_recharge;
            G_dailyLocalSurfaceRunoff[n] = surface_runoff;
            if(day==last_day_in_month+1)   {
                additionalOutIn.additionalOutputInput(n,20) =  G_soilWaterContent[n]; // is obsolete due to snow canopy and soil being in wghm states
                additionalOutIn.additionalOutputInput(n,22) =  G_canopyWaterContent[n];
                additionalOutIn.additionalOutputInput(n,24) =  G_snow[n];


            }

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

		}

		if ((2 == options.grid_store) || (4 == options.grid_store) || (6 == options.grid_store)) {

			// use the whole amount of prec for all parts off the cell (land/open water)
			G_monthlyPrecipitation(n,month) += dailyPrec;
      		G_monthlyTemperature(n,month) += dailyTempC;
     		G_monthlySunshine(n,month) += dailySunshine / 10;	// sunshine percentage output [in %]
      		G_monthlyExtRad(n,month) += ext_rad * conv_mmd_to_Wm2; // HMS radiation output adjustment
            G_monthlyShortUpRad(n,month) += refl_rad * conv_mmd_to_Wm2; // HMS radiation output adjustment
            G_monthlyLongDownRad(n,month) += long_wave_rad_in * conv_mmd_to_Wm2; // HMS radiation output adjustment
            G_monthlyLongUpRad(n,month) += long_wave_rad_out * conv_mmd_to_Wm2; // HMS radiation output adjustment
            G_monthlyShortDownRad(n,month) += solar_rad * conv_mmd_to_Wm2;	// convert from mm/day to W/m2 // HMS radiation output adjustment
			G_monthlyNetShortRad(n,month) += net_short_wave_rad * conv_mmd_to_Wm2; //HMS net_short_wave_rad output
      		G_monthlyNetLongRad(n,month) += net_long_wave_rad * conv_mmd_to_Wm2; //HMS net_long_wave_rad output
			G_monthlyNetRadiation(n,month) += net_rad * conv_mmd_to_Wm2;	// convert from mm/day to W/m2

			// open water pet, if cell consists solely of open water bodies
			G_monthlyOpenWaterPET(n,month) += dailyOpenWaterPET;
			// pet over land area (whole cell is considered to be land)
			G_monthlyPET(n,month) += dailyPET;

            if (geo.G_contcell[n]) {
                // consider pet and open water pet for land area and water fraction respectively
                G_monthlyTotalPET(n,month) += ( (routing.G_landfreq[n] * dailyPET + routing.G_fwaterfreq[n] * dailyOpenWaterPET)
                                                                               / geo.G_contfreq[n] );
                // grid cell aet; without canopy and snow evap
                G_monthlyAET(n,month) += ( dailyAET * landAreaFrac / geo.G_contfreq[n] );

                // for land AET in [mm]; accounts for land area fraction only (same as for pet and aet)
                G_monthlyLandAET(n,month) +=  ( dailyCellLandAET[n] * landAreaFrac / geo.G_contfreq[n] );
                G_monthlyLandAET_uncorr(n,month) +=  ( dailyCellLandAET_uncorr[n] * landAreaFrac / geo.G_contfreq[n] );

                G_monthlyLAI(n,month) += dailyLai;
                G_monthlyAlbedo(n,month) += (albedo * (landAreaFrac / 100.) + openWaterAlbedo * (1. - (landAreaFrac / 100.)));
                G_monthlyInterception(n,month)	+= (dailyCanopyEvapo * landAreaFrac	/ geo.G_contfreq[n]);

                // store canopy water content at end of month (consistent with water balance calculations)
                // *SE 2012-06* changed to monthly averages
                G_monthlycanopyWaterContent(n,month) += ( G_canopyWaterContent[n] * landAreaFrac / geo.G_contfreq[n] );
                G_monthlymaxcanopyWaterContent(n,month) += ( max_canopy_storage * landAreaFrac / geo.G_contfreq[n] );

                // snow parameters
                G_monthlySnowFall(n,month)	+= (dailySnowFall * landAreaFrac / geo.G_contfreq[n]);
                G_monthlyRainFall(n,month) += (dailyRainFall * landAreaFrac / geo.G_contfreq[n]);
                G_monthlySnowMelt(n,month)	+= (snowmelt * landAreaFrac / geo.G_contfreq[n]);
                G_monthlySnowEvap(n,month)	+= (dailySnowEvapo * landAreaFrac / geo.G_contfreq[n]);

                // store snow water equivalent at end of month (consistent with water balance calculations)
                // *SE 2012-06* changed to monthly averages
                G_monthlySWE(n,month) += ( G_snow[n] * landAreaFrac / geo.G_contfreq[n] );
                G_monthlySnowCoverAvg(n,month) += (dailySnowCoverFrac * landAreaFrac / geo.G_contfreq[n]);

                // soil water content
                // store soil water content at end of month (consistent with water balance calculations)
                // *SE 2012-06* changed to monthly averages
                G_monthlySoilWaterAvg(n,month) += ( G_soilWaterContent[n] * landAreaFrac / geo.G_contfreq[n] );

                // potential runoff from landcells, including the amount of gw-recharge
                G_monthlyRunoff(n,month) += (total_daily_runoff * landAreaFrac / geo.G_contfreq[n]);

                // urban runoff from landcells (accounts only for built-up area)
                G_monthlyUrbanRunoff(n,month) += immediate_runoff;

                // surface runoff from landcells
                G_monthlySurfaceRunoff(n,month) += (surface_runoff * landAreaFrac / geo.G_contfreq[n]);

                // goundwater recharge and gw-runoff into river network
                G_monthlyGwRecharge(n,month) += (daily_gw_recharge * landAreaFrac / geo.G_contfreq[n]);


                // calculate monthly storage per grid cell
                // here: from vertical water balance (first part)
                // in km³
                if (options.grid_store_TypeForStorages == 1) {
                    G_monthlyCanopyStorage(n,month) += G_canopyWaterContent[n] * geo.areaOfCellByArrayPos(n) / 1000000.0 * landAreaFrac / 100.0; // G_canopyWaterContent in mm, G_monthlyCanopyStorage in km³
                    G_monthlySnowStorage(n,month) += G_snow[n]	* geo.areaOfCellByArrayPos(n) / 1000000.0 * landAreaFrac / 100.0; // G_snow in mm, G_monthlySnowStorage in km³
                    G_monthlySoilStorage(n,month) += G_soilWaterContent[n] * geo.areaOfCellByArrayPos(n) / 1000000.0 * landAreaFrac / 100.0; // G_soilWaterContent in mm, G_monthlySoilStorage in km³
                } else if (options.grid_store_TypeForStorages == 0) {
                    G_monthlyCanopyStorage(n,month) = G_canopyWaterContent[n] * geo.areaOfCellByArrayPos(n) / 1000000.0 * landAreaFrac / 100.0; // G_canopyWaterContent in mm, G_monthlyCanopyStorage in km³
                    G_monthlySnowStorage(n,month) = G_snow[n]	* geo.areaOfCellByArrayPos(n) / 1000000.0 * landAreaFrac / 100.0; // G_snow in mm, G_monthlySnowStorage in km³
                    G_monthlySoilStorage(n,month) = G_soilWaterContent[n] * geo.areaOfCellByArrayPos(n) / 1000000.0 * landAreaFrac / 100.0; // G_soilWaterContent in mm, G_monthlySoilStorage in km³
                }

            }
            // 100% water cells (e.g. of Caspian Sea) have now G_fwaterfreq + G_landfreq = 0. The water balance is calculated in the outlet cell of large lakes; therefore a calculation of e.g. precipitation for this cells leads to inconsistent quantification of global fluxes. A division through 0 has to be avoided now (alternative: define it as not land cells in ng).
            else {
                G_monthlyTotalPET(n,month) = 0.;
                G_monthlyAET(n,month) = 0.;
                G_monthlyLandAET(n,month) = 0.;
                G_monthlyLandAET_uncorr(n,month) = 0.;
                G_monthlyInterception(n,month) = 0.;
                G_monthlycanopyWaterContent(n,month) = 0.;
                G_monthlymaxcanopyWaterContent(n,month) = 0.;
                G_monthlySnowFall(n,month) = 0.;
                G_monthlyRainFall(n,month) = 0.;
                G_monthlySnowMelt(n,month) = 0.;
                G_monthlySnowEvap(n,month) = 0.;
                G_monthlySWE(n,month) = 0.;
                G_monthlySnowCoverAvg(n,month) = 0.;
                G_monthlySoilWaterAvg(n,month) = 0.;
                G_monthlyRunoff(n,month) = 0.;
                G_monthlySurfaceRunoff(n,month) = 0.;
                G_monthlyGwRecharge(n,month) = 0.;
            }
        }


        // daily 31 storage calculation
        if ((3 == options.grid_store) || (4 == options.grid_store)) {
            if ((options.outSingleStoragesDaily || options.outTotalWaterInStoragesDaily_km3 || options.outTotalWaterInStoragesDaily_mm)) G_daily31CanopyStorage(n,day_in_month - 1) = G_canopyWaterContent[n] * geo.areaOfCellByArrayPos(n) / 1000000.0 * landAreaFrac / 100.0; // G_canopyWaterContent in mm, G_dailyCanopyStorage in km³
            if ((options.outSingleStoragesDaily || options.outTotalWaterInStoragesDaily_km3 || options.outTotalWaterInStoragesDaily_mm)) G_daily31SnowStorage(n,day_in_month - 1) = G_snow[n]	* geo.areaOfCellByArrayPos(n) / 1000000.0 * landAreaFrac / 100.0; // G_snow in mm, G_dailySnowStorage in km³
            if ((options.outSingleStoragesDaily || options.outTotalWaterInStoragesDaily_km3 || options.outTotalWaterInStoragesDaily_mm)) G_daily31SoilStorage(n,day_in_month - 1) = G_soilWaterContent[n] * geo.areaOfCellByArrayPos(n) / 1000000.0 * landAreaFrac / 100.0; // G_soilWaterContent in mm, G_dailySoilStorage in km³
        }

        if ((5 == options.grid_store) || (6 == options.grid_store)) {
            if ((options.outSingleStoragesDaily || options.outTotalWaterInStoragesDaily_km3 || options.outTotalWaterInStoragesDaily_mm)||(((options.scoutCanopyStor)||(options.scoutTWSkm3)||(options.scoutTWSmm))&&(2==options.day_store))) G_daily365CanopyStorage(n,day - 1) = G_canopyWaterContent[n] * geo.areaOfCellByArrayPos(n) / 1000000.0 * landAreaFrac / 100.0; // G_canopyWaterContent in mm, G_dailyCanopyStorage in km³
            if ((options.outSingleStoragesDaily || options.outTotalWaterInStoragesDaily_km3 || options.outTotalWaterInStoragesDaily_mm)||(((options.scoutSnowStor)||(options.scoutTWSkm3)||(options.scoutTWSmm))&&(2==options.day_store))) G_daily365SnowStorage(n,day - 1) = G_snow[n]	* geo.areaOfCellByArrayPos(n) / 1000000.0 * landAreaFrac / 100.0; // G_snow in mm, G_dailySnowStorage in km³
            if ((options.outSingleStoragesDaily || options.outTotalWaterInStoragesDaily_km3 || options.outTotalWaterInStoragesDaily_mm)||(((options.scoutSoilStor)||(options.scoutTWSkm3)||(options.scoutTWSmm))&&(2==options.day_store))) G_daily365SoilStorage(n,day - 1) = G_soilWaterContent[n] * geo.areaOfCellByArrayPos(n) / 1000000.0 * landAreaFrac / 100.0; // G_soilWaterContent in mm, G_dailySoilStorage in km³
        }


        // daily results
        // daily output option 31 start
		if ((3 == options.grid_store) || (4 == options.grid_store)) {

            // FP 2015-06 // Only continental cells
            // not 100% lake area
            // FP20161010N001 Replacement: geo.G_contfreq[n] instead (geo.G_landfreq[n] + geo.G_fwaterfreq[n])
            if (geo.G_contcell[n]) {

             // precipitation
                if (options.outPrecDaily) G_daily31Precip(n,day_in_month - 1) = dailyPrec;

             	// PET and interception
                if (options.outPETDaily)			G_daily31PET(n,day_in_month - 1) = dailyPET;

                if (options.outTotalPETDaily)		G_daily31TotalPET(n,day_in_month - 1) =
                                                                                        (routing.G_landfreq[n] * dailyPET + routing.G_fwaterfreq[n] * dailyOpenWaterPET)
                                                                                        / geo.G_contfreq[n];
                if (options.outLAIDaily)			G_daily31LAI(n,day_in_month - 1) = dailyLai;
                if (options.outLAIDaily)			G_daily31Kc(n,day_in_month - 1) = dailyKc;
                if (options.outAlbedoDaily)			G_daily31Albedo(n,day_in_month - 1) = albedo * (landAreaFrac / 100.) + openWaterAlbedo * (1. - (landAreaFrac / 100.));
                if (options.outInterceptionDaily)	G_daily31Interception(n,day_in_month - 1) = dailyCanopyEvapo * landAreaFrac / geo.G_contfreq[n];
                if (options.outCanopyWaterDaily)	G_daily31canopyWaterContent(n,day_in_month - 1) = G_canopyWaterContent[n] * landAreaFrac / geo.G_contfreq[n];
                if (options.outmaxCanopyWaterDaily)	G_daily31maxcanopyWaterContent(n,day_in_month - 1) = max_canopy_storage * landAreaFrac / geo.G_contfreq[n];

				// snow parameters
                if (options.outSnowFallDaily)		G_daily31SnowFall(n,day_in_month - 1) = dailySnowFall * landAreaFrac / geo.G_contfreq[n];
                if (options.outSnowMeltDaily)		G_daily31SnowMelt(n,day_in_month - 1) = snowmelt * landAreaFrac / geo.G_contfreq[n];
                if (options.outSWEDaily)			G_daily31SWE(n,day_in_month - 1) = G_snow[n] * landAreaFrac / geo.G_contfreq[n];
                if (options.outSnowCoverDaily)		G_daily31SnowCoverAvg(n,day_in_month - 1) = dailySnowCoverFrac * landAreaFrac / geo.G_contfreq[n];

				// soil water content
                if (options.outSoilWaterDaily)		G_daily31SoilWaterAvg(n,day_in_month - 1) = G_soilWaterContent[n] * landAreaFrac / geo.G_contfreq[n];

				// surface runoff from landcells
                if (options.outSurfaceRunoffDaily)	G_daily31SurfaceRunoff(n,day_in_month - 1) = surface_runoff * landAreaFrac / geo.G_contfreq[n];

			   	// goundwater recharge and gw-runoff into river network
                if (options.outGWRechargeDaily)		G_daily31GwRecharge(n,day_in_month - 1) = daily_gw_recharge * landAreaFrac / geo.G_contfreq[n];

				// radiation output adjustment - daily radiation output
                if(options.outExtRadiationDaily)		G_daily31ExtRad(n,day_in_month - 1) = ext_rad * conv_mmd_to_Wm2;
                if(options.outShortDownRadiationDaily)		G_daily31ShortDownRad(n,day_in_month - 1) = solar_rad * conv_mmd_to_Wm2;
                if(options.outShortUpRadiationDaily)		G_daily31ShortUpRad(n,day_in_month - 1) = refl_rad * conv_mmd_to_Wm2;
                if(options.outNetShortWaveRadiationDaily)		G_daily31NetShortRad(n,day_in_month - 1) = net_short_wave_rad * conv_mmd_to_Wm2;
                if(options.outLongDownRadiationDaily)		G_daily31LongDownRad(n,day_in_month - 1) = long_wave_rad_in * conv_mmd_to_Wm2;
                if(options.outLongUpRadiationDaily)		G_daily31LongUpRad(n,day_in_month - 1) = long_wave_rad_out * conv_mmd_to_Wm2;
                if(options.outNetLongWaveRadiationDaily)		G_daily31NetLongRad(n,day_in_month -  1) = net_long_wave_rad * conv_mmd_to_Wm2;
                if(options.outNetRadiationDaily)		G_daily31NetRadiation(n,day_in_month - 1) = net_rad * conv_mmd_to_Wm2;

            }
            // 100% water cells (e.g. of Caspian Sea) have now G_fwaterfreq + G_landfreq = 0. The water balance is calculated in the outlet cell of large lakes; therefore a calculation of e.g. precipitation for this cells leads to inconsistent quantification of global fluxes. A division through 0 has to be avoided now (alternative: define it as not land cells in ng).
            else {
                if (options.outTotalPETDaily)		G_daily31TotalPET(n,day_in_month - 1) = 0.;
                if (options.outInterceptionDaily)	G_daily31Interception(n,day_in_month - 1) = 0.;
                if (options.outCanopyWaterDaily)	G_daily31canopyWaterContent(n,day_in_month - 1) = 0.;
                if (options.outmaxCanopyWaterDaily)	G_daily31maxcanopyWaterContent(n,day_in_month - 1) = 0.;
                if (options.outSnowFallDaily)		G_daily31SnowFall(n,day_in_month - 1) = 0.;
                if (options.outSnowMeltDaily)		G_daily31SnowMelt(n,day_in_month - 1) = 0.;
                if (options.outSWEDaily)		G_daily31SWE(n,day_in_month - 1) = 0.;
                if (options.outSoilWaterDaily)		G_daily31SoilWaterAvg(n,day_in_month - 1) = 0.;
                if (options.outSurfaceRunoffDaily)	G_daily31SurfaceRunoff(n,day_in_month - 1) = 0.;
                if (options.outGWRechargeDaily)		G_daily31GwRecharge(n,day_in_month - 1) = 0.;
                if (options.outPrecDaily) G_daily31Precip(n,day_in_month - 1) = 0.;
                if (options.outLAIDaily)			G_daily31LAI(n,day_in_month - 1) = 0.;
                if (options.outLAIDaily)			G_daily31Kc(n,day_in_month - 1) = 0.;
                if (options.outAlbedoDaily)			G_daily31Albedo(n,day_in_month - 1) = 0.;
            }

        }

		// daily output option 365
		if ((5 == options.grid_store) || (6 == options.grid_store)) {
			// 100% water cells (e.g. of Caspian Sea) have now G_fwaterfreq + G_landfreq = 0. The water balance is calculated in the outlet cell of large lakes; therefore a calculation of e.g. precipitation for this cells leads to inconsistent quantification of global fluxes. A division through 0 has to be avoided now (alternative: define it as not land cells in ng).
			// not 100% lake area
			if (geo.G_contcell[n]) {
				// precipitation
				if ((options.outPrecDaily)||((options.scoutPrecip)&&(2==options.day_store)))		G_daily365Precip(n,day - 1) = dailyPrec;
                if ((options.outTemp)||((options.scoutPrecip)&&(2==options.day_store)))		G_daily365temp(n,day - 1) = dailyTempC;
                if ((options.outRedTemp)||((options.scoutPrecip)&&(2==options.day_store)&&(7==options.petOpt)))		G_daily365redtemp(n,day - 1) = dailyTempC_redu;


				// PET and interception
				if ((options.outPETDaily)||((options.scoutLandPET)&&(2==options.day_store)))		G_daily365PET(n,day - 1) = dailyPET;
				if ((options.outTotalPETDaily)||((options.scoutCellPET)&&(2==options.day_store)))	G_daily365TotalPET(n,day - 1) = (routing.G_landfreq[n] * dailyPET + routing.G_fwaterfreq[n] * dailyOpenWaterPET) / geo.G_contfreq[n];
				if ((options.outLAIDaily)||((options.scoutLAI)&&(2==options.day_store)))		G_daily365LAI(n,day - 1) = dailyLai;
				if ((options.outLAIDaily)||((options.scoutKc)&&(2==options.day_store)))		G_daily365Kc(n,day - 1) = dailyKc;
				if ((options.outAlbedoDaily)||((options.scoutAlbedo)&&(2==options.day_store)))		G_daily365Albedo(n,day - 1) = albedo * (landAreaFrac / 100.) + openWaterAlbedo * (1. - (landAreaFrac / 100.));
				if ((options.outInterceptionDaily)||((options.scoutInterception)&&(2==options.day_store)))	G_daily365Interception(n,day - 1) = dailyCanopyEvapo * landAreaFrac / geo.G_contfreq[n];
				if ((options.outCanopyWaterDaily)||((options.scoutCanopyWater)&&(2==options.day_store)))	G_daily365canopyWaterContent(n,day - 1) = G_canopyWaterContent[n] * landAreaFrac / geo.G_contfreq[n]; // [mm]
				if ((options.outmaxCanopyWaterDaily)||((options.scoutmaxCanopyWater)&&(2==options.day_store)))	G_daily365maxcanopyWaterContent(n,day - 1)	= max_canopy_storage * landAreaFrac / geo.G_contfreq[n];

				// snow parameters
				if ((options.outSnowFallDaily)||((options.scoutSnowfall)&&(2==options.day_store)))	G_daily365SnowFall(n,day - 1) = dailySnowFall * landAreaFrac / geo.G_contfreq[n];
				if ((options.outSnowMeltDaily)||((options.scoutSnowmelt)&&(2==options.day_store)))	G_daily365SnowMelt(n,day - 1) = snowmelt * landAreaFrac / geo.G_contfreq[n];
				if ((options.outSWEDaily)||((options.scoutSnowWater)&&(2==options.day_store)))		G_daily365SWE(n,day - 1) = G_snow[n] * landAreaFrac / geo.G_contfreq[n]; // [mm]
				if ((options.outSnowCoverDaily)||((options.scoutSnowCov)&&(2==options.day_store)))	G_daily365SnowCoverAvg(n,day - 1) = dailySnowCoverFrac * landAreaFrac / geo.G_contfreq[n];

				// soil water content
				if ((options.outSoilWaterDaily)||((options.scoutSoilWater)&&(2==options.day_store)))	G_daily365SoilWaterAvg(n,day - 1) = G_soilWaterContent[n] * landAreaFrac / geo.G_contfreq[n]; // [mm]

				// surface runoff from landcells
				if ((options.outSurfaceRunoffDaily)||((options.scoutSurfaceRunoff)&&(2==options.day_store)))	G_daily365SurfaceRunoff(n,day - 1) = surface_runoff * landAreaFrac / geo.G_contfreq[n];

				// groundwater recharge and gw-runoff into river network
				if ((options.outGWRechargeDaily)||((options.scoutGwRecharge)&&(2==options.day_store)))	G_daily365GwRecharge(n,day - 1) = daily_gw_recharge * landAreaFrac / geo.G_contfreq[n];

				// radiation output adjustment - daily radiation output
				if((options.outExtRadiationDaily)||((options.scoutExtRad)&&(2==options.day_store)))		G_daily365ExtRad(n,day-1) = ext_rad * conv_mmd_to_Wm2;
				if((options.outShortDownRadiationDaily)||((options.scoutShortDownRad)&&(2==options.day_store)))		G_daily365ShortDownRad(n,day-1) = solar_rad * conv_mmd_to_Wm2;
				if((options.outShortUpRadiationDaily)||((options.scoutShortUpRad)&&(2==options.day_store)))		G_daily365ShortUpRad(n,day-1) = refl_rad * conv_mmd_to_Wm2;
				if((options.outNetShortWaveRadiationDaily)||((options.scoutNetShortRad)&&(2==options.day_store)))		G_daily365NetShortRad(n,day-1) = net_short_wave_rad * conv_mmd_to_Wm2;
				if((options.outLongDownRadiationDaily)||((options.scoutLongDownRad)&&(2==options.day_store)))		G_daily365LongDownRad(n,day-1) = long_wave_rad_in * conv_mmd_to_Wm2;
				if((options.outLongUpRadiationDaily)||((options.scoutLongUpRad)&&(2==options.day_store)))		G_daily365LongUpRad(n,day-1) = long_wave_rad_out * conv_mmd_to_Wm2;
				if((options.outNetLongWaveRadiationDaily)||((options.scoutNetLongRad)&&(2==options.day_store)))		G_daily365NetLongRad(n,day-1) = net_long_wave_rad * conv_mmd_to_Wm2;
				if((options.outNetRadiationDaily)||((options.scoutNetRad)&&(2==options.day_store)))		G_daily365NetRadiation(n,day-1) = net_rad * conv_mmd_to_Wm2;
				if(options.scoutTemp&&(2==options.day_store)) 			G_daily365TempC(n,day-1) = dailyTempC;
			}
			else {
				if (options.outTotalPETDaily || (options.scoutCellPET&&(2==options.day_store)))	G_daily365TotalPET(n,day - 1) = 0.;
				if (options.outLAIDaily || (options.scoutLAI&&(2==options.day_store)))		G_daily365LAI(n,day - 1) = 0.;
				if (options.outLAIDaily || (options.scoutKc&&(2==options.day_store)))		G_daily365Kc(n,day - 1) = 0.;
				if (options.outAlbedoDaily || (options.scoutAlbedo&&(2==options.day_store)))		G_daily365Albedo(n,day - 1) = 0.;
				if (options.outInterceptionDaily || (options.scoutInterception&&(2==options.day_store)))	G_daily365Interception(n,day - 1) = 0.;
				if (options.outCanopyWaterDaily || (options.scoutCanopyWater&&(2==options.day_store)))	G_daily365canopyWaterContent(n,day - 1) = 0.;
				if (options.outmaxCanopyWaterDaily || (options.scoutmaxCanopyWater&&(2==options.day_store)))	G_daily365maxcanopyWaterContent(n,day - 1) = 0.;
				if (options.outSnowFallDaily || (options.scoutSnowfall&&(2==options.day_store)))	G_daily365SnowFall(n,day - 1) = 0.;
				if (options.outSnowMeltDaily || (options.scoutSnowmelt&&(2==options.day_store)))	G_daily365SnowMelt(n,day - 1) = 0.;
				if (options.outSWEDaily || (options.scoutSnowWater&&(2==options.day_store)))		G_daily365SWE(n,day - 1) = 0.;
				if (options.outSnowCoverDaily || (options.scoutSnowCov&&(2==options.day_store)))	G_daily365SnowCoverAvg(n,day - 1) = 0.;
				if (options.outSoilWaterDaily || (options.scoutSoilWater&&(2==options.day_store)))	G_daily365SoilWaterAvg(n,day - 1) = 0.;
				if (options.outSurfaceRunoffDaily || (options.scoutSurfaceRunoff&&(2==options.day_store)))	G_daily365SurfaceRunoff(n,day - 1) = 0.;
				if (options.outGWRechargeDaily || (options.scoutGwRecharge&&(2==options.day_store)))	G_daily365GwRecharge(n,day - 1) = 0.;
				if (options.outPETDaily || (options.scoutLandPET&&(2==options.day_store)))		G_daily365PET(n,day - 1) = 0.;
				if (options.outPrecDaily || (options.scoutPrecip&&(2==options.day_store)))		G_daily365Precip(n,day - 1) = 0.;
                if (options.outTemp || (options.scoutPrecip&&(2==options.day_store)))		G_daily365temp(n,day - 1) = 0.;
                if (options.outRedTemp || (options.scoutPrecip&&(2==options.day_store)&&(7==options.petOpt)))		G_daily365redtemp(n,day - 1) = 0.;
			}
		}
	}
}


void dailyWaterBalanceClass::init(const std::string inputDir, const short nBasins) {

	if ((2 == options.grid_store) || (4 == options.grid_store) || (6 == options.grid_store)) {
		G_monthlyPrecipitation.initialize();
		G_monthlyTemperature.initialize();
		G_monthlyNetRadiation.initialize();
		G_monthlyNetShortRad.initialize(); //HMS net_short_wave_rad output
		G_monthlyNetLongRad.initialize(); //HMS net_long_wave_rad output
		G_monthlyShortDownRad.initialize(); // HMS radiation output adjustment
		G_monthlyExtRad.initialize(); // HMS radiation output adjustment
		G_monthlyShortUpRad.initialize(); // HMS radiation output adjustment
		G_monthlyLongDownRad.initialize(); // HMS radiation output adjustment
		G_monthlyLongUpRad.initialize(); // HMS radiation output adjustment
		G_monthlySunshine.initialize(); // sunshine percentage output
		G_monthlyLAI.initialize();
		G_monthlyAlbedo.initialize();
		G_monthlyInterception.initialize();
		G_monthlycanopyWaterContent.initialize();
		G_monthlymaxcanopyWaterContent.initialize();
		G_monthlyAET.initialize();
		G_monthlyLandAET.initialize(); // land evap only
		G_monthlyLandAET_uncorr.initialize(); // land evap only
		G_monthlyPET.initialize();
		G_monthlyOpenWaterPET.initialize();
		G_monthlyTotalPET.initialize();
		G_monthlyRunoff.initialize();
		G_monthlyUrbanRunoff.initialize();
		G_monthlySurfaceRunoff.initialize();
		G_monthlyGwRecharge.initialize();
		G_monthlySnowCoverAvg.initialize();
		G_monthlySWE.initialize();
		G_monthlySnowFall.initialize();
		G_monthlyRainFall.initialize();
		G_monthlySnowMelt.initialize();
		G_monthlySnowEvap.initialize();
		G_monthlySoilWaterAvg.initialize();

	}
    //OE: should be here? monthly storage components
    G_monthlyCanopyStorage.initialize();
    G_monthlySnowStorage.initialize();
    G_monthlySoilStorage.initialize();
    G_monthly_totalWaterInStorages.initialize();
    G_monthly_totalWaterInStorages_mm.initialize();

    // daily output option 31 start
	if ((3 == options.grid_store) || (4 == options.grid_store)) {

		if (options.outSurfaceRunoffDaily)	G_daily31SurfaceRunoff.initialize();
		if (options.outGWRechargeDaily)		G_daily31GwRecharge.initialize();
		if (options.outSoilWaterDaily)		G_daily31SoilWaterAvg.initialize();
		if (options.outLAIDaily)			G_daily31LAI.initialize();
		if (options.outLAIDaily)			G_daily31Kc.initialize();
		if (options.outAlbedoDaily)			G_daily31Albedo.initialize();
		if (options.outInterceptionDaily)	G_daily31Interception.initialize();
		if (options.outCanopyWaterDaily)	G_daily31canopyWaterContent.initialize();
		if (options.outmaxCanopyWaterDaily)	G_daily31maxcanopyWaterContent.initialize();
		if (options.outPrecDaily)			G_daily31Precip.initialize();
		if (options.outPETDaily)			G_daily31PET.initialize();
		if (options.outTotalPETDaily)		G_daily31TotalPET.initialize();
		if (options.outSnowCoverDaily)		G_daily31SnowCoverAvg.initialize();
		if (options.outSWEDaily)			G_daily31SWE.initialize();
		if (options.outSnowFallDaily)		G_daily31SnowFall.initialize();
		if (options.outSnowMeltDaily)		G_daily31SnowMelt.initialize();

		if (options.outSingleStoragesDaily || options.outTotalWaterInStoragesDaily_km3 || options.outTotalWaterInStoragesDaily_mm)                G_daily31CanopyStorage.initialize();
        if (options.outSingleStoragesDaily || options.outTotalWaterInStoragesDaily_km3 || options.outTotalWaterInStoragesDaily_mm)                G_daily31SnowStorage.initialize();
        if (options.outSingleStoragesDaily || options.outTotalWaterInStoragesDaily_km3 || options.outTotalWaterInStoragesDaily_mm)                G_daily31SoilStorage.initialize();

	    // daily radiation output
        if(options.outExtRadiationDaily)			G_daily31ExtRad.initialize();
        if(options.outShortDownRadiationDaily)		G_daily31ShortDownRad.initialize();
        if(options.outShortUpRadiationDaily)		G_daily31ShortUpRad.initialize();
        if(options.outNetShortWaveRadiationDaily)	G_daily31NetShortRad.initialize();
        if(options.outLongDownRadiationDaily)		G_daily31LongDownRad.initialize();
        if(options.outLongUpRadiationDaily)			G_daily31LongUpRad.initialize();
        if(options.outNetLongWaveRadiationDaily)	G_daily31NetLongRad.initialize();
        if(options.outNetRadiationDaily)			G_daily31NetRadiation.initialize();
    }


	// daily output option 365 start
	if ((5 == options.grid_store) || (6 == options.grid_store)) {
		if ((options.outSurfaceRunoffDaily)||((options.scoutSurfaceRunoff)&&(2==options.day_store)))	G_daily365SurfaceRunoff.initialize();
		if ((options.outGWRechargeDaily)||((options.scoutGwRecharge)&&(2==options.day_store)))			G_daily365GwRecharge.initialize();
		if ((options.outSoilWaterDaily)||((options.scoutSoilWater)&&(2==options.day_store)))			G_daily365SoilWaterAvg.initialize();
		if ((options.outLAIDaily)||((options.scoutLAI)&&(2==options.day_store)))						G_daily365LAI.initialize();
		if ((options.outLAIDaily)||((options.scoutKc)&&(2==options.day_store)))							G_daily365Kc.initialize();
		if ((options.outAlbedoDaily)||((options.scoutAlbedo)&&(2==options.day_store)))					G_daily365Albedo.initialize();
		if ((options.outInterceptionDaily)||((options.scoutInterception)&&(2==options.day_store)))		G_daily365Interception.initialize();
		if ((options.outCanopyWaterDaily)||((options.scoutCanopyWater)&&(2==options.day_store)))		G_daily365canopyWaterContent.initialize();
		if ((options.outmaxCanopyWaterDaily)||((options.scoutmaxCanopyWater)&&(2==options.day_store)))	G_daily365maxcanopyWaterContent.initialize();
        if ((options.outPrecDaily)||((options.scoutPrecip)&&(2==options.day_store)))					G_daily365Precip.initialize();
        if ((options.outTemp)||((options.scoutPrecip)&&(2==options.day_store)))					        G_daily365temp.initialize();
        if ((options.outRedTemp)||((options.scoutPrecip)&&(2==options.day_store)&&(7==options.petOpt)))	 G_daily365redtemp.initialize();
		if ((options.outPETDaily)||((options.scoutLandPET)&&(2==options.day_store)))					G_daily365PET.initialize();
		if ((options.outTotalPETDaily)||((options.scoutCellPET)&&(2==options.day_store)))				G_daily365TotalPET.initialize();
		if ((options.outSnowCoverDaily)||((options.scoutSnowCov)&&(2==options.day_store)))				G_daily365SnowCoverAvg.initialize();
		if ((options.outSWEDaily)||((options.scoutSnowWater)&&(2==options.day_store)))					G_daily365SWE.initialize();
		if ((options.outSnowFallDaily)||((options.scoutSnowfall)&&(2==options.day_store)))				G_daily365SnowFall.initialize();
		if ((options.outSnowMeltDaily)||((options.scoutSnowmelt)&&(2==options.day_store)))				G_daily365SnowMelt.initialize();
        if ((options.outSingleStoragesDaily || options.outTotalWaterInStoragesDaily_km3 || options.outTotalWaterInStoragesDaily_mm)||(((options.scoutCanopyStor)||(options.scoutTWSkm3)||(options.scoutTWSmm))&&(2==options.day_store)))	G_daily365CanopyStorage.initialize();
        if ((options.outSingleStoragesDaily || options.outTotalWaterInStoragesDaily_km3 || options.outTotalWaterInStoragesDaily_mm)||(((options.scoutSnowStor)||(options.scoutTWSkm3)||(options.scoutTWSmm))&&(2==options.day_store)))      G_daily365SnowStorage.initialize();
        if ((options.outSingleStoragesDaily || options.outTotalWaterInStoragesDaily_km3 || options.outTotalWaterInStoragesDaily_mm)||(((options.scoutSoilStor)||(options.scoutTWSkm3)||(options.scoutTWSmm))&&(2==options.day_store)))      G_daily365SoilStorage.initialize();

		// daily radiation output
        if((options.outExtRadiationDaily)||((options.scoutExtRad)&&(2==options.day_store)))					G_daily365ExtRad.initialize();
        if((options.outShortDownRadiationDaily)||((options.scoutShortDownRad)&&(2==options.day_store)))		G_daily365ShortDownRad.initialize();
        if((options.outShortUpRadiationDaily)||((options.scoutShortUpRad)&&(2==options.day_store)))			G_daily365ShortUpRad.initialize();
        if((options.outNetShortWaveRadiationDaily)||((options.scoutNetShortRad)&&(2==options.day_store)))	G_daily365NetShortRad.initialize();
        if((options.outLongDownRadiationDaily)||((options.scoutLongDownRad)&&(2==options.day_store)))		G_daily365LongDownRad.initialize();
        if((options.outLongUpRadiationDaily)||((options.scoutLongUpRad)&&(2==options.day_store)))			G_daily365LongUpRad.initialize();
        if((options.outNetLongWaveRadiationDaily)||((options.scoutNetLongRad)&&(2==options.day_store)))		G_daily365NetLongRad.initialize();
        if((options.outNetRadiationDaily)||((options.scoutNetRad)&&(2==options.day_store)))					G_daily365NetRadiation.initialize();

        if ((options.scoutTemp&&(2==options.day_store)))		G_daily365TempC.initialize();
	}

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

	// albedo
	if (options.calc_albedo) {
		G_AlbedoCycle.read(options.input_dir + "/G_ALBEDO_CYCLE.12.UNF0");
	}
}

void dailyWaterBalanceClass::annualInit() {
	if ((2 == options.grid_store) || (4 == options.grid_store) || (6 == options.grid_store)) {
		for (int n = 0; n < ng; ++n) {
			for (short month = 0; month < 12; month++) {
				G_monthlyPrecipitation(n,month) = 0.;
				G_monthlyTemperature(n,month) = 0.;
				G_monthlySunshine(n,month) = 0.;   // sunshine percentage output
				G_monthlyNetRadiation(n,month) = 0.;
				G_monthlyNetShortRad(n,month) = 0.; // HMS net_short_wave_rad output
        		G_monthlyNetLongRad(n,month) = 0.; // HMS net_long_wave_rad output
                G_monthlyShortDownRad(n,month) = 0.; // HMS radiation output adjustment
				G_monthlyExtRad(n,month) = 0.; // HMS radiation output adjustment
                G_monthlyShortUpRad(n,month) = 0.; // HMS radiation output adjustment
				G_monthlyLongDownRad(n,month) = 0.; // HMS radiation output adjustment
				G_monthlyLongUpRad(n,month) = 0.; // HMS radiation output adjustment
				G_monthlyLAI(n,month) = 0.;
				G_monthlyAlbedo(n,month) = 0.;
				G_monthlyInterception(n,month) = 0.;
				G_monthlycanopyWaterContent(n,month) = 0.;
				G_monthlymaxcanopyWaterContent(n,month) = 0.;
				G_monthlyAET(n,month) = 0.;
				G_monthlyLandAET(n,month) = 0.;
				G_monthlyLandAET_uncorr(n,month) = 0.;
				G_monthlyPET(n,month) = 0.;
				G_monthlyOpenWaterPET(n,month) = 0.;
				G_monthlyTotalPET(n,month) = 0.;
				G_monthlyRunoff(n,month) = 0.;
				G_monthlyUrbanRunoff(n,month) = 0.;
				G_monthlySurfaceRunoff(n,month) = 0.;
				G_monthlyGwRecharge(n,month) = 0.;
				G_monthlySnowCoverAvg(n,month) = 0.;
				G_monthlySWE(n,month) = 0.;
				G_monthlySnowMelt(n,month) = 0.;
				G_monthlySnowEvap(n,month) = 0.;
				G_monthlySnowFall(n,month) = 0.;
                G_monthlyRainFall(n,month) = 0.;
				G_monthlySoilWaterAvg(n,month) = 0.;
				G_monthlyCanopyStorage(n,month) = 0.;
				G_monthlySnowStorage(n,month) = 0.;
				G_monthlySoilStorage(n,month) = 0.;
				G_monthly_totalWaterInStorages(n,month) = 0.;
				G_monthly_totalWaterInStorages_mm(n,month) = 0.;

			}
		}
	}
}

void dailyWaterBalanceClass::daily31outInit() {
	short d;
	int n;
	for (n = 0; n < ng; n++)
		for (d = 0; d < 31; d++) {
			if (options.outSurfaceRunoffDaily)
				G_daily31SurfaceRunoff(n,d) = 0.;
			if (options.outGWRechargeDaily)
				G_daily31GwRecharge(n,d) = 0.;
			if (options.outSoilWaterDaily)
				G_daily31SoilWaterAvg(n,d) = 0.;
			if (options.outLAIDaily)
				G_daily31LAI(n,d) = 0.;
			if (options.outLAIDaily)
				G_daily31Kc(n,d) = 0.;
			if (options.outAlbedoDaily)
				G_daily31Albedo(n,d) = 0.;
			if (options.outInterceptionDaily)
				G_daily31Interception(n,d) = 0.;
			if (options.outCanopyWaterDaily)
				G_daily31canopyWaterContent(n,d) = 0.;
			if (options.outmaxCanopyWaterDaily)
				G_daily31maxcanopyWaterContent(n,d) = 0.;
			if (options.outPrecDaily)
				G_daily31Precip(n,d) = 0.;
			if (options.outPETDaily)
				G_daily31PET(n,d) = 0.;
			if (options.outTotalPETDaily)
				G_daily31TotalPET(n,d) = 0.;
			if (options.outSnowCoverDaily)
				G_daily31SnowCoverAvg(n,d) = 0.;
			if (options.outSWEDaily)
				G_daily31SWE(n,d) = 0.;
			if (options.outSnowFallDaily)
				G_daily31SnowFall(n,d) = 0.;
			if (options.outSnowMeltDaily)
				G_daily31SnowMelt(n,d) = 0.;
			if (options.outSingleStoragesDaily || options.outTotalWaterInStoragesDaily_km3 || options.outTotalWaterInStoragesDaily_mm)
				G_daily31CanopyStorage(n,d) = 0.;
			if (options.outSingleStoragesDaily || options.outTotalWaterInStoragesDaily_km3 || options.outTotalWaterInStoragesDaily_mm)
				G_daily31SnowStorage(n,d) = 0.;
			if (options.outSingleStoragesDaily || options.outTotalWaterInStoragesDaily_km3 || options.outTotalWaterInStoragesDaily_mm)
				G_daily31SoilStorage(n,d) = 0.;

			// radiation output adjustment start - daily radiation output
			if(options.outExtRadiationDaily)
				G_daily31ExtRad(n,d)		= 0.;
			if(options.outShortDownRadiationDaily)
				G_daily31ShortDownRad(n,d)		= 0.;
			if(options.outShortUpRadiationDaily)
				G_daily31ShortUpRad(n,d)		= 0.;
			if(options.outNetShortWaveRadiationDaily)
				G_daily31NetShortRad(n,d)		= 0.;
			if(options.outLongDownRadiationDaily)
				G_daily31LongDownRad(n,d)		= 0.;
			if(options.outLongUpRadiationDaily)
				G_daily31LongUpRad(n,d)		= 0.;
			if(options.outNetLongWaveRadiationDaily)
				G_daily31NetLongRad(n,d)		= 0.;
			if(options.outNetRadiationDaily)
				G_daily31NetRadiation(n,d)		= 0.;
	}
}

void dailyWaterBalanceClass::daily365outInit() {

	short d;
	for (int n = 0; n < ng; n++)
		for (d = 0; d < 365; d++) {
			if ((options.outSurfaceRunoffDaily)||((options.scoutSurfaceRunoff)&&(2==options.day_store)))
				G_daily365SurfaceRunoff(n,d) = 0.;
			if ((options.outGWRechargeDaily)||((options.scoutGwRecharge)&&(2==options.day_store)))
				G_daily365GwRecharge(n,d) = 0.;
			if ((options.outSoilWaterDaily)||((options.scoutSoilWater)&&(2==options.day_store)))
				G_daily365SoilWaterAvg(n,d) = 0.;
			if ((options.outLAIDaily)||((options.scoutLAI)&&(2==options.day_store)))
				G_daily365LAI(n,d) = 0.;
			if ((options.outLAIDaily)||((options.scoutKc)&&(2==options.day_store)))
				G_daily365Kc(n,d) = 0.;
			if ((options.outAlbedoDaily)||((options.scoutAlbedo)&&(2==options.day_store)))
				G_daily365Albedo(n,d) = 0.;
			if ((options.outInterceptionDaily)||((options.scoutInterception)&&(2==options.day_store)))
				G_daily365Interception(n,d) = 0.;
			if ((options.outCanopyWaterDaily)||((options.scoutCanopyWater)&&(2==options.day_store)))
				G_daily365canopyWaterContent(n,d) = 0.;
			if ((options.outmaxCanopyWaterDaily)||((options.scoutmaxCanopyWater)&&(2==options.day_store)))
				G_daily365maxcanopyWaterContent(n,d) = 0.;
			if ((options.outPrecDaily)||((options.scoutPrecip)&&(2==options.day_store)))
				G_daily365Precip(n,d) = 0.;
            if ((options.outTemp)||((options.scoutPrecip)&&(2==options.day_store)))
                G_daily365temp(n,d) = 0.;
            if ((options.outRedTemp)||((options.scoutPrecip)&&(2==options.day_store)&&(7==options.petOpt)))
                G_daily365redtemp(n,d) = 0.;
			if ((options.outPETDaily)||((options.scoutLandPET)&&(2==options.day_store)))
				G_daily365PET(n,d) = 0.;
			if ((options.outTotalPETDaily)||((options.scoutCellPET)&&(2==options.day_store)))
				G_daily365TotalPET(n,d) = 0.;
			if ((options.outSnowCoverDaily)||((options.scoutSnowCov)&&(2==options.day_store)))
				G_daily365SnowCoverAvg(n,d) = 0.;
			if ((options.outSWEDaily)||((options.scoutSnowWater)&&(2==options.day_store)))
				G_daily365SWE(n,d) = 0.;
			if ((options.outSnowFallDaily)||((options.scoutSnowfall)&&(2==options.day_store)))
				G_daily365SnowFall(n,d) = 0.;
			if ((options.outSnowMeltDaily)||((options.scoutSnowmelt)&&(2==options.day_store)))
				G_daily365SnowMelt(n,d) = 0.;
			if ((options.outSingleStoragesDaily || options.outTotalWaterInStoragesDaily_km3 || options.outTotalWaterInStoragesDaily_mm)||(((options.scoutCanopyStor)||(options.scoutTWSkm3)||(options.scoutTWSmm))&&(2==options.day_store)))
				G_daily365CanopyStorage(n,d) = 0.;
			if ((options.outSingleStoragesDaily || options.outTotalWaterInStoragesDaily_km3 || options.outTotalWaterInStoragesDaily_mm)||(((options.scoutSnowStor)||(options.scoutTWSkm3)||(options.scoutTWSmm))&&(2==options.day_store)))
				G_daily365SnowStorage(n,d) = 0.;
			if ((options.outSingleStoragesDaily || options.outTotalWaterInStoragesDaily_km3 || options.outTotalWaterInStoragesDaily_mm)||(((options.scoutSoilStor)||(options.scoutTWSkm3)||(options.scoutTWSmm))&&(2==options.day_store)))
				G_daily365SoilStorage(n,d) = 0.;

			// radiation output adjustment - daily radiation output
			if((options.outExtRadiationDaily)||((options.scoutExtRad)&&(2==options.day_store)))
				G_daily365ExtRad(n,d)		= 0.;
			if((options.outShortDownRadiationDaily)||((options.scoutShortDownRad)&&(2==options.day_store)))
				G_daily365ShortDownRad(n,d)		= 0.;
			if((options.outShortUpRadiationDaily)||((options.scoutShortUpRad)&&(2==options.day_store)))
				G_daily365ShortUpRad(n,d)		= 0.;
			if((options.outNetShortWaveRadiationDaily)||((options.scoutNetShortRad)&&(2==options.day_store)))
				G_daily365NetShortRad(n,d)		= 0.;
			if((options.outLongDownRadiationDaily)||((options.scoutLongDownRad)&&(2==options.day_store)))
				G_daily365LongDownRad(n,d)		= 0.;
			if((options.outLongUpRadiationDaily)||((options.scoutLongUpRad)&&(2==options.day_store)))
				G_daily365LongUpRad(n,d)		= 0.;
			if((options.outNetLongWaveRadiationDaily)||((options.scoutNetLongRad)&&(2==options.day_store)))
				G_daily365NetLongRad(n,d)		= 0.;
			if((options.outNetRadiationDaily)||((options.scoutNetRad)&&(2==options.day_store)))
				G_daily365NetRadiation(n,d)		= 0.;
			if ((options.scoutTemp&&(2==options.day_store)))
				G_daily365TempC(n,d)		= 0.;
		}
}

void dailyWaterBalanceClass::setParameters() {
	canopyEvapoExp      = 0.66666666;
	runoffFracBuiltUp   = 0.5;

	openWaterAlbedo     = 0.08;
	kc_OpenWater        = 1.05;

	a_s = 0.25; // a_s: fraction of extraterrestrial radiation on overcast days
	b_s = 0.5; // a_s + b_s: fraction of extraterrestrial radiation on clear days

	a_c_arid = 1.35; // int-wave radiation coefficients for clear skies
	b_c_arid = -0.35; // sum has to be 1.0
	a_c_humid = 1.00;
	b_c_humid = 0.00;
}

void dailyWaterBalanceClass::setStoragesToZero() {
	for (int n = 0; n <= ng - 1; n++) {
		G_soilWaterContent[n]	= 0.;
		G_canopyWaterContent[n]	= 0.;
		G_snow[n]				= 0.;
		G_dailyStorageTransfer[n] = 0.;
	}
    for (int n = 0; n < ng; n++) {
		for (short elev = 0; elev < 101; elev++) {
			G_SnowInElevation(n,elev) = 0.;
		}
	}
}
// AEMS: new function to set storages:
void dailyWaterBalanceClass::setStorages(WghmStateFile &wghmState, SnowInElevationFile &snow_in_elevation, AdditionalOutputInputFile &additionalOutIn){
//    setStoragesToZero();
//    extern routingClass routing;
//    extern geoClass geo; // AEMS: added

    // AEMS: for all grid cells:
    for(int n=0;n<ng;n++)
    {
        // CH_laf:
        double landAreaFrac = routing.getLandAreaFrac(n);

        if(landAreaFrac <= 0.) {
            G_canopyWaterContent[n] = 0.;
            G_snow[n] = 0.;
            G_soilWaterContent[n] = 0.;
            for (short elev = 0; elev <= 100; elev++) {
                G_SnowInElevation(n, elev) = 0.;
            }
        }
        else {
            G_canopyWaterContent[n] = wghmState.cell(n).canopy(0) * geo.G_contfreq[n] / landAreaFrac;
            G_snow[n] = wghmState.cell(n).snow(0) * geo.G_contfreq[n] / landAreaFrac;
            G_soilWaterContent[n] = wghmState.cell(n).soil(0) * geo.G_contfreq[n] / landAreaFrac;
            for (short elev = 0; elev <= 100; elev++) {
                G_SnowInElevation(n, elev) = snow_in_elevation.snowInElevation(n,elev)*geo.G_contfreq[n]/landAreaFrac;
            }
        }        // AEMS: other storages will be set in routing.setStorages()
    }
}
void dailyWaterBalanceClass::readLCTdata(const std::string inputDir) {

	std::string filename = inputDir + "/LCT_22.DAT";

	// open LCT_22.DAT
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

void dailyWaterBalanceClass::writeAnnualGrids(const string output_dir, const short year) {

	// everything is stored as monthly grids:
	// nothing to do....

	// The routine is kept. Perhaps we want to store some of the 'storages'.
}

void dailyWaterBalanceClass::writeMonthlyGrids(const string outDir, const short year, const short month2store) {

	if (2 != options.grid_store && 4 != options.grid_store && 6	!= options.grid_store)
		return;

	extern routingClass routing;
	extern geoClass geo;

	char filename[250];
	double cellArea;
    double storages_routing = 0.;
    int n, month;

	if (options.outPrec) {
		G_monthlyPrecipitation.write(outDir + "/G_PRECIPITATION_" + to_string(year) + ".12.UNF0");
	}
	if (options.outInterception) {
		G_monthlyInterception.write(outDir + "/G_INTERCEPTION_" + to_string(year) + ".12.UNF0");
	}
	if (options.outPET) {
		G_monthlyPET.write(outDir + "/G_PET_" + to_string(year) + ".12.UNF0");
	}
	if (options.outTotalPET) {
		G_monthlyTotalPET.write(outDir + "/G_TOTAL_PET_" + to_string(year) + ".12.UNF0");
	}
	if (options.outAET) {
		G_monthlyAET.write(outDir + "/G_AET_" + to_string(year) + ".12.UNF0");
	}
	if (options.outCellAET) {
		G_monthlyLandAET.write(outDir + "/G_LAND_AET_" + to_string(year) + ".12.UNF0");
		G_monthlyLandAET_uncorr.write(outDir + "/G_LAND_AET_UNCORR_" + to_string(year) + ".12.UNF0");
	}
	if (options.outRunoff) {
		G_monthlyRunoff.write(outDir + "/G_RUNOFF_" + to_string(year) + ".12.UNF0");
	}
	if (options.outUrbanRunoff) {
		G_monthlyUrbanRunoff.write(outDir + "/G_URBAN_RUNOFF_" + to_string(year) + ".12.UNF0");
	}
	if (options.outSurfaceRunoff) {
		G_monthlySurfaceRunoff.write(outDir + "/G_SURFACE_RUNOFF_" + to_string(year) + ".12.UNF0");
	}
    if (options.outGWRecharge) {
		G_monthlyGwRecharge.write(outDir + "/G_GW_RECHARGE_mm_" + to_string(year) + ".12.UNF0");

        for (n = 0; n < ng; n++)
           for (short m = 0; m < 12; m++)
              G_monthlyArray(n,m) = G_monthlyGwRecharge(n,m)
                 * (geo.areaOfCellByArrayPos(n) * (geo.G_contfreq[n] / 100.)) / 1000000.0;
		G_monthlyArray.write(outDir + "/G_GW_RECHARGE_km3_" + to_string(year) + ".12.UNF0");
	}
	if (options.outOpenWaterPET) {
		G_monthlyOpenWaterPET.write(outDir + "/G_OPEN_WATER_PET_" + to_string(year) + ".12.UNF0");
	}

	int numberOfDaysInMonth[12] = { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };
	for (n = 0; n <= ng - 1; n++) {
		for (month = 0; month < 12; month++) {
			G_monthlySWE(n,month) /= numberOfDaysInMonth[month];
			G_monthlySnowCoverAvg(n,month) /= numberOfDaysInMonth[month];
			G_monthlySoilWaterAvg(n,month) /= numberOfDaysInMonth[month];
			G_monthlyNetRadiation(n,month) /= numberOfDaysInMonth[month]; // HMS radiation output adjustment
            G_monthlyShortDownRad(n,month) /= numberOfDaysInMonth[month];  // HMS radiation output adjustment
			G_monthlyExtRad(n,month) /= numberOfDaysInMonth[month];  // HMS radiation output adjustment
            G_monthlyShortUpRad(n,month) /= numberOfDaysInMonth[month];  // HMS radiation output adjustment
			G_monthlyNetShortRad(n,month) /= numberOfDaysInMonth[month];  // HMS radiation output adjustment
            G_monthlyLongDownRad(n,month) /= numberOfDaysInMonth[month];  // HMS radiation output adjustment
            G_monthlyLongUpRad(n,month) /= numberOfDaysInMonth[month];  // HMS radiation output adjustment
			G_monthlyNetLongRad(n,month) /= numberOfDaysInMonth[month];  // HMS radiation output adjustment
            G_monthlyTemperature(n,month) /= numberOfDaysInMonth[month];  // monthly temperature output
			G_monthlySunshine(n,month) /= numberOfDaysInMonth[month];  // monthly sunshine output
			G_monthlyLAI(n,month) /= numberOfDaysInMonth[month];
			G_monthlyAlbedo(n,month) /= numberOfDaysInMonth[month];
			G_monthlycanopyWaterContent(n,month) /= numberOfDaysInMonth[month];
			G_monthlymaxcanopyWaterContent(n,month) /= numberOfDaysInMonth[month];

            // GRACE application - divide by days in month to get monthly mean values
            if (options.grid_store_TypeForStorages == 1) {
				G_monthlyCanopyStorage(n,month) /= numberOfDaysInMonth[month];
				G_monthlySnowStorage(n,month) /= numberOfDaysInMonth[month];
				G_monthlySoilStorage(n,month) /= numberOfDaysInMonth[month];
			}

			// Storage summation restructured
			storages_routing = routing.G_monthlyLocLakeStorage(n,month) // all values are monthly means [km3]
				+ routing.G_monthlyLocWetlStorage(n,month)
				+ routing.G_monthlyGloLakeStorage(n,month)
				+ routing.G_monthlyGloWetlStorage(n,month)
				+ routing.G_monthlyRiverStorage(n,month)
				+ routing.G_monthlyGwStorage(n,month);


			if (options.resOpt == 1)
				storages_routing += routing.G_monthlyResStorage(n,month);

			if (options.glacierOpt == 1)
			    storages_routing += routing.G_monthlyGlacierStorage(n,month);

			if (options.grid_store_TypeForStorages == 1)
                storages_routing /= numberOfDaysInMonth[month];

			//add all waterstorages and subtract wateruse         // [km³]
			G_monthly_totalWaterInStorages(n,month) =
					G_monthlyCanopyStorage(n,month)
					+ G_monthlySnowStorage(n,month)
					+ G_monthlySoilStorage(n,month)
					+ storages_routing;

			//convert total storage from km³ to mm
			G_monthly_totalWaterInStorages_mm(n,month) =
					G_monthly_totalWaterInStorages(n,month)
                    / (geo.areaOfCellByArrayPos(n) * (geo.G_contfreq[n] / 100.)) * 1000000.;
        }
	}

	if (options.outSingleStorages) {
		// all storages are monthly means [km3]

		// TOTAL_STORAGES
		G_monthlyCanopyStorage.write(outDir + "/G_CANOPY_WATER_STORAGE_km3_" + to_string(year) + ".12.UNF0");

		for (n = 0; n < ng; n++)
			for (short m = 0; m < 12; m++)
				if (0 == geo.G_contcell[n]) // ocean grid cells within land mask (e.g. Caspian Sea)
					G_monthlyArray(n,m) = 0.;
				else
					G_monthlyArray(n,m) = G_monthlyCanopyStorage(n,m) // TOTAL_STORAGES
						/ (geo.areaOfCellByArrayPos(n) * (geo.G_contfreq[n] / 100.)) * 1000000;

		G_monthlyArray.write(outDir + "/G_CANOPY_WATER_STORAGE_mm_" + to_string(year) + ".12.UNF0");

		// TOTAL_STORAGES
		G_monthlySnowStorage.write(outDir + "/G_SNOW_WATER_STORAGE_km3_" + to_string(year) + ".12.UNF0");

		for (n = 0; n < ng; n++)
			for (short m = 0; m < 12; m++)
				if (0 == geo.G_contcell[n]) // ocean grid cells within land mask (e.g. Caspian Sea)
					G_monthlyArray(n,m) = 0.;
				else
					G_monthlyArray(n,m) = G_monthlySnowStorage(n,m) // TOTAL_STORAGES
						/ (geo.areaOfCellByArrayPos(n) * (geo.G_contfreq[n] / 100.)) * 1000000;

		G_monthlyArray.write(outDir + "/G_SNOW_WATER_STORAGE_mm_" + to_string(year) + ".12.UNF0");

		// TOTAL_STORAGES
		G_monthlySoilStorage.write(outDir + "/G_SOIL_WATER_STORAGE_km3_" + to_string(year) + ".12.UNF0");

		for (n = 0; n < ng; n++)
			for (short m = 0; m < 12; m++)
				if (0 == geo.G_contcell[n]) // ocean grid cells within land mask (e.g. Caspian Sea)
					G_monthlyArray(n,m) = 0.;
				else
					G_monthlyArray(n,m) = G_monthlySoilStorage(n,m) // TOTAL_STORAGES
						/ (geo.areaOfCellByArrayPos(n) * (geo.G_contfreq[n] / 100.)) * 1000000;

		G_monthlyArray.write(outDir + "/G_SOIL_WATER_STORAGE_mm_" + to_string(year) + ".12.UNF0");
	}

    // write total monthly storage to file [km3]
    if (options.outTotalWaterInStorages_km3) {
		G_monthly_totalWaterInStorages.write(outDir + "/G_TOTAL_STORAGES_km3_" + to_string(year) + ".12.UNF0");
    }

    // [mm]
	if (options.outTotalWaterInStorages_mm) {
		for (n = 0; n < ng; n++)
			for (short m = 0; m < 12; m++)
                if (0 == geo.G_contcell[n])
					// ocean grid cells within land mask (e.g. Caspian Sea)
                    G_monthlyArray(n,m) = 0.;
                else
                    G_monthlyArray(n,m) = G_monthly_totalWaterInStorages_mm(n,m);

		G_monthlyArray.write(outDir + "/G_TOTAL_STORAGES_mm_" + to_string(year) + ".12.UNF0");
	}

	if (options.outLAI) {
		G_monthlyLAI.write(outDir + "/G_LAI_" + to_string(year) + ".12.UNF0");
	}
	if (options.outAlbedo) {
		G_monthlyAlbedo.write(outDir + "/G_ALBEDO_" + to_string(year) + ".12.UNF0");
	}
	if (options.outCanopyWater) {
		G_monthlycanopyWaterContent.write(outDir + "/G_CANOPY_WATER_mm_" + to_string(year) + ".12.UNF0");
	}
	if (options.outmaxCanopyWater) {
		G_monthlymaxcanopyWaterContent.write(outDir + "/G_MAX_CANOPY_WATER_" + to_string(year) + ".12.UNF0");
	}
	if (options.outSoilWater) {
		G_monthlySoilWaterAvg.write(outDir + "/G_SOIL_WATER_" + to_string(year) + ".12.UNF0");
	}
	if (options.outSnowCover) {
		G_monthlySnowCoverAvg.write(outDir + "/G_SNOW_COVER_FRAC_" + to_string(month2store) + "_" + to_string(year) + ".12.UNF0");
	}
	if (options.outSWE) {
		G_monthlySWE.write(outDir + "/G_SNOW_WATER_EQUIV_" + to_string(year) + ".12.UNF0");
	}
	if (options.outSnowFall) {
		G_monthlySnowFall.write(outDir + "/G_SNOW_FALL_" + to_string(year) + ".12.UNF0");
	}
    if (options.outSnowFall) {
		G_monthlyRainFall.write(outDir + "/G_RAIN_FALL_" + to_string(year) + ".12.UNF0");
    }
    if (options.outSnowMelt) {
		G_monthlySnowMelt.write(outDir + "/G_SNOW_MELT_" + to_string(year) + ".12.UNF0");
	}
	if (options.outSnowEvap) {
		G_monthlySnowEvap.write(outDir + "/G_SNOW_EVAP_" + to_string(year) + ".12.UNF0");
	}
	if (options.outNetRadiation) {
		G_monthlyNetRadiation.write(outDir + "/G_NET_RADIATION_" + to_string(year) + ".12.UNF0");
	}
    if (options.outNetShortWaveRadiation) {
		// net_short_wave_rad output
		G_monthlyNetShortRad.write(outDir + "/G_NET_SHORTWAVE_RAD_" + to_string(year) + ".12.UNF0");
		G_monthlyExtRad.write(outDir + "/G_EXT_RAD_" + to_string(year) + ".12.UNF0");
		G_monthlyShortUpRad.write(outDir + "/G_SHORTWAVE_UP_RAD_" + to_string(year) + ".12.UNF0");
		G_monthlyShortDownRad.write(outDir + "/G_SHORTWAVE_DOWN_RAD_" + to_string(year) + ".12.UNF0");
	}
    if (options.outNetLongWaveRadiation) {
		// net_long_wave_rad output
		G_monthlyNetLongRad.write(outDir + "/G_NET_LONGWAVE_RAD_" + to_string(year) + ".12.UNF0");
		G_monthlyLongDownRad.write(outDir + "/G_LONGWAVE_DOWN_RAD_" + to_string(year) + ".12.UNF0");
		G_monthlyLongUpRad.write(outDir + "/G_LONGWAVE_UP_RAD_" + to_string(year) + ".12.UNF0");
	}
    if (options.outSunshine) {
		// monthly sunshine output in % (0 if not driven with CRU data)
		G_monthlySunshine.write(outDir + "/G_SUNSHINE_" + to_string(year) + ".12.UNF0");
    }
    if (options.outTemperature) {
		// monthly temperature output in ?C
		G_monthlyTemperature.write(outDir + "/G_TEMPERATURE_" + to_string(year) + ".12.UNF0");
    }
}

void dailyWaterBalanceClass::writeDaily31Grids(const string outDir, const int year, const int month) {

	if (3 != options.grid_store && 4 != options.grid_store)
		return;

	if (options.outSurfaceRunoffDaily) {
		G_daily31SurfaceRunoff.write(outDir + "/G_SURFACE_RUNOFF_" + to_string(year) + "_" + to_string(month+1) +".31.UNF0");
	}
	if (options.outGWRechargeDaily) {
		G_daily31GwRecharge.write(outDir + "/G_GW_RECHARGE_" + to_string(year) + "_" + to_string(month+1) +".31.UNF0");
	}
	if (options.outSoilWaterDaily) {
		G_daily31SoilWaterAvg.write(outDir + "/G_SOIL_WATER_" + to_string(year) + "_" + to_string(month+1) +".31.UNF0");
	}
	if (options.outLAIDaily) {
		G_daily31LAI.write(outDir + "/G_LAI_" + to_string(year) + "_" + to_string(month+1) +".31.UNF0");
	}
	if ((options.outLAIDaily) && (options.use_kc == 1)) {
		G_daily31Kc.write(outDir + "/G_KC_" + to_string(year) + "_" + to_string(month+1) +".31.UNF0");
	}
	if (options.outAlbedoDaily) {
		G_daily31Albedo.write(outDir + "/G_ALBEDO_" + to_string(year) + "_" + to_string(month+1) +".31.UNF0");
	}
	if (options.outInterceptionDaily) {
		G_daily31Interception.write(outDir + "/G_INTERCEPTION_" + to_string(year) + "_" + to_string(month+1) +".31.UNF0");
	}
	if (options.outCanopyWaterDaily) {
		G_daily31canopyWaterContent.write(outDir + "/G_CANOPY_WATER_mm_" + to_string(year) + "_" + to_string(month+1) +".31.UNF0");
	}
	if (options.outmaxCanopyWaterDaily) {
		G_daily31maxcanopyWaterContent.write(outDir + "/G_MAX_CANOPY_WATER_" + to_string(year) + "_" + to_string(month+1) +".31.UNF0");
	}
    if (options.outPrecDaily) {
		G_daily31Precip.write(outDir + "/G_PRECIPITATION_" + to_string(year) + "_" + to_string(month+1) +".31.UNF0");
    }
	if (options.outPETDaily) {
		G_daily31PET.write(outDir + "/G_PET_" + to_string(year) + "_" + to_string(month+1) +".31.UNF0");
	}
	if (options.outTotalPETDaily) {
		G_daily31TotalPET.write(outDir + "/G_TOTAL_PET_" + to_string(year) + "_" + to_string(month+1) +".31.UNF0");
	}
	if (options.outSnowCoverDaily) {
		G_daily31SnowCoverAvg.write(outDir + "/G_SNOW_COVER_FRAC_" + to_string(year) + "_" + to_string(month+1) +".31.UNF0");
	}
	if (options.outSWEDaily) {
		G_daily31SWE.write(outDir + "/G_SNOW_WATER_EQUIV_" + to_string(year) + "_" + to_string(month+1) +".31.UNF0");
	}
	if (options.outSnowFallDaily) {
		G_daily31SnowFall.write(outDir + "/G_SNOW_FALL_" + to_string(year) + "_" + to_string(month+1) +".31.UNF0");
	}
	if (options.outSnowMeltDaily) {
		G_daily31SnowMelt.write(outDir + "/G_SNOW_MELT_" + to_string(year) + "_" + to_string(month+1) +".31.UNF0");
	}
    if (options.outSingleStoragesDaily) {
		G_daily31CanopyStorage.write(outDir + "/G_CANOPY_WATER_STORAGE_km3_" + to_string(year) + "_" + to_string(month+1) +".31.UNF0");
    }
    if (options.outSingleStoragesDaily) {
		G_daily31SnowStorage.write(outDir + "/G_SNOW_WATER_STORAGE_km3_" + to_string(year) + "_" + to_string(month+1) +".31.UNF0");
    }
    if (options.outSingleStoragesDaily) {
		G_daily31SoilStorage.write(outDir + "/G_SOIL_WATER_STORAGE_km3_" + to_string(year) + "_" + to_string(month+1) +".31.UNF0");
    }
	// radiation output adjustment - daily rad output anyway if daily outputs are selected.
    if (options.outExtRadiationDaily) {
		G_daily31ExtRad.write(outDir + "/G_EXT_RAD_" + to_string(year) + "_" + to_string(month+1) +".31.UNF0");
  	}
    if (options.outShortDownRadiationDaily) {
		G_daily31ShortDownRad.write(outDir + "/G_SHORTWAVE_DOWN_RAD_" + to_string(year) + "_" + to_string(month+1) +".31.UNF0");
  	}
    if (options.outShortUpRadiationDaily) {
		G_daily31ShortUpRad.write(outDir + "/G_SHORTWAVE_UP_RAD_" + to_string(year) + "_" + to_string(month+1) +".31.UNF0");
  	}
	if (options.outNetShortWaveRadiationDaily) {		
		G_daily31NetShortRad.write(outDir + "/G_NET_SHORTWAVE_RAD_" + to_string(year) + "_" + to_string(month+1) +".31.UNF0");
  	}
    if (options.outLongDownRadiationDaily) {
		G_daily31LongDownRad.write(outDir + "/G_LONGWAVE_DOWN_RAD_" + to_string(year) + "_" + to_string(month+1) +".31.UNF0");
  	}
    if (options.outLongUpRadiationDaily) {
		G_daily31LongUpRad.write(outDir + "/G_LONGWAVE_UP_RAD_" + to_string(year) + "_" + to_string(month+1) +".31.UNF0");
  	}
	if (options.outNetLongWaveRadiationDaily) {		
		G_daily31NetLongRad.write(outDir + "/G_NET_LONG_WAVE_RAD_" + to_string(year) + "_" + to_string(month+1) +".31.UNF0");
  	}
	if (options.outNetRadiationDaily) {		
		G_daily31NetRadiation.write(outDir + "/G_NET_RADIATION_" + to_string(year) + "_" + to_string(month+1) +".31.UNF0");
	}
}

void dailyWaterBalanceClass::writeDaily365Grids(const string outDir, const int year) {

	if (5 != options.grid_store && 6 != options.grid_store)
		return;

    extern routingClass routing;
    extern geoClass geo;

	if (options.outSurfaceRunoffDaily) {
		G_daily365SurfaceRunoff.write(outDir + "/G_SURFACE_RUNOFF_" + to_string(year) + ".365.UNF0");
	}
	if (options.outGWRechargeDaily) {
		G_daily365GwRecharge.write(outDir + "/G_GW_RECHARGE_" + to_string(year) + ".365.UNF0");
	}
	if (options.outSoilWaterDaily) {
		G_daily365SoilWaterAvg.write(outDir + "/G_SOIL_WATER_" + to_string(year) + ".365.UNF0");
	}
	if (options.outLAIDaily) {
		G_daily365LAI.write(outDir + "/G_LAI_" + to_string(year) + ".365.UNF0");
	}
	if ((options.outLAIDaily) && (options.use_kc == 1)) {
		G_daily365Kc.write(outDir + "/G_KC_" + to_string(year) + ".365.UNF0");
	}
	if (options.outAlbedoDaily) {
		G_daily365Albedo.write(outDir + "/G_ALBEDO_" + to_string(year) + ".365.UNF0");
	}
	if (options.outInterceptionDaily) {
		G_daily365Interception.write(outDir + "/G_INTERCEPTION_" + to_string(year) + ".365.UNF0");
	}
	if (options.outCanopyWaterDaily) {
		G_daily365canopyWaterContent.write(outDir + "/G_CANOPY_WATER_mm_" + to_string(year) + ".365.UNF0");
	}
	if (options.outmaxCanopyWaterDaily) {
		G_daily365maxcanopyWaterContent.write(outDir + "/G_MAX_CANOPY_WATER_" + to_string(year) + ".365.UNF0");
	}
    if (options.outPrecDaily) {
		// daily precip output
		G_daily365Precip.write(outDir + "/G_PRECIPITATION_" + to_string(year) + ".365.UNF0");
	}
    if (options.outTemp) {
        // daily temperature output
        G_daily365temp.write(outDir + "/G_TEMPERATURE_" + to_string(year) + ".365.UNF0");
    }
    if (options.outRedTemp&&(7==options.petOpt)) {
        // daily reduced temperature output
        G_daily365redtemp.write(outDir + "/G_REDUCED_TEMPERATURE_" + to_string(year) + ".365.UNF0");
    }
	if (options.outPETDaily) {
		G_daily365PET.write(outDir + "/G_PET_" + to_string(year) + ".365.UNF0");
	}
	if (options.outTotalPETDaily) {
		G_daily365TotalPET.write(outDir + "/G_TOTAL_PET_" + to_string(year) + ".365.UNF0");
	}
	if (options.outSnowCoverDaily) {
		G_daily365SnowCoverAvg.write(outDir + "/G_SNOW_COVER_FRAC_" + to_string(year) + ".365.UNF0");
	}
	if (options.outSWEDaily) {
		G_daily365SWE.write(outDir + "/G_SNOW_WATER_EQUIV_" + to_string(year) + ".365.UNF0");
	}
	if (options.outSnowFallDaily) {
		G_daily365SnowFall.write(outDir + "/G_SNOW_FALL_" + to_string(year) + ".365.UNF0");
	}
	if (options.outSnowMeltDaily) {
		G_daily365SnowMelt.write(outDir + "/G_SNOW_MELT_" + to_string(year) + ".365.UNF0");
    }
    if (options.outSingleStoragesDaily) {
		G_daily365CanopyStorage.write(outDir + "/G_CANOPY_WATER_STORAGE_km3_" + to_string(year) + ".365.UNF0");
    }
    if (options.outSingleStoragesDaily) {
		G_daily365SnowStorage.write(outDir + "/G_SNOW_WATER_STORAGE_km3_" + to_string(year) + ".365.UNF0");
    }
    if (options.outSingleStoragesDaily) {
		G_daily365SoilStorage.write(outDir + "/G_SOIL_WATER_STORAGE_km3_" + to_string(year) + ".365.UNF0");
    }

	// radiation output adjustment - daily rad output anyway if daily outputs are selected.
    if (options.outExtRadiationDaily) {
		G_daily365ExtRad.write(outDir + "/G_EXT_RAD_" + to_string(year) + ".365.UNF0");
  	}
    if (options.outShortDownRadiationDaily) {
		G_daily365ShortDownRad.write(outDir + "/G_SHORTWAVE_DOWN_RAD_" + to_string(year) + ".365.UNF0");
  	}
    if (options.outShortUpRadiationDaily) {
		G_daily365ShortUpRad.write(outDir + "/G_SHORTWAVE_UP_RAD_" + to_string(year) + ".365.UNF0");
  	}
	if (options.outNetShortWaveRadiationDaily) {		
		G_daily365NetShortRad.write(outDir + "/G_NET_SHORTWAVE_RAD_" + to_string(year) + ".365.UNF0");
  	}
    if (options.outLongDownRadiationDaily) {
		G_daily365LongDownRad.write(outDir + "/G_LONGWAVE_DOWN_RAD_" + to_string(year) + ".365.UNF0");
  	}
    if (options.outLongUpRadiationDaily) {
		G_daily365LongUpRad.write(outDir + "/G_LONGWAVE_UP_RAD_" + to_string(year) + ".365.UNF0");
  	}
	if (options.outNetLongWaveRadiationDaily) {		
		G_daily365NetLongRad.write(outDir + "/G_NET_LONGWAVE_RAD_" + to_string(year) + ".365.UNF0");
  	}
	if (options.outNetRadiationDaily) {		
		G_daily365NetRadiation.write(outDir + "/G_NET_RADIATION_" + to_string(year) + ".365.UNF0");
	}
}

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

//    extern geoClass geo;

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
}
