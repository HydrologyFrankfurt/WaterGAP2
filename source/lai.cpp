
/***********************************************************************
*
*land use attributes are read in from LCT_22.DAT now (enhanced), LAI is
*calculated based on land use class, growing period is based on precipitation
*and not on PET anymore, all done by someone from CESR
*
*introduced crop coefficients as alternative for using different albedo when
*calculating PET (by Stephanie Eisner)
*
* see former changes at file lai.cpp.versioninfos.txt
*
***********************************************************************/
//#define store_grids

#include <iostream>
#include <fstream>
#include <cstdio>

#include "gridio.h"
#include "def.h"
#include "option.h"

#include "lai.h"

using namespace std;

extern optionClass options;

// #######################################
// Cell-specific calibration parameters
// FP
#include "calib_param.h"
extern calibParamClass calibParam; // JSON object and methods from watergap.cpp
// #######################################


laiClass::laiClass(void)
{
	LAIdataIsRead = false;
}

void laiClass::init(const char *inputDir, char *output_dir, const char *G_land_cover)
{
	// this initialization routine should be called once 
	// at the beginning of the simulation 
	// (for normal simulations or calibrations) 
	// or 
	// in case of a sensitivity analysis each time when 
	// starting with a new sample set.

	extern gridioClass gridIO;
	int n;  // Index over grid cells
	int i_lct;  // Index over land cover classes (replacing inconsistent use of n) // FP

	for (n = 0; n < ng; n++){
		G_days_since_start[n] = 0;
		G_GrowingStatus[n] = 0;
		G_PrecSum[n] = 0.;
	}

	char filename[250];
// not used any more in WG2.2	
//	float G_biomass[ng_land][9];
//	sprintf(filename, "%s/GBIOMASS_1970.9.UNF0", inputDir);	// from IMAGE
//	gridIO.readUnfFile(filename, 9 * ng_land, &G_biomass[0][0]);
// not used any more in WG2.2
	
	if (!LAIdataIsRead)
		readLAIdata(inputDir);

	// in case of a sensitivity analysis these values are modified.
	// we want to keep the original values.
// not used any more in WG2.2	
//	float specificLeafArea[nlct];	// [m2/kg dry mass]
//	float leafMassCorr[nlct];	// [-]
// not used any more in WG2.2	
	float LeafAreaIndex[nlct];	// [-]			WaterGAP2.2
	float decidiousPlantFraction[nlct];	// [-]
	float evergreensLAIreduction[nlct];	// [-]

	for (i_lct = 0; i_lct < nlct; i_lct++) {
// not used any more in WG2.2	
//		specificLeafArea[i_lct] = specificLeafAreaOrig[i_lct];
//		leafMassCorr[i_lct] = leafMassCorrOrig[i_lct];
// not used any more in WG2.2	
		/* WaterGAP2.2 */
		LeafAreaIndex[i_lct] = LeafAreaIndexOrig[i_lct]; // WaterGAP2.2
		decidiousPlantFraction[i_lct] = decidiousPlantFractionOrig[i_lct];
		evergreensLAIreduction[i_lct] = evergreensLAIreductionOrig[i_lct];
	}

// not used any more in WG2.2	
	// calculate cell-specific maximum LAI
//	for (n = 0; n < ng_land; n++) {
//		G_LAImax[n]
//			= specificLeafArea[G_land_cover[n] - 1]
//			* leafMassCorr[G_land_cover[n] - 1]
//			* (32.0 / 14000.0)	// convert gC/m2 to kg dry mass/m2 
//			* G_biomass[n][2];	// 3. column is leaf mass in gC/m2
//	}
// not used any more in WG2.2	

	/* WaterGAP2.2 */
	// this is the LAI in WaterGAP2.2 - not calculated from 
	// specific LAI * Biomass any more as in WaterGAP3.0
        //for (n = 0; n < ng_land; n++)
		for (n = 0; n < ng; n++) {
			// Cell-specific calibration parameters - Apply multiplier  // FP
			G_LAImax[n] = calibParam.getValue(M_LAI,n) * LeafAreaIndex[G_land_cover[n] - 1];
		}

        //for (n = ng_land; n < ng; n++)
        //	G_LAImax[n] = 0.;	// lake cells

	// these factors will be used to calculate cell-specific 
	// mimimum LAI. 
	// daily calculations should be kept as easy as possible.
	for (i_lct = 0; i_lct < nlct; i_lct++) {
		lai_factor_a[i_lct] = 0.1 * decidiousPlantFraction[i_lct];
		lai_factor_b[i_lct] = (1 - decidiousPlantFraction[i_lct]) * evergreensLAIreduction[i_lct];
	}

	if (options.outLAImax) { // new output options 2.2
		sprintf(filename, "%s/G_LAI_MAX.UNF0", output_dir);
		gridIO.writeUnfFile(filename, ng, G_LAImax);
	
		float G_LAImin[ng];
	
                //for (n = 0; n <= ng_land - 1; n++)
                for (n = 0; n <= ng - 1; n++)
			G_LAImin[n] = lai_factor_a[G_land_cover[n] - 1] + lai_factor_b[G_land_cover[n] - 1] * G_LAImax[n];
	
                //for (n = ng_land; n < ng; n++)
                //	G_LAImin[n] = 0;	// lake cells
	
		sprintf(filename, "%s/G_LAI_MIN.UNF0", output_dir);
		gridIO.writeUnfFile(filename, ng, G_LAImin);
	}

}
// end init (with calling readLAIdata())

void laiClass::readLAIdata(const char *inputDir)
{
	char filename[250];

	// open LAI.DAT 
// not used any more in WG2.2	
//	sprintf(filename, "%s/LAI.DAT", inputDir);
// not used any more in WG2.2	
	/* WaterGAP2.2 */
	sprintf(filename, "%s/LAI_22.DAT", inputDir);
	ifstream lai_file(filename);

	if (!lai_file) {
		cerr << "Can not open file " << filename << " for reading.\n";
		exit(-1);
	}
	// read commentlines indicated by #
	char string[250];

	while (!lai_file.eof() && lai_file.peek() == '#') {
		lai_file.getline(string, sizeof(string));
	}

	// read values
	short dummy_int;

	for (short i = 0; i <= nlct - 1; i++) {
		lai_file >> dummy_int;  // Reading index of land cover classes 1 ... nlct
		if ((dummy_int - 1) != i) {
			cerr << "Problem with " << filename;
			exit(-1);
		}
// not used any more in WG2.2	
//		lai_file >> specificLeafAreaOrig[i];	// [m2/kg dry mass]
//		lai_file >> leafMassCorrOrig[i];
// not used any more in WG2.2	
		/* WaterGAP2.2 */
		lai_file >> LeafAreaIndexOrig[i];	// [-]
		lai_file >> decidiousPlantFractionOrig[i];
		lai_file >> evergreensLAIreductionOrig[i];
		lai_file >> initialDays[i];
		lai_file >> kc_min[i];
		lai_file >> kc_max[i];
	}
	LAIdataIsRead = true;
	
}
// end readLAIdata() of file LAI_22.DAT

/* new version */
// use threshold of precipitation sum instead of dailyPET for decision of growing conditions
//double laiClass::getDailyLai(int n, char land_cover_type, bool arid_humid, double daily_temp_C,	 double daily_precipitation, double dailyPET)
double laiClass::getDailyLai(int n, char land_cover_type, bool arid, double daily_temp_C, double daily_precipitation)
{
	double LAImin = lai_factor_a[land_cover_type - 1] + lai_factor_b[land_cover_type - 1] * G_LAImax[n];
	double dailyLAI;
	
/* replace second condition with precipitation sum in growing/nogrowing function
	// in arid regions with distinct dry season growing period starts when precipitation equals half PET
	if (arid_humid) {
		if ((daily_temp_C > 5.) && ((2. * daily_precipitation) > dailyPET)) {
			// growing conditions
			dailyLAI = dailyLai_growing(&G_days_since_start[n], initialDays[land_cover_type - 1], &G_GrowingStatus[n], land_cover_type, LAImin, G_LAImax[n]);
			// for arid regions there is no initial phase
			//dailyLAI = dailyLai_growing(&G_days_since_start[n], 0, &G_GrowingStatus[n], LAImin, G_LAImax[n]);
		}
		else {
			// no growing conditions
			dailyLAI = dailyLai_nogrowing(&G_days_since_start[n], initialDays[land_cover_type - 1], &G_GrowingStatus[n], land_cover_type, LAImin, G_LAImax[n]);
			// for arid regions there is no initial phase
			//dailyLAI = dailyLai_nogrowing(&G_days_since_start[n], 0, &G_GrowingStatus[n], LAImin, G_LAImax[n]);
		}
	}
	// in humid areas second condition is not necessary
	else {
*/
		//if (daily_temp_C > 5.) {
		if (daily_temp_C > 8.) {
			// growing conditions
			dailyLAI = dailyLai_growing(&G_days_since_start[n], initialDays[land_cover_type - 1], &G_GrowingStatus[n], land_cover_type, arid, LAImin, G_LAImax[n], &G_PrecSum[n], daily_precipitation);
		}
		else {
			// no growing conditions
			dailyLAI = dailyLai_nogrowing(&G_days_since_start[n], initialDays[land_cover_type - 1], &G_GrowingStatus[n], land_cover_type, arid, LAImin, G_LAImax[n], &G_PrecSum[n], daily_precipitation);
		}
//	}
	return dailyLAI;
}

double laiClass::dailyLai_growing(char *days_since_start, char initialDays, char *GrowingStatus, char land_cover_type, bool arid, double LAImin, double LAImax, double *PrecSum, double prec)
{
	if (*GrowingStatus == 0) { // vegetation is not grown yet...
		// ...but initial days of growing season have been reached
		if (*days_since_start >= initialDays){
			(*days_since_start)++;
			(*PrecSum) += prec;
			// threshold of min precipitation has been met and growing season can start
			if(*PrecSum > 40.){
				// full vegetation development has been reached and status will switch
				if (*days_since_start >= initialDays + 30){
					*days_since_start = initialDays + 30;
					*GrowingStatus = 1;
				}
				// return either the LAI in sprout or the max LAI-value
				return (LAImin + (LAImax - LAImin) * (*days_since_start - initialDays) / 30.);
			}
			// sum of precipitation is not high enough to start growing season
			else {
				*days_since_start = initialDays;
				return (LAImin);
			}
		}
		// initial days of growing season have not been reached yet
		else {
			(*days_since_start)++;
			(*PrecSum) += prec;
			return (LAImin);
		}
	}
	else{// vegetation is already fully grown...
		if (*days_since_start <= 30) {
			// ...and we are in phase of senescence
			(*days_since_start)--;
			// if land_cover_type is '1' or '2' we have ervergreen plants and LAI will never clompletely degrade
			// therefore status will switch at once
			if (land_cover_type <= 2){
				*GrowingStatus = 0;
				//*PrecSum = 0.;
			}
			// LAI is completely degraded to min LAI and status will switch
			if (*days_since_start <= 0){
				*days_since_start = 0;
				*GrowingStatus = 0;
				*PrecSum = 0.;
			}
			return (LAImax - (LAImax - LAImin) * (30 - *days_since_start) / 30.);
		} else {
			// initial days for senescence phase have not been reached yet and we have growing conditions (again)
			// in arid regions we have no growing conditions if there is no rain 
			if ( (arid) && (prec<0.5))
				(*days_since_start)--;	
			// if conditions for LAI reduction should happen day after day, use this  equation (reset for initial days)
			else
				*days_since_start = 30 + initialDays;
			return (LAImax);
		}
	}
	
}

double laiClass::dailyLai_nogrowing(char *days_since_start, char initialDays, char *GrowingStatus, char land_cover_type, bool arid, double LAImin, double LAImax, double *PrecSum, double prec)
{

	if (*GrowingStatus == 0) { // vegetation is not fully grown yet 
		// initial days of growing season have been reached and plants will grow anyway
		//if (*days_since_start >= initialDays){
		if (*days_since_start > initialDays){
			(*days_since_start)++;
			(*PrecSum) += prec;
			// threshold of min precipitation has been met and growing season can start
			if(*PrecSum > 40.){
				if (*days_since_start >= initialDays + 30){
					*days_since_start = initialDays + 30;
					*GrowingStatus = 1;
				}
				return (LAImin + (LAImax - LAImin) * (*days_since_start - initialDays) / 30.);
			}
			// sum of precipitation is not high enough to start growing season
			else {
				*days_since_start = initialDays;
				return (LAImin);
			}
		}
		// no growing season
		else {
			(*PrecSum) += prec;
			return (LAImin);
		}
	}
	else{
		// we are in growing season but lai will be reduced now 
		if (*days_since_start <= 30) {
			(*days_since_start)--;
			// if land_cover_type is '1' or '2' we have ervergreen plants and LAI will never clompletely degrade
			// therefore status will switch at once
			//if (land_cover_type <= 2)
			//	*GrowingStatus = 0; 
			// LAI will completely degrade to min LAI and status will switch
			if (*days_since_start <= 0){
				*days_since_start = 0;
				*GrowingStatus = 0;
				*PrecSum = 0.;
			}
			return (LAImax - (LAImax - LAImin) * (30 - *days_since_start) / 30.);
		}
		// we are in growing season but the inital days for lai degrading have not been met 
		else {
			(*days_since_start)--;
			return (LAImax);
		}
	}
	
}
/* old version of daily lai calculation
double laiClass::getDailyLai(int n, char land_cover_type, bool aird_humid, double daily_temp_C,	 double daily_precipitation, double dailyPET)
{
	double LAImin;

	if ((daily_temp_C > 5.0) && ((2. * daily_precipitation) > dailyPET)) {
		// growing day
		if (30 == G_days_since_start[n]) {
			// previous day also was a growing day
			return G_LAImax[n];
		} else {
			if (G_days_since_start[n] >= 0) {
				// previous day was during the first 30 days of the growing period
				LAImin = lai_factor_a[land_cover_type - 1] + lai_factor_b[land_cover_type - 1] * G_LAImax[n];
				G_days_since_start[n]++;
				return (LAImin + (G_LAImax[n] - LAImin) * G_days_since_start[n] / 30.);
			} else {
				// if value is negative, previous day has been 
				// 'after' the end of the growing period
				// but now we are back in the growing period !!
				G_days_since_start[n] = 30 + G_days_since_start[n] + 1;
				LAImin = lai_factor_a[land_cover_type - 1] + lai_factor_b[land_cover_type - 1] * G_LAImax[n];
				return (LAImin + (G_LAImax[n] - LAImin) * G_days_since_start[n] / 30.);
			}
		}
	} else {
		// no growing day
		if (0 == G_days_since_start[n]) {
			// previous day also was no growing day
			return (lai_factor_a[land_cover_type - 1] + lai_factor_b[land_cover_type - 1] * G_LAImax[n]);
		} else {
			if (G_days_since_start[n] < 0) {
				// previous day was during the 30 days 
				// following the end of a growing period
				LAImin = lai_factor_a[land_cover_type - 1] + lai_factor_b[land_cover_type - 1] * G_LAImax[n];
				if (-29 == G_days_since_start[n]) {
					G_days_since_start[n] = 0;
					return LAImin;
				} else {
					G_days_since_start[n]--;
					return (G_LAImax[n] - (G_LAImax[n] - LAImin) * (-1. * G_days_since_start[n]) / 30.);
				}
			} else {
				// if value is positive, previous day has been 
				// within the first 30 days of a growing period
				// but now the period is already over
				G_days_since_start[n] = -30 + G_days_since_start[n] - 1;
				LAImin = lai_factor_a[land_cover_type - 1] + lai_factor_b[land_cover_type - 1] * G_LAImax[n];
				return (G_LAImax[n] - (G_LAImax[n] - LAImin) * (-1. * G_days_since_start[n]) / 30.);
			}
		}
	}
}
*/
double laiClass::getDailyKc(int n, char land_cover_type, double dailyLAI)
{
	double LAImin = lai_factor_a[land_cover_type - 1] + lai_factor_b[land_cover_type - 1] * G_LAImax[n];
	double dailyKc;
	
	if ((G_LAImax[n]-LAImin)==0.) dailyKc = kc_min[land_cover_type - 1];
	else	
		dailyKc = kc_min[land_cover_type - 1] + 
			(kc_max[land_cover_type - 1]-kc_min[land_cover_type - 1])* (dailyLAI-LAImin)/(G_LAImax[n]-LAImin);
		
	return dailyKc;
}
