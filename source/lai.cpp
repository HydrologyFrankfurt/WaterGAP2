
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

//#include <iostream>
//#include <fstream>
#include <cstdio>
#include "common.h"
#include "def.h"
//#include "option.h"
//#include "lai.h"
#include "globals.h"
#include "calib_param.h"

using namespace std;

//extern optionClass options;

// Cell-specific calibration parameters
#include "calib_param.h"
extern calibParamClass calibParam;


laiClass::laiClass(void)
{
	LAIdataIsRead = false;
}

void laiClass::init(const std::string inputDir, const std::string output_dir, const Grid<char> G_land_cover,AdditionalOutputInputFile &additionalOutIn,calibParamClass &calParam)
{
	// this initialization routine should be called once 
	// at the beginning of the simulation 
	// (for normal simulations or calibrations) 
	// or 
	// in case of a sensitivity analysis each time when 
	// starting with a new sample set.

	int n;  // Index over grid cells
	int i_lct;  // Index over land cover classes (replacing inconsistent use of n)

    // HMS if cda-mode, read in from additionalOutputInput, otherwise (normal run) set to 0
    if (additionalOutIn.additionalfilestatus == 0) { // no additionaloutin is used (normal run)
        for (n = 0; n < ng; n++){
            G_days_since_start[n] = 0;
            G_GrowingStatus[n] = 0;
            G_PrecSum[n] = 0.;
        }
    } else {
        for (n = 0; n < ng; n++) {
            G_days_since_start[n] = additionalOutIn.additionalOutputInput(n,0);
            G_GrowingStatus[n] = additionalOutIn.additionalOutputInput(n,1);
            G_PrecSum[n] = additionalOutIn.additionalOutputInput(n,2);
        }
    }
	
	if (!LAIdataIsRead)
		readLAIdata(inputDir);

	// in case of a sensitivity analysis these values are modified.
	// we want to keep the original values.
	float LeafAreaIndex[nlct];	// [-]
	float decidiousPlantFraction[nlct];	// [-]
	float evergreensLAIreduction[nlct];	// [-]

	for (i_lct = 0; i_lct < nlct; i_lct++) {
		LeafAreaIndex[i_lct] = LeafAreaIndexOrig[i_lct];
		decidiousPlantFraction[i_lct] = decidiousPlantFractionOrig[i_lct];
		evergreensLAIreduction[i_lct] = evergreensLAIreductionOrig[i_lct];
	}

	// this is the LAI in WaterGAP2.2 - not calculated from 
	// specific LAI * Biomass any more as in WaterGAP3.0
	for (n = 0; n < ng; n++) {
		// Cell-specific calibration parameters - Apply multiplier
		G_LAImax[n] = calParam.getValue(M_LAI,n) * LeafAreaIndex[G_land_cover[n] - 1];
	}

	// these factors will be used to calculate cell-specific 
	// mimimum LAI. 
	// daily calculations should be kept as easy as possible.
	for (i_lct = 0; i_lct < nlct; i_lct++) {
		lai_factor_a[i_lct] = 0.1 * decidiousPlantFraction[i_lct];
		lai_factor_b[i_lct] = (1 - decidiousPlantFraction[i_lct]) * evergreensLAIreduction[i_lct];
	}

	if (options.outLAImax) {
		G_LAImax.write(output_dir + "/G_LAI_MAX.UNF0");
	
		Grid<float> G_LAImin;
	
        for (n = 0; n <= ng - 1; n++)
			G_LAImin[n] = lai_factor_a[G_land_cover[n] - 1] + lai_factor_b[G_land_cover[n] - 1] * G_LAImax[n];

		G_LAImin.write(output_dir + "/G_LAI_MIN.UNF0");
	}

}

void laiClass::readLAIdata(const std::string inputDir)
{
	char filename[250];

	// open LAI.DAT 
	sprintf(filename, "%s/LAI_22.DAT", inputDir.c_str());
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

		lai_file >> LeafAreaIndexOrig[i];	// [-]
		lai_file >> decidiousPlantFractionOrig[i];
		lai_file >> evergreensLAIreductionOrig[i];
		lai_file >> initialDays[i];
		lai_file >> kc_min[i];
		lai_file >> kc_max[i];
	}
	LAIdataIsRead = true;
	
}

// use threshold of precipitation sum instead of dailyPET for decision of growing conditions
double laiClass::getDailyLai(int n, char land_cover_type, bool arid, double daily_temp_C, double daily_precipitation, const short day, const short last_day_in_month, AdditionalOutputInputFile &additionalOutIn)
{
	double LAImin = lai_factor_a[land_cover_type - 1] + lai_factor_b[land_cover_type - 1] * G_LAImax[n];
	double dailyLAI;
	
	if (daily_temp_C > 8.) {
		// growing conditions
		dailyLAI = dailyLai_growing(&G_days_since_start[n], initialDays[land_cover_type - 1], &G_GrowingStatus[n], land_cover_type, arid, LAImin, G_LAImax[n], &G_PrecSum[n], daily_precipitation);
	}
	else {
		// no growing conditions
		dailyLAI = dailyLai_nogrowing(&G_days_since_start[n], initialDays[land_cover_type - 1], &G_GrowingStatus[n], land_cover_type, arid, LAImin, G_LAImax[n], &G_PrecSum[n], daily_precipitation);
	}

    //if (additionalOutIn.additionalfilestatus ==1) { //HMS save only in case of C/DA runs
    // save G_PrecSum, G_GrowingStatus and G_days_since_start to additionalOutIn
    if (day == last_day_in_month + 1) // for the last day of a month
    {
        additionalOutIn.additionalOutputInput(n,0) = G_days_since_start[n];
        additionalOutIn.additionalOutputInput(n,1) = G_GrowingStatus[n];
        additionalOutIn.additionalOutputInput(n,2) = G_PrecSum[n];

    }

	return dailyLAI;
}

double laiClass::dailyLai_growing(int *days_since_start, short initialDays, int *GrowingStatus, char land_cover_type, bool arid, double LAImin, double LAImax, double *PrecSum, double prec)
{
	if (*GrowingStatus == 0) {

		// vegetation is not grown yet...
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
	else{

		// vegetation is already fully grown...
		if (*days_since_start <= 30) {
			// ...and we are in phase of senescence
			(*days_since_start)--;
			// if land_cover_type is '1' or '2' we have ervergreen plants and LAI will never clompletely degrade
			// therefore status will switch at once
			if (land_cover_type <= 2){
				*GrowingStatus = 0;
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

double laiClass::dailyLai_nogrowing(int *days_since_start, short initialDays, int *GrowingStatus, char land_cover_type, bool arid, double LAImin, double LAImax, double *PrecSum, double prec)
{
	if (*GrowingStatus == 0) {
		// vegetation is not fully grown yet
		// initial days of growing season have been reached and plants will grow anyway
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