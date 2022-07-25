
/***********************************************************************
*
*new attributes due to new land classes used and thrown out sensitivity
*analysis (Stephanie Eisner)
*
* see former changes at file s_max.cpp.versioninfos.txt
*
***********************************************************************/

#include <cstdio>
#include "gridio.h"
#include "def.h"
#include "option.h"

#include "s_max.h"
#include "daily.h"

// #######################################
// Cell-specific calibration parameters
// FP
#include "calib_param.h"
extern calibParamClass calibParam; // JSON object and methods from watergap.cpp
// #######################################

using namespace std;

soilWatCapClass::soilWatCapClass(void)
{
	parameterFlag = false;
}

void soilWatCapClass::setParameters()
{
	//const double rootingDepthOrig[nlct] = { 1.0, 1.0, 1.0, 1.0, 2.0, 1.0, 0.5, 1., 2.0, 2.0, 2.0, 2.0, 2.0, 1.0, 0.1, 1.5, 1.5, 2.0, 4.0, 1.0 };
	// CLC & GLCT maps combined
	const double rootingDepthOrig[nlct] = { 2.0, 4.0, 2.0, 2.0, 2.0, 1.0, 0.5, 1.5, 1.5, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.1, 1.0, 2.0 };

	// keep original values 
	for (short n = 0; n <= nlct - 1; n++)
		rootingDepth[n] = (float)rootingDepthOrig[n];

	parameterFlag = true;
}


void soilWatCapClass::createMaxSoilWaterCapacityGrid(char *input_dir, char *output_dir, char *G_land_cover)
{
	extern gridioClass gridIO;
	extern optionClass options;
	extern dailyWaterBalanceClass dailyWaterBalance;

	// parameters in dailyWaterBalanceClass read from 'LCT_22.DAT' 
	// if (!parameterFlag)
	setParameters();

	// total available soil water capacity in 1 m soil [mm] 
        // from BATJES (1996) (at old model versions)
        // since 2.2b it is from WISE-Database; UNF-file is renamed to TAWC
	//                                                       
	// cells which have no value in the database contain     
	// the value -9999, these are cells with ice-cover       
	// (Greenland and some islands)                           
	float G_TAWC[ng];
	float G_RootDepth[ng];
	char filename[250];

	sprintf(filename, "%s/G_TAWC.UNF0", input_dir);
	gridIO.readUnfFile(filename, ng, G_TAWC);

	// calculate maximum soil water capacity
	// by considering land cover specific rooting depth
	for (int n = 0; n < ng; n++) {
		// Cell-specific calibration parameters - Apply multiplier  // FP
		G_RootDepth[n] = calibParam.getValue(M_ROOT_D,n) * dailyWaterBalance.get_rootingDepth_lct(G_land_cover[n] - 1);
		if (G_TAWC[n] < 0)	// should only occur as NoData value (-9999)
			G_Smax[n] = -9999;
		else
			G_Smax[n] = G_TAWC[n] * G_RootDepth[n];
	}
	// write grid for rooting depth
	if (options.outRootingDepth) { // new output options 2.2
		sprintf(filename, "%s/G_ROOT_DEPTH.UNF0", output_dir);
		gridIO.writeUnfFile(filename, ng, G_RootDepth);
	}
	// write grid for maximum soil water capacity
	if (options.outmaxSoilWater) { // new output options 2.2
		sprintf(filename, "%s/G_Smax.UNF0", output_dir);
		gridIO.writeUnfFile(filename, ng, G_Smax);
	}
}
