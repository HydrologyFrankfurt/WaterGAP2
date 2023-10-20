
/***********************************************************************
*
*new attributes due to new land classes used and thrown out sensitivity
*analysis (Stephanie Eisner)
*
* see former changes at file s_max.cpp.versioninfos.txt
*
***********************************************************************/

#include <cstdio>
#include "grid.h"
#include "def.h"
#include "common.h"
#include "globals.h"
// Cell-specific calibration parameters
#include "calib_param.h"
//extern calibParamClass calibParam; // JSON object and methods from watergap.cpp

using namespace std;

soilWatCapClass::soilWatCapClass(void)
{
	parameterFlag = false;
}

void soilWatCapClass::setParameters()
{
	// CLC & GLCT maps combined
	const double rootingDepthOrig[nlct] = { 2.0, 4.0, 2.0, 2.0, 2.0, 1.0, 0.5, 1.5, 1.5, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.1, 1.0, 2.0 };

	// keep original values 
	for (short n = 0; n <= nlct - 1; n++)
		rootingDepth[n] = (float)rootingDepthOrig[n];

	parameterFlag = true;
}


void soilWatCapClass::createMaxSoilWaterCapacityGrid(const std::string input_dir, const std::string output_dir, const Grid<char> G_land_cover, calibParamClass calParam)//OE: calParam added
{
	// parameters in dailyWaterBalanceClass read from 'LCT_22.DAT'
	setParameters();

	// total available soil water capacity in 1 m soil [mm] 
        // from BATJES (1996) (at old model versions)
        // since 2.2b it is from WISE-Database; UNF-file is renamed to TAWC
	//                                                       
	// cells which have no value in the database contain     
	// the value -9999, these are cells with ice-cover       
	// (Greenland and some islands)                           
	Grid<float> G_TAWC;
	Grid<float> G_RootDepth;

	G_TAWC.read(input_dir + "/G_TAWC.UNF0");

	// calculate maximum soil water capacity
	// by considering land cover specific rooting depth
	for (int n = 0; n < ng; n++) {
		// Cell-specific calibration parameters - Apply multiplier  // FP
		G_RootDepth[n] = calParam.getValue(M_ROOT_D,n) * dailyWaterBalance.get_rootingDepth_lct(G_land_cover[n] - 1);
		if (G_TAWC[n] < 0)	// should only occur as NoData value (-9999)
			G_Smax[n] = -9999;
		else
			G_Smax[n] = G_TAWC[n] * G_RootDepth[n];
	}
	// write grid for rooting depth
	if (options.outRootingDepth) {
		G_RootDepth.write(output_dir + "/G_ROOT_DEPTH.UNF0");
	}
	// write grid for maximum soil water capacity
	if (options.outmaxSoilWater) {
		G_Smax.write(output_dir + "/G_Smax.UNF0");
	}
}
