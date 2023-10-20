/***********************************************************************
*
* see former changes at file gw_frac.cpp.versioninfos.txt
*
***********************************************************************/

//#include <iostream>
#include <cstdio>
#include <cmath>

#include "common.h"
#include "def.h"
//#include "option.h"
//#include "gw_frac.h"
#include "globals.h"
#include "calib_param.h"

using namespace std;

// Cell-specific calibration parameters
#include "calib_param.h"
extern calibParamClass calibParam;

inline short round_to_short(const float value)
{
	return (short) floor(value + 0.5);
}


groundwaterFactorClass::groundwaterFactorClass()
{
	//nothing to do
}


void groundwaterFactorClass::createGrids(const std::string inputDirectory, const std::string outputDirectory,calibParamClass &calParam)
{
	extern optionClass options;

	int n;

	// SLOPE FACTOR

    Grid<char> G_slope_class;	// slope class [10-70]
	G_slope_class.read(inputDirectory + "/G_SLOPE_CLASS.UNF1");

	// file has value 0 if cell is land in IMAGE2.2-landmask
	// but water in FAO soil map of the world
	// (data source for slope class)
    for (n = 0; n <= ng - 1; n++)
		if (0 == G_slope_class[n])
			G_slope_class[n] = 10;

    Grid<float> G_slope_factor;
	short slope_class_table[7] = { 10, 20, 30, 40, 50, 60, 70 };
	float slope_factor_table[7] = { 1.00, 0.95, 0.90, 0.75, 0.60, 0.30, 0.15 };

        for (n = 0; n <= ng - 1; n++) {
            // Enabling treatment of maximum slope class in G_SLOPE_CLASS.UNF1 (currently in file max. 69)
            // Maximum class index 6
            if (slope_class_table[6] == G_slope_class[n]) {
              G_slope_factor[n] = slope_factor_table[6];
            }
            // Class indices 0 - 5
            else {
                // Loop over class indices 0 - 5
                for (short i = 0; i <= 5; i++) {
                        // plain classes (implicit conversion/cast of char to short)
                        if (slope_class_table[i] == G_slope_class[n]) {
                                G_slope_factor[n] = slope_factor_table[i];
                                break;
                        } else if ((G_slope_class[n] > slope_class_table[i])  // interpolation with next class
                                           && (G_slope_class[n] < slope_class_table[i + 1])) {
                                // linear interpolation of values in slope_table
                                G_slope_factor[n] = slope_factor_table[i]
                                        + (((slope_factor_table[i + 1] - slope_factor_table[i])
                                                / (slope_class_table[i + 1] - slope_class_table[i]))
                                           * (G_slope_class[n] - slope_class_table[i]));
                                break;
                        }
                        // When breaks are not effective, the class boundaries are not respected
                        if (6 == i) {
                                cout << "ERROR groundwaterFactorClass::createGrids" << endl;
                                cout << "Class boundaries outside valid range in SLOPE-INPUT-FILE: " << inputDirectory << "/G_SLOPE_CLASS.UNF1" << endl;
                                cout << "at n: " << n << ", G_slope_class[n]: " << G_slope_class[n] << endl;
                                cerr << "ERROR groundwaterFactorClass::createGrids" << endl;
                                cerr << "Class boundaries outside valid range in SLOPE-INPUT-FILE: " << inputDirectory << "/G_SLOPE_CLASS.UNF1" << endl;
                                cerr << "at n: " << n << ", G_slope_class[n]: " << G_slope_class[n] << endl;
                        }
                }
            }
        }

        // Writing G_slope_factor to file
        if (options.outRGmax) {
            cout << "writing slope factors to file " << outputDirectory << "/G_SLOPE_FACTOR.UNF0" << endl;
	    G_slope_factor.write(outputDirectory + "/G_SLOPE_FACTOR.UNF0");
        }

	// TEXTURE FACTOR
	G_texture.read(inputDirectory + "/G_TEXTURE.UNF1");

	// values in file:
	// -2: landcell in IMAGE2.2 landmask but not
	//     in FAO soil map of the world
	// -1: land in IMAGE2.2 and FAO soil map, but no value available
	//  0: lakes, cells not considered (n>ng_land)
	//  1: cells with 100% galciers or rock
	//  2: cells with 100% salt lakes
	// 10-30: normal range of values

	Grid<float> G_texture_factor;
	Grid<float> G_Rgmax_f;

	// as 'float'. will be copied to 'short' array at the end.
	short texture_table[3] = { 10, 20, 30 };
	float Rgmax_table[3] = { 5., 3., 1.5 }; // standard value for Rgmax
	float Rgmax_tableWFD[3] = { 7., 4.5, 2.5 }; // values for WFD daily 

	if (0 == options.time_series){ // when WFD daily input is used, use modified Rgmax (BSc thesis Tögel)
		for (int i = 0; i < 3; i++)
			Rgmax_table[i] = Rgmax_tableWFD[i];
	}
	
	float texture_factor_table[3] = { 1, 0.95, 0.70 };
	for (n = 0; n <= ng - 1; n++) {

		// Efficiency enhancement through local copies instead of often used method getValue // FP
		double M_RG_MAX__at__n = calParam.getValue(M_RG_MAX,n);

		if ((G_texture[n] <= 0) || (2 == G_texture[n])) {
			// for those cases we assume that texture value is 20
			G_texture_factor[n] = 0.95;
			G_Rgmax_f[n] = 3;
			continue;
		}
		if (1 == G_texture[n]) {
			G_texture_factor[n] = 0;
			G_Rgmax_f[n] = 0;
			continue;
		}
		for (short i = 0; i <= 2; i++) {
            // plain classes (implicit conversion/cast of signed char to short)
			if (texture_table[i] == G_texture[n]) {
				G_texture_factor[n] = texture_factor_table[i];
				// Cell-specific calibration parameters - Apply multiplier  // FP
				G_Rgmax_f[n] = M_RG_MAX__at__n * Rgmax_table[i];
				break;
			} else if ((G_texture[n] > texture_table[i])
					   && (G_texture[n] < texture_table[i + 1])) {
				// linear interpolation of values in texture_table
				G_texture_factor[n] = texture_factor_table[i]
					+ (((texture_factor_table[i + 1] - texture_factor_table[i])
						/ (texture_table[i + 1] - texture_table[i]))
					   * (G_texture[n] - texture_table[i]));
				// Cell-specific calibration parameters - Apply multiplier  // FP
				G_Rgmax_f[n] = ( M_RG_MAX__at__n * Rgmax_table[i] )
					+ (
						( M_RG_MAX__at__n *
							(Rgmax_table[i + 1] - Rgmax_table[i])
								/ (texture_table[i + 1] - texture_table[i])
						)
					   * (G_texture[n] - texture_table[i])
					  );
				break;
			}

            // Corrected error treatment (value outside class limits) for G_TEXTURE.UNF1
            if (3 == i) {
                cout << "ERROR groundwaterFactorClass::createGrids" << endl;
                cout << "Class boundaries outside valid range in TEXTURE-INPUT-FILE: " << inputDirectory << "/G_TEXTURE.UNF1" << endl;
                cout << "at n: " << n << ", G_texture[n]: " << G_texture[n] << endl;
                cerr << "ERROR groundwaterFactorClass::createGrids" << endl;
                cerr << "Class boundaries outside valid range in TEXTURE-INPUT-FILE: " << inputDirectory << "/G_TEXTURE.UNF1" << endl;
                cerr << "at n: " << n << ", G_texture[n]: " << G_texture[n] << endl;
			}
		}
	}

    Grid<char> G_permaglac;	// permafrost and glacier [%]
	G_permaglac.read(inputDirectory + "/G_PERMAGLAC.UNF1");
	// related factor f_pg is: 1 - (G_permaglac[n]/100.0)

    Grid<char> G_aquifer_factor;	// aquifer factor [0-100]
	G_aquifer_factor.read(inputDirectory + "/G_AQ_FACTOR.UNF1");

	float slopeFactor, textureFactor, aquiferFactor;

	for (n = 0; n <= ng - 1; n++) {
		if (G_texture_factor[n] < 0)
			G_gwFactor[n] = -99;
		else {
			slopeFactor = G_slope_factor[n];
			textureFactor = G_texture_factor[n];
			aquiferFactor = (short) G_aquifer_factor[n] / 100.0;
			G_permaFactor[n] = 1. - (((float)G_permaglac[n] / 100.));

			// limit values to the range 0 - 1
			if (slopeFactor < 0)
				slopeFactor = 0;
			if (slopeFactor > 1)
				slopeFactor = 1;
			if (aquiferFactor < 0)
				aquiferFactor = 0;
			if (aquiferFactor > 1)
				aquiferFactor = 1;
			G_aquiferFactor[n] = aquiferFactor;
			if (textureFactor < 0)
				textureFactor = 0;
			if (textureFactor > 1)
				textureFactor = 1;
			G_textureFactor[n] = textureFactor;
			if (G_permaFactor[n] < 0)
				G_permaFactor[n] = 0;
			if (G_permaFactor[n] > 1)
				G_permaFactor[n] = 1;

			// Cell-specific calibration parameters - Apply multiplier  // FP
			G_gwFactor[n] = calParam.getValue(M_GW_F,n) * slopeFactor * textureFactor * aquiferFactor * G_permaFactor[n];

			// gw factor is set to 0.95 for cells with value > 1
			if (G_gwFactor[n] > 1.)
				G_gwFactor[n] = 0.95;


		}
	}

	for (n = 0; n < ng; n++) {
		if (G_Rgmax_f[n] >= 0)
			G_Rgmax[n] = round_to_short(G_Rgmax_f[n] * 100);
		else
			G_Rgmax[n] = -9999;
	}

    	// G_GW_FACT_CORR.UNF0 introduced in version WG2.2b (see documentation of WG2.2b)
	// read file G_GW_FACT_CORR.UNF0 including corrected GW_FACTOR in the Mississippi Embayment Regional Aquifer
	// Mississippi Embayment Regional Aquifer: fg = 0.1, all other cells: fg = -99
	G_gwFactorCorr.read(string(options.input_dir) + "/G_GW_FACTOR_CORR.UNF0");
	
	// write output
	if (options.outRGmax) {
		G_Rgmax.write(outputDirectory + "/G_RG_max.UNF2");
		G_textureFactor.write(outputDirectory + "/G_TEXTURE_FACTOR.UNF0");
		G_aquiferFactor.write(outputDirectory + "/G_AQUIFER_FACTOR.UNF0");
		G_permaFactor.write(outputDirectory + "/G_PERMA_FACTOR.UNF0");
	}
	
	if (options.outGWFactor) {
		// output file without correction for the Mississippi Embayment Regional Aquifer)
		G_gwFactor.write(outputDirectory + "/G_GW_FACTOR_UNCORR.UNF0");
	}
	
	// Correction of G_GW_FACTOR in the Mississippi Embayment Regional Aquifer
	// G_GW_FACTOR with corrected values is used in daily.cpp
	for (n = 0; n <= ng - 1; n++) {
		if (G_gwFactorCorr[n] > 0.) {

			//all other cells retain the G_gwFactor values computed above
			// Cell-specific calibration parameters - Apply multiplier
			// As the regional correction is considered to be permanent, the multiplier has to be applied also for these values
			G_gwFactor[n] = calParam.getValue(M_GW_F,n) * G_gwFactorCorr[n];

			// gw factor is set to 0.95 for cells with value > 1
			if (G_gwFactor[n] > 1.)
				G_gwFactor[n] = 0.95;
		}
	}

    	// write output (G_GW_FACTOR (corrected) is already written in permafrost.cpp)
    	if (options.outGWFactor) {
		// output file with correction for the Mississippi Embayment Regional Aquifer
		G_gwFactor.write(outputDirectory + "/G_GW_FACTOR_CORR.UNF0");
	}	
}

void groundwaterFactorClass::setgwFactor(int cell, float value)
{
        if (cell <0||cell>=ng) {
		cerr << "groundwaterFactorClass::setgwFactor cell number out of range!\n";
		exit(-1);
	}
	G_gwFactor[cell] = value;
}

void groundwaterFactorClass::setRgmax(int cell, short value)
{
        if (cell <0||cell>=ng) {
		cerr << "groundwaterFactorClass::setgwFactor cell number out of range!\n";
		exit(-1);
	}
	G_Rgmax[cell] = value;
}

void groundwaterFactorClass::setpermaFactor(int cell, float value)
{
        if (cell <0||cell>=ng) {
		cerr << "groundwaterFactorClass::setpermaFactor cell number outer range!\n";
		exit(-1);
	}
	G_permaFactor[cell] = value;
}

float groundwaterFactorClass::getgwFactor(int cell)
{
	return G_gwFactor[cell];
}

short groundwaterFactorClass::getRgmax(int cell)
{
	return G_Rgmax[cell];
}

float groundwaterFactorClass::getpermaFactor(int cell)
{
	return G_permaFactor[cell];
}

// these lines can be used to run the routine independent of WaterGAP
//void main(void) {
//  short G_Rgmax[ng_land];
//  float G_gw_factor[ng_land];
//
//  groundwater_fraction("input", "OUTPUT", G_Rgmax, G_gw_factor);
//}
