
/***********************************************************************
*
*New land cover map and new covert type classification (IGBP) by CESR, sources
*described in Verzano, K. (2009), Climate Change Impacts on Flood Related
*Hydrological Processes: Further Development and Application of a Global Scale
*Hydrological Model, 166 pp., University of Kassel.
*
* see former changes at file land.cpp.versioninfos.txt
*
***********************************************************************/
//#include <iostream>
//#include <fstream>
#include <cstdio>
#include "common.h"
#include "def.h"
#include "globals.h"
//#include "land.h"
using namespace std;

void landClass::init(const std::string input_dir)
{
	// built up area e.g. cities, fraction of cell area (0-1)
	G_built_up.read(input_dir + "/GBUILTUP.UNF0");

	// not used any more in WG2.2
	//
	//  read potential vegetation        
	//                                   
	// 6: Ice                            
	// 7: Tundra                         
	// 8: Wooded tundra                  
	// 9: Boreal forest                  
	// 10: Cool conifer forest            
	// 11: Temperate mixed forest         
	// 12: Temperate deciduous forest     
	// 13: Warm mixed forest             
	// 14: Steppe                        
	// 15: Hot desert                    
	// 16: Scrubland                     
	// 17: Savanna                       
	// 18: Tropical woodland             
	// 19: Tropical forest                 

	// the GPOTVEG_PROX file also contains vegetation for the lake areas
	// because we need it to decide whether the area is arid or humid
	// in the evapotransiration calculations
	// this file has been created with an ARC/VIEW function by just looking
	// at the value of the closest cell
	// the file is based on the IMAGE file GPOTVEG_1970.UNF1, which only contains
	// information for land areas

	//  read land cover types                       
	//                                              
	// 1: Agricultural land: Mainly cropland        
	// 2: Agricultural land: Mixed cropland/pasture 
	// 3: Agricultural land: Mainly pasture         
	// 4: Agricultural land: Marginal               
	// 5: Regrowth forest                           
	// 6: Ice                                       
	// 7: Tundra                                    
	// 8: Wooded tundra                             
	// 9: Boreal Forest                             
	// 10: Cool conifer forest                      
	// 11: Temp. mixed forest                       
	// 12: Temp. deciduous forest                   
	// 13: Warm mixed forest                        
	// 14: Grassland/Steppe                         
	// 15: Hot desert                               
	// 16: Scrubland                                
	// 17: Savanna                                  
	// 18: Tropical woodland                        
	// 19: Tropical forest                          
	// 20: Major cities                             
// not used any more in WG2.2	

	//  read land cover types (IGBP land cover legend + CLC (re-)classification) 
	//
	// 1: Evergreen needle leaf forest
	// 2: Evergreen broadleaf forest
	// 3: Deciduous needle leaf forest
	// 4: Deciduous broadleaf forest
	// 5: Mixed forest
	// 6: Closed Shrubland
	// 7: Open Shrubland
	// 8: Woody Savanna
	// 9: Savanna
	// 10: Grassland
	// 11: Permanent Wetland
	// 12: Cropland
	// 13: Urban and built up
	// 14: Cropland/ natural vegetation mosaik 
	// !!!new: 141: Cropland / permanent crops (CLC); implemented as class 18 ::new!!!
	// 15: Snow and Ice
	// 16: Barren or sparsely vegetated
	// 17: Water bodies (for WaterGAP lake cells only)

	/* WaterGAP2.2 */
	// Renamed landcover input (as GLCC is not any longer source but MODIS)
	G_landCover22.read(input_dir +"/G_LANDCOVER.UNF1");

	// if the cell is a land cell, then landcover is used
	// otherwise potential vegetation from IMAGE 2.2 is used
    for (int n = 0; n < ng; n++)
		G_landCover[n] = G_landCover22[n];
}

void landClass::readLandCoverMap(int actual_year, const std::string land_cover_dir)
{
	G_landCover22.read(string(land_cover_dir) + "/GLCT_" + to_string(actual_year) + ".UNF1");

    for (long n = 0; n < ng; n++)
		G_landCover[n] = G_landCover22[n];
}

void landClass::initPM(const std::string input_dir)
{
	G_altitude.read(input_dir + "/G_ALTITUDE.UNF2");
}
