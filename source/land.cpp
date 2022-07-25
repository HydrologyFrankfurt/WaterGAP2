
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
#include <iostream>
#include <fstream>
#include <cstdio>

#include "gridio.h"
#include "def.h"

#include "land.h"
using namespace std;

void landClass::init(const char *input_dir)
{
	extern gridioClass gridIO;
	char filename[250];

	// built up area e.g. cities, fraction of cell area (0-1)
	sprintf(filename, "%s/GBUILTUP.UNF0", input_dir);
        //gridIO.readUnfFile(filename, ng_land, G_built_up);
        gridIO.readUnfFile(filename, ng, G_built_up);

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

	//sprintf(filename, "%s/GPOTVEG_PROX.UNF1", input_dir);
	//gridIO.readUnfFile(filename, ng, &G_potVegetation[0]);
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

// not used any more in WG2.2	
//	sprintf(filename, "%s/GLCT_1970.UNF1", input_dir);
//	gridIO.readUnfFile(filename, ng_land, G_landCover22);
// not used any more in WG2.2	
	/* WaterGAP2.2 */
    sprintf(filename, "%s/G_LANDCOVER.UNF1", input_dir); // Renamed landcover input (as GLCC is not any longer source but MODIS)
        //gridIO.readUnfFile(filename, ng_land, G_landCover22);
        gridIO.readUnfFile(filename, ng, G_landCover22);

	// if the cell is a land cell, then landcover is used
	// otherwise potential vegetation from IMAGE 2.2 is used
	int n;

        //for (n = 0; n < ng_land; n++)
        for (n = 0; n < ng; n++)
		G_landCover[n] = G_landCover22[n];

	// add class "water bodies" to lake cells
        //for (n = ng_land; n < ng; n++)
                //G_landCover[n] = 17; // add class "water bodies" to lake cells
// not used any more in WG2.2	
	// add pot. vegatation to lake cells
	//for (n = ng_land; n < ng; n++)
	//	G_landCover[n] = G_potVegetation[n];
// not used any more in WG2.2	
}

void landClass::readLandCoverMap(int actual_year, const char *land_cover_dir)
{
	extern gridioClass gridIO;
	char filename[250];

	sprintf(filename, "%s/GLCT_%d.UNF1", land_cover_dir, actual_year);
        //gridIO.readUnfFile (filename, ng_land, G_landCover22);
        gridIO.readUnfFile (filename, ng, G_landCover22);

	long n;
	// add land use values from GLCT to land cells
        //for (n = 0; n < ng_land; n++)
        for (n = 0; n < ng; n++)
		G_landCover[n] = G_landCover22[n]; // WaterGAP3.1

	// add class "water bodies" to lake cells
        //for (n = ng_land; n < ng; n++)
        //	G_landCover[n] = 17; // add class "water bodies" to lake cells

}

void landClass::initPM(const char *input_dir)
{
	extern gridioClass gridIO;
	char filename[250];

	G_altitude = new short[ng];

	sprintf(filename, "%s/G_ALTITUDE.UNF2", input_dir);
	gridIO.readUnfFile(filename, ng, G_altitude);
}
