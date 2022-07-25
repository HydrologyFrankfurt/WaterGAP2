
/***********************************************************************
*
* see former changes at file land.h.versioninfos.txt
*
***********************************************************************/
#include "def.h"

class landClass {
  public:
	void init(const char *input_dir);
	void initPM(const char *input_dir);
	void readLandCoverMap(int actual_year, const char *land_cover_dir);

        //float G_built_up[ng_land];
        float G_built_up[ng];
// not used any more in WG2.2	
	char G_potVegetation[ng];
// not used any more in WG2.2	
	char G_landCover[ng];
        //char G_landCover22[ng_land];
        char G_landCover22[ng];

	short *G_altitude;	// only for Penman-Monteith (option: 2)
};
