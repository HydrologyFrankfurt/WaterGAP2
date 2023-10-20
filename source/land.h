
/***********************************************************************
*
* see former changes at file land.h.versioninfos.txt
*
***********************************************************************/
#include "def.h"
#include "grid.h"
#include <string>

class landClass {
  public:
	void init(const std::string input_dir);
	void initPM(const std::string input_dir);
	void readLandCoverMap(int actual_year, const std::string land_cover_dir);

        Grid<float> G_built_up;
	// not used any more in WG2.2

	Grid<char> G_potVegetation;
	// not used any more in WG2.2

	Grid<char> G_landCover;
        Grid<char> G_landCover22;

	Grid<short> G_altitude;	// only for Penman-Monteith (option: 2)
};
