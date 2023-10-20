#if !defined (_s_max_h_)
#define _s_max_h_
/***********************************************************************
*
* see former changes at file s_max.h.versioninfos.txt
*
***********************************************************************/
#include "def.h"
#include "grid.h"
#include <string>

class soilWatCapClass {
  public:
		soilWatCapClass(void);
		Grid<float> G_Smax;	// total available soil water capacity
		void createMaxSoilWaterCapacityGrid(std::string input_dir, std::string output_dir, Grid<char> G_land_cover,calibParamClass calParam);//OE: calParam added);
		void setParameters();

  private:
    float rootingDepth[nlct];
		bool parameterFlag;
};
#endif
