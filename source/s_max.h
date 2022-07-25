
/***********************************************************************
*
* see former changes at file s_max.h.versioninfos.txt
*
***********************************************************************/
#include "def.h"

class soilWatCapClass {
  public:
	soilWatCapClass(void);
	float G_Smax[ng];	// total available soil water capacity
	void createMaxSoilWaterCapacityGrid(char *input_dir, char *output_dir, char *G_land_cover);
    void setParameters();


  private:
    float rootingDepth[nlct];
	bool parameterFlag;
  	// variables which is read from 'LCT_22.DAT' now in dailyWaterBalanceClass as rootingDepth_lct
	// float rootingDepth[nlct];
};
