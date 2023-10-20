
/***********************************************************************
*
* see former changes at file calib_basins.h.versioninfos.txt
*
***********************************************************************/
#if !defined (_calib_basins_h_)
#define _calib_basins_h_

#include <vector>
#include <string>
#include "grid.h"

using namespace std;


class cbasinClass {
  public:
	vector<float> calibArea;
	vector<float> calibLandArea;
	vector<string> name;
	vector<int> cellNum;
	vector<float> gamma;
	vector<float> cellCorrFactor;
	vector<float> statCorrFactor;
	short numberOfCalibBasins;

	short prepare(string input_dir, string output_dir, string routing_dir);

	short getStationNumber(int cellNumber);
};
#endif
