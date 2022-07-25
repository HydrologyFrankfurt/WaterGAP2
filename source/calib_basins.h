
/***********************************************************************
*
* see former changes at file calib_basins.h.versioninfos.txt
*
***********************************************************************/
#if !defined (_calib_basins_h_)
#define _calib_basins_h_

class cbasinClass {
  public:
	float *calibArea;
	float *calibLandArea;
	char **name;
	int *cellNum;
	float *gamma;
	float *cellCorrFactor;
	float *statCorrFactor;
	short numberOfCalibBasins;

	short prepare(char *input_dir, char *output_dir, char *routing_dir);

	short getStationNumber(int cellNumber);
};
#endif
