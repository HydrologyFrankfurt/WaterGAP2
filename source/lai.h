
/***********************************************************************
*
* see former changes at file lai.h.versioninfos.txt
*
***********************************************************************/
#include "def.h"

class laiClass {

  private:
	// variables used in 'getDailyLAI'
	// values of these variables are calculated in 'init'
	float G_LAImax[ng];
	char G_days_since_start[ng]; 	// counter for days with growing conditions
	char G_GrowingStatus[ng]; 		// defines the status of vegetation development
	double G_PrecSum[ng]; 			// threshold for min precipitation sum to start growing season
	float lai_factor_a[nlct];
	float lai_factor_b[nlct];
	double dailyLai_growing(char *days_since_start, char initialDays, char *GrowingStatus, char land_cover_type, bool arid, double LAImin, double LAImax, double *PrecSum, double prec);
	double dailyLai_nogrowing(char *days_since_start, char initialDays, char *GrowingStatus, char land_cover_type, bool arid, double LAImin, double LAImax, double *PrecSum, double prec);

	void readLAIdata(const char *inputDir);
	bool LAIdataIsRead;
	// variables which are read from 'LAI.DAT'
// not used any more in WG2.2	
//	float specificLeafAreaOrig[nlct];	// [m2/kg dry mass]
//	float leafMassCorrOrig[nlct];	// [-]
// not used any more in WG2.2	
	float LeafAreaIndexOrig[nlct];			// [-]			WaterGAP2.2
	float decidiousPlantFractionOrig[nlct];	// [-]
	float evergreensLAIreductionOrig[nlct];	// [-]
	short initialDays[nlct];				// initial days to start/end growing season
	double kc_min[nlct];					// minimum Kc kc_ini for calculation of veg. specific pet
	double kc_max[nlct];					// maximum Kc kc_mid for calculation of veg. specific pet

	// variables for SA
	double LAImultiplier;
	double fdpError;	// fraction of decidious plants
	double rfepError;	// reduction factor for evergreen plants

  public:
	 laiClass(void);
	void init(const char *inputDir, char *output_dir, const char *G_land_cover);
	// use threshold of precipitation sum instead of dailyPET for decision of growing conditions
	//double getDailyLai(int n, char land_cover_type, bool arid_humid, double daily_temp_C, double precipitation, double dailyPET);
	double getDailyLai(int n, char land_cover_type, bool arid, double daily_temp_C, double precipitation);
	double getDailyKc(int n, char land_cover_type, double dailyLAI);
	void setSAParameters();

};
