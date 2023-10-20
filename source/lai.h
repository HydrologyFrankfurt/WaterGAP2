
/***********************************************************************
*
* see former changes at file lai.h.versioninfos.txt
*
***********************************************************************/
#include "def.h"
#include "grid.h"
#include "additionalOutputInputFile.h"
#include <string>

class laiClass {

  private:
	// variables used in 'getDailyLAI'
	// values of these variables are calculated in 'init'
	Grid<float> G_LAImax;
	Grid<int> G_days_since_start; 	// counter for days with growing conditions
	Grid<int> G_GrowingStatus; 		// defines the status of vegetation development
	Grid<> G_PrecSum; 			// threshold for min precipitation sum to start growing season
	float lai_factor_a[nlct];
	float lai_factor_b[nlct];
	double dailyLai_growing(int *days_since_start, short initialDays, int *GrowingStatus, char land_cover_type, bool arid, double LAImin, double LAImax, double *PrecSum, double prec);
	double dailyLai_nogrowing(int *days_since_start, short initialDays, int *GrowingStatus, char land_cover_type, bool arid, double LAImin, double LAImax, double *PrecSum, double prec);

	void readLAIdata(const std::string inputDir);
	bool LAIdataIsRead;

	// variables which are read from 'LAI.DAT'
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
	void init(const std::string inputDir, const std::string output_dir, const Grid<char> G_land_cover,AdditionalOutputInputFile &additionalOutIn,calibParamClass &calParam);

	// use threshold of precipitation sum instead of dailyPET for decision of growing conditions
	double getDailyLai(int n, char land_cover_type, bool arid, double daily_temp_C, double precipitation, const short day, const short last_day_in_month, AdditionalOutputInputFile &additionalOutIn);
	double getDailyKc(int n, char land_cover_type, double dailyLAI);
	void setSAParameters();

};
