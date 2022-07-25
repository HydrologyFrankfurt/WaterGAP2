
/***********************************************************************
*
* see former changes at file daily.h.versioninfos.txt
*
***********************************************************************/
#if !defined (_daily_h_)
#define _daily_h_

#include <vector>
#include "def.h"
using namespace std;

class dailyWaterBalanceClass {
  public:
	void calcNewDay(const short day,
					const short month, const short day_in_month, const short year, const int n);
	void init(const char *inputDir, const short nBasins);
	void setParameters();
    void annualInit();
    void daily31outInit();
        void daily365outInit();

	void setStoragesToZero();
	// Superbasin calculation currently obsolete, probably not working correctly in WaterGAP2.2b, do in postprocessing  // FP 2015-10
//	void storeDailyValuesToFile(const short year);
//	void storeAnnualValuesToFile(const short year);
	void writeAnnualGrids(const char *output_dir, const short year);
	void writeMonthlyGrids(const char *output_dir, const short year);
        void writeDaily31Grids(char *outputDir, const int year, const int month);
        void writeDaily365Grids(char *outputDir, const int year);
	// Superbasin calculation currently obsolete, probably not working correctly in WaterGAP2.2b, do in postprocessing  // FP 2015-10
	// Unused methods get...Vector
//	void getPETVector(vector < float >&Vector);
//	void getAETVector(vector < float >&Vector);
//	void getInterceptVector(vector < float >&Vector);
//	void getRunoffVector(vector < float >&Vector);
//	void getSoilSatVector(vector < float >&Vector);
//	void getSnowVector(vector < float >&Vector);
    double get_rootingDepth_lct(short lct) {return rootingDepth_lct[lct];}
    void set_rootingDepth_lct(short lct, float value) {rootingDepth_lct[lct]=value;}
	//void readElevationFile();
	//void readDailyTempVarFile(); // not used in 2.1f
	dailyWaterBalanceClass();
	~dailyWaterBalanceClass();
	
	// added for cell AET calculation (2.1f)
	double getDailyCellLandAET(const int cellNumber); 

        double G_gammaHBV[ng];	// calibration factor
        double G_cellCorrFact[ng];	// local correction factor

	// parameters for improved snow algorithm
        //short int G_Elevation[ng_land][101];	// [meters] ImageID (1-66663), Number of elevation >0, mean elevation, 100 x elevation
        short int G_Elevation[ng][101];	// [meters] ImageID (1-66869), Number of elevation >0, mean elevation, 100 x elevation
        double G_snow[ng];	// [mm]


	// actual water content down to the rooting depth [mm]
        double G_soilWaterContent[ng];

    // HMS 2014-06-04 reintroduced
    // daily values of precipitation minus PET
    // for cells with lakes
    // required for the routing component
    double G_lakeBalance[ng];

	// added for open water evaporation reduction (2.1f)
        double G_openWaterPET[ng];
        double G_openWaterPrec[ng];



        double G_groundwater[ng];
	// G_dailyLocalRunoff[ng] will be substituted by...
        double G_dailyLocalSurfaceRunoff[ng];
        double G_dailyLocalGWRunoff[ng];
		// CR 2015-08: G_dailyStorageTransfer introduced in 22b
		// If land area fraction becomes zero, soil, snow, and canopy storages from the last time step
		// are transferred to routing.cpp using G_dailyStorageTransfer[n] (in mm).
		// In routing.cpp, this value is multiplied by G_landAreaFracPrevTimestep[n]/100*cell_area.
		double G_dailyStorageTransfer[ng];

    // array for groundwater_recharge
    double G_dailyGwRecharge[ng];
    double (*G_monthlyGwRecharge)[12];
	// ... to get seperate values for surface water and sub-surface water reaching the river network
        double G_canopyWaterContent[ng];
    short getSingleCellNumber(int cellNumber); // HMS 2013-11-21 Single cell output
    //daily storages now public as they are used in routing.cpp
    double (*G_daily31CanopyStorage)[31];
    double (*G_daily31SnowStorage)[31];
    double (*G_daily31SoilStorage)[31];
    double (*G_daily2CanopyStorage)[2];
    double (*G_daily2SnowStorage)[2];
    double (*G_daily2SoilStorage)[2];
    double (*G_daily365CanopyStorage)[365];
    double (*G_daily365SnowStorage)[365];
    double (*G_daily365SoilStorage)[365];
    double (*G_daily365LAI)[365];
    double (*G_daily365Kc)[365];
    double (*G_daily365Albedo)[365];
    double (*G_daily365TempC)[365];  // degree celsius
    double (*G_daily365ExtRad)[365]; // W/m2
    double (*G_daily365ShortDownRad)[365]; // W/m2
    double (*G_daily365ShortUpRad)[365]; // W/m2
    double (*G_daily365LongDownRad)[365]; // W/m2
    double (*G_daily365LongUpRad)[365]; // W/m2
    double (*G_daily365NetLongRad)[365]; // W/m2
    double (*G_daily365NetShortRad)[365]; // W/m2
    double (*G_daily365NetRadiation)[365]; // W/m2
    double (*G_daily365Interception)[365];
    double (*G_daily365canopyWaterContent)[365];
    double (*G_daily365maxcanopyWaterContent)[365];
    double (*G_daily365Precip)[365];
    double (*G_daily365PET)[365];
    double (*G_daily365TotalPET)[365];
    double (*G_daily365SurfaceRunoff)[365];

    double (*G_daily365GwRecharge)[365];
    double (*G_daily365GwStorage)[365];
    double (*G_daily365SnowCoverAvg)[365];
    double (*G_daily365SWE)[365];
    double (*G_daily365SnowFall)[365];
    double (*G_daily365SnowMelt)[365];
    double (*G_daily365SoilWaterAvg)[365];

    // FP20161018N002 Reservoir operation start years
    float G_SnowInElevation[ng][101];	// [mm] array to store snow cover for 100 seperate elevations within G_snow cell

  private:
	// Superbasin calculation currently obsolete, probably not working correctly in WaterGAP2.2b, do in postprocessing  // FP 2015-10
//	void normalizeDailyResults();
	// indicator: are daily values already devided by basin area?
//	bool normalizedFlag;
	// number of specified basins
//	short nSpecBasins;

	double calc_ext_rad(const short day, const int n);

	// the following grids are only needed to store information for output files
	// they can be removed from the code without any influence on the simulation
	// (if we run into memory problems)
	// output binary files for monthly/yearly values
        double (*G_monthlyPrecipitation)[12];
        double (*G_monthlyTemperature)[12];
        double (*G_monthlyNetRadiation)[12];	// W/m2
        double (*G_monthlyNetShortRad)[12]; // W/m2
        double (*G_monthlyNetLongRad)[12]; // W/m2
        double (*G_monthlyShortDownRad)[12];	// W/m2
        double (*G_monthlyExtRad)[12];	// W/m2
        double (*G_monthlyShortUpRad)[12];	// W/m2
        double (*G_monthlyLongDownRad)[12];	// W/m2
        double (*G_monthlyLongUpRad)[12];	// W/m2
        double (*G_monthlySunshine)[12];	// sunshine percentage output
        double (*G_monthlyLAI)[12];
        double (*G_monthlyAlbedo)[12];
        double (*G_monthlyInterception)[12];
        double (*G_monthlycanopyWaterContent)[12];
        double (*G_monthlymaxcanopyWaterContent)[12];
        double (*G_monthlyAET)[12];
	// has been moved from routing.h
        double (*G_monthlyLandAET)[12]; // [mm]
        double (*G_monthlyLandAET_uncorr)[12]; // [mm]
        double (*G_monthlyPET)[12];
        double (*G_monthlyOpenWaterPET)[12];
        double (*G_monthlyTotalPET)[12];
        double (*G_monthlyRunoff)[12];
        double (*G_monthlyUrbanRunoff)[12];
        double (*G_monthlySurfaceRunoff)[12];
	//float (*G_monthlyGwRunoff)[12];
        //double (*G_monthlyGwRecharge)[12]; // now public to access in routing
        double (*G_monthlySnowCoverAvg)[12];
        double (*G_monthlySWE)[12];
        double (*G_monthlySnowFall)[12];
        double (*G_monthlyRainFall)[12]; // FP20161011N001
        double (*G_monthlySnowMelt)[12];
        double (*G_monthlySnowEvap)[12];
        double (*G_monthlySoilWaterAvg)[12];

	//monthly storages
	double (*G_monthlyCanopyStorage)[12];
	double (*G_monthlySnowStorage)[12];
	double (*G_monthlySoilStorage)[12];
	double (*G_monthly_totalWaterInStorages)[12];
	double (*G_monthly_totalWaterInStorages_mm)[12];
        double G_monthlyArray[ng][12]; // temporary array used to convert double to float values

  // daily output option 31 start
        double (*G_daily31LAI)[31];
        double (*G_daily31Kc)[31];
        double (*G_daily31Albedo)[31];
        double (*G_daily31ExtRad)[31]; // W/m2
    double (*G_daily31ShortDownRad)[31]; // W/m2
    double (*G_daily31ShortUpRad)[31]; // W/m2
    double (*G_daily31LongDownRad)[31]; // W/m2
    double (*G_daily31LongUpRad)[31]; // W/m2
    double (*G_daily31NetLongRad)[31]; // W/m2
    double (*G_daily31NetShortRad)[31]; // W/m2
    double (*G_daily31NetRadiation)[31]; // W/m2
        double (*G_daily31Interception)[31];
        double (*G_daily31canopyWaterContent)[31];
        double (*G_daily31maxcanopyWaterContent)[31];
    double (*G_daily31Precip)[31];
        double (*G_daily31PET)[31];
        double (*G_daily31TotalPET)[31];
        double (*G_daily31SurfaceRunoff)[31];
        double (*G_daily31GwRunoff)[31];
        double (*G_daily31GwRecharge)[31];
        double (*G_daily31GwStorage)[31];
        double (*G_daily31SnowCoverAvg)[31];
        double (*G_daily31SWE)[31];
        double (*G_daily31SnowFall)[31];
        double (*G_daily31SnowMelt)[31];
        double (*G_daily31SoilWaterAvg)[31];
  // daily output option 31 end
  
/*  // daily output option 365 start // now as public to get access from routing.cpp
	float (*G_daily365LAI)[365];
	float (*G_daily365Kc)[365];
	float (*G_daily365Albedo)[365];
  	float (*G_daily365TempC)[365];  // degree celsius
	float (*G_daily365ExtRad)[365]; // W/m2
  	float (*G_daily365GloRad)[365]; // W/m2
  	float (*G_daily365ReflRad)[365]; // W/m2
  	float (*G_daily365InLongRad)[365]; // W/m2
  	float (*G_daily365OutLongRad)[365]; // W/m2
  	float (*G_daily365NetLongRad)[365]; // W/m2
  	float (*G_daily365NetShortRad)[365]; // W/m2
  	float (*G_daily365NetRadiation)[365]; // W/m2	
	float (*G_daily365Interception)[365];
	float (*G_daily365canopyWaterContent)[365];
	float (*G_daily365maxcanopyWaterContent)[365];
    float (*G_daily365Precip)[365];
	float (*G_daily365PET)[365];
	float (*G_daily365TotalPET)[365];
	float (*G_daily365SurfaceRunoff)[365];
	float (*G_daily365GwRunoff)[365];
	float (*G_daily365GwRecharge)[365];
	float (*G_daily365GwStorage)[365];
	float (*G_daily365SnowCoverAvg)[365];
	float (*G_daily365SWE)[365];
	float (*G_daily365SnowFall)[365];
	float (*G_daily365SnowMelt)[365];
	float (*G_daily365SoilWaterAvg)[365];
  // daily output option 365 end */
	// Superbasin calculation currently obsolete, probably not working correctly in WaterGAP2.2b, do in postprocessing  // FP 2015-10
	// arrays to store daily values for superbasins
//        double (*dailyPrecSb)[365];
//        double (*dailyEffPrecSb)[365];
//        double (*dailySnowSb)[365];
//        double (*dailyAET_Sb)[365];
//        double (*dailyPET_Sb)[365];
//        double (*dailyRunoffSb)[365];
//        double (*dailySaturationSb)[365];
//        double (*dailySunshineSb)[365];
//        double (*dailyLaiSb)[365];
//        double (*dailyCanopyWaterSb)[365];
//        double (*dailyCanopyEvapoSb)[365];
//        double (*dailyTempSb)[365];
//        double *snowDaysSb;
//        double *precAsSnowSb;

//	double alphaArid;  // OBSOLETE now read from calibration parameter file // FP
//	double alphaHumid;  // OBSOLETE now read from calibration parameter file // FP
//	double maxCanopyStoragePerLAI;  // OBSOLETE now read from calibration parameter file // FP
	double canopyEvapoExp;
	double runoffFracBuiltUp;
//	double maxDailyPET_humid;  // OBSOLETE now read from calibration parameter file // FP
//	double maxDailyPET_arid;  // OBSOLETE now read from calibration parameter file // FP

	
	// variables which are read from 'LCT_22.DAT'
	void readLCTdata(const char *inputDir);
	bool LCTdataIsRead;
	float rootingDepth_lct[nlct];
	double albedoSnow_lct[nlct];
	double albedo_lct[nlct];
	double ddf_lct[nlct];
        double emissivity_lct[nlct]; // land use dependent emissivity

	// albedo for open water bodies and snow;
	double openWaterAlbedo;
	// pet-correction for open water bodies
	double kc_OpenWater;
	
	// added for innerannual variation in albedo
	//float (*G_AlbedoSoil);
	//float (*G_AlbedoCanopy);
	// use given albedo in the first step
	float (*G_AlbedoCycle)[12];
	
//	double snowMeltTemp;  // OBSOLETE now read from calibration parameter file // FP
//	double snowFreezeTemp;  // OBSOLETE now read from calibration parameter file // FP
	double a_s, b_s;
	double a_c_arid, b_c_arid, a_c_humid, b_c_humid;	// cloudiness coefficients


	// parameters for improved snow algorithm
        //float G_SnowInElevation[ng_land][101];	// [mm] array to store snow cover for 100 seperate elevations within G_snow cell
        // FP20161018N002 Reservoir operation start years: now public G_SnowInElevation
        //float G_SnowInElevation[ng][101];	// [mm] array to store snow cover for 100 seperate elevations within G_snow cell
        //int thresh_elev[ng_land];	// threshold elevation as level for minimum temperature in high elevated cells
        int thresh_elev[ng];	// threshold elevation as level for minimum temperature in high elevated cells
	
	// added for cell AET calculation (2.1f)
	double dailyCellLandAET[ng];
	double dailyCellLandAET_uncorr[ng];


};
#endif
