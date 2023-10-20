
/***********************************************************************
*
* see former changes at file daily.h.versioninfos.txt
*
***********************************************************************/
#if !defined (_daily_h_)
#define _daily_h_

#include <vector>
#include "def.h"
#include <string>
// AEMS:
#include "wghmStateFile.h"
#include "additionalOutputInputFile.h"
#include "grid.h"
#include "snowInElevationFile.h"
#include "calib_param.h"

using namespace std;

class dailyWaterBalanceClass {
  public:
	void calcNewDay(const short day, const short month, const short day_in_month, const short last_day_in_month, const short year, const int n,  WghmStateFile &wghmState, AdditionalOutputInputFile &additionalOutIn, SnowInElevationFile &snow_in_elevation, const short readinstatus, calibParamClass &calParam);//OE: calParam added
	void init(const string inputDir, const short nBasins);
	void setParameters();
        void setStorages(WghmStateFile &wghmState, SnowInElevationFile &snow_in_elevation,AdditionalOutputInputFile &additionalOutIn); //AEMS: new function to set storages
        void annualInit();
        void daily31outInit();
        void daily365outInit();

	void setStoragesToZero();
	void writeAnnualGrids(const string output_dir, const short year);
	void writeMonthlyGrids(const string output_dir, const short year, const short month2store);
        void writeDaily31Grids(const string outputDir, const int year, const int month);
        void writeDaily365Grids(const string outputDir, const int year);
        double get_rootingDepth_lct(short lct) {return rootingDepth_lct[lct];}
        void set_rootingDepth_lct(short lct, double value) {rootingDepth_lct[lct]=value;}
	dailyWaterBalanceClass();
	~dailyWaterBalanceClass();
	
	// added for cell AET calculation (2.1f)
	double getDailyCellLandAET(const int cellNumber); 

        Grid<> G_gammaHBV;	// calibration factor
        Grid<> G_cellCorrFact;	// local correction factor

	// parameters for improved snow algorithm
        VariableChannelGrid<101,short int> G_Elevation;	// [meters] ImageID (1-66869), Number of elevation >0, mean elevation, 100 x elevation
        Grid<> G_snow;	// [mm]


	// actual water content down to the rooting depth [mm]
        Grid<> G_soilWaterContent;

        // daily values of precipitation minus PET
        // for cells with lakes
        // required for the routing component
        Grid<> G_lakeBalance;

	// added for open water evaporation reduction (2.1f)
        Grid<> G_openWaterPET;
        Grid<> G_openWaterPrec;

        Grid<> G_groundwater;

	// G_dailyLocalRunoff[ng] will be substituted by...
        Grid<> G_dailyLocalSurfaceRunoff;
        Grid<> G_dailyLocalGWRunoff;

        // If land area fraction becomes zero, soil, snow, and canopy storages from the last time step
        // are transferred to routing.cpp using G_dailyStorageTransfer[n] (in mm).
        // In routing.cpp, this value is multiplied by G_landAreaFracPrevTimestep[n]/100*cell_area.
        Grid<> G_dailyStorageTransfer;

        // array for groundwater_recharge
        Grid<> G_dailyGwRecharge;
        MonthlyGrid<false> G_monthlyGwRecharge;

        // ... to get seperate values for surface water and sub-surface water reaching the river network
        Grid<> G_canopyWaterContent;

        short getSingleCellNumber(int cellNumber);

        Daily31Grid<false> G_daily31CanopyStorage;
        Daily31Grid<false> G_daily31SnowStorage;
        Daily31Grid<false> G_daily31SoilStorage;
        VariableChannelGrid<2> G_daily2CanopyStorage;
        VariableChannelGrid<2> G_daily2SnowStorage;
        VariableChannelGrid<2> G_daily2SoilStorage;
        DailyGrid<false> G_daily365CanopyStorage;
        DailyGrid<false> G_daily365SnowStorage;
        DailyGrid<false> G_daily365SoilStorage;
        DailyGrid<false> G_daily365LAI;
        DailyGrid<false> G_daily365Kc;
        DailyGrid<false> G_daily365Albedo;
        DailyGrid<false> G_daily365TempC;  // degree celsius
        DailyGrid<false> G_daily365ExtRad; // W/m2
        DailyGrid<false> G_daily365ShortDownRad; // W/m2
        DailyGrid<false> G_daily365ShortUpRad; // W/m2
        DailyGrid<false> G_daily365LongDownRad; // W/m2
        DailyGrid<false> G_daily365LongUpRad; // W/m2
        DailyGrid<false> G_daily365NetLongRad; // W/m2
        DailyGrid<false> G_daily365NetShortRad; // W/m2
        DailyGrid<false> G_daily365NetRadiation; // W/m2
        DailyGrid<false> G_daily365Interception;
        DailyGrid<false> G_daily365canopyWaterContent;
        DailyGrid<false> G_daily365maxcanopyWaterContent;
        DailyGrid<false> G_daily365Precip;
        DailyGrid<false> G_daily365temp;
        DailyGrid<false> G_daily365redtemp;
        DailyGrid<false> G_daily365PET;
        DailyGrid<false> G_daily365TotalPET;
        DailyGrid<false> G_daily365SurfaceRunoff;
        DailyGrid<false> G_daily365GwRecharge;
        DailyGrid<false> G_daily365SnowCoverAvg;
        DailyGrid<false> G_daily365SWE;
        DailyGrid<false> G_daily365SnowFall;
        DailyGrid<false> G_daily365SnowMelt;
        DailyGrid<false> G_daily365SoilWaterAvg;

        VariableChannelGrid<101,double> G_SnowInElevation;	// [mm] array to store snow cover for 100 seperate elevations within G_snow cell

  private:
	double calc_ext_rad(const short day, const int n);

	// the following grids are only needed to store information for output files
	// they can be removed from the code without any influence on the simulation
	// (if we run into memory problems)
	// output binary files for monthly/yearly values
        MonthlyGrid<false> G_monthlyPrecipitation;
        MonthlyGrid<false> G_monthlyTemperature;
        MonthlyGrid<false> G_monthlyNetRadiation;	// W/m2
        MonthlyGrid<false> G_monthlyNetShortRad; // W/m2
        MonthlyGrid<false> G_monthlyNetLongRad; // W/m2
        MonthlyGrid<false> G_monthlyShortDownRad;	// W/m2
        MonthlyGrid<false> G_monthlyExtRad;	// W/m2
        MonthlyGrid<false> G_monthlyShortUpRad;	// W/m2
        MonthlyGrid<false> G_monthlyLongDownRad;	// W/m2
        MonthlyGrid<false> G_monthlyLongUpRad;	// W/m2
        MonthlyGrid<false> G_monthlySunshine;	// sunshine percentage output
        MonthlyGrid<false> G_monthlyLAI;
        MonthlyGrid<false> G_monthlyAlbedo;
        MonthlyGrid<false> G_monthlyInterception;
        MonthlyGrid<false> G_monthlycanopyWaterContent;
        MonthlyGrid<false> G_monthlymaxcanopyWaterContent;
        MonthlyGrid<false> G_monthlyAET;
        MonthlyGrid<false> G_monthlyLandAET; // [mm]
        MonthlyGrid<false> G_monthlyLandAET_uncorr; // [mm]
        MonthlyGrid<false> G_monthlyPET;
        MonthlyGrid<false> G_monthlyOpenWaterPET;
        MonthlyGrid<false> G_monthlyTotalPET;
        MonthlyGrid<false> G_monthlyRunoff;
        MonthlyGrid<false> G_monthlyUrbanRunoff;
        MonthlyGrid<false> G_monthlySurfaceRunoff;
        MonthlyGrid<false> G_monthlySnowCoverAvg;
        MonthlyGrid<false> G_monthlySWE;
        MonthlyGrid<false> G_monthlySnowFall;
        MonthlyGrid<false> G_monthlyRainFall;
        MonthlyGrid<false> G_monthlySnowMelt;
        MonthlyGrid<false> G_monthlySnowEvap;
        MonthlyGrid<false> G_monthlySoilWaterAvg;

	//monthly storages
	MonthlyGrid<> G_monthlyCanopyStorage;
	MonthlyGrid<> G_monthlySnowStorage;
	MonthlyGrid<> G_monthlySoilStorage;
	MonthlyGrid<> G_monthly_totalWaterInStorages;
	MonthlyGrid<> G_monthly_totalWaterInStorages_mm;
        MonthlyGrid<> G_monthlyArray; // temporary array used to convert double to float values

        // daily output option 31 start
        Daily31Grid<false> G_daily31LAI;
        Daily31Grid<false> G_daily31Kc;
        Daily31Grid<false> G_daily31Albedo;
        Daily31Grid<false> G_daily31ExtRad; // W/m2
        Daily31Grid<false> G_daily31ShortDownRad; // W/m2
        Daily31Grid<false> G_daily31ShortUpRad; // W/m2
        Daily31Grid<false> G_daily31LongDownRad; // W/m2
        Daily31Grid<false> G_daily31LongUpRad; // W/m2
        Daily31Grid<false> G_daily31NetLongRad; // W/m2
        Daily31Grid<false> G_daily31NetShortRad; // W/m2
        Daily31Grid<false> G_daily31NetRadiation; // W/m2
        Daily31Grid<false> G_daily31Interception;
        Daily31Grid<false> G_daily31canopyWaterContent;
        Daily31Grid<false> G_daily31maxcanopyWaterContent;
        Daily31Grid<false> G_daily31Precip;
        Daily31Grid<false> G_daily31PET;
        Daily31Grid<false> G_daily31TotalPET;
        Daily31Grid<false> G_daily31SurfaceRunoff;
        Daily31Grid<false> G_daily31GwRecharge;
        Daily31Grid<false> G_daily31SnowCoverAvg;
        Daily31Grid<false> G_daily31SWE;
        Daily31Grid<false> G_daily31SnowFall;
        Daily31Grid<false> G_daily31SnowMelt;
        Daily31Grid<false> G_daily31SoilWaterAvg;

	double canopyEvapoExp;
	double runoffFracBuiltUp;
	
	// variables which are read from 'LCT_22.DAT'
	void readLCTdata(const string inputDir);
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

	// use given albedo in the first step
	MonthlyGrid<true,float> G_AlbedoCycle;

	double a_s, b_s;
	double a_c_arid, b_c_arid, a_c_humid, b_c_humid;	// cloudiness coefficients

        Grid<int> thresh_elev;	// threshold elevation as level for minimum temperature in high elevated cells
	
	// for cell AET calculation
	Grid<> dailyCellLandAET;
	Grid<> dailyCellLandAET_uncorr;
};
#endif
