
/***********************************************************************
*
* see former changes at file routing.h.versioninfos.txt
*
***********************************************************************/
#if !defined (_routing_h_)
#define _routing_h_

#include <vector>  
#include <fstream> 
#include "def.h"
#include "grid.h"
#include "configFile.h"
#include "wghmStateFile.h"
#include "additionalOutputInputFile.h"
#include "calib_param.h"
#include <array>
#include <string>
using namespace std;

class routingClass {
  public:
    const double minStorVol = 1.e-15;   // minimal storage volume (smaller volumes set to zero) to counter numerical inaccuracies
	routingClass();
	~routingClass();
    void init(const short nBasins, ConfigFile configFile, WghmStateFile &wghmState,
              AdditionalOutputInputFile &additionalOutIn);
	void setParameters();
        void setStorages(WghmStateFile &wghmState, AdditionalOutputInputFile &additionalOutIn); //AEMS: new subroutine to set storages according to the start values read
	void checkTimeStep();
    // Do NOT introduce nSpecBasins, as superbasin area calculation probably wrong & obsolete
    void annualInit(const short year, int start_month, AdditionalOutputInputFile &additionalOutIn);
    void daily31outInit();
    void daily365outInit();   // daily output option 365
    void initLakeDepthActive(calibParamClass &calParam);//OE: calParam added
    void initWetlDepthActive(calibParamClass &calParam);//OE: calParam added
	void setStoragesToZero();
    void startendOutInit();

    void initFractionStatus();
    void initFractionStatusAdditionalOI(AdditionalOutputInputFile &additionalOutIn);
    void setLakeWetlToMaximum(const short start_year);
    void routing(const short year, short day, short month, const short day_in_month, const short last_day_in_month,
                 WghmStateFile &wghmState, AdditionalOutputInputFile &additionalOutIn, const short readinstatus,calibParamClass &calParam);// OE: added calParam // Update to access current year, month, and day_in_month (or day) for checks, Felix Portmann 2015
    void updateLandAreaFrac(AdditionalOutputInputFile &additionalOutIn);
    void update_landarea_red_fac_PDAF(calibParamClass &calParam, AdditionalOutputInputFile &additionalOutIn);
	void annualWaterUsePostProcessing(short year, AdditionalOutputInputFile &additionalOutIn);
	double getActualStorage(int cellNumber, short selection);
	double getActualStorageRatio(int cellNumber, short selection);
    void writeAnnualReservoirCapacityandLandStorageChange(const string outputDir, const short year);
	double getAnnualDischarge(const int stationNumber);

	double getAnnualUpstreamInflow(const int stationNumber);
	double getAnnualSatisfWaterUse(const int stationNumber);
	void writeLakeStorageGrids(const string outputDir, int year);
    void writeLakeStorRatioGrids(const string outputDir, int year);
	void writeMonthlyGrids(const string outputDir, int year, int month2store);
    void writeTotalCellRunoffGrid(const string outputDir, int year); // for CFA calculation
    void writeDaily31Grids(const string outputDir, const int year, const int month);  // daily output option 31
    void writeDaily365Grids(const string outputDir, const int year);  // daily output option 365
    void writeStartendGrids(const string outputDir, const int year);
    void writeAnnualFreqGrids(const string outputDir, int year);
    void updateGloResPrevYear_pct();
    void storeSingleCellDailyValuesToFile(const short year);
    void createTestOutFile(const string outputDir);
    void closeTestOutFile();
    void createDailyFiles(const string outputDir);
	void closeDailyFiles();
	void appendDailyDischarge(const short actualYear, const short year);
	void appendDailyVelocity(const short actualYear, const short year); //------------for velocity file
	void appendDailyLakeWetl(const short actualYear, const short year);
	void createMonthlyDischargeFile(const string outputDir);
	void closeMonthlyDischargeFile();
	void appendMonthlyDischarge(const short actualYear, const short year);
	void createAnnualDischargeFile(const string outputDir);
	void closeAnnualDischargeFile();
	void appendAnnualDischarge(const short actualYear, const short year);
	void appendCalibInfoToAnnualDischargeFile(const double gamma);

	// routines to get information for SA
	void getAnnualDischargeVector(vector < float >&annualDischargeVector);
	void getMonthlyDischargeVector(short month, vector < float >&monthlyDischargeVector);
	void getMinMonthlyDischargeVector(vector < float >&monthlyDischargeVector);
	void getMinToMaxDischargeVector(vector < float >&monthlyDischargeVector);

    void calcNextDay_M(const short month);

    // reading dailyNUs und dailyNUg
    void dailyNUInit(const string input_dir, const short year,calibParamClass &calParam);//OE: calParam added

#ifdef _WATERGAP_CHECKS_GLOBAL_H
    // Check Felix Portmann 2015 - for invalid data
    void check_LT0__annualInit();
    void check_NaN__annualInit();
    void check_LT0__365TWS();
    void check_NaN__365TWS();
    void check_LT0__365AET();
    void check_NaN__365AET();

    void list_LT0__annualInit();
    void list_LT0__365TWS();
    void list_LT0__365AET();
    void list_NaN__annualInit();
    void list_NaN__365TWS();
    void list_NaN__365AET();
#endif

    Grid<> G_statCorrFact;	// correction factor at station locations

    Grid<> G_landAreaFrac;              // % of cell
    Grid<> G_landWaterExclGloLakAreaFrac;
    Grid<> G_landAreaFracNextTimestep;
    Grid<> G_landAreaFracPrevTimestep;

    short resYearCurrentToUse = 0;
    Grid<> G_landfreq;
    Grid<> G_fwaterfreq;
    Grid<> G_landfreq_prevyear;
    Grid<> G_fwaterfreq_prevyear;
    Grid<> G_glores_prevyear;
    Grid<> G_glores_change;
    Grid<> G_landfreq_change;
    Grid<> G_fwaterfreq_change;
    Grid<> G_landAreaFrac_change;
    Grid<> G_landStorageChange;

	// for glacier algorithm
	Grid<> G_glacAdaptedArea;
	Grid<> G_glacAreaFrac;
	Grid<> G_glacAreaFracPrevYear;
	Grid<> G_glacAreaFrac_change;
	Grid<> G_glacPrecip;

	MonthlyGrid<false> G_monthlyResStorage;
	MonthlyGrid<false> G_monthlyRiverStorage;
	MonthlyGrid<false> G_monthlyLocLakeStorage;
	MonthlyGrid<false> G_monthlyLocWetlStorage;
	MonthlyGrid<false> G_monthlyGloLakeStorage;
	MonthlyGrid<false> G_monthlyGloWetlStorage;
	MonthlyGrid<> G_monthlySwStorage;
    MonthlyGrid<false> G_monthlyMeanWaterTempRiver;
    MonthlyGrid<false> G_monthlyMeanWaterTempLocLake;
    MonthlyGrid<false> G_monthlyMeanWaterTempLocWetl;
    MonthlyGrid<false> G_monthlyMeanWaterTempGloLake;
    MonthlyGrid<false> G_monthlyMeanWaterTempReservoir;
    MonthlyGrid<false> G_monthlyMeanWaterTempGloWetl;
	MonthlyGrid<false> G_monthlyGlacierStorage;

    // daily storage
    Daily31Grid<false> G_daily31SurfStor;
    Daily31Grid<false> G_daily31GwStor;
    DailyGrid<false> G_daily365SurfStor;
    DailyGrid<false> G_daily365GwStor;

    //reduction factors now as arrays
    Grid<> G_locLakeAreaReductionFactor;
    Grid<> G_locWetlAreaReductionFactor;
    Grid<> G_gloWetlAreaReductionFactor;
    Grid<> G_riverAreaReductionFactor;
    Grid<> G_gloLakeEvapoReductionFactor;
    Grid<> G_gloResEvapoReductionFactor;

    // storage volume of lakes and wetlands [km3]
    Grid<> G_locLakeStorage;
	Grid<> G_locLakeMaxStorage;	//WG22b: required for Option "use_alloc" (water abstraction from second cell)
    Grid<> G_locWetlStorage;
    Grid<> G_gloLakeStorage;
	Grid<> G_gloLakeMaxStorage;	//WG22b: required for Option "use_alloc" (water abstraction from second cell)
    Grid<> G_gloWetlStorage;

    Grid<> G_gloResStorage; // for new reservoir algorithm
	Grid<> G_glacierStorage; // for glacier algorithm

	Grid<> G_gloResMaxStorage;	//WG22b: required for Option "use_alloc" (water abstraction from second cell)
    Grid<> G_riverStorage;			// [km3]
	Grid<> G_riverStorageMax;

	Grid<> G_actualUse; // used for new use_allocation (M.Hunger 2/2006)
    MonthlyGrid<false> G_monthlyGwStorage; // moved from private to public

    // matrices must be of type float
    MonthlyGrid<> G_dailyNUs;
    MonthlyGrid<> G_dailyNUg;

    MonthlyGrid<> G_dailyNUsAggregated;

	Grid<int> G_glwdUnit;

    // GLWDunit pattern with index to outflow grid cell Image22-Number (n = Image22nr - 1)
    Grid<int> G_outflow_cell_assignment;

    // Storage transfer from surface storages
    double canopyWaterContent_change_km3;
    double soilWaterContent_change_km3;
    double snowWaterContent_change_km3;
    double SnowInElevation_n_change_km3;
    double SnowInElevation_elev_change_km3;

    // Mean NUs for calculating mean demand
    Grid<> G_mean_NUs;
	
    // for calcNextDay_M, NUg/NUs at the respective day
	Grid<> G_dailydailyNUs;
	Grid<> G_dailydailyNUg; // grid with NAg converted to daily values based on input, is updated if there is unsatisfied NAs
    Grid<> G_reducedReturnFlow; // Grid with reduced return flow due to unsatisfied NAs from irrig see (updateNetAbstractionGW)
    Grid<> G_unsatisfiedNAsFromIrrig; // Grid with unsatisfied NAs from irrig which is associated with G_reducedReturnFlow
	Grid<> G_unsatisfiedNAsFromOtherSectors; // Grid with unsatisfied NAs from other sectors than irrig associated with G_reducedReturnFlow
    Grid<> G_reducedReturnFlowPrevYear; // Grid with reduced return flow due to unsatisfied NAs from irrig see (updateNetAbstractionGW) from previous year
    Grid<> G_unsatisfiedNAsFromIrrigPrevYear; // Grid with unsatisfied NAs from irrig which is associated with G_reducedReturnFlow from previous year
    Grid<> G_unsatisfiedNAsFromOtherSectorsPrevYear; // Grid with unsatisfied NAs from other sectors than irrig associated with G_reducedReturnFlow from previous year
	Grid<> G_dailydailyNUsAggregated;

    // Lakes, wetlands and reservoirs fractions or areas

    // Fractions of grid cell area (in percent of cell area)
    // needs to be public for listing characteristics of single cell out list
    Grid<> G_glo_lake;	// for STANDARD (G_GLOLAK.UNF0)
    Grid<> G_loc_lake;	// for STANDARD (G_LOCLAK.UNF0)
    Grid<> G_loc_res; // reservoir fractins of grid cell area

    Grid<> G_glo_wetland;       // % of cell
    Grid<> G_loc_wetland;       // % of cell

    Grid<> G_reg_lake; // regulated lake fraction (dynamic?)
    Grid<> G_glo_res; // reservoir fraction (dynamic)

    // Absolute areas in km2
    Grid<> G_lake_area; // km2
    Grid<> G_reservoir_area; // actual (yearly) reservoir area km2 (only selected reservoirs in operation in the current year)
    Grid<> G_reservoir_area_full; // maximum reservoir area km2 (if all reservoirs are in operation)

    // Indicator of status at model start (first model year):
    // - initialized with value 0 and
    // - changed to value 1 after first setting
    // (1) for updateGloResPrevYear
    short statusStarted_updateGloResPrevYear;
    // (2) for G_fwaterfreq & G_landfreq
    Grid<short> statusStarted_landfreq_fwaterfreq_landAreaFrac;

    // Indicator of status of landAreaFracNextTimestep
    // - initialized with value 0 and
    // - changed to value 1 after first setting
    Grid<short> statusStarted_landAreaFracNextTimestep;

    inline double getLandAreaFrac(int n) {
       if (0 == statusStarted_landAreaFracNextTimestep[n])  // First model day, to be used in daily.calcNewDay before changes are applied in routing.routing
           return G_landAreaFrac[n];
       else  // Following model days
           return G_landAreaFracNextTimestep[n];
    }

	// Cell-specific active lake depth (in km) as array to enable individual time series
	// Start with calibration parameter
	Grid<> G_lakeDepthActive;  // km

	// Cell-specific active wetland depth (in km) as array to enable individual time-series
	// Start with calibration parameter
	Grid<> G_wetlDepthActive;  // km

    // minimum percentage of land area fraction smaller than water fractions considered as error/strange (in percent of cell area)
    double limit_errorstrange_landAreaFrac_pct = -0.001;

    // maximum difference of snow water content change (in cell and in corresponding subgrid elevation levels) when introducing new reservoir on first day of the year
    double limit_errorstrange_snowWaterContent_change_km3 = 0.001;

  private:

    double updateNetAbstractionGW(int n);
  
    int i_nbc;
    int secondCell;

	//GFZ:GLOLAK begin - variable definition
	VariableChannelGrid<2,long,ng_glolakcells> big_lakes_cells;
	Grid<float,ng_glolakcells> big_lakes_area_frac;

	short checkTimeStep(double riverVelocity, short numberOfTimeSteps);
	
	double getRiverVelocity(const double RiverSlope, const double G_riverBottomWidth, const double Roughness,
                            double G_riverInflow, const int n, calibParamClass &calParam); // OE: added calParam

	double getNewRiverVelocity (const double RiverSlope, const double G_riverBottomWidth, const double Roughness,    // storage-based FV
                                double G_riverStorage, double G_riverLength, const int n,calibParamClass &calParam); // OE: added calParam

    Grid<> G_potCellRunoff;	// needed to calculate correction factors

	Grid<> G_riverInflow;			// needed for new routing algorithm in WG2.2; [km/d]
	Grid<> G_AnnualCellRunoff;			// needed to calculate correction factors

	// velocity
	Grid<> G_RiverSlope; 		// slope in [m/m]
	Grid<> G_RiverWidth; 		//[m]
	Grid<> G_Roughness;
	Grid<> G_riverLength; 		// length of river stretch; for routing in WG3.1
    Grid<> G_riverLengthCheck;
	Grid<> G_bankfull_flow;
	Grid<> G_RiverWidth_bf; 		//[m]  river top width at bankfull stage
	Grid<> G_RiverDepth_bf; 		//[m]  river top width at bankfull stage
	Grid<> G_riverBottomWidth;

    // river evaporation
    Grid<> G_riverAreaFrac;
    Grid<> G_maxRiverAreaFrac;
    Grid<> G_riverAreaFracNextTimestep_Frac;
    Grid<> G_riverAreaFrac_ChangeFrac;
	double riverDepth;
	double crossSectionalArea;

	VariableChannelGrid<8,int> G_neighbourCell; // used for new use_allocation (M.Hunger 2/2006)

	Grid<> G_totalUnsatisfiedUse;
	Grid<> G_fractreturngw_irrig;
	Grid<> G_dailyRemainingUse{}; // daily remaining use change, which is used in next time step to adapt NAg
    Grid<> G_withdrawalIrrigFromSwb{}; // withdrawal from irrigation assigned with G_dailyRemainingUse to adapt NAg m3 per day
    Grid<> G_consumptiveUseIrrigFromSwb{}; // consumptive use from irrigation assigned with G_dailyRemainingUse to adapt NAg m3 per day
	MonthlyGrid<> G_monthlyWUIrrigFromSwb{}; // daily water use irrigation, needed for updateNAg
    MonthlyGrid<> G_monthlyCUIrrigFromSwb{}; // daily consumptive use irrigation, needed for update NAg
    MonthlyGrid<> G_monthlySatisAllocatedUseinSecondCell{};
    MonthlyGrid<> G_monthlyRedistributeAllocUse{};
    MonthlyGrid<> G_monthlyRedistributeGLWDUse{}; // grid to redistribute use for output G_ACTUAL_WATER_CONSUMPTION_INCELLUSED which was satisfied via glwd aggregation
	Grid<> G_UnsatisfiedUsePrevYear;
	
    Grid<> G_AllocatedUse;// this array describes the remaining Use which is allocated to cell n originating from another cell, which water use couldn't be satisfied
    Grid<> G_unsatUseRiparian; // this array describes reallocated NAs which has first was allocated to outflowcell of a reservoir or global lake and then couldn't be satisfied and thus has been reallocated to original cell
	Grid<> G_daily_allocatedUseNextDay;
	Grid<> G_AllocUseToNeigborCell;
	Grid<> G_UnsatAllocUse{}; // amount of water which was allocated from cell n to another cell and couldn't be satisfied in the other cell
	Grid<> G_PrevUnsatAllocUse{};
    Grid<> G_PrevTotalUnstatisfiedUse{};
	Grid<> G_dailyAllocatedUse;
	Grid<> G_dailySatisAllocatedUseInSecondCell{};
	Grid<> G_daily_UnsatAllocUseNextDay;

	Grid<int> G_downstreamCell; 		// used for new use_allocation (M.Hunger 2/2006)
	Grid<int> G_routingCell; 			// used for new routing order in WaterGAP3.1 (F.Voss 3/2008)
    Grid<char> G_LDD;
	Grid<int> G_SecondCell;			// CR 2015-09: used to redistribute unsatisfied allocated water use over neighboring cells.
	Grid<int> G_CellPositionRoutOrder;

    // add consistent output for Precipitation (2.2) (e.g. precip only in global lake outflow cell)
    Grid<> G_GloLakePrecip;
    Grid<> G_ReservoirPrecip;
    Grid<> G_dailyRiverPrecip;
    MonthlyGrid<false> G_monthlyConsistentPrecip;
    Daily31Grid<false> G_daily31ConsistentPrecip;
    DailyGrid<false> G_daily365ConsistentPrecip;

    // Groundwater_km3
    Grid<> G_groundwaterStorage;
    // groundwater recharge below surface water bodies HMS 2013-11-22
    Grid<> gwr_loclak; // groundwater recharge below local lakes
    Grid<> gwr_glolak; // ...global lakes
    Grid<> gwr_locwet; // ...local wetlands
    Grid<> gwr_glowet; // ...global wetlands
    Grid<> gwr_res; // ...reservois
    Grid<> gwr_rivers; // ...rivers
    Grid<> G_gwr_swb; // ...below all surface waterbodies
    Grid<> G_fswbLandAreaFracNextTimestep; // fraction of surface water bodies without glo lakes and reservoirs, for adaption of landAreaFrac
    Grid<> G_fswbLandAreaFrac; // previous time step fraction of surface water bodies without glo lakes and reservoirs, for adaption of landAreaFrac
    Grid<> G_fswbLandAreaFrac_ChangePct; // change of fraction of surface water bodies without glo lakes and reservoirs, for adaption of landAreaFrac
    Grid<> G_fLocLake; //adapt for swb change and riverfractionchange
    Grid<> G_fLocWet; //adapt for swb change and riverfractionchange
    Grid<> G_fGloWet; //adapt for swb change and riverfractionchange
    Grid<> G_fGloLake; //do not adapt at all, just for output purposes
    Grid<> G_fswb; // fraction of surface waterbodies at the cell, can have values > 1 in outflow cells
    Grid<> G_fswbInit; // initial fraction of surface waterbodies at the cell including glo lakes and glo reservoirs (just for output)
    Grid<> G_fswb_catchment; // fraction of surface waterbodies at the cell including catchment area of swb (excluding glores and glolake)
    Grid<> G_gloWetlEvapoReductionFactor; // HMS test for a consistent evapo output

	// for cell AET calculation
	Grid<> G_dailyLocLakeEvapo; // [mm]
	Grid<> G_dailyGloLakeEvapo; // [mm]
    Grid<> G_dailyResEvapo; // [mm]
	Grid<> G_dailyLocWetlEvapo; // [mm]
	Grid<> G_dailyGloWetlEvapo; // [mm]
    Grid<> G_locwetlextent; // [km2]
    Grid<> G_glowetlextent; // [km2]
    Grid<> G_transportedVolume_for_AET; // [km3]

    // cell AET calculation
	MonthlyGrid<false> G_monthlyCellAET; // total cell evap [mm]
    MonthlyGrid<false> G_monthlyCellAETWCa; // total cell evap plus consumptive water use
    MonthlyGrid<false> G_monthlyOpenWaterEvap; // evap from open water bodies [mm]

    // for river evaporation
    Grid<> G_dailyRiverEvapo; // mm

    Grid<> G_daily_res_in;
    Grid<> G_daily_res_out;

	double lakeOutflowExp, wetlOutflowExp;

	double defaultRiverVelocity;
	short timeStepsPerDay, defaultTimeStepsPerDay;
	short nSpecBasins;	// number of specified basins
	double riverVelocity;     // needs to be corrected
	double maxRiverVelocity ;   // maximum river velocity
	int maxNumberOfTimeSteps;  // maximum number of time steps

	// reduction of open water PET
	double evapoReductionExp; 
	double evapoReductionExpReservoir; 
	
	// arrays to store daily values
	DailyGrid<false,double,-1> dailyRiverDischarge;
	DailyGrid<false,double,-1> dailyLocLakeStorage;
	DailyGrid<false,double,-1> dailyLocWetlStorage;
	DailyGrid<false,double,-1> dailyGloLakeStorage;
	DailyGrid<false,double,-1> dailyGloWetlStorage;
	DailyGrid<false,double,-1> dailyResStorage;
    DailyGrid<false> G_daily365TotalWaterInStorages_km3;
    DailyGrid<false> G_daily365TotalWaterInStorages_mm;
    DailyGrid<false> G_daily365GwRunoff;
    DailyGrid<false> G_daily365LandAreaFrac;

    // special grid with start/end value of each year (for water balance studies)
    VariableChannelGrid<2,double,false> G_startendTotalWaterInStorages_km3;

	DailyGrid<false,double,-1> dailyRiverVelocity;

	// grids with monthly values
	MonthlyGrid<false> G_monthlyRiverAvail;
	MonthlyGrid<false> G_monthlyRiverInUpstream;
	MonthlyGrid<false> G_monthlyCellRunoff;
    MonthlyGrid<false> G_monthlyPotCellRunoff;
    MonthlyGrid<false> G_monthlyPotCellRunoffDeficit;
	MonthlyGrid<false> G_monthlyCellSurfaceRunoff;
	MonthlyGrid<false> G_monthlyGwrunSurfrun;

	// for river velocity
	MonthlyGrid<false> G_monthlyVelocity; // velocity grid

	// total surface water storage (lakes/wetlands/rivers/canopies)
	MonthlyGrid<false> G_monthlySurfStor; // surface water storage
	MonthlyGrid<false> G_monthlySurfStor_mm; // surface water storage [mm]
    MonthlyGrid<false> G_monthlyMinRiverAvail;
    MonthlyGrid<false> G_monthlyMaxRiverAvail;
    MonthlyGrid<false> G_monthlyGwrSwb;
    MonthlyGrid<false> G_monthlyFswb;
    MonthlyGrid<false> G_monthlyLandAreaFrac;
    MonthlyGrid<false> G_monthlyRiverAreaFrac; // river evaporation
    MonthlyGrid<false> G_monthlyRiverPET; // river evaporation

    MonthlyGrid<false> G_monthlyGwRunoff;
    MonthlyGrid<false> G_monthlyLocWetlExtent;
    MonthlyGrid<false> G_monthlyGloWetlExtent;

    // grids with daily values
	Daily31Grid<false> G_daily31RiverAvail;
	Daily31Grid<false> G_daily31Velocity;
	Daily31Grid<false> G_daily31CellRunoff;
	Daily31Grid<false> G_daily31CellSurfaceRunoff;
    Daily31Grid<false> G_daily31LocLakeStor;
    Daily31Grid<false> G_daily31LocWetlStor;
    Daily31Grid<false> G_daily31GloLakeStor;
    Daily31Grid<false> G_daily31GloWetlStor;
    Daily31Grid<false> G_daily31RiverStor;
    Daily31Grid<false> G_daily31ResStor;
    Daily31Grid<false> G_daily31CellAET;
    Daily31Grid<false> G_daily31TotalWaterInStorages_km3;
    Daily31Grid<false> G_daily31TotalWaterInStorages_mm;
	Daily31Grid<false> G_daily31GwrSwb;
    Daily31Grid<false> G_daily31Fswb;
    Daily31Grid<false> G_daily31GwRunoff;
    Daily31Grid<false> G_daily31LandAreaFrac;
    Daily31Grid<false> G_daily31GwrunSurfrun;
    Daily31Grid<false> G_daily31CellAETWCa;
    Daily31Grid<false> G_daily31WaterTemp;
    Daily31Grid<false> G_daily31locLakeTemp;
    Daily31Grid<false> G_daily31locWetlTemp;
    Daily31Grid<false> G_daily31gloLakeTemp;
    Daily31Grid<false> G_daily31reservoirTemp;
    Daily31Grid<false> G_daily31gloWetlandTemp;

	DailyGrid<false> G_daily365RiverAvail;
	DailyGrid<false> G_daily365Velocity;
	DailyGrid<false> G_daily365CellRunoff;
    DailyGrid<false> G_daily365CellRunoff_mm;
	DailyGrid<false> G_daily365CellSurfaceRunoff;
    DailyGrid<false> G_daily365LocLakeStor;
    DailyGrid<false> G_daily365LocWetlStor;
    DailyGrid<false> G_daily365GloLakeStor;
    DailyGrid<false> G_daily365GloWetlStor;
    DailyGrid<false> G_daily365RiverStor;
    DailyGrid<false> G_daily365ResStor;
	DailyGrid<false> G_daily365GlacierStorage;
    DailyGrid<false> G_daily365GlacierStorage_mm;
	DailyGrid<false> G_daily365CellAET;
    DailyGrid<false> G_daily365GwrSwb;
    DailyGrid<false> G_daily365Fswb;
    DailyGrid<false> G_daily365GwrunSurfrun;
    DailyGrid<false> G_daily365CellAETWCa;

    // for reservoir algorithm
	// monthly output grids
	MonthlyGrid<false> G_monthlyResStorageRatio;
	MonthlyGrid<false> G_monthlyResInflow;
	MonthlyGrid<false> G_monthlyResOutflow;
	MonthlyGrid<false> G_monthlySatisfiedUse;
    MonthlyGrid<false> G_monthlyActualUse;
    MonthlyGrid<false> G_monthlyNUg; //  to finally allow an output of total consumptive water use
    MonthlyGrid<false> G_monthlyAllocatedUse;

	// for glacier algorithm
	MonthlyGrid<false> G_monthlyGlacierArea;
	MonthlyGrid<false> G_monthlyGlacierAreaFrac;
	MonthlyGrid<false> G_monthlyGlacierRunoff;
	MonthlyGrid<false> G_monthlyGlacierPrecip;

    // for reservoir algorithm
	// files which do not change in time
	Grid<char> G_res_type;    //reservoir type
	Grid<char> G_start_month; //start of operational year
	
    Grid<> G_mean_outflow;   //mean annual inflow to global lakes plus P - PET of global lake (derived from G_RIVER_AVAIL) (1971-2000)

    Grid<char> G_reg_lake_status;    //regulated lake status (1 == regulated lake)
    Grid<int> G_res_start_year;

    // always needed
    Grid<> G_stor_cap; // storage capacity (km3) (from published data) of reservoirs in operation in current year
    Grid<> G_actualStorageCapacity; // array to save storage capacity of current year
    Grid<> G_stor_cap_full;      // storage capacity (from published data) of reservoirs when all reservoirs are in operation

	Grid<> G_mean_demand;   //demand of reservoir cells, aggregated over downstream cells (max 5)
    Grid<> G_mean_cons_use; //mean annual consumptive use per cell (1971-2000)

	//allocation coefficient (< 1 if more than one reservoir is upstream), only for irrigation reservoirs
    VariableChannelGrid<reservoir_dsc> G_alloc_coeff;

	Grid<> K_release;
	Grid<> annual_release;


	ofstream dailyDischargeFile;
    ofstream testoutFile;
	ofstream dailyVelocityFile; // ++++++++++++++++++++++++++++++++++++++++++++++++ new flow velocity
	ofstream monthlyDischargeFile;
	ofstream annualDischargeFile;

	ofstream monthlyReservoirStorageFile;
	ofstream monthlyReservoirStorageRatioFile;
	ofstream monthlyReservoirInflowFile;
	ofstream monthlyReservoirOutflowFile;

	ofstream locLakeFile;
	ofstream gloLakeFile;
	ofstream locWetlFile;
	ofstream gloWetlFile;
	ofstream gloResFile;

    // Temporary Storage to hold Temperatures and Inflow volumes from previous day
    Grid<double> G_RiverWTempPreStep;
    Grid<double> G_RiverAvailPreStep;
    Grid<double> G_locLakeTempPrevStep;
    Grid<double> G_locWetlTempPrevStep;
    Grid<double> G_gloLakeTempPrevStep;
    Grid<double> G_reservoirTempPrevStep;
    Grid<double> G_gloWetlandTempPrevStep;
    Grid<double> SumWaterTempRiver;
    Grid<double> NumWaterTempRiver;
    Grid<double> SumWaterTempLocLake;
    Grid<double> NumWaterTempLocLake;
    Grid<double> SumWaterTempLocWetl;
    Grid<double> NumWaterTempLocWetl;
    Grid<double> SumWaterTempGloLake;
    Grid<double> NumWaterTempGloLake;
    Grid<double> SumWaterTempReservoir;
    Grid<double> NumWaterTempReservoir;
    Grid<double> SumWaterTempGloWetl;
    Grid<double> NumWaterTempGloWetl;

    // storage for all inflow cells
    VariableChannelGrid<9, int> G_inflow_cells;

    vector<double> calcWTempReturns{0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};

    // storage for ice thickness
    Grid<double>G_iceThicknesslocLake;
    Grid<double>G_iceThicknesslocWetland;
    Grid<double>G_iceThicknessGloLake;
    Grid<double>G_iceThicknessReservoir;
    Grid<double>G_iceThicknessGloWetland;
    Grid<double>G_iceThicknessRiver;

    // power plant withdrawal and consumption
    Grid<double> G_yearlyWW;
    Grid<double> G_yearlyWC;

    // temporary variables (for checking)
    double tmp_value_double_01;
    double tmp_value_double_02;
};
#endif
