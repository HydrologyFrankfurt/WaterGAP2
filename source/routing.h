
/***********************************************************************
*
* see former changes at file routing.h.versioninfos.txt
*
***********************************************************************/

#include <vector>  
#include <fstream> 
#include "def.h"
#include <array>
using namespace std;

class routingClass {
  public:
	routingClass();
	~routingClass();
	void init(const short nBasins);
	void setParameters();
	void checkTimeStep();
    // FP20161018N002 Reservoir operation start years  // FP20161018N002: Do NOT introduce nSpecBasins, as superbasin area calculation probably wrong & obsolete
    void annualInit(const short year);
    void daily31outInit();
        void daily365outInit();   // daily output option 365
	void initLakeDepthActive();
	void initWetlDepthActive();
	void setStoragesToZero();
        void startendOutInit();

    // FP20161018N002 Reservoir operation start years
    void initFractionStatus();
    void setLakeWetlToMaximum(const short start_year);
    void routing(const short year, short day, short month, const short day_in_month); // Update to access current year, month, and day_in_month (or day) for checks, Felix Portmann 2015
        void updateLandAreaFrac();
	void annualWaterUsePostProcessing(short year);
	double getActualStorage(int cellNumber, short selection);
	double getActualStorageRatio(int cellNumber, short selection);
    // FP20161018N002 Reservoir operation start years
    double writeAnnualReservoirCapacityandLandStorageChange(char *outputDir, const short year);
	double getAnnualDischarge(const int stationNumber);

	double getAnnualUpstreamInflow(const int stationNumber);
	double getAnnualSatisfWaterUse(const int stationNumber);
	void writeLakeStorageGrids(char *outputDir, int year);
    void writeLakeStorRatioGrids(char *outputDir, int year);
	void writeMonthlyGrids(char *outputDir, int year);
    void writeTotalCellRunoffGrid(char *outputDir, int year); // HMS 2014-06-04 reintroduced for CFA calculation
    void writeDaily31Grids(char *outputDir, const int year, const int month);  // daily output option 31
    void writeDaily365Grids(char *outputDir, const int year);  // daily output option 365
    void writeStartendGrids(char *outputDir, const int year);
	//void writeAnnualCellRunoffGrid(char *outputDir, int year);  // HMS 2014-06-04 commented out
    // FP20161018N002 Reservoir operation start years
    void writeAnnualFreqGrids(char *outputDir, int year);
    void updateGloResPrevYear_pct();
    void storeSingleCellDailyValuesToFile(const short year);	  // HMS 2013-11-21
void createTestOutFile(const char *outputDir); // HMS Testout
void closeTestOutFile(); // HMS Testout
    void createDailyFiles(const char *outputDir);
	void closeDailyFiles();
	void appendDailyDischarge(const short actualYear, const short year);
	void appendDailyVelocity(const short actualYear, const short year); //------------for velocity file
	void appendDailyLakeWetl(const short actualYear, const short year);
	void createMonthlyDischargeFile(const char *outputDir);
	void closeMonthlyDischargeFile();
	void appendMonthlyDischarge(const short actualYear, const short year);
	//void appendMonthlyVelocity(const short actualYear, const short year);  //------------for velocity file
	void createAnnualDischargeFile(const char *outputDir);
	void closeAnnualDischargeFile();
	void appendAnnualDischarge(const short actualYear, const short year);
	void appendCalibInfoToAnnualDischargeFile(const double gamma);

	// routines to get information for SA
	void getAnnualDischargeVector(vector < float >&annualDischargeVector);
	void getMonthlyDischargeVector(short month, vector < float >&monthlyDischargeVector);
	void getMinMonthlyDischargeVector(vector < float >&monthlyDischargeVector);
	void getMinToMaxDischargeVector(vector < float >&monthlyDischargeVector);

        // copied from water_use.h
        void calcNextDay_M(const short day);

        // reading dailyNUs und dailyNUg
        void dailyNUInit(const char *input_dir, const short year);

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

	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ new flow velocity
	//float  getRiverVelocity(const float G_riverVeloSlope,
	//const float transportedVolumePrevStep, int n, const int TimeStepRowPerDay);
	//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ new flow velocity:end

        double G_statCorrFact[ng];	// correction factor at station locations
	
	//short G_landAreaFrac[ng];	// % of cell
        double G_landAreaFrac[ng];
        double G_landWaterExclGloLakAreaFrac[ng];
        double G_landAreaFracNextTimestep[ng]; //HMS 2013-11-22 gwr_swb
        double G_landAreaFracPrevTimestep[ng]; //HMS 2015-04-29
    // FP20161018N002 Reservoir operation start years
    short resYearCurrentToUse = 0;
    double G_landfreq[ng];
    double G_fwaterfreq[ng];
    double G_landfreq_prevyear[ng];
    double G_fwaterfreq_prevyear[ng];
    double G_glores_prevyear[ng];
    double G_glores_change[ng];
    double G_landfreq_change[ng];
    double G_fwaterfreq_change[ng];
//    double G_landAreaFrac_prevyear[ng];  // replaced by G_landAreaFracPrevTimestep[n] in annualInit
    double G_landAreaFrac_change[ng];
    double G_landStorageChange[ng];

	double (*G_monthlyResStorage)[12];
	double (*G_monthlyRiverStorage)[12];
	double (*G_monthlyLocLakeStorage)[12];
	double (*G_monthlyLocWetlStorage)[12];
	double (*G_monthlyGloLakeStorage)[12];
	double (*G_monthlyGloWetlStorage)[12];
	double (*G_monthlySwStorage)[12];

    double (*G_dailySatisfiedUse)[365];
	double (*G_daily_lateral_water_storage)[365]; //HMS 2013-11-22 gwr_swb
    // daily storage
    double (*G_daily31SurfStor)[31];
    double (*G_daily31GwStor)[31];
    double (*G_daily365SurfStor)[365];
    double (*G_daily365GwStor)[365];

    //reduction factors now as arrays
    double G_locLakeAreaReductionFactor[ng];
    double G_locWetlAreaReductionFactor[ng];
    double G_gloWetlAreaReductionFactor[ng];
    double G_riverAreaReductionFactor[ng];
    double G_gloLakeEvapoReductionFactor[ng];
    double G_gloResEvapoReductionFactor[ng];

    // storage volume of lakes and wetlands [km3]
    double G_locLakeStorage[ng];
	double G_locLakeMaxStorage[ng];	//WG22b: required for Option "use_alloc" (water abstraction from second cell)
    double G_locWetlStorage[ng];
    double G_gloLakeStorage[ng];
	double G_gloLakeMaxStorage[ng];	//WG22b: required for Option "use_alloc" (water abstraction from second cell)
    double G_gloWetlStorage[ng];
    double *G_gloResStorage; // for new reservoir algorithm
	double G_gloResMaxStorage[ng];	//WG22b: required for Option "use_alloc" (water abstraction from second cell)
    double G_riverStorage[ng];			// [km3]
	double G_riverStorageMax[ng];

	double G_actualUse[ng]; // used for new use_allocation (M.Hunger 2/2006)
        double (*G_monthlyGwStorage)[12]; // moved from private to public

// FP 2015 removed, because unused
//	float *G_landAreaOfCell;
//	float *G_contAreaOfCell;
//	float *G_test_cellArea;

// FP 2015 removed, because unused
        // copied from water_use.h
//    double G_totalDailyUse[ng];	// this grid is used by the routing algorithm
//	double G_dailyWaterUse[ng][12];	// [km3/day] based on monthly averages
	
        // matrices must be of type float
        double G_dailyNUs[ng][12];
        double G_dailyNUg[ng][12];

    double G_dailyNUsAggregated[ng][12];
//    double G_dailyNUsRedistr[ng][12]; // FP 2015 removed, because unused
        //CR: G_fracAggNUs nicht mehr benötigt
    //double G_fracAggNUs[ng][12];
	int G_glwdUnit[ng];
    // FP20161018N002 Reservoir operation start years
    // GLWDunit pattern with index to outflow grid cell Image22-Number (n = Image22nr - 1)
    int G_outflow_cell_assignment[ng];
    // FP20161018N002 Reservoir operation start years
    // Storage transfer from surface storages
    double canopyWaterContent_change_km3;
    double soilWaterContent_change_km3;
    double snowWaterContent_change_km3;
    double SnowInElevation_n_change_km3;
    double SnowInElevation_elev_change_km3;

    // Mean NUs for calculating mean demand
    double G_mean_NUs[ng];
	
        // for calcNextDay_M, NUg/NUs at the respective day
	double G_dailydailyNUs[ng];
	double G_dailydailyNUg[ng];
	double G_dailydailyNUsAggregated[ng];
//	double G_dailydailyNUsRedistr[ng]; // FP 2015 removed, because unused

    // Lakes, wetlands and reservoirs fractions or areas

    // Fractions of grid cell area (in percent of cell area)
    // HMS 2013-11-21 needs to be public for listing characteristics of single cell out list
    double G_glo_lake[ng];	// for STANDARD (G_GLOLAK.UNF0)
    double G_loc_lake[ng];	// for STANDARD (G_LOCLAK.UNF0)

    //signed char G_glo_wetland[ng];	// % of cell
    double G_glo_wetland[ng];
    //signed char G_loc_wetland[ng];	// % of cell
    double G_loc_wetland[ng];
    // FP20161018N002 Reservoir operation start years
    double G_reg_lake[ng]; // regulated lake fraction (dynamic?)
    double G_glo_res[ng]; // reservoir fraction (dynamic)
    // Absolute areas in km2
    // int G_lake_area[ng];	// km2     // obsolete
    // int G_reservoir_area[ng];	// km2    // obsolete
    double G_lake_area[ng]; // km2
    double G_reservoir_area[ng]; // actual (yearly) reservoir area km2 (only selected reservoirs in operation in the current year)
    double G_reservoir_area_full[ng]; // maximum reservoir area km2 (if all reservoirs are in operation)

    // Indicator of status at model start (first model year):
    // - initialized with value 0 and
    // - changed to value 1 after first setting
    // (1) for updateGloResPrevYear
    short statusStarted_updateGloResPrevYear;
    // (2) for G_fwaterfreq & G_landfreq
    short statusStarted_landfreq_fwaterfreq_landAreaFrac[ng];

    // Indicator of status of landAreaFracNextTimestep
    // - initialized with value 0 and
    // - changed to value 1 after first setting
    short statusStarted_landAreaFracNextTimestep[ng];

    inline double getLandAreaFrac(int n) { // HMS 2013-11-22 landAreaFrac at 1st day or not? gwr_swb
       if (0 == statusStarted_landAreaFracNextTimestep[n])  // First model day, to be used in daily.calcNewDay before changes are applied in routing.routing
           return G_landAreaFrac[n];
       else  // Following model days
           return G_landAreaFracNextTimestep[n];
    }

	// Cell-specific active lake depth (in km) as array to enable individual time series  // FP
	// Start with calibration parameter
	double G_lakeDepthActive[ng];  // km

	// Cell-specific active wetland depth (in km) as array to enable individual time-series  // FP
	// Start with calibration parameter
	double G_wetlDepthActive[ng];  // km

//	// Superbasin calculation currently obsolete, probably not working correctly in WaterGAP2.2b, do in postprocessing  // FP 2015-10
//	inline double getLandAreaOfBasin(int n) {
//		return basinLandArea[n];
//	}

    // minimum percentage of land area fraction smaller than water fractions considered as error/strange (in percent of cell area)
    double limit_errorstrange_landAreaFrac_pct = -0.001;

    // maximum difference of snow water content change (in cell and in corresponding subgrid elevation levels) when introducing new reservoir on first day of the year
    double limit_errorstrange_snowWaterContent_change_km3 = 0.001;

  private:

    double updateNetAbstractionGW(int n, int month);
  
  int i_nbc;
  int secondCell;

	//GFZ:GLOLAK begin - variable definition
	long big_lakes_cells[2][ng_glolakcells];
	float big_lakes_area_frac[ng_glolakcells];
	//float G_tempGloLakeStorage[ng_glolakcells]; // ET: als locale Variable in writeMonthlyGrids()
	//float G_tempResStorage[ng_glolakcells];     // ET: als locale Variable in writeMonthlyGrids()
	//GFZ:GLOLAK end

	short checkTimeStep(double riverVelocity, short numberOfTimeSteps);
	
	double getRiverVelocity(const double RiverSlope, const double G_riverBottomWidth, const double Roughness,
							double G_riverInflow, const int n);

	double getNewRiverVelocity (const double RiverSlope, const double G_riverBottomWidth, const double Roughness,    // DC 2018-02 storage-based FV
								double G_riverStorage, double G_riverLength, const int n);

    double G_potCellRunoff[ng];	// needed to calculate correction factors // HMS 2014-06-04 reintroduced

	double G_riverInflow[ng];			// needed for new routing algorithm in WG2.2; [km/d]
	double G_AnnualCellRunoff[ng];			// needed to calculate correction factors

	// amount of water transported through river segments
	// (= availability)
	// float G_annualRiverAvail[ng];

	// velocity
	double G_RiverSlope[ng]; 		// slope in [m/m]
	double G_RiverWidth[ng]; 		//[m] 
	double G_Roughness[ng]; 
	double G_riverLength[ng]; 		// length of river stretch; for routing in WG3.1
	double G_bankfull_flow[ng];
	double G_RiverWidth_bf[ng]; 		//[m]  river top width at bankfull stage
	double G_RiverDepth_bf[ng]; 		//[m]  river top width at bankfull stage
	double G_riverBottomWidth[ng];
	// velocity

    // river evaporation
    double G_riverAreaFrac[ng];
    double G_maxRiverAreaFrac[ng];
    double G_riverAreaFracNextTimestep_Frac[ng];
    double G_riverAreaFrac_ChangeFrac[ng];
	double riverDepth;
	double crossSectionalArea;

	// lakes, wetlands and reservoirs areas

	//float G_glo_lake[ng];	// for STANDARD (G_GLOLAK.UNF0)
	//float G_loc_lake[ng];	// for STANDARD (G_LOCLAK.UNF0)

	//signed char G_glo_wetland[ng];	// % of cell
	//float G_glo_wetland[ng];
	//signed char G_loc_wetland[ng];	// % of cell
	//float G_loc_wetland[ng];
	//int G_lake_area[ng];	// km2
	//int G_reservoir_area[ng];	// km2

	int G_neighbourCell[ng][8]; // used for new use_allocation (M.Hunger 2/2006)

	// cellnumber (1-66896)
    // double G_satisfiedUse[ng];   // CR 2015-09: A distinction between satisfied use and actual use is not possible in 22b.
    using gridarray = std::array<double, ng>;
    using monthlygridarray = std::array<std::array<double, 12>, ng>;
	double G_totalUnsatisfiedUse[ng];
	double G_fractreturngw_irrig[ng];
	gridarray G_dailyRemainingUse{}; // daily remaining use, which is used in next time step to adapt NAg
	monthlygridarray G_monthlyWUIrrig{}; // daily water use irrigation, needed for updateNAg
    monthlygridarray G_monthlyCUIrrig{}; // daily consumptive use irrigation, needed for update NAg
    monthlygridarray G_monthlySatisAllocatedUseinSecondCell{};
    monthlygridarray G_monthlyRedistributeAllocUse{};
	double G_UnsatisfiedUsePrevYear[ng];
	
    //CR: introduced in WG22b
    double G_AllocatedUse[ng];
    double G_unsatUseRiparian[ng];
    double G_remainingUse[ng];
    double G_remainingUseRes[ng];
    double G_remainingUseGloLake[ng];
    double G_remainingUseGloLakeRedistr[ng];
    double G_remainingUseResRedistr[ng];	
	double G_totalDesiredUseCell[ng];
	double G_daily_allocatedUseNextDay[ng];
	double G_AllocUseToNeigborCell[ng];
	gridarray G_UnsatAllocUse{};
	gridarray G_PrevUnsatAllocUse{};
    gridarray G_PrevTotalUnstatisfiedUse{};
	double G_dailyAllocatedUse[ng];
	gridarray G_dailySatisAllocatedUseInSecondCell{};
	double G_daily_UnsatAllocUseNextDay[ng];

	// arrays used to detect use stisfaction from neighbour cells
	// may be calculated and written to file for test purposes
	/*
	double G_extractedFromNeighb[ng];
	double G_deliveredToNeighb[ng];
	*/

	//int G_upstreamCells[ng][9];	// because of the routing algorithm not used any more
	int G_downstreamCell[ng]; 		// used for new use_allocation (M.Hunger 2/2006)
	int G_routingCell[ng]; 			// used for new routing order in WaterGAP3.1 (F.Voss 3/2008)
//	char G_LDD[ng];             	// CR 2015-05: used in WG22 to set 'transportedVolume' inland sinks to zero.
									// In WG22b, 'transportedVolume' is computed for inland sinks; it can be interpreted as water resource or additional AET.
	int G_SecondCell[ng];			// CR 2015-09: used to redistribute unsatisfied allocated water use over neighboring cells.
	int G_CellPositionRoutOrder[ng];
	
	//double G_landAreaFrac[ng];	// % of cell
	// Superbasin calculation currently obsolete, probably not working correctly in WaterGAP2.2b, do in postprocessing  // FP 2015-10
//	double *basinLandArea;

	// distance between cells for all directions
	// float G_celldistance[9][360];
	// double verticDist;	// vertical distance is independent of latitude

        // add consistent output for Precipitation (2.2) (e.g. precip only in global lake outflow cell)
        double G_GloLakePrecip[ng];
        double G_ReservoirPrecip[ng];
        double G_dailyRiverPrecip[ng];
        double (*G_monthlyConsistentPrecip)[12];
        double (*G_daily31ConsistentPrecip)[31];
        double (*G_daily365ConsistentPrecip)[365];

		// Groundwater_km3
        double G_groundwaterStorage[ng];
        // groundwater recharge below surface water bodies HMS 2013-11-22
        double gwr_loclak[ng]; // groundwater recharge below local lakes
        double gwr_glolak[ng]; // ...global lakes
        double gwr_locwet[ng]; // ...local wetlands
        double gwr_glowet[ng]; // ...global wetlands
        double gwr_res[ng]; // ...reservois
        double gwr_rivers[ng]; // ...rivers
        double G_gwr_swb[ng]; // ...below all surface waterbodies
        double G_fswbLandAreaFracNextTimestep[ng]; // fraction of surface water bodies without glo lakes and reservoirs, for adaption of landAreaFrac
        double G_fswbLandAreaFrac[ng]; // previous time step fraction of surface water bodies without glo lakes and reservoirs, for adaption of landAreaFrac
        double G_fswbLandAreaFrac_ChangePct[ng]; // change of fraction of surface water bodies without glo lakes and reservoirs, for adaption of landAreaFrac
        double G_fLocLake[ng]; //adapt for swb change and riverfractionchange
        double G_fLocWet[ng]; //adapt for swb change and riverfractionchange
        double G_fGloWet[ng]; //adapt for swb change and riverfractionchange
        double G_fGloLake[ng]; //do not adapt at all, just for output purposes
        double G_fswb[ng]; // fraction of surface waterbodies at the cell, can have values > 1 in outflow cells
        double G_fswbInit[ng]; // initial fraction of surface waterbodies at the cell including glo lakes and glo reservoirs (just for output)
        double G_fswb_catchment[ng]; // fraction of surface waterbodies at the cell including catchment area of swb (excluding glores and glolake)
        double G_gloWetlEvapoReductionFactor[ng]; // HMS test for a consistent evapo output
	// added for cell AET calculation (2.1f)
	double G_dailyLocLakeEvapo[ng]; // [mm]
	double G_dailyGloLakeEvapo[ng]; // [mm]
    double G_dailyResEvapo[ng]; // [mm]
	double G_dailyLocWetlEvapo[ng]; // [mm]
	double G_dailyGloWetlEvapo[ng]; // [mm]
        double G_locwetlextent[ng]; // [km2]
        double G_glowetlextent[ng]; // [km2]
    // added for river evaporation
    double G_dailyRiverEvapo[ng]; // mm
//    double G_locLakeEvapoRedFactor[ng]; //CR: only used for testing version 22b; remove later

    double G_daily_res_in[ng];
    double G_daily_res_out[ng];

//	double lakeDepth, wetlandDepth;  // OBSOLETE now read from calibration parameter file // FP
	double lakeOutflowExp, wetlOutflowExp;
//	double glo_storageFactor, loc_storageFactor;	//double storageFactor;  // OBSOLETE now read from calibration parameter file // FP
//	double k_g;  // OBSOLETE now read from calibration parameter file // FP

	double defaultRiverVelocity;
	//double vStep;	// velocity relative to time step [km per time step]
	short timeStepsPerDay, defaultTimeStepsPerDay;
	short nSpecBasins;	// number of specified basins
	//int TimeStepRowPerDay[ng]; // time steps for each row, dependend on latitude
	double riverVelocity;     // needs to be corrected
	double maxRiverVelocity ;   // maximum river velocity
	int maxNumberOfTimeSteps;  // maximum number of time steps

	// reduction of open water PET (2.1f)
	double evapoReductionExp; 
	double evapoReductionExpReservoir; 
	
	// arrays to store daily values
	double (*dailyRiverDischarge)[365];
	double (*dailyLocLakeStorage)[365];
	double (*dailyLocWetlStorage)[365];
	double (*dailyGloLakeStorage)[365];
	double (*dailyGloWetlStorage)[365];
	double (*dailyResStorage)[365];
    double (*G_daily365TotalWaterInStorages_km3)[365];
    double (*G_daily365TotalWaterInStorages_mm)[365];
    double (*G_daily365GwRunoff)[365];
    double (*G_daily365LandAreaFrac)[365];
    // special grid with start/end value of each year (for water balance studies)
    double (*G_startendTotalWaterInStorages_km3)[2];

	double (*dailyRiverVelocity)[365]; // +++++++++++++++++++++++++++ new flow velocity

	// grids with monthly values
	double (*G_monthlyRiverAvail)[12];
	double (*G_monthlyRiverInUpstream)[12];
	double (*G_monthlyCellRunoff)[12];
    double (*G_monthlyPotCellRunoff)[12];
    double (*G_monthlyPotCellRunoffDeficit)[12];
	double (*G_monthlyCellSurfaceRunoff)[12];
	// added for cell AET calculation (2.1f)
	double (*G_monthlyCellAET)[12]; // [mm]
        double (*G_monthlyCellAETWCa)[12];
        double (*G_monthlyOpenWaterEvap)[12]; // [mm]
	// added for river velocity
	double (*G_monthlyVelocity)[12];
	// total surface water storage (lakes/wetlands/rivers/canopies)
	double (*G_monthlySurfStor)[12];
	double (*G_monthlySurfStor_mm)[12];
        double (*G_monthlyMinRiverAvail)[12];
        double (*G_monthlyMaxRiverAvail)[12];
        double (*G_monthlyGwrSwb)[12];
        double (*G_monthlyFswb)[12];
        double (*G_monthlyLandAreaFrac)[12];
        double (*G_monthlyRiverAreaFrac)[12]; // river evaporation
        double (*G_monthlyRiverPET)[12]; // river evaporation

        double (*G_monthlyGwRunoff)[12];
        double (*G_monthlyLocWetlExtent)[12];
        double (*G_monthlyGloWetlExtent)[12];

        // grids with daily values
	double (*G_daily31RiverAvail)[31]; 
	double (*G_daily31Velocity)[31];
	double (*G_daily31CellRunoff)[31];
	double (*G_daily31CellSurfaceRunoff)[31];
    double (*G_daily31LocLakeStor)[31];
    double (*G_daily31LocWetlStor)[31];
    double (*G_daily31GloLakeStor)[31];
    double (*G_daily31GloWetlStor)[31];
    double (*G_daily31RiverStor)[31];
    double (*G_daily31ResStor)[31];
    double (*G_daily31CellAET)[31];
    double (*G_daily31TotalWaterInStorages_km3)[31];
    double (*G_daily31TotalWaterInStorages_mm)[31];
	double (*G_daily31GwrSwb)[31];
    double (*G_daily31Fswb)[31];
    double (*G_daily31GwRunoff)[31];
    double (*G_daily31LandAreaFrac)[31];
    double (*G_daily31GwrunSurfrun)[31];
    double (*G_daily31CellAETWCa)[31];

	double (*G_daily365RiverAvail)[365]; 
	double (*G_daily365Velocity)[365];
	double (*G_daily365CellRunoff)[365];
    double (*G_daily365CellRunoff_mm)[365];
	double (*G_daily365CellSurfaceRunoff)[365];
    double (*G_daily365LocLakeStor)[365];
    double (*G_daily365LocWetlStor)[365];
    double (*G_daily365GloLakeStor)[365];
    double (*G_daily365GloWetlStor)[365];
    double (*G_daily365RiverStor)[365];
    double (*G_daily365ResStor)[365];
	double (*G_daily365CellAET)[365];      
    double (*G_daily365GwrSwb)[365];
    double (*G_daily365Fswb)[365];
    double (*G_daily365GwrunSurfrun)[365];
    double (*G_daily365CellAETWCa)[365];
        //double (*G_daily365SurfStor)[365];
		
//CR 2015-09 TEST Optionen 19, 20, 28
// CR ##########################################################################		
//	double (*G_daily365ActualUse)[365];
//	double (*G_daily365AllocUse)[365];
//	double (*G_daily365AllocUseNextDay)[365];
//	double (*G_daily365TotalUnsatUse)[365];
// 	double (*G_daily365UnsatAllocUseNextDay)[365];
// CR ##########################################################################		

        // for new reservoir algorithm
	// monthly output grids
	double (*G_monthlyResStorageRatio)[12];
	double (*G_monthlyResInflow)[12];
	double (*G_monthlyResOutflow)[12];
	double (*G_monthlySatisfiedUse)[12];
        double (*G_monthlyActualUse)[12];
        double (*G_monthlyNUg)[12]; // HMS 2017-02-03 to finally allow an output of total consumptive water use
	// CR 15-08-19 implemented for test output:
    double (*G_monthlyAllocatedUse)[12]; // FP 2015, not used anymore, at the moment, implement if needed again

        // for new reservoir algorithm
	//files which do not change in time
	char (*G_res_type);    //reservoir type
	char (*G_start_month); //start of operational year
	
        double (*G_mean_outflow);   //mean annual inflow to global lakes plus P - PET of global lake (derived from G_RIVER_AVAIL) (1971-2000)
    // FP20161018N002 Reservoir operation start years
    //char (*G_reg_lake_status);    //regulated lake status (1 == regulated lake)
    char G_reg_lake_status[ng];    //regulated lake status (1 == regulated lake)
    int (*G_res_start_year);
    // always needed
    double G_stor_cap[ng]; // storage capacity (km3) (from published data) of reservoirs in operation in current year
    double G_actualStorageCapacity[ng]; // array to save storage capacity of current year
    double (*G_stor_cap_full);      // storage capacity (from published data) of reservoirs when all reservoirs are in operation

	double (*G_mean_demand);   //demand of reservoir cells, aggregated over downstream cells (max 5)
        double (*G_mean_cons_use); //mean annual consumptive use per cell (1971-2000)
	
	// ET: use of G_downstreamCell
	// int    G_dsc_res[ng][reservoir_dsc]; //downstream cells of reservoirs (max reservoir_dsc, 0 = fill value), only for irrigation reservoirs
	
	double (*G_alloc_coeff)[reservoir_dsc];//allocation coefficient (< 1 if more than one reservoir is upstream), only for irrigation reservoirs
	// ET: wird nicht gebraucht double G_temp_array[ng][reservoir_dsc];
	double (*K_release);          
	double (*annual_release);     


	ofstream dailyDischargeFile;
        ofstream testoutFile; // HMS testout
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

    // temporary variables (for checking)
    double tmp_value_double_01;
    double tmp_value_double_02;
};
