// modified watergap.cpp version
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <sys/time.h>
#include "def.h"

#include "integrateWGHM.h" // OE 2018
#include "common.h" // OE 2018
#include "globals.h" // OE 2018

#include "calib_param.h"

// Definitions for global checks (all grid cells):
// variables (chk_year, chk_month, chk_day_in_month), classes/functions, exceptions
// Checks Felix Portmann 2015
#ifdef _WATERGAP_CHECKS_GLOBAL_H
#include "watergap_checks_global.h"
#endif

//Check for pregenerated git commit hash in gitcommit.h (script get_git_commit_hash.sh to execute
// before build in the IDE of choice)
#if __has_include("gitcommit.h")
#include "gitcommit.h"
#else
static const char* GIT_COMMIT_HASH = "Git version information not generated.";
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

// Conversion factors  // FP
#include "conversion.h"


void integrate_wghm_ (const char* s, ConfigFile *& configFile, WghmStateFile *& wghmState,calibParamClass *& calParam, AdditionalOutputInputFile *& additionalOutIn, SnowInElevationFile *& snow_in_elevation,long *step, long *total_steps,long *year, long *month,const char* s2) // OE 2018
{
    bool store_states_only_at_simulation_end = true;
    bool store_states = false;
    bool first_day_after_PDAF  = true;
    try{

        // directories as strings
        std::string str_filename_global; // e.g. "/home/temp140/nobackup/portmann/WaterGAP2_2/JSON/JSON_WG_param2jsonfile"
        std::string str_filename_configFile; // new configuration file that unifies TIME.DAT and DATA.DIR
        std::string str_pathfilename_global; // e.g. "/home/temp140/nobackup/portmann/WaterGAP2_2/JSON/JSON_WG_param2jsonfile"
        // #######################################
        
        std::cout <<"CHECK: initstate "<< wghmState->cell(0).tws(0) << std::endl;
        //std::cout <<"CHECK: additionalfilestatus "<< additionalOutIn->additionalfilestatus << std::endl;
        std::string configName(s);
        std::string progName(s2);//OE
        std::cout << "progName " << progName << std::endl; //OE
        //std::cout << std::fixed;
        std::cout << std::setprecision(16); // to set the setprecision for double
        int yr;int mo;int lastStep;int stepping;//OE
        yr=*year;mo=*month;lastStep=*total_steps;stepping=*step;//OE
        configFile=new ConfigFile(configName,yr,mo,progName); //OE
        if (progName == "OL") {
            std::cout << "CHECK: for the year <" << configFile->startYear << ">" << "and month <" << configFile->startMonth << ">"<< std::endl; //OE
        }
        else {
            std::cout << "CHECK: for the year <" << yr << ">" << "and month <" << mo << ">"<< std::endl; //OE
                        std::cout << "CHECK: for the year <" << configFile->startYear << ">" << "and month <" << configFile->startMonth << ">"<< std::endl; //OE
        }
        
        // For compatibility reasons (32/64-bit processors)
        // all "long" variables were changed to "int".
        // Assure that int size is 4 byte. Exit if not.
        if (sizeof(int) != 4) {
            cerr << "Size of 'int' (" << sizeof(int)
            << " bytes) is too small on this machine! 4 bytes expected. Exiting...\n"
            << endl;
            exit(-1);
        }

#ifdef _WATERGAP_CHECKS_GLOBAL_H
        // Check Felix Portmann 2015: date (year/month/day_in_month) to check for invalid data
        // global checks

        // start in year 1901 (from 1901-01-01)
        chk_year = 1901;
        chk_month = 0; // start with 0
        chk_day_in_month = 1; // start with 1
        // NO TEST (with all functions included): specify impossible year 3000
        //chk_year = 3000;
#endif

        // endif Check Felix Portmann 2015

        char filename[350];
        char filename_no[350];
        char string[350];
        short readinstatus = 1; // read in values at first model day in case it is needed (set to 0 after first model day)
        short day, day_in_month, month, actual_year; // counter for time steps
        short start_month, end_month;
        // number_of_days_in_month[12] moved here and first_day_in_month[12] introduced
        const short number_of_days_in_month[12] = { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };
        const short first_day_in_month[12] = {0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334};
        const short last_day_in_month[12] = {30, 58, 89, 119, 150, 180, 211, 242, 272, 303, 333, 364}; // AEMS: last day needed to save snow in elevation
        short totalInitYears, remainingInitYears;
        short nSpecBasins; // number of specified basins (in STATIONS.DAT)
        short nSingleCells; // HMS 2013-11-21 number of secified single cells (in SINGLECELLS.DAT)
        const int anzThreads = 8;

#ifdef _OPENMP
        omp_set_num_threads(anzThreads);
        cout << anzThreads << " Threads" << endl;
#endif

        // print information on simulation time
        // WaterGAP version number and compilation time.
        getISOdateTime(string, sizeof(string));
        cout << string << endl;
        cout << "This is WaterGAP 2.2e from" << endl;
        cout << __DATE__ << ' ' << __TIME__ << endl;
        cout << "Commithash: " << GIT_COMMIT_HASH << endl;
        cout << "--------------------------" << endl;

        //
        // initialize classes
        // OE:auskommentiert
        //calParam->callCheckJson(); //OE: checked

        options.init(*configFile);
        options.createModelSettingsFile(options.output_dir);
        options.closeModelSettingsFile();

        if (1 == options.rout_prepare)
            prepare_routing_files(options.input_dir, options.routing_dir, 1, options.resOpt);
        geo.init(options.input_dir, options.fileEndianType);

        // init climate with local data
        climate.init();
        if ((options.time_series == 1) || (options.time_series == 3)){
            climateYear.init();
            cout << "climateYear.init()" << endl;
        }
        // init glacier algorithm with local data
        if (options.glacierOpt == 1) {
            glacierYear.init();
            cout << "glacierYear.init()" << endl;
        }
        // init clcl-class (Koeppen-map)
        if (options.clclOpt)
            clcl.init();

        switch (options.petOpt) {
            case 0:         // default climate variables
                break;
            case 1:         // default climate variables + windspeed + tmin
                climate.initPM();
                break;
            case 2:         // default climate variables + windspeed + tmin + altitude
                climate.initPM();
                land.initPM(options.input_dir);
                break;
            case 3:        // default climate variables + temperature range
                climate.initHG();
                break;
            case 4:         // default climate variables + windspeed + vapour pressure
                climate.initPM();
                break;
            case 5:         // default climate variables + windspeed + vapour pressure
                climate.initPM();
                break;
            case 6:         // default climate variables + vapour pressure
                climate.initPM();
                break;
            case 7:        // default climate variables + temperature reduction factor
                climate.initMD();
                break;
        }

        nSpecBasins = cbasin.prepare(options.input_dir, options.output_dir, options.routing_dir);
        cout << "Number of specified basins: " << nSpecBasins << endl;
        
        // Initialize cell water balance
        // Global parameters at dailyWaterBalance.setParameters, not cell-specific calibration parameters
        dailyWaterBalance.init(options.input_dir, nSpecBasins);
        G_sbasin.read(options.output_dir + "/G_CALIB_BASIN.UNF2");

        directUpstSt.findDirectStations(options.routing_dir);
        if (options.outDirectUpStations) { // new output options 2.1f
            sprintf(filename, "%s/DIRECT_UPSTREAM_STATIONS.DAT", options.output_dir.c_str());
            sprintf(filename_no, "%s/NO_DIRECT_UPSTREAM_STATIONS.DAT", options.output_dir.c_str());
            directUpstSt.writeListToFile(filename, filename_no);
        }
        allUpstSt.findAllStations(options.routing_dir);
        if (options.outAllUpStations) { // new output options 2.1f
            sprintf(filename, "%s/ALL_UPSTREAM_STATIONS.DAT", options.output_dir.c_str());
            sprintf(filename_no, "%s/NO_ALL_UPSTREAM_STATIONS.DAT",    options.output_dir.c_str());
            allUpstSt.writeListToFile(filename, filename_no);
        }

        // if this is a normal simulation run, the model is only run once.
        // in case of calibration or sensitivity analysis the mode is run in a loop
        bool lastSimLoopFlag = true;

        // if calibration if combined with the senistivity analysis
        // the sampled parameters should only be set during the first loop
        bool firstSimLoopFlag = true;

        // in case of a calibration run, this flag indicates if the model is
        // still inside the calibration loop, or if it is already the test
        // simulation
        bool testRunFlag = false;

        //
        // initialize calibration
        short calStatNum; // if zero, this is no calibration run

        // variables which are only required during calibration runs
        short upSt;
        Grid<int> G_outflow_cell;
        float gamma;

        calStatNum = calibGamma.checkStations(&gamma); // zero (0) if file 'C_STATION.DAT' does not exist
        
        // this is also used as a flag to indicate if this
        // is a calibration run or not.
        // calStatNum = 0: no calibration, otherwise calibrate the model
        if (calStatNum > 0) {
            lastSimLoopFlag = false;
            cout << "Calibration: Station " << calStatNum << " selected.\n";
            calibGamma.init();
        }
        geo.calibrun = calStatNum;
        // FP & HMS Writing of masks only in case of no calibration run, i.e. stations (geo.calibrun == 0)
        // FP20161010N002: remove internal selection "if (calibrun == 0)"; now selection before calling "if (geo.calibrun == 0)"
        if (geo.calibrun == 0) {
            geo.writeContCellMask(options.output_dir); // write mask of continental cells
            geo.writeGridAreas(options.output_dir); // write characteristic areas
            geo.writeContFreq(options.output_dir); // write continent area as percentag of grid cell area  // FP20161010N001
        }
        // end if (calibrun == 0)

        const int raindaysRandomInit = -421112;

        //
        // read grid files which are independent of time
        //

        land.init(options.input_dir); // requires GCRC_2

        //-------GROUNDWATER RECHARGE (read (semi)arid cells)-----
        G_aindex.read(options.input_dir + "/G_ARID_HUMID.UNF2");
        G_LDD.read(options.routing_dir + "/G_LDD_2.UNF1");

        //================== ELEVATION =====================================
        dailyWaterBalance.G_Elevation.read(options.input_dir + "/G_ELEV_RANGE.101.UNF2");
        //=============== END ELEVATION ====================================

        lai.init(options.input_dir, options.output_dir, land.G_landCover, *additionalOutIn,*calParam); //OE: *calParam added


        //routing.createTestOutFile(options.output_dir);
        // Initializing elements necessary for routing,
        // Global parameters in routing.setParameters, but not cell-specific calibration parameters
        routing.init(nSpecBasins, *configFile, *wghmState, *additionalOutIn);

        if (calStatNum > 0) {
            
            G_outflow_cell.read(options.routing_dir + "/G_OUTFLC.UNF4");
            
        }
        
        if (1 == options.day_store)
            routing.createDailyFiles(options.output_dir);
        
        if (2 == options.day_store) {    // HMS 2013-11-21 single cell output
            nSingleCells = options.createSingleCellFileList(options.input_dir, options.output_dir);
            
            cout << "Number of single cells where daily output is created: " << nSingleCells << endl;
        }
        
        routing.createMonthlyDischargeFile(options.output_dir);
        routing.createAnnualDischargeFile(options.output_dir);
        
        
        // /*  // START EXCLUDE & END LOOP
        // Initialize cell-specific active depth of lakes and wetlands (typically from calibration parameters)  // FP
        routing.initLakeDepthActive(*calParam);//OE: *calParam added
        routing.initWetlDepthActive(*calParam);//OE: *calParam added

        // start of the loop for calibration or sensitivity analysis
        // if a normal simulation is executed the following do-while-loop is
        // only executed once
        do {
            if (additionalOutIn->additionalfilestatus == 0) { // to avoid setting initial fractions during CDA monthly runs.

                routing.initFractionStatus(); //avoid taking last landareafrac vom previous calibration loop by setting statusStarted_updateGloResPrevYear and statusStarted_landfreq_fwaterfreq_landAreaFrac to 0.
            }
            if (additionalOutIn->additionalfilestatus == 1) { // if initial CDA monthly runs to read in fractions from additionalOutputInput.
                routing.initFractionStatusAdditionalOI(*additionalOutIn); //avoid taking last landareafrac vom previous calibration loop by setting statusStarted_updateGloResPrevYear and statusStarted_landfreq_fwaterfreq_landAreaFrac to 0.
            }
            if (configFile->additionalfile.empty()) {
                std::cout << "storages are not read in from additionalOutIn but set to 0" << std::endl;
                routing.setStoragesToZero();
            } else {
                routing.setStorages(*wghmState, *additionalOutIn);
            }
            //routing.setStoragesToZero();
            if (configFile->startvaluefile.empty()) {
                std::cout << "lakes and wetland storage not read in but set to max" << std::endl;
                routing.setLakeWetlToMaximum(options.start_year);
            }
            // CH_laf:
            if((first_day_after_PDAF) && (readinstatus == 1) && (additionalOutIn->additionalfilestatus == 1)){
                //TEST
                routing.annualInit(options.start_year, configFile->startMonth, *additionalOutIn); // TODO test purpose. Arrays are already here created (unecessary), hence better split adding arrays and annual init.
                routing.update_landarea_red_fac_PDAF(*calParam, *additionalOutIn);
                first_day_after_PDAF = false;
            }
            // mark cells where calculations should be done
            if (1 == options.basin)
                for (int n = 0; n < ng; ++n)
                    G_toBeCalculated[n] = 1;
            else
                for (int n = 0; n < ng; ++n) {
                    if (G_sbasin[n] > 0)
                        G_toBeCalculated[n] = 1;
                    else if (G_sbasin[n] < 0)
                        G_toBeCalculated[n] = -1;
                    else
                        G_toBeCalculated[n] = 0;
                }

            // Filling of dailyWaterBalance.G_gammaHBV (P_GAMRUN_C)
            // Read from / write to UNF0-file G_GAMMA_HBV.UNF0 only during calibration runs  // FP
            if (calStatNum > 0) {
                dailyWaterBalance.G_gammaHBV.read(options.input_dir + "/G_GAMMA_HBV.UNF0");
            }
            // end if (calStatNum > 0)
            else {
                // Normal runs: transfer values from JSON-object (from calibration parameter file) // FP
                for (int n = 0; n < ng; n++) {
                    dailyWaterBalance.G_gammaHBV[n] = calParam->getValue(P_GAMRUN_C, n);
                }
            }
            // end else filling G_gammaHBV (P_GAMRUN_C)

            // Filling of dailyWaterBalance.G_cellCorrFact (P_CFA)
            // Filling of routing.G_statCorrFact (P_CFS)
            // Read from G_CORR_FACTOR.UNF0 & G_STAT_CORR.UNF0 only during calibration runs
            if (calStatNum > 0) {
                if (!testRunFlag) {
                    dailyWaterBalance.G_cellCorrFact.read("G_CORR_FACTOR.UNF0");
                    // make sure, correction factors in the calibration basin are set to
                    // 1, during calibration (not for the test run)
                    for (int n = 0; n < ng; ++n)
                        if (G_sbasin[n] == calStatNum) {
                            dailyWaterBalance.G_cellCorrFact[n] = 1.0;
                            routing.G_statCorrFact[n] = 1.0;
                        }
                }
            }
            // end if (calStatNum > 0)
            else {
                for (int n = 0; n < ng; n++) {
                    dailyWaterBalance.G_cellCorrFact[n] = calParam->getValue(P_CFA, n);
                    routing.G_statCorrFact[n] = calParam->getValue(P_CFS, n);
                }
            }
            // end else filling G_cellCorrFact (P_CFA) and G_statCorrFact (P_CFS)


            //replace gammas (if some are given in the STATIONS-file)
            bool gammaInStationList = false;

            for (short n = 1; n <= cbasin.numberOfCalibBasins; n++) {
                if (cbasin.gamma[n - 1] > 0) {
                    gammaInStationList = true;
                    allUpstSt.replaceGridValues(n, cbasin.gamma[n - 1], dailyWaterBalance.G_gammaHBV, cbasin.gamma);
                    if (cbasin.statCorrFactor[n - 1] > 0) {
                        routing.G_statCorrFact[cbasin.cellNum[n - 1] - 1] = cbasin.statCorrFactor[n - 1];

                        // do not set the value for the basin which should
                        // be calibrated
                        if (n == calStatNum)
                            routing.G_statCorrFact[cbasin.cellNum[n - 1] - 1] = 1.0;
                    }
                }
            }

            // Write potentially changed values to G_GAMMA_HBV_NEW.UNF0 & G_STAT_CORR.UNF0
            // at the end of calibration runs
            if (gammaInStationList && (0 == calStatNum)) {
                // if values of gamma and correction factors are given in
                // STATIONS.DAT, a grid which includes the new gammas is
                // stored to the working directory.
                // This grid might be needed for further runs.
                // Important: a consistent G_CORR_FACTOR.UNF0 has to be generated
                // separately (should be result of the calibration run)!
                cout << "Writing update of gamma and cell correction factors from file STATIONS.DAT to files G_GAMMA_HBV_NEW.UNF0 and G_STAT_CORR.UNF0" << endl;
                
                dailyWaterBalance.G_gammaHBV.write("G_GAMMA_HBV_NEW.UNF0");
                routing.G_statCorrFact.write("G_STAT_CORR.UNF0");
            }

            if (configFile->additionalfile.empty()) { //TODO: subject for improved if clause. Might be better to write similar to the if clause above but also dont hurt if it is written like this
                std::cout << "storages are not read in from additionalOutIn but set to 0" << std::endl;
                dailyWaterBalance.setStoragesToZero();
            } else {
                dailyWaterBalance.setStorages(*wghmState, *snow_in_elevation, *additionalOutIn);
            }


            if (calStatNum > 0) {
                cout << "Calibration: New value of gamma: " << gamma << endl;

                for (int i = 0; i < ng; ++i) {
                    dailyWaterBalance.G_gammaHBV[i] = gamma;
                    //dailyWaterBalance.G_cellCorrFact[i] = cellCorrFactor;
                }

                // replace gammas of upstream stations (which are already available)
                // and mark cells where calculations are necessary
                for (int n = 0; n < ng; ++n) {
                    if (G_sbasin[n] == calStatNum)
                        G_toBeCalculated[n] = 1;
                    else
                        G_toBeCalculated[n] = 0;
                }
                for (int n = 1; n <= allUpstSt.getNumberOfUpstreamStations(
                                                                           calStatNum); ++n) {
                    upSt = allUpstSt.getUpstreamStation(calStatNum, n);
                    if (cbasin.gamma[upSt - 1] > 0) {
                        cout << "Replacing gamma: " << cbasin.gamma[upSt - 1] << " for basin " << upSt << endl;
                        for (int m = 0; m < ng; ++m)
                            if (G_sbasin[m] == upSt) {
                                dailyWaterBalance.G_gammaHBV[m] = cbasin.gamma[upSt    - 1];
                                G_toBeCalculated[m] = 1;
                            }
                    } else {
                        cerr << "Gamma not available\n";
                        exit(-1);
                    }
                }
                G_toBeCalculated[G_outflow_cell[cbasin.cellNum[calStatNum - 1] - 1]    - 1] = -1;
            } //End: if (calStatNum > 0)
            if (0 != options.timeStepCheckFlag)
                routing.checkTimeStep(); // depends on G_toBeCalculated[]
            
            if (calStatNum > 0) {
                // all these files have already been created above,
                // but should be overwritten, if this is a calibration run
                if (1 == options.day_store) {
                    routing.closeDailyFiles();
                    routing.createDailyFiles(options.output_dir);
                }
                //routing.closeTestOutFile(); // HMS Testout
                routing.closeMonthlyDischargeFile();
                routing.createMonthlyDischargeFile(options.output_dir);
                
                // file with annual information should not be overwritten
                // but gamma values should be stored as additional information
                routing.appendCalibInfoToAnnualDischargeFile(gamma);
            }


            // grids with max. soil water capacity and groundwater factors
            maxSoilWaterCap.createMaxSoilWaterCapacityGrid(options.input_dir, options.output_dir, land.G_landCover, *calParam);//OE: *calParam added
            GW.createGrids(options.input_dir, options.output_dir,*calParam); //OE: added *calParam

            //
            // FP20161018N003 Enable reading input data from climate land mask
            // read climate grid files which do not change in time
            // this is for calculations of climate normal period only
            // HERE CHANGE ALSO FOR OTHER LOCATIONS
            // Read data from normal land mask
            if (0 == options.climate_spatial_resolution) {
                climate.read_climate_longtermAvg();
            }
            // Read data from climate land mask
            else if (1 == options.climate_spatial_resolution) {
                climate.read_climate_longtermAvg_withConvert();
            }

            //===============KOEPPEN CLIMATE CLASSIFICATION===========================
            if (options.clclOpt) {
                
                ifstream station_file(options.input_dir + "/G_KOEPPEN.UNF2");
                
                if (!station_file) {
                    
                    cout << "Calculate KOEPPEN map... " << endl;
                    
                    // new time series input option
                    // read long term mean data
                    // monthly values [0.01 oC]
                    climate.G_temperature.read(options.climate_dir + "/GTEMP_1971_2000.12.UNF2");
                    
                    // read monthly precipitation [mm]
                    climate.G_precipitation.read(options.climate_dir + "/GPREC_1971_2000.12.UNF2");
                    
                    // calculate climate zones
                    clcl.calcKoep();
                    // assign alpha values
                    clcl.alphaKoep();
                    
                    // write KOEPPEN map (grid)
                    if (options.outClcl) {
                        clcl.cls.write(options.output_dir + "/G_KOEPPEN.UNF2");
                    }
                } else {
                    station_file.close();
                    cout << "Read KOEPPEN map from UNF-file... " << endl;
                    // read climate zones from file
                    clcl.cls.read(options.input_dir + "/G_KOEPPEN.UNF2");
                    
                    // assign alpha values
                    clcl.alphaKoep();
                }
                
            }
            //===============KOEPPEN CLIMATE CLASSIFICATION===========================

            // start of loop for simulation years
            // calculate from start to end year
            short end_year = -1, start_year = -1;
            
            start_year = options.start_year;
            end_year = options.end_year;
            
            // for variable land cover map
            int landCoverYear = 0;
            if (options.landCoverOpt) {
                land.readLandCoverMap(options.landCoverYears[landCoverYear], options.land_cover_dir);
                landCoverYear++;
            }

            // for glacier algorithm
            if (options.glacierOpt == 1)
                glacierYear.read_initial_glacier_area(options.start_year);

            totalInitYears = options.init_years;

#ifdef _WATERGAP_CHECKS_GLOBAL_H
            cout << " CHECK totalInitYears (from options): " << totalInitYears << endl; // Check Felix Portmann 2015
#endif

            // Loop over years
            for (actual_year = start_year; actual_year <= end_year; actual_year += options.time_step) {
                cout << "Year of simulation: " << actual_year << endl;
                // for variable land cover map
                if (options.landCoverOpt && // use variable land cover map
                    landCoverYear < options.landCoverYearsInList && // last map is not read yet
                    actual_year >= options.landCoverYearsStart[landCoverYear]) { // actual_year use next map

                    land.readLandCoverMap(options.landCoverYears[landCoverYear], options.land_cover_dir);
                    landCoverYear++;
                }


                if (options.glacierOpt == 1) {
                    glacierYear.read_glacier_data_daily_per_year(actual_year);
                    glacierYear.read_glacier_data_yearly_per_year(actual_year);
                }

                // and daily climate data from yearly files
                // new time series input option start - add .365
                if (options.time_series == 1) {
                    if (0 == options.climate_spatial_resolution) {
                        cout << "Begin read_climate_data_daily_per_year(" << actual_year << ")" << endl;
                        climateYear.read_climate_data_daily_per_year(actual_year);
                        cout << "End read_climate_data_daily_per_year(" << actual_year << ")" << endl;
                    } else {
                        cerr << "reading in .365 UNF at climate landmask is not yet implemented." << endl;
                    }
                }

                //===============Permafrost Calculations===========================

                // with daily climate data for permafrost calculating we need monthly values for temperature and precipitation
                if (options.permaOpt) {
                    for (int month = 0; month < 12; ++month) {
                        // read daily temperature and precipitation from .31 files
                        if (options.time_series == 0)
                            if (0==options.climate_spatial_resolution) {
                                climate.read_climate_data_daily(month+1, actual_year);
                            } else {
                                climate.read_climate_data_daily_withConvert(month+1, actual_year);
                            }
                        //or split daily data from .365 files
                            else if (options.time_series == 1)
                                climateYear.split_climate_data_from_year_to_daily(first_day_in_month[month], number_of_days_in_month[month]);
                        
                        double dummy_temp, dummy_prec;
                        for (int n = 0; n < ng; ++n) {
                            dummy_temp = 0.; dummy_prec = 0.;
                            for (int day_in_month = 1;     day_in_month <= number_of_days_in_month[month]; day_in_month++) {
                                dummy_temp   += climate.G_temperature_d(n,day_in_month);
                                dummy_prec += climate.G_precipitation_d(n,day_in_month);
                            }
                            climate.G_temperature(n,month)  = (short)(dummy_temp/number_of_days_in_month[month] *100.);
                            climate.G_precipitation(n,month)= (short)(dummy_prec);
                        }
                    }
                }

                if (options.permaOpt) {

                    permaClass permafrost;

#pragma omp parallel for  \
shared(actual_year, cout, permafrost)

                    for (int n = 0; n < ng; ++n) {
#ifdef _OPENMP
                        if (omp_get_num_threads() != anzThreads) {
                            cout << "Falsche Threadzahl" << endl;
                            exit(1);
                        }
#endif
                        permafrost.calcFrostNumber(actual_year, n);
                    }

                    // write yearly grids
                    permafrost.writeYearlyGrids(options.output_dir, actual_year);

                }
                //===============Permafrost Calculations===========================

                // Loop over initialization years (from totalInitYears = options.init_years)
                for (remainingInitYears = totalInitYears; remainingInitYears >= 0; remainingInitYears--) {
                    // 'remainingInitYears' is different from 0
                    // only for the first simulation year

#ifdef _WATERGAP_CHECKS_GLOBAL_H
                    // Check Felix Portmann 2015 - for invalid data
                    cout << " Loop initial years: CHECK totalInitYears: " << totalInitYears << " remainingInitYears: " << remainingInitYears << endl;
#endif

                    if (remainingInitYears > 0) {
                        cout << "  Year of initialization phase: " << remainingInitYears << endl;
                    }

                    dailyWaterBalance.annualInit();
                    routing.annualInit(actual_year, configFile->startMonth, *additionalOutIn); // Do NOT introduce nSpecBasins, as superbasin area calculation probably wrong & obsolete
                    
                    // read dailyNUs und dailyNUg or reallocated values
                    if (((2 == options.subtract_use)  || (3 == options.subtract_use)) && (0 == options.time_series)) {
                        routing.dailyNUInit(options.water_use_dir, actual_year, *calParam);//OE: *calParam added
                        //cout << "Initialization of dailyNUInit for year " << actual_year << " successfully done." << endl;
                    }

#ifdef _WATERGAP_CHECKS_GLOBAL_H
                    // Check Felix Portmann 2015 - for invalid data
                    // (1) At the start of the current year(s)
                    if (actual_year >= chk_year) {
                        // (1.1) Every calculation year (at start)
                        cout << "Check initialized monthly data arrays at start of the current calculation year " << actual_year << endl;
                        routing.check_LT0__annualInit();
                        routing.check_NaN__annualInit();
                    }
#endif

                    //
                    // loop for daily calculations
                    //
                    day = 0; // initialisation for parameter day

                    // AEMS: count number of days since begin of the year
                    if (actual_year == configFile->startYear) {
                        for(int m = 1; m<configFile->startMonth;m++) {
                            day += number_of_days_in_month[m-1];
                        }
                    }

                    //write startendvalue of TWS
                    if (options.outTotalWaterInStoragesStartEndDaily_km3) {
                        //   daily.startendOutInit();
                        routing.startendOutInit();
                    }

                    //reformulated for daily output options (.31, .365)
                    if ((5 == options.grid_store) || (6 == options.grid_store)) {
                        dailyWaterBalance.daily365outInit();
                        routing.daily365outInit();

#ifdef _WATERGAP_CHECKS_GLOBAL_H
                        // Check Felix Portmann 2015 - for invalid data
                        if (actual_year >= chk_year) {
                            // (1.2) Every calculation year (at start)
                            cout << "Check initialized daily (365) data arrays at start of the current calculation year " << actual_year << endl;
                            routing.check_LT0__365TWS();
                            routing.check_NaN__365TWS();
                            routing.check_LT0__365AET();
                            routing.check_NaN__365AET();
                        }
#endif

                    }
                    //for (month = 0; month < 12; ++month) { // HMS old code (run for a whole year)
                    end_month = configFile->endMonth;
                    start_month = configFile->startMonth;
                    const int normal_year{11};
                    month = 0;
                    int a = 0;
                    if (start_year == end_year){
                        start_month = configFile->startMonth ;
                        end_month = configFile->endMonth;

                    }

                    else if(actual_year == start_year){
                        //First iteration of year loop
                        end_month = 12;
                        start_month = configFile->startMonth ;
                    }
                    else if(actual_year == end_year){
                        end_month = configFile->endMonth;
                        start_month = 01;
                    }
                    else if(actual_year  < end_year && actual_year > start_year){
                        end_month = 12;
                        start_month = 01;
                    }

                    // while(month <= normal_year) {
                    for (month = start_month-1; month < end_month; month++) {
                        //reformulated for daily output options (.31, .365)
                        if ((3 == options.grid_store) || (4 == options.grid_store)) {
                            dailyWaterBalance.daily31outInit();
                            routing.daily31outInit();
                        }
                        //
                        // read grid files for daily climate data which change in time
                        //

                        // new time series input option start
                        if (options.time_series == 0) {
                            if (0==options.climate_spatial_resolution)
                                climate.read_climate_data_daily(month + 1, actual_year);
                            else
                                climate.read_climate_data_daily_withConvert(month + 1, actual_year);
                        }
                        if (options.time_series == 1) {
                            climateYear.split_climate_data_from_year_to_daily(day, number_of_days_in_month[month]);
                        }
                        // new time series input option end

                        // AEMS: Reset monthly storages:
                        wghmState->resetCells(number_of_days_in_month[month]);

                        //// CH_laf:
                        //if((first_day_after_PDAF) && (readinstatus == 1) && (additionalOutIn->additionalfilestatus == 1)){
                        //    routing.update_landarea_red_fac_PDAF(*calParam, *additionalOutIn);
                        //    first_day_after_PDAF = false;
                        //}

                        for (day_in_month = 1; day_in_month <= number_of_days_in_month[month]; day_in_month++) {

                            if (store_states_only_at_simulation_end) {

                                if ((actual_year == end_year) && (month+1 == end_month) && (day_in_month == number_of_days_in_month[month])) {
                                    store_states = true;
                                }
                            }
                            day++; // [1..365]

                            //cout <<"Year, Month, Day : "<< actual_year<<", "<<month << ", " << day << endl;

#pragma omp parallel for  \
shared(day, month, day_in_month, actual_year, cout, dailyWaterBalance)

                            for (int n = 0; n < ng; ++n) {
                                // Only execute for continental grid cells
                                if (geo.G_contcell[n]) {
#ifdef _OPENMP
                                    if (omp_get_num_threads() != anzThreads) {
                                        cout << "Falsche Threadzahl" << endl;
                                        exit(1);
                                    }
#endif
                                    dailyWaterBalance.calcNewDay(day, month, day_in_month, last_day_in_month[month], actual_year, n, *wghmState, *additionalOutIn, *snow_in_elevation,readinstatus,*calParam); //OE: added *calParam

                                }
                                // endif geo.G_contcell[n]
                            }
                            //endfor loop over grid cells

#ifdef _WATERGAP_CHECKS_GLOBAL_H
                            // Check Felix Portmann 2015 - for invalid data
                            // Messaging subsequent checking
                            if ( (actual_year >= chk_year) && (month >= chk_month) && (day_in_month >= chk_day_in_month) ) {
                                cout << "Checking: routing of year: " << actual_year << " day: " << day << " month: " << month << " day_in_month: " << day_in_month << endl;
                            }
#endif

                            if ((2 == options.subtract_use) || (3 == options.subtract_use)) { // enhanced to reallocated water use
                                routing.calcNextDay_M(month);
                            }
                            routing.routing(actual_year, day, month, day_in_month, last_day_in_month[month], *wghmState, *additionalOutIn,readinstatus, *calParam); //OE: *calParam added   // routing through river network
                            routing.updateLandAreaFrac(*additionalOutIn);
                            

                            //AEMS: Save daily to wghmState
                            if (progName == "OL" && store_states) {
                                if(!configFile->outputdailyfile.empty()) {
                                    std::stringstream ss;
                                    ss << "_" << actual_year << "-" << std::setw(2)<<std::setfill('0')<<month+1<<"-"<<std::setw(2)<<std::setfill('0')<<day_in_month;
                                    std::string date = ss.str();

                                    std::string fn = configFile->outputdailyfile;
                                    std::string::size_type pos = fn.rfind('.');
                                    if(pos == std::string::npos)
                                        fn = fn + date;
                                    else
                                        fn.insert(pos,date);

                                    std::cout << "Save day to <"<<fn<<"> ... "<< std::flush;
                                    //wghmState.saveDay(fn, day_in_month-1); // HMS: not sure why day_in_month-1 is used here.
                                    wghmState->saveDay(fn, day_in_month-1);

                                    std::cout <<"done."<< std::endl << std::flush;

                                }
                            }
                            if (readinstatus == 1) // make sure that only once per model run states/additionalinout etc are read in
                                readinstatus = 0;
                        } // end of loop for day_in_month
                        
                        
                        
                       // HG (2020/07) added this part to apply landAreaFrac after updating it
                        for (int n = 0; n < ng; ++n){

                            double landAreaFrac = routing.getLandAreaFrac(n);

                            for (short elev = 0; elev <= 100; elev++){
                                if (landAreaFrac ==0.){
                                    snow_in_elevation->snowInElevation(n, elev)=0.;
                                }
                                else{
                                    snow_in_elevation->snowInElevation(n, elev) = dailyWaterBalance.G_SnowInElevation(n,elev) * landAreaFrac / geo.G_contfreq[n];
                                }
                            }
                            for (day_in_month = 1; day_in_month <= number_of_days_in_month[month]; day_in_month++) {
                                wghmState->cell(n).canopy(day_in_month-1) = dailyWaterBalance.G_canopyWaterContent[n] * landAreaFrac / geo.G_contfreq[n];
                                wghmState->cell(n).snow(day_in_month-1) = dailyWaterBalance.G_snow[n] * landAreaFrac / geo.G_contfreq[n];
                                wghmState->cell(n).soil(day_in_month-1) = dailyWaterBalance.G_soilWaterContent[n] * landAreaFrac / geo.G_contfreq[n];
                            }
                        }
                        // HG (2020/07) end

                        

                        // OE 2018
                        if (progName == "OL" && store_states) {
                            // AEMS: Save monthly mean and last day to wghmState
                            if (!configFile->outputmeanfile.empty())
                            {
                                std::stringstream ss;
                                ss << "_" << actual_year <<"-"<<std::setw(2)<<std::setfill('0')<<month+1;//<<"-"<<std::setw(1)<<std::setfill('0')<<day_in_month;
                                std::string date = ss.str();
                                std::string fn = configFile->outputmeanfile;
                                std::string::size_type pos = fn.rfind('.');
                                if(pos == std::string::npos)
                                    fn = fn + date;
                                else
                                    fn.insert(pos,date);

                                std::cout << "Save monthly mean to <"<<fn<<"> ... "<< std::flush;
                                wghmState->saveMean(fn);
                                std::cout <<"done."<< std::endl << std::flush;

                            }


                            if (!configFile->outputlastdayfile.empty())
                            {
                                std::string fn = configFile->outputlastdayfile;
                                std::cout << "Save last day to <"<<fn<<"> ... "<< std::flush;
                                //wghmState.saveDay(fn, number_of_days_in_month[month]-1); // HMS I dont understand why number_of_days.. is used here, should be 0, or?
                                wghmState->saveDay(fn, number_of_days_in_month[month] - 1);
                                std::cout << "number_of_days_in_month[month] - 1 <"<<number_of_days_in_month[month] - 1<<"> ... "<< std::flush;
                                std::cout <<"done."<< std::endl << std::flush;
                            }


                            if (!configFile->outputadditionalfile.empty()) // AEMS
                            {

                                std::string fn = configFile->outputadditionalfile;
                                std::cout << "Save G_days_since_start, G_growingStatus and G_PrecSum of last day to <"<<fn<<"> ... "<< std::flush;
                                additionalOutIn->save(fn);
                                std::cout <<"done."<< std::endl << std::flush;
                            }

                            if (!configFile->outputsnowlastdayfile.empty()) // AEMS
                            {
                                std::string fn = configFile->outputsnowlastdayfile;
                                std::cout << "Save snow in elevation of last day to <"<<fn<<"> ... "<< std::flush;
                                snow_in_elevation->save(fn);
                                std::cout <<"done."<< std::endl << std::flush;
                            }
                        }

                        // In the following methods, it has to be guaranteed that
                        //  - only continental grid cells are processed
                        //  - or zero values are still maintained from initialization
                        //      or written, if necessary
                        if ((actual_year >= options.evalStartYear) && (remainingInitYears == 0)) {

                            // daily output option 31 start
                            if ((3 == options.grid_store) || (4 == options.grid_store)) {
                                routing.writeDaily31Grids(options.output_dir, actual_year, month);
                                dailyWaterBalance.writeDaily31Grids(options.output_dir, actual_year, month);
                            }
                            // daily output option 31 end
                        }
                        //month++;
                    } // end of loop for month


                    // Transfer land area fraction of reservoirs (in percent) of current year in array of previous year for use in the next year
                    if (month==12)  // should only be done at end of year
                        routing.updateGloResPrevYear_pct();


                    // Write daily output only with start of evaluation
                    //  and outside initialization years
                    if ((actual_year >= options.evalStartYear) && (remainingInitYears == 0)) {
                        // daily output option 365 start
                        if ((5 == options.grid_store) || (6 == options.grid_store) || (2 == options.day_store)) {
                            routing.writeDaily365Grids(options.output_dir, actual_year);
                            dailyWaterBalance.writeDaily365Grids(options.output_dir, actual_year);
                        }
                        // daily output option 365 end
                        if (options.outTotalWaterInStoragesStartEndDaily_km3)
                            routing.writeStartendGrids(options.output_dir, actual_year);
                    }
                    // end if after start of evalaution and after initialization


                    if (0 != options.subtract_use) {
                        //routing.annualWaterUsePostProcessing(actual_year); //HMS not needed (see below)
                        if(month==12) // only at the end of a year
                        {
                            routing.annualWaterUsePostProcessing(actual_year, *additionalOutIn); // AEMS: additionalOutIn added
                            //AEMS: the outputadditionalfile has to be updated
                            if (progName == "OL" && store_states) {
                                if (!configFile->outputadditionalfile.empty()) // AEMS
                                {
                                    std::string fn = configFile->outputadditionalfile;
                                    std::cout << "Overwrite <"<<fn<<"> ... "<< std::flush;
                                    additionalOutIn->save(fn);
                                    std::cout <<"done."<< std::endl << std::flush;
                                }
                            }
                        }
                    }

                    // write daily river discharge information into file
                    if (options.day_store == 1) {
                        routing.appendDailyDischarge(actual_year, remainingInitYears);
                        routing.appendDailyLakeWetl(actual_year, remainingInitYears);
                        routing.appendDailyVelocity(actual_year, remainingInitYears); //++++++++ new flow velocity
                    }
                    routing.appendMonthlyDischarge(actual_year, remainingInitYears);
                    routing.appendAnnualDischarge(actual_year, remainingInitYears);

                    if (calStatNum > 0) {
                        if (actual_year >= options.evalStartYear) {
                            // pass simulation results to the calibration object
                            calibGamma.setRunoff(actual_year, routing.getAnnualDischarge(calStatNum));
                            calibGamma.setWaterUse(actual_year, routing.getAnnualSatisfWaterUse(calStatNum));
                            calibGamma.setUpstInflow(actual_year, routing.getAnnualUpstreamInflow(calStatNum));
                        }
                    }


#ifdef _WATERGAP_CHECKS_GLOBAL_H
                    // START Check Felix Portmann 2015 - for invalid data
                    // (2) At the end of the current year(s)
                    if (actual_year >= chk_year) {

                        // START (2.1) Every calculation year (before writing)
                        cout << "Check data at the end of the current calculation year " << actual_year << " (before writing)" << endl;

                        // Check monthly arrays used/initialized in "annualInit", e.g. AET, cell runoff
                        // monthly data - check & list
                        routing.list_LT0__annualInit(); // (selected variables)
                        routing.list_NaN__annualInit();
                        // monthly data - check & throw an exception
                        routing.check_LT0__annualInit();
                        routing.check_NaN__annualInit();

                        // Check daily output arrays used/initialized in "annualInit", e.g. total storage, cell AET
                        // daily data - check & list OR throw an exception
                        if ((5 == options.grid_store) || (6 == options.grid_store)) {
                            // can be negative = exclude
                            // only include listing if necessary
                            // routing.list_LT0__365TWS();
                            routing.list_LT0__365AET(); // only a small number of cells are negative (because discharge is larger than precipitation)

                            routing.list_NaN__365TWS();
                            routing.list_NaN__365AET();

                            // can be negative = exclude
                            // routing.check_LT0__365TWS();
                            // routing.check_LT0__365AET();
                            routing.check_NaN__365TWS();
                            routing.check_NaN__365AET();
                        }
                    }
                    // endif if (year >= chk_year)
                    // END (2.1) Every calculation year (before writing)
#endif

                }    // end of loop for remainingInitYear


                totalInitYears = 0;

                if (2 == options.day_store)
                    routing.storeSingleCellDailyValuesToFile(actual_year);
                if (actual_year >= options.evalStartYear) {

                    // write output grids
                    if (options.grid_store != 0) {
                        if ((1 == options.grid_store) && (options.outSingleStorages)) {
                            routing.writeLakeStorRatioGrids(options.output_dir, actual_year);
                            //dailyWaterBalance.writeAnnualGrids(options.output_dir, actual_year);
                            routing.writeLakeStorageGrids(options.output_dir, actual_year);
                        }
                        dailyWaterBalance.writeMonthlyGrids(options.output_dir, actual_year, configFile->startMonth); // OE:2add: options.snowCoverFrac_dir, configFile->startMonth);
                        routing.writeMonthlyGrids(options.output_dir, actual_year, configFile->startMonth);//OE:2add:, options.discharge_dir, options.reanPrecip_dir, configFile->startMonth
                    }

                    if ((options.grid_store != 0) && (calStatNum > 0)) {
                        // these grids are required to calculate
                        // the cell related correction factor
                        routing.writeTotalCellRunoffGrid(options.output_dir, actual_year);
                    }

                    // FP20161018N002 Reservoir operation start years
                    if ((0 != options.grid_store) && (options.outLandWaterFractions)) {
                        //write out GFREQ, GFREQW etc., fractions as defined in routing.annualInit at start of every year
                        routing.writeAnnualFreqGrids(options.output_dir, actual_year);
                    }

#ifdef _WATERGAP_CHECKS_GLOBAL_H
                    // Check Felix Portmann 2015 - for invalid data
                    // START (2.2) Check evaluation year (after writing)
                    // Checks, Felix Portmann 2015
                    if (actual_year >= chk_year) {

                        cout << "Check data at the end of the current evaluation year " << actual_year << " (after writing)" << endl;

                        // Check monthly arrays used/initialized in "annualInit", e.g. AET, cell runoff
                        // monthly data - check & list
                        routing.list_LT0__annualInit(); // (selected variables)
                        routing.list_NaN__annualInit();
                        // monthly data - check & throw an exception
                        routing.check_LT0__annualInit();
                        routing.check_NaN__annualInit();

                        // Check daily output arrays used/initialized in "annualInit", e.g. total storage, cell AET
                        // daily data - check & list OR throw an exception
                        if ((5 == options.grid_store) || (6 == options.grid_store)) {
                            // can be negative = exclude
                            // only include listing if necessary
                            // routing.list_LT0__365TWS();
                            routing.list_LT0__365AET(); // only a small number of cells are negative (because discharge is larger than precipitation)

                            routing.list_NaN__365TWS();
                            routing.list_NaN__365AET();

                            // can be negative = exclude
                            // routing.check_LT0__365TWS();
                            // routing.check_LT0__365AET();
                            routing.check_NaN__365TWS();
                            routing.check_NaN__365AET();
                        }

                    }
                    // endif if (year >= chk_year)
                    // END (2.2) Evaluation year (after writing)
                    // END Check Felix Portmann 2015
#endif

                }
                // endif (actual_year >= options.evalStartYear)
            }
            //End: for (actual_year = start_year; actual_year <= end_year; actual_year += options.time_step)

            if (calStatNum > 0) {
                // this is a calibration run
                firstSimLoopFlag = false;

                if (testRunFlag) {
                    // if this was the test run, we are ready after writing the results
                    cout << "now, search for CFS " << endl;
                    calibGamma.writeCorrFactors(gamma, calibGamma.cellCorrFactorInd);
                    calibGamma.writeCalibStatus(calibGamma.getCalibStatus());
                    lastSimLoopFlag = true;
                } else {
                    // continue with model execution if gamma has not reached any limit
                    // otherwise this has been the last run
                    float gamma_old;

                    gamma_old = gamma;
                    cout << "now, I search for new gamma " << endl;
                    gamma = calibGamma.findNewGamma(gamma);
                    if (gamma < 0) {
                        // If return value of gamma was -99 the calibration
                        // is finished. Another run with the previous value of
                        // gamma is done to test the behaviour of the model.
                        testRunFlag = true;
                        gamma = gamma_old;
                    }
                }
            }
            // endif (calStatNum > 0)

        } while (!lastSimLoopFlag);
        // end of calibration loop
        // */ // END EXCLUDE & END LOOP
        /*** Clean up ***/
        climate.cleanup();
        // Clean memory for annual data
        if (options.time_series == 1){
            climateYear.cleanup();
            cout << "climateYear.cleanup()" << endl;
        }
        if ( options.clclOpt){
            clcl.cleanup();
            //    cout << "clcl.cleanup()" << endl;
        }
        if ( options.glacierOpt){
            glacierYear.cleanup();
            cout << "glacierYear.cleanup()" << endl;
        }
        if (progName == "OL") {
            // OE 2018
            cout << month <<", " << actual_year << endl;
            cout << configFile->endMonth <<", " << configFile->endYear << endl;
            if (month == configFile->endMonth && actual_year-1 == configFile->endYear) {
                delete wghmState;
                wghmState=NULL;
                std::cout << "die Adresse vom freigegebenen Speicher fuer wghmState ist " << wghmState << ">" << std::endl;
                delete calParam;
                calParam=NULL;
                std::cout << "die Adresse vom freigegebenen Speicher fuer calParam ist " << calParam << ">" << std::endl;
                delete additionalOutIn;
                additionalOutIn=NULL;
                std::cout << "die Adresse vom freigegebenen Speicher fuer additionalOutIn ist " <<     additionalOutIn << ">" << std::endl;
                delete snow_in_elevation;
                snow_in_elevation=NULL;
                std::cout << "die Adresse vom freigegebenen Speicher fuer snow_in_elevation ist " << snow_in_elevation << ">" << std::endl;
                delete configFile;
                configFile=NULL;
                std::cout << "die Adresse vom freigegebenen Speicher fuer configFile ist " << configFile << ">" << std::endl;
            }
        }
    }
    catch(std::exception &e)
    {
        throw(Exception(std::string("In watergap::main():\n")+e.what()));
    }
    
    return;
}

void integrate_wghm(const char* s, ConfigFile *& configFile,WghmStateFile *& wghmState, calibParamClass *& calParam, AdditionalOutputInputFile *& additionalOutIn, SnowInElevationFile *& snow_in_elevation,long *step, long *total_steps,long *year, long *month,const char* s2){
    
    integrate_wghm_(s,configFile,wghmState,calParam, additionalOutIn,snow_in_elevation,step,total_steps,year, month,s2); //keine funkition wird hier definiert, sondern nur aufgerufen
}



