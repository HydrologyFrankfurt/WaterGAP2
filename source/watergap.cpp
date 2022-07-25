/**********************************************************************
* Modification (MH20160601N001):
* The arguments of the program was not working properly. During this
* modification the bug is fixed temporarily (in a rather orthodox manner).
*
* This modification is required for the parallel implementation of model
* calibration using Borg-MOEA 1.8.
*
* Modified by: H.M. Mehedi Hasan, Dated on: June 01, 2016
***********************************************************************/
/**********************************************************************
* Modifications by: Felix T. Portmann (FP), Institute of Physical Geography, Goethe University Frankfurt am Main
*
* FP20160915N001:
* Arguments of program from update MH20160601N001:
* Now in error message also arguments of single cells output are included
* (watergap.cpp, options.cpp)
*
* FP20160915N002: scoutSurfaceRunoff -> options.scoutCellSRunoff
* Correction of initialization of G_daily365CellSurfaceRunoff: options.scoutCellSRunoff instead of erroneous scoutSurfaceRunoff
* (routing.cpp)
*
* FP20160926N001: changed for daily output .31 and .365 arrays: replaced "+=" by consistent "="
* for several variables e.g. river availability, transported volume, river velocity
* (routing.cpp)
*
* FP20161005N001: double G_landAreaFracPrevTimestep__at__n = routing.G_landAreaFracPrevTimestep[n];
* (daily.cpp)
*
* FP20161010N001 Replacement: geo.G_contfreq[n] instead (geo.G_landfreq[n] + geo.G_fwaterfreq[n])
* Class geo.writeContFreq for writing resulting grid to file GCONTFREQ.UNF0
* (geo.cpp/h, daily.cpp, routing.cpp, watergap.cpp)
*
* FP20161010N002: remove internal selection "if (calibrun == 0)"; now selection before calling "if (geo.calibrun == 0)"
* (geo.cpp)
*
* FP20161011N001: Introduce monthly rainfall output (as monthly snowfall output is also written)
* (daily.cpp)
*
* FP20161012N001 Enabling treatment of maximum slope class in G_SLOPE_CLASS.UNF1 (currently in file max. 69)
* FP20161012N002 Corrected error treatment (for values outside class limits) for G_TEXTURE.UNF1
* (gw_frac.cpp)
*
* FP20161018N001 Evaluation of arguments since first argument
* FP20161018N002 Reservoir operation start years
* FP20161018N003 Enable reading input data from climate land mask
* FP20161018N004 Net use / net abstraction for reservoir algorithm: new filename
* (climate.cpp, def.h, geo.cpp, option.h/cpp, routing.h/cpp, watergap.cpp)
***********************************************************************/

/***********************************************************************
 *Since 2.1f (around the year 2006) the development of WaterGAP (version with
 *spatial resolution of 0.5°) at CESR and IPG resulted in different model
 *versions and features. At the WaterGAP coordination meeting in 2011 it was
 *decided to create a harmonized model version of WaterGAP 2 integrating all
 *current developments. This model version (WaterGAP 2.2 refvers) is the basis
 *for further code development. The old headers (version control infos) were
 *saved in a textfile and replaced by comments on changes since 2.1f. Names of
 *persons who contributed to the code changes were added in order to get the
 *right person in touch for questions. Where possible, references were added.
 *Within the most of the .cpp headers, comments were added (where the changes
 *occur). Here, a short summary is done. More detailled documentation is in
 *preparation.
 *
 * WaterGAP 2.2 refvers has a lot of new features, including: handling of daily
 * climate input (as days per month and as days per year); variable flow
 * velocity including bankfull flow; water use can be taken out of surface
 * water or ground water sources (water_use.cpp is the basis for the
 * preprocessing tool GWSWUSE); a new land cover classification scheme (IGBP)
 * and data basis (combination of GLCT & Corine) as well as revised land use
 * attributes are used; arid/humid definition can be calculated at Köppen
 * climates; several new radiation (e.g. land use class based emissivity value
 * for calculating outgoing long wave radiation) and PET options (e.g. crop
 * coefficients) were included; distribution of permafrost can be calculated;
 * the growing period is based on precipitation and not on PET anymore;
 * precipitation correction is done after Adam & Lettenmaier; some user friedly
 * things were added (e.g. MODEL_SETTINGS_OUT.DAT), the original sensitivity
 * analysis was deleted, several input files were updated and many output files
 * were added.
 *
 *WaterGAP 2.2 refvers was calibrated within CESR at 1xxx GRDC stations using
 *WATCH Forcing Data as climate input; IPG used 1323 GRDC stations and CRU 3.1
 *climate input (with GPCC v6 precip). The calibration period was set to
 *1971-2000 (if discharge data were available at this time span).
 *
 * Stephanie Eisner & Hannes Müller Schmied, Kassel and Frankfurt in June 2012
 *
 * see former changes at file watergap.cpp.versioninfos.txt
 *
 ***********************************************************************/
// WaterGAP
// Water - Global Assessment and Prognosis
//
// original developed by Frank Kaspar, Kassel, 1997-2000
// Center for Environmental Systems Reserach
// University of Kassel
// kaspar@usf.uni-kassel.de

// since that the following people were involved in code development and
//provided model input files (alphabetical order): Adam, Linda; Aus der Beek,
//Tim; Döll, Petra; Eisner, Stephanie; Fiedler (now Wilkinson), Kristina; Flörke,
//Martina; Hoffmann-Dobrev, Heike; Hunger, Martin; Müller Schmied, Hannes = HMS;
//Portmann, Felix = FP; Schneider, Christof; Schulze (now Verzano), Kerstin; Siebert,
//Stefan; Teichner (now Kynast), Ellen; Voß, Frank; Weiß, Martina; Zhang, Jing (tbc.)

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>

#ifdef _OPENMP
#include <omp.h>
#endif
#include <sys/time.h>
using namespace std;

//#include "mmgr.h"

#include "gridio.h"
#include "daily.h"
#include "rout_prepare.h"
#include "calib_basins.h"
#include "upstream_stations.h"
#include "routing.h"
#include "option.h"
#include "gw_frac.h"
#include "s_max.h"
#include "lai.h"
#include "dailyprec.h"
#include "timestring.h"
#include "calibration.h"
#include "climate.h"
#include "climateYear.h"
#include "geo.h"
#include "land.h"
#include "permafrost.h"
#include "clcl.h"
#include "def.h"

// Definitions for global checks (all grid cells):
// variables (chk_year, chk_month, chk_day_in_month), classes/functions, exceptions
// Checks Felix Portmann 2015
#ifdef _WATERGAP_CHECKS_GLOBAL_H
#include "watergap_checks_global.h"
#endif

cbasinClass cbasin;
optionClass options;
laiClass lai;
dailyPrecNrdClass dailyPrecNrd;
routingClass routing;
climateClass climate;
climateYearClass climateYear;
upstreamStationClass directUpstSt;
upstreamStationClass allUpstSt;
calibGammaClass calibGamma;
dailyWaterBalanceClass dailyWaterBalance;
gridioClass gridIO;
geoClass geo;
landClass land;
soilWatCapClass maxSoilWaterCap;
groundwaterFactorClass GW;
clclClass clcl;

// variables which are used as externals
// in other routines

// #######################################
// Section to treat calibration parameters
// FP 2015-10

#include "calib_param.h"
calibParamClass calibParam;

// Conversion factors  // FP
#include "conversion.h"

//
// Important settings for file with calibration parameters
// e.g. for selected descriptors (see also code where other descriptors are set)
//
//int n_descriptors = 13; // number of descriptors, e.g. datetime (string & number (e.g. int))
//int n_ordinators = 2; // number of ordinators, e..g arcid, gcrc (number (int))
//int n_parameters = 26; // number of calibration parameters (arrays of number (double))

//short n_descriptors = 13; // number of descriptors, e.g. datetime (string & number (e.g. int))
//short n_ordinators = 2; // number of ordinators, e..g arcid, gcrc (number (int))
//short n_parameters = 26; // number of calibration parameters (arrays of number (double))

// directories as strings
std::string str_inputdir_global; // e.g. "/home/temp140/nobackup/portmann/WaterGAP2_2/JSON/JSON_WG_param2jsonfile"
std::string str_filename_global; // e.g. "/home/temp140/nobackup/portmann/WaterGAP2_2/JSON/JSON_WG_param2jsonfile"
std::string str_filename_calibParam; // e.g. "/home/temp140/nobackup/portmann/WaterGAP2_2/JSON/JSON_WG_param2jsonfile"
std::string str_pathfilename_global; // e.g. "/home/temp140/nobackup/portmann/WaterGAP2_2/JSON/JSON_WG_param2jsonfile"
// #######################################


// input grids
//---------ARID GROUNDWATER recharge calculations--------
short G_aindex[ng];

short G_toBeCalculated[ng]; // 0: no calculation
short G_SingleCelltoBeCalculated[ng]; // HMS single cell output 2013-11-21

signed short G_sbasin[ng]; // superbasins


// for measuring wallclock execution time of watergap
double timediff(struct timeval tv2, struct timeval tv1) {
    return (double) (tv2.tv_sec - tv1.tv_sec) + ((double) (tv2.tv_usec	- tv1.tv_usec) / 1000000.0);
}
// Structures for holding current time for gettimeofday (sys/time.h)
// for calculating difference in time with timediff()
struct timeval tv1_global, tv2_global;

// Measuring CPU time
// http://stackoverflow.com/questions/6129029/type-of-clocks-per-sec/6129107#6129107
// 2015-11-04
#include <time.h>
// CPU cycles
clock_t clockt1_global, clockt2_global;
// Function to cakcukate CPU time in seconds
double CPUtimediff_sec(clock_t clockt_2, clock_t clockt_1) {
	return ((double) (clockt_2 - clockt_1)) / CLOCKS_PER_SEC;
}

// grids for temporal local use, e.g. for transforming data and subsequent writing
// Felix Portmann 2015
float G_float[ng];

// main WaterGAP program
int main(int argc, char* argv[]) {

    // For compatibility reasons (32/64-bit processors)
    // all "long" variables were changed to "int".
    // Assure that int size is 4 byte. Exit if not.
    if (sizeof(int) != 4) {
        cerr << "Size of 'int' (" << sizeof(int)
             << " bytes) is too small on this machine! 4 bytes expected. Exiting...\n"
             << endl;
        exit(-1);
    }

	//Check arguments
    //	cout << "argv[0]: '" << argv[0] << "'" << endl;
    //	cout << "argv[1]: '" << argv[1] << "'" << endl;
	// Analyzing parameters as written by AEMS of University of Bonn for WaterGAP


    // Modification Number: MH20160601N001 (part-1)

    /*
    //stopped code segment: start
    switch(argc)
    	{
    		// valid situation (1)
    		case 2:
                // cout << "argv[0]: '" << argv[0] << "'" << endl;
                // cout << "argv[1]: '" << argv[1] << "'" << endl;
    			str_filename_calibParam = argv[1];
                cout << "Working with calibration parameter file in INPUT directory: " << str_filename_calibParam << endl;
    			break;

            // error
    		default:
    			cerr << "USAGE: " << argv[0] << " <string inputfile_calibration_parameters_cellspecific.json e.g. calib_parameters_initial_cell.json> "<< endl;
    			return -1;
    	}
	//end of the stopped segment
	*/

    // inserted code block (start)
    bool succeed = false;
    // Looking for -p argument if found use the given filename for parameters.json
    for (int i = 1; i < argc - 1; i++){
        if (!strcmp(argv[i], "-p")){
            str_filename_calibParam = argv[i+1];
            succeed = true;
            break;
        }
    }
    // If there is no -p argument given standard filename is assumed
    if (!succeed){
            str_filename_calibParam = "parameters.json";
        }

    cout << "Working with calibration parameter file in INPUT directory: " << str_filename_calibParam << endl;
	// new code block (end)

	// end of part-1 of MH20160601N001

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

    char filename[250];
    char filename_no[250];
    char string[250];

    short day, day_in_month, month, actual_year; // counter for time steps
    // number_of_days_in_month[12] moved here and first_day_in_month[12] introduced
    const short number_of_days_in_month[12] = { 31, 28, 31, 30, 31,	30, 31, 31, 30, 31, 30, 31 };
    const short first_day_in_month[12] = {0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334};
    short totalInitYears, remainingInitYears;
    short nSpecBasins;	// number of specified basins (in STATIONS.DAT)
    short nSingleCells; // HMS 2013-11-21 number of secified single cells (in SINGLECELLS.DAT)
    const int anzThreads = 8;

#ifdef _OPENMP
    omp_set_num_threads(anzThreads);
    cout << anzThreads << " Threads" << endl;
#endif
	// Start time
	gettimeofday(&tv1_global, NULL);  // wallclock time
	clockt1_global = clock();  // CPU cycles

    // print information on simulation time
    // WaterGAP version number and compilation time.
    getISOdateTime(string, sizeof(string));
    cout << string << endl;
    cout << "This is WaterGAP 2.2b from" << endl;
    cout << __DATE__ << ' ' << __TIME__ << endl;
    cout << "--------------------------" << endl;

    //
    // initialize classes
    //

    //read OPTIONS.DAT
	//(or the specified file in command line) // 2nd item after parameter file
//	options.init(argc - 1, &argv[1]);  // original (wrong: cannot read switch + filename correctly)
//	options.init(argc - 2, &argv[2]);  // "original style" when argv[1] is name of parameterfile (wrong: cannot read switch + filename correctly)

	// Modification Number: MH20160601N001 (part-2)

	/*
	//stopped code block: start
	options.init(0, &argv[1]); // Just ignore any arguments, second argument may be void
    //end of stopped code block
    */

    //new code (start)
    // FP20161018N001 Evaluation of arguments since first argument
    // Do not extract again the parameter file, as this has been done before, possibly without option switch "-p")
    if (argc > 2) options.init(argc, argv);
    else options.init(0, NULL);
    //end of new code

    // end of part-2 of MH20160601N001



	// compose string of input directory path
	str_inputdir_global = options.input_dir;
//	cout << "CHECK str_inputdir_global '" << str_inputdir_global << "'" << endl;




	// Define and read calibration parameters
	cout << "Define and read calibration parameters from JSON file" << endl;

	// Establish vectors of names of
	//	calibration parameters (as used in file, in order to query them)
	// cout << "Define calibration parameters - calibParam.defineVNamesCalibParamJson()" << endl;
	calibParam.defineVNamesCalibParamJson();

	//	cout << "Read file - calibParam.readJson(str_pathfilename)" << endl;
	// Includes a global check whether a JSON object has been generated (program continuation) or not (program exit)
	str_pathfilename_global = str_inputdir_global + "/" + str_filename_calibParam; // path + file
//	cout << "CHECK str_pathfilename_global '" << str_pathfilename_global << "'" << endl;
	calibParam.readJson(str_pathfilename_global);

//	cout << "If desired, check read-in & parsed JSON object and print first array values - calibParam.callCheckJson()" << endl;
//	calibParam.callCheckJson();

//	cout << "If desired, check iteration time on all values of all parameters in JSON object- calibParam.checkJsonSpeed()" << endl;
//	calibParam.checkJsonSpeed();

	//	cout << "If desired, print values of selected grid cells - calibParam.printValueAllParam()" << endl;
	//	calibParam.printValueAllParam(0);
	//	calibParam.printValueAllParam(1);
	//	calibParam.printValueAllParam(ng-1);
	// CHECK Call of method
//		for (int n = 0; n < 2; n++) {
//			cout << "n = " << n << endl;
//			cout << calibParam.getValue(P_GAMRUN_C,n) << endl;
//			cout << calibParam.getValue(P_CFA,n) << endl;
//			cout << calibParam.getValue(P_CFS,n) << endl;
//			cout << calibParam.getValue(M_ROOT_D,n) << endl;
//			cout << calibParam.getValue(M_RIVRGH_C,n) << endl;
//			cout << calibParam.getValue(P_LAK_D,n) << endl;
//			cout << calibParam.getValue(P_WET_D,n) << endl;
//			cout << calibParam.getValue(P_SWOUTF_C,n) << endl;
//			cout << calibParam.getValue(M_EVAREDEX,n) << endl;
//			cout << calibParam.getValue(M_NETRAD,n) << endl;
//			cout << calibParam.getValue(P_PTC_HUM,n) << endl;
//			cout << calibParam.getValue(P_PTC_ARI,n) << endl;
//			cout << calibParam.getValue(P_PET_MXDY,n) << endl;
//			cout << calibParam.getValue(P_MCWH,n) << endl;
//			cout << calibParam.getValue(M_LAI,n) << endl;
//			cout << calibParam.getValue(P_T_SNOWFZ,n) << endl;
//			cout << calibParam.getValue(P_T_SNOWMT,n) << endl;
//			cout << calibParam.getValue(M_DEGDAY_F,n) << endl;
//			cout << calibParam.getValue(P_T_GRADNT,n) << endl;
//			cout << calibParam.getValue(M_GW_F,n) << endl;
//			cout << calibParam.getValue(M_RG_MAX,n) << endl;
//			cout << calibParam.getValue(P_PCRITGWA,n) << endl;
//			cout << calibParam.getValue(P_GWOUTF_C,n) << endl;
//			cout << calibParam.getValue(M_NETABSSW,n) << endl;
//			cout << calibParam.getValue(M_NETABSGW,n) << endl;
//			cout << calibParam.getValue(M_PREC,n) << endl;
//		}
//	}

	// /*  // START EXCLUDE 1
    options.createModelSettingsFile(options.output_dir);
    options.closeModelSettingsFile();

    cout << "Endian type: " << gridIO.getEndianType() << '\n';
    gridIO.setFileEndianType(options.fileEndianType);

    if (1 == options.rout_prepare)
        prepare_routing_files(options.input_dir, options.routing_dir, 1, options.resOpt);
    geo.init(options.input_dir, options.fileEndianType);

    // init climate with local data
    climate.init();
    if ((options.time_series == 1) || (options.time_series == 3)){
        climateYear.init();
        cout << "climateYear.init()" << endl;
    }

    // init clcl-class (Koeppen-map)
    if (options.clclOpt)
        clcl.init();

    if (options.petOpt) {

        if (options.petOpt == 3) {
            climate.initHG();
        } else {
            // Penman-Monteith and ... requires windspeed and vapour pressure
            climate.initPM();
            if (options.petOpt == 2)
                land.initPM(options.input_dir);
        }

    }

    nSpecBasins = cbasin.prepare(options.input_dir, options.output_dir,	options.routing_dir);
    cout << "Number of specified basins: " << nSpecBasins << endl;
	// Initialize cell water balance
	// Global parameters at dailyWaterBalance.setParameters, not cell-specific calibration parameters
	dailyWaterBalance.init(options.input_dir, nSpecBasins);
    sprintf(filename, "%s/G_CALIB_BASIN.UNF2", options.output_dir);
    gridIO.readUnfFile(filename, ng, G_sbasin);

    directUpstSt.findDirectStations(options.routing_dir);
    if (options.outDirectUpStations) { // new output options 2.1f
        sprintf(filename, "%s/DIRECT_UPSTREAM_STATIONS.DAT", options.output_dir);
        sprintf(filename_no, "%s/NO_DIRECT_UPSTREAM_STATIONS.DAT", options.output_dir);
        directUpstSt.writeListToFile(filename, filename_no);
    }
    allUpstSt.findAllStations(options.routing_dir);
    if (options.outAllUpStations) { // new output options 2.1f
        sprintf(filename, "%s/ALL_UPSTREAM_STATIONS.DAT", options.output_dir);
        sprintf(filename_no, "%s/NO_ALL_UPSTREAM_STATIONS.DAT",	options.output_dir);
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
    int *G_outflow_cell = 0;
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

    dailyPrecNrd.init(raindaysRandomInit, options.output_dir);

    // read correction factors for precipitation (Adam and Lettenmaier)
    // (percentage)
    if (options.prec_correct == 1) {
        climate.readPrecAdjustFactors(options.input_dir, options.climate_dir);
    }

    if (6 == options.time_series) { // new time series input option
        climate.createTempAvg(options.climate_dir, options.evalStartYear, options.end_year);
        climate.createPrecAvg(options.climate_dir, options.evalStartYear, options.end_year);
        if (1 == options.prec_correct)
            climate.adjustPrecipitation();
        if (0 == options.cloud)
            climate.createSunshineAvg(options.climate_dir, options.evalStartYear, options.end_year);
        if (1 == options.cloud)
            climate.createCloudAvg(options.climate_dir, options.evalStartYear, options.end_year);
        if (2 == options.cloud)
            climate.createRadiationAvg(options.climate_dir,	options.evalStartYear, options.end_year);
        climate.createRaindaysAvg(options.climate_dir, options.evalStartYear, options.end_year);
    }
    //
    // read grid files which are independent of time
    //

    land.init(options.input_dir); // requires GCRC_2

    //-------GROUNDWATER RECHARGE (read (semi)arid cells)-----

    //dailyWaterBalance.readAridCells();
    sprintf(filename, "%s/G_ARID_HUMID.UNF2", options.input_dir); // renamed from AINDEX
    gridIO.readUnfFile(filename, ng, G_aindex);
    //cerr << "Unable to open file:" << filename << endl;
    //cout << " finished reading arid cells " << endl;

    //================== ELEVATION =====================================
    //dailyWaterBalance.readElevationFile();
    sprintf(filename, "%s/G_ELEV_RANGE.101.UNF2", options.input_dir);
    //gridIO.readUnfFile(filename, ng_land * 101,	&dailyWaterBalance.G_Elevation[0][0]);
    gridIO.readUnfFile(filename, ng * 101,	&dailyWaterBalance.G_Elevation[0][0]);
    //=============== END ELEVATION ====================================

    lai.init(options.input_dir, options.output_dir, land.G_landCover);


    //routing.createTestOutFile(options.output_dir);
	// Initializing elements necessary for routing,
	// Global parameters in routing.setParameters, but not cell-specific calibration parameters
    routing.init(nSpecBasins);

    if (calStatNum > 0) {
        G_outflow_cell = new int[ng];

        sprintf(filename, "%s/G_OUTFLC.UNF4", options.routing_dir);
        gridIO.readUnfFile(filename, ng, G_outflow_cell);

        if (6 == options.time_series) // new time series input option
            options.end_year = options.evalStartYear;
    }

    if (1 == options.day_store)
        routing.createDailyFiles(options.output_dir);

    if (2 == options.day_store) {    // HMS 2013-11-21 single cell output
        nSingleCells = options.createSingleCellFileList(options.input_dir, options.output_dir);

        cout << "Number of single cells where daily output is created: " << nSingleCells << endl;
    }

    routing.createMonthlyDischargeFile(options.output_dir);
    routing.createAnnualDischargeFile(options.output_dir);
//*/  // END EXCLUDE 1

// /*  // START EXCLUDE & END LOOP
	// Initialize cell-specific active depth of lakes and wetlands (typically from calibration parameters)  // FP
	routing.initLakeDepthActive();
	routing.initWetlDepthActive();

    // start of the loop for calibration or sensitivity analysis
    // if a normal simulation is executed the following do-while-loop is
    // only executed once
    do { //while (!lastSimLoopFlag)
        // nb=waterbasins(options.input_dir);
        routing.initFractionStatus(); //avoid taking last landareafrac vom previous calibration loop by setting statusStarted_updateGloResPrevYear and statusStarted_landfreq_fwaterfreq_landAreaFrac to 0.
        routing.setStoragesToZero();
        routing.setLakeWetlToMaximum(options.start_year);

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
			sprintf(filename, "%s/G_GAMMA_HBV.UNF0", options.input_dir);
			gridIO.readUnfFile(filename, ng, dailyWaterBalance.G_gammaHBV);
		}
		// end if (calStatNum > 0)
		else {
			// Normal runs: transfer values from JSON-object (from calibration parameter file) // FP
			for (int n = 0; n < ng; n++) {
					dailyWaterBalance.G_gammaHBV[n] = calibParam.getValue(P_GAMRUN_C,n);
			}
		}
		// end else filling G_gammaHBV (P_GAMRUN_C)

		// Filling of dailyWaterBalance.G_cellCorrFact (P_CFA)
		// Filling of routing.G_statCorrFact (P_CFS)
		// Read from G_CORR_FACTOR.UNF0 & G_STAT_CORR.UNF0 only during calibration runs
        if (calStatNum > 0) {
            if (!testRunFlag) {
                sprintf(filename, "G_CORR_FACTOR.UNF0");
                gridIO.readUnfFile(filename, ng, dailyWaterBalance.G_cellCorrFact);
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
			// Obsolete reading of files
			//			sprintf(filename, "%s/G_CORR_FACTOR.UNF0", options.input_dir);
			//			gridIO.readUnfFile(filename, ng, dailyWaterBalance.G_cellCorrFact);
			//			sprintf(filename, "%s/G_STAT_CORR.UNF0", options.input_dir);
			//			gridIO.readUnfFile(filename, ng, routing.G_statCorrFact);
			// Normal runs: transfer values from JSON-object (from calibration parameter file) // FP
			for (int n = 0; n < ng; n++) {
					dailyWaterBalance.G_cellCorrFact[n] = calibParam.getValue(P_CFA,n);
					routing.G_statCorrFact[n] = calibParam.getValue(P_CFS,n);
			}
		}
		// end else filling G_cellCorrFact (P_CFA) and G_statCorrFact (P_CFS)


        //replace gammas (if some are given in the STATIONS-file)
        bool gammaInStationList = false;

        for (short n = 1; n <= cbasin.numberOfCalibBasins; n++) {
            if (cbasin.gamma[n - 1] > 0) {
                gammaInStationList = true;
                allUpstSt.replaceGridValues(n, cbasin.gamma[n - 1],	dailyWaterBalance.G_gammaHBV, cbasin.gamma);
                //allUpstSt.replaceGridValues(n,
                //              cbasin.cellCorrFactor[n-1],
                //              dailyWaterBalance.G_cellCorrFact,
                //              cbasin.cellCorrFactor);
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
            sprintf(filename, "G_GAMMA_HBV_NEW.UNF0");
            gridIO.writeUnfFile(filename, ng, dailyWaterBalance.G_gammaHBV);
            sprintf(filename, "G_STAT_CORR.UNF0");
            gridIO.writeUnfFile(filename, ng, routing.G_statCorrFact);
        }

        dailyWaterBalance.setStoragesToZero();

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
                            //cout << "Cell: " << m << endl;
                            dailyWaterBalance.G_gammaHBV[m] = cbasin.gamma[upSt	- 1];
                            //dailyWaterBalance.G_cellCorrFact[m]
                            //= cbasin.cellCorrFactor[upSt-1];
                            G_toBeCalculated[m] = 1;
                        }
                } else {
                    cerr << "Gamma not available\n";
                    exit(-1);
                }
            }
            G_toBeCalculated[G_outflow_cell[cbasin.cellNum[calStatNum - 1] - 1]	- 1] = -1;
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
        maxSoilWaterCap.createMaxSoilWaterCapacityGrid(options.input_dir, options.output_dir, land.G_landCover);
        GW.createGrids(options.input_dir, options.output_dir);

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
        else if (1 == options.climate_spatial_resolution){
            climate.read_climate_longtermAvg_withConvert();
        }

        //===============KOEPPEN CLIMATE CLASSIFICATION===========================
        if (options.clclOpt) {

            sprintf(filename, "%s/G_KOEPPEN.UNF2", options.input_dir);
            ifstream station_file(filename);

            if (!station_file) {

                cout << "Calculate KOEPPEN map... " << endl;

                if (options.time_series != 5) { // new time series input option
                    // read long term mean data
                    // monthly values [0.01 oC]
                    sprintf(filename, "%s/GTEMP_1971_2000.12.UNF2",	options.climate_dir);
                    gridIO.readUnfFile(filename, ng * 12, &climate.G_temperature[0][0]);

                    // read monthly precipitation [mm]
                    sprintf(filename, "%s/GPREC_1971_2000.12.UNF2",	options.climate_dir);
                    gridIO.readUnfFile(filename, ng * 12, &climate.G_precipitation[0][0]);
                    if ((1 == options.prec_correct) && ((4 == options.time_series) ||
                                                        (5 == options.time_series) || (6 == options.time_series))) // new time series input option prevent precip correction if daily precip is used.
                        climate.adjustPrecipitation();
                }

                // calculate climate zones
                clcl.calcKoep();
                // assign alpha values
                clcl.alphaKoep();

                // write KOEPPEN map (grid)
                if (options.outClcl) {
                    char filename[250];
                    sprintf(filename, "%s/G_KOEPPEN.UNF2", options.output_dir);
                    gridIO.writeUnfFile(filename, ng, clcl.cls);
                }
            } else {
                station_file.close();
                cout << "Read KOEPPEN map from UNF-file... " << endl;
                // read climate zones from file
                gridIO.readUnfFile(filename, ng, &clcl.cls[0]);
                // assign alpha values
                clcl.alphaKoep();
            }

        }
        //===============KOEPPEN CLIMATE CLASSIFICATION===========================

        // start of loop for simulation years
        // calculate from start to end year
        short end_year = -1, start_year = -1;

        if (options.time_series <= 4) { // new time series input option
            start_year = options.start_year;
            end_year = options.end_year;
        }
        if (options.time_series > 4) { // new time series input option
            end_year = options.evalStartYear;
            start_year = options.evalStartYear;
        }

        // for variable land cover map
        int landCoverYear = 0;
        if (options.landCoverOpt) {
            land.readLandCoverMap(options.landCoverYears[landCoverYear], options.land_cover_dir);
            landCoverYear++;
        }

        totalInitYears = options.init_years;

#ifdef _WATERGAP_CHECKS_GLOBAL_H
        cout << " CHECK totalInitYears (from options): " << totalInitYears << endl; // Check Felix Portmann 2015
#endif

        // Loop over years
        for (actual_year = start_year; actual_year <= end_year; actual_year	+= options.time_step) {
            cout << "Year of simulation: " << actual_year << endl;
            // for variable land cover map
            if (options.landCoverOpt && // use variable land cover map
                    landCoverYear < options.landCoverYearsInList && // last map is not read yet
                    actual_year >= options.landCoverYearsStart[landCoverYear]){ // actual_year use next map

                land.readLandCoverMap(options.landCoverYears[landCoverYear], options.land_cover_dir);
                landCoverYear++;
            }

            //
            // read grid files for monthly climate data which change in time
            //
            if ((2 == options.time_series) || (3 == options.time_series) || (4	== options.time_series)) { // new time series input option
                if (0 == options.climate_spatial_resolution)
                    climate.read_climate_data_monthly(actual_year, 12);
                else
                    climate.read_climate_data_monthly_withConvert(actual_year, 12);
            }


            // and daily climate data from yearly files
            // new time series input option start - add .365
            if ((options.time_series == 1) || (options.time_series == 3)) {
                if (0 == options.climate_spatial_resolution) {
                    cout << "Begin read_climate_data_daily_per_year(" << actual_year << ", 365)" << endl;
                    climateYear.read_climate_data_daily_per_year(actual_year, 365);
                    cout << "End read_climate_data_daily_per_year(" << actual_year << ", 365)" << endl;
                } else {
                    cerr << "reading in .365 UNF at climate landmask is not yet implemented." << endl;
                    //cout << "Begin read_climate_data_daily_per_year(" << actual_year << ", 365)" << endl;
                    //climateYear.read_climate_data_daily_per_year_withConvert(actual_year, 365);
                    //cout << "End read_climate_data_daily_per_year(" << actual_year << ", 365)" << endl;
                }
            }
            // new time series input option end - add .365

            //===============Permafrost Calculations===========================

            // with daily climate data for permafrost calculating we need monthly values for temperature and precipitation
            if (options.time_series <= 3 && options.permaOpt) {
                for (int month = 0; month < 12; ++month) {
                    // read daily temperature and precipitation from .31 files
                    if (options.time_series == 0 || options.time_series == 2)
                        if (0==options.climate_spatial_resolution) {
                            climate.read_climate_data_daily(month+1, actual_year, 31);
                        } else {
                            climate.read_climate_data_daily_withConvert(month+1, actual_year, 31);
                        }
                    //or split daily data from .365 files
                    else if (options.time_series == 1 || options.time_series == 3)
                        climateYear.split_climate_data_from_year_to_daily(first_day_in_month[month], number_of_days_in_month[month]);

                    double dummy_temp, dummy_prec;
                    for (int n = 0; n < ng; ++n) {
                        dummy_temp = 0.; dummy_prec = 0.;
                        for (int day_in_month = 1; 	day_in_month <= number_of_days_in_month[month]; day_in_month++) {
                            dummy_temp   += climate.G_temperature_d[n][day_in_month];
                            dummy_prec += climate.G_precipitation_d[n][day_in_month];
                        } // for(day_in_month)
                        climate.G_temperature[n][month]  = (short)(dummy_temp/number_of_days_in_month[month] *100.);
                        climate.G_precipitation[n][month]= (short)(dummy_prec);
                    } // for(n)
                }//for(month)
            }

            if (options.permaOpt) {

                permaClass permafrost;

#pragma omp parallel for  \
    shared(actual_year,cout,permafrost)

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

            } // end of permaOpt
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

                // read dailyNUs und dailyNUg or reallocated values
                if (((2 == options.subtract_use)  || (3 == options.subtract_use)) && ((4 == options.time_series) || (0 == options.time_series))) {
                    routing.dailyNUInit(options.water_use_dir, actual_year);
                    //cout << "Initialization of dailyNUInit for year " << actual_year << " successfully done." << endl;
                }

                dailyWaterBalance.annualInit();
                routing.annualInit(actual_year); // FP20161018N002: Do NOT introduce nSpecBasins, as superbasin area calculation probably wrong & obsolete

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
                for (month = 0; month < 12; ++month) {
                    //reformulated for daily output options (.31, .365)
                    if ((3 == options.grid_store) || (4 == options.grid_store)){
                        dailyWaterBalance.daily31outInit();
                        routing.daily31outInit();
                    }
                    //
                    // read grid files for daily climate data which change in time
                    //

                    // new time series input option start
                    if ((options.time_series == 0) || (options.time_series == 2)) {
                        if (0==options.climate_spatial_resolution)
                            climate.read_climate_data_daily(month + 1, actual_year,	31);
                        else
                            climate.read_climate_data_daily_withConvert(month + 1, actual_year, 31);
                    }
                    if ((options.time_series == 1) || (options.time_series == 3)) {
                        climateYear.split_climate_data_from_year_to_daily(day, number_of_days_in_month[month]);
                    }
                    // new time series input option end

                    for (day_in_month = 1; day_in_month	<= number_of_days_in_month[month]; day_in_month++) {

                        day++; // [1..365]

                        //		cout <<"Year, Month, Day : "<< actual_year<<", "<<month << ", " << day << endl;

#pragma omp parallel for  \
    shared(day,month,day_in_month,actual_year,cout,dailyWaterBalance)
                        for (int n = 0; n < ng; ++n) {
                            // FP 2015-06
                            // Only execute for continental grid cells
                            if (geo.G_contcell[n]) {
#ifdef _OPENMP
                                if (omp_get_num_threads() != anzThreads) {
                                    cout << "Falsche Threadzahl" << endl;
                                    exit(1);
                                }
#endif
                                dailyWaterBalance.calcNewDay(day, month, day_in_month, actual_year, n);

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

                        if ((2 == options.subtract_use) || (3 == options.subtract_use)){ // enhanced to reallocated water use
                            routing.calcNextDay_M(day);
                        }
                        routing.routing(actual_year, day, month, day_in_month); // routing through river network
                        routing.updateLandAreaFrac();
                    } // end of loop for day_in_month

                    // FP 2015-06
                    // In the following methods, it has to be guaranteed that
                    //  - only continental grid cells are processed
                    //  - or zero values are still maintained from initialization
                    //      or written, if necessary
                    if ((actual_year >= options.evalStartYear) && (remainingInitYears == 0)) {

                        // daily output option 31 start
                        if ((3 == options.grid_store) || (4	== options.grid_store)) {
                            routing.writeDaily31Grids(options.output_dir, actual_year, month);
                            dailyWaterBalance.writeDaily31Grids(options.output_dir, actual_year, month);
                        }
                        // daily output option 31 end
                    }

                } // end of loop for month


                // FP20161018N002 Reservoir operation start years
                // Transfer land area fraction of reservoirs (in percent) of current year in array of previous year for use in the next year
                routing.updateGloResPrevYear_pct();


                // FP20161018N002 Reservoir operation start years
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

                if (0 != options.subtract_use)
                    routing.annualWaterUsePostProcessing(actual_year);

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
                        calibGamma.setWaterUse(actual_year,	routing.getAnnualSatisfWaterUse(calStatNum));
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

            }	// end of loop for remainingInitYear

            totalInitYears = 0;

			// Superbasin calculation currently obsolete, probably not working correctly in WaterGAP2.2b, do in postprocessing  // FP 2015-10
//            if (1 == options.day_store)
//                dailyWaterBalance.storeDailyValuesToFile(actual_year);
            if (2 == options.day_store)
                routing.storeSingleCellDailyValuesToFile(actual_year); // HMS 2013-11-21 write values into file
			// Superbasin calculation currently obsolete, probably not working correctly in WaterGAP2.2b, do in postprocessing  // FP 2015-10
//			if (options.outSuperbasinClimate)  // new output options 2.1f  // Currently obsolete, probably not working correctly in WaterGAP2.2b, do in postprocessing  // FP 2015
//				dailyWaterBalance.storeAnnualValuesToFile(actual_year);  // Currently obsolete, probably not working correctly in WaterGAP2.2b, do in postprocessing  // FP 2015

            if (actual_year >= options.evalStartYear) {

                // write output grids
                if (options.grid_store != 0) {
                    if ((1 == options.grid_store) && (options.outSingleStorages)) {
						routing.writeLakeStorRatioGrids(options.output_dir, actual_year);
                        //dailyWaterBalance.writeAnnualGrids(options.output_dir, actual_year);
                        routing.writeLakeStorageGrids(options.output_dir, actual_year);
                    }
					dailyWaterBalance.writeMonthlyGrids(options.output_dir, actual_year);
                    routing.writeMonthlyGrids(options.output_dir, actual_year);

                }

                if ((options.grid_store != 0) && (calStatNum > 0)) {
                    // these grids are required to calculate
                    // the cell related correction factor
                    routing.writeTotalCellRunoffGrid(options.output_dir, actual_year);
                }

                // FP20161018N002 Reservoir operation start years
                if ((0 != options.grid_store) && (options.outLandWaterFractions )) {
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
                calibGamma.writeCorrFactors(gamma,calibGamma.cellCorrFactorInd);
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
    //	cout << "climate.cleanup()" << endl;
    // Clean memory for annual data
    if ((options.time_series == 1) || (options.time_series == 3)){
        climateYear.cleanup();
        cout << "climateYear.cleanup()" << endl;
    }
    if ( options.clclOpt){
        clcl.cleanup();
        //	cout << "clcl.cleanup()" << endl;
    }

	// Final time
	gettimeofday(&tv2_global, NULL);  // wallclock time
	clockt2_global = clock();  // CPU cycles
	cout << "Wallclock Time (seconds): = " << timediff(tv2_global, tv1_global) << endl;
	cout << "CPU Time (seconds): = " << CPUtimediff_sec(clockt2_global, clockt1_global) << endl;

}
// end main
