/***********************************************************************
 *
 *variable flow velocity, developed by Kerstin Verzano (former Schulze) and
 *Frank Voss as well as restructured routing scheme. This is one of the main
 *changes of WaterGAP 2.2 refvers compared to 2.1 versions. This is described
 *in the paper Verzano, K., I. Bärlund, M. Flörke, B. Lehner, E. Kynast, F.
 *Voß, and J. Alcamo (2012), Modeling variable river flow velocity on
 *continental scale: Current situation and climate change impacts in Europe,
 *Journal of Hydrology, 424-425, 238-251, doi:10.1016/j.jhydrol.2012.01.005.
 *[online] Available from:
 *http://linkinghub.elsevier.com/retrieve/pii/S0022169412000248 as well as in
 *the thesis Verzano, K. (2009), Climate Change Impacts on Flood Related
 *Hydrological Processes: Further Development and Application of a Global Scale
 *Hydrological Model, 166 pp., University of Kassel. Within that framework, the
 *bankfull flow concept was integrated by Cristof Schneider and described in the
 *paper Schneider, C., M. Flörke, S. Eisner, and F. Voss (2011), Large scale
 *modelling of bankfull flow: An example for Europe, Journal of Hydrology,
 *408(3-4), 235-245, doi:10.1016/j.jhydrol.2011.08.004. [online] Available
 *from: http://linkinghub.elsevier.com/retrieve/pii/S0022169411005300
 *
 *water use is now divided into water uptake from groundwater and surface
 *water. This was coded by Heike Hoffmann-Dobrev, a tool to calculate the net
 *uses of groundwater and surface water (GWSWUSE) was developed from Felix
 *Portmann and Heike Hoffman-Dobrev. The concept is described in the paper
 *Döll, P., Hoffmann-Dobrev, H., Portmann, F.T., Siebert, S., Eicker, A.,
 *Rodell, M., Strassberg, G., Scanlon, B.R. (2012): Impact of water withdrawals
 *from groundwater and surface water on continental water storage variations,
 *Journal of Geodynamics (in press).
 *
 *calculation of groundwater storages was moved from daily.cpp to routing.cpp
 *by Heike Hoffmann-Dobrev; storage otuputs were restrctured and enhanced by
 *Heike Hoffmann-Dobrev, Linda Adam, Stephanie Eisner and Hannes Müller
 *Schmied); several new outputs are integrated (by Ellen Kyast (former
 *Teichert) Frank Voss, Stephanie Eisner (e.g. max and min monhtly river
 *availbility), Hannes Müller Schmied (e.g. restructuring .31 and .365 output)
 *
 *new reservoir algorithm from Kristina Fiedler, described at the paper Döll,
 *P., K. Fiedler, and J. Zhang (2009), Global-scale analysis of river flow
 *alterations due to water withdrawals and reservoirs, Hydrology and Earth
 *System Sciences, 13(12), 2413-2432, doi:10.5194/hess-13-2413-2009. [online]
 *Available from: http://www.hydrol-earth-syst-sci.net/13/2413/2009/ as well as
 *work on that by Ellen Kynast (former Teichert)
 *
 *handling of .365 climate input by Heike Hoffmann-Dobrev
 *
 *litle inconsistencies due to numeric problems when reading in float values of
 *surface water body input files were catched throug a adapted storage number
 *of lake/wetland fraction error message (Linda Adam, Hannes Müller Schmied)
 *
 *all mean input (eg. mean temperature) were changed from 1961-1990 to
 *1971-2000 because this period is the focus period of calibration (Hannes
 *Müller Schmied)
 *
 *the allocated water use of surface water as preprocessing step for
 *calibration (where water use from neighbour cell is not possible at catchment
 *borders) is now G_ACTUAL_NAs (Stephanie Eisner, Hannes Müller Schmied)
 *
 *important: it was decided to do the calibration with station correction but
 *normal runs without station correction e.g. when the model time period
 *differs from calibration period. This can be switched on and of with the new
 *option 24 (Hannes Müller Schmied).
 *
 * see former changes at file routing.cpp.versioninfos.txt
 *
 ***********************************************************************/
//#define NDEBUG

// Felix Portmann (FP) 2015 & 2016 several changes
//  Saved in UTF-8 format
//  Inserted "// end"-statements
//      - for classes
//      - for/if/else statements (mostly class/method "routing")
//  Replacing all "++n" by "n++" (filling / plain treatment procedures)
//  Use of geo.G_contcell[n] (geo.cpp)
//      - With query (geo.G_landfreq[n] + geo.G_fwaterfreq[n]) < limit_contcell_pct)
//      instead of query ((geo.G_landfreq[n] + geo.G_fwaterfreq[n]) == 0) (or "== 0."
//      - For correctly identifying contintal cells vs. ocean cells of Caspian Sea
//  Moving some brackets "}" to correct position in class/method "routing"
//      for G_tobeCalculated[n]
//  Inserted checks functions for testing "less than zero" and "Not a Number" (NaN)
//  Using cell-specific Calibration Parameters
// FP20161010N001 Replacement: geo.G_contfreq[n] instead (geo.G_landfreq[n] + geo.G_fwaterfreq[n])
// FP20161018N002 Reservoir operation start years

#include <cstdio>
#include <cmath>
#include <cassert>
#include <fstream>
#include <iomanip>
#include <vector>

#include "daily.h"
#include "calib_basins.h"
#include "upstream_stations.h"
#include "option.h"
#include "gridio.h"
#include "geo.h"
#include "timestring.h"
#include "def.h"
#include "calibration.h"
#include "routing.h"

using namespace std;

extern optionClass options;
extern gridioClass gridIO;
extern geoClass geo;
extern cbasinClass cbasin;
// FP20161018N002 Reservoir operation start years
extern dailyWaterBalanceClass dailyWaterBalance;

// #######################################
// Cell-specific calibration parameters
// FP
#include "calib_param.h"

extern calibParamClass calibParam; // JSON object and methods from watergap.cpp
// #######################################
extern calibGammaClass calibGamma;
// Conversion factors  // FP
#include "conversion.h"


const unsigned char numberOfDaysInMonth[12] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};

template<class T>
inline int round(const T value) {
    return (int) floor(value + 0.5);
}

// Definitions for global checks (all grid cells):
// use of extern variables (chk_year, chk_month, chk_day_in_month), classes/functions, exceptions
// Felix Portmann 2015
#ifdef _WATERGAP_CHECKS_GLOBAL_H
#include "routing_checks_global.h"
#endif

routingClass::routingClass() {
    defaultRiverVelocity = 86.4; // [km/d] = 1 m/s
    defaultTimeStepsPerDay = 1;
}
// end routingClass

void routingClass::setParameters() {
//	lakeDepth = 0.005; // [km]  // OBSOLETE now read from calibration parameter file // FP
//	wetlandDepth = 0.002; // [km]  // OBSOLETE now read from calibration parameter file // FP
    lakeOutflowExp = 1.5;
    wetlOutflowExp = 2.5;
//	k_g = 0.01;	// for groundwater storage	// [1/d]  // OBSOLETE now read from calibration parameter file // FP
//calibParam.getValue(P_SWOUTF_C,n)
//	glo_storageFactor = 0.0125;  // [1/d]  // ADAPT to new value 0.01 for comparison of outputs! // OBSOLETE now read from calibration parameter file // FP
//	loc_storageFactor = 0.01;    // [1/d]  // OBSOLETE now read from calibration parameter file // FP

    evapoReductionExp = 3.32193;
    evapoReductionExpReservoir = 2.81383;
}
// end setParameters

void routingClass::init(const short nBasins) {
    // - reads the files with routing parameters and
    //   grids with lake and wetland locations
    // - calculates land area fraction
    // - reads files with flow direction, etc.
    // - dynamic memory allocation

    // requirements:
    // - G_sbasin has been read in 'main'
    char filename[250];
    // GFZ:GLOLAK begin - initialise and read variables of biglakefile (from Kristina)
    if (0 != ng_glolakcells) {
        for (int n = 0; n < ng_glolakcells; n++) {
            big_lakes_area_frac[n] = 0.0;
            big_lakes_cells[0][n] = -99;
            big_lakes_cells[1][n] = -99;
        }
        sprintf(filename, "%s/wg_big_lakes_area_frac.txt", options.input_dir);
        ifstream biglakefile(filename);
        if (biglakefile) {
            int n = 0;
            // float data; // FP20161018N002 Reservoir operation start years: uncomment unused variable
            while (biglakefile && (n < ng_glolakcells)) {
                biglakefile >> big_lakes_cells[0][n];
                biglakefile >> big_lakes_cells[1][n];
                biglakefile >> big_lakes_area_frac[n];
                n++;
            }

            biglakefile.close();
            cout << "Allocate GloLakeStorage for " << ng_glolakcells << " lakecells in wg_big_lakes_area_frac.txt."
                 << endl;
        } else {
            cout << "ERROR: wg_big_lakes_area_frac.txt not read." << endl;
            cout << "Using default values (no allocation of lakestorage on all lakecells)." << endl;
        }
    }
    // GFZ:GLOLAK end

    setParameters();

    // read parameters from ROUTING.DAT
    double v = defaultRiverVelocity;
    short s = defaultTimeStepsPerDay;

    try {
        sprintf(filename, "ROUTING.DAT");
        ifstream routingOptionFile(filename);

        if (!routingOptionFile) {
            throw routingOptionFile.rdstate();
        } else {
            char string[250];

            // read commentlines indicated by #
            while (!routingOptionFile.eof() && routingOptionFile.peek() == '#') {
                routingOptionFile.getline(string, sizeof(string));
            }

            routingOptionFile >> v; // river velocity
            if (!routingOptionFile.good())
                throw routingOptionFile.rdstate();
            else
                defaultRiverVelocity = v;

            routingOptionFile >> s;
            if (!routingOptionFile.good())
                throw routingOptionFile.rdstate();
            else
                defaultTimeStepsPerDay = s;

        }
    } catch (ios::iostate &status) {
        cout << filename << " not read. Using default values.\n";
        // cout << " Status: " << status << endl;
    }

    // defaults for velocity and timesteps
    riverVelocity = defaultRiverVelocity; // [km/day]
    timeStepsPerDay = defaultTimeStepsPerDay;

    cout << "Routing options: riverVelocity = " << riverVelocity << " km/day; "
         << timeStepsPerDay << " time steps per day." << endl;

    // read lakes, wetlands and reservoirs
    // and remove nodata values

    // fractions of grid cell area
    sprintf(filename, "%s/G_GLOLAK.UNF0", options.input_dir); // FP20161018N002 Global lakes WITHOUT regulated lakes
    gridIO.readUnfFile(filename, ng, &G_glo_lake[0]);
    sprintf(filename, "%s/G_LOCLAK.UNF0", options.input_dir);
    gridIO.readUnfFile(filename, ng, &G_loc_lake[0]);
    sprintf(filename, "%s/G_GLOWET.UNF0", options.input_dir);
    gridIO.readUnfFile(filename, ng, &G_glo_wetland[0]);
    sprintf(filename, "%s/G_LOCWET.UNF0", options.input_dir);
    gridIO.readUnfFile(filename, ng, &G_loc_wetland[0]);
    // FP20161018N002 Reservoir operation start years
    // read all areas in km2 into arrays of data type double
    sprintf(filename, "%s/G_LAKAREA.UNF0", options.input_dir);
    gridIO.readUnfFile(filename, ng, &G_lake_area[0]);
    sprintf(filename, "%s/G_RESAREA.UNF0", options.input_dir);
    gridIO.readUnfFile(filename, ng,
                       &G_reservoir_area_full[0]);  // reservoir AND regulated lake areas when all reservoirs/reg.lakes are in operation
    // fraction of grid cell area
    sprintf(filename, "%s/G_REGLAKE.UNF0",
            options.input_dir);  // regulated lake fraction when all reg.lakes are in operation
    gridIO.readUnfFile(filename, ng, &G_reg_lake[0]);
    //read regulated lake status (1 == lake is regulated lake), area included in G_RESAREA.UNF0
    sprintf(filename, "%s/G_REG_LAKE.UNF1", options.input_dir);
    gridIO.readUnfFile(filename, ng, &G_reg_lake_status[0]); //char

    // and remove nodata values
    for (int n = 0; n < ng; n++) {

        if (G_glo_lake[n] < 0.) G_glo_lake[n] = 0.;
        if (G_glo_wetland[n] < 0.) G_glo_wetland[n] = 0.;
        if (G_loc_lake[n] < 0.) G_loc_lake[n] = 0.;
        if (G_loc_wetland[n] < 0.) G_loc_wetland[n] = 0.;
        if (G_lake_area[n] < 0.) G_lake_area[n] = 0.;
        // FP20161018N002 Reservoir operation start years
        if (G_reservoir_area_full[n] < 0.) G_reservoir_area_full[n] = 0.;
        //G_reservoir_area[n] = 0;
        if (G_reg_lake[n] < 0.) G_reg_lake[n] = 0.;

        // initial values for fractions
        // (1) water areas (global lakes + local lakes + reservoirs (including regulated lakes)
        // (2) land areas (including wetlands)
        G_fwaterfreq[n] = 0.;
        G_landfreq[n] = 0.;

        // HMS 20170202 this should be done for each run during a calibration and not only once.
        // therefore it is moved to routing.initFractionStatus.
        /*// initialize G_fswb and related fractions
            G_fLocLake[n] = G_loc_lake[n] / 100.;
            G_fLocWet[n] = G_loc_wetland[n] / 100.;
            G_fGloLake[n] = G_glo_lake[n] / 100.; // do not adapt at all, just for output purposes
            //G_fGloLakeOutflowCell[n] =  (double)G_lake_area[n] / (geo.areaOfCellByArrayPos(n)); // can reach values > 1 due to relocalization into outflow cell.
            //G_fGloResOutflowCell[n] = (double)G_reservoir_area[n] / (geo.areaOfCellByArrayPos(n)); // can reach values > 1 due to relocalization into outflow cell.
            G_fGloWet[n] = G_glo_wetland[n] / 100.;

            //G_fswbInit[n] = G_fLocLake[n] + G_fLocWet[n] + G_fGloLakeOutflowCell[n] + G_fGloResOutflowCell[n] + G_fGloWet[n]; // for output purposes
            G_fswbInit[n] = G_fLocLake[n] + G_fLocWet[n] + G_fGloWet[n]; // for output purposes and for specifying catchment area (for fractional routing), outflow cells are not considered anymore
            G_fswbLandAreaFrac[n] = G_fswbInit[n]; // for changing LandAreaFrac only
            G_fswbLandAreaFracNextTimestep[n] = G_fswbLandAreaFrac[n];  // initialize for first time step
*/
        //define status at first occurrence / start = 0:
        // statusStarted_landAreaFracNextTimestep[n] = 0; //obsolete as covered in initFractionStatus
        // statusStarted_landfreq_fwaterfreq_landAreaFrac[n] = 0; //obsolete as covered in initFractionStatus

    }
    // end loop over grid cells

    // GLWDunit pattern with index to outflow grid cell Image22-Number (n = Image22nr - 1)
    sprintf(filename, "%s/G_OUTFLOW_CELL_ASSIGNMENT.UNF4", options.input_dir);
    gridIO.readUnfFile(filename, ng, &G_outflow_cell_assignment[0]); //int or short


    // read cell number of downstream cell;
    // added for new use_allocation algorithm, M.Hunger 2/2006,
    // and also for the routing algorithm (F.Voss 4/2008)
    sprintf(filename, "%s/G_OUTFLC.UNF4", options.routing_dir);
    gridIO.readUnfFile(filename, ng, G_downstreamCell);

    // new reservoir algorithm start

    // FP20161018N002 Reservoir operation start years

    //define status at model start / first occurence, e.g. when reading for annual reservoir fraction arrays
    // statusStarted_updateGloResPrevYear = 0; //obsolete as covered in initFractionStatus




    // Status of algorithm depends on
    // whether naturalized runs are performed or whether reservoirs (incl. regulated lakes) are to be used
    if (options.antNatOpt == 0) { // anthropogenic runs: consider reservoirs and regulated lakes as reservoirs
        cout
                       << "routing.init: if options.antNatOpt == 0, anthropogenic run: consider reservoirs and regulated lakes as reservoirs"
                       << endl;

        //select reservoir algorithm yes/no

        // Adding to lake area depending on
        // (1) reservoir algorithm option
        // (2) status of the reservoir / regulated lake: type, mean outflow
        // Calculate here, otherwise, areas would be added each year.

        // Anthropogenic run
        //read files which do not change in time
        // Use reservoir algorithm for reservoirs/reg.lakes
        if (options.resOpt == 1) {
            cout << "routing.init: if options.resOpt == 1, initialze new reservoir algorithm" << endl;

            // Dynamic arrays are deleted in destructor class ~routingClass()
            G_gloResStorage = new double[ng];
            G_res_type = new char[ng];
            G_start_month = new char[ng];
            G_mean_outflow = new double[ng];
            // FP20161018N002 Reservoir operation start years
            // G_res_start_year & G_stor_cap_full (instead G_stor_cap)
            G_res_start_year = new int[ng];
            G_stor_cap_full = new double[ng];

            G_mean_demand = new double[ng];
            G_mean_cons_use = new double[ng];
            G_alloc_coeff = new double[ng][reservoir_dsc];
            K_release = new double[ng];
            annual_release = new double[ng];

            //read reservoir_type
            sprintf(filename, "%s/G_RES_TYPE.UNF1", options.input_dir);
            gridIO.readUnfFile(filename, ng, &G_res_type[0]); //char

            // FP20161018N002 Reservoir operation start years
            //read reservoir operation start year
            sprintf(filename, "%s/G_START_YEAR.UNF4", options.input_dir);
            gridIO.readUnfFile(filename, ng, &G_res_start_year[0]); //int

            //read file -> start month of operational year
            sprintf(filename, "%s/G_START_MONTH.UNF1", options.routing_dir);
            gridIO.readUnfFile(filename, ng, &G_start_month[0]);//char

            //read mean annual outflow of reservoirs (long term mean [km3/month], calculated from G_RIVER_AVAIL), was G_MEAN_INFLOW before but is in fact inflow to reservoir plus P-PET of reservoir
            sprintf(filename, "%s/G_MEAN_OUTFLOW.UNF0", options.input_dir);
            gridIO.readUnfFile(filename, ng, &G_mean_outflow[0]);

            //read file -> reservoir volume from published data (GLWD1), storage capacity [km3]
            sprintf(filename, "%s/G_STORAGE_CAPACITY.UNF0", options.input_dir);
            gridIO.readUnfFile(filename, ng, &G_stor_cap_full[0]); // FP20161018N002 Reservoir operation start years

            //read file -> mean annual net use of surface water per cell (e.g. mean value of 1971 to 2000)
            // net use is given in m3/yr
            // FP20161018N002 Reservoir operation start years:
            // Specify name with averaging period of full year start and end in options
            sprintf(filename, "%s/G_NUs_%d_%d.UNF0", options.input_dir, options.resNUsMeanYearFirst,
                    options.resNUsMeanYearLast); // e.g. G_NUs_7100.UNF0
            gridIO.readUnfFile(filename, ng, &G_mean_NUs[0]);

            //read file -> allocation coefficient
            sprintf(filename, "%s/G_ALLOC_COEFF.%d.UNF0", options.routing_dir, reservoir_dsc); // is 5
            gridIO.readUnfFile(filename, ng * reservoir_dsc, &G_alloc_coeff[0][0]);

                // FP20161018N002 Reservoir operation start years:
                // Remark: Calculate this only once (1x)
                // Change in code: G_reservoir_area_full instead G_reservoir_area
                // Use all reservoirs and regulated lakes
                //calculate MEAN demand of downstream area
                for (int n = 0; n < ng; n++) {
                    G_mean_demand[n] = 0; //init
                    G_mean_demand[n] += G_mean_NUs[n]; //water demand of cell itself
                    G_glores_prevyear[n] = 0.;
                    G_glores_change[n] = 0.;
                    short i = 0;
                    int downstreamCell = G_downstreamCell[n];

                // max reservoir_dsc cells downstream or to the basin outlet or area up to next reservoir
                while (i < reservoir_dsc && downstreamCell > 0
                       && G_reservoir_area_full[downstreamCell] <= 0) {
                    //add demand of downstream cells, m3/yr
                    G_mean_demand[n] += G_mean_NUs[downstreamCell - 1] * G_alloc_coeff[n][i++];

                    // next downstream cell
                    downstreamCell = G_downstreamCell[downstreamCell - 1];
                }
                // end while loop
            }
            // end loop over grid cells

            // FP20161018N002 Reservoir operation start years
            // Remark: Calculate this only once (1x)
            // Change in code: G_stor_cap_full instead G_stor_cap
            // Calculate G_mean_outflow[n] & G_mean_demand[n]
            for (int n = 0; n < ng; n++) {
                if (G_stor_cap_full[n] < 0.)
                    G_stor_cap_full[n] = 0.;
                K_release[n] = 0.1; //HMS 2015-03-02 Experience from ISIMIP, leads to a reduction of initial outflow before start month is reached

                // FP20161018N002 Reservoir operation start years
                // No change in code, only 2 remarks:
                // (1) Calculate this only once (1x)
                // (2) Otherwise values are increasing until nan!
                // G_mean_outflow is long term mean annual value ([km3/month]); required number is mean annual value
                G_mean_outflow[n] = G_mean_outflow[n] * 12. * 1000000000. /
                                    31536000.; // 31536000=(365 * 24 * 60 * 60); //km3/month -> km3/yr -> m3/s
                //G_mean_demand[n] = G_mean_demand[n] * 1000000000. / 31536000.; // 31536000=(365 * 24 * 60 * 60); //km3/yr -> m3/s
                G_mean_demand[n] = G_mean_demand[n] /
                                   31536000.; // 31536000=(365 * 24 * 60 * 60); //km3/yr -> m3/s // G_mean_demand is in m�/yr as G_NUs_6190 is in km�/yr
            }
            // end loop over grid cells

        }
            // endif (options.resOpt == 1)
            // FP20161018N002 Reservoir operation start years
            // Anthropogenic run
            // Reservoir algorithm is not used, reservoirs and regulated lakes are handled as global lakes
        else if (options.resOpt == 0) { // reservoirs = lakes
            cout
                           << "routing.init: else if options.resOpt == 0, reservoirs and regulated lakes are handled as global lakes"
                           << endl;

            // If reservoir algorithm option is switched off, read in reservoir fraction (in percent) from reference year,
            // e.g. when all reservoirs are in operation
            // Preprocessing needed: Copy e.g. the last year G_RES_2014.UNF0 (ISIMIP21a: G_RES_2012.UNF0)
            sprintf(filename, "%s/G_RES/G_RES_FRAC.UNF0", options.input_dir);
            gridIO.readUnfFile(filename, ng, &G_glo_res[0]);

            for (int n = 0; n < ng; n++) {
                // treat area
                // for regulated lakes and reservoirs (in G_reservoir_area_full[n]), add area to lake area
                if (G_reservoir_area_full[n] > 0.) {
                    G_lake_area[n] += G_reservoir_area_full[n]; //add reservoir area to lake area
                    G_reservoir_area_full[n] = 0.; //set area to zero (inhibits further accounting elsewhere)
                }
                // treat land cover fractions
                //G_glo_lake[n] += G_reg_lake[n]; // add reg_lake fraction to glo_lake fraction // not needed, as G_glo_res includes G_reg_lake
                G_reg_lake[n] = 0.;  //set fraction to zero (inhibits further accounting elsewhere)
                G_glo_lake[n] += G_glo_res[n]; // add reservoir & regulated lake fraction to glo_lake fraction
                G_glo_res[n] = 0.;  //set fraction to zero (inhibits further accounting elsewhere)
                G_glores_prevyear[n] = 0.;  //set fraction to zero (inhibits further accounting elsewhere)
            }
            // end loop over grid cells
        }
        // end if (options.resOpt == 0)
    } // end if (options.antNatOpt == 0)
    if ((options.antNatOpt == 1) && (options.resOpt ==1)){
        cout << "Since options.antNatOpt ==1, options.resOpt cannot be 1: using options.resOpt ==1" << endl;
        options.resOpt = 0;

    }
    if ((options.antNatOpt == 1) && (options.resOpt == 0)) { //naturalized run
        cout << "routing.init: if options.antNatOpt == 1, naturalized run: using regulated lakes as lakes" << endl;

        for (int n = 0; n < ng; n++) {
            // treat land cover fractions
            G_glo_res[n] = 0.; //no reservoir fractions are used at all  // set fraction to zero (inhibits further accounting elsewhere)
            G_glores_prevyear[n] = G_glo_res[n]; //initialize with percent fraction of reservoirs of start or reference year
            G_glo_lake[n] += G_reg_lake[n]; //global lake fraction is increased by regulated lake fraction only (NOT reservoirs)
            G_reg_lake[n] = 0.; // set fraction to zero (inhibits further accounting elsewhere)
            // treat area
            if (G_reg_lake_status[n] == 1) {
                G_lake_area[n] += G_reservoir_area_full[n]; //add reservoir area to lake area in case of regulated lakes
                G_reservoir_area_full[n] = 0.; //set area to zero (inhibits further accounting elsewhere)  // new FP vs. ISIMIP of HMS
            }
        }
        // end loop over grid cells

    }
    // end else if (options.antNatOpt == 1)

    // FP20161018N002 Reservoir operation start years
    // Counting number of reservoirs and regulated lakes
    long count = 0;
    for (int n = 0; n < ng; n++) {
        if (G_reservoir_area_full[n] > 0.) count++;
    }
    cout << "routing.init: nr. of reservoirs/regulated lakes detected (area > 0)  : " << count << endl;

    //new reservoir algorithm end

    // set correction factor for upstream locations to 1 (= no correction) use statcorrOpt == 1 during calibration runs (here, station correction should occur). During normal runs, this line should be normally 0 (runs without station correction).
    if (1 == options.statcorrOpt) {
        for (int n = 0; n < ng; n++)
            G_statCorrFact[n] = 1.;
    }

    //before calculating landareaFrac, initial riverFrac needs to be calculated, based on bankfull flow and default flow velocity.

    // read stretch length for each river segment
    // (added for new routing-algorithm, F.Voss 4/2008)
    sprintf(filename, "%s/G_RIVER_LENGTH.UNF0", options.routing_dir);
    gridIO.readUnfFile(filename, ng, G_riverLength);

#ifdef _WATERGAP_CHECKS_GLOBAL_H
                                                                                                                            // Check Felix Portmann 2015 - for invalid data
            double limit_zero_RiverLength = 0.00000000000000001;
            long n_cells_zero__file = 0;
            long n_cells_zero__transformed = 0;
            printf("Check for G_riverLength[n] < %e , i.e. zero\n", limit_zero_RiverLength);
            cout << "(1) as read from from file (if any, possibly none)" << filename << endl;
            for (int n = 0; n < ng; n++) {
                if (G_riverLength[n] < limit_zero_RiverLength) {
                    printf("G_riverLength[n] zero %e, geo.G_contcell[n] %d, at grid cell n: %d, geo.G_landfreq[n]: %e, geo.G_fwaterfreq[n]: %e \n", G_riverLength[n], geo.G_contcell[n], n, geo.G_landfreq[n], geo.G_fwaterfreq[n]);
                    n_cells_zero__file += 1;
                }
            }
            printf("(Data read from file) Number of cells with zero percentage of continental area: %d, limit %e\n", n_cells_zero__file, limit_zero_RiverLength);
#endif

    // Introduced in WG22b: Reduce river length in coastal cells in proportion to continental area fraction:
    for (int n = 0; n < ng; n++) {
        G_riverLength[n] *= geo.G_contfreq[n] / 100.;
    }

#ifdef _WATERGAP_CHECKS_GLOBAL_H
                                                                                                                            // Check Felix Portmann 2015 - for invalid data
            cout << "(2) WG22b reduced G_riverLength[n] 'in coastal cells'" << endl;
            for (int n = 0; n < ng; n++) {
                if (G_riverLength[n] < limit_zero_RiverLength) {
                    printf("G_riverLength[n] zero %e, geo.G_contcell[n] %d, at grid cell n: %d, geo.G_landfreq[n]: %e, geo.G_fwaterfreq[n]: %e \n", G_riverLength[n], geo.G_contcell[n], n, geo.G_landfreq[n], geo.G_fwaterfreq[n]);
                    n_cells_zero__transformed += 1;
                }
            }
            printf("(Transformed data) Number of cells with zero percentage of continental area: %d, limit %e\n", n_cells_zero__transformed, limit_zero_RiverLength);
#endif

    // velocity -------------------------------------------------------
    //if (options.riverveloOpt == 1) {
    sprintf(filename, "%s/G_RIVERSLOPE.UNF0", options.routing_dir);
    gridIO.readUnfFile(filename, ng, G_RiverSlope);

    sprintf(filename, "%s/G_ROUGHNESS.UNF0", options.input_dir);
    gridIO.readUnfFile(filename, ng, G_Roughness);

    sprintf(filename, "%s/G_BANKFULL.UNF0", options.input_dir);
    gridIO.readUnfFile(filename, ng, G_bankfull_flow); //[m/sec]

    for (int n = 0; n < ng; n++) {
        // to avoid negative bottom width
        if (G_bankfull_flow[n] < 0.05)
            G_bankfull_flow[n] = 0.05;

        G_RiverWidth_bf[n] = 2.71 * pow(G_bankfull_flow[n], 0.557); //[m]
        G_RiverDepth_bf[n] = 0.349 * pow(G_bankfull_flow[n], 0.341); //[m]
        // calculate river bottom width asssuming a trapezoidal channes with 2/1 run to rise ratio; [m]
        G_riverBottomWidth[n] = G_RiverWidth_bf[n] - 2.0 * 2.0 * G_RiverDepth_bf[n];
        //CR: G_riverStorageMax[n] is required to compute riverEvapoReductionFactor introduced in WG22b.
        G_riverStorageMax[n] = G_riverLength[n] * 0.5 * G_RiverDepth_bf[n] / 1000. *
                               (G_riverBottomWidth[n] / 1000. + G_RiverWidth_bf[n] / 1000.); // [km3]
        //set initial value for G_riverAreaFrac to 0 (later updated)
        G_riverAreaFrac[n] = 0.;

    }
    //}
    // End velocity

    // FP20161018N002 Reservoir operation start years:
    // G_landAreaFrac[n] and G_landWaterExclGloLakAreaFrac[n] now calculated in routing.annualInit

    // set initialization of G_landAreaFracNextTimestep to -99 for the first model day
    // to calculate new special value at the end of the first model day
    for (int n = 0; n < ng; n++) {
        G_landAreaFracNextTimestep[n] = -99.0;
        G_landAreaFracNextTimestep[n] = 0.; // Set to zero to prevent undesired results
        G_landAreaFrac[n] = 0.; // Set to zero to prevent undesired results
        G_landAreaFracPrevTimestep[n] = 0.; // Set to zero to prevent undesired results
    }

    // Superbasin calculation currently obsolete, probably not working correctly in WaterGAP2.2b, do in postprocessing  // FP 2015-10
    // calculate land area of basins
//		extern signed short G_sbasin[ng];  // STRANGE: SEEMS NOT TO BE RECOGNIZED by Qtcreator as to be from watergap.cpp)
//	basinLandArea = new double[nBasins];

//	for (int n = 0; n < nBasins; n++)
//		basinLandArea[n] = 0.;
//	for (int n = 0; n < ng; n++) {
//		if (G_sbasin[n] > 0) {
//			basinLandArea[G_sbasin[n] - 1] += G_landAreaFrac[n] * geo.areaOfCellByArrayPos(n) / 100.0;
//		}
//	}

    // set G_GloLakePrecip and G_ReservoirPrecip to 0

    for (int n = 0; n < ng; n++) {
        G_GloLakePrecip[n] = 0.;
        G_ReservoirPrecip[n] = 0.;
        G_dailyRiverPrecip[n] = 0.;
    }

    // read image numbers of neighbouring cells of each cell
    // (added for new use_allocation algorithm, M.Hunger 2/2006)
    if (options.use_alloc > 0) {
        sprintf(filename, "%s/G_NEIGHBOUR_CELLS.8.UNF4", options.routing_dir);
        gridIO.readUnfFile(filename, ng * 8, &G_neighbourCell[0][0]);
    }

    // read ordered routing number of each cell
    // and assign array position for ordered routing scheme
    // (added for new routing-algorithm, F.Voss 4/2008)
    int *G_routOrder;
    G_routOrder = new int[ng];
    sprintf(filename, "%s/G_ROUT_ORDER.UNF4", options.routing_dir);
    gridIO.readUnfFile(filename, ng, G_routOrder);
#pragma omp parallel for \
        shared(G_routingCell, G_routOrder)
    for (int i = 0; i < ng; ++i)
        G_routingCell[G_routOrder[i] - 1] = i + 1;
    delete[] G_routOrder;
    G_routOrder = NULL;

    // read frgi to adapt NAg if necessary
    if (options.subtract_use > 0) {
        sprintf(filename, "%s/G_FRACTRETURNGW_IRRIG.UNF0", options.input_dir);
        gridIO.readUnfFile(filename, ng, G_fractreturngw_irrig);
    }

//	// read stretch length for each river segment
//	// (added for new routing-algorithm, F.Voss 4/2008)
//	sprintf(filename, "%s/G_RIVER_LENGTH.UNF0", options.routing_dir);
//	gridIO.readUnfFile(filename, ng, G_riverLength);

//	// velocity -------------------------------------------------------
//        //if (options.riverveloOpt == 1) {
//		sprintf(filename, "%s/G_RIVERSLOPE.UNF0", options.routing_dir);
//		gridIO.readUnfFile(filename, ng, G_RiverSlope);

//		sprintf(filename, "%s/G_ROUGHNESS.UNF0", options.input_dir);
//		gridIO.readUnfFile(filename, ng, G_Roughness);

//		sprintf(filename, "%s/G_BANKFULL.UNF0", options.input_dir);
//		gridIO.readUnfFile(filename, ng, G_bankfull_flow); //[m/sec]

//		for (int n = 0; n < ng; n++) {
//			// to avoid negative bottom width
//			if (G_bankfull_flow[n] < 0.05)
//				G_bankfull_flow[n] = 0.05;

//			G_RiverWidth_bf[n] = 2.71 * pow(G_bankfull_flow[n], 0.557); //[m]
//			G_RiverDepth_bf[n] = 0.349 * pow(G_bankfull_flow[n], 0.341); //[m]
//			// calculate river bottom width asssuming a trapezoidal channes with 2/1 run to rise ratio; [m]
//			G_riverBottomWidth[n] = G_RiverWidth_bf[n] - 2.0 * 2.0 * G_RiverDepth_bf[n];
//			//CR: G_riverStorageMax[n] is required to compute riverEvapoReductionFactor introduced in WG22b.
//			G_riverStorageMax[n] = G_riverLength[n] * 0.5 * G_RiverDepth_bf[n]/1000. * (G_riverBottomWidth[n]/1000. + G_RiverWidth_bf[n]/1000.); // [km3]
//		}
//        //}
//	// Ende velocity

    // local drain direction (= flow direction map)
//	sprintf(filename, "%s/G_LDD_2.UNF1", options.routing_dir);
//	gridIO.readUnfFile(filename, ng, G_LDD);                        // CR 2015-05-12: required in WG22 (and previous versions) to set 'transportedVolume' to zero

    // number of superbasins varies.
    // therefore size of the arrays is allocated dynamically.
    nSpecBasins = nBasins; // number of specified basins

    cout << " : ##### routing.init:  nSpecBasins= " << nSpecBasins << endl;

    dailyRiverDischarge = new double[nSpecBasins][365];
    // initialize daily storages
    if (options.day_store == 1) {
        dailyLocLakeStorage = new double[nSpecBasins][365];
        dailyLocWetlStorage = new double[nSpecBasins][365];
        dailyGloLakeStorage = new double[nSpecBasins][365];
        dailyGloWetlStorage = new double[nSpecBasins][365];
        dailyRiverVelocity = new double[nSpecBasins][365];   //-------------new for daily velocity file------------------
        if (options.resOpt == 1)
            dailyResStorage = new double[nSpecBasins][365];
    }
    // endif options.day_store

    // initialize monthly storages
    if ((2 == options.grid_store) || (4 == options.grid_store) || (6 == options.grid_store)) {
        G_monthlyConsistentPrecip = new double[ng][12];
        G_monthlyCellRunoff = new double[ng][12];
        G_monthlyCellSurfaceRunoff = new double[ng][12];
        G_monthlyRiverAvail = new double[ng][12];
        G_monthlyPotCellRunoff = new double[ng][12];
        G_monthlyPotCellRunoffDeficit = new double[ng][12];
        if (options.outRiverInUpstream) G_monthlyRiverInUpstream = new double[ng][12]; // total inflow from upstream

        // added for cell AET calculation (2.1f)
        G_monthlyCellAET = new double[ng][12]; // total cell evap
        G_monthlyCellAETWCa = new double[ng][12]; // total cell evap plus consumptive water use
        G_monthlyOpenWaterEvap = new double[ng][12]; // evap from open water bodies
        G_monthlyVelocity = new double[ng][12]; // velocity grid
        G_monthlySurfStor = new double[ng][12]; // surface water storage
        G_monthlySurfStor_mm = new double[ng][12]; // surface water storage [mm]

        G_monthlyMinRiverAvail = new double[ng][12];
        G_monthlyMaxRiverAvail = new double[ng][12];

        G_monthlyGwrSwb = new double[ng][12];
        G_monthlyFswb = new double[ng][12];
        G_monthlyLandAreaFrac = new double[ng][12];
        G_monthlyRiverAreaFrac = new double[ng][12];
        G_monthlyRiverPET = new double[ng][12];
        G_monthlyGwRunoff = new double[ng][12];
        G_monthlyLocWetlExtent = new double[ng][12];
        G_monthlyGloWetlExtent = new double[ng][12];

        //define arrays for storage output
        if (options.resOpt == 1) {
            G_monthlyResStorage = new double[ng][12];
            G_monthlyResStorageRatio = new double[ng][12];
            G_monthlyResInflow = new double[ng][12];
            G_monthlyResOutflow = new double[ng][12];
        }
        G_monthlyRiverStorage = new double[ng][12];
        G_monthlyLocLakeStorage = new double[ng][12];
        G_monthlyLocWetlStorage = new double[ng][12];
        G_monthlyGloLakeStorage = new double[ng][12];
        G_monthlyGloWetlStorage = new double[ng][12];
        G_monthlySatisfiedUse = new double[ng][12];
        G_monthlyGwStorage = new double[ng][12];
        // G_monthlySwStorage		= new double[ng][12]; // not used at the moment, implement if needed
        G_monthlyActualUse = new double[ng][12];
        G_monthlyNUg = new double[ng][12];
        // CR 15-08-19 used for test output
        G_monthlyAllocatedUse = new double[ng][12]; // FP 2015, not used anymore, at the moment, implement if needed again
        //KF *** end
    }
    // endif options (initialize monthly storages)

    // daily output option 31 start
    if ((3 == options.grid_store) || (4 == options.grid_store)) {
        if (options.outPrecDaily) G_daily31ConsistentPrecip = new double[ng][31];
        if (options.outCellRunoffDaily) G_daily31CellRunoff = new double[ng][31];
        if (options.outGWRunoffDaily) G_daily31GwRunoff = new double[ng][31];
        if (options.outCellSurfaceDaily) G_daily31CellSurfaceRunoff = new double[ng][31];
        if (options.outRiverAvailDaily) G_daily31RiverAvail = new double[ng][31];
        if (options.outRiverVeloDaily) G_daily31Velocity = new double[ng][31];
        if (options.outCellAETDaily) G_daily31CellAET = new double[ng][31];
        if (options.outSurfStorDaily) G_daily31SurfStor = new double[ng][31];
        if (options.outSingleStoragesDaily) G_daily31LocLakeStor = new double[ng][31];
        if (options.outSingleStoragesDaily) G_daily31LocWetlStor = new double[ng][31];
        if (options.outSingleStoragesDaily) G_daily31GloLakeStor = new double[ng][31];
        if (options.outSingleStoragesDaily) G_daily31GloWetlStor = new double[ng][31];
        if (options.outSingleStoragesDaily) G_daily31RiverStor = new double[ng][31];
        if (options.outSingleStoragesDaily) G_daily31ResStor = new double[ng][31];
        if (options.outSingleStoragesDaily || options.outGWStorageDaily) G_daily31GwStor = new double[ng][31];
        if (options.outTotalWaterInStoragesDaily_km3 ||
            options.outTotalWaterInStoragesDaily_mm)
            G_daily31TotalWaterInStorages_km3 = new double[ng][31];
        if (options.outTotalWaterInStoragesDaily_mm) G_daily31TotalWaterInStorages_mm = new double[ng][31];
        if (options.outGwrSwbDaily) G_daily31GwrSwb = new double[ng][31];
        if (options.outFswbDaily) G_daily31Fswb = new double[ng][31];
        if (options.outLandAreaFracDaily) G_daily31LandAreaFrac = new double[ng][31];
        if (options.outGwrunSurfrunDaily) G_daily31GwrunSurfrun = new double[ng][31];
        if (options.outCellAETWCaDaily) G_daily31CellAETWCa = new double[ng][31];

    }
    // endif options (initialize output daily 31 storages)
    // daily output option 31 end
    //HMS first and last tws value per year
    if (options.outTotalWaterInStoragesStartEndDaily_km3) G_startendTotalWaterInStorages_km3 = new double[ng][2];
    // daily output option 365 start
    if ((5 == options.grid_store) || (6 == options.grid_store)) {

//CR #############################################################################################
//	G_daily365ActualUse				= new double[ng][365];
//	G_daily365AllocUse				= new double[ng][365];
//	G_daily365AllocUseNextDay		= new double[ng][365];
//	G_daily365TotalUnsatUse			= new double[ng][365];
//	G_daily365UnsatAllocUseNextDay	= new double[ng][365];
//CR #############################################################################################

        if ((options.outPrecDaily) || ((options.scoutcPrecip) && (2 == options.day_store)))
            G_daily365ConsistentPrecip = new double[ng][365];

        if ((options.outCellRunoffDaily) ||
            (((options.scoutCellRunoffkm3) || (options.scoutCellRunoffmm)) && (2 == options.day_store)))
            G_daily365CellRunoff = new double[ng][365];
        if ((options.outCellRunoffDaily) ||
            ((options.scoutCellRunoffmm) && (2 == options.day_store)))
            G_daily365CellRunoff_mm = new double[ng][365];
        if ((options.outCellSurfaceDaily) ||
            ((options.scoutCellSRunoff) && (2 == options.day_store)))
            G_daily365CellSurfaceRunoff = new double[ng][365]; // FP20160915N002
        if ((options.outRiverAvailDaily) || ((options.scoutQ) && (2 == options.day_store)))
            G_daily365RiverAvail = new double[ng][365];
        if ((options.outRiverVeloDaily) || ((options.scoutFlowVelo) && (2 == options.day_store)))
            G_daily365Velocity = new double[ng][365];
        if ((options.outGWRunoffDaily) || ((options.scoutGwRunoff) && (2 == options.day_store)))
            G_daily365GwRunoff = new double[ng][365];
        if ((options.outLandAreaFracDaily) ||
            ((options.scoutLandAreaFrac) && (2 == options.day_store)))
            G_daily365LandAreaFrac = new double[ng][365];
        if ((options.outCellAETDaily) || ((options.scoutCellAET) && (2 == options.day_store)))
            G_daily365CellAET = new double[ng][365];
        if ((options.outSurfStorDaily) || ((options.scoutSurfaceStor) && (2 == options.day_store)))
            G_daily365SurfStor = new double[ng][365];
        if ((options.outSingleStoragesDaily) ||
            ((options.scoutLocLake) && (2 == options.day_store)))
            G_daily365LocLakeStor = new double[ng][365];
        if ((options.outSingleStoragesDaily) ||
            ((options.scoutLocWet) && (2 == options.day_store)))
            G_daily365LocWetlStor = new double[ng][365];
        if ((options.outSingleStoragesDaily) ||
            ((options.scoutGloLake) && (2 == options.day_store)))
            G_daily365GloLakeStor = new double[ng][365];
        if ((options.outSingleStoragesDaily) ||
            ((options.scoutGloWet) && (2 == options.day_store)))
            G_daily365GloWetlStor = new double[ng][365];
        if ((options.outSingleStoragesDaily) || ((options.scoutRiver) && (2 == options.day_store)))
            G_daily365RiverStor = new double[ng][365];
        if ((options.outSingleStoragesDaily) || ((options.scoutReservoir) && (2 == options.day_store)))
            G_daily365ResStor = new double[ng][365];
        if ((options.outSingleStoragesDaily || options.outGWStorageDaily) ||
            ((options.scoutGwStor) && (2 == options.day_store)))
            G_daily365GwStor = new double[ng][365];
        if ((options.outTotalWaterInStoragesDaily_km3 || options.outTotalWaterInStoragesDaily_mm) ||
            (((options.scoutTWSkm3) || (options.scoutTWSmm)) && (2 == options.day_store)))
            G_daily365TotalWaterInStorages_km3 = new double[ng][365];
        if ((options.outTotalWaterInStoragesDaily_mm) ||
            ((options.scoutTWSmm) && (2 == options.day_store)))
            G_daily365TotalWaterInStorages_mm = new double[ng][365];
        if ((options.outGwrSwbDaily) || ((options.scoutGwrSwb) && (2 == options.day_store)))
            G_daily365GwrSwb = new double[ng][365];
        if ((options.outFswbDaily) || ((options.scoutFswb) && (2 == options.day_store)))
            G_daily365Fswb = new double[ng][365];
        if (options.outGwrunSurfrunDaily) G_daily365GwrunSurfrun = new double[ng][365];
        if (options.outCellAETWCaDaily) G_daily365CellAETWCa = new double[ng][365];
    }
    // endif options (initialize output daily 365 storages)
    // daily output option 365 end

    setStoragesToZero();

}
// end init

void routingClass::initFractionStatus() {
    // this subroutine is necessary to reset the fractions to its initial values
    // in each first year during a calibration loop
    // otherwise, fractions are used from the previous loop and calibration fails.

    for (int n = 0; n < ng; n++) {
        statusStarted_landfreq_fwaterfreq_landAreaFrac[n] = 0;
        statusStarted_landAreaFracNextTimestep[n] = 0;

        // initialize G_fswb and related fractions
        G_fLocLake[n] = G_loc_lake[n] / 100.;
        G_fLocWet[n] = G_loc_wetland[n] / 100.;
        G_fGloLake[n] = G_glo_lake[n] / 100.; // do not adapt at all, just for output purposes
        G_fGloWet[n] = G_glo_wetland[n] / 100.;
        G_fswbInit[n] = G_fLocLake[n] + G_fLocWet[n] +
                        G_fGloWet[n]; // for output purposes and for specifying catchment area (for fractional routing), outflow cells are not considered anymore
        G_fswbLandAreaFrac[n] = G_fswbInit[n]; // for changing LandAreaFrac only
        G_fswbLandAreaFracNextTimestep[n] = G_fswbLandAreaFrac[n];  // initialize for first time step

    }
    statusStarted_updateGloResPrevYear = 0;
}

void routingClass::setStoragesToZero() {
    // this subroutine concerns those storages
    // which have to be set to zero at the beginning of
    // a normal simulation, or at the beginning of each
    // pass of a calibration or sensitivity loop
    //
    // additional storages have to be set to zero at the beginning
    // of each simulation year

    // set storage value of lakes, wetlands and rivers to zero
    for (int n = 0; n < ng; n++) {
        G_locLakeStorage[n] = 0.;
        G_locWetlStorage[n] = 0.;
        G_gloLakeStorage[n] = 0.;
        G_gloWetlStorage[n] = 0.;
        // FP20161018N002 Reservoir operation start years
        G_landStorageChange[n] = 0.;
        //initialize release coefficient for reservoirs
        if ((options.antNatOpt == 0) && (options.resOpt == 1)) {
            G_gloResStorage[n] = 0.; // for new reservoir algorithm
        }

        if ((options.riverveloOpt == 0) || (options.riverveloOpt == 1)) {
            G_riverStorage[n] = 0.;
        }
        else {
            G_riverStorage[n] = G_riverStorageMax[n];
        }

        G_groundwaterStorage[n] = 0.;
        // G_totalUnsatisfiedUse is copied to G_UnsatisfiedUsePrevYear
        // at the beginning of each year.
        // Therefore it has to be set to zero before doing this
        // for the first time.
        G_totalUnsatisfiedUse[n] = 0.;        // should be initialized in any case

        // Variables introduced for Option "aggrNUsGloLakResOpt":
        // The amount of unsatisfied use in riparian cells of global lakes or reservoirs
        // can be satisfied from river storage in the next time step.
        // This is done independently from Option "delayedUseSatisfaction".
        G_unsatUseRiparian[n] = 0.;
        G_remainingUse[n] = 0.;
        G_remainingUseRes[n] = 0.;
        G_remainingUseGloLake[n] = 0.;
        G_remainingUseGloLakeRedistr[n] = 0.;
        G_remainingUseResRedistr[n] = 0.;
        G_totalDesiredUseCell[n] = 0.;


        // New approach in 2.2b for Option "use_alloc":
        // If the remainingUse cannot be satisfied in a cell, the value is added to the remainingUse of a second cell.
        // Depending on the routing order, the water balance of the second cell is computed in the same or the next timestep.
        // This way, NAs is consistently included in the storage equations of lakes and the river.
        G_AllocatedUse[n] = 0.;
        G_daily_allocatedUseNextDay[n] = 0.;
        G_dailyAllocatedUse[n] = 0.;
        G_UnsatAllocUse[n] = 0.;
        G_daily_UnsatAllocUseNextDay[n] = 0.;
        G_AllocUseToNeigborCell[n] = 0.;
        G_SecondCell[n] = 0;    // int

    }
}
// end setStoragesToZero
// reading dailyNUs und dailyNUg once a year in m3/month

void routingClass::dailyNUInit(const char *input_dir, const short new_year) {
    char filename[250];
    short month;

    // if this is a calibration run, allocated surface water and adapted net abstraction from groundwater should be
    // used.
    if (3 == options.subtract_use) {
        sprintf(filename, "%s/G_ACTUAL_NAS_%d.12.UNF0", input_dir, new_year);
        gridIO.readUnfFile(filename, ng * 12, &G_dailyNUs[0][0]);

        sprintf(filename, "%s/G_ACTUAL_NAG_%d.12.UNF0", input_dir, new_year);
        gridIO.readUnfFile(filename, ng * 12, &G_dailyNUg[0][0]);
    } else { // normal runs
        sprintf(filename, "%s/G_NETUSE_SW_m3_%d.12.UNF0", input_dir, new_year);
        gridIO.readUnfFile(filename, ng * 12, &G_dailyNUs[0][0]);

        sprintf(filename, "%s/G_NETUSE_GW_m3_%d.12.UNF0", input_dir, new_year);
        gridIO.readUnfFile(filename, ng * 12, &G_dailyNUg[0][0]);
    }

    //cout << "read file " << input_dir << "/G_NETUSE_GW__m3_" << new_year << ".12.UNF0" << endl;

    sprintf(filename, "%s/G_IRRIG_WITHDRAWAL_USE_m3_%d.12.UNF0", input_dir, new_year);
    gridIO.readUnfFile(filename, ng * 12, &G_monthlyWUIrrig[0][0]);

    sprintf(filename, "%s/G_IRRIG_CONS_USE_m3_%d.12.UNF0", input_dir, new_year);
    gridIO.readUnfFile(filename, ng * 12, &G_monthlyCUIrrig[0][0]);

    // Cell-specific calibration parameters - Apply multiplier  // FP
    for (int n_local = 0; n_local < ng; n_local++) {
        for (month = 0; month < 12; month++) {
            G_dailyNUs[n_local][month] = calibParam.getValue(M_NETABSSW, n_local) * G_dailyNUs[n_local][month];
            G_dailyNUg[n_local][month] = calibParam.getValue(M_NETABSGW, n_local) * G_dailyNUg[n_local][month];
        }
    }
    // end loop over grid cells


    // read file -> read GLWD units of global lakes and reservoirs (Option 28)
    sprintf(filename, "%s/GLWDunits.UNF4", options.input_dir);
    gridIO.readUnfFile(filename, ng, &G_glwdUnit[0]);

    // CR 2015: Calculation of array G_dailyNUsAggregated[n][month] for outflow and riparian cells of global lakes/reservoirs and all other cells.
    // outflow cells: all positive values of dailyNUs are aggregated over riparian cells and are added to the value of the outflow cell
    // negative values in riparian cells: G_dailyNUsAggregated[n][month] = dailyNUs[n][month] (initial value)
    // positive values in riparian cells (no outflow cell): G_dailyNUsAggregated[n][month] = 0 (since value is added to outflow cell)
    // all other cells (neither riparian nor outflow cell): G_dailyNUsAggregated[n][month] = dailyNUs[n][month] (initial value)
    if (1 == options.aggrNUsGloLakResOpt) {

        //  cout << "options.aggrNUsGloLakResOpt == 1: Aggregation of NUs over riparian cells of global lakes/reservoirs" << endl;

        for (int n = 0; n < ng; n++) {
            for (month = 0; month < 12; month++) {
                G_dailyNUsAggregated[n][month] = G_dailyNUs[n][month];
            }
        }

        // Aggregation of NUs values (except negative values) over GLWD units
        // (GLWD units include all riparian cells of a global lake or reservoir)

        for (month = 0; month < 12; month++) {

            for (int n = 0; n < ng; n++) {

                if (G_lake_area[n] > 0. || G_reservoir_area[n] > 0.) {        // if cell is an outflow cell

                    for (int i = 0; i < ng; i++) {

                        if (n == i) {        // to avoid counting the NUs value of the outflow cell twice
                            continue;
                        }


                        if (G_glwdUnit[n] == G_glwdUnit[i]) {

                            if (G_dailyNUs[i][month] > 0.) {

                                G_dailyNUsAggregated[n][month] += G_dailyNUsAggregated[i][month];
                                G_dailyNUsAggregated[i][month] = 0.;

                            } //else (negative value): G_dailyNUsAggregated remains G_dailyNUs
                        }
                    }
                }
            }
        }

    } // end if (1 == options.aggrNUsGloLakResOpt)

// calculate daily values, if this is a calibration run, G_dailyNUs is also in m3
    for (month = 0; month <= 11; month++) {
        for (int n = 0; n <= ng - 1; n++) {
            G_dailyNUs[n][month] = G_dailyNUs[n][month] / (1000000000. * (double) numberOfDaysInMonth[month]);
            G_dailyNUg[n][month] = G_dailyNUg[n][month] / (1000000000. * (double) numberOfDaysInMonth[month]);
            if (1 == options.aggrNUsGloLakResOpt) {
                G_dailyNUsAggregated[n][month] =
                               G_dailyNUsAggregated[n][month] / (1000000000. * (double) numberOfDaysInMonth[month]);
            }
        }
    }
}
// end dailyNUInit

void routingClass::annualInit(const short year) {
    // FP20161018N002: Do NOT introduce nSpecBasins as local variable nBasins, as superbasin area calculation probably wrong & obsolete
    short m;  // index for months
    int n;  // index for grid cells
    char filename[250];

    // Part 1: Initialize always needed arrays
    for (n = 0; n < ng; n++) {

        // store the values of the grid, so that we can compare
        // at the end of the year if unsatisfied use of the previous year
        // has been satisfied with water availability from the current year
        G_UnsatisfiedUsePrevYear[n] = G_totalUnsatisfiedUse[n];
        // G_dailyRemainingUse[n] = 0.;	// CR 2015-05-11: required for old approach for Option "use_alloc" in WG22
        // G_satisfiedUse[n] = 0.;      // CR 2015-09: A distinction between satisfied use and actual use is not possible in 22b.
        G_potCellRunoff[n] = 0.;        // HMS 2014-06-04 reintroduced

        //		G_AnnualCellRunoff[n] 	= 0.; // HMS 2014-06-04 commented out
        //G_annualRiverAvail[n] = 0.;
        G_actualUse[n] = 0.; // used for new use_allocation (M.Hunger 2/2006)
        // FP20161018N002 Reservoir operation start years
        // Reservoir and regulated lakes area and storage capacity of current year
        G_reservoir_area[n] = 0.;
        G_stor_cap[n] = 0.;

    }
    // end loop over grid cells

    //
    // FP20161018N002 Reservoir operation start years
    //
    // Part 2: Fill/initialize arrays for the current year
    // e.g. G_reservoir_area, G_stor_cap
    //
    // Status of algorithm depends on
    // whether naturalized runs are performed or whether reservoirs (incl. regulated lakes) are to be used
    // Here, only anthrogenic runs are of interest, as then year-specific values have to be considered
    if (options.antNatOpt == 0) { // anthropogenic runs: consider reservoirs and regulated lakes as reservoirs
        //     cout << "routing.annualInit Part 2: if options.antNatOpt == 0, anthropogenic run: consider reservoirs and regulated lakes as reservoirs" << endl;

        // Anthropogenic run
        // Use reservoirs and regulated lakes
        if (options.resOpt == 1) {

            // Always for regulated lakes, reservoir_area and storage capacity needs to be considered
            for (n = 0; n < ng; n++) {
                // 1. Do not consider another time land cover fraction of regulated lakes
                //      because they are already in G_glo_res (that should be used in subsequent treatment, if concerned)
                //      But as G_reservoir_area is initialized with zero every year, the transfer has to be done every year
                //      from G_reservoir_area_full as read from the respective file (that always should contain regulated lakes, too)
                // 2. Treat area and storage capacity
                if (G_reg_lake_status[n] == 1) {
                    //reservoir area is needed anyhow from regulated lakes, otherwise they cannot be calculated
                    //in the reservoir algorithm
                    G_reservoir_area[n] = G_reservoir_area_full[n];
                    //Previous selection, currently discarted:
                    //but storage capacity is only needed when regulated lake gets operated
                    //not true anymore, as we let calculate regulated lakes with stor cap also if it is not yet regulated.
                    //if (((options.resYearOpt == 1) && (year >= G_res_start_year[n]))
                    //        ||((options.resYearOpt == 0) && (options.resYearReference >= G_res_start_year[n])))
                    G_stor_cap[n] = G_stor_cap_full[n];
                }
                // end if G_reg_lake_status[n] == 1
            }
            // end loop over grid cells

            // Define desired year of reservoir data (fraction, area, storage capacity)
            // Depending on reservoir start year option yes/no
            if (options.resYearOpt == 1) { // if reservoir year option is switched on

                // Limit years to valid domain specified by
                // options.resYearFirstToUse & options.resYearLastToUse
                // Maximum = last year
                if (year > options.resYearLastToUse) {
                    resYearCurrentToUse = options.resYearLastToUse;
                }
                    // Minimum = first year
                else if (year < options.resYearFirstToUse) {
                    resYearCurrentToUse = options.resYearFirstToUse;
                }
                    // Within domain: Current year
                else {
                    resYearCurrentToUse = year;
                }

            }
                // end if resYearOpt == 1
            else if (options.resYearOpt == 0) { //resYearOpt == 0

                // Using reference year
                resYearCurrentToUse = options.resYearReference;  // ISIMIP21a: G_RES_2000.UNF0

            }
            // end else if resYearOpt == 0

            // Read desired reservoir and regulated lake fractions
            // within domain of years, e.g. current, min, max, reference year
            // ATTENTION: Fractions should be consistent to condition resYearCurrentToUse >= G_res_start_year[n]
            sprintf(filename, "%s/G_RES/G_RES_%d.UNF0", options.input_dir, resYearCurrentToUse);
            gridIO.readUnfFile(filename, ng, &G_glo_res[0]);


            // Transfer the values at start of the first model year to array of previous
            if (0 == statusStarted_updateGloResPrevYear) {
                routingClass::updateGloResPrevYear_pct();  // copying values of current year
                statusStarted_updateGloResPrevYear = 1; // indicate for the next iteration that first year has been treated / model has started
            }
            // end if

            // Treat area and storage capacity (outflow cells)
            for (n = 0; n < ng; n++) {
                // (1.1) update G_reservoir_area and G_stor_cap if reservoir is in operation, e.g. since reference year
                if (resYearCurrentToUse >= G_res_start_year[n]) {
                    G_reservoir_area[n] = G_reservoir_area_full[n];
                    G_stor_cap[n] = G_stor_cap_full[n];
                }
                // (1.2) Define current storage to be the maximum storage capacity of the selected year (current, max, min, reference)
                G_actualStorageCapacity[n] = G_stor_cap[n];
            }
            // end loop over grid cells


            // In case of erroneous attributes of a reservoir, treat it as a global lake:
            // Unknown reservoir type & negative mean outflow
            // Only look at conditions of current year
            for (n = 0; n < ng; n++) {
                if ((G_reservoir_area[n] > 0.) && (resYearCurrentToUse >= G_res_start_year[n])) {

                    if (G_res_type[n] == 0) {//only IF reservoir type is UNKNOWN (should not be the case anyhow)
                        cerr << "ERROR: routing.annuaInit: Year " << year
                             << " - Reservoir management: G_res_type is unknown at n: " << n << endl;
                        cerr << "Reservoir will be treated as global lake!" << endl;
                        cerr << "Adding " << G_reservoir_area_full[n] << " km2 to global lake area " << G_lake_area[n]
                             << endl;
                        G_lake_area[n] += G_reservoir_area_full[n]; //use reservoir area as lake area, or sum up reservoir and lake area
                        G_reservoir_area[n] = 0.; //set reservoir area to zero (inhibits further accounting elsewhere)
                        G_reservoir_area_full[n] = 0.; //set full reservoir area to zero (inhibits further accounting in next year)
                        cerr << "Final global lake area " << G_lake_area[n] << endl;
                    } else if (G_mean_outflow[n] <= 0.) {
                        cerr << "ERROR: routing.annuaInit: Year " << year
                             << " - Reservoir management: G_mean_outflow in reservoir is " << G_mean_outflow[n]
                             << " i.e. less than or equal to zero at n: " << n << endl;
                        cerr << "Reservoir will be treated as global lake!" << endl;
                        cerr << "Adding " << G_reservoir_area_full[n] << " km2 to global lake area " << G_lake_area[n]
                             << endl;
                        G_lake_area[n] += G_reservoir_area_full[n]; //use reservoir area as lake area, or sum up reservoir and lake area
                        G_reservoir_area[n] = 0.; //set reservoir area to zero (inhibits further accounting elsewhere)
                        G_reservoir_area_full[n] = 0.; //set full reservoir area to zero (inhibits further accounting in next year)
                        cerr << "Final global lake area " << G_lake_area[n] << endl;

                    }

                }
                // end if option conditions
            }
            // end loop over grid cells

        }
        // end if options.resOpt == 1

    }
    //end calculate only for anthropogenic runs (options.antNatOpt == 0)


    // Count the number of cells with reservoir area
    short count = 0;
    for (n = 0; n < ng; n++) {
        if (G_reservoir_area[n] > 0.) count++;
    }
    cout << "routing.annualInit: " << count << " reservoirs (area > 0) are considered in this year " << year << endl;

    //
    // FP20161018N002 Reservoir operation start years
    //


    // Part 3: Always used arrays that depend on G_glo_res (Unit: percent fraction of grid cell area)
    for (n = 0; n < ng; n++) {
        // Now, calculate land area fraction (instead of once at model initialisation at routing.init because G_glo_res can be dynamic)

        // (1) Land area fraction excluding lakes and reservoirs/regulated lakes (G_glo_res includes G_reg_lake)
        // fixed G_fGloLake = G_glo_lake/100 and newly read-in G_glo_res
        G_landWaterExclGloLakAreaFrac[n] = /*(double)*/ (geo.G_contfreq[n]
                                                         - G_glo_lake[n] - G_glo_res[n]);
        if (G_landWaterExclGloLakAreaFrac[n] < limit_errorstrange_landAreaFrac_pct) {
            cerr << "Warning: routing.annualInit: Year " << year
                 << " - Strange number for lake/wetland fraction 'G_landWaterExclGloLakAreaFrac' for cell " << n << ": "
                 << G_landWaterExclGloLakAreaFrac[n] << ", limit: " << limit_errorstrange_landAreaFrac_pct << endl;
            cerr << "Warning: routing.annualInit: Year " << year << " - Set values < 0 to zero" << endl;
        }
        if (G_landWaterExclGloLakAreaFrac[n] < 0.)
            G_landWaterExclGloLakAreaFrac[n] = 0.;

        // (2)
        // (2.1) Fraction of water and land: G_fwaterfreq and G_landfreq (model start and updatedevery year)
        // and associated changes
        //G_fwaterfreq is sum of G_glo_lake, G_loc_lake, and G_glo_res (includes G_reg_lake)
        //G_landfreq is G_cont - G_fwaterfreq, contains wetlands (is used to calculate PET in daily.cpp)
        // ATTENTION: Any change to previous year is assumed to be an increase in G_glo_res

        // (2.2) G_landAreaFrac for the first year = at model start
        // For first year

        if (0 == statusStarted_landfreq_fwaterfreq_landAreaFrac[n]) {

            // Water and area fractions (in percent)
            // fixed G_fLocLake from first reading in routing.init: G_loc_lake/100
            // fixed G_fGloLake = G_glo_lake/100
            // as read-in G_glo_res
            G_fwaterfreq[n] = G_glo_lake[n] + G_loc_lake[n] + G_glo_res[n];
            G_landfreq[n] = geo.G_contfreq[n] - G_fwaterfreq[n];


            G_fwaterfreq_change[n] = 0.;
            G_landfreq_change[n] = 0.;

            // Initial land area fraction (in percent)
            G_landAreaFrac[n] = /*(double)*/ (geo.G_contfreq[n]
                                              - (G_glo_lake[n] + G_glo_wetland[n] + G_loc_lake[n] + G_loc_wetland[n] +
                                                 G_glo_res[n]));

            /*    if (n == 24026) {
                cout << "++++ first day: G_fwaterfreq[n] " << G_fwaterfreq[n] << endl;
                cout << "++++ first day: G_landfreq[n] " << G_landfreq[n] << endl;
                cout << "++++ first day: geo.G_contfreq[n] " << geo.G_contfreq[n] << endl;
                cout << "++++ first day: G_fwaterfreq[n] " << G_fwaterfreq[n] << endl;
                cout << "++++ first day: G_landAreaFrac[n] " << G_landAreaFrac[n] << endl;
            }*/
            // Treat undesired values
            if (G_landAreaFrac[n] < limit_errorstrange_landAreaFrac_pct) {
                cerr << "Warning: routing.annualInit: Year " << year
                     << " - Strange number for land area fraction 'G_landAreaFrac' for cell " << n << ": "
                     << G_landAreaFrac[n] << ", limit: " << limit_errorstrange_landAreaFrac_pct << endl;
                cerr << "Warning: routing.annualInit: Year " << year << " - Set values < 0 to zero" << endl;
            }
            if (G_landAreaFrac[n] < 0.)
                G_landAreaFrac[n] = 0.;


            // update indicator
            statusStarted_landfreq_fwaterfreq_landAreaFrac[n] = 1;

        }
            // For all years besides the first year, transfer data from previous year and recalculate
        else {

            G_fwaterfreq_prevyear[n] = G_fwaterfreq[n];
            G_landfreq_prevyear[n] = G_landfreq[n];

            // Consider fractional routing fractions of local lakes (as done in previous step, besides for the first year)
            // variable G_fLocLake from previous time step
            // fixed G_fGloLake = G_glo_lake/100
            // as read-in G_glo_res
            G_fwaterfreq[n] = G_glo_lake[n] + (G_fLocLake[n] * 100.) + G_glo_res[n];
            G_landfreq[n] = geo.G_contfreq[n] - G_fwaterfreq[n];

            G_fwaterfreq_change[n] = G_fwaterfreq[n] - G_fwaterfreq_prevyear[n];
            G_landfreq_change[n] = G_landfreq[n] - G_landfreq_prevyear[n];

        }
        // end else



    }
    // end loop over grid cells


    // Part 4: Needed only in case of anthropogenic runs using reservoirs
    if (options.antNatOpt == 0) { // anthropogenic runs: consider reservoirs and regulated lakes as reservoirs
        //    cout << "routing.annualInit Part 4: if options.antNatOpt == 0, anthropogenic run: consider reservoirs and regulated lakes as reservoirs" << endl;

        // Anthropogenic run
        // Use reservoirs and regulated lakes
        if (options.resOpt == 1) {

            // HERE: Adapt code in order to add water from reduced land fraction to reservoir volume
            // include also G_riverAreaFrac[n]
            for (n = 0; n < ng; n++) {

                //check if change arrays etc. are already initialized
                //0 == statusStarted_modelStartConditions[n] (initial condition) at routing.init, i.e. start of program execution
                /*     if (0 == statusStarted_modelStartConditions[n]) {
                 G_glores_prevyear[n] = G_glo_res[n]; //initialize with percent fraction of reservoirs of start or reference year
                 // Force other changes to zero
                 G_glores_change[n] = 0.;
                 G_fwaterfreq_change[n] = 0.;
                 G_landfreq_change[n] = 0.;
                 G_landAreaFrac_change[n] = 0.;
                 statusStarted_modelStartConditions[n] = 1; // indicate for the next iteration that first year has been treated / model has started
              }
              else {
                 G_glores_change[n] = G_glo_res[n] - G_glores_prevyear[n];
                 G_fwaterfreq_change[n] = G_fwaterfreq[n] - G_fwaterfreq_prevyear[n];
                 G_landfreq_change[n] = G_landfreq[n] - G_landfreq_prevyear[n];

              }
        */

                // Calculate (possible) increase in reservoir area at start of the year
                // indicating operation start of a new reservoir
                G_glores_change[n] = G_glo_res[n] - G_glores_prevyear[n];


                // Perform changes in water storages only if a new reservoir started operation
                if (G_glores_change[n] > 0.) {

                    // Subtract additional reservoir fraction from land area fraction
                    // (1) Land area fraction including wetlands, local lakes, rivers
                    // day 0 (current day, in current year: year / Jan / 01): G_glo_res[n]
                    // From routing.updateLandAreaFrac:
                    // day-1 (previous day, in previous year: year-1 / Dec / 31): G_landAreaFrac[n]
                    // day-2 (day before previous day, in previous year: year-1 / Dec / 30): G_landAreaFracPrevTimestep[n]
                //    if (n == 22473) {
                //        cout << "YEAR " << year << endl;
                //        cout << "G_landAreaFrac before adapt " << G_landAreaFrac[n] << endl;
                //        cout << "G_glo_res[n] before adapt " << G_glo_res[n] << endl;
                //    }
                    if (G_glores_prevyear[n] > 0.) { // in case there was already a reservoir fraction, adapt only the changes.
                        G_landAreaFrac[n] = G_landAreaFrac[n] - G_glores_change[n];
                    }
                    else { //reservoir fraction increased from 0.
                        G_landAreaFrac[n] = G_landAreaFrac[n] - G_glo_res[n];
                    }
                //    if (n == 22473) {
                //        cout << "G_landAreaFrac after adapt " << G_landAreaFrac[n] << endl;
                //    }
                    // Treat undesired values
                    if (G_landAreaFrac[n] < limit_errorstrange_landAreaFrac_pct) {
                        cerr
                                       << "Warning: routing.annualInit: (options.antNatOpt == 0) & (options.resOpt == 1) & (G_glores_change[n] > 0.) - Strange number for land area fraction 'G_landAreaFrac' for cell "
                                       << n << ": " << G_landAreaFrac[n] << endl;
                    }
                    if (G_landAreaFrac[n] < 0.)
                        G_landAreaFrac[n] = 0.;


                    // START ISIMIP2b POSSIBLE uncomment 2016.
                    // ATTENTION: POSSIBLY For ISIMIP2B, GLWD units consistent to outflow cells are missing.
                    // Therefore, the calculation of net change and subsequent transfer of storage to reservoir storage is commented out
                    // Scaling of storage as usual in daily.calcNewDay
                    //

                    // Calculate net change
                    // Changes of day 0 from wetland / loc.lakes/rivers (fractional routing) + new reservoir
                    // with respect to day-1
                    G_landAreaFrac_change[n] = G_landAreaFrac[n] - G_landAreaFracPrevTimestep[n];
                    /*  if (n == 40812) {
                        cout << "G_landAreaFrac_change[n] " << G_landAreaFrac_change[n] << endl;
                        cout << "G_landAreaFrac[n] " << G_landAreaFrac[n] << endl;
                        cout << "G_landAreaFracPrevTimestep[n] " << G_landAreaFracPrevTimestep[n] << endl;
                    }*/
                    // In case that net change is not negative,
                    // assume that it is a result of fractional routing only, reset G_landAreaFrac to value of day-2
                    // and let daily.calcNewDay execute as usual, implicitly applying no changes (from fractional routing
                    if (G_landAreaFrac_change[n] >= 0.) {
                        G_landAreaFrac[n] = G_landAreaFracPrevTimestep[n];
                        G_landAreaFrac_change[n] = 0.;
                    }
                        // Net change is negative (reservoir and/or added: Adapt storage
                    else {

                        // Calculate absolute storage volume in km3 from reduced ("lost") land surface fraction
                        canopyWaterContent_change_km3 =
                                       dailyWaterBalance.G_canopyWaterContent[n] * geo.areaOfCellByArrayPos(n) /
                                       1000000.0 * ((-G_landAreaFrac_change[n]) / 100.0);
                        soilWaterContent_change_km3 =
                                       dailyWaterBalance.G_soilWaterContent[n] * geo.areaOfCellByArrayPos(n) /
                                       1000000.0 * ((-G_landAreaFrac_change[n]) / 100.0);
                        snowWaterContent_change_km3 =
                                       dailyWaterBalance.G_snow[n] * geo.areaOfCellByArrayPos(n) / 1000000.0 *
                                       ((-G_landAreaFrac_change[n]) / 100.0);
                        SnowInElevation_n_change_km3 = 0.;
                        for (short elev = 0; elev < 101; elev++) {
                            SnowInElevation_elev_change_km3 = dailyWaterBalance.G_SnowInElevation[n][elev] / 100.0 *
                                                              geo.areaOfCellByArrayPos(n) / 1000000.0 *
                                                              ((-G_landAreaFrac_change[n]) / 100.0);
                            SnowInElevation_n_change_km3 += SnowInElevation_elev_change_km3;
                        }
                        if (((SnowInElevation_n_change_km3 - snowWaterContent_change_km3) >
                             limit_errorstrange_snowWaterContent_change_km3)
                            || ((snowWaterContent_change_km3 - SnowInElevation_n_change_km3) >
                                limit_errorstrange_snowWaterContent_change_km3)) {
                            cerr
                                           << "Warning: routing.annualInit: (options.antNatOpt == 0) & (options.resOpt == 1) & (G_glores_change[n] > 0.) - Exceeded limit of difference between 'SnowInElevation_n_change_km3' and 'snowWaterContent_change_km3' for cell "
                                           << n << ": " << SnowInElevation_n_change_km3 << ", "
                                           << snowWaterContent_change_km3 << ", limit: "
                                           << limit_errorstrange_snowWaterContent_change_km3 << endl;
                        }
                        //dailyWaterBalance.G_dailyStorageTransfer[n] = 0.;

                        // Add storage to reservoir storage volume located in outflow cell, as identified via GLWDunit grid

                        G_gloResStorage[G_outflow_cell_assignment[n] - 1] +=
                                       canopyWaterContent_change_km3 + soilWaterContent_change_km3 +
                                       snowWaterContent_change_km3; //

                        // ATTENTION: dailyWaterBalance.G_canopyWaterContent etc. are in units of Millimeter (mm)
                        // Through the reduction of land surface fraction, water storage is implicitly reduced
                        // To ensure that no further scaling is done within daily.calcNewDay, set equal previous value (day-1) to current newly determined uvalue (day 0)
                        /*  if (n == 40812) {
                            cout << "G_landAreaFracPrevTimestep[n] " << G_landAreaFracPrevTimestep[n] << endl;
                            cout << "G_landAreaFrac[n] " << G_landAreaFrac[n] << endl;
                        }*/
                        G_landAreaFracPrevTimestep[n] = G_landAreaFrac[n];
                        /*   if (n == 40812) {
                            cout << "G_landAreaFracPrevTimestep[n] " << G_landAreaFracPrevTimestep[n] << endl;
                            cout << "G_landAreaFrac[n] " << G_landAreaFrac[n] << endl;
                        }*/
                    }
                    // end else

                    //
                    // END ISIMIP2b uncomment


                }
                // end if (G_glores_change[n] > 0.)
            }
            // end loop over grid cells

        }
        // end if options.resOpt == 1

    }
    //end calculate only for anthropogenic runs (options.antNatOpt == 0)
// FP20161018N002: Do NOT introduce nSpecBasins as local variable nBasins, as superbasin area calculation probably wrong & obsolete
//    // Superbasin calculation currently obsolete, probably not working correctly in WaterGAP2.2b, do in postprocessing
//    // calculate land area of basins
//    extern signed short G_sbasin[ng];
//    basinLandArea = new double[nBasins];
//
//    for (int n = 0; n < nBasins; ++n)
//        basinLandArea[n] = 0.;
//    for (int n = 0; n < ng; ++n) {
//        if (G_sbasin[n] > 0) {
//            basinLandArea[G_sbasin[n] - 1] += G_landAreaFrac[n] * geo.areaOfCellByArrayPos(n) / 100.0;
//        }
//    }
//    // end loop over grid cells

    // END ISIMIP21a insert ersion


    // Initialize arrays for writing output
    if ((2 == options.grid_store) || (4 == options.grid_store) || (6 == options.grid_store)) {
        for (n = 0; n < ng; n++)
            for (m = 0; m < 12; m++) {
                G_monthlyConsistentPrecip[n][m] = 0.;

                G_monthlyRiverAvail[n][m] = 0.;
                G_monthlyPotCellRunoff[n][m] = 0.;
                G_monthlyPotCellRunoffDeficit[n][m] = 0.;
                G_monthlyCellRunoff[n][m] = 0.;
                G_monthlyCellSurfaceRunoff[n][m] = 0.;
                // added for cell AET calculation (2.1f)
                G_monthlyCellAET[n][m] = 0.;
                G_monthlyCellAETWCa[n][m] = 0.;
                G_monthlyOpenWaterEvap[n][m] = 0.;
                G_monthlyVelocity[n][m] = 0.;
                G_monthlySurfStor[n][m] = 0.;
                G_monthlySurfStor_mm[n][m] = 0.;
                G_monthlyMinRiverAvail[n][m] = 0.;
                G_monthlyMaxRiverAvail[n][m] = 0.;
                G_monthlyGwrSwb[n][m] = 0.;
                G_monthlyLandAreaFrac[n][m] = 0.;
                G_monthlyFswb[n][m] = 0.;
                G_monthlyGwRunoff[n][m] = 0.;
                G_monthlyLocWetlExtent[n][m] = 0.;
                G_monthlyGloWetlExtent[n][m] = 0.;
                if (options.outRiverPET == 1) G_monthlyRiverAreaFrac[n][m] = 0.;
                if (options.outRiverPET == 1) G_monthlyRiverPET[n][m] = 0.;

                if (options.resOpt == 1) {
                    G_monthlyResStorage[n][m] = 0.;
                    G_monthlyResStorageRatio[n][m] = 0.;
                    G_monthlyResInflow[n][m] = 0.;
                    G_monthlyResOutflow[n][m] = 0.;
                }

                if (options.outRiverInUpstream)
                    G_monthlyRiverInUpstream[n][m] = 0.;
                G_monthlyRiverStorage[n][m] = 0.;
                G_monthlyGwStorage[n][m] = 0.; // moved from daily.cpp
                G_monthlyLocLakeStorage[n][m] = 0.;
                G_monthlyLocWetlStorage[n][m] = 0.;
                G_monthlyGloLakeStorage[n][m] = 0.;
                G_monthlyGloWetlStorage[n][m] = 0.;
                G_monthlySatisfiedUse[n][m] = 0.;
                G_monthlyActualUse[n][m] = 0.;
                G_monthlyNUg[n][m] = 0.;
                // CR 150819 test output
                G_monthlyAllocatedUse[n][m] = 0.; // FP 2015, not used anymore, at the moment, implement if needed again
                G_monthlySatisAllocatedUseinSecondCell[n][m] = 0.;
            }
        // end loop over months


    }
    // end if options (for writing grid output)

}
// end annualInit

// HMS 2015-03-02 I do not understand why this is called monthlyInit - should it be named better daily31outInit? (same in daily.cpp)
// (FP 2015: perhaps naming is because the data arrays contain data from one month)
// daily output option 31 start
void routingClass::daily31outInit() {
    short d;
    int n;

    //if ((3 == options.grid_store) || (4 == options.grid_store)) {
    for (n = 0; n < ng; n++)
        for (d = 0; d < 31; d++) {
            if (options.outPrecDaily) G_daily31ConsistentPrecip[n][d] = 0.;

            if (options.outCellRunoffDaily) G_daily31CellRunoff[n][d] = 0.;
            if (options.outCellSurfaceDaily) G_daily31CellSurfaceRunoff[n][d] = 0.;
            if (options.outGWRunoffDaily) G_daily31GwRunoff[n][d] = 0.;
            if (options.outRiverAvailDaily) G_daily31RiverAvail[n][d] = 0.;
            if (options.outRiverVeloDaily) G_daily31Velocity[n][d] = 0.;
            if (options.outCellAETDaily) G_daily31CellAET[n][d] = 0.;
            if (options.outSurfStorDaily) G_daily31SurfStor[n][d] = 0.;
            if (options.outSingleStoragesDaily) G_daily31LocLakeStor[n][d] = 0.;
            if (options.outSingleStoragesDaily) G_daily31LocWetlStor[n][d] = 0.;
            if (options.outSingleStoragesDaily) G_daily31GloLakeStor[n][d] = 0.;
            if (options.outSingleStoragesDaily) G_daily31GloWetlStor[n][d] = 0.;
            if (options.outSingleStoragesDaily) G_daily31RiverStor[n][d] = 0.;
            if (options.outSingleStoragesDaily) G_daily31ResStor[n][d] = 0.;
            if (options.outSingleStoragesDaily || options.outGWStorageDaily) G_daily31GwStor[n][d] = 0.;
            if (options.outTotalWaterInStoragesDaily_km3 ||
                options.outTotalWaterInStoragesDaily_mm)
                G_daily31TotalWaterInStorages_km3[n][d] = 0.;
            if (options.outTotalWaterInStoragesDaily_mm) G_daily31TotalWaterInStorages_mm[n][d] = 0.;
            if (options.outGwrSwbDaily) G_daily31GwrSwb[n][d] = 0.;
            if (options.outFswbDaily) G_daily31Fswb[n][d] = 0.;
            if (options.outLandAreaFracDaily) G_daily31LandAreaFrac[n][d] = 0.;
            if (options.outGwrunSurfrunDaily) G_daily31GwrunSurfrun[n][d] = 0.;
            if (options.outCellAETWCaDaily) G_daily31CellAETWCa[n][d] = 0.;
        }
    // endfor loop over 31 days
}
//}
// end daily31outInit
// daily output option 31 end

// HMS 2017_02_21 to store daily TWS value for first and last day in year.
void routingClass::startendOutInit() {
    short d;
    int n;

    for (n = 0; n < ng; n++)
        for (d = 0; d < 1; d++) {
            G_startendTotalWaterInStorages_km3[n][d] = 0.;
        }
}

// daily output option 365 start
void routingClass::daily365outInit() {  // HMS 2015-03-02 where is this class called (not in watergap.cpp)
    short d;
    int n;

    //if ((5 == options.grid_store) || (6 == options.grid_store)) {
    for (n = 0; n < ng; n++)
        for (d = 0; d < 365; d++) {

//CR #############################################################################################
//	G_daily365ActualUse[n][d]				= 0.;
//	G_daily365AllocUse[n][d]				= 0.;
//	G_daily365AllocUseNextDay[n][d]			= 0.;
//	G_daily365TotalUnsatUse[n][d]			= 0.;
//	G_daily365UnsatAllocUseNextDay[n][d]	= 0.;
//CR #############################################################################################

            if ((options.outPrecDaily) ||
                ((options.scoutcPrecip) && (2 == options.day_store)))
                G_daily365ConsistentPrecip[n][d] = 0.;
            if ((options.outCellRunoffDaily) ||
                (((options.scoutCellRunoffkm3) || (options.scoutCellRunoffmm)) && (2 == options.day_store)))
                G_daily365CellRunoff[n][d] = 0.;
            if ((options.outCellRunoffDaily) ||
                ((options.scoutCellRunoffmm) && (2 == options.day_store)))
                G_daily365CellRunoff_mm[n][d] = 0.;
            if ((options.outCellSurfaceDaily) ||
                ((options.scoutCellSRunoff) && (2 == options.day_store)))
                G_daily365CellSurfaceRunoff[n][d] = 0.; // FP20160915N002
            if ((options.outLandAreaFracDaily) ||
                ((options.scoutLandAreaFrac) && (2 == options.day_store)))
                G_daily365LandAreaFrac[n][d] = 0.;
            if ((options.outGWRunoffDaily) || ((options.scoutGwRunoff) && (2 == options.day_store)))
                G_daily365GwRunoff[n][d] = 0.;
            if ((options.outRiverAvailDaily) || ((options.scoutQ) && (2 == options.day_store)))
                G_daily365RiverAvail[n][d] = 0.;
            if ((options.outRiverVeloDaily) || ((options.scoutFlowVelo) && (2 == options.day_store)))
                G_daily365Velocity[n][d] = 0.;
            if ((options.outCellAETDaily) || ((options.scoutCellAET) && (2 == options.day_store)))
                G_daily365CellAET[n][d] = 0.;
            if ((options.outSurfStorDaily) ||
                ((options.scoutSurfaceStor) && (2 == options.day_store)))
                G_daily365SurfStor[n][d] = 0.;
            if ((options.outSingleStoragesDaily) ||
                ((options.scoutLocLake) && (2 == options.day_store)))
                G_daily365LocLakeStor[n][d] = 0.;
            if ((options.outSingleStoragesDaily) ||
                ((options.scoutLocWet) && (2 == options.day_store)))
                G_daily365LocWetlStor[n][d] = 0.;
            if ((options.outSingleStoragesDaily) ||
                ((options.scoutGloLake) && (2 == options.day_store)))
                G_daily365GloLakeStor[n][d] = 0.;
            if ((options.outSingleStoragesDaily) ||
                ((options.scoutGloWet) && (2 == options.day_store)))
                G_daily365GloWetlStor[n][d] = 0.;
            if ((options.outSingleStoragesDaily) ||
                ((options.scoutRiver) && (2 == options.day_store)))
                G_daily365RiverStor[n][d] = 0.;
            if ((options.outSingleStoragesDaily) ||
                ((options.scoutReservoir) && (2 == options.day_store)))
                G_daily365ResStor[n][d] = 0.;
            if ((options.outSingleStoragesDaily || options.outGWStorageDaily) ||
                ((options.scoutGwStor) && (2 == options.day_store)))
                G_daily365GwStor[n][d] = 0.;
            if ((options.outTotalWaterInStoragesDaily_km3 || options.outTotalWaterInStoragesDaily_mm) ||
                (((options.scoutTWSkm3) || (options.scoutTWSmm)) && (2 == options.day_store)))
                G_daily365TotalWaterInStorages_km3[n][d] = 0.;
            if ((options.outTotalWaterInStoragesDaily_mm) ||
                ((options.scoutTWSmm) && (2 == options.day_store)))
                G_daily365TotalWaterInStorages_mm[n][d] = 0.;
            if ((options.outGwrSwbDaily) || ((options.scoutGwrSwb) && (2 == options.day_store)))
                G_daily365GwrSwb[n][d] = 0.;
            if ((options.outFswbDaily) || ((options.scoutFswb) && (2 == options.day_store))) G_daily365Fswb[n][d] = 0.;
            if (options.outGwrunSurfrunDaily) G_daily365GwrunSurfrun[n][d] = 0.;
            if (options.outCellAETWCaDaily) G_daily365CellAETWCa[n][d] = 0.;
        }
    // endfor loop over 365 days


}
//}
// end daily365outInit
// daily output option 365 end

void routingClass::routing(const short year, short day, short month, const short day_in_month) {
    // routing through river network
    extern dailyWaterBalanceClass dailyWaterBalance;

    extern short G_toBeCalculated[ng];

    const short int first_day_in_month[12] = {1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335};

    double G_localRunoff[ng];
    double G_localRunoffIntoRiver[ng]; // for options.fractionalRouting
    double G_localGWRunoff[ng]; // needed for different values for sub-/surface runoff
    double G_localGWRunoffIntoRiver[ng]; // for options.fractionalRouting
    // initialize
    for (int n = 0; n < ng; n++) {
        G_localRunoff[n] = 0.;
        G_localRunoffIntoRiver[n] = 0.;
        G_localGWRunoff[n] = 0.;
        G_localGWRunoffIntoRiver[n] = 0.;
    }
    double CellRunoff = 0.; // actual cell runoff generated within the cell
    double cellArea = 0., transportedVolume = 0., dailyRiverEvapo = 0.;
    double inflow = 0., inflowUpstream = 0., totalInflow = 0., outflow = 0., maxStorage = 0.;
    double K = 0.; // for new routing-equation
    short b = 0; // counters for nSpecBasins (avoid dynamic memory allocation)
    double potCellRunoff = 0.;

    double Kswbgw = 10.; // constant groundwater recharge below surface water bodies in mm/d --> 10 mm per day as gwr

    double fswbFracCorr = 0.; // correction for swb (loclak, locwet, glowet) if riverAreaFrac can not be satisfied
    double riverAreaFracDeficit = 0.; //to correct river area fraction if necessary
    extern short G_aindex[ng]; // to get info about arid cells in a first approach


    // open water PET reduction (2.1f)
    double locLakeAreaReductionFactor = 0., gloLakeEvapoReductionFactor = 0.,
                   locWetlAreaReductionFactor = 0., gloWetlAreaReductionFactor = 0.;
    double gloResEvapoReductionFactor = 0., riverAreaReductionFactor = 0.;
    //CR: introduced in 2.2b; openwaterPET, gwr, and RemainingUse are included in one storage equation.
    //If they cannot be fully satisfied, they are reduced by the ReductionFactor = PETgwrMax/PETgwr.
    double PETgwr = 0., PETgwrMax = 0., PETgwrRedFactor = 0., PETgwrRemUseRedFactor = 0.;
    double groundwater_runoff_km3 = 0.;
    double G_groundwaterStoragePrevStep = 0.;
    // WG22b: GWR beneath SWB [km3] is included in the GW storage equation of the next timestep
    double netGWin = 0.;
    double dailyLocalSurfaceRunoff = 0.;    // in km3

    // set arrays with daily values to zero
    // this has to be done already here because
    // these arrays are summed up with finer time steps later
    // Reservoir operation start years
    // Homogenized index n to general b for "nSpecBasins"
    for (b = 0; b < nSpecBasins; b++) {
        dailyRiverDischarge[b][day - 1] = 0.;
        if (options.day_store == 1) {
            dailyLocLakeStorage[b][day - 1] = 0.;
            dailyLocWetlStorage[b][day - 1] = 0.;
            dailyGloLakeStorage[b][day - 1] = 0.;
            dailyGloWetlStorage[b][day - 1] = 0.;
            dailyRiverVelocity[b][day - 1] = 0.;    // velocity
            if (options.resOpt == 1)
                dailyResStorage[b][day - 1] = 0.;
        }
        // end if options.day_store == 1
    }
    // end loop over nSpecBasins


#pragma omp parallel for private() \
                shared(G_localRunoff, G_localGWRunoff, G_dailyLocalSurfaceRunoff, G_dailyStorageTransfer, dailyWaterBalance, geo, options, month)

    // computation of GW storage and baseflow for options.aridareaOpt == 0
    if (0 == options.aridareaOpt) {
        // Loop over grid cells
        for (int n = 0; n < ng; n++) {

            // Only execute for continental grid cells
            if (geo.G_contcell[n]) {
                // Efficiency enhancement through local copies instead of often used method getValue // FP
                double P_GWOUTF_C__at__n_aridareaOpt = calibParam.getValue(P_GWOUTF_C, n);

                // transform runoff into desired units
                cellArea = geo.areaOfCellByArrayPos(n);

                G_groundwaterStoragePrevStep = G_groundwaterStorage[n];

                // WG22b: ODE (ordinary differential equation) is applied for groundwater storage: dS/dt = GWR - NAg - k*S  (k*S = groundwater_runoff_km3)

                netGWin = dailyWaterBalance.G_dailyGwRecharge[n] * cellArea * (G_landAreaFrac[n] / 100.) /
                          1000000.;        // mm -> km3
                // If G_landAreaFrac[n] == 0 in this timestep, G_dailyGwRecharge[n] has been set to 0 in daily.cpp.

                // G_dailydailyNUg is adapted in this timestep based on the remaining use of last time step
                if (options.subtract_use > 0) {
                    if (G_dailyRemainingUse[n] != 0) {
                        G_dailydailyNUg[n] = updateNetAbstractionGW(n, month);
                    }
                    netGWin -= G_dailydailyNUg[n];
                }

                // Analytical solution of dS/dt = netGWR - k*S
                // k_g in [1/d]
                // Cell-specific calibration parameters - Use parameter  // FP

                G_groundwaterStorage[n] = G_groundwaterStoragePrevStep * exp(-1. * P_GWOUTF_C__at__n_aridareaOpt) +
                                          (1. / P_GWOUTF_C__at__n_aridareaOpt) * netGWin *
                                          (1. - exp(-1. * P_GWOUTF_C__at__n_aridareaOpt));

                groundwater_runoff_km3 = G_groundwaterStoragePrevStep - G_groundwaterStorage[n] + netGWin;

                if (groundwater_runoff_km3 <=
                    0.) {    // different differential equation (ODE) applies: dS/dt = GWR - NAg (without k*S) -> S(t) = S(t-1) + netGWin

                    groundwater_runoff_km3 = 0.;

                    G_groundwaterStorage[n] = G_groundwaterStoragePrevStep + netGWin;
                }
                if ((2 == options.grid_store) || (4 == options.grid_store) || (6 == options.grid_store))
                    G_monthlyGwRunoff[n][month] += groundwater_runoff_km3;
                //if daily 365 output and gwout (same for 31; same for land area fraction)
                //G_daily365GwRunoff[n][day - 1] = groundwater_runoff_km3;

                // CR 2015-08: If landAreaFraction == 0., the remaining canopy, snow, and soil storage from the previous timestep becomes surface runoff (see daily.cpp).
                if (G_landAreaFrac[n] == 0.) {
                    dailyLocalSurfaceRunoff = dailyWaterBalance.G_dailyStorageTransfer[n] * cellArea / 1000000. *
                                              G_landAreaFracPrevTimestep[n] / 100.;
                } else {
                    dailyLocalSurfaceRunoff = dailyWaterBalance.G_dailyLocalSurfaceRunoff[n] * cellArea / 1000000. *
                                              G_landAreaFrac[n] / 100.;
                }


                if (1 ==
                    options.fractionalRoutingOpt) { // only a fraction (depending on amount of surface water bodies) of surface and groundwater runoff should go to the local lakes later, rest should go directly to the river
                    G_fswb_catchment[n] = G_fswbInit[n] *
                                          20.;                  // G_fswb_catchment[n] = G_fswb[n] * 5.; // add a catchment area to surface water bodies (specified as x times the size of swb)...
                    if (G_fswb_catchment[n] > 1.) // ... and prevent that G_fswb is greater then 1
                        G_fswb_catchment[n] = 1.;
                    if (G_aindex[n] == 1) {// in semi-arid/arid areas, all groundwater reaches the river directly
                        G_localGWRunoffIntoRiver[n] = groundwater_runoff_km3;
                    } else // in humid areas, groundwater is also routed through surface water bodies
                        G_localGWRunoffIntoRiver[n] = (1. - G_fswb_catchment[n]) * groundwater_runoff_km3;

                    if (G_aindex[n] == 0) // route groundwater only through surface water bodies in humid regions
                        G_localGWRunoff[n] = G_fswb_catchment[n] * groundwater_runoff_km3;

                    else
                        G_localGWRunoff[n] = 0.; // in arid regions, all groundwater flows directly to the river

                    // Firstly, amount of local runoff flowing directly into the river is calculated, then local gw runoff is added for the rest term.
                    //G_localRunoffIntoRiver[n] = (1. - G_fswb_catchment[n]) * dailyWaterBalance.G_dailyLocalSurfaceRunoff[n]
                    //        * cellArea / 1000000. * G_landAreaFrac[n] / 100.;
                    G_localRunoffIntoRiver[n] = (1. - G_fswb_catchment[n]) *
                                                dailyLocalSurfaceRunoff;    // dailyLocalSurfaceRunoff in km3

                    G_localRunoff[n] = G_fswb_catchment[n] * dailyLocalSurfaceRunoff + G_localGWRunoff[n];


                } else {    //normal version (without fractional routing)
                    // G_dailyLocalRunoff is devided in G_dailyLocalGWRunoff and G_dailyLocalSurfaceRunoff
                    // values are for land area only

                    G_localGWRunoff[n] = groundwater_runoff_km3;

                    G_localRunoff[n] = dailyLocalSurfaceRunoff + G_localGWRunoff[n];

                }

                // river storage is calculated in km3 instead of mm
                // this makes transport between cells much more easy.
                // therefore area/1000000.0
                //
                // WG22b: timeStepsPerDay removed from equations.

                // the division by timeStepsPerDay is done*/
                // so that this array contains the value which has
                // to be added during each time step of
                // the routing algorithm.
                //

                // dailyWaterBalance.G_lakeBalance[n] /= timeStepsPerDay; // HMS 2014-06-04 G_lakeBalance reintroduced
                // added for open water PET reduction (2.1f)
                // not needed anymore (as timeStepsPerDay is now 1 (analytical solution of ODE)
                // dailyWaterBalance.G_openWaterPrec[n] /= (double) timeStepsPerDay;
                // dailyWaterBalance.G_openWaterPET[n] /= (double) timeStepsPerDay;

                // to be changed for fractionalrouting option (normal routing without G_localRunffIntoRiver and G_localRGWRunoffIntoRiver
                // reintroduced for CFA calculation
                potCellRunoff = (G_localRunoff[n] + G_localRunoffIntoRiver[n] +
                                 G_localGWRunoffIntoRiver[n]    // is in km3
                                 + ((dailyWaterBalance.G_lakeBalance[n]    // is in mm
                                     / 1000000.0)    // mm -> km3
                                    * (((double) G_lake_area[n] + (double) G_reservoir_area[n])    // km2
                                       + ((geo.areaOfCellByArrayPos(n) / 100.0)    // % -> km2
                                          * (G_loc_lake[n] + G_loc_wetland[n] + G_glo_wetland[n])))));    // [%]

                G_potCellRunoff[n] += potCellRunoff;
                if (2 == options.grid_store) {
                    G_monthlyPotCellRunoff[n][month] += potCellRunoff;
                }




            }
            // endif (geo.G_contcell[n]) [only execute for continental cells]

        } // end loop over all cells [n]

    } // end of computation of GW storage and baseflow for options.aridareaOpt == 0

#ifdef _WATERGAP_CHECKS_GLOBAL_H
                                                                                                                            // Check Felix Portmann 2015 - for invalid data
        if ( (year >= chk_year) && (month >= chk_month) && (day_in_month >= chk_day_in_month) ) {

            check_array_double_LT0_YYMMDD(G_localRunoff, ng, year, month, day_in_month);
            check_array_double_LT0_YYMMDD(G_localRunoffIntoRiver, ng, year, month, day_in_month);
            check_array_double_LT0_YYMMDD(G_localGWRunoff, ng, year, month, day_in_month);
            check_array_double_LT0_YYMMDD(G_localGWRunoffIntoRiver, ng, year, month, day_in_month);

            check_array_double_NaN_YYMMDD(G_localRunoff, ng, year, month, day_in_month);
            check_array_double_NaN_YYMMDD(G_localRunoffIntoRiver, ng, year, month, day_in_month);
            check_array_double_NaN_YYMMDD(G_localGWRunoff, ng, year, month, day_in_month);
            check_array_double_NaN_YYMMDD(G_potCellRunoff, ng, year, month, day_in_month);

            check_array_double_NaN_YYMMDD(G_localGWRunoffIntoRiver, ng, year, month, day_in_month);

        }
        // endif checks
#endif

    for (int n = 0; n < ng; n++) {
        cellArea = geo.areaOfCellByArrayPos(n);

        // added for cell AET calculation (2.1f)
        // set daily variables to zero
        G_dailyLocLakeEvapo[n] = 0.;
        G_dailyGloLakeEvapo[n] = 0.;
        G_dailyResEvapo[n] = 0.;
        G_dailyLocWetlEvapo[n] = 0.;
        G_dailyGloWetlEvapo[n] = 0.;
        G_dailyRiverEvapo[n] = 0.;
        // added for new routing algorithm
        G_riverInflow[n] = 0.;
        G_daily_res_in[n] = 0.;
        G_daily_res_out[n] = 0.;
        G_locwetlextent[n] = 0.;
        G_glowetlextent[n] = 0.;


        // WG22b (option 'use_alloc'): MaxStorage of local lakes and global lakes
        // required to identify the second cell with the maximum 'storageSum'
        if (options.use_alloc > 0) {

            G_locLakeMaxStorage[n] = 0.;
            G_gloLakeMaxStorage[n] = 0.;

            if (G_loc_lake[n] > 0.) {
                // Cell-specific active lake depth - Use array (initialized from cell-specific parameter)
                G_locLakeMaxStorage[n] = ((G_loc_lake[n]) / 100.) * cellArea * G_lakeDepthActive[n];
            }

            if (G_lake_area[n] > 0.) {
                // Cell-specific active lake depth - Use array (initialized from cell-specific parameter)
                G_gloLakeMaxStorage[n] = ((double) G_lake_area[n]) * G_lakeDepthActive[n];
            }
        }
        //	end if (options.use_alloc > 0)
    }
    // end loop over grid cells (initialize)


    double daily_totalUnsatisfiedUse = 0.;

    double remainingUse = 0., dailyUse = 0., totalDesiredUse = 0., totalNeighbourStorage = 0., totalDesiredUseCell = 0.;
    double daily_allocatedUseNextDay = 0., daily_UnsatAllocUseNextDay = 0., dailyActualUse = 0., demandWithoutAllocUse = 0., unsatAllocUse = 0., unsatUse = 0.;
    //introduced in WG2.2b to implement NUs into storage equations:
    double PETgwrRemUseMax = 0., PETgwrRemUse = 0., RiverEvapoRemUse = 0., RivEvapoRemUseMax = 0.;
    //introduced for the new approach in WG2.2b to abstract 50% of NUs from each SWB, if both SWB exist:
    double remainingUseGloLake = 0., remainingUseGloLakeStart = 0., remainingUseRes = 0., remainingUseResStart = 0., remainingUseRiverStart = 0.;


#pragma omp parallel for  \
            private(dailyUse, totalDesiredUse, remainingUse) \
            shared(G_toBeCalculated, options)

    // loop over timeStepsPerDay (currently fixed to 1)
    for (short st = 1; st <= timeStepsPerDay; st++) {

        /* routing order from cell to cell: no parallelization allowed */

        // routingCell: (WG2.2) added for ordered routing scheme
        for (int routingCell = 0; routingCell < ng; routingCell++) {

            // added for new routing order in Wg2.2 // moved some lines up
            int n = (G_routingCell[routingCell] - 1);

            // Efficiency enhancement through local copies instead of often used method getValue
            double P_GWOUTF_C__at__n_routing = calibParam.getValue(P_GWOUTF_C, n);
            double M_EVAREDEX__at__n_routing = calibParam.getValue(M_EVAREDEX, n);

            if (options.use_alloc > 0)
                G_CellPositionRoutOrder[n] = routingCell;


            // initialize values
            remainingUse = 0.;
            dailyUse = 0.;
            totalDesiredUse = 0.;
            totalNeighbourStorage = 0.;
            totalDesiredUseCell = 0.;
            PETgwrRemUseMax = 0.;
            PETgwrRemUse = 0.;
            RiverEvapoRemUse = 0.;
            RivEvapoRemUseMax = 0.;
            remainingUseGloLake = 0.;
            remainingUseGloLakeStart = 0.;
            remainingUseRes = 0.;
            remainingUseResStart = 0.;
            remainingUseRiverStart = 0.;
            dailyActualUse = 0.;
            daily_allocatedUseNextDay = 0.;
            daily_UnsatAllocUseNextDay = 0.;
            demandWithoutAllocUse = 0.;
            unsatAllocUse = 0.;
            unsatUse = 0.;

            daily_totalUnsatisfiedUse = 0.;


            // Only execute for continental grid cells
            if (geo.G_contcell[n]) {

                cellArea = geo.areaOfCellByArrayPos(n);


                // If landAreaFraction == 0., the remaining canopy, snow, and soil storage from the previous timestep
                // becomes surface runoff (see daily.cpp).
                if (G_landAreaFrac[n] == 0.) {
                    dailyLocalSurfaceRunoff = dailyWaterBalance.G_dailyStorageTransfer[n] * cellArea / 1000000. *
                                              G_landAreaFracPrevTimestep[n] / 100.;
                } else {
                    dailyLocalSurfaceRunoff = dailyWaterBalance.G_dailyLocalSurfaceRunoff[n] * cellArea / 1000000. *
                                              G_landAreaFrac[n] / 100.;
                }

                // If GWR under SWB / fractionalRouting is allowed,
                // GW storage is computed before local lakes in humid cells (baseflow directed into local lake)
                // and before the river in (semi-)arid cells (baseflow directed into river)

                // Calculate fswb*Rs and (a-fswb)*Rs for arid cells for Options fractionalRouting and gwr_swb:

                if ((1 == options.aridareaOpt) && (1 == G_aindex[n])) {    // (semi-)arid cells

                    G_fswb_catchment[n] = G_fswbInit[n] *
                                          20.;    // G_fswb_catchment[n] = G_fswb[n] * 5.; // add a catchment area to surface water bodies (specified as x times the size of swb)...

                    if (G_fswb_catchment[n] > 1.) // ... and prevent that G_fswb is greater then 1
                        G_fswb_catchment[n] = 1.;


                    G_localRunoffIntoRiver[n] = (1. - G_fswb_catchment[n]) *
                                                dailyLocalSurfaceRunoff;    // dailyLocalSurfaceRunoff in km3

                    G_localRunoff[n] = G_fswb_catchment[n] * dailyLocalSurfaceRunoff;

                }

                // Calculate baseflow, runoff into local lakes and direct runoff into river for humid cells:

                if ((1 == options.aridareaOpt) && (0 == G_aindex[n])) {    // humid cells only

                    G_groundwaterStoragePrevStep = G_groundwaterStorage[n];

                    // WG22b: ODE (ordinary differential equation) is applied for groundwater storage:
                    // dS/dt = GWR - NAg - k*S  (k*S = groundwater_runoff_km3)

                    // humid cells -> gwr_swb = 0
                    netGWin = (double) dailyWaterBalance.G_dailyGwRecharge[n] * cellArea *
                              ((double) G_landAreaFrac[n] / 100.) / 1000000.;        // mm -> km3

                    // G_dailydailyNUg is adapted in this time step based on the remaining use of last time step
                    if (options.subtract_use > 0) {
                        if (G_dailyRemainingUse[n] != 0) {
                            G_dailydailyNUg[n] = updateNetAbstractionGW(n, month);
                        }
                        netGWin -= G_dailydailyNUg[n];
                    }

                    // Analytical solution of dS/dt = netGWR - k*S
                    // k_g in [1/d]
                    // Cell-specific calibration parameters - Use parameter  // FP

                    G_groundwaterStorage[n] = G_groundwaterStoragePrevStep * exp(-1. * P_GWOUTF_C__at__n_routing) +
                                              (1. / P_GWOUTF_C__at__n_routing) * netGWin *
                                              (1. - exp(-1. * P_GWOUTF_C__at__n_routing));

                    groundwater_runoff_km3 = G_groundwaterStoragePrevStep - G_groundwaterStorage[n] + netGWin;

                    if (groundwater_runoff_km3 <=
                        0.) {    // different differential equation (ODE) applies:
                        // dS/dt = GWR - NAg (without k*S) -> S(t) = S(t-1) + netGWin

                        groundwater_runoff_km3 = 0.;

                        G_groundwaterStorage[n] = G_groundwaterStoragePrevStep + netGWin;
                    }
                    if ((2 == options.grid_store) || (4 == options.grid_store) || (6 == options.grid_store))
                        G_monthlyGwRunoff[n][month] += groundwater_runoff_km3;
                    //if daily 365 output and gwout (same for 31; same for land area fraction)
                    //G_daily365GwRunoff[n][day - 1] = groundwater_runoff_km3;

                    if (1 ==
                        options.fractionalRoutingOpt) {    // only a fraction (depending on amount of surface water bodies) of surface and groundwater runoff should go to the local lakes later, rest should go directly to the river

                        G_fswb_catchment[n] = G_fswbInit[n] *
                                              20.;                  // G_fswb_catchment[n] = G_fswb[n] * 5.; // add a catchment area to surface water bodies (specified as x times the size of swb)...
                        if (G_fswb_catchment[n] > 1.) // ... and prevent that G_fswb is greater then 1
                            G_fswb_catchment[n] = 1.;

                        if (G_aindex[n] ==
                            1) {                                        // in semi-arid/arid areas, all groundwater reaches the river directly
                            G_localGWRunoffIntoRiver[n] = groundwater_runoff_km3;
                        } else // in humid areas, groundwater is also routed through surface water bodies
                            G_localGWRunoffIntoRiver[n] = (1. - G_fswb_catchment[n]) * groundwater_runoff_km3;

                        if (G_aindex[n] == 0) // route groundwater only through surface water bodies in humid regions
                            G_localGWRunoff[n] = G_fswb_catchment[n] * groundwater_runoff_km3;

                        else
                            G_localGWRunoff[n] = 0.; // in arid regions, all groundwater flows directly to the river


                        G_localRunoffIntoRiver[n] = (1. - G_fswb_catchment[n]) *
                                                    dailyLocalSurfaceRunoff;    // dailyLocalSurfaceRunoff in km3

                        G_localRunoff[n] = (G_fswb_catchment[n] * dailyLocalSurfaceRunoff) + G_localGWRunoff[n];

                    } else {    // normal version (without fractional routing)
                        // G_dailyLocalRunoff is devided in G_dailyLocalGWRunoff and G_dailyLocalSurfaceRunoff
                        // values are for land area only
                        if ((G_aindex[n] == 1) && (options.aridareaOpt ==
                                                   1))    // in semi-arid/arid areas, all groundwater reaches the river directly
                            G_localGWRunoffIntoRiver[n] = groundwater_runoff_km3;
                        else
                            G_localGWRunoff[n] = groundwater_runoff_km3;

                        G_localRunoff[n] = dailyLocalSurfaceRunoff + G_localGWRunoff[n];
                    }


                    // to be changed for fractionalrouting option (normal routing without G_localRunffIntoRiver and
                    // G_localRGWRunoffIntoRiver
                    // reintroduced for CFA calculation
                    potCellRunoff = (G_localRunoff[n] + G_localRunoffIntoRiver[n] +
                                     G_localGWRunoffIntoRiver[n]    // is in km3
                                     + ((dailyWaterBalance.G_lakeBalance[n]    // is in mm
                                         / 1000000.0)    // mm -> km3
                                        * (((double) G_lake_area[n] + (double) G_reservoir_area[n])    // km2
                                           + ((geo.areaOfCellByArrayPos(n) / 100.0)    // % -> km2
                                              * (G_loc_lake[n] + G_loc_wetland[n] + G_glo_wetland[n])))));    // [%]

                    G_potCellRunoff[n] += potCellRunoff;

                    if (2 == options.grid_store) {
                        G_monthlyPotCellRunoff[n][month] += potCellRunoff;
                    }

                }    // end of computation of GW storage and baseflow for options.aridareaOpt == 1 in humid cells



                if (options.subtract_use > 0) {

                    if (1 == G_toBeCalculated[n]) {

                        // estimation of dailyUse and totalDesiredUseCell
                        // option for satisfaction of NUs of riparian cells of global lakes and reservoirs
                        if (options.aggrNUsGloLakResOpt == 1) {

                            if (G_reservoir_area[n] > 0. ||
                                G_lake_area[n] > 0.) {    // outflow cell of global lake and/or reservoir

                                dailyUse = G_dailydailyNUsAggregated[n];    // [km3/day]

                                // totalDesiredUseCell (introduced for option 'aggrNUsGloLakResOpt'):
                                // required to compute G_satisfiedUse and G_actualUse of outflow and riparian cells of global lakes/reservoirs
                                // based on the initial G_dailydailyNUs[n] (not aggregated over riparian cells)
                                totalDesiredUseCell = G_dailydailyNUs[n];

                            } else {    // no outflow cell (either riparian cell or other cell)

                                if (G_glwdUnit[n] > 0) {    // riparian cell, but no outflow cell

                                    // G_unsatUseRiparian[n] is set to 0 for outflow cells (see computation of water balance for global lakes and reservoirs)
                                    // For outflow cells,
                                    dailyUse = G_dailydailyNUsAggregated[n] + G_unsatUseRiparian[n];
                                    totalDesiredUseCell = G_dailydailyNUs[n] +
                                                          G_unsatUseRiparian[n]; // only required to compute G_satisfiedUse and G_actualUse for option 'aggrNUsGloLakResOpt' (see above)
                                } else {    // neither outflow nor riparian cell
                                    dailyUse = G_dailydailyNUs[n];
                                    totalDesiredUseCell = G_dailydailyNUs[n];  // only required to redistribute remainingUse over riparian cells for option 'aggrNUsGloLakResOpt' (see above)
                                }
                            }

                        } else {  // if options.aggrNUsGloLakResOpt == 0:

                            // G_dailydailyNUs (surface water abstractions) can be negative (i.e. return flow to surface water bodies)
                            dailyUse = G_dailydailyNUs[n]; // [km3/day]
                            totalDesiredUseCell = G_dailydailyNUs[n]; // only required to redistribute remainingUse over riparian cells for option 'aggrNUsGloLakResOpt' (see above)

                        }
                        // end estimation of dailyUse and totalDesiredUseCell

                        totalDesiredUse = dailyUse;


                        if (options.delayedUseSatisfaction == 1) {
                            totalDesiredUse += G_totalUnsatisfiedUse[n];
                            totalDesiredUseCell += G_totalUnsatisfiedUse[n];// only required to redistribute remainingUse over riparian cells for option 'aggrNUsGloLakResOpt' (see above)
                            // G_PrevTotalUnstatisifiedUse[n] is needed for comparing it to G_totalUnstatisfiedUse to
                            // adapt NAg based on this difference
                            G_PrevTotalUnstatisfiedUse[n] = G_totalUnsatisfiedUse[n];
                            G_totalUnsatisfiedUse[n] = 0;

                        }

                        if (options.use_alloc > 0) {

                            // New approach in WG22b:
                            // If a second cell exists, all unsatisfied use of cell n is added to the remainingUse
                            // in the second cell (G_AllocatedUse[secondCell]) and is satisfied in the same
                            // or the next time step, depending on the routing order.
                            // If there is no second cell and "options.delayedUseSatisfaction == 1":
                            // The remainingUse is satisfied in the next timestep (G_totalUnsatisfiedUse[n]).
                            // Thus, G_AllocatedUse[n] and G_totalUnsatisfiedUse[n] are included in the storage
                            // equations of surface water bodies like dailydailyNUs.
                            // G_UnsatAllocUse[n] is the fraction of allocated use that could not be satisfied in a nbc
                            // in the last timestep.

                            totalDesiredUse += G_UnsatAllocUse[n];
                            totalDesiredUseCell += G_UnsatAllocUse[n];
                            // G_PrevUnsatAllocUse is needed for comparing it to G_UnsatAllocUse[n] to adapt NAg based
                            // on this difference
                            G_PrevUnsatAllocUse[n] = G_UnsatAllocUse[n];
                            G_UnsatAllocUse[n] = 0.;

                            demandWithoutAllocUse = totalDesiredUse;

                            totalDesiredUse += G_AllocatedUse[n];
                            totalDesiredUseCell += G_AllocatedUse[n]; // only required to compute G_satisfiedUse and G_actualUse (see above)

                            // If n is identified as a "second cell" again during this time step,
                            // G_AllocatedUse[n] will be satisfied in the next timestep.
                            // G_dailyAllocatedUse[n] only used as output
                            G_dailyAllocatedUse[n] = G_AllocatedUse[n];
                            G_AllocatedUse[n] = 0.;
                            // daily_allocatedUseNextDay only used as output
                            daily_allocatedUseNextDay = G_daily_allocatedUseNextDay[n];
                            G_daily_allocatedUseNextDay[n] = 0.;
                            // daily_UnsatAllocUseNextDay only used as output
                            daily_UnsatAllocUseNextDay = G_daily_UnsatAllocUseNextDay[n];
                            G_daily_UnsatAllocUseNextDay[n] = 0.;

                        }
                        // end if (options.use_alloc > 0)

                        // 'remainingUse' contains always the amount of water that
                        // has not been satisfied (if "options.use_alloc" > 0 and/or "options.delayedUseSatisfaction" > 0)
                        // At this point, remainingUse can get negative!
                        // G_dailydailyNUs (surface water abstractions) can be negative (i.e. return flow to surface water bodies)
                        remainingUse = totalDesiredUse;
                    }
                    // end if (1 == G_toBeCalculated[n])

                } else {    // if (options.subtract_use == 0)
                    remainingUse = 0.;
                }

                // end "options.subtract_use"

                // added for cell AET calculation (2.2)
                // added for new routing algorithm WG2.2
                // store value for surface water storages in previous time step

                double G_riverStoragePrevStep = 0.;
                double G_locLakeStoragePrevStep = 0.;
                double G_locWetlStoragePrevStep = 0.;
                double G_gloLakeStoragePrevStep = 0.;
                double G_gloResStoragePrevStep = 0.;
                double G_gloWetlStoragePrevStep = 0.;

                if (0 != G_toBeCalculated[n]) {

                    // G_localRunoffPerStep contributes to the inflow of lake/wetland system
                    inflow = G_localRunoff[n];
                    cellArea = geo.areaOfCellByArrayPos(n);

                    // Efficiency enhancement through local copies instead of often used method getValue // FP
                    double P_SWOUTF_C__at__n = calibParam.getValue(P_SWOUTF_C, n);

                    // -------------- Routing through local lake ----------------------------------------------------
                    if (G_loc_lake[n] > 0.) {

                        // assign correct storage value which was calculated in the last timestep or maxStorage at model day one
                        // Note: G_locLakeStorage[n] at this point can differ from the output written at the end of the last time step,
                        // since water abstractions from a second cell (Option "use_alloc") are carried out after writing grids.
                        G_locLakeStoragePrevStep = G_locLakeStorage[n];

                        // maximum storage capacity
                        // Cell-specific active lake depth - Use array (initialized from cell-specific parameter)  // FP
                        maxStorage = ((G_loc_lake[n]) / 100.) * cellArea * G_lakeDepthActive[n]; // [km3]
                        G_locLakeMaxStorage[n] = maxStorage; // G_locLakeMaxStorage is used when NUs is abstracted from locLakeStorage at the end of the routing loop (in cell n and in second cell)

                        // This factor was added to reduce PET as a function of actual lake storage (2.1f).
                        // Without reduction PET would lead to a continuous decline of lake level in some cases ((semi)arid regions)
                        // use reduction factor from previous time step (or the initial value)
                        locLakeAreaReductionFactor = G_locLakeAreaReductionFactor[n];


                        //reintroduced CFA

                        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                        // CFA approach for surface water bodies from version 2.1f
                        //
                        // Water balance of lake / reservoir:
                        // dS/dt = CFA*(P - PET) + Qin - NAs - gwr_swb
                        //
                        // (P - PET) is considered as the runoff from the respective SWB. Therefore, only (P - PET) is multiplied by CFA
                        // and only P and PET are taken into account to compute PET_corrected.
                        // (different from soil water balance, where R equals (P - AET - dS/dt), and dS/dt is taken into account to compute PET_corrected.
                        //
                        // SWB: CFA affects R, AET and dS/dt.
                        // Land water balance (daily.cpp): CFA affects only R and AET.
                        // (as in version 2.1f, but code implementation changed in 2.2b)
                        //
                        // Approach to compute PET_corrected:
                        //
                        //	equation 1) P - PET = R
                        // 	equation 2) P - PET_corr = R * CFA
                        // 						   R = (P - PET_corr)/CFA
                        //
                        // 	equation 2) in equation 1):
                        // 	P - PET = (P - PET_corr)/CFA			| *CFA
                        // 	(P - PET) * CFA = P - PET_corr
                        //
                        // 	CFA*P - CFA * PET = P - PET_corr
                        // 	PET_corr = P - CFA*P + CFA * PET
                        //
                        // Final equation in routing.cpp:
                        //
                        // 	PET_corr = (1 - CFA) * P + CFA * PET
                        //
                        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
                        //
                        // compare daily.cpp, line 1217:
                        // 	AET_corr = dS/dt*(CFA - 1) - P*(CFA - 1) + CFA * AET_uncorr
                        //
                        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

                        G_dailyLocLakeEvapo[n] =
                                       ((1.0 - dailyWaterBalance.G_cellCorrFact[n])                        // [mm]
                                        * dailyWaterBalance.G_openWaterPrec[n] * locLakeAreaReductionFactor)
                                       + (dailyWaterBalance.G_cellCorrFact[n]
                                          * (dailyWaterBalance.G_openWaterPET[n] * locLakeAreaReductionFactor));


                        totalInflow = inflow + (dailyWaterBalance.G_openWaterPrec[n] * locLakeAreaReductionFactor)
                                               * (cellArea / 1000000.) * (G_loc_lake[n] / 100.);    //[mm]->[km3]


                        // gwr below surface water bodies
                        if ((1 == options.aridareaOpt) && (1 == G_aindex[n])) {
                            gwr_loclak[n] = Kswbgw * locLakeAreaReductionFactor * G_loc_lake[n] / 100. /
                                            (geo.G_contfreq[n] /
                                             100.); // gwr below local lakes in mm/d with respect to continental area

                        } else
                            gwr_loclak[n] = 0.;    // -> to avoid repeating the following equations for options.aridareaOpt==0

                        PETgwr = G_dailyLocLakeEvapo[n] * (cellArea / 1000000.) *
                                 (G_loc_lake[n] / 100.)                            //[mm]->[km3]
                                 + gwr_loclak[n] * cellArea * (geo.G_contfreq[n] / 100.) / 1000000.;    //[mm -> km3]

                        // Maximum amount of (openwaterPET + gwr_loclak) to ensure that G_locLakeStorage is always larger than -maxStorage:
                        PETgwrMax = G_locLakeStoragePrevStep + maxStorage + totalInflow;

                        if (PETgwr > PETgwrMax) {    // reduction of PET and gwr_locLak required

                            G_locLakeStorage[n] = (-1.) * maxStorage;
                            gwr_loclak[n] *= PETgwrMax / PETgwr;    //[mm]
                            G_dailyLocLakeEvapo[n] *= PETgwrMax / PETgwr;    //[mm]

                        } else { // no reduction of PET and gwr_loclak required

                            G_locLakeStorage[n] = G_locLakeStoragePrevStep + totalInflow - PETgwr;
                        }

                        // WG22b: outflow is a function of storage from previous time step (outflow = f(St-1)). Until WG22: outflow = f(S(t))
                        if (G_locLakeStoragePrevStep > 0.) {
                            // Cell-specific calibration parameters - Use parameter  // FP
                            // loc_storageFactor (now: only one factor for surface water P_SWOUTF_C) in [1/d]
                            outflow = P_SWOUTF_C__at__n * G_locLakeStoragePrevStep
                                      * pow((G_locLakeStoragePrevStep / maxStorage), lakeOutflowExp);

                            if (G_locLakeStorage[n] <=
                                0.) {    // Even if S(t-1) is positive, G_locLakeStorage[n] can become negative, if PET+gwr > totalInflow.
                                outflow = 0;                    // In this case, the outflow is set to zero.
                            } else {
                                if (outflow > G_locLakeStorage[n])
                                    outflow = G_locLakeStorage[n];
                            }

                        } else // if G_locLakeStoragePrevStep <= 0, no outflow should occur
                            outflow = 0.;

                        // New approach in WG22b: In inland sinks, outflow from SWB flows into the river.
                        // The computed transported volume /(= outflow of the river) represents the water resource of this cell and thus should not be zero (although
                        // it cannot be assigned to a downstream cell). This value can be interpreted as either additional AET or river availability of the inland sink.
                        /*if (G_LDD[n] == -1) {
                            if (((G_glo_lake[n] + G_glo_wetland[n]) == 0) && (G_loc_wetland[n] == 0)) {
                                outflow = 0.;
                            }
                        }*/

                        G_locLakeStorage[n] -= outflow;

                        // reduce G_locLakeStorage[n] to maximum storage capacity
                        if (G_locLakeStorage[n] > maxStorage) {
                            outflow += (G_locLakeStorage[n] - maxStorage);
                            G_locLakeStorage[n] = maxStorage;
                        }

                        // If (PETgwr > PETgwrMax), PET is reduced by PETgwrRedFactor. If a reduction is not required, PETgwrRedFactor is set to 1.
                        // HMS 2014-09-26 if PETgwr is 0 (because gwr_swb is 0 because reductionfactor is 0 because storage is empty),
                        // division through 0 --> has to be catched, because otherwise G_CELL_RUNOFF_TOTAL and G_CELL_AET gets nan.
                        if (PETgwr == 0.)
                            PETgwrRedFactor = 1.;
                        else
                            PETgwrRedFactor = PETgwrMax / PETgwr;

                        if (PETgwrMax > PETgwr) {// -> reduction of PET not required
                            PETgwrRedFactor = 1.;
                        }

                        if (2 == options.grid_store)
                            // sum up the negative fraction of the potential cell runoff which has not been realized
                            // due to reduction of potential evaporation
                            G_monthlyPotCellRunoffDeficit[n][month]
                                           += (0.0 - dailyWaterBalance.G_openWaterPET[n] * PETgwrRedFactor
                                                     * (1.0 - locLakeAreaReductionFactor)
                                                     * dailyWaterBalance.G_cellCorrFact[n]
                                                     * (cellArea / 1000000.0)
                                                     * (G_loc_lake[n] / 100.0));


                        inflow = outflow;

                        // Area reduction factor is now calculated for the next day.
                        // Cell-specific calibration parameters - Apply multiplier
                        G_locLakeAreaReductionFactor[n] = 1. -
                                                          pow(fabs(G_locLakeStorage[n] - maxStorage)
                                                              / (2. * maxStorage),
                                                              (M_EVAREDEX__at__n_routing * evapoReductionExp));

                        if (G_locLakeAreaReductionFactor[n] <
                            0.)
                            G_locLakeAreaReductionFactor[n] = 0.;

                        if (G_locLakeAreaReductionFactor[n] > 1.)  // (could occur due to numerical inaccuracies)
                            G_locLakeAreaReductionFactor[n] = 1.;

                    }
                    // outflow of local lakes is inflow to local wetlands


                    // transport through local wetlands
                    if (G_loc_wetland[n] > 0.) {

                        // assign correct storage value which was calculated in the last timestep or maxStorage at model day one
                        G_locWetlStoragePrevStep = G_locWetlStorage[n];

                        // maximum storage capacity
                        // Cell-specific active wetland depth - Use array (initialized from cell-specific parameter)
                        maxStorage = ((G_loc_wetland[n]) / 100.) * cellArea * G_wetlDepthActive[n];

                        // use reduction factor from previous time step (or the initial value)
                        locWetlAreaReductionFactor = G_locWetlAreaReductionFactor[n];

                        // G_dailyLocWetlEvapo[n]: see local lakes for explanation

                        G_dailyLocWetlEvapo[n] =
                                       ((1.0 - dailyWaterBalance.G_cellCorrFact[n])                        // [mm]
                                        * dailyWaterBalance.G_openWaterPrec[n] * locWetlAreaReductionFactor)
                                       + (dailyWaterBalance.G_cellCorrFact[n]
                                          * (dailyWaterBalance.G_openWaterPET[n] * locWetlAreaReductionFactor));

                        totalInflow = inflow + (dailyWaterBalance.G_openWaterPrec[n] * locWetlAreaReductionFactor
                                                * (cellArea / 1000000.) * (G_loc_wetland[n] / 100.));



                        G_locwetlextent[n] = locWetlAreaReductionFactor * cellArea * (G_loc_wetland[n] / 100.);
                        // gwr below surface water bodies
                        if ((1 == options.aridareaOpt) && (1 == G_aindex[n])) {
                            gwr_locwet[n] = Kswbgw * locWetlAreaReductionFactor * G_loc_wetland[n] / 100. /
                                            (geo.G_contfreq[n] /
                                             100.); // gwr below local wetlands in mm/d with respect to continental area

                        } else
                            gwr_locwet[n] = 0.; // -> to avoid repeating the following equations for options.aridareaOpt==0

                        PETgwr = G_dailyLocWetlEvapo[n] * (cellArea / 1000000.) *
                                 (G_loc_wetland[n] / 100.)                        //[mm -> km3]
                                 + gwr_locwet[n] * cellArea * (geo.G_contfreq[n] / 100.) / 1000000.;    //[mm -> km3]

                        // Maximum amount of (openwaterPET + gwr_locwet) to ensure that G_locWetlStorage does not fall below zero:
                        PETgwrMax = G_locWetlStoragePrevStep + totalInflow;

                        if (PETgwr > PETgwrMax) {

                            G_locWetlStorage[n] = 0.;
                            gwr_locwet[n] *= PETgwrMax / PETgwr;    // [mm]
                            G_dailyLocWetlEvapo[n] *= PETgwrMax / PETgwr;    // [mm]

                        } else {
                            G_locWetlStorage[n] = G_locWetlStoragePrevStep + totalInflow - PETgwr;
                        }

                        if (G_locWetlStorage[n] > 0.) {
                            // Cell-specific calibration parameters - Use parameter  // FP
                            // loc_storageFactor (now: only one factor for surface water P_SWOUTF_C) in [1/d]
                            outflow = P_SWOUTF_C__at__n * G_locWetlStorage[n]
                                      * pow((G_locWetlStorage[n] / maxStorage), wetlOutflowExp);

                            if (outflow > G_locWetlStorage[n]) {    // this could happen if S(t-1) >> S(t)
                                outflow = G_locWetlStorage[n];
                            }

                        } else
                            outflow = 0.;

                        // New approach in WG22b: In inland sinks, outflow from SWB flows into the river.
                        // The computed transported volume /(= outflow of the river) represents the water resource of this cell and thus should not be zero (although
                        // it cannot be assigned to a downstream cell). This value can be interpreted as either additional AET or river availability of the inland sink.

                        // substract outflow from storage
                        G_locWetlStorage[n] -= outflow;

                        // reduce G_locWetlStorage to maximum storage capacity

                        if (G_locWetlStorage[n] > maxStorage) {
                            outflow += (G_locWetlStorage[n] - maxStorage);
                            G_locWetlStorage[n] = maxStorage;
                        }

                        // If (PETgwr > PETgwrMax), PET is reduced by PETgwrRedFactor. If a reduction is not required, PETgwrRedFactor is set to 1.
                        // if PETgwr is 0 (because gwr_swb is 0 because reductionfactor is 0 because storage is empty),
                        // division through 0 --> has to be catched, because otherwise G_CELL_RUNOFF_TOTAL and G_CELL_AET gets nan.
                        if (PETgwr ==
                            0.)
                            PETgwrRedFactor = 1.;
                        else
                            PETgwrRedFactor = PETgwrMax / PETgwr;

                        if (PETgwrMax > PETgwr)    // -> reduction of PET not required
                            PETgwrRedFactor = 1.;

                        // reintroduced CFA
                        if (2 == options.grid_store)
                            // sum up the negative fraction of the potential cell runoff which has not been realized
                            // due to reduction of potential evaporation
                            G_monthlyPotCellRunoffDeficit[n][month]
                                           += (0.0 - dailyWaterBalance.G_openWaterPET[n] * PETgwrRedFactor
                                                     * (1.0 - locWetlAreaReductionFactor)
                                                     * dailyWaterBalance.G_cellCorrFact[n]
                                                     * (cellArea / 1000000.0)
                                                     * (G_loc_wetland[n] / 100.0));

                        inflow = outflow;

                        // Area reduction factor is now calculated for the next day.
                        // Cell-specific calibration parameters - Apply multiplier
                        G_locWetlAreaReductionFactor[n] = 1. -
                                                          pow(fabs(G_locWetlStorage[n] - maxStorage)
                                                              / (maxStorage),
                                                              (M_EVAREDEX__at__n_routing * evapoReductionExp));


                        if (G_locWetlAreaReductionFactor[n] < 0.)    // (could occur due to numerical inaccuracies)
                            G_locWetlAreaReductionFactor[n] = 0.;

                        if (G_locWetlAreaReductionFactor[n] > 1.)  // (could occur due to numerical inaccuracies)
                            G_locWetlAreaReductionFactor[n] = 1.;

                    }
                    // end if (G_loc_wetland[n] > 0.)

                    // water that comes out of the system of local lakes/wetlands
                    // is inflow into the river and is being routed through global lakes, reservoirs and global wetlands


                    inflow += G_riverInflow[n];
                    inflowUpstream = G_riverInflow[n];



                    // transport through global lakes
                    // input is defined as the inflow from the local wetlands in addition to the upstream flow

                    // new reservoir algorithm
                    if (G_lake_area[n] > 0.) //cell is a natural lake, or the reservoir type is unknown
                    {

                        //assign correct storage value which was calculated in the last timestep or maxStorage at model day one
                        G_gloLakeStoragePrevStep = G_gloLakeStorage[n];

                        // maximum storage capacity
                        // Cell-specific active lake depth - Use array (initialized from cell-specific parameter)  // FP
                        maxStorage = ((double) G_lake_area[n]) * G_lakeDepthActive[n];
                        G_gloLakeMaxStorage[n] = maxStorage; // CR: Introduced in WG22b: G_gloLakeMaxStorage is used when NUs is abstracted from G_gloLakeMaxStorage of a second (neighboring) cell.
                        //  maxStG_locWetlStorage[n]orage = (((double) G_glo_lake[n]) / 100.) * cellArea * lakeDepth;	// [km3]

                        // This factor was added to reduce PET as a function of actual lake storage (2.1f).
                        // Without reduction PET would lead to a continuous decline of lake level in some cases ((semi)arid regions)
                        //HMS use reduction factor from previous time step (or the initial value)
                        gloLakeEvapoReductionFactor = G_gloLakeEvapoReductionFactor[n];

                        // G_dailyGloLakeEvapo[n]: see local lakes for explanation

                        G_dailyGloLakeEvapo[n] =
                                       ((1.0 - dailyWaterBalance.G_cellCorrFact[n])                            // [mm]
                                        * dailyWaterBalance.G_openWaterPrec[n])
                                       + (dailyWaterBalance.G_cellCorrFact[n]
                                          * dailyWaterBalance.G_openWaterPET[n] * gloLakeEvapoReductionFactor);

                        // reintroduced CFA in calculation (!no distinguishing for aridareaOpt!)
                        totalInflow = inflow + (dailyWaterBalance.G_openWaterPrec[n]
                                                * ((double) G_lake_area[n] / 1000000.));

                        G_GloLakePrecip[n] = dailyWaterBalance.G_openWaterPrec[n] * ((double) G_lake_area[n]) /
                                             1000000.; // km3

                        // gwr below surface water bodies
                        if ((1 == options.aridareaOpt) && (1 == G_aindex[n])) {

                            gwr_glolak[n] = Kswbgw * gloLakeEvapoReductionFactor * ((double) G_lake_area[n] /
                                                                                    (cellArea * (geo.G_contfreq[n] /
                                                                                                 100.))); // gwr below global lakes in mm/d

                        } else {
                            gwr_glolak[n] = 0.;    // to avoid repeating the following equations for options.aridareaOpt==0
                        }


                        // remainingUse is evenly distributed between global reservoirs and global lakes (if existing) within a grid cell
                        // If remainingUse is negative (i.e. return flow to surface water bodies, from negative surface water abstractions G_dailydailyNUs)
                        // then it is added to the storage
                        if (G_reservoir_area[n] >
                            0.) {    // If the cell also contains a reservoir, 50% of NAs is withdrawn from each of the SWB. If "remainingUseGloLake" cannot be fully satisfied, it is added to the 50% NAs of the reservoir.
                            remainingUseGloLake = 0.5 * remainingUse;
                        } else {
                            remainingUseGloLake = remainingUse;
                        }

                        // Required in 22b to compute actual use for Option 28 "aggrNUsGloLakResOpt"
                        remainingUseGloLakeStart = remainingUseGloLake;


                        PETgwrRemUse = G_dailyGloLakeEvapo[n] * ((double) G_lake_area[n] /
                                                                 1000000.)                                    //[mm]->[km3]
                                       + gwr_glolak[n] * cellArea * (geo.G_contfreq[n] / 100.) / 1000000.    //[mm]->[km3]
                                       + remainingUseGloLake;

                        // Maximum amount of water available to satisfy PET, gwr and NAs,
                        // until St = -Smax (outflow = 0 if storage <= 0).
                        PETgwrRemUseMax = totalInflow + maxStorage + G_gloLakeStoragePrevStep;

                        if (PETgwrRemUse >
                            PETgwrRemUseMax) {    // if (PETgwrRemUse >= PETgwrRemUseMax): PETgwrRemUse could be zero -> division through zero.

                            G_gloLakeStorage[n] = (-1.) * maxStorage;
                            outflow = 0.;
                            G_dailyGloLakeEvapo[n] *=
                                           PETgwrRemUseMax / PETgwrRemUse;                                // [mm]
                            gwr_glolak[n] *= PETgwrRemUseMax / PETgwrRemUse;                                // [mm]

                            if (remainingUseGloLake > 0.) {

                                remainingUseGloLake -= remainingUseGloLake * PETgwrRemUseMax / PETgwrRemUse;    //[km3]

                            } else {
                                remainingUseGloLake = 0.;
                            }


                        }
                            // endif (PETgwrRemUse > PETgwrRemUseMax)
                        else {    //if (PETgwrRemUse <= PETgwrRemUseMax):

                            // This ODE (ordinary differential equation) applies for S >= 0 only:
                            // dS/dt = totalInflow - PET - gwr - NAs - k*S
                            // If the outflow becomes zero, S is computed again using a different ODE.

                            // Cell-specific calibration parameters - Use parameter
                            // glo_storageFactor (now: only one factor for surface water P_SWOUTF_C) in [1/d]
                            G_gloLakeStorage[n] = G_gloLakeStoragePrevStep * exp(-1. * P_SWOUTF_C__at__n)
                                                  + (1. / P_SWOUTF_C__at__n) * (totalInflow - PETgwrRemUse)
                                                  * (1. - exp(-1. * P_SWOUTF_C__at__n));

                            outflow = totalInflow + G_gloLakeStoragePrevStep - G_gloLakeStorage[n] - PETgwrRemUse;

                            // reduce G_gloLakeStorage[n] to maximum storage capacity
                            if (G_gloLakeStorage[n] > maxStorage) {
                                outflow += (G_gloLakeStorage[n] - maxStorage);
                                G_gloLakeStorage[n] = maxStorage;
                            }

                            if (outflow < 0.) {

                                // Different ODE (ordinary differential equation) applies for S < 0:
                                // dS/dt = totalInflow - PET - gwr - NAs  (without -k*S)
                                // S(t) = S(t-1) + inflow - PET - gwr - NAs

                                outflow = 0.;

                                // Storage equation without outflow:
                                G_gloLakeStorage[n] = G_gloLakeStoragePrevStep + totalInflow - PETgwrRemUse;


                            }
                            // endif (outflow < 0.)

                            remainingUseGloLake = 0.;
                            // Reduction of gwr_glolak[n] and evaporation not required in this case.


                        }
                        // endelse (PETgwrRemUse <= PETgwrRemUseMax)

                        // If (PETgwrRemUse > PETgwrRemUseMax), PET is reduced by PETgwrRemUseRedFactor.
                        // If a reduction is not required, PETgwrRemUseRedFactor is set to 1.
                        if (PETgwrRemUse == 0.) {
                            PETgwrRemUseRedFactor = 1.;
                        } else {
                            PETgwrRemUseRedFactor = PETgwrRemUseMax / PETgwrRemUse;
                        }
                        // endelse

                        if (PETgwrRemUseMax > PETgwrRemUse) {    // -> reduction of PET not required
                            PETgwrRemUseRedFactor = 1.;
                        }
                        // endif


                        if (2 == options.grid_store) {
                            // sum up the negative fraction of the potential cell runoff which has not been realized
                            // due to reduction of potential evaporation
                            G_monthlyPotCellRunoffDeficit[n][month]
                                           += (0.0 -
                                               (dailyWaterBalance.G_openWaterPET[n] * PETgwrRemUseRedFactor
                                                * (1.0 - gloLakeEvapoReductionFactor)
                                                * dailyWaterBalance.G_cellCorrFact[n]
                                                * G_lake_area[n]
                                                / 1000000.0));
                        }

                        inflow = outflow;

                        // introduced for 22b, Option 28 "aggrNUsGloLakResOpt"
                        dailyActualUse = remainingUseGloLakeStart - remainingUseGloLake;

                        // 'remainingUseGloLake' is redistributed over the riparian cells (and the outflow cell if the
                        // initial 'dailydailyNUs' was larger than zero) (see documentation WG22b)
                        // If a reservoir exists, 'remainingUseGloLake' is added to 'remainingUseRes', which is
                        // redistributed over the riparian cells at the end of the reservoir loop.
                        if (0 == G_reservoir_area[n]) {

                            if (options.aggrNUsGloLakResOpt == 1) {

                                G_remainingUse[n] = remainingUse;
                                G_remainingUseGloLake[n] = remainingUseGloLake;
                                G_remainingUseGloLakeRedistr[n] = remainingUseGloLake;
                                G_totalDesiredUseCell[n] = totalDesiredUseCell;

                                for (int i = 0; i < ng; i++) {

                                    if (G_glwdUnit[n] == G_glwdUnit[i]) {

                                        if (G_remainingUseGloLakeRedistr[n] > 0.) {

                                            if (G_totalDesiredUseCell[n] <=
                                                0.) {        // different computation of G_unsatUseRiparian[i], if G_totalDesiredUseCell[n] of the OUTFLOW CELL is negative - see WG22b documentation

                                                // (G_dailydailyNUsAggregated[n] - G_dailydailyNUs[n]) always larger than zero in this loop

                                                if (G_dailydailyNUs[i] <= 0.)
                                                    G_unsatUseRiparian[i] = 0.;

                                                else
                                                    G_unsatUseRiparian[i] = (G_dailydailyNUs[i] /
                                                                             (G_dailydailyNUsAggregated[n] -
                                                                              G_dailydailyNUs[n])) *
                                                                            G_remainingUseGloLakeRedistr[n];
                                                // G_unsatUseRiparian is also computed for the outflow cell (when i == n), but it is only used for riparian cells, not for outflow cells.

                                                G_remainingUseGloLake[n] = 0.;
                                                G_unsatUseRiparian[n] = 0.;

                                            } else {    // G_totalDesiredUseCell[n] in outflow cell > 0 (G_totalDesiredUseCell[n] includes G_dailydailyNUs and (depending on selected options) allocated use from nbc and totalUnsatisfiedUse[n] (from previous day)

                                                //	At this point, G_remainingUse[n] is always > 0. (within loop "G_remainingUseGloLakeRedistr[n] > 0.")
                                                G_unsatUseRiparian[i] = G_dailydailyNUs[i] / G_remainingUse[n] *
                                                                        G_remainingUseGloLakeRedistr[n];

                                                if (G_dailydailyNUs[i] < 0.)
                                                    G_unsatUseRiparian[i] = 0.;

                                                G_remainingUseGloLake[n] *= (G_totalDesiredUseCell[n] /
                                                                             G_remainingUse[n]);    // G_remainingUseGloLake[n] to be satisfied from river storage

                                                G_unsatUseRiparian[n] = 0.;
                                                // CR: Required since G_unsatUseRiparian[n] is added to 'remainingUse' at the beginning of each time step for all cells with GLWD units.
                                            }

                                        } else {    // if G_remainingUseGloLake[n] == 0.:
                                            G_unsatUseRiparian[i] = 0.;
                                            // required? -> already covered by G_unsatUseRiparian[i] = 0?
                                            G_unsatUseRiparian[n] = 0.;
                                            G_remainingUseGloLake[n] = 0.;
                                        }

                                    }    // end of loop over 'i' with matching GLWD unit
                                }
                                // end of loop over grid cells 'i'

                                remainingUse = G_remainingUseGloLake[n];

                            }   // end if (options.aggrNUsGloLakResOpt == 1)
                            else { // Option 'aggrNUsGloLakResOpt' == 0 and no reservoir in this cell
                                remainingUse = remainingUseGloLake;
                            }
                            // end else

                        }  // end if (0 == G_reservoir_area[n])
                        // if (G_reservoir_area[n] > 0.): 'remainingUseGloLake' is added to 'remainingUseRes' at the beginning of the reservoir loop


                        //Evaporation reduction factor is now calculated for the next day.
                        // Cell-specific calibration parameters - Apply multiplier  // FP
                        G_gloLakeEvapoReductionFactor[n] = 1. -
                                                           pow(fabs(G_gloLakeStorage[n] - maxStorage)
                                                               / (2. * maxStorage),
                                                               (M_EVAREDEX__at__n_routing * evapoReductionExp));


                        if (G_gloLakeEvapoReductionFactor[n] < 0.)    // (could occur due to numerical inaccuracies)
                            G_gloLakeEvapoReductionFactor[n] = 0.;

                        if (G_gloLakeEvapoReductionFactor[n] > 1.)  // (could occur due to numerical inaccuracies)
                            G_gloLakeEvapoReductionFactor[n] = 1.;

                    }
                    //end  if (G_lake_area[n] > 0)


                    // This variable should be defined when calculating total inflow to reservoirs.
                    //res_inflow = inflow - G_riverInflow[n];


                    // transport through reservoir if (options.resOpt == 1)
                    if (G_reservoir_area[n] > 0.) { //cell is a reservoir cell with a well known type

                        double c_ratio;
                        double prov_rel; //provisional release
                        double release;
                        double dailyUse;
                        double monthlyUse;

                        // define storage capacity to mean annual outflow ratio (c_ratio)
                        // G_mean_outflow:  m3/s,  G_stor_cap:  km3
                        // 31536000=(365 * 24 * 60 * 60)
                        c_ratio = G_stor_cap[n] / (G_mean_outflow[n] * 31536000. / 1000000000.);
                        // Reservoir operation start years
                        // maxStorage = G_reservoir_area[n] * 0.01; //10 m depth for active storage volume
                        // not calculated here because for global lake handling of regulated lakes before operation year we must seperate calculation.
                        maxStorage = G_stor_cap[n];

                        G_gloResMaxStorage[n] = maxStorage;    // required for Option "use_alloc" (water abstraction from second cell)

                        G_daily_res_in[n] = inflow; // [km3/day]

                        //assign correct storage value which was calculated in the last timestep or maxStorage at model day one
                        G_gloResStoragePrevStep = G_gloResStorage[n];

                        //same as gloLakeEvapo
                        // use reduction factor from previous time step (or the initial value)
                        gloResEvapoReductionFactor = G_gloResEvapoReductionFactor[n];

                        // G_dailyResEvapo[n]: see local lakes for explanation

                        G_dailyResEvapo[n] =
                                       ((1.0 - dailyWaterBalance.G_cellCorrFact[n])                            // [mm]
                                        * dailyWaterBalance.G_openWaterPrec[n])
                                       + (dailyWaterBalance.G_cellCorrFact[n]
                                          * (dailyWaterBalance.G_openWaterPET[n] * gloResEvapoReductionFactor));

                        //inflow from upstream PLUS lake water balance
                        // reintroduced CFA in calculation (!no distinguishing for aridareaOpt!)
                        totalInflow = inflow + (dailyWaterBalance.G_openWaterPrec[n]
                                                * ((double) G_reservoir_area[n] / 1000000.));

                        G_ReservoirPrecip[n] = dailyWaterBalance.G_openWaterPrec[n] * ((double) G_reservoir_area[n]) /
                                               1000000.; // km3

                        // gwr below surface water bodies
                        if ((1 == options.aridareaOpt) && (1 == G_aindex[n])) {
                            gwr_res[n] = Kswbgw * gloResEvapoReductionFactor * ((double) G_reservoir_area[n] /
                                                                                (cellArea * (geo.G_contfreq[n] /
                                                                                             100.))); // gwr below reservoirs in mm/d per continental area

                        } else {
                            gwr_res[n] = 0.;    // -> to avoid repeating the following equations for options.aridareaOpt==0
                        }

                        // remainingUse is evenly distributed between global reservoirs and global lakes (if existing) within a grid cell
                        // If remainingUse is negative (i.e. return flow to surface water bodies, from negative surface water abstractions G_dailydailyNUs)
                        //      then it is added to the storage
                        if (G_lake_area[n] >
                            0.) // If the cell also contains a global lake, 50% of NAs has already been subtracted from global lake storage.
                            remainingUseRes = 0.5 * remainingUse + remainingUseGloLake;

                        else
                            remainingUseRes = remainingUse;

                        // Required in 22b to compute actual use for Option 28 "aggrNUsGloLakResOpt"
                        remainingUseResStart = remainingUseRes;

                        // PET and gwr can occur until St = 0, while NAs is only withdrawn until St = 0.1*stor_capacity.
                        // Therefore, in case of reservoirs, the three variables cannot be equally satisfied.
                        // First, PET and gwr are satisfied. Then, if St > 0.1*stor_capacity, 'remainingUseRes' is subtracted from storage.

                        PETgwr = G_dailyResEvapo[n] *
                                 ((double) G_reservoir_area[n] / 1000000.) //[mm]->[km3]
                                 + gwr_res[n] * cellArea * (geo.G_contfreq[n] / 100.) / 1000000.;    //[mm]->[km3]

                        //WG22b: Water abstractions from reservoirs are reduced to 'remainingUseMax' to ensure that G_gloResStorage does not fall below 10% of the storage capacity.
                        PETgwrMax = G_gloResStoragePrevStep + totalInflow;

                        if (PETgwr > PETgwrMax) {

                            G_gloResStorage[n] = G_gloResStoragePrevStep + totalInflow - PETgwrMax;
                            gwr_res[n] *= PETgwrMax / PETgwr;
                            G_dailyResEvapo[n] *= PETgwrMax / PETgwr;

                        } else {
                            G_gloResStorage[n] = G_gloResStoragePrevStep + totalInflow - PETgwr;
                        }

                        // Treatment of NAs (net abstraction of surface water): add or substract from available storage
                        if (options.subtract_use > 0) {

                            // Negative values of NAs (return flow to surface water) are added (signs "- plus - equals +") to the storage:
                            if (remainingUseRes < 0.) {
                                G_gloResStorage[n] -= remainingUseRes;
                                remainingUseRes = 0.;
                            }
                                // Subtract remainingUseRes only if G_gloResStorage[n] > (0.1 * G_stor_cap[n]):
                            else {    // remainingUseRes >=0.
                                if (G_gloResStorage[n] > (0.1 * G_stor_cap[n])) {

                                    if (remainingUseRes < (G_gloResStorage[n] - (G_stor_cap[n] * 0.1))) {
                                        G_gloResStorage[n] -= remainingUseRes;
                                        remainingUseRes = 0.;
                                    }
                                        // When remainingUseRes is larger or equal to available storage
                                        // then subtract the available storage amount until 10% stor_cap is reached
                                    else {
                                        remainingUseRes -= (G_gloResStorage[n] - (G_stor_cap[n] * 0.1));
                                        G_gloResStorage[n] = G_stor_cap[n] * 0.1;
                                    }

                                }
                                // endif (gloResStorage > 10% stor_cap)
                                // else (gloResStorage <= 10% stor_cap): no abstraction allowed
                            }
                            // endelse remainingUseRes >=0.

                        }
                        // endif (options.subtract_use > 0)


                        // Different from local lakes and local wetlands, the outflow (release) of reservoirs is
                        // computed based on S(t) and not S(t-1).

                        // Calculate release (outflow) of reservoir:
                        // set rules at the beginning of the operational year! (only once!)
                        if (month == G_start_month[n] - 1) {
                            if (day == first_day_in_month[month]) {
                                // if (st == 1) {	// st always 1; not required anymore
                                    //1. calculate release coefficient for the actual year:
                                    //reduce release coefficent in this year to refill storage volume in reservoir
                                    if (G_gloResStorage[n] < (G_stor_cap[n] * 0.1))
                                        K_release[n] = 0.1;
                                    else
                                        K_release[n] = G_gloResStorage[n] / (maxStorage * 0.85);

                                //2 calculate annual release
                                annual_release[n] = K_release[n] * G_mean_outflow[n];
                                //} // st always 1; not required anymore
                            }
                        }

                        if ((G_res_type[n] + 0) == 1) {// (irrigation reservoir)
                            // calculate monthly demand of downstream area
                            dailyUse = 0.0;

                            // read in net surface water use
                            dailyUse = G_dailydailyNUs[n];// [km3/day]

                            short i = 0;
                            int downstreamCell = G_downstreamCell[n];
                            // max reservoir_dsc cells downstream or to the basin outlet or area up to next reservoir
                            while (i < reservoir_dsc && downstreamCell > 0 && G_reservoir_area[downstreamCell] <= 0) {
                                dailyUse += G_dailydailyNUs[downstreamCell - 1] * G_alloc_coeff[n][i++];

                                // next downstream cell
                                downstreamCell = G_downstreamCell[downstreamCell - 1];
                            }

                            monthlyUse = dailyUse * numberOfDaysInMonth[month]; // [km3/month]
                            monthlyUse = monthlyUse * 1000000000. / (numberOfDaysInMonth[month] * 86400.); //[m3/s]

                            if (G_mean_demand[n] >= 0.5 * G_mean_outflow[n])
                                //provisional monthly release = i_mean /2. * (1. + monthly_demand_sum / d_mean_annual)
                                prov_rel = G_mean_outflow[n] / 2. * (1. + monthlyUse / G_mean_demand[n]); //[m3/s]
                            else
                                //provisional monthly release = i_mean + monthly_demand_sum - d_mean)
                                prov_rel = G_mean_outflow[n] + monthlyUse - G_mean_demand[n]; //[m3/s]
                        } else if ((G_res_type[n] + 0) == 2) {//(non-irrigation)
                            prov_rel = G_mean_outflow[n]; //[m3/s]
                        } else
                            cerr << "unknown reservoir type in gcrcNumber " << n + 1 << '\n';

                        // calculate release [m3/s]
                        if (c_ratio >= 0.5) {
                            release = K_release[n] * prov_rel;
                        } else { //monthly_release = (c_ratio/0,5)^2 * provisional monthly release + (1-(c_ratio/0,5)^2) * monthly_inflow (just daily inflow?)
                            release = ((4. * c_ratio * c_ratio) * K_release[n] * prov_rel)
                                      //+ ((1.0 - ((4.*c_ratio*c_ratio))) * inflow * 1000000000./(2.*3600.));
                                      + ((1.0 - ((4. * c_ratio * c_ratio))) * inflow * 1000000000. / (24. * 3600.));
                        }

                        //new outflow
                        //reservoir storage volume should not be less than 10% of maximum!)...
                        if (G_gloResStorage[n] >= (G_stor_cap[n] * 0.1))
                            outflow = release * (24. * 3600.) / 1000000000.; //m3/s  ->  km3/routing time step
                            //...otherwise outflow will be reduced!
                            //(outflow should not be stopped, because we should serve ecosystem demands)
                        else
                            outflow = 0.1 * release * (24. * 3600.) / 1000000000.;

                        if (outflow <
                            0.) // to prevent that outflow (= river availability) gets negative in rare cases when net uses are negative downstream of a irrigation reservoir
                            outflow = 0.;

                        // substract outflow from storage
                        G_gloResStorage[n] -= outflow;

                        // reduce G_gloResStorage to maximum storage capacity
                        if (G_gloResStorage[n] > maxStorage) {
                            outflow += (G_gloResStorage[n] - maxStorage);
                            G_gloResStorage[n] = maxStorage;
                        }

                        // if G_gloResStorage is below '0' there is no outflow anymore
                        // 'gloResEvapoReductionFactor' should prevent that G_gloResStorage falls below zero.
                        if (G_gloResStorage[n] < 0.) {
                            outflow += G_gloResStorage[n];
                            G_gloResStorage[n] = 0.;
                        }

                        //daily outflow from reservoirs
                        //G_daily_res_out = (outflow * 1000000000.); // [m3]
                        G_daily_res_out[n] = outflow; // [km3/day]

                        // If (PETgwr > PETgwrMax), PET is reduced by PETgwrRedFactor.
                        // If a reduction is not required, PETgwrRedFactor is set to 1.
                        if (PETgwr == 0.)
                            PETgwrRemUseRedFactor = 1.;
                        else
                            PETgwrRemUseRedFactor = PETgwrMax / PETgwr;

                        if (PETgwrMax > PETgwr) {    // -> reduction of PET not required
                            PETgwrRemUseRedFactor = 1.;
                        }


                        if (2 == options.grid_store)
                            // sum up the negative fraction of the potential cell runoff which has not been realized
                            // due to reduction of potential evaporation
                            G_monthlyPotCellRunoffDeficit[n][month]
                                           += (0.0 - (dailyWaterBalance.G_openWaterPET[n] * PETgwrRemUseRedFactor // ???
                                                      * (1.0 - gloResEvapoReductionFactor)
                                                      * dailyWaterBalance.G_cellCorrFact[n]
                                                      * G_reservoir_area[n]
                                                      / 1000000.0));

                        // reintroduced CFA

                        inflow = outflow; // outflow of reservoir is inflow to global wetlands

                        //Redistribute remainingUseRes over all riparian cells

#ifdef _WATERGAP_CHECKS_GLOBAL_H
                                                                                                                                                // Check Felix Portmann 2015 - for invalid data
                        if ( (year >= chk_year) && (month >= chk_month) && (day_in_month >= chk_day_in_month) ) {
                            check_value_double_LT0(inflow);
                            check_value_double_LT0(outflow);
                        }
#endif


                        // introduced for 22b, Option 28 "aggrNUsGloLakResOpt"
                        dailyActualUse += remainingUseResStart - remainingUseRes;

                        // remainingUseRes is redistributed over the riparian cells (and the outflow cell if the initial
                        // dailydailyNUs was larger than zero):
                        if (options.aggrNUsGloLakResOpt == 1) {

                            G_remainingUse[n] = remainingUse;
                            G_remainingUseRes[n] = remainingUseRes;
                            G_remainingUseResRedistr[n] = remainingUseRes;
                            G_totalDesiredUseCell[n] = totalDesiredUseCell;

                            for (int i = 0; i < ng; i++) {

                                if (G_glwdUnit[n] == G_glwdUnit[i]) {

                                    if (G_remainingUseResRedistr[n] > 0.) {

                                        if (G_totalDesiredUseCell[n] <=
                                            0.) {        // different computation of G_unsatUseRiparian[i], if G_totalDesiredUseCell[n] of the OUTFLOW CELL is negative - see WG22b documentation

                                            if (G_dailydailyNUs[i] <= 0.) {

                                                G_unsatUseRiparian[i] = 0.;

                                            } else {

                                                G_unsatUseRiparian[i] = (G_dailydailyNUs[i] /
                                                                         (G_dailydailyNUsAggregated[n] -
                                                                          G_dailydailyNUs[n])) *
                                                                        G_remainingUseResRedistr[n]; // (G_remainingUse[n] - G_dailydailyNUs[n]) always larger than zero in this loop
                                                // G_unsatUseRiparian is also computed for outflow cells (when i == n), but it is not used in the code for outflow cells.
                                                G_remainingUseRes[n] = 0.;
                                                G_unsatUseRiparian[n] = 0.;
                                            }

                                        } else {    // G_totalDesiredUseCell[n] in outflow cell > 0 (G_totalDesiredUseCell[n] includes G_dailydailyNUs and (depending on selected options) allocated use from nbc and totalUnsatisfiedUse[n] (from previous day)

                                            G_unsatUseRiparian[i] = G_dailydailyNUs[i] / G_remainingUse[n] *
                                                                    G_remainingUseResRedistr[n]; // G_AllocatedUse OR G_totalUnsatisfiedUse are also distributed over riparian cells.

                                            if (G_dailydailyNUs[i] < 0.)
                                                G_unsatUseRiparian[i] = 0.;

                                            G_remainingUseRes[n] *= (G_totalDesiredUseCell[n] /
                                                                     G_remainingUse[n]);    // remainingUse to be satisfied from river storage

                                            G_unsatUseRiparian[n] = 0.;        // CR: Required since G_unsatUseRiparian[n] is added to 'remainingUse' at the beginning of each time step for all cells with GLWD units.

                                        }
                                        //  end if (G_remainingUseResRedistr[n] > 0.)
                                    } else {
                                        G_unsatUseRiparian[i] = 0.;
                                        // required? -> already covered by G_unsatUseRiparian[i] = 0?
                                        G_unsatUseRiparian[n] = 0.;
                                        G_remainingUseRes[n] = 0.;
                                    }
                                    // end if (G_remainingUseResRedistr[n] == 0.)

                                }
                                // end of loop over 'i' with matching GLWD unit

                            }
                            // end of loop over 'i'

                            remainingUse = G_remainingUseRes[n];

                            // end if (options.aggrNUsGloLakResOpt == 1)
                        } else {
                            remainingUse = remainingUseRes;
                        }

                        //Evaporation reduction factor is now calculated for the next day.
                        G_gloResEvapoReductionFactor[n] = 1. - pow(fabs(G_gloResStorage[n] - maxStorage)
                                                        / maxStorage, evapoReductionExpReservoir);

                        if (G_gloResEvapoReductionFactor[n] < 0.)    // (could occur due to numerical inaccuracies)
                            G_gloResEvapoReductionFactor[n] = 0.;

                        if (G_gloResEvapoReductionFactor[n] > 1.)  // (could occur due to numerical inaccuracies)
                            G_gloResEvapoReductionFactor[n] = 1.;

                    }
                    // end if (G_reservoir_area[n] > 0)
                    // reservoir algorithm end



                    // outflow of a global lake is inflow to global wetlands


                    // ------------------  Routing through global wetlands --------------------------------------------
                    if (G_glo_wetland[n] > 0) {

                        // assign correct storage value (for t) which was calculated in the last timestep or maxStorage at model day one
                        G_gloWetlStoragePrevStep = G_gloWetlStorage[n];

                        // maximum storage capacity
                        // Cell-specific active wetland depth - Use array (initialized from cell-specific parameter)
                        maxStorage = ((G_glo_wetland[n]) / 100.) * cellArea * G_wetlDepthActive[n];

                        // This factor was added to reduce PET as a function of actual wetland storage (2.1f).
                        // use reduction factor from previous time step (or the initial value)
                        gloWetlAreaReductionFactor = G_gloWetlAreaReductionFactor[n];


                        // G_dailyGloWetlEvapo[n]: see local lakes for explanation

                        G_dailyGloWetlEvapo[n] =
                                       ((1.0 - dailyWaterBalance.G_cellCorrFact[n])                            // [mm]
                                        * (dailyWaterBalance.G_openWaterPrec[n] * gloWetlAreaReductionFactor))
                                       + (dailyWaterBalance.G_cellCorrFact[n]
                                          * (dailyWaterBalance.G_openWaterPET[n] * gloWetlAreaReductionFactor));


                        // calculate inflow only in dependence of reduction factor to avoid inconsistent counting of precipitation
                        //   if (options.aridareaOpt == 1) // same handling of CFA like in local lakes, no further comments here
                        totalInflow = inflow + (dailyWaterBalance.G_openWaterPrec[n] * gloWetlAreaReductionFactor
                                                * (cellArea / 1000000.) * (G_glo_wetland[n] / 100.));    //[mm]->[km3]

                        //					else // no gwr_swb means no reduction-factor based variation of land area fraction
                        //						totalInflow = inflow + ((double)dailyWaterBalance.G_openWaterPrec[n]
                        //												* (double)dailyWaterBalance.G_cellCorrFact[n]//[mm]->[km3]
                        //                                                * (cellArea / 1000000.) * ((double) G_glo_wetland[n] / 100.));


                        G_glowetlextent[n] = gloWetlAreaReductionFactor * cellArea * (G_glo_wetland[n] / 100.);
                        // HMS 2013-11-22 gwr below surface water bodies
                        if ((1 == options.aridareaOpt) && (1 == G_aindex[n])) {
                            gwr_glowet[n] = Kswbgw * gloWetlAreaReductionFactor * G_glo_wetland[n] / 100. /
                                            (geo.G_contfreq[n] / 100.); // gwr below global wetlands in mm/d

                        } else
                            gwr_glowet[n] = 0.; // -> to avoid repeating the following equations for options.aridareaOpt==0

                        PETgwr = G_dailyGloWetlEvapo[n] * (cellArea / 1000000.) *
                                 ((G_glo_wetland[n]) / 100.)                        //[mm]->[km3]
                                 + gwr_glowet[n] * cellArea * (geo.G_contfreq[n] / 100.) / 1000000.;    //[mm -> km3]

                        // New approach in WG22b: In inland sinks, outflow from SWB flows into the river.
                        // The computed transported volume (= outflow of the river) represents the water resource of this cell and thus should not be zero (although
                        // it cannot be assigned to a downstream cell). This value can be interpreted as either additional AET or river availability of the inland sink.

                        // There are two options to compute PETgwrMax: Based on the analytical solution of 1) dS/dt=In-Out-k*S; 2) dS/dt=In-Out
                        // In the first equation, the outflow is appr. 0.6% of S(t-1), although the storage is empty at the end of the timestep.
                        // In this case, it seems more realistic if the remaining water evaporates or infiltrates (second equation).

                        // PETgwrMax = totalInflow + (glo_storageFactor * G_gloWetlStoragePrevStep * exp(-1. * glo_storageFactor)) / (1. - exp(-1 * glo_storageFactor));

                        PETgwrMax = totalInflow + G_gloWetlStoragePrevStep;

                        if (PETgwr > PETgwrMax) {    // reduction of PET and gwr_glowet required

                            G_gloWetlStorage[n] = 0.;
                            outflow = 0.;
                            gwr_glowet[n] *= PETgwrMax / PETgwr;        // [mm]
                            G_dailyGloWetlEvapo[n] *= PETgwrMax / PETgwr;        // [mm]

                        } else {

                            // Cell-specific calibration parameters - Use parameter
                            // glo_storageFactor (now: only one factor for surface water P_SWOUTF_C) in [1/d]
                            G_gloWetlStorage[n] = G_gloWetlStoragePrevStep * exp(-1. * P_SWOUTF_C__at__n)
                                                  + (1. / P_SWOUTF_C__at__n) * (totalInflow - PETgwr)
                                                    * (1. - exp(-1. * P_SWOUTF_C__at__n));

                            outflow = totalInflow + G_gloWetlStoragePrevStep - G_gloWetlStorage[n] - PETgwr;
                        }

#ifdef _WATERGAP_CHECKS_GLOBAL_H
                                                                                                                                                // Check Felix Portmann 2015 - for invalid data
                            if ( (year >= chk_year) && (month >= chk_month) && (day_in_month >= chk_day_in_month) ) {
                                check_value_double_LT0(inflow);
								check_value_double_LT0(P_SWOUTF_C__at__n);
                                check_value_double_LT0(totalInflow);
                                check_value_double_LT0(G_gloWetlStoragePrevStep);
                                check_value_double_LT0(G_gloWetlStorage[n]);
                                check_value_double_LT0(PETgwr);
                                check_value_double_LT0(outflow);

                                check_value_double_LT0(maxStorage);
                            }
#endif

                        if (G_gloWetlStorage[n] > maxStorage) {
                            outflow += (G_gloWetlStorage[n] - maxStorage);
                            G_gloWetlStorage[n] = maxStorage;
                        }


                        // If (PETgwr > PETgwrMax), PET is reduced by PETgwrRedFactor. If a reduction is not required,
                        // PETgwrRedFactor is set to 1.
                        // if PETgwr is 0 (because gwr_swb is 0 because reductionfactor is 0 because storage is empty),
                        // division through 0 --> has to be catched, because otherwise G_CELL_RUNOFF_TOTAL
                        // and G_CELL_AET gets nan.
                        if (PETgwr == 0.)
                            PETgwrRedFactor = 1.;
                        else
                            PETgwrRedFactor = PETgwrMax / PETgwr;

                        if (PETgwrMax > PETgwr)    // -> reduction of PET not required
                            PETgwrRedFactor = 1.;

                        if (2 == options.grid_store)
                            // sum up the negative fraction of the potential cell runoff which has not been realized
                            // due to reduction of potential evaporation
                            G_monthlyPotCellRunoffDeficit[n][month]
                                           += (0.0 - dailyWaterBalance.G_openWaterPET[n] * PETgwrRedFactor
                                                     * (1.0 - gloWetlAreaReductionFactor)
                                                     * dailyWaterBalance.G_cellCorrFact[n]
                                                     * (cellArea / 1000000.0)
                                                     * (G_glo_wetland[n] / 100.0));


                        inflow = outflow;

                        // Area reduction factor is now calculated for the next day.
                        // Cell-specific calibration parameters - Apply multiplier  // FP
                        G_gloWetlAreaReductionFactor[n] = 1. -
                                                          pow(fabs(G_gloWetlStorage[n] - maxStorage)
                                                              / maxStorage,
                                                              (M_EVAREDEX__at__n_routing * evapoReductionExp));

                        if (G_gloWetlAreaReductionFactor[n] < 0.)  // (could occur due to numerical inaccuracies)
                            G_gloWetlAreaReductionFactor[n] = 0.;

                        if (G_gloWetlAreaReductionFactor[n] > 1.)  // (could occur due to numerical inaccuracies)
                            G_gloWetlAreaReductionFactor[n] = 1.;

                    }    // end if (G_glo_wetland[n] > 0)



                    // If GWR under SWB / fractionalRouting is allowed,
                    // GW storage is computed before local lakes in humid cells (baseflow directed into local lake)
                    // and before the river in (semi-)arid cells (baseflow directed into river)

                    if ((1 == options.aridareaOpt) && (1 == G_aindex[n])) {        // (semi-)arid cells

                        // sum up total groundwater recharge [mm/d] from surface water bodies and calculate fraction
                        // of actual surface waterbody area
                        // gwr_swb is added to the groundwater storage in the next timestep
                        G_gwr_swb[n] = gwr_loclak[n] + gwr_glolak[n] + gwr_locwet[n] + gwr_glowet[n] +
                                       gwr_res[n]; // sum up gwr but do not use gwr_rivers because rivers are the final storage [mm/d]
                        G_groundwaterStoragePrevStep = G_groundwaterStorage[n];

                        // WG22b: ODE (ordinary differential equation) is applied for groundwater storage:
                        // dS/dt = GWR - NAg - k*S  (k*S = groundwater_runoff_km3)

                        // Groundwater recharge below surface water bodies from previous timestep:
                        netGWin = G_gwr_swb[n] * cellArea * (geo.G_contfreq[n] / 100.) /
                                  1000000.                        // mm -> km3
                                  + (double) dailyWaterBalance.G_dailyGwRecharge[n] * cellArea *
                                    ((double) G_landAreaFrac[n] / 100.) / 1000000.;    // mm -> km3

                        // G_dailydailyNUg is adapted in this timestep based on the remaining use of last time step
                        if (options.subtract_use > 0) {
                            if (G_dailyRemainingUse[n] != 0) {
                                G_dailydailyNUg[n] = updateNetAbstractionGW(n, month);
                            }
                            netGWin -= G_dailydailyNUg[n];
                        }

                        // Analytical solution of dS/dt = netGWR - k*S
                        // k_g (now: P_GWOUTF_C) in [1/d]
                        // Cell-specific calibration parameters - Use parameter

                        G_groundwaterStorage[n] = G_groundwaterStoragePrevStep * exp(-1. * P_GWOUTF_C__at__n_routing) +
                                                  (1. / P_GWOUTF_C__at__n_routing) * netGWin *
                                                  (1. - exp(-1. * P_GWOUTF_C__at__n_routing));

                        groundwater_runoff_km3 = G_groundwaterStoragePrevStep - G_groundwaterStorage[n] + netGWin;

                        if (groundwater_runoff_km3 <=
                            0.) {    // different differential equation (ODE) applies:
                                     // dS/dt = GWR - NAg (without k*S) -> S(t) = S(t-1) + netGWin

                            groundwater_runoff_km3 = 0.;

                            G_groundwaterStorage[n] = G_groundwaterStoragePrevStep + netGWin;
                        }
                        if ((2 == options.grid_store) || (4 == options.grid_store) || (6 == options.grid_store))
                            G_monthlyGwRunoff[n][month] += groundwater_runoff_km3;
                        // if daily 365 output and gwout (same for 31; same for land area fraction)
                        // G_daily365GwRunoff[n][day - 1] = groundwater_runoff_km3;

                        if (1 ==
                            options.fractionalRoutingOpt) { // only a fraction (depending on initial amount of surface water bodies) of surface and groundwater runoff should go to the local lakes later, rest should go directly to the river
                            G_fswb_catchment[n] = G_fswbInit[n] *
                                                  20.;                  // G_fswb_catchment[n] = G_fswb[n] * 5.; // add a catchment area to surface water bodies (specified as x times the size of swb)...
                            if (G_fswb_catchment[n] > 1.) // ... and prevent that G_fswb is greater then 1
                                G_fswb_catchment[n] = 1.;
                            if (G_aindex[n] ==
                                1) {// in semi-arid/arid areas, all groundwater reaches the river directly
                                G_localGWRunoffIntoRiver[n] = groundwater_runoff_km3;
                            } else {// in humid areas, groundwater is also routed through surface water bodies
                                G_localGWRunoffIntoRiver[n] = (1. - G_fswb_catchment[n]) * groundwater_runoff_km3;
                            }

                            if (G_aindex[n] ==
                                0) // route groundwater only through surface water bodies in humid regions
                                G_localGWRunoff[n] = G_fswb_catchment[n] * groundwater_runoff_km3;

                            else
                                G_localGWRunoff[n] = 0.; // in arid regions, all groundwater flows directly to the river

                            // Firstly, amount of local runoff flowing directly into the river is calculated, then local gw runoff is added for the rest term.

                        }

                        // to be changed for fractionalrouting option (normal routing without G_localRunffIntoRiver and
                        // G_localRGWRunoffIntoRiver
                        // reintroduced for CFA calculation
                        potCellRunoff = (G_localRunoff[n] + G_localRunoffIntoRiver[n] +
                                         G_localGWRunoffIntoRiver[n]    // is in km3
                                         + ((dailyWaterBalance.G_lakeBalance[n]    // is in mm
                                             / 1000000.0)    // mm ->*/ km
                                            * (((double) G_lake_area[n] + (double) G_reservoir_area[n])    // km2
                                               + ((geo.areaOfCellByArrayPos(n) / 100.0)    // % -> km2
                                                  * (G_loc_lake[n] + G_loc_wetland[n] + G_glo_wetland[n])))));    // [%]

                        G_potCellRunoff[n] += potCellRunoff;

                        if (2 == options.grid_store) {
                            G_monthlyPotCellRunoff[n][month] += potCellRunoff;
                        }

                    }
                    // end of computation of GW storage and baseflow for options.aridareaOpt == 1 && (semi-)arid cells



                    // calculate the amount of water which flows out of local/global lakes/wetlands
                    // and contributes to the actual cell runoff within the cell
                    // Different from WG22, cell runoff is computed as (transportedVolume - inflow from upstream cells) in WG2.2b .
                    // if (1 == options.fractionalRoutingOpt) // reduce G_riverInflow by the amount of water which flows not through the storage cascade
                    // CellRunoff = (inflow - G_riverInflow[n] + G_localRunoffIntoRiver[n] + G_localGWRunoffIntoRiver[n]);

                    // else  // normal version (without fractional routing)
                    // CellRunoff = (inflow - G_riverInflow[n]);

                    // if ((0 == options.fractionalRoutingOpt) && (G_aindex[n] == 1) && (options.aridareaOpt == 1)) // in arid regions with gwr_swb groundwater should run directly into the river
                    // CellRunoff = (inflow - G_riverInflow[n] + G_localGWRunoffIntoRiver[n]);

                    // outflow from global wetlands is inflow to the river system
                    G_riverInflow[n] = inflow;



                    // At this point G_riverInflow is the total inflow to the river system:
                    // runoff->loca_lakes->local_wetlands->global_lakes->global_wetlands->river_network

                    if (1 ==
                        options.fractionalRoutingOpt) { // if option is activated, add the water which is flowing directly to the river
                        G_riverInflow[n] += G_localRunoffIntoRiver[n];
                        G_riverInflow[n] += G_localGWRunoffIntoRiver[n];
                    }
                    if ((0 == options.fractionalRoutingOpt) && (G_aindex[n] == 1) && (options.aridareaOpt == 1)) {
                        G_riverInflow[n] += G_localGWRunoffIntoRiver[n];
                    }
                    // riverVelocity in [km/d]

                    if (options.riverveloOpt == 0) {
                        riverVelocity = defaultRiverVelocity;
                    }
                    else if (options.riverveloOpt == 1) {
                        // Cell-specific calibration parameters - Apply multiplier M_RIVRGH_C within method getRiverVelocity // FP
                        riverVelocity = getRiverVelocity(G_RiverSlope[n], G_riverBottomWidth[n], G_Roughness[n],
                                                         G_riverInflow[n], n);
                    }
                    else if (options.riverveloOpt == 2) {
                        riverVelocity = getNewRiverVelocity(G_RiverSlope[n], G_riverBottomWidth[n], G_Roughness[n],
                                                            G_riverStorage[n], G_riverLength[n], n);
                    }

                    // WG22b: To be consistent with global lakes/wetlands/groundwater (storages with ord. diff.
                    // equation), K is expressed in 1/d.
                    K = riverVelocity / G_riverLength[n]; // [(km/d)/km] = [1/d]
                    // K = G_riverLength[n] / riverVelocity; // [km / (km/d)] = [d]
                    // K = G_celldistance[geo.G_row[n] - 1][2] / riverVelocity; // [km / (km/d)] = [d]


                    // assign correct storage value which was calculated in the last timestep or maxStorage at model day one
                    G_riverStoragePrevStep = G_riverStorage[n];


                    // Computation of river evaporation prior to river storage equation to implement G_dailyRiverEvapo
                    // into storage equation.
                    // River evaporation start
                    if (options.riverEvapoOpt == 1) {

                        // use reduction factor from previous time step (or the initial value)
                        //riverAreaReductionFactor = G_riverAreaReductionFactor[n];	// riverAreaReductionFactor is not used anymore

                        G_riverAreaFrac[n] = G_riverAreaFracNextTimestep_Frac[n];

                        //HMS both calculated later if riverareafrac is finally calculated
                        G_dailyRiverEvapo[n] = ((1.0 - dailyWaterBalance.G_cellCorrFact[n]) *
                                                (dailyWaterBalance.G_openWaterPrec[n])
                                                + (dailyWaterBalance.G_cellCorrFact[n] *
                                                   dailyWaterBalance.G_openWaterPET[n]))
                                               * G_riverAreaFrac[n] / 100. * cellArea / 1000000.;

                        //dailyRiverEvapo = G_dailyRiverEvapo[n] * cellArea *(((double)geo.G_landfreq[n] + (double)geo.G_fwaterfreq[n]) / 100. / 1000000.); //evaporation from rivers in km3
                        //                  if (transportedVolume > dailyRiverEvapo)
                        //                  transportedVolume -= G_dailyRiverEvapo[n] * cellArea *(((double)geo.G_landfreq[n] + (double)geo.G_fwaterfreq[n]) / 100. / 1000000.); // subtract dailyRiverevapo from discharge
                        G_dailyRiverPrecip[n] =
                                       dailyWaterBalance.G_openWaterPrec[n] * G_riverAreaFrac[n] / 100. * cellArea /
                                       1000000.; // precip on river fraction in km3

                    } else {// (no evaporation from rivers)
                        G_dailyRiverPrecip[n] = 0.;
                        G_dailyRiverEvapo[n] = 0.;
                    }

                    G_riverInflow[n] += G_dailyRiverPrecip[n];

                    // At this point, RiverEvapoRemUse can get negative because remainingUse can be negative!
                    RiverEvapoRemUse = remainingUse +
                                       G_dailyRiverEvapo[n]; // [km3] // remainingUse = 0, if options.subtract_use == 0

                    // introduced for 22b, Option 28 "aggrNUsGloLakResOpt"
                    remainingUseRiverStart = remainingUse;



                    // WG22b: Unit of K changed from [d] to [1/d]
                    // Different from global wetlands, RivEvapoRemUseMax is computed based on the ODE dS/dt=In-Out-k*S.
                    // (outflow is assumed to occur even if S = 0 at the end of the time step.
                    // See comment in global wetland loop)
                    RivEvapoRemUseMax = G_riverInflow[n] + (K * G_riverStoragePrevStep * exp(-1. * K))
                                                           / (1. - exp(-1. * K));

                    if (RiverEvapoRemUse > RivEvapoRemUseMax) {

                        G_riverStorage[n] = 0.;

                        transportedVolume = G_riverInflow[n] + G_riverStoragePrevStep - RivEvapoRemUseMax;

                        if (remainingUse > 0.)
                            remainingUse -= remainingUse * RivEvapoRemUseMax / RiverEvapoRemUse;
                        else
                            remainingUse = 0.;

                        G_dailyRiverEvapo[n] *= RivEvapoRemUseMax / RiverEvapoRemUse; // [km3]

                    } else {


                        // WG22b: Unit of K changed from [d] to [1/d]
                        G_riverStorage[n] = G_riverStoragePrevStep * exp(-1. * K)
                                            + (1. / K) * (G_riverInflow[n] - RiverEvapoRemUse)
                                              * (1. - exp(-1. * K));

                        transportedVolume =
                                       G_riverInflow[n] + G_riverStoragePrevStep - G_riverStorage[n] - RiverEvapoRemUse;

                        remainingUse = 0.;

                    }


                    if (transportedVolume <
                        0.) {    // (due to numerical inaccuracies)
                        cout << "ERROR transported volume < 0 in cell: " << n << "routing.routing: year-month-day "
                             << year << "-" << month << "-" << day_in_month << " - julian day: " << day << " with: "
                             << transportedVolume << endl;

                        cout << "transportedVolume = G_riverInflow[n]: " << G_riverInflow[n]
                             << " + G_riverStoragePrevStep: " << G_riverStoragePrevStep << endl;
                        cout << " - G_riverStorage[n]: " << G_riverStorage[n] << " - RiverEvapoRemUse: "
                             << RiverEvapoRemUse << endl;
                        cout << "RivEvapoRemUseMax: " << RivEvapoRemUseMax << ", remainingUse (set to zero): "
                             << remainingUse << endl;
                        cout << "K: " << K << ", riverVelocity: " << riverVelocity << ", G_riverLength[n]: "
                             << G_riverLength[n] << endl;
                        cout << "G_riverAreaFrac[n]: " << G_riverAreaFrac[n] << ", G_dailyRiverEvapo[n]: "
                             << G_dailyRiverEvapo[n] << ", G_dailyRiverPrecip[n]: " << G_dailyRiverPrecip[n] << endl;


                        // transportedVolume = 0.; // set to 0 in those cases.

                    }


#ifdef _WATERGAP_CHECKS_GLOBAL_H
                                                                                                                                            // Check Felix Portmann 2015 - for invalid data
                    // Check when each grid cell is accessed
                    // From this point on, remainingUse should be positive or zero!
                    if ( (year >= chk_year) && (month >= chk_month) && (day_in_month >= chk_day_in_month) ) {
                        check_value_double_LT0(RivEvapoRemUseMax);
                        check_value_double_LT0(transportedVolume);
                        check_value_double_LT0(remainingUse);
                        check_value_double_LT0(remainingUseRes);
                        check_value_double_LT0(riverVelocity);
                        check_value_double_LT0(G_riverStoragePrevStep);
                        check_value_double_LT0(K);

                        check_value_double_NaN(RiverEvapoRemUse);
                        check_value_double_NaN(RivEvapoRemUseMax);
                        check_value_double_NaN(transportedVolume);
                        check_value_double_NaN(remainingUse);
                        check_value_double_NaN(remainingUseRes);
                        check_value_double_NaN(riverVelocity);
                        check_value_double_NaN(G_riverStoragePrevStep);
                        check_value_double_NaN(K);
                    }
                    // end Check Felix Portmann 2015
#endif

                    // Calculate river width (based on storage at the end of current time step) according to new river
                    // flow velocity approach (not implemented for river velocity yet,
                    // but documented in documentation for WG22b)

                    // In WG22, G_RiverWidth was a function of "transportedVolume" and could therefore drop to zero.
                    // G_RiverWidth[n] = 2.71 * pow(transportedVolume * 1/0.0000864, 0.557); // calculate G_RiverWidth in m after Allen et al. 1994; transportedVolume has to be converted from km3/d into m3/s (1/((60*60*24)/1000000000)
                    // In WG22b, the minimum river width is the river bottom width meaning that river evaporation could occur even if the river depth is zero.
                    // Therefore, G_riverAreaFrac[n] is reduced by the G_riverAreaReductionFactor.

                    if (options.riverEvapoOpt == 1) {

                        crossSectionalArea = G_riverStorage[n] / G_riverLength[n];    // [km3/km = km2]
                        // G_riverBottomWidth[n]: m -> km
                        riverDepth = -G_riverBottomWidth[n] / (4. * 1000.) +
                                     sqrt(G_riverBottomWidth[n] / 1000. * G_riverBottomWidth[n] / (16. * 1000.) +
                                          0.5 * crossSectionalArea); // [km]
                        G_RiverWidth[n] = G_riverBottomWidth[n] / 1000. + 4. * riverDepth; //[km]

                        if (G_RiverWidth[n] > G_RiverWidth_bf[n] / 1000.) //G_RiverWidth_bf[n] in m
                            G_RiverWidth[n] = G_RiverWidth_bf[n] / 1000.;


                        // G_riverAreaReductionFactor and G_riverAreaFracNextTimestep are now calculated for the next day.
                        // Cell-specific calibration parameters - Apply multiplier  // FP
                        G_riverAreaReductionFactor[n] = 1. -
                                                        pow(fabs(G_riverStorage[n] - G_riverStorageMax[n])
                                                            / G_riverStorageMax[n],
                                                            (M_EVAREDEX__at__n_routing * evapoReductionExp));

                        if (G_riverAreaReductionFactor[n] < 0.)  // (could occur due to numerical inaccuracies)
                            G_riverAreaReductionFactor[n] = 0.;

                        if (G_riverAreaReductionFactor[n] > 1.)  // (could occur due to numerical inaccuracies)
                            G_riverAreaReductionFactor[n] = 1.;


                        // Only continental cells
                        if (geo.G_contcell[n]) {
                            G_riverAreaFracNextTimestep_Frac[n] =
                                           G_riverAreaReductionFactor[n] * G_riverLength[n] * G_RiverWidth[n] * 100. /
                                           cellArea;
                        }
                            // 100% water cells (of Caspian Sea) have now G_fwaterfreq + G_landfreq = 0. The water balance is calculated in the outlet cell of large lakes; therefore a calculation of e.g. precipitation for this cells leads to inconsistent quantification of global fluxes. A division through 0 has to be avoided now (alternative: define it as not land cells in ng).
                            // to avoid not defined G_riverAreaFrac at Caspian Sea
                        else {
                            G_riverAreaFracNextTimestep_Frac[n] = 0.;
                        }

                        // Calculate the necessary percent change toward the next time step
                        G_riverAreaFrac_ChangeFrac[n] = G_riverAreaFracNextTimestep_Frac[n] - G_riverAreaFrac[n];

                    }

                    // WG22b: Different from global lakes and reservoirs, remaining use (if > 0) is subtracted
                    // from local lake storage at the end of routing.
                    if (1 <= options.subtract_use) {

                        if ((G_loc_lake[n] > 0.) && (remainingUse > 0.)) {

                            // WG22b: Water can be abstracted as long as locLakeStorage exceeds -maxStorage
                            // (and not zero, as it was the case in all previous versions).
                            if (G_locLakeStorage[n] > (-1.) * G_locLakeMaxStorage[n]) {
                                if (remainingUse < (G_locLakeMaxStorage[n] + G_locLakeStorage[n])) {
                                    G_locLakeStorage[n] -= remainingUse;
                                    remainingUse = 0.;
                                } else {
                                    remainingUse -= (G_locLakeMaxStorage[n] + G_locLakeStorage[n]);
                                    G_locLakeStorage[n] = (-1.) * G_locLakeMaxStorage[n];
                                }
                            }

                            // Area reduction factor for the next day is calculated again based on the reduced local
                            // lake storage.
                            // Cell-specific calibration parameters - Apply multiplier  // FP
                            G_locLakeAreaReductionFactor[n] = 1. -
                                                              pow(fabs(G_locLakeStorage[n] - G_locLakeMaxStorage[n])
                                                                  / (2. * G_locLakeMaxStorage[n]),
                                                                  (M_EVAREDEX__at__n_routing * evapoReductionExp));

                            if (G_locLakeAreaReductionFactor[n] <
                                0.)
                                G_locLakeAreaReductionFactor[n] = 0.;

                            if (G_locLakeAreaReductionFactor[n] > 1.)  // (could occur due to numerical inaccuracies)
                                G_locLakeAreaReductionFactor[n] = 1.;

                        }
                        // end of abstraction from local lakes ((G_loc_lake[n] > 0.) && (remainingUse > 0.))

                        // introduced for 22b, Option 28 "aggrNUsGloLakResOpt"
                        dailyActualUse += remainingUseRiverStart - remainingUse;



                        // satisfied use of the current day summed up for the whole year
                        // G_satisfiedUse[n] += totalDesiredUse - remainingUse;	// WG22
                        // for new Use allocation (M.Hunger 2/2006)
                        // all actual satisfied water uses
                        // G_actualUse[n] += totalDesiredUse - remainingUse;	// WG22

                        // WG22b: REDISTRIBUTE UNSATISFIED ALLOCATED USE OVER RESPECTIVE NEIGHBORING CELLS:
                        // Rationale: First priority: Satisfy water demand of cell n (from water storage in cell n).
                        // Second priority: Satisfy water demand allocated from neighboring cell(s) (from water storage
                        // in cell n).

                        if ((options.use_alloc > 0) && (G_dailyAllocatedUse[n] > 0.)) {
                        //if (G_dailyAllocatedUse[n] > 0.) {    // only possible if (options.use_alloc > 0)

                            if (remainingUse > 0.) {

                                if (demandWithoutAllocUse > 0.) {

                                    unsatUse = remainingUse;

                                    remainingUse = demandWithoutAllocUse - dailyActualUse;    // remainingUse of cell n


                                    if (remainingUse <
                                        0.)    // since "dailyActualUse" can exceed "demandWithoutAllocUse"
                                        remainingUse = 0.;

                                    unsatAllocUse = unsatUse -
                                                    remainingUse;    // total unsatisfiedUse (including allocated use)
                                                                     // minus remainingUse of cell n

                                } else {    // if (demandWithoutAllocUse <= 0.) -> negative G_dailydailyNUs

                                    // The initial water demand of cell n (demandWithoutAllocUse) is <= 0.
                                    // In this case, the remainingUse at this point only stems from the respective
                                    // neighboring cells.
                                    // The remainingUse of cell n is set to zero.
                                    unsatAllocUse = remainingUse;    // 	remainingUse to be redistributed over respective nbc

                                    remainingUse = 0.;    // remainingUse of cell n

                                }

                            }
                                // end (remainingUse > 0.)

                            else {    // if (remainingUse == 0)

                                unsatAllocUse = 0.;
                            }


                            for (i_nbc = 0; i_nbc < ng; i_nbc++) {

                                if (G_SecondCell[i_nbc] == n) {
                                    // Scale G_UnsatAllocUse[i_nbc] to just the allocated use coming from cell i_nbc
                                    G_UnsatAllocUse[i_nbc] = G_AllocUseToNeigborCell[i_nbc] * unsatAllocUse /
                                                             G_dailyAllocatedUse[n];    // At this point, G_dailyAllocatedUse[n] is always > 0.
                                    if (options.delayedUseSatisfaction == 0){
                                        G_dailyRemainingUse[i_nbc] = G_UnsatAllocUse[i_nbc];
                                    }else{
                                        // The adaption of NAg shall be done on a daily basis, hence the difference from
                                        // prev timestep is calculated. It can happen that secondCell is used in one
                                        // timestep and not in another thus both has to be subtracted.
                                        G_dailyRemainingUse[i_nbc] = (G_UnsatAllocUse[i_nbc] -
                                                                      G_PrevUnsatAllocUse[i_nbc] -
                                                                      G_PrevTotalUnstatisfiedUse[i_nbc]);
                                    }

                                    // In order to have information about the total satisfied water use in the
                                    // cell where it was demanded an array is created which contains information
                                    // about which saturated allocated use should be subtracted and which added where
                                    // to have this as output option.
                                    double sat_alloc_use = (G_AllocUseToNeigborCell[i_nbc] -
                                                                         G_UnsatAllocUse[i_nbc]);
                                    G_monthlyRedistributeAllocUse[i_nbc][month] += sat_alloc_use;
                                    G_monthlyRedistributeAllocUse[n][month] -= sat_alloc_use;

                                    // If "G_UnsatAllocUse[i_nbc]" will be added to the "remainingUse" of "i_nbc" the
                                    // next day, this amount of unsatisfied use must be taken into account to check the
                                    // following condition:
                                    // sum(actual use) + sum(unsat. use) = sum(NAs)

                                    if (options.delayedUseSatisfaction == 0) {

                                        G_totalUnsatisfiedUse[i_nbc] += G_UnsatAllocUse[i_nbc];
                                        G_UnsatAllocUse[i_nbc] = 0.;

                                    } else {    // options.delayedUseSatisfaction == 1:

                                        // used as output only:
                                        if (G_CellPositionRoutOrder[i_nbc] < G_CellPositionRoutOrder[n])
                                            G_daily_UnsatAllocUseNextDay[i_nbc] = G_UnsatAllocUse[i_nbc];

                                        else
                                            G_daily_UnsatAllocUseNextDay[i_nbc] = 0.;

                                        // WG22b: As in previous versions, subtraction of NUs can be delayed for up to one year.
                                        if (day == 365 &&
                                            (G_CellPositionRoutOrder[i_nbc] < G_CellPositionRoutOrder[n])) {

                                            G_totalUnsatisfiedUse[i_nbc] = G_UnsatAllocUse[i_nbc];
                                            G_UnsatAllocUse[i_nbc] = 0.;

                                        }

                                    }
                                    // end if options.delayedUseSatisfaction == 1

                                }
                                // end if (G_SecondCell[i_nbc] == n)

                            }
                            // end for (i_nbc = 0; i_nbc < ng; i_nbc++)

                        }
                        // end if (G_dailyAllocatedUse[n] > 0.): 'REDISTRIBUTE UNSATISFIED ALLOCATED USE OVER RESPECTIVE NEIGHBORING CELLS'


                        G_actualUse[n] += dailyActualUse;    // required to write annual value of G_actualUse[n]

                        // calculate satisfied allocated use in second cell for output
                        G_dailySatisAllocatedUseInSecondCell[n] = G_dailyAllocatedUse[n] - unsatAllocUse;


                        // monthly satisfied use: satisfied uses from this cell only
                        // monthly actual use: all satisfied uses (see below)
                        if ((2 == options.grid_store) || (4 == options.grid_store) || (6 == options.grid_store)) {
                            //    G_monthlySatisfiedUse[n][month] += totalDesiredUseCell - remainingUse;	// CR 2015-09: A distinction between satisfied use and actual use is not possible in 22b.
                            G_monthlyActualUse[n][month] += dailyActualUse;
                            G_monthlyNUg[n][month] += G_dailydailyNUg[n];
                            // Different from WG22, G_monthlyAllocatedUse does not necessarily have to be fully satisfied in the neighboring cell (nbc).
                            // It is only the water demand allocated to a nbc (without taking into account the current river and SWB storage in the nbc).
                            G_monthlyAllocatedUse[n][month] += G_dailyAllocatedUse[n];
                            G_monthlySatisAllocatedUseinSecondCell[n][month] += G_dailySatisAllocatedUseInSecondCell[n];
                        }


                        // ALLOCATION OF REMAINING USE TO SECOND CELL
                        // At this point, 'remainingUse' only includes unsatisfied water demand of the cell itself.
                        // Unsatisfied allocated use has been redistributed over the respective neighboring cells before.

                        if (options.use_alloc > 0) {

                            if (remainingUse > 0.) {    // Identification of secondCell

                                if (G_toBeCalculated[n] == 1) {
                                    // a second cell is allowed to satisfy the remaining demand
                                    short i;
                                    int nbc;
                                    double storageSum;

                                    // choose the neighbouring cell with the highest
                                    // amount of river and lake storage
                                    // which is not upstream of cell n
                                    totalNeighbourStorage = 0.0;
                                    secondCell = -1; // standard setting "no downstream cell"

                                    for (i = 0; i < 8; i++) {

                                        nbc = G_neighbourCell[n][i] - 1; // number of neighbour cell
                                        // For non-continent cells (Caspian Sea in SLM) G_contcell[nbc] == 0,
                                        // set nbc to unaccepted value "-1", to exclude this cell from further calculations here
                                        if (0 == geo.G_contcell[nbc]) {
                                            nbc = -1;
                                        }

                                        if (nbc < 0) {
                                            continue; //  neighbouring cell does not exist
                                        } else {
                                            if (G_toBeCalculated[nbc] == 1) {
                                                // " + G_locLakeMaxStorage[nbc]" and " + G_gloLakeMaxStorage[nbc]" inserted,
                                                // since in WG22b local/global lake storage can decrease to -maxStorage.
                                                // Otherwise, 'storageSum' could become negative, and the initial value of secondCell = -1 would not be updated (secondCell = nbc).

                                                storageSum = G_riverStorage[nbc]
                                                             + (G_locLakeStorage[nbc] + G_locLakeMaxStorage[nbc])
                                                             + (G_gloLakeStorage[nbc] + G_gloLakeMaxStorage[nbc]);

                                                if (options.resOpt == 1) {

                                                    storageSum += G_gloResStorage[nbc];

                                                }

                                                if (storageSum > totalNeighbourStorage) {

                                                    // check if neighbour cell is a direct upstream cell
                                                    if ((G_downstreamCell[nbc] - 1) != n) {
                                                        // set second cell only if it is not a direct
                                                        // upstream cell of current cell n
                                                        totalNeighbourStorage = storageSum;
                                                        secondCell = nbc;

                                                    }
                                                }

                                            }
                                            // end 'if (G_toBeCalculated[nbc] == 1)' (inner)

                                        }
                                        // endelse neighbouring cell exists

                                    }
                                    // endfor loop over eight cells (G_neighbourCell)

                                }
                                // end 'if (G_toBeCalculated[n] == 1)' (outer)

                                // WG22b: As in previous versions, subtraction of NUs can be delayed for up to one year.
                                if (day == 365 && (G_CellPositionRoutOrder[secondCell] < G_CellPositionRoutOrder[n]))
                                    secondCell = -1;

                            }
                                // end 'if (remainingUse > 0.)'	(Identification of secondCell)
                            else {
                                secondCell = -1;

                            }
                            // end 'else: remainingUse == 0.'


                            if (secondCell < 0) {    // remainingUse can be >= 0 at this point

                                G_AllocUseToNeigborCell[n] = 0.;

                                if (options.delayedUseSatisfaction == 1) {
                                    // The adaption of NAg shall be done on a daily basis, hence the difference from
                                    // prev timestep is calculated. It can happen that secondCell is used in one
                                    // timestep and not in another thus both has to be subtracted.
                                    G_dailyRemainingUse[n] = (remainingUse -
                                                              G_PrevTotalUnstatisfiedUse[n] -
                                                              G_PrevUnsatAllocUse[n]);
                                    G_totalUnsatisfiedUse[n] = remainingUse;

                                } else {    // if options.delayedUseSatisfaction == 0: G_totalUnsatisfiedUse[n] is only stored as output variable
                                    G_totalUnsatisfiedUse[n] += remainingUse;
                                    // alloc usage but no second cell identified
                                    // G_dailyRemainingUse is assigned to  adapt NAg in next timestep
                                    G_dailyRemainingUse[n] = remainingUse;
                                }


                            }
                                // end if (secondCell < 0)

                            else {    // secondCell > 0: remainingUse is always > 0. in this case

                                // G_AllocatedUse[secondCell] is increased (+= instead of =),
                                // as it may be identified as a "second cell" for other cells in the same timestep.
                                // When the remainingUse of the second cell is computed at the beginning of each
                                // timestep, G_AllocatedUse is set to zero directly after adding the value to the
                                // remainingUse.

                                G_AllocatedUse[secondCell] += remainingUse;

                                // WG22b: G_AllocUseToNeigborCell and G_SecondCell used to redistribute the unsatisfied
                                // fractions of allocated use to the original cells:

                                // If G_AllocatedUse cannot be fully satisfied in the second cell, the remainin
                                // allocated use is redistributed over the neighboring cells using
                                // G_AllocUseToNeigborCell[n] and G_SecondCell[n].
                                G_AllocUseToNeigborCell[n] = remainingUse;

                                // WG22b: G_daily_allocatedUseNextDay only used for output G_ALLOC_USE_NEXT_DAY_%d.365.UNF0:
                                //	If G_routingCell[n] > G_routingCell[secondCell],
                                //	remainingUse allocated to a nbc will be satisfied in the next timestep.
                                //	To check if water use components are written correctly, this unsatisfied allocated use
                                //	must be taken into account.
                                if (G_CellPositionRoutOrder[secondCell] <
                                    G_CellPositionRoutOrder[n])    // (Only G_AllocatedUse[secondCell] of the next time step is taken into account)
                                    G_daily_allocatedUseNextDay[secondCell] += remainingUse;

                                if (options.delayedUseSatisfaction == 1) {
                                    G_totalUnsatisfiedUse[n] = 0.;
                                    // Otherwise, remainingUse of cell n would be increased by G_totalUnsatisfiedUse[n] from the previous timestep again.
                                }
                                // else (options.delayedUseSatisfaction == 0):
                                // G_totalUnsatisfiedUse[n] is not changed in this timestep.

                            }
                            // end else (secondCell > 0): remainingUse is always > 0. in this case

                            G_SecondCell[n] = secondCell;

                        }
                            // end 'if (options.use_alloc > 0)': 'ALLOCATION OF REMAINING USE TO SECOND CELL'

                        else {    // (options.use_alloc == 0)	(remainingUse can be >= 0 at this point)

                            if (options.delayedUseSatisfaction == 1) {
                                // NAg adaption happens on a daily basis thus it is compared to the
                                // G_totalUnsatisfiedUse from last time step before assining remainingUse.
                                G_dailyRemainingUse[n] = remainingUse - G_totalUnsatisfiedUse[n];
                                G_totalUnsatisfiedUse[n] = remainingUse;
                            }
                            else {
                                // In this case, 'G_totalUnsatisfiedUse' is not added to 'remainingUse' at the beginning,
                                // but the sum of 'G_totalUnsatisfiedUse' is stored as an output variable.
                                G_totalUnsatisfiedUse[n] += remainingUse;
                                //NAg adaption based on remainingUse
                                G_dailyRemainingUse[n] = remainingUse;
                            }
                        }
                        // end else(options.use_alloc == 0)

                    }
                    // end if options.subtract_use >= 1

                    // Computation of SWB and land area fractions at the end of 'routing'.

                    //                 riverAvail = transportedVolume; // CR 2015-08-04 not used in WG22b

                    // WG22b: The river availability represents the water resource of a cell (including inland sinks).
                    // Thus, the outflow of SWB should be taken into account in inland sinks as well.

                    if (1 ==
                        options.statcorrOpt) // if station correction is used, multiply discharge with station correction (e.g. for calibration).
                        transportedVolume *= G_statCorrFact[n];

                    // WG22b: New approach to compute CellRunoff:
                    CellRunoff = (transportedVolume - inflowUpstream);


                    //-------------------------------------for daily velocity file----------------------
                    if (options.day_store == 1) {
                        for (b = 0; b <= nSpecBasins - 1; b++) {
                            if (n == cbasin.cellNum[b] - 1) {
                                //dailyRiverVelocity[b][day-1] = riverVelocity * timeStepsPerDay / 86.4;	// //km/timeStep --> km/day --> m/sec ((60*60*24)/1000)
                                dailyRiverVelocity[b][day - 1] += (riverVelocity /
                                                                   86.4); // //km/d --> m/s ((60*60*24)/1000)
                                break;
                            }
                        }
                    }

                    // ... and put the water into the downstream cell
                    if ((G_downstreamCell[n] - 1) >= 0) {
                        if (0 != G_toBeCalculated[(G_downstreamCell[n] - 1)])
                            G_riverInflow[(G_downstreamCell[n] - 1)] += transportedVolume;
                    }

                    // annualCellRunoff over the whole year
                    G_AnnualCellRunoff[n] += CellRunoff;

                    // monthly output (possibly combined with daily output) start
                    if ((2 == options.grid_store) || (4 == options.grid_store) || (6 == options.grid_store)) {
                        // G_monthlyCellRunoff should be defined as the inflow from global lakes,
                        // when calculating total inflow to reservoirs.
                        // G_monthlyCellRunoff[n][month] 			+= res_inflow;

                        // use precip on land, wetlands and local lakes and global lakes and reservoirs if they occur ans all only on continental area in km3
                        // Only continental cells
                        // each other cell that contains continental area
                        // ATTENTION/POTENTIAL ERROR: G_dailyRiverPrecip[n] is not taken into account in calculationg G_landWaterExclGloLakAreaFrac[n]


                        if (geo.G_contcell[n]) {
                            G_monthlyConsistentPrecip[n][month] += (
                                           (dailyWaterBalance.G_openWaterPrec[n] * G_landWaterExclGloLakAreaFrac[n] /
                                            geo.G_contfreq[n]) *
                                           ((geo.areaOfCellByArrayPos(n) * (geo.G_contfreq[n] / 100.)) / 1000000.) +
                                           G_GloLakePrecip[n] + G_ReservoirPrecip[n] + G_dailyRiverPrecip[n]);
                        }
                            // 100% water cells (of Caspian Sea) have now G_fwaterfreq + G_landfreq = 0. The water balance is calculated in the outlet cell of large lakes; therefore a calculation of e.g. precipitation for this cells leads to inconsistent quantification of global fluxes. A division through 0 has to be avoided now (alternative: define it as not land cells in ng).
                            // inner Caspian lake cell (not belonging to a continent)
                        else {
                            G_monthlyConsistentPrecip[n][month] = 0.;
                        }

                        // add cell runoff to monthly cell runoff
                        G_monthlyCellRunoff[n][month] += CellRunoff;

                        // Only continental cells
                        // each other cell that contains continental area
                        if (geo.G_contcell[n]) {
                            G_monthlyCellSurfaceRunoff[n][month] += ((CellRunoff - G_localGWRunoff[n])
                                                                     / (cellArea * (geo.G_contfreq[n] / 100.)) *
                                                                     1000000.);
                        }
                            // 100% water cells (of Caspian Sea) have now G_fwaterfreq + G_landfreq = 0. The water balance is calculated in the outlet cell of large lakes; therefore a calculation of e.g. precipitation for this cells leads to inconsistent quantification of global fluxes. A division through 0 has to be avoided now (alternative: define it as not land cells in ng).
                            // inner Caspian lake cell (not belonging to a continent)
                        else {
                            G_monthlyCellSurfaceRunoff[n][month] = 0.;
                        }

                        // add transportedVolume to monthly cell runoff
                        G_monthlyRiverAvail[n][month] += transportedVolume;

                        if (options.outRiverInUpstream)
                            G_monthlyRiverInUpstream[n][month] += inflowUpstream;


                        // add riverVelocity to monthlyVelocity
                        G_monthlyVelocity[n][month] += (riverVelocity / 86.4); // velocity grid [m/s]
                        // for those cells which are defined as inland sinks and do have considerable lake/wetland area (>'0'),
                        // river availability should be updated, otherwise river availability on basin scale could be '0'

                        // WG22b: The river availability represents the water resource of a cell.
                        // Thus, the outflow of SWB should be taken into account in inland sinks as well.

                        // monthly min and max river availability output start
                        //most certainly to be rationalised
                        if (options.outMinMaxRiverAvail) {
                            if (day_in_month == 1) {
                                G_monthlyMinRiverAvail[n][month] = transportedVolume;
                                G_monthlyMaxRiverAvail[n][month] = transportedVolume;
                            } else {
                                if (transportedVolume > G_monthlyMaxRiverAvail[n][month])
                                    G_monthlyMaxRiverAvail[n][month] = transportedVolume;
                                if (transportedVolume < G_monthlyMinRiverAvail[n][month])
                                    G_monthlyMinRiverAvail[n][month] = transportedVolume;
                            }
                        }
                        // monthly min and max river availability output end
                        if (options.outRiverPET == 1) G_monthlyRiverAreaFrac[n][month] += G_riverAreaFrac[n];
                        if (options.outRiverPET == 1)
                            G_monthlyRiverPET[n][month] += G_dailyRiverEvapo[n] * 1000000. /
                                                           (G_riverAreaFrac[n] / 100. * cellArea); //km3 > mm
                        //output for reservoir and storages start
                        if ((options.antNatOpt == 0) && (options.resOpt == 1)) {
                            G_monthlyResInflow[n][month] += G_daily_res_in[n]; // total inflow in month [m3/timeStepsPerDay] -> [m3/month]
                            G_monthlyResOutflow[n][month] += G_daily_res_out[n]; // total outflow in month [m3/timeStepsPerDay] -> [m3/month]

                            if (options.grid_store_TypeForStorages == 1)
                                G_monthlySurfStor[n][month] +=
                                               (G_locLakeStorage[n]
                                                + G_locWetlStorage[n]
                                                + G_gloLakeStorage[n]
                                                + G_gloWetlStorage[n]
                                                + G_riverStorage[n]
                                                + G_gloResStorage[n]);

                        }
                            // endif options.resOpt
                        else {
                            // formerly "else if (options.resOpt == 0)"
                            // mean value of total surface water storage in grid cell (lakes/wetlands/rivers) [km3] (converted to mm for optional output in a later step)
                            if (options.grid_store_TypeForStorages == 1)
                                G_monthlySurfStor[n][month] +=
                                               (G_locLakeStorage[n]
                                                + G_locWetlStorage[n]
                                                + G_gloLakeStorage[n]
                                                + G_gloWetlStorage[n]
                                                + G_riverStorage[n]);
                            // / (cellArea * (geo.G_contfreq[n] / 100.)) * 1000000. / timeStepsPerDay );
                        }

                        if (options.outGwrSwb == 1)
                            G_monthlyGwrSwb[n][month] += G_gwr_swb[n];

                        // if (options.outFswb == 1)
                        //     G_monthlyFswb[n][month] += G_fswb[n];
                        G_monthlyLocWetlExtent[n][month] += G_locwetlextent[n];
                        G_monthlyGloWetlExtent[n][month] += G_glowetlextent[n];


                    }
                    // endif options ((2 == options.grid_store) || (4 == options.grid_store) || (6 == options.grid_store))
                    // monthly output (possibly combined with daily output) endstart

                    // daily output option 31 start
                    if ((3 == options.grid_store) || (4 == options.grid_store)) {

                        if (options.outPrecDaily) {
                            // Only continental cells
                            if (geo.G_contcell[n]) {
                                G_daily31ConsistentPrecip[n][day_in_month - 1] = (
                                               (dailyWaterBalance.G_openWaterPrec[n] *
                                                G_landWaterExclGloLakAreaFrac[n] / geo.G_contfreq[n]) *
                                               ((geo.areaOfCellByArrayPos(n) * (geo.G_contfreq[n] / 100.)) / 1000000.) +
                                               G_GloLakePrecip[n] + G_ReservoirPrecip[n] + G_dailyRiverPrecip[n]);
                            }
                                // 100% water cells (of Caspian Sea) have now G_fwaterfreq + G_landfreq = 0. The water balance is calculated in the outlet cell of large lakes; therefore a calculation of e.g. precipitation for this cells leads to inconsistent quantification of global fluxes. A division through 0 has to be avoided now (alternative: define it as not land cells in ng).
                            else {
                                G_daily31ConsistentPrecip[n][day_in_month - 1] = 0.;
                            }
                        } // endif options.outPrecDaily

                        if (options.outCellRunoffDaily)
                            G_daily31CellRunoff[n][day_in_month - 1] = CellRunoff;

                        if (options.outCellSurfaceDaily) {
                            //  Only continental cells
                            if (geo.G_contcell[n]) {
                                G_daily31CellSurfaceRunoff[n][day_in_month - 1] =
                                               (CellRunoff - G_localGWRunoff[n])
                                               / (cellArea * (geo.G_contfreq[n] / 100.)) * 1000000.;
                            }
                                // 100% water cells (of Caspian Sea) have now G_fwaterfreq + G_landfreq = 0. The water balance is calculated in the outlet cell of large lakes; therefore a calculation of e.g. precipitation for this cells leads to inconsistent quantification of global fluxes. A division through 0 has to be avoided now (alternative: define it as not land cells in ng).
                            else {
                                G_daily31CellSurfaceRunoff[n][day_in_month - 1] = 0.;
                            }
                        } // endif options.outCellSurfaceDaily

                        if (options.outRiverAvailDaily) {
                            G_daily31RiverAvail[n][day_in_month - 1] = transportedVolume;  // FP20160926N001

                            // for those cells which are defined as inland sinks and do have considerable lake/wetland area (>'0'),
                            // river availability should be updated, otherwise river availability on basin scale could be '0'

                            // WG22b: The river availability represents the water resource of a cell (including inland sinks).
                            // Thus, the outflow of SWB should be taken into account in inland sinks as well.
                        } // endif options.outRiverAvailDaily

                        if (options.outRiverVeloDaily)
                            G_daily31Velocity[n][day_in_month - 1] = (riverVelocity /
                                                                      86.4); // velocity grid [m/s]

                        if (options.outGwrSwbDaily)
                            G_daily31GwrSwb[n][day_in_month -
                                               1] = G_gwr_swb[n]; // daily gwr below surface water bodies [mm]
                        if (options.outGWRunoffDaily)
                            // daily fraction of surface water bodies [mm]
                            G_daily31GwRunoff[n][day_in_month -
                                                 1] = groundwater_runoff_km3; // daily fraction of surface water bodies [mm]
                        if (options.outGwrunSurfrunDaily)
                            G_daily31GwrunSurfrun[n][day_in_month - 1] =
                                           groundwater_runoff_km3 + dailyLocalSurfaceRunoff; // in km3 day-1
                        if (options.outLandAreaFracDaily)
                            G_daily31LandAreaFrac[n][day_in_month - 1] = G_landAreaFrac[n];
                        if (options.outCellAETWCaDaily) // add NAg and actual use SW and later AET
                            G_daily31CellAETWCa[n][day_in_month - 1] =
                                           dailyActualUse + G_dailydailyNUg[n]; //already in km3 day-1
                    }
                    // endif option ((3 == options.grid_store) || (4 == options.grid_store))
                    // daily output option 31 end

                    // daily output option 365 start

                    if ((5 == options.grid_store) || (6 == options.grid_store)) {

                        if (options.outPrecDaily || ((options.scoutcPrecip) && (2 == options.day_store))) {
                            // FP 2015-06 // Only continental cells
                            if (geo.G_contcell[n]) {
                                G_daily365ConsistentPrecip[n][day - 1] = ((dailyWaterBalance.G_openWaterPrec[n] *
                                                                           G_landWaterExclGloLakAreaFrac[n] /
                                                                           geo.G_contfreq[n]) *
                                                                          ((geo.areaOfCellByArrayPos(n) *
                                                                            (geo.G_contfreq[n] / 100.)) / 1000000.) +
                                                                          G_GloLakePrecip[n] + G_ReservoirPrecip[n] +
                                                                          G_dailyRiverPrecip[n]);

                            }
                                // 100% water cells (of Caspian Sea) have now G_fwaterfreq + G_landfreq = 0. The water balance is calculated in the outlet cell of large lakes; therefore a calculation of e.g. precipitation for this cells leads to inconsistent quantification of global fluxes. A division through 0 has to be avoided now (alternative: define it as not land cells in ng).
                            else {
                                G_daily365ConsistentPrecip[n][day - 1] = 0.;
                            }
                        }

                        if ((options.outGwrSwbDaily) || ((options.scoutGwrSwb) && (2 == options.day_store)))
                            G_daily365GwrSwb[n][day - 1] = G_gwr_swb[n];
                        //if ((options.outFswbDaily)||((options.scoutFswb)&&(2==options.day_store)))
                        //    G_daily365Fswb[n][day - 1] = G_fswb[n];
                        if ((options.outGWRunoffDaily) || ((options.scoutGwRunoff) && (2 == options.day_store)))
                            //HMS 2017-02-03 the following line was included but makes no sense, or?
                            //G_daily365GwRunoff[n][day - 1] = G_groundwaterStorage[n];
                            G_daily365GwRunoff[n][day - 1] = groundwater_runoff_km3;
                        if (options.outGwrunSurfrunDaily)
                            G_daily365GwrunSurfrun[n][day - 1] = groundwater_runoff_km3 + dailyLocalSurfaceRunoff;
                        if ((options.outLandAreaFracDaily) || ((options.scoutLandAreaFrac) && (2 == options.day_store)))
                            G_daily365LandAreaFrac[n][day - 1] = G_landAreaFrac[n];

                        if (options.outCellRunoffDaily ||
                            (((options.scoutCellRunoffkm3) || (options.scoutCellRunoffmm)) &&
                             (2 == options.day_store))) {
                            G_daily365CellRunoff[n][day - 1] = CellRunoff;  // FP20160926N001
                            // FP 2015-06 // Only continental cells
                            if (geo.G_contcell[n]) {
                                G_daily365CellRunoff_mm[n][day - 1] =
                                               CellRunoff / (cellArea * (geo.G_contfreq[n] / 100.)) *
                                               1000000.;  // FP20160926N001
                            }
                                // 100% water cells (of Caspian Sea) have now G_fwaterfreq + G_landfreq = 0. The water balance is calculated in the outlet cell of large lakes; therefore a calculation of e.g. precipitation for this cells leads to inconsistent quantification of global fluxes. A division through 0 has to be avoided now (alternative: define it as not land cells in ng).
                            else {
                                G_daily365CellRunoff_mm[n][day - 1] = 0.;
                            }
                        }
                        // endif options

                        if (options.outCellSurfaceDaily ||
                            ((options.scoutCellSRunoff) && (2 == options.day_store))) { // FP20160915N002
                            // Only continental cells
                            // 100% water cells (of Caspian Sea) have now G_fwaterfreq + G_landfreq = 0. The water balance is calculated in the outlet cell of large lakes; therefore a calculation of e.g. precipitation for this cells leads to inconsistent quantification of global fluxes. A division through 0 has to be avoided now (alternative: define it as not land cells in ng).
                            if (geo.G_contcell[n]) {
                                G_daily365CellSurfaceRunoff[n][day - 1] = (CellRunoff - G_localGWRunoff[n]) /
                                                                          (cellArea * (geo.G_contfreq[n] / 100.)) *
                                                                          1000000.;  // FP20160926N001
                            }
                                // 100% water cells (of Caspian Sea) have now G_fwaterfreq + G_landfreq = 0. The water balance is calculated in the outlet cell of large lakes; therefore a calculation of e.g. precipitation for this cells leads to inconsistent quantification of global fluxes. A division through 0 has to be avoided now (alternative: define it as not land cells in ng).
                            else {
                                G_daily365CellSurfaceRunoff[n][day - 1] = 0.;
                            }
                        }
                        // endif options

                        if (options.outRiverAvailDaily || ((options.scoutQ) && (2 == options.day_store))) {
                            G_daily365RiverAvail[n][day - 1] = transportedVolume;  // FP20160926N001
                            // for those cells which are defined as inland sinks and do have considerable lake/wetland area (>'0'),
                            // river availability should be updated, otherwise river availability on basin scale could be '0'
                            // WG22b: The river availability represents the water resource of a cell (including inland sinks).
                            // Thus, the outflow of SWB should be taken into account in inland sinks as well.
                            /*if (G_LDD[n] == -1) {
                                if ((G_glo_wetland[n] > 0) || (G_lake_area[n] > 0))
                                    G_daily365RiverAvail[n][day - 1] += (inflowUpstream - transportedVolume);
                            }*/
                        }
                        if (options.outCellAETWCaDaily)
                            G_daily365CellAETWCa[n][day - 1] =
                                           dailyActualUse + G_dailydailyNUg[n]; //already in km3 day-1
                        // endif options

                        if ((options.outRiverVeloDaily) || ((options.scoutFlowVelo) && (2 == options.day_store)))
                            G_daily365Velocity[n][day - 1] = (riverVelocity /
                                                              86.4); // velocity grid [m/s]  // FP20160926N001

                    }
                    // endif option
                    // daily output option 365 end

                    for (b = 0; b <= nSpecBasins - 1; b++) {
                        if (n == cbasin.cellNum[b] - 1) {
                            dailyRiverDischarge[b][day - 1] += transportedVolume;
                            break;
                        }
                    }
                    // endfor loop over nSpecBasins

                }
                // for next loop in timeStepsPerDay

                G_riverInflow[n] = 0.;

            }
            // endif (geo.G_contcell[n]) [only execute for continental cells]

        } // end of loop for gridcells (routing)

    } // end of loop for routing (loop for routing time steps)


    // Calculate total actual Evapotranspiration (land and open water areas); consistent with corrected runoff (2.1f)
    // Store storage values and end of timestep
    if (options.grid_store > 1) {

        double dailyTotalAET;
        double dailyTotalAET_km3;
        double dailyOpenWaterEvap;
        double dailyOpenWaterEvap_km3;

        for (int n = 0; n < ng; n++) {
            if (0 != G_toBeCalculated[n]) {
                dailyTotalAET = 0.;
                dailyTotalAET_km3 = 0.;
                dailyOpenWaterEvap = 0.;
                dailyOpenWaterEvap_km3 = 0.;
                cellArea = geo.areaOfCellByArrayPos(n);


                // first store only land area evapotranspiration
                dailyTotalAET = dailyWaterBalance.getDailyCellLandAET(n);

                // then sum up land area evapotranspiration and evaporation
                // from open water areas
                // all areas in [km2]
                // resulting data unit: mm * km2


                // G_dailyLocLakeEvapo[n] etc. = openWaterPET * area reduction factor -> SWB area already reduced;
                // therefore, G_loc_lake[n] etc. is used in the following equation:
                dailyOpenWaterEvap_km3 =    // in km3
                               (G_dailyLocLakeEvapo[n] / 1000000. * cellArea * G_loc_lake[n] / 100.)
                               + (G_dailyLocWetlEvapo[n] / 1000000. * cellArea * G_loc_wetland[n] / 100.)
                               + (G_dailyGloLakeEvapo[n] / 1000000. * G_lake_area[n])
                               + (G_dailyResEvapo[n] / 1000000. * G_reservoir_area[n])
                               + (G_dailyGloWetlEvapo[n] / 1000000. * cellArea * G_glo_wetland[n] / 100.);

                if (options.riverEvapoOpt == 1) {

                    // dailyOpenWaterEvap now dailyOpenWaterEvap_km3 in km3
                    // POSSIBLE ERROR data unit of dailyOpenWaterEvap: mm * km2 does not correspond to km3 of G_dailyRiverEvapo[n];
                    //G_dailyRiverEvapo[n] wird oben in km3 berechnet und mit G_riverAreaFrac[n] multipliziert
                    //dailyOpenWaterEvap += G_dailyRiverEvapo[n] * cellArea * G_riverAreaFrac[n] / 100.;
                    dailyOpenWaterEvap_km3 += G_dailyRiverEvapo[n];        // G_dailyRiverEvapo[n] in km3
                }

                dailyTotalAET_km3 = (dailyTotalAET / 1000000. * cellArea * G_landAreaFrac[n] /
                                     100.)        // dailyTotalAET: mm -> km3
                                    + dailyOpenWaterEvap_km3;

                // monthly output start
                if ((2 == options.grid_store) || (4 == options.grid_store) || (6 == options.grid_store)) {
                    // FP 2015 // Only continental cells
                    if (geo.G_contcell[n]) {

                        // AET in mm
                        G_monthlyCellAET[n][month] += (dailyTotalAET_km3 * 1000000. /
                                                       ((cellArea * geo.G_contfreq[n] / 100.)));    // [mm]
                        G_monthlyCellAETWCa[n][month] += dailyTotalAET_km3;    // [km3]

                        // open water evaporation from surface storages (without rivers)
                        G_monthlyOpenWaterEvap[n][month] += (dailyOpenWaterEvap_km3 * 1000000. /
                                                             ((cellArea * geo.G_contfreq[n] / 100.)));    // [mm]

                    }
                        // 100% water cells (of Caspian Sea) have now G_fwaterfreq + G_landfreq = 0. The water balance is calculated in the outlet cell of large lakes; therefore a calculation of e.g. precipitation for this cells leads to inconsistent quantification of global fluxes. A division through 0 has to be avoided now (alternative: define it as not land cells in ng).
                    else {
                        G_monthlyCellAET[n][month] = 0.;
                        G_monthlyCellAETWCa[n][month] = 0.;
                        G_monthlyOpenWaterEvap[n][month] = 0.;
                    }
                }
                // monthly output end

                // daily output option 31 start
                if ((3 == options.grid_store) || (4 == options.grid_store)) {
                    if (options.outCellAETDaily) {
                        // FP 2015-06 // Only continental cells
                        if (geo.G_contcell[n]) {
                            // AET in mm
                            G_daily31CellAET[n][day_in_month - 1] = (dailyTotalAET_km3 * 1000000. /
                                                                     ((cellArea * geo.G_contfreq[n] / 100.)));
                        }
                            // 100% water cells (of Caspian Sea) have now G_fwaterfreq + G_landfreq = 0. The water balance is calculated in the outlet cell of large lakes; therefore a calculation of e.g. precipitation for this cells leads to inconsistent quantification of global fluxes. A division through 0 has to be avoided now (alternative: define it as not land cells in ng).
                        else {
                            G_daily31CellAET[n][day_in_month - 1] = 0.;
                        }
                    } // end if options.outCellAETDaily
                    if (options.outCellAETWCaDaily) {
                        if (geo.G_contcell[n]) {
                            // AETWCa in km3

                            G_daily31CellAETWCa[n][day_in_month - 1] += dailyTotalAET_km3;
                        }
                            // 100% water cells (of Caspian Sea) have now G_fwaterfreq + G_landfreq = 0. The water balance is calculated in the outlet cell of large lakes; therefore a calculation of e.g. precipitation for this cells leads to inconsistent quantification of global fluxes. A division through 0 has to be avoided now (alternative: define it as not land cells in ng).
                        else {
                            G_daily31CellAETWCa[n][day_in_month - 1] = 0.;
                        }
                    } // end if options.outCellAETWCaDaily
                }
                // daily output option 31 end

                // daily output option 365 start
                if ((5 == options.grid_store) || (6 == options.grid_store)) {

                    if ((options.outCellAETDaily) || ((options.scoutCellAET) && (2 == options.day_store))) {

                        // FP 2015-06 // Only continental cells
                        if (geo.G_contcell[n]) {
                            // AET in mm
                            G_daily365CellAET[n][day - 1] = (dailyTotalAET_km3 * 1000000. /
                                                             ((cellArea * geo.G_contfreq[n] / 100.)));
                        }
                            // 100% water cells (of Caspian Sea) have now G_fwaterfreq + G_landfreq = 0. The water balance is calculated in the outlet cell of large lakes; therefore a calculation of e.g. precipitation for this cells leads to inconsistent quantification of global fluxes. A division through 0 has to be avoided now (alternative: define it as not land cells in ng).
                        else {
                            G_daily365CellAET[n][day - 1] = 0.;
                        }

                    } // endif options.outCellAETDaily
                    if (options.outCellAETWCaDaily) {

                        if (geo.G_contcell[n]) {
                            // AET in km3
                            G_daily365CellAETWCa[n][day - 1] += dailyTotalAET_km3;
                        }
                            // 100% water cells (of Caspian Sea) have now G_fwaterfreq + G_landfreq = 0. The water balance is calculated in the outlet cell of large lakes; therefore a calculation of e.g. precipitation for this cells leads to inconsistent quantification of global fluxes. A division through 0 has to be avoided now (alternative: define it as not land cells in ng).
                        else {
                            G_daily365CellAETWCa[n][day - 1] = 0.;
                        }

                    } // endif options.outCellAETDaily
                }
                // daily output option 365 end

            }
            // endif G_toBeCalculated[n]
        }
        // end loop over grid cells
    }
    // end of: if (options.grid_store > 1)


    //calculate total monthly storage per grid cell
    // here: from lateral water balance (second part)
    // in km3
    if ((2 == options.grid_store) || (4 == options.grid_store) || (6 == options.grid_store)) {
        for (int n = 0; n < ng; n++) {
            if (0 != G_toBeCalculated[n]) {

                //output for reservoir and storages start
                if (options.resOpt == 1) {
                    if (G_reservoir_area[n] > 0.) {
                        G_monthlyResInflow[n][month] += G_daily_res_in[n]; // total inflow in month [km3/day] -> [km3/month]
                        G_monthlyResOutflow[n][month] += G_daily_res_out[n]; // total outflow in month [km3/day] -> [km3/month]
                    } else {
                        G_monthlyResInflow[n][month] = 0.; // no reservoir area --> no Inflow (anyhow, output is not yet implemented)
                        G_monthlyResOutflow[n][month] = 0.; // no reservoir area --> no Outflow
                    }
                }
                // endif options.resOpt

                // output if
                // G_monthlyResStorage, G_monthlyLocLakeStorage, G_monthlyGloLakeStorage, G_monthlyRiverStorage,
                // G_monthlySurfStor, G_monthlyLocWetlStorage, G_monthlyGloWetlStorage are values at end of month
                if (options.grid_store_TypeForStorages == 0) {
                    if (options.resOpt == 1) {

                        G_monthlySurfStor[n][month] =
                                       (G_locLakeStorage[n]
                                        + G_locWetlStorage[n]
                                        + G_gloLakeStorage[n]
                                        + G_gloWetlStorage[n]
                                        + G_riverStorage[n]
                                        + G_gloResStorage[n]);
                        // / (cellArea * (geo.G_contfreq[n] / 100.)) * 1000000.);
                    } // endif options.resOpt
                    else {
                        // total surface water storage in grid cell (lakes/wetlands/rivers) // [km3]
                        G_monthlySurfStor[n][month] =
                                       (G_locLakeStorage[n]
                                        + G_locWetlStorage[n]
                                        + G_gloLakeStorage[n]
                                        + G_gloWetlStorage[n]
                                        + G_riverStorage[n]);
                        // / (cellArea * (geo.G_contfreq[n] / 100.)) * 1000000. );
                    }
                }
                // endif (options.grid_store_TypeForStorages == 0)

                //FP20161213N001 Probably double-counting here, as same procedure already present in loop over time steps
                // if mean monthly surface storage value is of interest
                if (options.grid_store_TypeForStorages == 1) {
                    if ((options.antNatOpt == 0) && (options.resOpt == 1)) {
                        G_monthlySurfStor[n][month] +=
                                       (G_locLakeStorage[n]
                                        + G_locWetlStorage[n]
                                        + G_gloLakeStorage[n]
                                        + G_gloWetlStorage[n]
                                        + G_riverStorage[n]
                                        + G_gloResStorage[n]);


                    } else { //if (options.resOpt == 0){ HMS 2017-06-08 obsolete as either it is ant + resOpt == 1 or any other.
                        // mean value of total surface water storage in grid cell (lakes/wetlands/rivers) [km3] (converted to mm for optional output in a later step)
                            G_monthlySurfStor[n][month] +=
                                           (G_locLakeStorage[n]
                                            + G_locWetlStorage[n]
                                            + G_gloLakeStorage[n]
                                            + G_gloWetlStorage[n]
                                            + G_riverStorage[n]);
                        // / (cellArea * (geo.G_contfreq[n] / 100.)) * 1000000. / timeStepsPerDay );
                    }
                }
                // endif (options.grid_store_TypeForStorages == 1)
//FP20161213N001 Item 2 start
                //output for natural lakes and wetlands
                if (options.grid_store_TypeForStorages == 0) {
                    G_monthlyLocLakeStorage[n][month] = G_locLakeStorage[n];
                    G_monthlyGloLakeStorage[n][month] = G_gloLakeStorage[n];
                    G_monthlyRiverStorage[n][month] = G_riverStorage[n];
                    G_monthlyLocWetlStorage[n][month] = G_locWetlStorage[n];
                    G_monthlyGloWetlStorage[n][month] = G_gloWetlStorage[n];
                    G_monthlyGwStorage[n][month] = G_groundwaterStorage[n];

                    if ((options.antNatOpt == 0) && (options.resOpt == 1))
                        G_monthlyResStorage[n][month] = G_gloResStorage[n];

                } else if (options.grid_store_TypeForStorages == 1) {
                    G_monthlyLocLakeStorage[n][month] += G_locLakeStorage[n];
                    G_monthlyGloLakeStorage[n][month] += G_gloLakeStorage[n];
                    G_monthlyRiverStorage[n][month] += G_riverStorage[n];
                    G_monthlyLocWetlStorage[n][month] += G_locWetlStorage[n];
                    G_monthlyGloWetlStorage[n][month] += G_gloWetlStorage[n];
                    G_monthlyGwStorage[n][month] += G_groundwaterStorage[n];

                    if ((options.antNatOpt == 0) && (options.resOpt == 1))
                        G_monthlyResStorage[n][month] += G_gloResStorage[n];
                }

//FP20161213N001 Item 2 end
            }
            // end if G_toBeCalculated[n]
        }
        // end loop over grid cells
    }
    // endif option


    if ((3 == options.grid_store) || (4 == options.grid_store)) { // .31 gw storage output calculation in km3
        // calculation G_monthlyGwStorage formerly in daily.cpp, has to be calculated after NUg is subtracted
        for (int n = 0; n < ng; n++) {
            if (0 != G_toBeCalculated[n]) {
                if (options.outGWStorageDaily || options.outSingleStoragesDaily)
                    G_daily31GwStor[n][day_in_month - 1] += G_groundwaterStorage[n];
            } // end if G_toBeCalculated[n]
        }
        // end loop over grid cells
    }
    // endif option


    if ((5 == options.grid_store) || (6 == options.grid_store)) {  // .365 gw storage output - calculation in km3
        for (int n = 0; n < ng; n++) {
            if (0 != G_toBeCalculated[n]) {
                if ((options.outGWStorageDaily || options.outSingleStoragesDaily) ||
                    ((options.scoutGwStor) && (2 == options.day_store)))
                    G_daily365GwStor[n][day - 1] = G_groundwaterStorage[n];
            }
            // end if G_toBeCalculated[n]
        }
        // end loop over grid cells
    }
    // endif option


    // daily 31 storage calculation
    if (((3 == options.grid_store) || (4 == options.grid_store)) &&
        ((options.outTotalWaterInStoragesDaily_km3) || (options.outTotalWaterInStoragesDaily_mm))) {
        for (int n = 0; n < ng; n++) {
            if (0 != G_toBeCalculated[n]) {

                if (options.outSingleStoragesDaily) {
                    G_daily31LocLakeStor[n][day_in_month - 1] = G_locLakeStorage[n];
                    G_daily31LocWetlStor[n][day_in_month - 1] = G_locWetlStorage[n];
                    G_daily31GloLakeStor[n][day_in_month - 1] = G_gloLakeStorage[n];
                    G_daily31GloWetlStor[n][day_in_month - 1] = G_gloWetlStorage[n];
                    G_daily31RiverStor[n][day_in_month - 1] = G_riverStorage[n];
                    if (options.resOpt == 1)
                        G_daily31ResStor[n][day_in_month - 1] = G_gloResStorage[n];
                }
                // endif options.outSurfStorDaily

                // daily total water storage output (calculate here after subtracting water use from groundwater)
                // calculate daily 31 storage per grid cell
                // here: get from vertical water balance (first part)
                // in km³
                double storage_daily31 = 0.;
                storage_daily31 = dailyWaterBalance.G_daily31CanopyStorage[n][day_in_month - 1] // [km3]
                                  + dailyWaterBalance.G_daily31SnowStorage[n][day_in_month - 1]
                                  + dailyWaterBalance.G_daily31SoilStorage[n][day_in_month - 1];

                // here: add the storages from routing (and later groundwater) (second part)
                G_daily31TotalWaterInStorages_km3[n][day_in_month - 1] =
                               storage_daily31
                               + G_locLakeStorage[n]
                               + G_locWetlStorage[n]
                               + G_gloLakeStorage[n]
                               + G_gloWetlStorage[n]
                               + G_riverStorage[n]
                               + G_groundwaterStorage[n];
                if (options.resOpt) {
                    G_daily31TotalWaterInStorages_km3[n][day_in_month - 1] += G_gloResStorage[n];
                }

                // FP 2015-06 // Only continental cells
                if (geo.G_contcell[n]) {
                    // convert total storage from km3 to mm
                    G_daily31TotalWaterInStorages_mm[n][day_in_month - 1] =
                                   G_daily31TotalWaterInStorages_km3[n][day_in_month - 1]
                                   / (geo.areaOfCellByArrayPos(n) * (geo.G_contfreq[n] / 100.)) * 1000000.;
                }
                    // 100% water cells (of Caspian Sea) have now G_fwaterfreq + G_landfreq = 0. The water balance is calculated in the outlet cell of large lakes; therefore a calculation of e.g. precipitation for this cells leads to inconsistent quantification of global fluxes. A division through 0 has to be avoided now (alternative: define it as not land cells in ng).
                else {
                    G_daily31TotalWaterInStorages_mm[n][day_in_month - 1] = 0.;
                }

                // daily surface water storage km3
                if (options.outSurfStorDaily) {
                    if (options.resOpt == 1)
                        G_daily31SurfStor[n][day_in_month - 1] +=
                                       (G_locLakeStorage[n]
                                        + G_locWetlStorage[n]
                                        + G_gloLakeStorage[n]
                                        + G_gloWetlStorage[n]
                                        + G_riverStorage[n]
                                        + G_gloResStorage[n]);
                        // 	/ (cellArea * (geo.G_contfreq[n] / 100.)) * 1000000. / timeStepsPerDay );
                    else
                        G_daily31SurfStor[n][day_in_month - 1] += // [km3]
                                       (G_locLakeStorage[n]
                                        + G_locWetlStorage[n]
                                        + G_gloLakeStorage[n]
                                        + G_gloWetlStorage[n]
                                        + G_riverStorage[n]);
                    //	/ (cellArea * (geo.G_contfreq[n] / 100.)) * 1000000. / timeStepsPerDay );
                }
                // endif options.outSurfStorDaily

            }
            // endif G_toBeCalculated[n]

        }
        // end loop over grid cells

    }
    // endif option


    if (((5 == options.grid_store) || (6 == options.grid_store)) &&
        ((options.outTotalWaterInStoragesDaily_km3) || (options.outTotalWaterInStoragesDaily_mm) ||
         (((options.scoutTWSkm3) || (options.scoutTWSmm)) && (2 == options.day_store)))) {
        for (int n = 0; n < ng; n++) {
            if (0 != G_toBeCalculated[n]) {

                if ((options.outSingleStoragesDaily) || ((options.scoutLocLake) && (2 == options.day_store)))
                    G_daily365LocLakeStor[n][day - 1] = G_locLakeStorage[n];
                if ((options.outSingleStoragesDaily) || ((options.scoutLocWet) && (2 == options.day_store)))
                    G_daily365LocWetlStor[n][day - 1] = G_locWetlStorage[n];
                if ((options.outSingleStoragesDaily) || ((options.scoutGloLake) && (2 == options.day_store)))
                    G_daily365GloLakeStor[n][day - 1] = G_gloLakeStorage[n];
                if ((options.outSingleStoragesDaily) || ((options.scoutGloWet) && (2 == options.day_store)))
                    G_daily365GloWetlStor[n][day - 1] = G_gloWetlStorage[n];
                if ((options.outSingleStoragesDaily) || ((options.scoutRiver) && (2 == options.day_store)))
                    G_daily365RiverStor[n][day - 1] = G_riverStorage[n];
                if ((options.outSingleStoragesDaily) || ((options.scoutReservoir) && (2 == options.day_store))) {
                    if (options.resOpt == 1)
                        G_daily365ResStor[n][day - 1] = G_gloResStorage[n];
                }

                // daily total water storage output (calculate here after subtracting water use from groundwater)
                // calculate daily 365 storage per grid cell
                // here: get from vertical water balance (first part)
                // in km³

                double storage_daily365 = 0.;
                storage_daily365 = dailyWaterBalance.G_daily365CanopyStorage[n][day - 1] // [km3]
                                   + dailyWaterBalance.G_daily365SnowStorage[n][day - 1]
                                   + dailyWaterBalance.G_daily365SoilStorage[n][day - 1];

                // here: add the storages from routing (and later groundwater) (second part)
                G_daily365TotalWaterInStorages_km3[n][day - 1] =
                               storage_daily365
                               + G_locLakeStorage[n]
                               + G_locWetlStorage[n]
                               + G_gloLakeStorage[n]
                               + G_gloWetlStorage[n]
                               + G_riverStorage[n]
                               + G_groundwaterStorage[n];
                if (options.resOpt || (((options.scoutTWSkm3) || (options.scoutTWSmm)) && (2 == options.day_store)))
                    G_daily365TotalWaterInStorages_km3[n][day - 1] += G_gloResStorage[n];


                if (options.outTotalWaterInStoragesDaily_mm || options.scoutTWSmm) {
                    // FP 2015-06 // Only continental cells
                    if (geo.G_contcell[n]) {
                        // convert total storages from km3 to mm
                        G_daily365TotalWaterInStorages_mm[n][day - 1] =
                                       G_daily365TotalWaterInStorages_km3[n][day - 1]
                                       / (geo.areaOfCellByArrayPos(n) * (geo.G_contfreq[n] / 100.)) * 1000000.;
                    }
                        // 100% water cells (of Caspian Sea) have now G_fwaterfreq + G_landfreq = 0. The water balance is calculated in the outlet cell of large lakes; therefore a calculation of e.g. precipitation for this cells leads to inconsistent quantification of global fluxes. A division through 0 has to be avoided now (alternative: define it as not land cells in ng).
                    else {
                        G_daily365TotalWaterInStorages_mm[n][day - 1] = 0.;
                    }
                }
                // endif options

                // daily surface water storage km3
                if ((options.outSurfStorDaily) || ((options.scoutSurfaceStor) && (2 == options.day_store))) {
                    if (options.resOpt == 1)
                        G_daily365SurfStor[n][day - 1] =
                                       (G_locLakeStorage[n]
                                        + G_locWetlStorage[n]
                                        + G_gloLakeStorage[n]
                                        + G_gloWetlStorage[n]
                                        + G_riverStorage[n]
                                        + G_gloResStorage[n]);
                        // / (cellArea * (geo.G_contfreq[n] / 100.)) * 1000000. / timeStepsPerDay );

                    else
                        G_daily365SurfStor[n][day - 1] = // [km3]
                                       (G_locLakeStorage[n]
                                        + G_locWetlStorage[n]
                                        + G_gloLakeStorage[n]
                                        + G_gloWetlStorage[n]
                                        + G_riverStorage[n]);
                    // / (cellArea * (geo.G_contfreq[n] / 100.)) * 1000000. / timeStepsPerDay );
                }
                // endif options

            }
            // endif (0 != G_toBeCalculated[n])
        }
        // end loop over grid cells
    }
    // endif option

    //HMS startend TWS output
    if (options.outTotalWaterInStoragesStartEndDaily_km3) {
        for (int n = 0; n < ng; n++) {
            if (0 != G_toBeCalculated[n]) {
                if (day == 1) { // calculate for first day in year
                    G_startendTotalWaterInStorages_km3[n][0] =
                                   dailyWaterBalance.G_canopyWaterContent[n] * geo.areaOfCellByArrayPos(n) / 1000000.0 *
                                   G_landAreaFrac[n] / 100.0 // [km3]
                                   + dailyWaterBalance.G_snow[n] * geo.areaOfCellByArrayPos(n) / 1000000.0 *
                                     G_landAreaFrac[n] / 100.0
                                   + dailyWaterBalance.G_soilWaterContent[n] * geo.areaOfCellByArrayPos(n) / 1000000.0 *
                                     G_landAreaFrac[n] / 100.0
                                   + G_locLakeStorage[n]
                                   + G_locWetlStorage[n]
                                   + G_gloLakeStorage[n]
                                   + G_gloWetlStorage[n]
                                   + G_riverStorage[n]
                                   + G_groundwaterStorage[n];
                    if (options.resOpt)
                        G_startendTotalWaterInStorages_km3[n][0] += G_gloResStorage[n];
                }
                if (day == 365) { // calculate for last day in year
                    G_startendTotalWaterInStorages_km3[n][1] =
                                   dailyWaterBalance.G_canopyWaterContent[n] * geo.areaOfCellByArrayPos(n) / 1000000.0 *
                                   G_landAreaFrac[n] / 100.0 // [km3]
                                   + dailyWaterBalance.G_snow[n] * geo.areaOfCellByArrayPos(n) / 1000000.0 *
                                     G_landAreaFrac[n] / 100.0
                                   + dailyWaterBalance.G_soilWaterContent[n] * geo.areaOfCellByArrayPos(n) / 1000000.0 *
                                     G_landAreaFrac[n] / 100.0
                                   + G_locLakeStorage[n]
                                   + G_locWetlStorage[n]
                                   + G_gloLakeStorage[n]
                                   + G_gloWetlStorage[n]
                                   + G_riverStorage[n]
                                   + G_groundwaterStorage[n];
                    if (options.resOpt)
                        G_startendTotalWaterInStorages_km3[n][1] += G_gloResStorage[n];
                }
            }
        }
    }


    if (options.day_store == 1) {
        for (short b = 0; b < nSpecBasins; b++) {
            dailyLocLakeStorage[b][day - 1] = getActualStorageRatio(cbasin.cellNum[b],
                                                                    1);  // ATTENTION: method using only cell-specific calibration parameter P_LAK_D  // FP
            dailyLocWetlStorage[b][day - 1] = getActualStorageRatio(cbasin.cellNum[b], 2);
            dailyGloLakeStorage[b][day - 1] = getActualStorageRatio(cbasin.cellNum[b],
                                                                    3);  // ATTENTION: method using only cell-specific calibration parameter P_LAK_D  // FP
            dailyGloWetlStorage[b][day - 1] = getActualStorageRatio(cbasin.cellNum[b], 4);
            if (options.resOpt == 1)
                dailyResStorage[b][day - 1] = getActualStorageRatio(cbasin.cellNum[b], 5);
        }
    }
    // endif option (options.day_store == 1)


    // Computation of SWB and land area fractions for the next time step
    for (int n = 0; n < ng; n++) {

        //if ((1 == options.fractionalRoutingOpt) || (options.aridareaOpt == 1) || (options.outFswb == 1)){

        // a second cell is allowed to satisfy the remaining demand

        // HMS 2013-11-22 multiply actual size of each swb with specific area in % and divide through continental area in % to get actual fraction of surface water bodies for the next time step
        // G_lake_area and G_reservoir_area are in km2 and have therefore considered differently to get the fraction of the continental area of those
        // to avoid numerical problems (sometimes the area in % gets nan), each swb type is summed up seperately.
        // for arid and for humid cells


        // FP20161018N002 Reservoir operation start years
        // Initial values are now assigned in routing.init
        if ((G_loc_lake[n] > 0.) && (G_locLakeAreaReductionFactor[n] > 0.)) {
            G_fLocLake[n] = (G_locLakeAreaReductionFactor[n] * G_loc_lake[n] / 100.);
        } else {
            G_fLocLake[n] = 0.; // to prevent that a single fraction values of a previous grid cell is used (if the conditions of the if clause is not true), the fractions will be 0 if the if clause is not true.
        }

        if ((G_loc_wetland[n] > 0.) && (G_locWetlAreaReductionFactor[n] > 0.)) {
            G_fLocWet[n] = (G_locWetlAreaReductionFactor[n] * G_loc_wetland[n] / 100.);
        } else {
            G_fLocWet[n] = 0.;
        }

        //if (((double)G_lake_area[n] > 0.) && (gloLakeEvapoReductionFactor > 0.))
        //    G_fGloLakeOutflowCell[n] = (gloLakeEvapoReductionFactor * ((double)G_lake_area[n] / (geo.areaOfCellByArrayPos(n) ))); // can reach values > 1 due to relocalization into outflow cell.
        //else
        //    G_fGloLakeOutflowCell[n] = 0.;

        //if (((double)G_reservoir_area[n] > 0.) && (gloResEvapoReductionFactor > 0.))
        //    G_fGloResOutflowCell[n] = (gloResEvapoReductionFactor * ((double)G_reservoir_area[n] / (geo.areaOfCellByArrayPos(n) ))); // can reach values > 1 due to relocalization into outflow cell.
        //else
        //    G_fGloResOutflowCell[n] = 0.;

        if ((G_glo_wetland[n] > 0.) && (G_gloWetlAreaReductionFactor[n] > 0.)) {
            G_fGloWet[n] = (G_gloWetlAreaReductionFactor[n] * G_glo_wetland[n] / 100.);
        } else {
            G_fGloWet[n] = 0.;
        }

        // FP20161018N002 Reservoir operation start years
        // calculate changes in G_fswb (initialization is already done in routing.init)

        // Old fraction from previously calculated time step
        G_fswbLandAreaFrac[n] = G_fswbLandAreaFracNextTimestep[n];  // In first time step identical values from routing.init
        // New current fraction
        G_fswbLandAreaFracNextTimestep[n] = G_fLocLake[n] + G_fLocWet[n] + G_fGloWet[n];
        G_fswbLandAreaFrac_ChangePct[n] = G_fswbLandAreaFracNextTimestep[n] * 100. - G_fswbLandAreaFrac[n] *
                                                                                     100.; // G_fwsb is in fraction, G_landAreaFrac (where G_fswbChange is calculated) is in %


        // adapt changes in landareafrac due to changes in river area fraction
        if (options.riverEvapoOpt == 1) {

            G_maxRiverAreaFrac[n] = geo.G_contfreq[n] / 100. -
                                    G_fGloLake[n]; //river area fraction can not be larger than "available" land area fraction

            // FP comment: inland sink?? and 100% glo lake
            if (((G_lake_area[n] > 0.) || (G_reservoir_area[n] > 0.)) &&
                (G_fGloLake[n] == 1.)) { // for outlet cells containing 100% global lakes, no river evapo at all
                G_riverAreaFracNextTimestep_Frac[n] = 0.;
                G_riverAreaFrac_ChangeFrac[n] = 0.;
                G_riverAreaReductionFactor[n] = 0.;
                G_dailyRiverEvapo[n] = 0.;
            } else { // no inland sink and no 100% glo lake
                //if (G_riverAreaFracChange[n] <= G_maxRiverAreaFrac[n]) //river fraction can reach max value
                if (G_riverAreaFracNextTimestep_Frac[n] <= G_maxRiverAreaFrac[n]) {//river fraction can reach max value

                    if (G_fswbLandAreaFracNextTimestep[n] >
                        (G_maxRiverAreaFrac[n] - G_riverAreaFracNextTimestep_Frac[n])) {
                        //if (G_riverAreaFracChange[n] > (G_maxRiverAreaFrac[n] - G_fswbLandAreaFracNextTimestep[n])) // not enough land fraction is available for river fraction, thus reduce swb fraction

                        fswbFracCorr = (G_maxRiverAreaFrac[n] - G_riverAreaFracNextTimestep_Frac[n]) /
                                       G_fswbLandAreaFracNextTimestep[n];
                        //fswbFracCorr = 1 - (G_fswbLandAreaFracNextTimestep[n] + G_riverAreaFracChange[n]) / G_maxRiverAreaFrac[n];
                        G_locLakeAreaReductionFactor[n] *= fswbFracCorr;
                        G_locWetlAreaReductionFactor[n] *= fswbFracCorr;
                        G_gloWetlAreaReductionFactor[n] *= fswbFracCorr;

                        //fLocLak = G_fLocLake[n] + G_fLocLake[n] * fswbFracCorr;
                        //fLocWet = G_fLocWet[n] + G_fLocWet[n] * fswbFracCorr;
                        //fGloWet = G_fGloWet[n] + G_fGloWet[n] * fswbFracCorr;
                        if ((G_fLocLake[n] > 0.) && (G_locLakeAreaReductionFactor[n] > 0.)) {
                            //G_locLakeAreaReductionFactor[n] = G_locLakeAreaReductionFactor[n] * G_fLocLake[n] / ((double)G_loc_lake[n]/100.); // adapt areareductionfactor
                            G_fLocLake[n] = (G_locLakeAreaReductionFactor[n] * G_loc_lake[n] /
                                             100.); // and finalize swb fraction
                        } else {
                            G_locLakeAreaReductionFactor[n] = 0.;
                            G_fLocLake[n] = 0.;
                        }

                        if ((G_fLocWet[n] > 0.) && (G_locWetlAreaReductionFactor[n] > 0.)) {
                            //G_locWetlAreaReductionFactor[n] = G_locWetlAreaReductionFactor[n] * G_fLocWet[n] / ((double)G_loc_wetland[n]/100.); // adapt areareductionfactor
                            G_fLocWet[n] = (G_locWetlAreaReductionFactor[n] * G_loc_wetland[n] /
                                            100.); // and finalize swb fraction
                        } else {
                            G_locWetlAreaReductionFactor[n] = 0.;
                            G_fLocWet[n] = 0.;
                        }

                        if ((G_fGloWet[n] > 0.) && (G_gloWetlAreaReductionFactor[n] > 0.)) {
                            //G_gloWetlAreaReductionFactor[n] = G_gloWetlAreaReductionFactor[n] * G_fGloWet[n] / ((double)G_glo_wetland[n]/100.); // adapt areareductionfactor
                            G_fGloWet[n] = (G_gloWetlAreaReductionFactor[n] * G_glo_wetland[n] /
                                            100.); // and finalize swb fraction
                        } else {
                            G_gloWetlAreaReductionFactor[n] = 0.;
                            G_fGloWet[n] = 0.;
                        }

                    }
                    // endif (G_fswbLandAreaFracNextTimestep[n] > (G_maxRiverAreaFrac[n] - G_riverAreaFracNextTimestep[n]))

                }
                    // endif //river fraction can reach max value
                else { //river fraction cannot reach max value, thus has to be limited as well as river evaporation
                    //  cout << "riverfrac cannot reach max at n: " << n << " and day: " << day << " glolak: " << G_glo_lake[n] << " riverfrac: " << G_riverAreaFrac[n] << endl;
                    riverAreaFracDeficit = G_riverAreaFracNextTimestep_Frac[n] -
                                           G_maxRiverAreaFrac[n]; //calculate missing fraction
                    G_riverAreaFrac_ChangeFrac[n] -= riverAreaFracDeficit; //reduce fraction change
                    G_riverAreaReductionFactor[n] *= G_maxRiverAreaFrac[n] /
                                                     G_riverAreaFracNextTimestep_Frac[n];//adapt (reduce) riverAreaReductionFactor
                    G_riverAreaFracNextTimestep_Frac[n] = G_maxRiverAreaFrac[n]; //set fraction to max fraction

                    G_fLocLake[n] = 0.;
                    G_fLocWet[n] = 0.;
                    G_fGloWet[n] = 0.;
                }

            }
            // endelse (no inland sink and no 100% glo lake)
            //G_dailyRiverEvapo[n] = (double)dailyWaterBalance.G_openWaterPET[n] * G_riverAreaFrac[n] / 100.; // Evaporation from rivers in mm of cell area
            //dailyRiverEvapo = G_dailyRiverEvapo[n] * cellArea *(((double)geo.G_landfreq[n] + (double)geo.G_fwaterfreq[n]) / 100. / 1000000.); //evaporation from rivers in km3

            G_fswbLandAreaFracNextTimestep[n] = G_fLocLake[n] + G_fLocWet[n] + G_fGloWet[n];
            G_fswbLandAreaFrac_ChangePct[n] = G_fswbLandAreaFracNextTimestep[n] * 100. - G_fswbLandAreaFrac[n] *
                                                                                         100.; // G_fwsb is in fraction, G_landAreaFrac (where G_fswbChange is calculated) is in %

        }
            // endif (options.riverEvapoOpt == 1)
        else {
            //no evaporation from rivers
            G_riverAreaFracNextTimestep_Frac[n] = 0.;
            G_riverAreaFrac_ChangeFrac[n] = 0.;
        }
        // end treatment of river area (for evaporation)


        //calculate final G_fswb (should be now max. continental frac)  // FP HERE possibly reservoir fraction is missing?? // HMS: I agree, implemented
        G_fswb[n] = G_fLocLake[n] + G_fLocWet[n] + G_fGloLake[n] + G_fGloWet[n] + G_glo_res[n] +
                    G_riverAreaFracNextTimestep_Frac[n]; // for output purposes only!


        //}
        //else {//no fswb change
        //	G_fswbLandAreaFracChange[n] = 0.;
        //}

        //
        //now, adapt land area fraction for changes in fswb and riverareafrac
        // First model day special condition: subtract changes from G_landAreaFrac[n]
        if (0 == statusStarted_landAreaFracNextTimestep[n]) {
            G_landAreaFracNextTimestep[n] = G_landAreaFrac[n] -
                                            (G_fswbLandAreaFrac_ChangePct[n] + (G_riverAreaFrac_ChangeFrac[n] * 100.));
            statusStarted_landAreaFracNextTimestep[n] = 1; // indicate for the next iteration that first day has been treated
            /*  if (n == 40812) {
                    cout << "first model day special condition! " << endl;
                    cout << "G_landAreaFracNextTimestep[n] " << G_landAreaFracNextTimestep[n] << endl;
                    cout << "G_landAreaFrac[n] " << G_landAreaFrac[n] << endl;
                    cout << "G_fswbLandAreaFrac_ChangePct[n] " << G_fswbLandAreaFrac_ChangePct[n] << endl;
                }*/
        }
            // Following model days: only apply changes
        else {
            //on exception: reservoir fraction is added in this year, consider G_glores_change once.
            if (G_glores_change[n] > 0.) {
                G_landAreaFracNextTimestep[n] -= G_glores_change[n];
                G_glores_change[n] = 0.;
            }
            G_landAreaFracNextTimestep[n] -= (G_fswbLandAreaFrac_ChangePct[n] + (G_riverAreaFrac_ChangeFrac[n] * 100.));
            /*  if (n == 40812) {
                    cout << "any other model day condition! " << endl;
                    cout << "G_landAreaFracNextTimestep[n] " << G_landAreaFracNextTimestep[n] << endl;
                    cout << "G_fswbLandAreaFrac_ChangePct[n] " << G_fswbLandAreaFrac_ChangePct[n] << endl;
                }*/
        }



        // Felix Portmann FP & Petra Döll 2015
        // TEMPORARILY REMOVE - POSSIBLY REINTRODUCE
        // Introduce land area only when surface water body area has decreased so much that a limit of cell area
        // (limit_zero_landAreaFrac_pct, e.g. 1% of cell area) is land area
//            if ( (0. == G_landAreaFrac[n])
//                 && ( ((geo.G_landfreq[n]+geo.G_fwaterfreq[n]) - (G_fswb[n]*100.)) >= limit_zero_landAreaFrac_pct )
//                    ) {
//                G_landAreaFracNextTimestep[n] = ((geo.G_landfreq[n]+geo.G_fwaterfreq[n]) - (G_fswb[n]*100.));
//            }

        // FP & CR 2015-06
        // Control not allowed range of values: Set to zero if calculated land fraction is less than a limit
        // limit currently set to zero
        // e.g. also for Caspian Sea cells
        if (G_landAreaFracNextTimestep[n] < 0.) {
            G_landAreaFracNextTimestep[n] = 0.;
        }

        //else // ((0 == options.fractionalRoutingOpt) && (options.aridareaOpt == 0))
        //G_landAreaFracNewDay[n] = -99.0;

        // HMS 2014-10-28 moved from lines around 920 to here so that area is consistent for daily.cpp calculations
        //if (G_landAreaFracNewDay[n] > -99.0) // to prevent the initialisation with at the first time step
        //G_landAreaFrac[n] = G_landAreaFracNewDay[n];
        // Sum values for monthly averaging // only with specific writing options // Felix Portmann 2015-06
        if ((2 == options.grid_store) || (4 == options.grid_store) || (6 == options.grid_store)) {
            G_monthlyLandAreaFrac[n][month] += G_landAreaFrac[n];
        }

    }
    // endfor loop over grid cells
    // end of computation of SWB and land area fractions for the next day

    //HMS 2015-03-02 now write out final G_fswb
    for (int n = 0; n < ng; n++) {
        if ((2 == options.grid_store) || (4 == options.grid_store) || (6 == options.grid_store)) {
            if (options.outFswb == 1)
                G_monthlyFswb[n][month] += G_fswb[n];
        }
        if ((3 == options.grid_store) || (4 == options.grid_store)) {
            if (options.outFswbDaily)
                G_daily31Fswb[n][day_in_month - 1] = G_fswb[n]; // daily fraction of surface water bodies [mm]
        }
        if ((5 == options.grid_store) || (6 == options.grid_store)) {
            if ((options.outFswbDaily) || ((options.scoutFswb) && (2 == options.day_store)))
                G_daily365Fswb[n][day - 1] = G_fswb[n];
        }
    }
    // endfor loop over grid cells

#ifdef _WATERGAP_CHECKS_GLOBAL_H
                                                                                                                            // Check Felix Portmann 2015 - for invalid data
            if ( (year >= chk_year) && (month >= chk_month) && (day_in_month >= chk_day_in_month) ) {


                // Checks for values less than zero LT0

                // grids 1
                check_array_double_LT0_YYMMDD(G_riverInflow, ng, year, month, day_in_month);
                check_array_double_LT0_YYMMDD(G_riverStorage, ng, year, month, day_in_month);
                check_array_double_LT0_YYMMDD(G_riverLength, ng, year, month, day_in_month);
                check_array_double_LT0_YYMMDD(G_riverAreaFrac, ng, year, month, day_in_month);
                check_array_double_LT0_YYMMDD(G_dailyRiverEvapo, ng, year, month, day_in_month);
                check_array_double_LT0_YYMMDD(G_dailyRiverPrecip, ng, year, month, day_in_month);

                // grids 2
                check_array_double_LT0_YYMMDD(G_localRunoff, ng, year, month, day_in_month);
                check_array_double_LT0_YYMMDD(G_localRunoffIntoRiver, ng, year, month, day_in_month);
                check_array_double_LT0_YYMMDD(G_localGWRunoff, ng, year, month, day_in_month);
                check_array_double_LT0_YYMMDD(G_localGWRunoffIntoRiver, ng, year, month, day_in_month);
                check_array_short_LT0_YYMMDD(G_aindex, ng, year, month, day_in_month);

                // fractions
                // land area
                check_array_double_LT0_YYMMDD(G_landAreaFrac, ng, year, month, day_in_month);
                check_array_double_LT0_YYMMDD(G_landWaterExclGloLakAreaFrac, ng, year, month, day_in_month);
                check_array_double_LT0_YYMMDD(G_landAreaFracNextTimestep, ng, year, month, day_in_month);
                check_array_double_LT0_YYMMDD(G_landAreaFracPrevTimestep, ng, year, month, day_in_month);

                // river evaporation
                check_array_double_LT0_YYMMDD(G_riverAreaFrac, ng, year, month, day_in_month);
                check_array_double_LT0_YYMMDD(G_maxRiverAreaFrac, ng, year, month, day_in_month);
                check_array_double_LT0_YYMMDD(G_riverAreaFracNextTimestep_Frac, ng, year, month, day_in_month);
        //        check_array_double_LT0_YYMMDD(G_riverAreaFracChange, ng, year, month, day_in_month);

                // Surface Water Bodies swb
                check_array_double_LT0_YYMMDD(G_fswbLandAreaFracNextTimestep, ng, year, month, day_in_month);
                check_array_double_LT0_YYMMDD(G_fswbLandAreaFracInit, ng, year, month, day_in_month);
                check_array_double_LT0_YYMMDD(G_fswbLandAreaFrac, ng, year, month, day_in_month);
        //        check_array_double_LT0_YYMMDD(G_fswbLandAreaFracChange, ng, year, month, day_in_month);
                check_array_double_LT0_YYMMDD(G_fLocLake, ng, year, month, day_in_month);
                check_array_double_LT0_YYMMDD(G_fLocWet, ng, year, month, day_in_month);
                check_array_double_LT0_YYMMDD(G_fGloWet, ng, year, month, day_in_month);
                check_array_double_LT0_YYMMDD(G_fGloLake, ng, year, month, day_in_month);
                check_array_double_LT0_YYMMDD(G_fswb, ng, year, month, day_in_month);
                check_array_double_LT0_YYMMDD(G_fswbInit, ng, year, month, day_in_month);
                check_array_short_LT0_YYMMDD(G_fractionInit, ng, year, month, day_in_month);
                check_array_double_LT0_YYMMDD(G_fswb_catchment, ng, year, month, day_in_month);

                // scalars
                check_value_double_LT0(cellArea);
                check_value_double_LT0(CellRunoff);
                check_value_double_LT0(transportedVolume);
                check_value_double_LT0(dailyRiverEvapo);
                check_value_double_LT0(inflow);
                check_value_double_LT0(inflowUpstream);
                check_value_double_LT0(totalInflow);
                check_value_double_LT0(outflow);
                check_value_double_LT0(maxStorage);
                check_value_double_LT0(K);
                // check_value_double_LT0(potCellRunoff); // Update 2015-06-08
                check_value_double_LT0(riverAvail);
                check_value_double_LT0(Kswbgw);
                check_value_double_LT0(fswbFracCorr);
                check_value_double_LT0(riverAreaFracDeficit);
                check_value_double_LT0(fLocLak);
                check_value_double_LT0(fLocWet);
                check_value_double_LT0(fGloWet);
                check_value_double_LT0(locLakeAreaReductionFactor);
                check_value_double_LT0(gloLakeEvapoReductionFactor);
                check_value_double_LT0(locWetlAreaReductionFactor);
                check_value_double_LT0(gloWetlAreaReductionFactor);
                check_value_double_LT0(gloResEvapoReductionFactor);
                check_value_double_LT0(riverAreaReductionFactor);
                check_value_double_LT0(PETgwr);
                check_value_double_LT0(PETgwrMax);
                check_value_double_LT0(PETgwrRedFactor);
                check_value_double_LT0(PETgwrRemUseRedFactor);
                check_value_double_LT0(gloResEvapoReductionFactor);
                check_value_double_LT0(riverAreaReductionFactor);
                check_value_double_LT0(groundwater_runoff_km3);
                // check_value_double_LT0(G_groundwaterStoragePrevStep); // G_groundwaterStoragePrevStep may get negative because of netGWin // Felix Portmann 2015-06
                // additionally:
                check_value_double_LT0(b); // counter variable
                check_value_double_LT0(Kswbgw); // fixed value


                // Check for NaN

                // grids 1
                check_array_double_NaN_YYMMDD(G_riverInflow, ng, year, month, day_in_month);
                check_array_double_NaN_YYMMDD(G_riverStorage, ng, year, month, day_in_month);
                check_array_double_NaN_YYMMDD(G_riverLength, ng, year, month, day_in_month);
                check_array_double_NaN_YYMMDD(G_riverAreaFrac, ng, year, month, day_in_month);
                check_array_double_NaN_YYMMDD(G_dailyRiverEvapo, ng, year, month, day_in_month);
                check_array_double_NaN_YYMMDD(G_dailyRiverPrecip, ng, year, month, day_in_month);

                // grids 2
                check_array_double_NaN_YYMMDD(G_localRunoff, ng, year, month, day_in_month);
                check_array_double_NaN_YYMMDD(G_localRunoffIntoRiver, ng, year, month, day_in_month);
                check_array_double_NaN_YYMMDD(G_localGWRunoff, ng, year, month, day_in_month);
                check_array_double_NaN_YYMMDD(G_localGWRunoffIntoRiver, ng, year, month, day_in_month);
                check_array_short_NaN_YYMMDD(G_aindex, ng, year, month, day_in_month);

                // fractions
                // land area
                check_array_double_NaN_YYMMDD(G_landAreaFrac, ng, year, month, day_in_month);
                check_array_double_NaN_YYMMDD(G_landWaterExclGloLakAreaFrac, ng, year, month, day_in_month);
                check_array_double_NaN_YYMMDD(G_landAreaFracNextTimestep, ng, year, month, day_in_month);
                check_array_double_NaN_YYMMDD(G_landAreaFracPrevTimestep, ng, year, month, day_in_month);

                // river evaporation
                check_array_double_NaN_YYMMDD(G_riverAreaFrac, ng, year, month, day_in_month);
                check_array_double_NaN_YYMMDD(G_maxRiverAreaFrac, ng, year, month, day_in_month);
                check_array_double_NaN_YYMMDD(G_riverAreaFracNextTimestep_Frac, ng, year, month, day_in_month);
                check_array_double_NaN_YYMMDD(G_riverAreaFrac_ChangeFrac, ng, year, month, day_in_month);

                // Surface Water Bodies swb
                check_array_double_NaN_YYMMDD(G_fswbLandAreaFracNextTimestep, ng, year, month, day_in_month);
                check_array_double_NaN_YYMMDD(G_fswbLandAreaFracInit, ng, year, month, day_in_month);
                check_array_double_NaN_YYMMDD(G_fswbLandAreaFrac, ng, year, month, day_in_month);
                check_array_double_NaN_YYMMDD(G_fswbLandAreaFrac_ChangePct, ng, year, month, day_in_month);
                check_array_double_NaN_YYMMDD(G_fLocLake, ng, year, month, day_in_month);
                check_array_double_NaN_YYMMDD(G_fLocWet, ng, year, month, day_in_month);
                check_array_double_NaN_YYMMDD(G_fGloWet, ng, year, month, day_in_month);
                check_array_double_NaN_YYMMDD(G_fGloLake, ng, year, month, day_in_month);
                check_array_double_NaN_YYMMDD(G_fswb, ng, year, month, day_in_month);
                check_array_double_NaN_YYMMDD(G_fswbInit, ng, year, month, day_in_month);
                check_array_short_NaN_YYMMDD(G_fractionInit, ng, year, month, day_in_month);
                check_array_double_NaN_YYMMDD(G_fswb_catchment, ng, year, month, day_in_month);


                // Scalars
                check_value_double_NaN(cellArea);
                check_value_double_NaN(CellRunoff);
                check_value_double_NaN(transportedVolume);
                check_value_double_NaN(dailyRiverEvapo);
                check_value_double_NaN(inflow);
                check_value_double_NaN(inflowUpstream);
                check_value_double_NaN(totalInflow);
                check_value_double_NaN(outflow);
                check_value_double_NaN(maxStorage);
                check_value_double_NaN(K);
                check_value_double_NaN(potCellRunoff);
                check_value_double_NaN(riverAvail);
                check_value_double_NaN(Kswbgw);
                check_value_double_NaN(fswbFracCorr);
                check_value_double_NaN(riverAreaFracDeficit);
                check_value_double_NaN(fLocLak);
                check_value_double_NaN(fLocWet);
                check_value_double_NaN(fGloWet);
                check_value_double_NaN(locLakeAreaReductionFactor);
                check_value_double_NaN(gloLakeEvapoReductionFactor);
                check_value_double_NaN(locWetlAreaReductionFactor);
                check_value_double_NaN(gloWetlAreaReductionFactor);
                check_value_double_NaN(gloResEvapoReductionFactor);
                check_value_double_NaN(riverAreaReductionFactor);
                check_value_double_NaN(PETgwr);
                check_value_double_NaN(PETgwrMax);
                check_value_double_NaN(PETgwrRedFactor);
                check_value_double_NaN(PETgwrRemUseRedFactor);
                check_value_double_NaN(gloResEvapoReductionFactor);
                check_value_double_NaN(riverAreaReductionFactor);
                check_value_double_NaN(groundwater_runoff_km3);
                check_value_double_NaN(G_groundwaterStoragePrevStep);

                // additionally:
                check_value_double_NaN(b); // counter variable
                check_value_double_NaN(Kswbgw); // fixed value
                check_value_double_NaN(netGWin); // possibly negative value (ONLY NaN testing)

            }
            // endif check date reached // Check Felix Portmann 2015
#endif

}
// end routing
// end of routing


void routingClass::annualWaterUsePostProcessing(short year) {
    if (options.subtract_use >= 1) {
        char filename[250];
        int n;

        // Has there been a satisfaction of the water use of the previous year?
        // We want to store a grid with the amount of use which is still
        // unsatisfied from the previous year for control purposes.
        // After this year the unsatisfied use of the previous year will
        // not be considered any more. Therefore it is removed from
        // the grid 'G_totalUnsatisfiedUse'.
        for (n = 0; n <= ng - 1; n++) {
            if (G_totalUnsatisfiedUse[n] < G_UnsatisfiedUsePrevYear[n]) {
                // This means: all the desired uses of the current year have been
                // satisfied and a fraction of the previous year also.
                G_UnsatisfiedUsePrevYear[n] = G_totalUnsatisfiedUse[n];
                G_totalUnsatisfiedUse[n] = 0.;
            } else {
                // This means: we have even more unsatisfied use than at the end
                // of the previous year.
                // The grid of the year before has not to be changed,
                // but the grid of the current year should only contain the values
                // of the current year.
                G_totalUnsatisfiedUse[n] -= G_UnsatisfiedUsePrevYear[n];
            }
        }

        if ((options.grid_store > 0) && (year >= options.evalStartYear)) {
            double G_array[ng];

            if (options.outUnsatUseSWPrev) { // new output options 2.1f
                for (n = 0; n < ng; n++)
                    G_array[n] = G_UnsatisfiedUsePrevYear[n]; // conversion double -> float
                sprintf(filename, "%s/G_UNSAT_USE_SW_PREV_%d.UNF0", options.output_dir, year);
                gridIO.writeUnfFile(filename, ng, &G_array[0]);
            }
            if (options.outUnsatUseSW) { // new output options 2.1f
                for (n = 0; n < ng; n++)
                    G_array[n] = (float) G_totalUnsatisfiedUse[n]; // conversion double -> float
                sprintf(filename, "%s/G_UNSAT_USE_SW_%d.UNF0", options.output_dir, year);
                gridIO.writeUnfFile(filename, ng, &G_array[0]);
            }
            if (options.outActualNAS) { // new output options 2.1f
                for (n = 0; n < ng; n++) {
                    G_actualUse[n] = G_actualUse[n] *
                                     1000000000.; // conversion km3/month --> m3/month to avoid loss of small numbers during casting operation.
                    G_array[n] = G_actualUse[n]; // conversion double -> float
                }
                sprintf(filename, "%s/G_ACTUAL_NAS_%d.UNF0", options.output_dir, year);
                gridIO.writeUnfFile(filename, ng, &G_array[0]);
            }
            // detect use satisfaction from neighbour cells
            /*
			 for (n = 0; n < ng; n++)
			 G_array[n] = (float) G_extractedFromNeighb[n]; // conversion double -> float
			 sprintf(filename, "%s/G_EXTRACTED_FROM_NEIGHB_%d.UNF0", options.output_dir, year);
			 gridIO.writeUnfFile(filename, ng, &G_array[0]);
			 for (n = 0; n < ng; n++)
			 G_array[n] = (float) G_deliveredToNeighb[n]; // conversion double -> float
			 sprintf(filename, "%s/G_DELIVERED_TO_NEIGHB_%d.UNF0", options.output_dir, year);
			 gridIO.writeUnfFile(filename, ng, &G_array[0]);
			 */

        }
        // Note: ACTUAL_USE_[YEAR] + UNSAT_USE_[YEAR] does not equal CONS_WATER_USE_[YEAR]
        // as unsatisfied use from the previous year (YEAR-1) might have partially
        // been satisfied in YEAR
        // ACTUAL_USE_[YEAR] + UNSAT_USE_[YEAR] = CONS_WATER_USE[YEAR] + (UNSAT_USE_[YEAR-1] - UNSAT_USE_PREV_[YEAR])

    }

}
// end annualWaterUsePostProcessing

void routingClass::checkTimeStep() {
    // this is the public function
    // should be called after G_toBeCaluclated has been generated
    // and after initialization of this class
    defaultTimeStepsPerDay = checkTimeStep(defaultRiverVelocity, defaultTimeStepsPerDay);
    // velocity
    riverVelocity = defaultRiverVelocity; // [km/day]
}
// end checkTimeStep

void routingClass::updateLandAreaFrac() {
    //LandAreaFracNextTimestep is already calculated, now G_LandAreaFrac has to be updated to be used consistently in daily and routing for the next timestep
    for (int n = 0; n < ng; n++) {
        G_landAreaFracPrevTimestep[n] = G_landAreaFrac[n];
        G_landAreaFrac[n] = G_landAreaFracNextTimestep[n];
    }
}
// end updateLandAreaFrac

// FP20161018N002 Reservoir operation start years
// Transfer land area fraction of reservoirs (in percent) of current year in array of previous year for use in the next year
void routingClass::updateGloResPrevYear_pct() {
    for (int n = 0; n < ng; n++) {
        G_glores_prevyear[n] = G_glo_res[n];
    }
}
// end updateGloResPrevYear
/**
 *
 * @brief This class updates the groundwater net abstraction due to not statisfied SW usage
 *
 * This method corrects the groundwater net abstraction in case of not fullfilled net abstraction from surface water
 * bodies(swb). GWSWUSE is assuming that all water abstraction from swb is fullfilled and calculates gw recharge
 * based on this assumption. Thus groundwater recharge can be overestimated. Hence the unstatisfied surface water
 * use from previous time step is compared to irrigation usage from previous time step. The GW recharge of current
 * time step is then adapted, assuming the remaining use is solely irrigation use.
 * @param n
 * @param month
 */
double routingClass::updateNetAbstractionGW(int n, int month) {

    // Only update NAg if WUsi is existent
    if (G_monthlyWUIrrig[n][month] > 0) {
        double daily_WUIrrig = G_monthlyWUIrrig[n][month] / (1000000000. * (double) numberOfDaysInMonth[month]);
        double daily_CUIrrig = G_monthlyCUIrrig[n][month] / (1000000000. * (double) numberOfDaysInMonth[month]);
        double eff = daily_CUIrrig / daily_WUIrrig;
        double frgi = G_fractreturngw_irrig[n];
        double factor = 1 - (1 - frgi) * (1 - eff);
        double wusi_new = 1 / factor * (daily_WUIrrig * factor - G_dailyRemainingUse[n]);
        double NAgnew = G_dailydailyNUg[n] - (frgi * (1 - eff) * (wusi_new - daily_WUIrrig));
        return NAgnew;
    }
    else{
        return G_dailydailyNUg[n];
    }
}

short routingClass::checkTimeStep(double riverVelocity, short numberOfTimeSteps) {
    extern short G_toBeCalculated[ng];

    // find most nothern and southern row
    short minRow = 1000;
    short maxRow = 0;

    for (int n = 0; n < ng; n++) {
        if (G_toBeCalculated[n] != 0) {
            if (geo.G_row[n] < minRow)
                minRow = geo.G_row[n];
            if (geo.G_row[n] > maxRow)
                maxRow = geo.G_row[n];
        }
    }

    // transfer to northern hemisphere
    if (minRow > 180)
        minRow = 360 - minRow + 1;
    if (maxRow > 180)
        maxRow = 360 - maxRow + 1;

    // find most northern row
    if (maxRow < minRow)
        minRow = maxRow;

    // calculate horizontal extend of cell
    const double pi = 3.141592653589793;
    const double earthRadius = 6371.211; // [km]
    double lat, horizDist;

    lat = 0.25 + minRow * 0.5;
    horizDist = 2 * pi * earthRadius * sin(lat * pi / 180.0) * 0.5 / 360.0;

    //cout << "Horizontal extend of smallest cell: " << horizDist << endl;

    short newNumberOfTimeSteps;

    newNumberOfTimeSteps = numberOfTimeSteps;
    if (riverVelocity / (horizDist * numberOfTimeSteps) > 0.75) {
        cout << "Time step to great: " << numberOfTimeSteps << " steps per day\n";
        newNumberOfTimeSteps = (short) ceil(riverVelocity / (0.6 * horizDist));
        cout << "Adjusted to: " << newNumberOfTimeSteps << endl;
    }
    return newNumberOfTimeSteps;
}
// end checkTimeStep

// Set (initialize) cell-specific active lake depth (in km) as array to enable individual time-series  // FP
// Start with calibration parameter
void routingClass::initLakeDepthActive() {

    // Initialize array with calibration parameter values
    // Data unit of depth has to be converted (m to km)
    for (int n = 0; n < ng; n++) {
        G_lakeDepthActive[n] = (calibParam.getValue(P_LAK_D, n) *
                                CONV_FACTOR__LENGTH__M_TO_KM); // [km]  // OBSOLETE now read from calibration parameter file // FP
    }

}
// end setLakeDepthActive

// Set (initialize) cell-specific active wetland depth (in km) as array to enable individual time-series  // FP
// Start with calibration parameter
void routingClass::initWetlDepthActive() {

    // Initialize array with calibration parameter values
    // Data unit of depth has to be converted (m to km)
    for (int n = 0; n < ng; n++) {
        G_wetlDepthActive[n] = (calibParam.getValue(P_WET_D, n) *
                                CONV_FACTOR__LENGTH__M_TO_KM); // [km]  // OBSOLETE now read from calibration parameter file // FP
    }

}
// end setWetlDepthActive


void routingClass::setLakeWetlToMaximum(const short start_year) {
    // execute at start of each calibration year
    // or at first model year // FP20161018N002: any possible specified year
    for (int n = 0; n <= ng - 1; n++) {

        // Cell-specific active lake depth - Use array (initialized from cell-specific parameter)  // FP
        // Cell-specific active wetland depth - Use array (initialized from cell-specific parameter)  // FP
        G_locLakeStorage[n] = (G_loc_lake[n] / 100.)
                              * geo.areaOfCellByArrayPos(n) * G_lakeDepthActive[n];
        G_locWetlStorage[n] = (G_loc_wetland[n] / 100.)
                              * geo.areaOfCellByArrayPos(n) * G_wetlDepthActive[n];

        G_gloLakeStorage[n] = G_lake_area[n] *
                              G_lakeDepthActive[n]; // added for new reservoir algorithm  (remark: int G_lake_area, in full km2)

        G_gloWetlStorage[n] = (G_glo_wetland[n] / 100.)
                              * geo.areaOfCellByArrayPos(n) * G_wetlDepthActive[n];

        // Use reservoir algorithm with reservoirs and regulated lakes
        if ((options.antNatOpt == 0) && (options.resOpt == 1)) {
//			G_gloResStorage[n] = G_stor_cap[n] * 0.85; //for new reservoir algorithm (vol. from published data) // Obsolete

            // Under certain conditions set storage capacity to 85% of full capacity

            //set to storage capacity if no operational year and start year not later than reference year (e.g. 2000)
            // for both, reg. lakes and reservoirs
            if ((options.resYearOpt == 0) && (options.resYearReference >= G_res_start_year[n]) &&
                (G_stor_cap_full[n] > -99)) {
                G_gloResStorage[n] = G_stor_cap_full[n] * 0.85; //for new reservoir algorithm (vol. from published data)
            }
            //set reservoir storage to capacity for reservoirs and regulated lakes which are already operating at model start year
            if ((options.resYearOpt == 1) && (start_year >= G_res_start_year[n]) && (G_stor_cap_full[n] > -99)) {
                G_gloResStorage[n] = G_stor_cap_full[n] * 0.85; //for new reservoir algorithm (vol. from published data)
            }

// HERE CHECK IF CORRECT / NO DOUBLE COUNTING
            // In case of regulated lake: storage = product of area * lake_depth
            //increase reservoir storage if regulated lake is not yet operating but handle it like if it would be a global lake
            if (((options.resYearOpt == 1) && (G_reg_lake_status[n] == 1) && (start_year < G_res_start_year[n]))
                || ((options.resYearOpt == 0) && (G_reg_lake_status[n] == 1) &&
                    (options.resYearReference < G_res_start_year[n]))) {
                G_gloResStorage[n] += G_reservoir_area_full[n] * G_lakeDepthActive[n];
            }
            //increase global lake storage if regulated lake and resOpt == 0
            //if ((options.resOpt == 0) && (G_reg_lake_status[n] == 1)) {
            //    G_gloLakeStorage[n] += G_reservoir_area_full[n] * G_lakeDepthActive[n];
            //}
        }
        // end if (options.resOpt == 1)

        //initialize area reduction factors
        // initial value at start of first model year or at each calibration year)
        G_locLakeAreaReductionFactor[n] = 1.;
        G_locWetlAreaReductionFactor[n] = 1.;
        G_gloLakeEvapoReductionFactor[n] = 1.;
        G_gloWetlAreaReductionFactor[n] = 1.;
        G_gloResEvapoReductionFactor[n] = 1.;
        if (options.riverEvapoOpt == 1) {
            G_riverAreaReductionFactor[n] = 0.5; //corresponds to around 20% of bankfull storage for the first model day
            // river area fraction (as fraction of cell area)
            // with river width of bankfull flow [unit: m]
            // ATTENTION:
            // Area adds to water area!
            // Make sure that total area does not exceed cell area!
            // E.g. subtract from surface water bodies (e.g. lakes, wetlands)
            G_riverAreaFracNextTimestep_Frac[n] =
                           G_riverAreaReductionFactor[n] * G_riverLength[n] * G_RiverWidth_bf[n] / 1000. /
                           geo.areaOfCellByArrayPos(n);
        } else {
            //no evaporation from rivers
            G_riverAreaReductionFactor[n] = 0.;
            G_riverAreaFracNextTimestep_Frac[n] = 0.;
            G_riverAreaFrac_ChangeFrac[n] = 0.;
        }
        // end treatment of river area (for evaporation)

    }
    // end loop over grid cells

}
// end setLakeWetlToMaximum

double routingClass::getActualStorage(int cellNumber, short selection) {
    switch (selection) {
        case '1':
            return G_locLakeStorage[cellNumber - 1];
        case '2':
            return G_locWetlStorage[cellNumber - 1];
        case '3':
            return G_gloLakeStorage[cellNumber - 1];
        case '4':
            return G_gloWetlStorage[cellNumber - 1];
            // ET: for new reservoir algorithm
        case '5':
            return G_gloResStorage[cellNumber - 1];
        default:
            cerr << "Your selection is not available\n";
            return -9999;
    }
}
// end getActualStorage

double routingClass::getActualStorageRatio(int cellNumber, short selection) {
    double maxStorage;

    switch (selection) {
        case 1: {
            if (G_loc_lake[cellNumber - 1] <= 0)
                return -9999;
            else {
                // Cell-specific active lake depth - Use array (initialized from cell-specific parameter)  // FP
                maxStorage = (G_loc_lake[cellNumber - 1] / 100.0)
                             * geo.areaOfCell(cellNumber) * G_lakeDepthActive[cellNumber - 1];
                return G_locLakeStorage[cellNumber - 1] / maxStorage;
            }
        }
        case 2: {
            if (G_loc_wetland[cellNumber - 1] <= 0)
                return -9999;
            else {
                // Cell-specific active wetland depth - Use array (initialized from cell-specific parameter)  // FP
                maxStorage = (G_loc_wetland[cellNumber - 1] / 100.0)
                             * geo.areaOfCell(cellNumber) * G_wetlDepthActive[cellNumber - 1];
                return G_locWetlStorage[cellNumber - 1] / maxStorage;
            }
        }
        case 3: {
///			if (G_glo_lake[cellNumber - 1] <= 0)
            if (G_lake_area[cellNumber - 1] <= 0) // ET: for new reservoir algorithm
                return -9999;
            else {

///				maxStorage = (G_glo_lake[cellNumber - 1] / 100.0)
///					* geo.areaOfCell(cellNumber) * lakeDepth;
                // Cell-specific active lake depth - Use array (initialized from cell-specific parameter)  // FP
                maxStorage = G_lake_area[cellNumber - 1] * G_lakeDepthActive[cellNumber - 1];

                return G_gloLakeStorage[cellNumber - 1] / maxStorage;
            }
        }
        case 4: {
            if (G_glo_wetland[cellNumber - 1] <= 0)
                return -9999;
            else {
                // Cell-specific active wetland depth - Use array (initialized from cell-specific parameter)  // FP
                maxStorage = (G_glo_wetland[cellNumber - 1] / 100.0)
                             * geo.areaOfCell(cellNumber) * G_wetlDepthActive[cellNumber - 1];
                return G_gloWetlStorage[cellNumber - 1] / maxStorage;
            }
        }
            // for new reservoir algorithm
        case 5: {
            if (G_reservoir_area[cellNumber - 1] <= 0)
                return -9999;
            else {

                maxStorage = G_stor_cap[cellNumber - 1] * 0.85; //85 % of storage capacity (from published volume data)

                return G_gloResStorage[cellNumber - 1] / maxStorage;
            }
        }

        default:
            cerr << "Selection not available\n";
            return -9999;
    }
}
// end getActualStorageRatio


// FP20161018N002 Reservoir operation start years
double routingClass::writeAnnualReservoirCapacityandLandStorageChange(char *outputDir, const short year) {

    extern optionClass options;

    char filename[250];

    // Reservoir capacity of current year
    if (options.outResCapacity) {
        int n;
        float G_array[ng];
        for (n = 0; n < ng; n++)
            G_array[n] = (float) G_actualStorageCapacity[n];
        sprintf(filename, "%s/G_RESERVOIR_STORAGE_CAPACITY_km3_%d.UNF0", outputDir, year);
        gridIO.writeUnfFile(filename, ng, &G_array[0]);
    }
    // end if option

    // Land storage change
    if (options.outLandStorageChange) {

        int n;
        float G_array[ng];

        for (n = 0; n < ng; n++)
            G_array[n] = (float) G_landStorageChange[n];
        sprintf(filename, "%s/G_LAND_STORAGE_CHANGE_km3_%d.UNF0", outputDir, year);
        gridIO.writeUnfFile(filename, ng, &G_array[0]);

        for (n = 0; n < ng; n++) { // for output in mm
            if (geo.G_contfreq[n] == 0) {
                G_array[n] = 0.;
            } else {
                G_array[n] = (float) G_landStorageChange[n]// conversion double -> float
                             / (geo.areaOfCellByArrayPos(n) * (geo.G_contfreq[n] / 100.)) * 1000000;
            }
        }
        // endfor loop over grid cells
        sprintf(filename, "%s/G_LAND_STORAGE_CHANGE_mm_%d.UNF0", outputDir, year);
        gridIO.writeUnfFile(filename, ng, &G_array[0]);
        //Now, set back G_landStorageChange to 0 (to avoid that this value occurs now in each year)
        for (n = 0; n < ng; n++) {
            if (G_landStorageChange[n] != 0.)
                G_landStorageChange[n] = 0.;
        }

    }
    // end if option

}
// end writeAnnualReservoirCapacityandLandStorageChange


void
routingClass::storeSingleCellDailyValuesToFile(const short year)   // HMS 2013-11-21 Storage of single cell values begin
{
    //extern cbasinClass cbasin;

    // store daily values of single cells defined in SINGLECELLS.DAT
    FILE *file_ptr;
    char filename[250];
    extern optionClass options;
    extern dailyWaterBalanceClass dailyWaterBalance;
    int actualSingleCellNumber;
    short calibcounter = calibGamma.getCallCounter();


    for (short i = 0; i < options.numberOfSingleCells; ++i) {
        actualSingleCellNumber = options.cellnr[i] - 1;
        sprintf(filename, "%s/SINGLE_CELL_DAILY_VALUES_%d_%d.%s", options.output_dir, year, calibcounter,
                options.cellname[i]);
        file_ptr = fopen(filename, "w");
        fprintf(file_ptr, "# Daily values of single cell '%s':\n", options.cellname[i]);
        fprintf(file_ptr, "Day\t");
        if (options.scoutTemp) fprintf(file_ptr, "Temperature\t");
        if (options.scoutExtRad) fprintf(file_ptr, "ExtRad\t");
        if (options.scoutShortDownRad) fprintf(file_ptr, "ShortDownRad\t");
        if (options.scoutAlbedo) fprintf(file_ptr, "Albedo\t");
        if (options.scoutShortUpRad) fprintf(file_ptr, "ShortUpRad\t");
        if (options.scoutNetShortRad) fprintf(file_ptr, "NetShortRad\t");
        if (options.scoutLongDownRad) fprintf(file_ptr, "LongDownRad\t");
        if (options.scoutLongUpRad) fprintf(file_ptr, "LongUpRad\t");
        if (options.scoutNetLongRad) fprintf(file_ptr, "NetLongRad\t");
        if (options.scoutNetRad) fprintf(file_ptr, "NetRad\t");
        if (options.scoutLAI) fprintf(file_ptr, "LAI\t");
        if (options.scoutKc) fprintf(file_ptr, "Kc\t");
        if (options.scoutLandPET) fprintf(file_ptr, "LandPET\t");
        if (options.scoutCellPET) fprintf(file_ptr, "CellPET\t");
        if (options.scoutPrecip) fprintf(file_ptr, "Precip\t");
        if (options.scoutcPrecip) fprintf(file_ptr, "cPrecip\t");
        if (options.scoutCanopyWater) fprintf(file_ptr, "CanopyWater\t");
        if (options.scoutmaxCanopyWater) fprintf(file_ptr, "maxCanopyWater\t");
        if (options.scoutInterception) fprintf(file_ptr, "Interception\t");
        if (options.scoutSnowfall) fprintf(file_ptr, "Snowfall\t");
        if (options.scoutSnowCov) fprintf(file_ptr, "SnowCov\t");
        if (options.scoutSnowmelt) fprintf(file_ptr, "Snowmelt\t");
        if (options.scoutSnowWater) fprintf(file_ptr, "SnowWater\t");
        if (options.scoutSoilWater) fprintf(file_ptr, "SoilWater\t");
        if (options.scoutSurfaceRunoff) fprintf(file_ptr, "SurfaceRunoff\t");
        if (options.scoutGwRunoff) fprintf(file_ptr, "GwRunoff\t");
        if (options.scoutGwRecharge) fprintf(file_ptr, "GwRecharge\t");
        if (options.scoutCellAET) fprintf(file_ptr, "CellAET\t");
        if (options.scoutCellRunoffkm3) fprintf(file_ptr, "CellRunoffkm3\t");
        if (options.scoutCellRunoffmm) fprintf(file_ptr, "CellRunoffmm\t");
        if (options.scoutCellSRunoff) fprintf(file_ptr, "CellSRunoff\t");
        if (options.scoutQ) fprintf(file_ptr, "Q\t");
        if (options.scoutFlowVelo) fprintf(file_ptr, "FlowVelo\t");
        if (options.scoutLocLake) fprintf(file_ptr, "LocLake\t");
        if (options.scoutLocWet) fprintf(file_ptr, "LocWet\t");
        if (options.scoutGloLake) fprintf(file_ptr, "GloLake\t");
        if (options.scoutReservoir) fprintf(file_ptr, "Reservoir\t");
        if (options.scoutGloWet) fprintf(file_ptr, "GloWet\t");
        if (options.scoutRiver) fprintf(file_ptr, "River\t");
        if (options.scoutSurfaceStor) fprintf(file_ptr, "SurfaceStor\t");
        if (options.scoutTWSkm3) fprintf(file_ptr, "TWSkm3\t");
        if (options.scoutTWSmm) fprintf(file_ptr, "TWSmm\t");
        if (options.scoutGwStor) fprintf(file_ptr, "GwStor\t");
        if (options.scoutCanopyStor) fprintf(file_ptr, "CanopyStor\t");
        if (options.scoutSnowStor) fprintf(file_ptr, "SnowStor\t");
        if (options.scoutSoilStor) fprintf(file_ptr, "SoilStor\t");
        if (options.scoutGwrSwb) fprintf(file_ptr, "GwrSwb\t");
        if (options.scoutFswb) fprintf(file_ptr, "Fswb\t");
        if (options.scoutLandAreaFrac) fprintf(file_ptr, "LandAreaFrac\t");
        fprintf(file_ptr, "Day\n");

        for (short day = 0; day <= 364; day++) {
            fprintf(file_ptr, "%3d\t", day + 1);

//CR 2015-09 TEST Optionen 19, 20, 28
// CR ###################################################################################################
//			fprintf(file_ptr, "%10.6f\t", G_daily365ActualUse[actualSingleCellNumber][day]);
//			fprintf(file_ptr, "%10.6f\t", G_daily365AllocUse[actualSingleCellNumber][day]);
//			fprintf(file_ptr, "%10.6f\t", G_daily365AllocUseNextDay[actualSingleCellNumber][day]);
//			fprintf(file_ptr, "%10.6f\t", G_daily365TotalUnsatUse[actualSingleCellNumber][day]);
//			fprintf(file_ptr, "%10.6f\t", G_daily365UnsatAllocUseNextDay[actualSingleCellNumber][day]);
// CR ###################################################################################################

            if (options.scoutTemp)
                fprintf(file_ptr, "%10.6f\t", dailyWaterBalance.G_daily365TempC[actualSingleCellNumber][day]);
            if (options.scoutExtRad)
                fprintf(file_ptr, "%10.6f\t", dailyWaterBalance.G_daily365ExtRad[actualSingleCellNumber][day]);
            if (options.scoutShortDownRad)
                fprintf(file_ptr, "%10.6f\t", dailyWaterBalance.G_daily365ShortDownRad[actualSingleCellNumber][day]);
            if (options.scoutAlbedo)
                fprintf(file_ptr, "%10.6f\t", dailyWaterBalance.G_daily365Albedo[actualSingleCellNumber][day]);
            if (options.scoutShortUpRad)
                fprintf(file_ptr, "%10.6f\t", dailyWaterBalance.G_daily365ShortUpRad[actualSingleCellNumber][day]);
            if (options.scoutNetShortRad)
                fprintf(file_ptr, "%10.6f\t", dailyWaterBalance.G_daily365NetShortRad[actualSingleCellNumber][day]);
            if (options.scoutLongDownRad)
                fprintf(file_ptr, "%10.6f\t", dailyWaterBalance.G_daily365LongDownRad[actualSingleCellNumber][day]);
            if (options.scoutLongUpRad)
                fprintf(file_ptr, "%10.6f\t", dailyWaterBalance.G_daily365LongUpRad[actualSingleCellNumber][day]);
            if (options.scoutNetLongRad)
                fprintf(file_ptr, "%10.6f\t", dailyWaterBalance.G_daily365NetLongRad[actualSingleCellNumber][day]);
            if (options.scoutNetRad)
                fprintf(file_ptr, "%10.6f\t", dailyWaterBalance.G_daily365NetRadiation[actualSingleCellNumber][day]);
            if (options.scoutLAI)
                fprintf(file_ptr, "%10.6f\t", dailyWaterBalance.G_daily365LAI[actualSingleCellNumber][day]);
            if (options.scoutKc)
                fprintf(file_ptr, "%10.6f\t", dailyWaterBalance.G_daily365Kc[actualSingleCellNumber][day]);
            if (options.scoutLandPET)
                fprintf(file_ptr, "%10.6f\t", dailyWaterBalance.G_daily365PET[actualSingleCellNumber][day]);
            if (options.scoutCellPET)
                fprintf(file_ptr, "%10.6f\t", dailyWaterBalance.G_daily365TotalPET[actualSingleCellNumber][day]);
            if (options.scoutPrecip)
                fprintf(file_ptr, "%10.6f\t", dailyWaterBalance.G_daily365Precip[actualSingleCellNumber][day]);
            if (options.scoutcPrecip)
                fprintf(file_ptr, "%10.6f\t", G_daily365ConsistentPrecip[actualSingleCellNumber][day]);
            if (options.scoutCanopyWater)
                fprintf(file_ptr, "%10.6f\t",
                        dailyWaterBalance.G_daily365canopyWaterContent[actualSingleCellNumber][day]);
            if (options.scoutmaxCanopyWater)
                fprintf(file_ptr, "%10.6f\t",
                        dailyWaterBalance.G_daily365maxcanopyWaterContent[actualSingleCellNumber][day]);
            if (options.scoutInterception)
                fprintf(file_ptr, "%10.6f\t", dailyWaterBalance.G_daily365Interception[actualSingleCellNumber][day]);
            if (options.scoutSnowfall)
                fprintf(file_ptr, "%10.6f\t", dailyWaterBalance.G_daily365SnowFall[actualSingleCellNumber][day]);
            if (options.scoutSnowCov)
                fprintf(file_ptr, "%10.6f\t", dailyWaterBalance.G_daily365SnowCoverAvg[actualSingleCellNumber][day]);
            if (options.scoutSnowmelt)
                fprintf(file_ptr, "%10.6f\t", dailyWaterBalance.G_daily365SnowMelt[actualSingleCellNumber][day]);
            if (options.scoutSnowWater)
                fprintf(file_ptr, "%10.6f\t", dailyWaterBalance.G_daily365SWE[actualSingleCellNumber][day]);
            if (options.scoutSoilWater)
                fprintf(file_ptr, "%10.6f\t", dailyWaterBalance.G_daily365SoilWaterAvg[actualSingleCellNumber][day]);
            if (options.scoutSurfaceRunoff)
                fprintf(file_ptr, "%10.6f\t", dailyWaterBalance.G_daily365SurfaceRunoff[actualSingleCellNumber][day]);
            if (options.scoutGwRunoff) fprintf(file_ptr, "%10.6f\t", G_daily365GwRunoff[actualSingleCellNumber][day]);
            if (options.scoutGwRecharge)
                fprintf(file_ptr, "%10.6f\t", dailyWaterBalance.G_daily365GwRecharge[actualSingleCellNumber][day]);
            if (options.scoutCellAET) fprintf(file_ptr, "%10.6f\t", G_daily365CellAET[actualSingleCellNumber][day]);
            if (options.scoutCellRunoffkm3)
                fprintf(file_ptr, "%10.6f\t", G_daily365CellRunoff[actualSingleCellNumber][day]);
            if (options.scoutCellRunoffmm)
                fprintf(file_ptr, "%10.6f\t", G_daily365CellRunoff_mm[actualSingleCellNumber][day]);
            if (options.scoutCellSRunoff)
                fprintf(file_ptr, "%10.6f\t", G_daily365CellSurfaceRunoff[actualSingleCellNumber][day]);
            if (options.scoutQ) fprintf(file_ptr, "%10.6f\t", G_daily365RiverAvail[actualSingleCellNumber][day]);
            if (options.scoutFlowVelo) fprintf(file_ptr, "%10.6f\t", G_daily365Velocity[actualSingleCellNumber][day]);
            if (options.scoutLocLake) fprintf(file_ptr, "%10.6f\t", G_daily365LocLakeStor[actualSingleCellNumber][day]);
            if (options.scoutLocWet) fprintf(file_ptr, "%10.6f\t", G_daily365LocWetlStor[actualSingleCellNumber][day]);
            if (options.scoutGloLake) fprintf(file_ptr, "%10.6f\t", G_daily365GloLakeStor[actualSingleCellNumber][day]);
            if (options.scoutReservoir) fprintf(file_ptr, "%10.6f\t", G_daily365ResStor[actualSingleCellNumber][day]);
            if (options.scoutGloWet) fprintf(file_ptr, "%10.6f\t", G_daily365GloWetlStor[actualSingleCellNumber][day]);
            if (options.scoutRiver) fprintf(file_ptr, "%10.6f\t", G_daily365RiverStor[actualSingleCellNumber][day]);
            if (options.scoutSurfaceStor)
                fprintf(file_ptr, "%10.6f\t", G_daily365SurfStor[actualSingleCellNumber][day]);
            if (options.scoutTWSkm3)
                fprintf(file_ptr, "%10.6f\t", G_daily365TotalWaterInStorages_km3[actualSingleCellNumber][day]);
            if (options.scoutTWSmm)
                fprintf(file_ptr, "%10.6f\t", G_daily365TotalWaterInStorages_mm[actualSingleCellNumber][day]);
            if (options.scoutGwStor) fprintf(file_ptr, "%10.6f\t", G_daily365GwStor[actualSingleCellNumber][day]);
            if (options.scoutCanopyStor)
                fprintf(file_ptr, "%10.6f\t", dailyWaterBalance.G_daily365CanopyStorage[actualSingleCellNumber][day]);
            if (options.scoutSnowStor)
                fprintf(file_ptr, "%10.6f\t", dailyWaterBalance.G_daily365SnowStorage[actualSingleCellNumber][day]);
            if (options.scoutSoilStor)
                fprintf(file_ptr, "%10.6f\t", dailyWaterBalance.G_daily365SoilStorage[actualSingleCellNumber][day]);
            if (options.scoutGwrSwb) fprintf(file_ptr, "%10.6f\t", G_daily365GwrSwb[actualSingleCellNumber][day]);
            if (options.scoutFswb) fprintf(file_ptr, "%10.6f\t", G_daily365Fswb[actualSingleCellNumber][day]);
            if (options.scoutLandAreaFrac)
                fprintf(file_ptr, "%10.6f\t", G_daily365LandAreaFrac[actualSingleCellNumber][day]);
            fprintf(file_ptr, "%3d\n", day + 1);
        }
        fclose(file_ptr);
    }
}
// end storeSingleCellDailyValuesToFile
// HMS 2013-11-21 Storage of single cell values end

void routingClass::writeLakeStorageGrids(char *outputDir, int year) {
    // values in [km3]  // yearly outputs (1 == options.grid_store)
    char filename[250];
    int n;
    double G_array[ng];

    for (n = 0; n < ng; n++)
        G_array[n] = G_locLakeStorage[n]; // conversion double -> float
    sprintf(filename, "%s/G_LOC_LAKE_STORAGE_km3_%d.UNF0", outputDir, year);
    gridIO.writeUnfFile(filename, ng, &G_array[0]);

    for (n = 0; n < ng; n++)
        G_array[n] = G_locWetlStorage[n]; // conversion double -> float
    sprintf(filename, "%s/G_LOC_WETL_STORAGE_km3_%d.UNF0", outputDir, year);
    gridIO.writeUnfFile(filename, ng, &G_array[0]);

    for (n = 0; n < ng; n++)
        G_array[n] = G_gloLakeStorage[n]; // conversion double -> float
    sprintf(filename, "%s/G_GLO_LAKE_STORAGE_km3_%d.UNF0", outputDir, year);
    gridIO.writeUnfFile(filename, ng, &G_array[0]);

    for (n = 0; n < ng; n++)
        G_array[n] = G_gloWetlStorage[n]; // conversion double -> float
    sprintf(filename, "%s/G_GLO_WETL_STORAGE_km3_%d.UNF0", outputDir, year);
    gridIO.writeUnfFile(filename, ng, &G_array[0]);

    // for new reservoir algorithm
    if (options.resOpt == 1) {
        for (n = 0; n < ng; n++)
            G_array[n] = G_gloResStorage[n]; // conversion double -> float
        sprintf(filename, "%s/G_RES_STORAGE_km3_%d.UNF0", outputDir, year);
        gridIO.writeUnfFile(filename, ng, &G_array[0]);
    }


}
// end writeLakeStorageGrids

void routingClass::writeLakeStorRatioGrids(char *outputDir, int year) {
    // range of values: 1 is maximum   // yearly outputs (1 == options.grid_store)
    char filename[250];
    float G_storageRatio[ng];
    int n;

    // local lakes
    for (n = 0; n <= ng - 1; n++) {
        G_storageRatio[n] = G_locLakeStorage[n] / ((G_loc_lake[n] / 100.0)
                                                   * geo.areaOfCellByArrayPos(n) * G_lakeDepthActive[n]);
        // term in brackets is equal to maxStorage
    }
    sprintf(filename, "%s/G_LOC_LAKE_STOR_RATIO_%d.UNF0", outputDir, year);
    gridIO.writeUnfFile(filename, ng, G_storageRatio);

    // local wetlands
    for (int n = 0; n <= ng - 1; n++) {
        // Cell-specific active wetland depth - Use array (initialized from cell-specific parameter)  // FP
        G_storageRatio[n] = G_locWetlStorage[n] / ((G_loc_wetland[n] / 100.0)
                                                   * geo.areaOfCellByArrayPos(n) * G_wetlDepthActive[n]);
    }
    sprintf(filename, "%s/G_LOC_WETL_STOR_RATIO_%d.UNF0", outputDir, year);
    gridIO.writeUnfFile(filename, ng, G_storageRatio);

    // global lakes
    for (int n = 0; n <= ng - 1; n++) {
        G_storageRatio[n] = G_gloLakeStorage[n]
                            // for new reservoir algorithm			/ ((G_lake_area[n] + G_reservoir_area[n])* lakeDepth);
                            // Cell-specific active lake depth - Use array (initialized from cell-specific parameter)  // FP
                            / ((double) G_lake_area[n] * G_lakeDepthActive[n]);

        ///		G_storageRatio[n] = G_gloLakeStorage[n]
        ///			/ ((G_glo_lake[n] / 100.0)
        ///			   * geo.areaOfCellByArrayPos(n) * lakeDepth);

    }
    sprintf(filename, "%s/G_GLO_LAKE_STOR_RATIO_%d.UNF0", outputDir, year);
    gridIO.writeUnfFile(filename, ng, G_storageRatio);

    // global wetlands
    for (int n = 0; n <= ng - 1; n++) {
        // Cell-specific active wetland depth - Use array (initialized from cell-specific parameter)  // FP
        G_storageRatio[n] = G_gloWetlStorage[n] / ((G_glo_wetland[n] / 100.0)
                                                   * geo.areaOfCellByArrayPos(n) * G_wetlDepthActive[n]);
    }
    sprintf(filename, "%s/G_GLO_WETL_STOR_RATIO_%d.UNF0", outputDir, year);
    gridIO.writeUnfFile(filename, ng, G_storageRatio);

    // for new reservoir algorithm
    // reservoirs
    if (options.resOpt == 1) {
        for (int n = 0; n <= ng - 1; n++) {
            G_storageRatio[n] = G_gloResStorage[n] /
                                (G_stor_cap[n] * 0.85); //85 % of storage capacity (from published volume data)
        }
        sprintf(filename, "%s/G_RES_STOR_RATIO_%d.UNF0", outputDir, year);
        gridIO.writeUnfFile(filename, ng, G_storageRatio);
    }
}
// end writeLakeStorRatioGrids

void routingClass::writeMonthlyGrids(char *outputDir, int year) {
    if (2 != options.grid_store && 4 != options.grid_store && 6 != options.grid_store) return;

    char filename[250];
    int n;
    short m;
    double G_monthlyArray[ng][12]; // used to write float array to file
    double numberOfDaysInMonth[12] = {31., 28., 31., 30., 31., 30., 31., 31., 30., 31., 30., 31.};

    if (options.outRiverAvail) { // new output options 2.1f
        for (n = 0; n < ng; n++)
            for (m = 0; m < 12; m++)
                G_monthlyArray[n][m] = G_monthlyRiverAvail[n][m]; // conversion double -> float
        sprintf(filename, "%s/G_RIVER_AVAIL_%d.12.UNF0", outputDir, year);
        gridIO.writeUnfFile(filename, ng * 12, &G_monthlyArray[0][0]);
    }
    // endif options.outRiverAvail

    if (options.outPrec) { // new output options 2.2 Precip which is really used by WaterGAP
        for (n = 0; n < ng; n++)
            for (m = 0; m < 12; m++)
                G_monthlyArray[n][m] = G_monthlyConsistentPrecip[n][m];
        sprintf(filename, "%s/G_CONSISTENT_PRECIPITATION_km3_%d.12.UNF0", outputDir, year);
        gridIO.writeUnfFile(filename, ng * 12, &G_monthlyArray[0][0]);
    }
    // endif options.outPrec

    //new output option min/max daily discharge
    if (options.outMinMaxRiverAvail) {
        for (n = 0; n < ng; n++)
            for (m = 0; m < 12; m++)
                G_monthlyArray[n][m] = G_monthlyMinRiverAvail[n][m]; // conversion double -> float
        sprintf(filename, "%s/G_MIN_RIVER_AVAIL_%d.12.UNF0", outputDir, year);
        gridIO.writeUnfFile(filename, ng * 12, &G_monthlyArray[0][0]);

        for (n = 0; n < ng; n++)
            for (m = 0; m < 12; m++)
                G_monthlyArray[n][m] = G_monthlyMaxRiverAvail[n][m]; // conversion double -> float
        sprintf(filename, "%s/G_MAX_RIVER_AVAIL_%d.12.UNF0", outputDir, year);
        gridIO.writeUnfFile(filename, ng * 12, &G_monthlyArray[0][0]);
    }
    // endif options.outMinMaxRiverAvail

    if (options.outRiverPET) {
        for (n = 0; n < ng; n++)
            for (m = 0; m < 12; m++)
                G_monthlyArray[n][m] =
                               G_monthlyRiverAreaFrac[n][m] / numberOfDaysInMonth[m]; // conversion double -> float
        sprintf(filename, "%s/G_RIVER_FRACTIONS_%d.12.UNF0", outputDir, year);
        gridIO.writeUnfFile(filename, ng * 12, &G_monthlyArray[0][0]);

        for (n = 0; n < ng; n++)
            for (m = 0; m < 12; m++)
                G_monthlyArray[n][m] = G_monthlyRiverPET[n][m]; // conversion double -> float
        sprintf(filename, "%s/G_RIVER_PET_%d.12.UNF0", outputDir, year);
        gridIO.writeUnfFile(filename, ng * 12, &G_monthlyArray[0][0]);

    }
    // endif options.outRiverPET

    if (options.outRiverInUpstream) { // new output options 2.2
        for (n = 0; n < ng; n++)
            for (m = 0; m < 12; m++)
                G_monthlyArray[n][m] = G_monthlyRiverInUpstream[n][m]; // conversion double -> float
        sprintf(filename, "%s/G_RIVER_IN_UPSTREAM_%d.12.UNF0", outputDir, year);
        gridIO.writeUnfFile(filename, ng * 12, &G_monthlyArray[0][0]);
    }
    // endif options.outRiverInUpstream

    if (options.outRiverVelo) { // new output options 2.2
        for (n = 0; n < ng; n++)
            for (m = 0; m < 12; m++)
                G_monthlyArray[n][m] = G_monthlyVelocity[n][m] / numberOfDaysInMonth[m]; // conversion double -> float
        sprintf(filename, "%s/G_RIVER_VELOCITY_%d.12.UNF0", outputDir, year);
        gridIO.writeUnfFile(filename, ng * 12, &G_monthlyArray[0][0]);
    }
    // endif options.outCellSurface

    if (options.outLocWetlExt) {
        for (n = 0; n < ng; n++)
            for (m = 0; m < 12; m++)
                G_monthlyArray[n][m] = G_monthlyLocWetlExtent[n][m] / numberOfDaysInMonth[m];
        sprintf(filename, "%s/G_LOC_WETL_EXTENT_km2_%d.12.UNF0", outputDir, year);
        gridIO.writeUnfFile(filename, ng * 12, &G_monthlyArray[0][0]);
    }
    if (options.outGloWetlExt) {
        for (n = 0; n < ng; n++)
            for (m = 0; m < 12; m++)
                G_monthlyArray[n][m] = G_monthlyGloWetlExtent[n][m] / numberOfDaysInMonth[m];
        sprintf(filename, "%s/G_GLO_WETL_EXTENT_km2_%d.12.UNF0", outputDir, year);
        gridIO.writeUnfFile(filename, ng * 12, &G_monthlyArray[0][0]);
    }
    if (options.outCellRunoff) { // new output options 2.2
        for (n = 0; n < ng; n++)
            for (m = 0; m < 12; m++)
                G_monthlyArray[n][m] = G_monthlyCellRunoff[n][m]; // conversion double -> float
        sprintf(filename, "%s/G_CELL_RUNOFF_%d.12.UNF0", outputDir, year);
        gridIO.writeUnfFile(filename, ng * 12, &G_monthlyArray[0][0]);

        for (n = 0; n < ng; n++) { // total water resources in mm/month
            for (m = 0; m < 12; m++) {
                // FP 2015-06 // Only continental cells
                if (geo.G_contcell[n]) {
                    G_monthlyArray[n][m] = G_monthlyCellRunoff[n][m] // conversion double -> float
                                           / (geo.areaOfCellByArrayPos(n) * (geo.G_contfreq[n] / 100.)) * 1000000;
                }
                    // 100% water cells (of Caspian Sea) have now G_fwaterfreq + G_landfreq = 0. The water balance is calculated in the outlet cell of large lakes; therefore a calculation of e.g. precipitation for this cells leads to inconsistent quantification of global fluxes. A division through 0 has to be avoided now (alternative: define it as not land cells in ng).
                else {
                    G_monthlyArray[n][m] = 0.;
                }
            }
            //endfor loop over months
        }
        // endfor loop over grid cells
        sprintf(filename, "%s/G_CELL_RUNOFF_mm_%d.12.UNF0", outputDir, year);
        gridIO.writeUnfFile(filename, ng * 12, &G_monthlyArray[0][0]);

        //if (options.outCellRunoffTotal) { // new output options 2.2b to be implemented as output option
        // dynamic memory allocation
        // static allocation has lead to 'segmentation fault'
        double (*const G_totalCellRunoff)[12] = new double[ng][12];
        for (n = 0; n < ng; n++)
            for (m = 0; m < 12; m++) {
                G_monthlyPotCellRunoff[n][m] -= G_monthlyPotCellRunoffDeficit[n][m];    // contains negative values
                G_monthlyArray[n][m] = G_monthlyPotCellRunoff[n][m]; // conversion double -> float
            }
        sprintf(filename, "%s/G_CELL_RUNOFF_TOTAL_%d.12.UNF0", outputDir, year);
        gridIO.writeUnfFile(filename, ng * 12, &G_monthlyArray[0][0]);
        delete[]G_totalCellRunoff;
        //}

    }
    // endif options.outCellRunoff

    if (options.outCellSurface) { // new output options 2.2
        for (n = 0; n < ng; n++)
            for (m = 0; m < 12; m++)
                G_monthlyArray[n][m] = G_monthlyCellSurfaceRunoff[n][m]; // conversion double -> float
        sprintf(filename, "%s/G_CELL_SURFACE_RUNOFF_%d.12.UNF0", outputDir, year);
        gridIO.writeUnfFile(filename, ng * 12, &G_monthlyArray[0][0]);
    }
    // endif options.outCellSurface

    if (options.outSurfStor) { // new output options 2.2

        if (options.grid_store_TypeForStorages == 0) {

            for (n = 0; n < ng; n++)
                for (m = 0; m < 12; m++)
                    G_monthlyArray[n][m] = G_monthlySurfStor[n][m]; // / numberOfDaysInMonth[m]; // conversion double -> float
            sprintf(filename, "%s/G_SURFACE_WATER_STORAGE_km3_%d.12.UNF0", outputDir, year);
            gridIO.writeUnfFile(filename, ng * 12, &G_monthlyArray[0][0]);

            for (n = 0; n < ng; n++) {
                for (m = 0; m < 12; m++) {
                    // FP 2015-06 // Only continental cells
                    if (geo.G_contcell[n]) {
                        G_monthlyArray[n][m] =
                                       (G_monthlySurfStor[n][m]// / numberOfDaysInMonth[m]; // conversion double -> float
                                        / (geo.areaOfCellByArrayPos(n) * geo.G_contfreq[n] / 100.)) * 1000000;
                    }
                        // 100% water cells (of Caspian Sea) have now G_fwaterfreq + G_landfreq = 0. The water balance is calculated in the outlet cell of large lakes; therefore a calculation of e.g. precipitation for this cells leads to inconsistent quantification of global fluxes. A division through 0 has to be avoided now (alternative: define it as not land cells in ng).
                    else {
                        G_monthlyArray[n][m] = 0.;
                    }
                }
                //endfor loop over months
            }
            // endfor loop over grid cells
            sprintf(filename, "%s/G_SURFACE_WATER_STORAGE_mm_%d.12.UNF0", outputDir, year);
            gridIO.writeUnfFile(filename, ng * 12, &G_monthlyArray[0][0]);

        }
            // endif (options.grid_store_TypeForStorages == 0)
        else {
            for (n = 0; n < ng; n++)
                for (m = 0; m < 12; m++)
                    G_monthlyArray[n][m] =
                                   G_monthlySurfStor[n][m] / numberOfDaysInMonth[m]; // conversion double -> float
            sprintf(filename, "%s/G_SURFACE_WATER_STORAGE_MEAN_km3_%d.12.UNF0", outputDir, year);
            gridIO.writeUnfFile(filename, ng * 12, &G_monthlyArray[0][0]);

            for (n = 0; n < ng; n++) {
                for (m = 0; m < 12; m++) {
                    // FP 2015-06 // Only continental cells
                    if (geo.G_contcell[n]) {
                        G_monthlyArray[n][m] =
                                       (G_monthlySurfStor[n][m] / numberOfDaysInMonth[m] // conversion double -> float
                                        / (geo.areaOfCellByArrayPos(n) * geo.G_contfreq[n] / 100.)) * 1000000;
                    }
                        // 100% water cells (of Caspian Sea) have now G_fwaterfreq + G_landfreq = 0. The water balance is calculated in the outlet cell of large lakes; therefore a calculation of e.g. precipitation for this cells leads to inconsistent quantification of global fluxes. A division through 0 has to be avoided now (alternative: define it as not land cells in ng).
                    else {
                        G_monthlyArray[n][m] = 0.;
                    }
                }
                //endfor loop over months
            }
            // endfor loop over grid cells
            sprintf(filename, "%s/G_SURFACE_WATER_STORAGE_MEAN_mm_%d.12.UNF0", outputDir, year);
            gridIO.writeUnfFile(filename, ng * 12, &G_monthlyArray[0][0]);

        }
        // endelse (i.e. options.grid_store_TypeForStorages != 0)

    }
    // endif options.outSurfStor

    if (options.outCellAET) { // new output options 2.1f
        for (n = 0; n < ng; n++)
            for (m = 0; m < 12; m++)
                G_monthlyArray[n][m] = G_monthlyCellAET[n][m]; // conversion double -> float
        sprintf(filename, "%s/G_CELL_AET_%d.12.UNF0", outputDir, year);
        gridIO.writeUnfFile(filename, ng * 12, &G_monthlyArray[0][0]);

    }
    // endif options.outCellAET


    if (options.outOpenWaterEvap) { // new output options 2.2
        for (n = 0; n < ng; n++)
            for (m = 0; m < 12; m++)
                G_monthlyArray[n][m] = G_monthlyOpenWaterEvap[n][m]; // conversion double -> float
        sprintf(filename, "%s/G_OPEN_WATER_EVAP_%d.12.UNF0", outputDir, year);
        gridIO.writeUnfFile(filename, ng * 12, &G_monthlyArray[0][0]);

    }
    // endif options.outOpenWaterEvap

    if (options.outGwrSwb) { // monthly sum of groundwater recharge below surface water bodies (with respect to continental area) [mm/month]
        for (n = 0; n < ng; n++)
            for (m = 0; m < 12; m++)
                G_monthlyArray[n][m] = G_monthlyGwrSwb[n][m];
        sprintf(filename, "%s/G_GWR_SURFACE_WATER_BODIES_mm_%d.12.UNF0", outputDir, year);
        gridIO.writeUnfFile(filename, ng * 12, &G_monthlyArray[0][0]);

        for (n = 0; n < ng; n++)
            for (short m = 0; m < 12; m++) {
                G_monthlyArray[n][m] = G_monthlyGwrSwb[n][m]
                                       * (geo.areaOfCellByArrayPos(n) * (geo.G_contfreq[n] / 100.)) / 1000000.;
            }
        sprintf(filename, "%s/G_GWR_SURFACE_WATER_BODIES_km3_%d.12.UNF0", outputDir, year);
        gridIO.writeUnfFile(filename, ng * 12, &G_monthlyArray[0][0]);
    }
    // endif options.outGwrSwb

    if (options.outTotalGWR) {
        for (n = 0; n < ng; n++)
            for (m = 0; m < 12; m++)
                G_monthlyArray[n][m] = G_monthlyGwrSwb[n][m] + dailyWaterBalance.G_monthlyGwRecharge[n][m];
        sprintf(filename, "%s/G_TOTAL_GW_RECHARGE_mm_%d.12.UNF0", outputDir, year);
        gridIO.writeUnfFile(filename, ng * 12, &G_monthlyArray[0][0]);

        for (n = 0; n < ng; n++)
            for (m = 0; m < 12; m++)
                G_monthlyArray[n][m] = (G_monthlyGwrSwb[n][m] + dailyWaterBalance.G_monthlyGwRecharge[n][m])
                                       * (geo.areaOfCellByArrayPos(n) * (geo.G_contfreq[n] / 100.)) / 1000000.;
        sprintf(filename, "%s/G_TOTAL_GW_RECHARGE_km3_%d.12.UNF0", outputDir, year);
        gridIO.writeUnfFile(filename, ng * 12, &G_monthlyArray[0][0]);
    }

    if (options.outFswb) { // mean monthly fraction of surface water bodies
        for (n = 0; n < ng; n++)
            for (m = 0; m < 12; m++) {
                G_monthlyArray[n][m] = G_monthlyFswb[n][m] / numberOfDaysInMonth[m];
            }
        sprintf(filename, "%s/G_FRACTION_SURFACE_WATER_BODIES_%d.12.UNF0", outputDir, year);
        gridIO.writeUnfFile(filename, ng * 12, &G_monthlyArray[0][0]);
    }
    // endif options.outFswb
    if (options.outLandAreaFrac) {
        // Write land area fractions
        for (n = 0; n < ng; n++)
            for (m = 0; m < 12; m++)
                G_monthlyArray[n][m] = G_monthlyLandAreaFrac[n][m] / numberOfDaysInMonth[m];
        sprintf(filename, "%s/G_LAND_AREA_FRACTIONS_%d.12.UNF0", outputDir, year);
        gridIO.writeUnfFile(filename, ng * 12, &G_monthlyArray[0][0]);
    }

    // GFZ:GLOLAK begin - Allocation of GloLakeStorage of Outflowcell to all Lakecells
    //for lakes and reservoirs
    if (0 != ng_glolakcells) {

        float G_tempGloLakeStorage[ng_glolakcells];
        //float G_tempResStorage[ng_glolakcells];

        long ofcell, acell;
        for (short month = 0; month <= 11; month++) {
            //Initialisation of the temporary Glolakstorage-array
            for (n = 0; n < ng_glolakcells; n++) {
                G_tempGloLakeStorage[n] = 0.0;
                //G_tempResStorage[n] = 0; //added for new reservoir storage
            }
            //computation of glolakstorage for each lakecell and saving in temporary array
            for (n = 0; n < ng_glolakcells; n++) {
                ofcell = big_lakes_cells[0][n];
                acell = big_lakes_cells[1][n];
                if (ofcell != -99)
                    G_tempGloLakeStorage[n] += G_monthlyGloLakeStorage[acell - 1][month] * big_lakes_area_frac[n];
                //G_tempResStorage[n] += 	G_monthlyResStorage[acell-1][month] * big_lakes_area_frac[n];
            }
            //Rewriting of global array for Glolakstorage
            for (n = 0; n < ng_glolakcells; n++) {
                ofcell = big_lakes_cells[0][n];
                G_monthlyGloLakeStorage[ofcell - 1][month] = G_tempGloLakeStorage[n];
                //G_monthlyResStorage[ofcell-1][month] = G_tempResStorage[n];
            }
        }
    }
    // endif (0 != ng_glolakcells)
    //GFZ:GLOLAK end

    // monthly storage grids
    if (options.grid_store_TypeForStorages == 1) {
        for (n = 0; n <= ng - 1; n++) {
            for (short month = 0; month <= 11; month++) {
                G_monthlyGwStorage[n][month] /= numberOfDaysInMonth[month]; // moved from daily.cpp
                G_monthlyLocLakeStorage[n][month] /= numberOfDaysInMonth[month];
                G_monthlyLocWetlStorage[n][month] /= numberOfDaysInMonth[month];
                G_monthlyGloLakeStorage[n][month] /= numberOfDaysInMonth[month];
                G_monthlyGloWetlStorage[n][month] /= numberOfDaysInMonth[month];
                G_monthlyRiverStorage[n][month] /= numberOfDaysInMonth[month];
                if (options.resOpt == 1)
                    G_monthlyResStorage[n][month] /= numberOfDaysInMonth[month];
            }
        } //for(n)
    } // endif options.grid_store_TypeForStorages

    if (options.outGWRunoff) { // moved from daily.cpp
        for (n = 0; n < ng; n++)
            for (m = 0; m < 12; m++)
                G_monthlyArray[n][m] = G_monthlyGwRunoff[n][m];
        sprintf(filename, "%s/G_GW_RUNOFF_km3_%d.12.UNF0", outputDir, year);
        gridIO.writeUnfFile(filename, ng * 12, &G_monthlyArray[0][0]);

        for (n = 0; n < ng; n++) { // for output in mm
            for (m = 0; m < 12; m++) {
                // FP 2015-06 // Only continental cells
                if (geo.G_contcell[n]) {
                    G_monthlyArray[n][m] = G_monthlyGwRunoff[n][m] // conversion double -> float
                                           / (geo.areaOfCellByArrayPos(n) * (geo.G_contfreq[n] / 100.)) * 1000000.0;
                }
                    // 100% water cells (of Caspian Sea) have now G_fwaterfreq + G_landfreq = 0. The water balance is calculated in the outlet cell of large lakes; therefore a calculation of e.g. precipitation for this cells leads to inconsistent quantification of global fluxes. A division through 0 has to be avoided now (alternative: define it as not land cells in ng).
                else {
                    G_monthlyArray[n][m] = 0.;
                }
            }
            //endfor loop over months
        }
        // endfor loop over grid cells
        sprintf(filename, "%s/G_GW_RUNOFF_mm_%d.12.UNF0", outputDir, year);
        gridIO.writeUnfFile(filename, ng * 12, &G_monthlyArray[0][0]);
    }
    // endif options.outGWRunoff
    //write output of monthly storage (components)

    if (options.outLocLakeStorage || options.outSingleStorages) { // moved from daily.cpp
        for (n = 0; n < ng; n++)
            for (m = 0; m < 12; m++)
                G_monthlyArray[n][m] = G_monthlyGwStorage[n][m];
        sprintf(filename, "%s/G_GROUND_WATER_STORAGE_km3_%d.12.UNF0", outputDir, year);
        gridIO.writeUnfFile(filename, ng * 12, &G_monthlyArray[0][0]);


        for (n = 0; n < ng; n++) { // for output in mm
            for (m = 0; m < 12; m++) {
                // FP 2015-06 // Only continental cells
                if (geo.G_contcell[n]) {
                    G_monthlyArray[n][m] =
                                   G_monthlyGwStorage[n][m]// / numberOfDaysInMonth[m]; // conversion double -> float
                                   / (geo.areaOfCellByArrayPos(n) * (geo.G_contfreq[n] / 100.)) * 1000000.0;
                }
                    // 100% water cells (of Caspian Sea) have now G_fwaterfreq + G_landfreq = 0. The water balance is calculated in the outlet cell of large lakes; therefore a calculation of e.g. precipitation for this cells leads to inconsistent quantification of global fluxes. A division through 0 has to be avoided now (alternative: define it as not land cells in ng).
                else {
                    G_monthlyArray[n][m] = 0.;
                }
            }
            //endfor loop over months
        }
        // endfor loop over grid cells

        sprintf(filename, "%s/G_GROUND_WATER_STORAGE_mm_%d.12.UNF0", outputDir, year);
        gridIO.writeUnfFile(filename, ng * 12, &G_monthlyArray[0][0]);

    }
    // endif options
    if (options.outLocLakeStorage || options.outSingleStorages) { // neue OUTPUT-Option
        for (n = 0; n < ng; n++)
            for (m = 0; m < 12; m++)
                G_monthlyArray[n][m] = G_monthlyLocLakeStorage[n][m]; // LOCAL LAKES - conversion double -> float
        sprintf(filename, "%s/G_LOC_LAKE_STORAGE_km3_%d.12.UNF0", outputDir, year);
        gridIO.writeUnfFile(filename, ng * 12, &G_monthlyArray[0][0]);

        for (n = 0; n < ng; n++) { // for output in mm
            for (m = 0; m < 12; m++) {
                // FP 2015-06 // Only continental cells
                if (geo.G_contcell[n]) {
                    G_monthlyArray[n][m] = G_monthlyLocLakeStorage[n][m] // LOCAL LAKES - conversion double -> float
                                           / (geo.areaOfCellByArrayPos(n) * (geo.G_contfreq[n] / 100.)) * 1000000;
                }
                    // 100% water cells (of Caspian Sea) have now G_fwaterfreq + G_landfreq = 0. The water balance is calculated in the outlet cell of large lakes; therefore a calculation of e.g. precipitation for this cells leads to inconsistent quantification of global fluxes. A division through 0 has to be avoided now (alternative: define it as not land cells in ng).
                else {
                    G_monthlyArray[n][m] = 0.;
                }
            }
            //endfor loop over months
        }
        // endfor loop over grid cells
        sprintf(filename, "%s/G_LOC_LAKE_STORAGE_mm_%d.12.UNF0", outputDir, year);
        gridIO.writeUnfFile(filename, ng * 12, &G_monthlyArray[0][0]);
    }
    // endif options

    if (options.outLocWetlStorage || options.outSingleStorages) { // neue OUTPUT-Option
        for (n = 0; n < ng; n++)
            for (m = 0; m < 12; m++)
                G_monthlyArray[n][m] = G_monthlyLocWetlStorage[n][m]; // LOCAL WETLANDS - conversion double -> float
        sprintf(filename, "%s/G_LOC_WETL_STORAGE_km3_%d.12.UNF0", outputDir, year);
        gridIO.writeUnfFile(filename, ng * 12, &G_monthlyArray[0][0]);

        for (n = 0; n < ng; n++) { // for output in mm
            for (m = 0; m < 12; m++) {
                // FP 2015-06 // Only continental cells
                if (geo.G_contcell[n]) {
                    G_monthlyArray[n][m] = G_monthlyLocWetlStorage[n][m] // LOCAL WETLANDS - conversion double -> float
                                           / (geo.areaOfCellByArrayPos(n) * (geo.G_contfreq[n] / 100.)) * 1000000;
                }
                    // 100% water cells (of Caspian Sea) have now G_fwaterfreq + G_landfreq = 0. The water balance is calculated in the outlet cell of large lakes; therefore a calculation of e.g. precipitation for this cells leads to inconsistent quantification of global fluxes. A division through 0 has to be avoided now (alternative: define it as not land cells in ng).
                else {
                    G_monthlyArray[n][m] = 0.;
                }
            }
            //endfor loop over months
        }
        // endfor loop over grid cells
        sprintf(filename, "%s/G_LOC_WETL_STORAGE_mm_%d.12.UNF0", outputDir, year);
        gridIO.writeUnfFile(filename, ng * 12, &G_monthlyArray[0][0]);
    }
    // endif options

    if (options.outGloLakeStorage || options.outSingleStorages) { // neue OUTPUT-Option
        for (n = 0; n < ng; n++)
            for (m = 0; m < 12; m++)
                G_monthlyArray[n][m] = G_monthlyGloLakeStorage[n][m]; // GLOBAL LAKES - conversion double -> float
        sprintf(filename, "%s/G_GLO_LAKE_STORAGE_km3_%d.12.UNF0", outputDir, year);
        gridIO.writeUnfFile(filename, ng * 12, &G_monthlyArray[0][0]);

        for (n = 0; n < ng; n++) { // for output in mm
            for (m = 0; m < 12; m++) {
                // FP 2015-06 // Only continental cells
                if (geo.G_contcell[n]) {
                    G_monthlyArray[n][m] = G_monthlyGloLakeStorage[n][m] // GLOBAL LAKES - conversion double -> float
                                           / (geo.areaOfCellByArrayPos(n) * (geo.G_contfreq[n] / 100.)) * 1000000;
                }
                    // 100% water cells (of Caspian Sea) have now G_fwaterfreq + G_landfreq = 0. The water balance is calculated in the outlet cell of large lakes; therefore a calculation of e.g. precipitation for this cells leads to inconsistent quantification of global fluxes. A division through 0 has to be avoided now (alternative: define it as not land cells in ng).
                else {
                    G_monthlyArray[n][m] = 0.;
                }
            }
            //endfor loop over months
        }
        // endfor loop over grid cells
        sprintf(filename, "%s/G_GLO_LAKE_STORAGE_mm_%d.12.UNF0", outputDir, year);
        gridIO.writeUnfFile(filename, ng * 12, &G_monthlyArray[0][0]);
    }
    // endif option

    if ((options.outResStorage || options.outSingleStorages) && options.resOpt == 1) { // neue OUTPUT-Option
        for (n = 0; n < ng; n++)
            for (m = 0; m < 12; m++)
                G_monthlyArray[n][m] = G_monthlyResStorage[n][m]; // RESERVOIRS - conversion double -> float

        if (options.grid_store_TypeForStorages == 0)
            sprintf(filename, "%s/G_RES_STORAGE_km3_%d.12.UNF0", outputDir, year);
        else
            sprintf(filename, "%s/G_RES_STORAGE_MEAN_km3_%d.12.UNF0", outputDir, year);

        gridIO.writeUnfFile(filename, ng * 12, &G_monthlyArray[0][0]);


        for (n = 0; n < ng; n++) { // for output in mm
            for (m = 0; m < 12; m++) {
                // FP 2015-06 // Only continental cells
                if (geo.G_contcell[n]) {
                    G_monthlyArray[n][m] = G_monthlyResStorage[n][m] // RESERVOIRS - conversion double -> float
                                           / (geo.areaOfCellByArrayPos(n) * (geo.G_contfreq[n] / 100.)) * 1000000;
                }
                    // 100% water cells (of Caspian Sea) have now G_fwaterfreq + G_landfreq = 0. The water balance is calculated in the outlet cell of large lakes; therefore a calculation of e.g. precipitation for this cells leads to inconsistent quantification of global fluxes. A division through 0 has to be avoided now (alternative: define it as not land cells in ng).
                else {
                    G_monthlyArray[n][m] = 0.;
                }
            }
            //endfor loop over months
        }
        // endfor loop over grid cells

        if (options.grid_store_TypeForStorages == 0)
            sprintf(filename, "%s/G_RES_STORAGE_mm_%d.12.UNF0", outputDir, year);
        else
            sprintf(filename, "%s/G_RES_STORAGE_MEAN_mm_%d.12.UNF0", outputDir, year);

        gridIO.writeUnfFile(filename, ng * 12, &G_monthlyArray[0][0]);

    }
    // endif options

    if ((options.outResStorage || options.outSingleStorages) && options.resOpt == 1) { // neue OUTPUT-Option
        for (n = 0; n < ng; n++)
            for (m = 0; m < 12; m++)
                G_monthlyArray[n][m] = G_monthlyResOutflow[n][m] /
                                       1000000000.; // RESERVOIRS - conversion double -> float; [m3/month] -> [km3/month]
        sprintf(filename, "%s/G_RES_OUT_km3_%d.12.UNF0", outputDir, year);
        gridIO.writeUnfFile(filename, ng * 12, &G_monthlyArray[0][0]);
    }
    // endif options

    if (options.outGloWetlStorage || options.outSingleStorages) { // neue OUTPUT-Option
        for (n = 0; n < ng; n++)
            for (m = 0; m < 12; m++)
                G_monthlyArray[n][m] = G_monthlyGloWetlStorage[n][m]; // GLOBAL WETLANDS - conversion double -> float
        sprintf(filename, "%s/G_GLO_WETL_STORAGE_km3_%d.12.UNF0", outputDir, year);
        gridIO.writeUnfFile(filename, ng * 12, &G_monthlyArray[0][0]);

        for (n = 0; n < ng; n++) { // for output in mm
            for (m = 0; m < 12; m++) {
                // FP 2015-06 // Only continental cells
                if (geo.G_contcell[n]) {
                    G_monthlyArray[n][m] = G_monthlyGloWetlStorage[n][m] // GLOBAL WETLANDS - conversion double -> float
                                           / (geo.areaOfCellByArrayPos(n) * (geo.G_contfreq[n] / 100.)) * 1000000;
                }
                    // 100% water cells (of Caspian Sea) have now G_fwaterfreq + G_landfreq = 0. The water balance is calculated in the outlet cell of large lakes; therefore a calculation of e.g. precipitation for this cells leads to inconsistent quantification of global fluxes. A division through 0 has to be avoided now (alternative: define it as not land cells in ng).
                else {
                    G_monthlyArray[n][m] = 0.;
                }
            }
            //endfor loop over months
        }
        // endfor loop over grid cells
        sprintf(filename, "%s/G_GLO_WETL_STORAGE_mm_%d.12.UNF0", outputDir, year);
        gridIO.writeUnfFile(filename, ng * 12, &G_monthlyArray[0][0]);

    }
    // endif options

    if (options.outRiverStorage || options.outSingleStorages) { // neue OUTPUT-Option
        for (n = 0; n < ng; n++)
            for (m = 0; m < 12; m++)
                G_monthlyArray[n][m] = G_monthlyRiverStorage[n][m]; // RIVER STORAGE - conversion double -> float
        sprintf(filename, "%s/G_RIVER_STORAGE_km3_%d.12.UNF0", outputDir, year);
        gridIO.writeUnfFile(filename, ng * 12, &G_monthlyArray[0][0]);

        for (n = 0; n < ng; n++) { // for output in mm
            for (m = 0; m < 12; m++) {
                // FP 2015-06 // Only continental cells
                if (geo.G_contcell[n]) {
                    G_monthlyArray[n][m] = G_monthlyRiverStorage[n][m] // RIVER STORAGE - conversion double -> float
                                           / (geo.areaOfCellByArrayPos(n) * (geo.G_contfreq[n] / 100.)) * 1000000;
                }
                    // 100% water cells (of Caspian Sea) have now G_fwaterfreq + G_landfreq = 0. The water balance is calculated in the outlet cell of large lakes; therefore a calculation of e.g. precipitation for this cells leads to inconsistent quantification of global fluxes. A division through 0 has to be avoided now (alternative: define it as not land cells in ng).
                else {
                    G_monthlyArray[n][m] = 0.;
                }
            }
            //endfor loop over months
        }
        // endfor loop over grid cells
        sprintf(filename, "%s/G_RIVER_STORAGE_mm_%d.12.UNF0", outputDir, year);
        gridIO.writeUnfFile(filename, ng * 12, &G_monthlyArray[0][0]);

    }
    // endif options


    if (options.outActualNAS) { // neue OUTPUT-Option
        for (n = 0; n < ng; n++)
            for (m = 0; m < 12; m++) {
                //G_monthlyActualUse[n][m] = G_monthlyActualUse[n][m] * 1000000000; // conversion km3/month --> m3/month to avoid loss of small numbers during casting operation.
                G_monthlyArray[n][m] = G_monthlyActualUse[n][m] * 1000000000; // ACTUAL USE - conversion double -> float
            }
        sprintf(filename, "%s/G_ACTUAL_NAS_%d.12.UNF0", outputDir, year);
        gridIO.writeUnfFile(filename, ng * 12, &G_monthlyArray[0][0]);
    }

    if (options.outActualNAG) {
        for (n = 0; n < ng; n++)
            for (m = 0; m < 12; m++) {
                G_monthlyArray[n][m] = G_monthlyNUg[n][m] * 1000000000; // ACTUAL USE - conversion double -> float
            }
        sprintf(filename, "%s/G_ACTUAL_NAG_%d.12.UNF0", outputDir, year);
        gridIO.writeUnfFile(filename, ng * 12, &G_monthlyArray[0][0]);

    }

    if (options.outSatisAllocUsein2ndCell) {
        for (n = 0; n < ng; n++)
            for (m = 0; m < 12; m++) {
            G_monthlyArray[n][m] = G_monthlySatisAllocatedUseinSecondCell[n][m] ;
        }
        sprintf(filename, "%s/G_SATIS_ALLOC_USE_IN_2NDCELL_%d.12.UNF0", outputDir, year);
        gridIO.writeUnfFile(filename, ng * 12, &G_monthlyArray[0][0]);
    }

    if (options.outWCa) {
        for (n = 0; n < ng; n++) {
            for (m = 0; m < 12; m++) {
                G_monthlyArray[n][m] = (G_monthlyActualUse[n][m] + G_monthlyNUg[n][m]) * 1000000000;
            }
        }
        sprintf(filename, "%s/G_ACTUAL_WATER_CONSUMPTION_%d.12.UNF0", outputDir, year);
        gridIO.writeUnfFile(filename, ng * 12, &G_monthlyArray[0][0]);
    }

    if (options.outWCaInCellUsed) {
        for (n = 0; n < ng; n++) {
            for (m = 0; m < 12; m++) {
                G_monthlyArray[n][m] = (G_monthlyRedistributeAllocUse[n][m] + G_monthlyActualUse[n][m] + G_monthlyNUg[n][m]) * 1000000000;
            }
        }
        sprintf(filename, "%s/G_ACTUAL_WATER_CONSUMPTION_INCELLUSED_%d.12.UNF0", outputDir, year);
        gridIO.writeUnfFile(filename, ng * 12, &G_monthlyArray[0][0]);
    }

    if (options.outCellAETWCa) { // new output options 2.2b [km3/month]
        //    if (options.outActualNAS) { // if set to 1, G_monthlyActualUse already in m3/month (see above)
        //         for (n = 0; n < ng; n++) {
        //                 for (m = 0; m < 12; m++) {
        //                     G_monthlyArray[n][m] = G_monthlyCellAET[n][m] + G_monthlyActualUse[n][m] / 1000000000 + G_monthlyNUg[n][m] ; // conversion double -> float
        //                 }
        //         }
        //     } else { // G_monthlyActualUse is in km3/month
        for (n = 0; n < ng; n++) {
            for (m = 0; m < 12; m++) {
                G_monthlyArray[n][m] = G_monthlyCellAETWCa[n][m] + G_monthlyActualUse[n][m] +
                                       G_monthlyNUg[n][m]; // conversion double -> float
            }
        }
        //     }
        sprintf(filename, "%s/G_CELLAET_CONSUSE_km3_%d.12.UNF0", outputDir, year);
        gridIO.writeUnfFile(filename, ng * 12, &G_monthlyArray[0][0]);
    }
// endif options.outCellAET

    // FP 2015, not used anymore, at the moment, implement if needed again
    // CR 15-08-19 G_AllocatedUse = test output: demand allocated from neighboring cell(s)
    if (options.outAllocUsein2ndCell) { // new output options 2.1f
        for (n = 0; n < ng; n++)
            for (m = 0; m < 12; m++)
                G_monthlyArray[n][m] = G_monthlyAllocatedUse[n][m]; // conversion double -> float
        sprintf(filename, "%s/G_ALLOC_USE_IN_2NDCELL_%d.12.UNF0", outputDir, year);
        gridIO.writeUnfFile(filename, ng * 12, &G_monthlyArray[0][0]);
    }
    // endif options.outAllocUse

}
// end writeMonthlyGrids

// daily output option 31 start
void routingClass::writeDaily31Grids(char *outputDir, const int year, const int month) {
    if (3 != options.grid_store && 4 != options.grid_store) return;

    char filename[250];
    double G_dailyArray[ng][31];

    if (options.outCellRunoffDaily) { // new output options 2.2
        for (int i = 0; i < ng; ++i) {
            for (int j = 0; j < 31; ++j)
                G_dailyArray[i][j] = G_daily31CellRunoff[i][j];
        }
        sprintf(filename, "%s/G_CELL_RUNOFF_%d_%d.31.UNF0", outputDir, year, month + 1);
        gridIO.writeUnfFile(filename, ng * 31, &G_dailyArray[0][0]);
    }
    if (options.outGWRunoffDaily) { // new output options 2.2
        for (int i = 0; i < ng; ++i) {
            for (int j = 0; j < 31; ++j)
                G_dailyArray[i][j] = G_daily31GwRunoff[i][j];
        }
        sprintf(filename, "%s/G_GW_RUNOFF_%d_%d.31.UNF0", outputDir, year, month + 1);
        gridIO.writeUnfFile(filename, ng * 31, &G_dailyArray[0][0]);
    }
    if (options.outGwrunSurfrunDaily) {
        for (int i = 0; i < ng; ++i) {
            for (int j = 0; j < 31; ++j)
                G_dailyArray[i][j] = G_daily31GwrunSurfrun[i][j];
        }
        sprintf(filename, "%s/G_GWRUN_SURFRUN_%d_%d.31.UNF0", outputDir, year, month + 1);
        gridIO.writeUnfFile(filename, ng * 31, &G_dailyArray[0][0]);
    }

    if (options.outPrecDaily) { // new output options 2.2
        for (int i = 0; i < ng; ++i) {
            for (int j = 0; j < 31; ++j)
                G_dailyArray[i][j] = G_daily31ConsistentPrecip[i][j];
        }
        sprintf(filename, "%s/G_CONSISTENT PRECIPITATION_km3_%d_%d.31.UNF0", outputDir, year, month + 1);
        gridIO.writeUnfFile(filename, ng * 31, &G_dailyArray[0][0]);
    }
    if (options.outCellSurfaceDaily) { // new output options 2.2
        for (int i = 0; i < ng; ++i) {
            for (int j = 0; j < 31; ++j)
                G_dailyArray[i][j] = G_daily31CellSurfaceRunoff[i][j];
        }
        sprintf(filename, "%s/G_CELL_SURFACE_RUNOFF_%d_%d.31.UNF0", outputDir, year, month + 1);
        gridIO.writeUnfFile(filename, ng * 31, &G_dailyArray[0][0]);
    }
    if (options.outRiverAvailDaily) { // new output options 2.2
        for (int i = 0; i < ng; ++i) {
            for (int j = 0; j < 31; ++j)
                G_dailyArray[i][j] = G_daily31RiverAvail[i][j];
        }
        sprintf(filename, "%s/G_RIVER_AVAIL_%d_%d.31.UNF0", outputDir, year, month + 1);
        gridIO.writeUnfFile(filename, ng * 31, &G_dailyArray[0][0]);
    }
    if (options.outRiverVeloDaily) { // new output options 2.2
        for (int i = 0; i < ng; ++i) {
            for (int j = 0; j < 31; ++j)
                G_dailyArray[i][j] = G_daily31Velocity[i][j]; // conversion double -> float
        }
        sprintf(filename, "%s/G_RIVER_VELOCITY_%d_%d.31.UNF0", outputDir, year, month + 1);
        gridIO.writeUnfFile(filename, ng * 31, &G_dailyArray[0][0]);
    }
    if (options.outCellAETDaily) { // new output options 2.2
        for (int i = 0; i < ng; ++i) {
            for (int j = 0; j < 31; ++j)
                G_dailyArray[i][j] = G_daily31CellAET[i][j];
        }
        sprintf(filename, "%s/G_CELL_AET_%d_%d.31.UNF0", outputDir, year, month + 1);
        gridIO.writeUnfFile(filename, ng * 31, &G_dailyArray[0][0]);
    }
    if (options.outCellAETWCaDaily) { // new output options 2.2
        for (int i = 0; i < ng; ++i) {
            for (int j = 0; j < 31; ++j)
                G_dailyArray[i][j] = G_daily31CellAETWCa[i][j];
        }
        sprintf(filename, "%s/G_CELLAET_CONSUSE_km3_%d_%d.31.UNF0", outputDir, year, month + 1);
        gridIO.writeUnfFile(filename, ng * 31, &G_dailyArray[0][0]);
    }
    if (options.outSurfStorDaily) { // new output options 2.2
        for (int i = 0; i < ng; ++i) {
            for (int j = 0; j < 31; ++j)
                G_dailyArray[i][j] = G_daily31SurfStor[i][j];
        }
        sprintf(filename, "%s/G_SURFACE_WATER_STORAGE_km3_%d_%d.31.UNF0", outputDir, year, month + 1);
        gridIO.writeUnfFile(filename, ng * 31, &G_dailyArray[0][0]);
    }
    if (options.outSingleStoragesDaily) { // new output options 2.2
        for (int i = 0; i < ng; ++i) {
            for (int j = 0; j < 31; ++j)
                G_dailyArray[i][j] = G_daily31LocLakeStor[i][j];
        }
        sprintf(filename, "%s/G_LOC_LAKE_STORAGE_km3_%d_%d.31.UNF0", outputDir, year, month + 1);
        gridIO.writeUnfFile(filename, ng * 31, &G_dailyArray[0][0]);
    }
    if (options.outSingleStoragesDaily) { // new output options 2.2
        for (int i = 0; i < ng; ++i) {
            for (int j = 0; j < 31; ++j)
                G_dailyArray[i][j] = G_daily31LocWetlStor[i][j];
        }
        sprintf(filename, "%s/G_LOC_WETL_STORAGE_km3_%d_%d.31.UNF0", outputDir, year, month + 1);
        gridIO.writeUnfFile(filename, ng * 31, &G_dailyArray[0][0]);
    }
    if (options.outSingleStoragesDaily) { // new output options 2.2
        for (int i = 0; i < ng; ++i) {
            for (int j = 0; j < 31; ++j)
                G_dailyArray[i][j] = G_daily31GloLakeStor[i][j];
        }
        sprintf(filename, "%s/G_GLO_LAKE_STORAGE_km3_%d_%d.31.UNF0", outputDir, year, month + 1);
        gridIO.writeUnfFile(filename, ng * 31, &G_dailyArray[0][0]);
    }
    if (options.outSingleStoragesDaily) { // new output options 2.2
        for (int i = 0; i < ng; ++i) {
            for (int j = 0; j < 31; ++j)
                G_dailyArray[i][j] = G_daily31GloWetlStor[i][j];
        }
        sprintf(filename, "%s/G_GLO_WETL_STORAGE_km3_%d_%d.31.UNF0", outputDir, year, month + 1);
        gridIO.writeUnfFile(filename, ng * 31, &G_dailyArray[0][0]);
    }
    if (options.outSingleStoragesDaily) { // new output options 2.2
        for (int i = 0; i < ng; ++i) {
            for (int j = 0; j < 31; ++j)
                G_dailyArray[i][j] = G_daily31RiverStor[i][j];
        }
        sprintf(filename, "%s/G_RIVER_STORAGE_km3_%d_%d.31.UNF0", outputDir, year, month + 1);
        gridIO.writeUnfFile(filename, ng * 31, &G_dailyArray[0][0]);
    }
    if (options.outSingleStoragesDaily) { // new output options 2.2
        for (int i = 0; i < ng; ++i) {
            for (int j = 0; j < 31; ++j)
                G_dailyArray[i][j] = G_daily31ResStor[i][j];
        }
        sprintf(filename, "%s/G_RESERVOIR_STORAGE_km3_%d_%d.31.UNF0", outputDir, year, month + 1);
        gridIO.writeUnfFile(filename, ng * 31, &G_dailyArray[0][0]);
    }
    if (options.outSingleStoragesDaily || options.outGWStorageDaily) { // new output options 2.2
        for (int i = 0; i < ng; ++i) {
            for (int j = 0; j < 31; ++j)
                G_dailyArray[i][j] = G_daily31GwStor[i][j];
        }
        sprintf(filename, "%s/G_GROUND_WATER_STORAGE_km3_%d_%d.31.UNF0", outputDir, year, month + 1);
        gridIO.writeUnfFile(filename, ng * 31, &G_dailyArray[0][0]);
    }
    if (options.outTotalWaterInStoragesDaily_km3) { // new output options 2.2
        for (int i = 0; i < ng; ++i) {
            for (int j = 0; j < 31; ++j)
                G_dailyArray[i][j] = G_daily31TotalWaterInStorages_km3[i][j];
        }
        sprintf(filename, "%s/G_TOTAL_STORAGES_km3_%d_%d.31.UNF0", outputDir, year, month + 1);
        gridIO.writeUnfFile(filename, ng * 31, &G_dailyArray[0][0]);
    }
    if (options.outTotalWaterInStoragesDaily_mm) { // new output options 2.2
        for (int i = 0; i < ng; ++i) {
            for (int j = 0; j < 31; ++j)
                G_dailyArray[i][j] = G_daily31TotalWaterInStorages_mm[i][j];
        }
        sprintf(filename, "%s/G_TOTAL_STORAGES_mm_%d_%d.31.UNF0", outputDir, year, month + 1);
        gridIO.writeUnfFile(filename, ng * 31, &G_dailyArray[0][0]);
    }
    if (options.outGwrSwbDaily) { // groundwater recharge below surface waterbodies (with respect to continental area) [mm/month]
        for (int i = 0; i < ng; ++i) {
            for (int j = 0; j < 31; ++j)
                G_dailyArray[i][j] = G_daily31GwrSwb[i][j];
        }
        sprintf(filename, "%s/G_GWR_SURFACE_WATER_BODIES_mm_%d_%d.31.UNF0", outputDir, year);
        gridIO.writeUnfFile(filename, ng * 31, &G_dailyArray[0][0]);
    }
    if (options.outFswbDaily) { // fraction of surface water bodies
        for (int i = 0; i < ng; ++i) {
            for (int j = 0; j < 31; ++j)
                G_dailyArray[i][j] = G_daily31Fswb[i][j];
        }
        sprintf(filename, "%s/G_FRACTION_SURFACE_WATER_BODIES_%d_%d.31.UNF0", outputDir, year);
        gridIO.writeUnfFile(filename, ng * 31, &G_dailyArray[0][0]);
    }
    if (options.outLandAreaFracDaily) { // fraction of surface water bodies
        for (int i = 0; i < ng; ++i) {
            for (int j = 0; j < 31; ++j)
                G_dailyArray[i][j] = G_daily31LandAreaFrac[i][j];
        }
        sprintf(filename, "%s/G_LAND_AREA_FRACTION_%d_%d.31.UNF0", outputDir, year);
        gridIO.writeUnfFile(filename, ng * 31, &G_dailyArray[0][0]);
    }
}
// end writeDaily31Grids
// daily output option 31 end

// daily output option 365 start
void routingClass::writeDaily365Grids(char *outputDir, const int year) {
    if (5 != options.grid_store && 6 != options.grid_store)
        return;

    char filename[250];


    if (options.outCellRunoffDaily) { // new output options 2.2
        sprintf(filename, "%s/G_CELL_RUNOFF_%d.365.UNF0", outputDir, year);
        gridIO.writeUnfFile(filename, ng * 365, &G_daily365CellRunoff[0][0]);
    }
    if (options.outCellRunoffDaily) { // new output options 2.2
        sprintf(filename, "%s/G_CELL_RUNOFF_mm_%d.365.UNF0", outputDir, year);
        gridIO.writeUnfFile(filename, ng * 365, &G_daily365CellRunoff_mm[0][0]);
    }
    if (options.outLandAreaFracDaily) { // new output options 2.2
        sprintf(filename, "%s/G_LAND_AREA_FRACTION_%d.365.UNF0", outputDir, year);
        gridIO.writeUnfFile(filename, ng * 365, &G_daily365LandAreaFrac[0][0]);
    }
    if (options.outGWRunoffDaily) { // new output options 2.2
        sprintf(filename, "%s/G_GW_RUNOFF_%d.365.UNF0", outputDir, year);
        gridIO.writeUnfFile(filename, ng * 365, &G_daily365GwRunoff[0][0]);
    }
    if (options.outGwrunSurfrunDaily) { // new output options 2.2b ISIMIP2b
        sprintf(filename, "%s/G_GWRUN_SURFRUN_%d.365.UNF0", outputDir, year);
        gridIO.writeUnfFile(filename, ng * 365, &G_daily365GwrunSurfrun[0][0]);
    }
    if (options.outCellAETWCaDaily) { // new output options 2.2b ISIMIP2b
        sprintf(filename, "%s/G_CELLAET_CONSUSE_km3_%d.365.UNF0", outputDir, year);
        gridIO.writeUnfFile(filename, ng * 365, &G_daily365CellAETWCa[0][0]);
    }
    if (options.outPrecDaily) { // new output options 2.2
        sprintf(filename, "%s/G_CONSISTENT_PRECIPITATION_km3_%d.365.UNF0", outputDir, year);
        gridIO.writeUnfFile(filename, ng * 365, &G_daily365ConsistentPrecip[0][0]);
    }
    if (options.outCellSurfaceDaily) { // new output options 2.2
        sprintf(filename, "%s/G_CELL_SURFACE_RUNOFF_%d.365.UNF0", outputDir, year);
        gridIO.writeUnfFile(filename, ng * 365, &G_daily365CellSurfaceRunoff[0][0]);
    }
    if (options.outRiverAvailDaily) { // new output options 2.2
        sprintf(filename, "%s/G_RIVER_AVAIL_%d.365.UNF0", outputDir, year);
        gridIO.writeUnfFile(filename, ng * 365, &G_daily365RiverAvail[0][0]);
    }
    if (options.outRiverVeloDaily) { // new output options 2.2
        sprintf(filename, "%s/G_RIVER_VELOCITY_%d.365.UNF0", outputDir, year);
        gridIO.writeUnfFile(filename, ng * 365, &G_daily365Velocity[0][0]);
    }
    if (options.outCellAETDaily) { // new output options 2.2
        sprintf(filename, "%s/G_CELL_AET_%d.365.UNF0", outputDir, year);
        gridIO.writeUnfFile(filename, ng * 365, &G_daily365CellAET[0][0]);
    }
    if (options.outSurfStorDaily) { // new output options 2.2
        sprintf(filename, "%s/G_SURFACE_WATER_STORAGE_km3_%d.365.UNF0", outputDir, year);
        gridIO.writeUnfFile(filename, ng * 365, &G_daily365SurfStor[0][0]);
    }
    if (options.outSingleStoragesDaily) { // new output options 2.2
        sprintf(filename, "%s/G_LOC_LAKE_STORAGE_km3_%d.365.UNF0", outputDir, year);
        gridIO.writeUnfFile(filename, ng * 365, &G_daily365LocLakeStor[0][0]);
    }
    if (options.outSingleStoragesDaily) { // new output options 2.2
        sprintf(filename, "%s/G_LOC_WETL_STORAGE_km3_%d.365.UNF0", outputDir, year);
        gridIO.writeUnfFile(filename, ng * 365, &G_daily365LocWetlStor[0][0]);
    }
    if (options.outSingleStoragesDaily) { // new output options 2.2
        sprintf(filename, "%s/G_GLO_LAKE_STORAGE_km3_%d.365.UNF0", outputDir, year);
        gridIO.writeUnfFile(filename, ng * 365, &G_daily365GloLakeStor[0][0]);
    }
    if (options.outSingleStoragesDaily) { // new output options 2.2
        sprintf(filename, "%s/G_GLO_WETL_STORAGE_km3_%d.365.UNF0", outputDir, year);
        gridIO.writeUnfFile(filename, ng * 365, &G_daily365GloWetlStor[0][0]);
    }
    if (options.outSingleStoragesDaily) { // new output options 2.2
        sprintf(filename, "%s/G_RIVER_STORAGE_km3_%d.365.UNF0", outputDir, year);
        gridIO.writeUnfFile(filename, ng * 365, &G_daily365RiverStor[0][0]);
    }
    if (options.outSingleStoragesDaily) { // new output options 2.2
        sprintf(filename, "%s/G_RESERVOIR_STORAGE_km3_%d.365.UNF0", outputDir, year);
        gridIO.writeUnfFile(filename, ng * 365, &G_daily365ResStor[0][0]);
    }
    if (options.outSingleStoragesDaily || options.outGWStorageDaily) { // new output options 2.2
        sprintf(filename, "%s/G_GROUND_WATER_STORAGE_km3_%d.365.UNF0", outputDir, year);
        gridIO.writeUnfFile(filename, ng * 365, &G_daily365GwStor[0][0]);
    }
    if (options.outTotalWaterInStoragesDaily_km3) { // new output options 2.2
        sprintf(filename, "%s/G_TOTAL_STORAGES_km3_%d.365.UNF0", outputDir, year);
        gridIO.writeUnfFile(filename, ng * 365, &G_daily365TotalWaterInStorages_km3[0][0]);
    }
    if (options.outTotalWaterInStoragesDaily_mm) { // new output options 2.2
        sprintf(filename, "%s/G_TOTAL_STORAGES_mm_%d.365.UNF0", outputDir, year);
        gridIO.writeUnfFile(filename, ng * 365, &G_daily365TotalWaterInStorages_mm[0][0]);
    }
    if (options.outGwrSwbDaily) { // groundwater recharge below surface waterbodies (with respect to continental area) [mm/month]
        sprintf(filename, "%s/G_GWR_SURFACE_WATER_BODIES_mm_%d.365.UNF0", outputDir, year);
        gridIO.writeUnfFile(filename, ng * 365, &G_daily365GwrSwb[0][0]);
    }
    if (options.outFswbDaily) { // fraction of surface water bodies
        sprintf(filename, "%s/G_FRACTION_SURFACE_WATER_BODIES_%d.365.UNF0", outputDir, year);
        gridIO.writeUnfFile(filename, ng * 365, &G_daily365Fswb[0][0]);
    }
}
// end writeDaily365Grids
// daily output option 365 end

void routingClass::writeStartendGrids(char *outputDir, const int year) {
    char filename[250];

    sprintf(filename, "%s/G_TOTAL_STORAGES_STARTEND_km3_%d.2.UNF0", outputDir, year);
    gridIO.writeUnfFile(filename, ng * 2, &G_startendTotalWaterInStorages_km3[0][0]);
}

void routingClass::writeTotalCellRunoffGrid(char *outputDir, int year) {
    // used for calculation of correction factor during calibration
    char filename[250];
    int n;
    double G_array[ng];

    if (options.outPotCellRunoffAnnual) { // new output options 2.1f
        for (n = 0; n < ng; n++)
            G_array[n] = G_potCellRunoff[n];
        sprintf(filename, "%s/G_POT_CELL_RUNOFF_%d.UNF0", outputDir, year);
        gridIO.writeUnfFile(filename, ng, &G_array[0]);
    }
}
// end writeTotalCellRunoffGrid


// FP20161018N002 Reservoir operation start years
void routingClass::writeAnnualFreqGrids(char *outputDir, int year) {
    // yearly output for GFREQ, GFREQW and its changes (unit: percent of grid cell area)
    char filename[250];
    int n;
    float G_array[ng];

    if (options.outLandWaterFractions) { // new output options 2.2 ISIMIP
        for (n = 0; n < ng; n++)
            G_array[n] = G_fwaterfreq[n];
        sprintf(filename, "%s/GFREQW_%d.UNF0", outputDir, year);
        gridIO.writeUnfFile(filename, ng, &G_array[0]);

        for (n = 0; n < ng; n++)
            G_array[n] = G_fwaterfreq_change[n];
        sprintf(filename, "%s/GFREQW_change_%d.UNF0", outputDir, year);
        gridIO.writeUnfFile(filename, ng, &G_array[0]);

        for (n = 0; n < ng; n++)
            G_array[n] = G_landfreq[n];
        sprintf(filename, "%s/GFREQ_%d.UNF0", outputDir, year);
        gridIO.writeUnfFile(filename, ng, &G_array[0]);

        for (n = 0; n < ng; n++)
            G_array[n] = G_landfreq_change[n];
        sprintf(filename, "%s/GFREQ_change_%d.UNF0", outputDir, year);
        gridIO.writeUnfFile(filename, ng, &G_array[0]);

        for (n = 0; n < ng; n++)
            G_array[n] = G_landAreaFrac[n];
        sprintf(filename, "%s/GLAND_AREA_FRAC_%d.UNF0", outputDir, year);
        gridIO.writeUnfFile(filename, ng, &G_array[0]);

        for (n = 0; n < ng; n++)
            G_array[n] = G_landAreaFrac_change[n];
        sprintf(filename, "%s/GLAND_AREA_FRAC_change_%d.UNF0", outputDir, year);
        gridIO.writeUnfFile(filename, ng, &G_array[0]);
    }

}
// end writeAnnualFreqGrids

void createAsciiFile(char *filename, ofstream &streamName) {
    streamName.open(filename);
    if (!streamName) {
        cerr << "Can not open file " << filename << " for writing" << endl;
        exit(-1);
    }
    streamName << "# " << getTimeString();
}
// end createAsciiFile

void routingClass::createDailyFiles(const char *outputDir) {

    extern cbasinClass cbasin;

    char filename[250];
    if (options.outStationDischargeDaily) { // new output options 2.1f
        sprintf(filename, "%s/STATION_DISCHARGE_DAILY.OUT", outputDir);
        createAsciiFile(filename, dailyDischargeFile);
        dailyDischargeFile << "# Year\tDay\t";
        for (int n = 0; n <= cbasin.numberOfCalibBasins - 1; n++)
            dailyDischargeFile << cbasin.name[n] << '\t';
        dailyDischargeFile << endl;
    }
    // +++ daily velo file ++++++++++++++++++++++++++++++++++++++++++++++++ new flow velocity
    if (options.outStationVelocityDaily) { // new output options 2.2
        sprintf(filename, "%s/STATION_VELOCITY_DAILY.OUT", outputDir);
        createAsciiFile(filename, dailyVelocityFile);
        dailyVelocityFile << "# Year\tDay\t";
        for (int n = 0; n <= cbasin.numberOfCalibBasins - 1; n++)
            dailyVelocityFile << cbasin.name[n] << '\t';
        dailyVelocityFile << endl;
    }
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ new flow velocity: end

    if (options.outLocLakeStorageDaily) { // new output options 2.1f
        sprintf(filename, "%s/LOC_LAKE_STORAGE.OUT", outputDir);
        createAsciiFile(filename, locLakeFile);
        for (int n = 0; n <= cbasin.numberOfCalibBasins - 1; n++)
            locLakeFile << cbasin.name[n] << '\t';
        locLakeFile << endl;
    }
    if (options.outLocWetlStorageDaily) { // new output options 2.1f
        sprintf(filename, "%s/LOC_WETL_STORAGE.OUT", outputDir);
        createAsciiFile(filename, locWetlFile);
        for (int n = 0; n <= cbasin.numberOfCalibBasins - 1; n++)
            locWetlFile << cbasin.name[n] << '\t';
        locWetlFile << endl;
    }
    if (options.outGloLakeStorageDaily) { // new output options 2.1f
        sprintf(filename, "%s/GLO_LAKE_STORAGE.OUT", outputDir);
        createAsciiFile(filename, gloLakeFile);
        for (int n = 0; n <= cbasin.numberOfCalibBasins - 1; n++)
            gloLakeFile << cbasin.name[n] << '\t';
        gloLakeFile << endl;
    }

    // new output options.outResStorageDaily
    if (options.outResStorageDaily && (options.resOpt == 1)) { // new output options 2.1f
        sprintf(filename, "%s/RESERVOIR_STORAGE.OUT", outputDir);
        createAsciiFile(filename, gloResFile);
        for (int n = 0; n <= cbasin.numberOfCalibBasins - 1; n++)
            gloResFile << cbasin.name[n] << '\t';
        gloResFile << endl;

    }
    if (options.outGloWetlStorage) { // new output options 2.1f
        sprintf(filename, "%s/GLO_WETL_STORAGE.OUT", outputDir);
        createAsciiFile(filename, gloWetlFile);
        for (int n = 0; n <= cbasin.numberOfCalibBasins - 1; n++)
            gloWetlFile << cbasin.name[n] << '\t';
        gloWetlFile << endl;
    }
}
// end createDailyFiles

void routingClass::closeDailyFiles() {
    if (options.outStationDischargeDaily) dailyDischargeFile.close();
    if (options.outStationVelocityDaily) dailyVelocityFile.close(); //+++ new flow velocity
    if (options.outLocLakeStorageDaily) locLakeFile.close();
    if (options.outLocWetlStorageDaily) locWetlFile.close();
    if (options.outGloLakeStorageDaily) gloLakeFile.close();
    if (options.outResStorageDaily && (options.resOpt == 1))
        gloResFile.close(); // new output options.outResStorageDaily
    if (options.outGloWetlStorageDaily) gloWetlFile.close();
}
// end closeDailyFiles

//-----------------------Velo-File----------------------------------------------------
void routingClass::appendDailyVelocity(const short actualYear, const short remainingInitYears) {
    // write daily river velocity information into file
    if (options.outStationVelocityDaily) {
        for (short d = 0; d <= 364; d++) {
            dailyVelocityFile << setw(4) << actualYear << '\t' << setw(3) << d + 1 << '\t';
            for (short b = 0; b < nSpecBasins; ++b)
                dailyVelocityFile << setiosflags(ios::fixed | ios::showpoint)
                                  << setprecision(6) << setw(13)
                                  << dailyRiverVelocity[b][d] << '\t';
            dailyVelocityFile << endl;
        }
    }
}
// end appendDailyVelocity
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ new flow velocity


void routingClass::appendDailyDischarge(const short actualYear, const short remainingInitYears) {
    // write daily river discharge information into file
    if (options.outStationDischargeDaily) { // new output options 2.1f
        for (short d = 0; d <= 364; d++) {
            dailyDischargeFile << setw(4) << actualYear << '\t' << setw(3) << d + 1 << '\t';
            for (short b = 0; b < nSpecBasins; ++b)
                dailyDischargeFile << setiosflags(ios::fixed | ios::showpoint)
                                   << setprecision(6) << setw(13)
                                   << dailyRiverDischarge[b][d] << '\t';
            dailyDischargeFile << endl;
        }
    }
}
// end appendDailyDischarge

void routingClass::appendDailyLakeWetl(const short actualYear, const short remainingInitYears) {
    if (options.outLocLakeStorageDaily) { // new output options 2.1f
        // store daily values of lake and wetland water content
        for (short d = 0; d <= 364; d++) {
            locLakeFile << setw(4) << actualYear << '\t' << setw(3) << d + 1 << '\t';
            for (short b = 0; b < nSpecBasins; ++b)
                locLakeFile << setiosflags(ios::fixed | ios::showpoint)
                            << setprecision(6) << setw(13)
                            << dailyLocLakeStorage[b][d] << '\t';
            locLakeFile << endl;
        }
    }
    if (options.outLocWetlStorageDaily) { // new output options 2.1f
        for (short d = 0; d <= 364; d++) {
            locWetlFile << setw(4) << actualYear << '\t' << setw(3) << d + 1 << '\t';
            for (short b = 0; b < nSpecBasins; ++b)
                locWetlFile << setiosflags(ios::fixed | ios::showpoint)
                            << setprecision(6) << setw(13)
                            << dailyLocWetlStorage[b][d] << '\t';
            locWetlFile << endl;
        }
    }
    if (options.outGloLakeStorageDaily) { // new output options 2.1f
        for (short d = 0; d <= 364; d++) {
            gloLakeFile << setw(4) << actualYear << '\t' << setw(3) << d + 1 << '\t';
            for (short b = 0; b < nSpecBasins; ++b)
                gloLakeFile << setiosflags(ios::fixed | ios::showpoint)
                            << setprecision(6) << setw(13)
                            << dailyGloLakeStorage[b][d] << '\t';
            gloLakeFile << endl;
        }
    }

    // new output options.outResStorageDaily
    if (options.outResStorageDaily && (options.resOpt == 1)) { // new output options 2.1f
        for (short d = 0; d <= 364; d++) {
            gloResFile << setw(4) << actualYear << '\t' << setw(3) << d + 1 << '\t';
            for (short b = 0; b < nSpecBasins; ++b)
                gloResFile << setiosflags(ios::fixed | ios::showpoint)
                           << setprecision(6) << setw(13) << dailyResStorage[b][d]
                           << '\t';
            gloResFile << endl;
        }
    }

    if (options.outGloWetlStorageDaily) { // new output options 2.1f
        for (short d = 0; d <= 364; d++) {
            gloWetlFile << setw(4) << actualYear << '\t' << setw(3) << d + 1 << '\t';
            for (short b = 0; b < nSpecBasins; ++b)
                gloWetlFile << setiosflags(ios::fixed | ios::showpoint)
                            << setprecision(6) << setw(13)
                            << dailyGloWetlStorage[b][d] << '\t';
            gloWetlFile << endl;
        }
    }
}
// end appendDailyLakeWetl

void routingClass::createMonthlyDischargeFile(const char *outputDir) {
    char filename[250];
    if (options.outStationDischargeMonthly) { // new output options 2.1f

        sprintf(filename, "%s/STATION_DISCHARGE_MONTHLY.OUT", outputDir);
        createAsciiFile(filename, monthlyDischargeFile);
        monthlyDischargeFile
                       << "Stationname\t"; //HMS 2017-03-16: changed from "# Station name" to one single word in order to easy handling of output file e.g. in R-scripts.
        for (int n = 0; n <= cbasin.numberOfCalibBasins - 1; n++)
            monthlyDischargeFile << cbasin.name[n] << '\t';
        monthlyDischargeFile << endl;
    }

    // new output option (outResvoirMonthly)
    if (options.outResvoirMonthly && // output option for this 4 files
        options.resOpt == 1 && // New reservoir algorithm (N. Hanasaki) will be used (*)
        ((2 == options.grid_store) || (4 == options.grid_store) ||
         (6 == options.grid_store)) // Store grid files, including monthly values
                   ) {
        // create reservoir storage file (ascii)
        sprintf(filename, "%s/RESERVOIR_STORAGE_MONTHLY.OUT", outputDir);
        createAsciiFile(filename, monthlyReservoirStorageFile);
        monthlyReservoirStorageFile << "# Year" << '\t' << "Month" << '\t';
        for (int n = 0; n < ng; n++) {
            if (G_reservoir_area[n] > 0.0) {
                monthlyReservoirStorageFile << n + 1 << '\t';
            }
        }
        monthlyReservoirStorageFile << endl;

        // create reservoir storage ratio file (ascii)
        sprintf(filename, "%s/RESERVOIR_STORAGE_RATIO_MONTHLY.OUT", outputDir);
        createAsciiFile(filename, monthlyReservoirStorageRatioFile);
        monthlyReservoirStorageRatioFile << "# Year" << '\t' << "Month" << '\t';
        for (int n = 0; n < ng; n++) {
            if (G_reservoir_area[n] > 0.0) {
                monthlyReservoirStorageRatioFile << n + 1 << '\t';
            }
        }
        monthlyReservoirStorageRatioFile << endl;

        // create reservoir inflow file (ascii)
        sprintf(filename, "%s/RESERVOIR_INFLOW_MONTHLY.OUT", outputDir);
        createAsciiFile(filename, monthlyReservoirInflowFile);
        monthlyReservoirInflowFile << "# Year" << '\t' << "Month" << '\t';
        for (int n = 0; n < ng; n++) {
            if (G_reservoir_area[n] > 0.0) {
                monthlyReservoirInflowFile << n + 1 << '\t';
            }
        }
        monthlyReservoirInflowFile << endl;

        // create reservoir outflow file (ascii)
        sprintf(filename, "%s/RESERVOIR_OUTFLOW_MONTHLY.OUT", outputDir);
        createAsciiFile(filename, monthlyReservoirOutflowFile);
        monthlyReservoirOutflowFile << "# Year" << '\t' << "Month" << '\t';
        for (int n = 0; n < ng; n++) {
            if (G_reservoir_area[n] > 0.0) {
                monthlyReservoirOutflowFile << n + 1 << '\t';
            }
        }
        monthlyReservoirOutflowFile << endl;
    }

}
// end createMonthlyDischargeFile

void routingClass::closeMonthlyDischargeFile() {
    if (options.outStationDischargeMonthly)
        monthlyDischargeFile.close();

    if (options.outResvoirMonthly && // output option for this 4 files
        options.resOpt == 1 && // New reservoir algorithm (N. Hanasaki) will be used (*)
        ((2 == options.grid_store) || (4 == options.grid_store) ||
         (6 == options.grid_store)) // Store grid files, including monthly values
                   ) {
        monthlyReservoirStorageFile.close();
        monthlyReservoirStorageRatioFile.close();
        monthlyReservoirInflowFile.close();
        monthlyReservoirOutflowFile.close();
    }
}
// end closeMonthlyDischargeFile

void routingClass::appendMonthlyDischarge(const short actualYear, const short remainingInitYears) {
    char number_of_days_in_month[12] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};

    if (options.outStationDischargeMonthly) { // new output options 2.1f
        short day = 0;

        for (short i = 0; i <= 11; i++) {
            if (remainingInitYears != 0)
                monthlyDischargeFile << "# ";
            monthlyDischargeFile << setw(4) << (actualYear - 1901) * 12 + i + 1 << '\t';
            for (short b = 0; b < nSpecBasins; ++b) {
                double monthlyDischarge = 0;
                for (short d = 0; d <= number_of_days_in_month[i] - 1; d++) {
                    monthlyDischarge += dailyRiverDischarge[b][day + d];
                }
                monthlyDischargeFile
                               << setiosflags(ios::fixed | ios::showpoint)
                               << setprecision(6) << setw(13) << monthlyDischarge << '\t';
            }
            monthlyDischargeFile << endl;
            day += number_of_days_in_month[i];
        }
    }

    if (options.outResvoirMonthly && // output option for this 4 files
        options.resOpt == 1 && // New reservoir algorithm (N. Hanasaki) will be used
        ((2 == options.grid_store) || (4 == options.grid_store) ||
         (6 == options.grid_store)) // Store grid files, including monthly values
                   ) {
        for (short m = 0; m < 12; m++) {
            if (remainingInitYears != 0) {
                monthlyReservoirStorageFile << "# ";
                monthlyReservoirStorageRatioFile << "# ";
                monthlyReservoirInflowFile << "# ";
                monthlyReservoirOutflowFile << "# ";
            }
            monthlyReservoirStorageFile << setw(4) << actualYear << '\t';
            monthlyReservoirStorageFile << setw(4) << m + 1 << '\t';

            monthlyReservoirStorageRatioFile << setw(4) << actualYear << '\t';
            monthlyReservoirStorageRatioFile << setw(4) << m + 1 << '\t';

            monthlyReservoirInflowFile << setw(4) << actualYear << '\t';
            monthlyReservoirInflowFile << setw(4) << m + 1 << '\t';

            monthlyReservoirOutflowFile << setw(4) << actualYear << '\t';
            monthlyReservoirOutflowFile << setw(4) << m + 1 << '\t';

            for (int n = 0; n < ng; n++) {
                //G_monthlyResStorage[n][m] /= number_of_days_in_month[m];

                if (G_reservoir_area[n] > 0.0) {
                    if (options.grid_store_TypeForStorages == 0)
                        G_monthlyResStorageRatio[n][m] = G_monthlyResStorage[n][m] / G_stor_cap[n] * 0.85;
                    else
                        G_monthlyResStorageRatio[n][m] =
                                       G_monthlyResStorage[n][m] / number_of_days_in_month[m] / G_stor_cap[n] * 0.85;

                    // (in routing, timestep flow values [m3/timestep] are summed up to calculate total monthly values here [m3/month]!)
                    // [m3/s]
                    G_monthlyResInflow[n][m] = G_monthlyResInflow[n][m] / (number_of_days_in_month[m] * 86400.);
                    G_monthlyResOutflow[n][m] = G_monthlyResOutflow[n][m] / (number_of_days_in_month[m] * 86400.);

                    double tmp_G_monthlyResStorage;
                    if (options.grid_store_TypeForStorages == 0) tmp_G_monthlyResStorage = G_monthlyResStorage[n][m];
                    else tmp_G_monthlyResStorage = G_monthlyResStorage[n][m] / number_of_days_in_month[m];

                    monthlyReservoirStorageFile << setiosflags(ios::fixed | ios::showpoint)
                                                << setprecision(6) << setw(13) << tmp_G_monthlyResStorage << '\t';
                    monthlyReservoirStorageRatioFile << setiosflags(ios::fixed | ios::showpoint)
                                                     << setprecision(6) << setw(13) << G_monthlyResStorageRatio[n][m]
                                                     << '\t';
                    monthlyReservoirInflowFile << setiosflags(ios::fixed | ios::showpoint)
                                               << setprecision(6) << setw(13) << (G_monthlyResInflow[n][m]) << '\t';
                    monthlyReservoirOutflowFile << setiosflags(ios::fixed | ios::showpoint)
                                                << setprecision(6) << setw(13) << (G_monthlyResOutflow[n][m]) << '\t';
                } else
                    G_monthlyResStorageRatio[n][m] = -99.99;

            }
            monthlyReservoirStorageFile << endl;
            monthlyReservoirStorageRatioFile << endl;
            monthlyReservoirInflowFile << endl;
            monthlyReservoirOutflowFile << endl;
        }
    }
}
// end appendMonthlyDischarge

void routingClass::createAnnualDischargeFile(const char *outputDir) {
    if (options.outStationDischargeAnnual) { // new output options 2.1f
        char filename[250];

        sprintf(filename, "%s/STATION_DISCHARGE_ANNUAL.OUT", outputDir);
        createAsciiFile(filename, annualDischargeFile);
        annualDischargeFile
                       << "Stationname\t"; //HMS 2017-03-16: changed from "# Station name" to one single word in order to easy handling of output file e.g. in R-scripts.
        for (int n = 0; n <= cbasin.numberOfCalibBasins - 1; n++)
            annualDischargeFile << cbasin.name[n] << '\t';
        annualDischargeFile << endl;
    }
}
// end createAnnualDischargeFile

void routingClass::closeAnnualDischargeFile() {
    if (options.outStationDischargeAnnual) annualDischargeFile.close();
}
// end closeAnnualDischargeFile

void routingClass::appendAnnualDischarge(const short actualYear, const short remainingInitYears) {
    if (options.outStationDischargeAnnual) { // new output options 2.1f
        if (remainingInitYears != 0)
            annualDischargeFile << "# ";
        annualDischargeFile << setw(4) << actualYear << '\t';

        for (short b = 0; b < nSpecBasins; ++b) {
            double annualDischarge = 0;
            for (short i = 0; i <= 364; i++) {
                annualDischarge += dailyRiverDischarge[b][i];
            }
            annualDischargeFile << setiosflags(ios::fixed | ios::showpoint)
                                << setprecision(6) << setw(13) << annualDischarge << '\t';
        }
        annualDischargeFile << endl;
    }
}
// end appendAnnualDischarge

void routingClass::appendCalibInfoToAnnualDischargeFile(const double gamma) {
    if (options.outStationDischargeAnnual) { // new output options 2.1f
        annualDischargeFile << "# Next step of calibration." << endl;
        annualDischargeFile << "# Gamma = " << gamma << endl;
    }
}
// end appendCalibInfoToAnnualDischargeFile

double routingClass::getAnnualDischarge(const int stationNumber) {
    double annualDischarge = 0;

    for (short i = 0; i <= 364; i++) {
        annualDischarge += dailyRiverDischarge[stationNumber - 1][i];
    }
    return annualDischarge;
}
// end getAnnualDischarge

double routingClass::getAnnualUpstreamInflow(const int stationNumber) {
    // calculates annual inflow from upstream basins

    extern upstreamStationClass directUpstSt;
    int upSt;
    short i;
    double stationAnnualDischarge = 0, totalAnnualDischarge = 0;

    for (int n = 1; n <= directUpstSt.getNumberOfUpstreamStations(stationNumber); n++) {
        stationAnnualDischarge = 0;
        upSt = directUpstSt.getUpstreamStation(stationNumber, n);
        for (i = 0; i <= 364; i++) {
            stationAnnualDischarge += dailyRiverDischarge[upSt - 1][i];
        }
        totalAnnualDischarge += stationAnnualDischarge;
    }
    return totalAnnualDischarge;
}
// end getAnnualUpstreamInflow

double routingClass::getAnnualSatisfWaterUse(const int stationNumber) {
    //extern waterUseClass waterUse;
    extern signed short G_sbasin[ng];

    double satisfiedBasinUse = 0;

    for (int n = 0; n <= ng - 1; n++) {
        if (stationNumber == G_sbasin[n]) {
            //satisfiedBasinUse += G_satisfiedUse[n];
            // CR 2015-09: A distinction between satisfied use and actual use is not possible in 22b.
            satisfiedBasinUse += G_actualUse[n];
        }
    }

    return satisfiedBasinUse;
}
// end getAnnualSatisfWaterUse

double routingClass::getRiverVelocity(const double RiverSlope, const double G_riverBottomWidth,
                                      const double Roughness, double G_riverInflow, const int n) {
    double river_velocity;

    double riverDepth;
    double wettedPerimeter;
    double crossSectionalArea;
    double hydraulicRad;
    double incoming_discharge;

    incoming_discharge = (G_riverInflow * 1000. * 1000. * 1000.) / (60. * 60. * 24.); //[km3/day]  --> [m3/sec]

    riverDepth = 0.349 * pow(incoming_discharge, 0.341); //[m]

    // trapezoidal channel shape with channel sides 2:1 run to rise ratio
    crossSectionalArea = riverDepth * (2.0 * riverDepth + G_riverBottomWidth); // trapezoidal river shape
    wettedPerimeter = G_riverBottomWidth + 2.0 * riverDepth * sqrt(5.0); // sqrt(1+2^2)
    //crossSectionalArea = riverDepth * (riverDepth + G_riverBottomWidth); // trapezoidal river shape
    //wettedPerimeter = G_riverBottomWidth + 2.0 * riverDepth * sqrt(2.0); // sqrt(1+2^2)

    hydraulicRad = crossSectionalArea / wettedPerimeter;

    // calculate riverVelocity
    // Cell-specific calibration parameters - Apply multiplier  // FP
    river_velocity = 1. / (calibParam.getValue(M_RIVRGH_C, n) * Roughness) * pow(hydraulicRad, (2. / 3.)) *
                     pow(RiverSlope, 0.5); //[m/sec]
    river_velocity = river_velocity * 86.4; //m/sec -->km/day (60*60*24)/1000

    // return value in [km/day]; if value is below 1cm/day then set limit

    if (river_velocity < 0.00001)
        return 0.00001;
    else
        return river_velocity;
    //return G_riverInflow;
    //return RiverSlope;
    //return G_riverBottomWidth;
    //return Roughness;

}
// end getRiverVelocity

double routingClass::getNewRiverVelocity(const double RiverSlope, const double G_riverBottomWidth,
                                         const double Roughness, double G_riverStorage, double G_riverLength, const int n) {
    double river_velocity;

    double riverDepth;
    double wettedPerimeter;
    double crossSectionalArea;
    double hydraulicRad;
    double incoming_discharge;

    // trapezoidal channel shape with channel sides 2:1 run to rise ratio
    crossSectionalArea = (G_riverStorage / G_riverLength) * 1000000.;	// km3/km *10^6 -> m2

    // G_riverBottomWidth in m
    riverDepth = - G_riverBottomWidth/4.0 + sqrt(G_riverBottomWidth * G_riverBottomWidth/16.0 + 0.5 * crossSectionalArea);

    // Version 22 and 22b:
    wettedPerimeter = G_riverBottomWidth + 2.0 * riverDepth * sqrt(5.0); // sqrt(1+2^2)

    hydraulicRad = crossSectionalArea / wettedPerimeter;

    // calculate riverVelocity
    // Cell-specific calibration parameters - Apply multiplier  // FP
    river_velocity = 1. / (calibParam.getValue(M_RIVRGH_C,n) * Roughness) * pow(hydraulicRad, (2. / 3.)) * pow(RiverSlope, 0.5); //[m/sec]
    river_velocity = river_velocity * 86.4; //m/sec -->km/day (60*60*24)/1000

    // return value in [km/day]; if value is below 1cm/day then set limit

    if (river_velocity < 0.00001)
        return 0.00001;
    else
        return river_velocity;
    //return G_riverInflow;
    //return RiverSlope;
    //return G_riverBottomWidth;
    //return Roughness;
}
// end getNewRiverVelocity

void routingClass::getAnnualDischargeVector(vector<float>&annualDischargeVector) {
    for (short n = 0; n < nSpecBasins; n++) {
        float annualDischarge = 0;
        for (short i = 0; i <= 364; i++) {
            annualDischarge += dailyRiverDischarge[n][i];
        }
        annualDischarge *= 1000000000 / (double) (365 * 24 * 60 * 60); // km3/a -> m3/s
        annualDischargeVector.push_back(annualDischarge);
    }
}
// end getAnnualDischargeVector

void routingClass::getMonthlyDischargeVector(short month, vector<float> &monthlyDischargeVector) {
    short day = 0;
    char number_of_days_in_month[12] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
    for (short i = 0; i <= 11; i++) {
        if (i == month) {
            for (short b = 0; b < nSpecBasins; ++b) {
                float monthlyDischarge = 0;
                for (short d = 0; d < number_of_days_in_month[i]; d++) {
                    monthlyDischarge += dailyRiverDischarge[b][day + d];
                }
                monthlyDischarge *= 1000000000 / (float) (24 * 60 * 60 * number_of_days_in_month[i]);
                monthlyDischargeVector.push_back(monthlyDischarge);
            }
            return;
        }
        day += number_of_days_in_month[i];
    }
}
// end getMonthlyDischargeVector

void routingClass::getMinMonthlyDischargeVector(vector<float> &monthlyDischargeVector) {
    char number_of_days_in_month[12] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
    for (short b = 0; b < nSpecBasins; ++b) {
        float minValue = 99E29;
        short day = 0;
        for (short i = 0; i <= 11; i++) {
            float monthlyDischarge = 0;
            for (short d = 0; d <= number_of_days_in_month[i] - 1; d++) {
                monthlyDischarge += dailyRiverDischarge[b][day + d];
            }
            monthlyDischarge *= 1000000000 / (float) (24 * 60 * 60 * number_of_days_in_month[i]);
            if (monthlyDischarge < minValue)
                minValue = monthlyDischarge;
            day += number_of_days_in_month[i];
        }
        monthlyDischargeVector.push_back(minValue);
    }
}
// end getMinMonthlyDischargeVector

void routingClass::getMinToMaxDischargeVector(vector<float> &monthlyDischargeVector) {
    char number_of_days_in_month[12] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
    for (short b = 0; b < nSpecBasins; ++b) {
        float minValue = 99E29;
        float maxValue = -99E29;
        short day = 0;
        for (short i = 0; i <= 11; i++) {
            float monthlyDischarge = 0;
            for (short d = 0; d <= number_of_days_in_month[i] - 1; d++) {
                monthlyDischarge += dailyRiverDischarge[b][day + d];
            }
            monthlyDischarge *= 1000000000 / (float) (24 * 60 * 60 * number_of_days_in_month[i]);
            if (monthlyDischarge < minValue)
                minValue = monthlyDischarge;
            if (monthlyDischarge > maxValue)
                maxValue = monthlyDischarge;
            day += number_of_days_in_month[i];
        }
        monthlyDischargeVector.push_back(minValue / maxValue);
    }
}
// end getMinToMaxDischargeVector

// Implementation of class method calcNextDay_M(...)
void routingClass::calcNextDay_M(const short new_day) {
    static short month = 11;
    static short old_day;
    short old_month;

    const short int last_day_per_month[12]
                   = {31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365};
    old_month = month;

    // does the day beint to the next year already ?
    if (new_day == 1) {
        month = 0;
        old_day = 0;
    }
    // does the day beint to the next month ?
    while (new_day > last_day_per_month[month])
        month++;

    //cout << new_day << ' ' << month << endl;

    if (new_day != old_day + 1) {
        cerr << "Day counter mismatch!\n";
        exit(-1);
    }

    for (int n = 0; n <= ng - 1; n++) {
        //G_totalDailyUse[n] = G_dailyWaterUse[n][month];
        G_dailydailyNUs[n] = G_dailyNUs[n][month];
        G_dailydailyNUg[n] = G_dailyNUg[n][month];

        if (1 == options.aggrNUsGloLakResOpt) {
            G_dailydailyNUsAggregated[n] = G_dailyNUsAggregated[n][month];
        }}
    old_day = new_day;
}
// end calcNextDay_M

routingClass::~routingClass() {
    delete[] dailyRiverDischarge;
    delete[] G_monthlySwStorage;
    delete[] G_monthlyGwStorage;
    delete[] G_monthlyGwRunoff;
    delete[] G_monthlyLocWetlExtent;
    delete[] G_monthlyGloWetlExtent;

    if (options.day_store == 1) {
        delete[] dailyRiverVelocity; // ++++++++++++++++++++++++++++++++++++++++++++++++ new flow velocity
        delete[] dailyLocLakeStorage;
        delete[] dailyLocWetlStorage;
        delete[] dailyGloLakeStorage;
        if (options.resOpt == 1) delete[] dailyResStorage;
        delete[] dailyGloWetlStorage;
    }

    if (options.resOpt == 1) {
        delete[] G_gloResStorage;
        G_gloResStorage = NULL;
        delete[] G_res_type;
        G_res_type = NULL;
        delete[] G_start_month;
        G_start_month = NULL;
        delete[] G_mean_outflow;
        G_mean_outflow = NULL;
        // FP20161018N002 Reservoir operation start years
        // G_res_start_year & G_stor_cap_full (instead G_stor_cap)
        delete[] G_res_start_year;
        G_res_start_year = NULL;
        delete[] G_stor_cap_full;
        G_stor_cap_full = NULL;
        delete[] G_mean_demand;
        G_mean_demand = NULL;
        delete[] G_mean_cons_use;
        G_mean_cons_use = NULL;
        delete[] G_alloc_coeff;
        G_alloc_coeff = NULL;
        delete[] K_release;
        K_release = NULL;
        delete[] annual_release;
        annual_release = NULL;
    }

    if ((2 == options.grid_store) || (4 == options.grid_store) || (6 == options.grid_store)) {
        delete[] G_monthlyRiverAvail;
        G_monthlyRiverAvail = NULL;
        delete[] G_monthlyConsistentPrecip;
        G_monthlyConsistentPrecip = NULL;
        delete[] G_monthlyPotCellRunoff;
        G_monthlyPotCellRunoff = NULL;
        delete[] G_monthlyPotCellRunoffDeficit;
        G_monthlyPotCellRunoffDeficit = NULL;

        if (options.outRiverInUpstream)
            delete[] G_monthlyRiverInUpstream;
        G_monthlyRiverInUpstream = NULL;
        delete[] G_monthlyCellRunoff;
        G_monthlyCellRunoff = NULL;
        delete[] G_monthlyCellSurfaceRunoff;
        G_monthlyCellSurfaceRunoff = NULL;
        // added for cell AET calculation (2.1f)
        delete[]G_monthlyCellAET;
        G_monthlyCellAET = NULL;
        delete[]G_monthlyCellAETWCa;
        G_monthlyCellAETWCa = NULL;
        delete[]G_monthlyOpenWaterEvap;
        G_monthlyOpenWaterEvap = NULL;

        delete[] G_monthlyVelocity;
        G_monthlyVelocity = NULL;
        delete[] G_monthlySurfStor;
        G_monthlySurfStor = NULL;
        delete[] G_monthlySurfStor_mm;
        G_monthlySurfStor_mm = NULL;

        delete[] G_monthlyMinRiverAvail;
        G_monthlyMinRiverAvail = NULL;
        delete[] G_monthlyMaxRiverAvail;
        G_monthlyMaxRiverAvail = NULL;

        delete[] G_monthlyGwrSwb;
        G_monthlyGwrSwb = NULL;
        delete[] G_monthlyFswb;
        G_monthlyFswb = NULL;
        if (options.outRiverPET) delete[] G_monthlyRiverAreaFrac;
        G_monthlyRiverAreaFrac = NULL;
        if (options.outRiverPET) delete[] G_monthlyRiverPET;
        G_monthlyRiverPET = NULL;

        //define arrays for new output
        if (options.resOpt == 1) {
            delete[] G_monthlyResStorage;
            G_monthlyResStorage = NULL;
            delete[] G_monthlyResStorageRatio;
            G_monthlyResStorageRatio = NULL;
            delete[] G_monthlyResInflow;
            G_monthlyResInflow = NULL;
            delete[] G_monthlyResOutflow;
            G_monthlyResOutflow = NULL;
        }
        delete[] G_monthlyRiverStorage;
        G_monthlyRiverStorage = NULL;
        delete[] G_monthlyLocLakeStorage;
        G_monthlyLocLakeStorage = NULL;
        delete[] G_monthlyLocWetlStorage;
        G_monthlyLocWetlStorage = NULL;
        delete[] G_monthlyGloLakeStorage;
        G_monthlyGloLakeStorage = NULL;
        delete[] G_monthlyGloWetlStorage;
        G_monthlyGloWetlStorage = NULL;
        delete[] G_monthlySatisfiedUse;
        G_monthlySatisfiedUse = NULL;
        delete[] G_monthlyActualUse;
        G_monthlyActualUse = NULL;
        delete[] G_monthlyNUg;
        G_monthlyNUg = NULL;
        // FP 2015, not used anymore, at the moment, implement if needed again
        //CR 15-08-19 G_AllocatedUse = test output: demand allocated from neighboring cell(s)
        delete[] G_monthlyAllocatedUse;
        G_monthlyAllocatedUse = NULL;
    }

    // daily output option 31 start
    if ((3 == options.grid_store) || (4 == options.grid_store)) {
        if (options.outPrecDaily) {
            delete[] G_daily31ConsistentPrecip;
            G_daily31ConsistentPrecip = NULL;
        }
        if (options.outCellRunoffDaily) {
            delete[] G_daily31CellRunoff;
            G_daily31CellRunoff = NULL;
        }
        if (options.outGWRunoffDaily) {
            delete[] G_daily31GwRunoff;
            G_daily31GwRunoff = NULL;
        }
        if (options.outCellSurfaceDaily) {
            delete[] G_daily31CellSurfaceRunoff;
            G_daily31CellSurfaceRunoff = NULL;
        }
        if (options.outRiverAvailDaily) {
            delete[] G_daily31RiverAvail;
            G_daily31RiverAvail = NULL;
        }
        if (options.outRiverVeloDaily) {
            delete[] G_daily31Velocity;
            G_daily31Velocity = NULL;
        }
        if (options.outCellAETDaily) {
            delete[] G_daily31CellAET;
            G_daily31CellAET = NULL;
        }
        if (options.outCellAETWCaDaily) {
            delete[] G_daily31CellAETWCa;
            G_daily31CellAETWCa = NULL;
        }
        if (options.outSurfStorDaily) {
            delete[] G_daily31SurfStor;
            G_daily31SurfStor = NULL;
        }
        if (options.outSingleStoragesDaily) {
            delete[] G_daily31LocLakeStor;
            G_daily31LocLakeStor = NULL;
        }
        if (options.outSingleStoragesDaily) {
            delete[] G_daily31LocWetlStor;
            G_daily31LocWetlStor = NULL;
        }
        if (options.outSingleStoragesDaily) {
            delete[] G_daily31GloLakeStor;
            G_daily31GloLakeStor = NULL;
        }
        if (options.outSingleStoragesDaily) {
            delete[] G_daily31GloWetlStor;
            G_daily31GloWetlStor = NULL;
        }
        if (options.outSingleStoragesDaily) {
            delete[] G_daily31RiverStor;
            G_daily31RiverStor = NULL;
        }
        if (options.outSingleStoragesDaily) {
            delete[] G_daily31ResStor;
            G_daily31ResStor = NULL;
        }
        if (options.outSingleStoragesDaily || options.outGWStorageDaily) {
            delete[] G_daily31GwStor;
            G_daily31GwStor = NULL;
        }
        if (options.outTotalWaterInStoragesDaily_km3 || options.outTotalWaterInStoragesDaily_mm)
            delete[] G_daily31TotalWaterInStorages_km3;
        G_daily31TotalWaterInStorages_km3 = NULL;
        if (options.outTotalWaterInStoragesDaily_mm) delete[] G_daily31TotalWaterInStorages_mm;
        G_daily31TotalWaterInStorages_mm = NULL;
        if (options.outGwrSwbDaily) {
            delete[] G_daily31GwrSwb;
            G_daily31GwrSwb = NULL;
        }
        if (options.outFswbDaily) {
            delete[] G_daily31Fswb;
            G_daily31Fswb = NULL;
        }
        if (options.outLandAreaFracDaily) {
            delete[] G_daily31LandAreaFrac;
            G_daily31LandAreaFrac = NULL;
        }
        if (options.outGwrunSurfrunDaily) {
            delete[] G_daily31GwrunSurfrun;
            G_daily31GwrunSurfrun = NULL;
        }
    }
    // daily output option 31 end

    // daily output option 365 start
    if ((5 == options.grid_store) || (6 == options.grid_store)) {
        if ((options.outPrecDaily) || ((options.scoutcPrecip) && (2 == options.day_store))) {
            delete[] G_daily365ConsistentPrecip;
            G_daily365ConsistentPrecip = NULL;
        }
        if ((options.outCellRunoffDaily) || ((options.scoutCellRunoffkm3) || (options.scoutCellRunoffmm) && (2 ==
                                                                                                             options.day_store))) {
            delete[] G_daily365CellRunoff;
            G_daily365CellRunoff = NULL;
        }
        if ((options.outCellRunoffDaily) || ((options.scoutCellRunoffmm) && (2 == options.day_store))) {
            delete[] G_daily365CellRunoff_mm;
            G_daily365CellRunoff = NULL;
        }
        if ((options.outCellSurfaceDaily) || ((options.scoutCellSRunoff) && (2 == options.day_store))) {
            delete[] G_daily365CellSurfaceRunoff;
            G_daily365CellSurfaceRunoff = NULL;
        } // FP20160915N002
        // FP20161018N002 Reservoir operation start years: outLandAreaFracDaily instead outLandAreaFrac
        if ((options.outLandAreaFracDaily) || ((options.scoutLandAreaFrac) && (2 == options.day_store))) {
            delete[] G_daily365LandAreaFrac;
            G_daily365LandAreaFrac = NULL;
        }
        if ((options.outGWRunoffDaily) || ((options.scoutGwRunoff) && (2 == options.day_store))) {
            delete[] G_daily365GwRunoff;
            G_daily365GwRunoff = NULL;
        }
        if ((options.outRiverAvailDaily) || ((options.scoutQ) && (2 == options.day_store))) {
            delete[] G_daily365RiverAvail;
            G_daily365RiverAvail = NULL;
        }
        if ((options.outRiverVeloDaily) || ((options.scoutFlowVelo) && (2 == options.day_store))) {
            delete[] G_daily365Velocity;
            G_daily365Velocity = NULL;
        }
        if ((options.outCellAETDaily) || ((options.scoutCellAET) && (2 == options.day_store))) {
            delete[] G_daily365CellAET;
            G_daily365CellAET = NULL;
        }
        if ((options.outSurfStorDaily) || ((options.scoutSurfaceStor) && (2 == options.day_store))) {
            delete[] G_daily365SurfStor;
            G_daily365SurfStor = NULL;
        }
        if ((options.outSingleStoragesDaily) || ((options.scoutLocLake) && (2 == options.day_store))) {
            delete[] G_daily365LocLakeStor;
            G_daily365LocLakeStor = NULL;
        }
        if ((options.outSingleStoragesDaily) || ((options.scoutLocWet) && (2 == options.day_store))) {
            delete[] G_daily365LocWetlStor;
            G_daily365LocWetlStor = NULL;
        }
        if ((options.outSingleStoragesDaily) || ((options.scoutGloLake) && (2 == options.day_store))) {
            delete[] G_daily365GloLakeStor;
            G_daily365GloLakeStor = NULL;
        }
        if ((options.outSingleStoragesDaily) || ((options.scoutGloWet) && (2 == options.day_store))) {
            delete[] G_daily365GloWetlStor;
            G_daily365GloWetlStor = NULL;
        }
        if ((options.outSingleStoragesDaily) || ((options.scoutRiver) && (2 == options.day_store))) {
            delete[] G_daily365RiverStor;
            G_daily365RiverStor = NULL;
        }
        if ((options.outSingleStoragesDaily) || ((options.scoutReservoir) && (2 == options.day_store))) {
            delete[] G_daily365ResStor;
            G_daily365ResStor = NULL;
        }
        if ((options.outSingleStoragesDaily || options.outGWStorageDaily) ||
            ((options.scoutGwStor) && (2 == options.day_store))) {
            delete[] G_daily365GwStor;
            G_daily365GwStor = NULL;
        }
        if ((options.outTotalWaterInStoragesDaily_km3 || options.outTotalWaterInStoragesDaily_mm) ||
            ((options.scoutTWSkm3) || (options.scoutTWSmm) && (2 == options.day_store))) {
            delete[] G_daily365TotalWaterInStorages_km3;
            G_daily365TotalWaterInStorages_km3 = NULL;
        }
        if ((options.outTotalWaterInStoragesDaily_mm) || ((options.scoutTWSmm) && (2 == options.day_store))) {
            delete[] G_daily365TotalWaterInStorages_mm;
            G_daily365TotalWaterInStorages_mm = NULL;
        }
        if ((options.outGwrSwbDaily) || ((options.scoutGwrSwb) && (2 == options.day_store))) {
            delete[] G_daily365GwrSwb;
            G_daily365GwrSwb = NULL;
        }
        if ((options.outFswbDaily) || ((options.scoutFswb) && (2 == options.day_store))) {
            delete[] G_daily365Fswb;
            G_daily365Fswb = NULL;
        }
        if (options.outGwrunSurfrunDaily) {
            delete[] G_daily365GwrunSurfrun;
            G_daily365GwrunSurfrun = NULL;
        }
        if (options.outCellAETWCaDaily) {
            delete[] G_daily365CellAETWCa;
            G_daily365CellAETWCa = NULL;
        }
    }
    // daily output option 365 end
    monthlyDischargeFile.close();
    annualDischargeFile.close();
    if (options.outStationDischargeDaily) dailyDischargeFile.close();
    if (options.outStationVelocityDaily) dailyVelocityFile.close(); //+++ new flow velocity
    if (options.outLocLakeStorageDaily) locLakeFile.close();
    if (options.outLocWetlStorageDaily) locWetlFile.close();
    if (options.outGloLakeStorageDaily) gloLakeFile.close();
    if (options.outResStorageDaily && (options.resOpt == 1))
        gloResFile.close(); // neue OUTPUT-Option
    if (options.outGloWetlStorageDaily) gloWetlFile.close();
}
// end ~routingClass

