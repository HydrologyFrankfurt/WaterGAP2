
/***********************************************************************
*
*many new options and output options by Frank Voss, Stephanie Eisner, Martina
*Weiß, Kristina Fiedler, Heike Hoffmann-Dobrev and Hannes Müller Schmied;
*option settings as an output file by Hannes Müller Schmied
*
* see former changes at file option.cpp.versioninfos.txt
*
***********************************************************************/

#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include "def.h"
#include "timestring.h"  // HMS MODEL_SETTINGS output
#include "stack.h" // HMS 2013-11-21 SINGLECELL output
#include "geo.h" // HMS 2013-11-21 SINGLECELL output
#include "routing.h" // HMS 2013-11-21 SINGLECELL output
#include "land.h" // HMS 2013-11-21 SINGLECELL output

#include "option.h"
using namespace std;

void optionClass::init(int option_c,  char* option_v[])
{
	FILE *file_ptr;
	filename_data		= new char[50];
	filename_time		= new char[50];
	filename_opt		= new char[50];
	filename_out_opt	= new char[50];
	filename_stations	= new char[50];
	filename_singlecells	= new char[50];	         // HMS 2013-11-21 for outputs at predefined IMAGE-Nr
    filename_scout_opt      = new char[50];            // HMS 2013-11-21 output options for single cell out
    char option[6];      // to process options, e.g. "-d", "-t", "-o", "-r", "-s", "-c", "-u"
  	char optionValue[50]; // option value (filenames)
	char string[250];

	// ==== for variable land cover map
	landCoverYearsInList = 0;
	landCoverYears       = NULL;
	landCoverYearsStart  = NULL;
	landCoverYearsEnd    = NULL;

	// default files
	strcpy(filename_data, "DATA.DIR");
	strcpy(filename_time,"TIME.DAT");
	strcpy(filename_opt,"OPTIONS.DAT");
	strcpy(filename_out_opt,"OUTPUT_OPTIONS.DAT");
	strcpy(filename_stations,"STATIONS.DAT");
	strcpy(filename_singlecells,"SINGLECELLS.DAT"); // HMS 2013-11-21 for outputs at predefined IMAGE-Nr
    strcpy(filename_scout_opt,"SINGLECELL_OUTPUT_OPTIONS.DAT");
        //----- check command line for filenames
	// MH20160601N001 Workaround for correct optional reading of option/settings file names as arguments
  	if (option_c) {
        //  Modification Number: MH20160601N001 (part-3)

        //stopped code block : start
        /*
        //cout << "CHECK option_c: '" << option_c << "'" << endl;
    	for (int i=0; i<option_c; i++) {
      		// seperate option and option value
  			strcpy(option, option_v[i]);
            // cout << "CHECK option: '" << option << "'" << endl;
            // cout << "CHECK option_v[i]: '" << option_v[i] << "'" << endl;
			option[2]=0;                                // cut string at 3rd position
  			strcpy(optionValue, option_v[i]+2);
            // cout << "CHECK option_v[i]+2: '" << option_v[i]+2 << "'" << endl;
            // cout << "CHECK option: '" << option << "'" << endl;
            // cout << "CHECK optionValue: '" << optionValue << "'" << endl;
		  	// check DATA.DIR file

        */
 		// end of stopped block

 		// new code block: start

 		// end of new code
        // FP20161018N001 Evaluation of arguments since first argument
        // Do not extract again the parameter file, as this has been done before, possibly without option switch "-p")
        for (int i = 1; i < option_c - 1; i++) {
            strcpy(option, option_v[i]);
            strcpy(optionValue, option_v[i+1]);
 		// end of part-3 of MH20160601N001

 			if (!strcmp(option, "-d")){
				strncpy(filename_data, optionValue, 50);
 			}
		  	// check TIME.DAT file
  			if (!strcmp(option, "-t")){
				strncpy(filename_time, optionValue, 50);
 			}
 		  	// check OPTIONS.DAT file
 			if (!strcmp(option, "-o")){
				strncpy(filename_opt, optionValue, 50);
 			}
 		  	// check OUTPUT_OPTIONS.DAT file
 			if (!strcmp(option, "-r")){
 				strncpy(filename_out_opt, optionValue, 50);
 			}
 		  	// check STATIONS.DAT file
 			if (!strcmp(option, "-s")){
 				strncpy(filename_stations, optionValue, 50);
 			}
 		  	// check SINGLECELLS.DAT file
 			if (!strcmp(option, "-c")){
 				strncpy(filename_singlecells, optionValue, 50);
                        }
                        // check SINGLECELL_OUTPUT_OPTIONS.DAT file
                        if (!strcmp(option, "-u")){
                                strncpy(filename_scout_opt, optionValue, 50);
                        }
		}
  	}

	// the file DATA.DIR contains path names of
	// the directories for input, output,
	// climate input and routing files
	file_ptr = fopen(filename_data, "r");
	if (file_ptr != NULL) {
		printf("Reading %s:\n", filename_data);
		fscanf(file_ptr, "%s", input_dir);
		fscanf(file_ptr, "%s", output_dir);
		fscanf(file_ptr, "%s", climate_dir);
		fscanf(file_ptr, "%s", routing_dir);
		fscanf(file_ptr, "%s", water_use_dir);
		fscanf(file_ptr, "%s", land_cover_dir);
		fclose(file_ptr);
	} else {
		cerr << "Unable to open file:" << filename_data << endl;
		exit(-1);
	}


	// set default options
	// which will be used, when OPTIONS.DAT is not found
	basin 			= 1;
	fileEndianType 	= 2;
	grid_store 		= 0;
	grid_store_TypeForStorages = 0;
	day_store 		= 0;
	cloud 			= 0;
	raindays 		= 0;
	prec_correct 	= 0;
	intercept 		= 0;
	calc_albedo		= 0;
	rout_prepare 	= 1;
	timeStepCheckFlag 		= 1;
	riverveloOpt 	= 0;
	time_series 	= 0;
	subtract_use 	= 0;
	use_alloc 		= 0;
	delayedUseSatisfaction 	= 0;
	petOpt 			= 0;
	clclOpt 		= 0;
	permaOpt 		= 0;
        resOpt   		= 0;
	landCoverOpt    = 0;
	use_kc          = 0;
        statcorrOpt     = 0;
        aridareaOpt     = 0;
        fractionalRoutingOpt     = 0;
        riverEvapoOpt   = 0;
		aggrNUsGloLakResOpt = 1;
		gammaHBV_CorrOpt = 1;
    // FP20161018N003 Enable reading input data from climate land mask
    climate_spatial_resolution = 0;
    // FP20161018N002 Reservoir operation start years
    resYearOpt      = 0;
    resYearReference      = 2000;
    resYearFirstToUse     = 1900;
    resYearLastToUse      = 2010;
    antNatOpt      = 0;
    // FP20161018N004 Net use / net abstraction for reservoir algorithm: new filename
    resNUsMeanYearFirst = 1971;
    resNUsMeanYearLast = 2000;

	// read OPTIONS.DAT
	file_ptr = fopen(filename_opt, "r");
	if (file_ptr != NULL) {
		printf("Reading %s:\n", filename_opt);
		int i = 0, n;

		while (fgets(string, 250, file_ptr))
			if (0 == strncmp(string, "Value:", 6)) {
				i++;
				n = atoi(&string[7]);
				if (1 == i) {
					fileEndianType = n;
					switch (n) {
					case 1:
						printf("   x86 input files.\n");
						break;
					case 2:
						printf("   sparc input files.\n");
						break;
					}
				}
				if (2 == i) {
					basin = n;
					switch (n) {
					case 0:
						printf("   Only superbasins will be calculated.\n");
						break;
					case 1:
						printf("   All waterbasins will be calculated.\n");
						break;
					}
				}
				if (3 == i) {
					grid_store = n;
					switch (n) {
					case 0:
						printf("   Grid files will not be stored.\n");
						break;
					case 1:
						printf("   Grid files will be stored.\n");
						break;
					case 2:
						printf("   Grid files will be stored (including monthly grids).\n");
						break;
					case 3:
						printf("   Grid files will be stored (including daily grids as .31 output for each month).\n");
						break;
					case 4:
						printf("   Grid files will be stored (including daily (as .31 output for each month) and monthly grids).\n");
						break;
					case 5:
						printf("   Grid files will be stored (as .365 output for one year).\n");
						break;
					case 6:
						printf("   Grid files will be stored (including daily (as .365 output for one year) and monthly grids).\n");
						break;
					}
				}
				if (4 == i) {
					grid_store_TypeForStorages = n;
					switch (n) {
					case 0:
						printf("   storage/reservoir volumes at end of time step (*)\n");
						break;
					case 1:
						printf("   mean storage/reservoir volumes (e.g. for worldqual input).\n");
						break;
					}
				}
				if (5 == i) {
					day_store = n;
					switch (n) {
					case 0:
						printf("   Daily values will not be stored.\n");
						break;
					case 1:
						printf("   Daily values will be stored (discharge and water balance for stations defined in STATIONS.DAT).\n");
						break;
					case 2:
						printf("   Daily values will be stored (for IMAGE-Nrs defined in SINGLECELLS.DAT).\n");
						break;
					}
				}

				if (6 == i) {
					time_series = n;
					switch (n) {
					case 0:
						printf("   Daily time series (.31) according to 'TIME.DAT'.\n");
						break;
					case 1:
						printf("   Daily time series (.365) according to 'TIME.DAT'.\n");
						break;
					case 2:
						printf("   Daily precip time series (.31) but monthly time series according to 'TIME.DAT'.\n");
						break;
					case 3:
						printf("   Daily precip time series (.365) but monthly time series according to 'TIME.DAT'.\n");
						break;
					case 4:
						printf("   Monthly time series according to 'TIME.DAT'.\n");
						break;
					case 5:
						printf("   Averages for climate normal period.\n");
						break;
					case 6:
						printf("   Averages for period given in 'TIME.DAT'.\n");
						break;
					}
				}
				if (7 == i) {
					cloud = n;
					switch (n) {
					case 0:
						printf("   Sunshine percentage will be used.\n");
						break;
					case 1:
						printf("   Cloudiness will be used.\n");
						break;
					case 2:
						printf("   Observed incoming short wave radiation and net long wave radiation will be used.\n");
						break;
					case 3:
						printf("   Observed incoming short wave radiation and incoming long wave radiation will be used.\n");
						break;
					case 4:
						printf("   Observed radiation (only incoming short-wave) will be used.\n");
						break;
					}
				}
				if (8 == i) {
					raindays = n;
					switch (n) {
					case 0:
						printf("   Daily rain linearly interpolated from monthly values.\n");
						break;
					case 1:
						printf("   'Number of rain days' will be considered for daily rain.\n");
						break;
					}
				}
				if (9 == i) {
					prec_correct = n;
					switch (n) {
					case 0:
						printf("   Precipitation will not be corrected.\n");
						break;
					case 1:
						printf("   Precipitation will be corrected after Adam & Lettenmaier (usage recommended only for precip input which is not yet corrected, e.g. GPCC).\n");
						break;
					}
				}
				if (10 == i) {
					intercept = n;
					switch (n) {
					case 0:
						printf("   Interception will not be considered.\n");
						break;
					case 1:
						printf("   Interception will be considered.\n");
						break;
					}
				}
				if (11 == i) {
					calc_albedo = n;
					switch (n) {
					case 0:
						printf("   Use constant albedo all over the year.\n");
						break;
					case 1:
						printf("   Use monthly variation in albedo.\n");
						break;
					}
				}
				if (12 == i) {
					petOpt = n;
					switch (n) {
					case 0:
						printf("   PET: Priestley-Taylor.\n");
						break;
					case 1:
						printf("   PET: Penman-Monteith.\n");
						break;
					case 2:
						printf("   PET: Penman-Monteith with altitude-dep. pressure and e_s with temp. range.\n");
						break;
					case 3: 		//mw-H
						printf("   PET: Hargreaves \n");
						break;
					case 4: 		//mw-T
						printf("   PET: combination of Turc and Doorenbos&Pruitt \n");
						break;
					case 5: 		//mw-KP
						printf("   PET: Kimberly-Penman \n");
						break;
					case 6: 		//mw-PT humid/arid ueber RH
						printf("   PET: PT (arid areas based on RH) \n");
						break;
					}
				}
				if (13 == i) {
					use_kc = n;
					switch (n) {
					case 0:
						printf("   Crop coefficients will not be used for estimation of PET.\n");
						break;
					case 1:
						printf("   Crop coefficients will be used for estimation of PET.\n");
						break;
					}
				}
				if (14 == i) {
					landCoverOpt = n;
					switch (n) {
					case 0:
						printf("   const land cover map will be used.\n");
						break;
					case 1:
						printf("   variable land cover map will be used.\n");
						break;
					}
				}
				if (15 == i) {
					rout_prepare = n;
					switch (n) {
					case 0:
						printf("   Old routing files will be used.\n");
						break;
					case 1:
						printf("   Routing files will be prepared again.\n");
						break;
					}
				}
				if (16 == i) {
					timeStepCheckFlag = n;
					switch (n) {
					case 0:
						printf("   Routing time step will not be checked.\n");
						break;
					case 1:
						printf("   Routing time step will be checked.\n");
						break;
					}
				}
				if (17 == i) {
					riverveloOpt = n;
					switch (n) {
					case 0:
						printf("   Use constant flow velocity (according to ROUTING.DAT).\n");
						break;
					case 1:
						printf("   Bankfull-based flow velocity will be calculated.\n");
						break;
                    case 2:
                        printf("   Storage-based flow velocity will be calculated.\n");
                        break;
					}
				}
				if (18 == i) {
					subtract_use = n;
					switch (n) {
					case 0:
						printf("   Water use will not be considered.\n");
						break;
					case 1:
						printf
							("   Water use will be subtracted from river/lake storages (Daily irrig. values).\n");
						break;
					case 2:
						printf
							("   Water use will be subtracted from river/lake storages (Monthly irrig. values).\n");
						break;
					case 3:
						printf
							("   Water use will be subtracted from river/lake storages (Reallocated values).\n");
						break;
					}
				}
				if (19 == i) {
					use_alloc = n;
					switch (n) {
					case 0:
						printf("   Water use is satisfied only by the cells themselfs.\n");
						break;
					case 1:
						printf("   A second cell is allowed (from neighbouring cells).\n");
						break;
					}
				}
				if (20 == i) {
					delayedUseSatisfaction = n;
					switch (n) {
					case 0:
						printf("   Delayed satisfaction of water use not allowed.\n");
						break;
					case 1:
						printf("   Delayed satisfaction of water use allowed.\n");
						break;
					}
				}
				if (21 == i) {
					clclOpt = n;
					switch (n) {
					case 0:
						printf("   Use arid/humid definition based on G_ARID_HUMID.UNF2.\n");
						break;
					case 1:
						printf("   Use KOEPPEN map for arid/humid definition.\n");
						break;
					}
				}
				if (22 == i) {
					permaOpt = n;
					switch (n) {
					case 0:
						printf("   Areal extend of permafrost will not be calculated.\n");
						break;
					case 1:
						printf("   Areal extend of permafrost will be calculated.\n");
						break;
					}
				}

				if (23 == i) {
					resOpt = n;
					switch (n) {
					case 0:
						printf("   Reservoirs will be treated as global lakes.\n");
						break;
					case 1:
						printf("   New reservoir algorithm will be used.\n");
						break;
					}
				}

                                if (24 == i) {
                                        statcorrOpt = n;
                                        switch (n) {
                                        case 0:
                                                printf("   Station correction will be switched off.\n");
                                                break;
                                        case 1:
                                                printf("   Station correction will be used.\n");
                                                break;
                                        }
                                }

                                if (25 == i) {
                                        aridareaOpt = n;
                                        switch (n) {
                                        case 0:
                                                printf("   No groundwater recharge below surface water bodies in arid / semi-arid regions.\n");
                                                break;
                                        case 1:
                                                printf("   Additional groundwater recharge below surface water bodies in arid / semi-arid regions.\n");
                                                break;
                                        }
                                }

                                if (26 == i) {
                                        fractionalRoutingOpt = n;
                                        switch (n) {
                                        case 0:
                                                printf("   No modified routing of groundwater outflow and fast runoff.\n");
                                                break;
                                        case 1:
                                                printf("   Groundwater outflow and fast runoff will be routed in dependance on current size of surface water bodies.\n");
                                                break;
                                        }
                                }
								if (27 == i) {
                                        riverEvapoOpt = n;
                                        switch (n) {
                                        case 0:
                                                printf("   No evaporation from rivers.\n");
                                                break;
                                        case 1:
                                                printf("   Evaporation from rivers is allowed.\n");
                                                break;
                                        }
                                }
								if (28 == i) {
                                        aggrNUsGloLakResOpt = n;
                                        switch (n) {
                                        case 0:
                                                printf("   NUs of riparian cells of global lakes/reservoirs is not taken into account.\n");
                                                break;
                                        case 1:
                                                printf("   NUs of riparian cells of global lakes/reservoirs is taken into account.\n");
                                                break;
                                        }
								}
								if (29 == i) {
                                        gammaHBV_CorrOpt = n;
                                        switch (n) {
                                        case 0:
                                                printf("   G_gammaHBV is not corrected.\n");
                                                break;
                                        case 1:
                                                printf("   G_gammaHBV is corrected.\n");
                                                break;
                                        }
								}

                                // FP20161018N003 Enable reading input data from climate land mask
                                if (30 == i) {
                                    climate_spatial_resolution = n;
                                    switch (n) {
                                    case 0:
                                        printf("   Climate will be read from traditional WaterGAP land mask (SLM or WLM).\n");
                                        break;
                                    case 1:
                                        printf("   Climate will be read from climate land mask.\n");
                                        break;
                                    }
                                }

                                // FP20161018N001 Reservoir start operation years
                                if (31 == i) {
                                        resYearOpt = n;
                                        switch (n) {
                                        case 0:
                                                printf("   Reservoir operation start year is not considered.\n");
                                                break;
                                        case 1:
                                                printf("   Reservoirs are becoming active in their start operation year (*)\n");
                                                break;
                                        }
                                }

                                if (32 == i) {
                                        resYearReference = n;
                                        printf("   If reservoir operation start year is not considered (resYearOpt == 0), then use reservoir status from reference year %d).\n", resYearReference);
                                }

                                if (33 == i) {
                                        resYearFirstToUse = n;
                                        printf("   If reservoir operation start year is considered (resYearOpt == 1), use this year as first status of reservoirs: %d).\n", resYearFirstToUse);
                                }

                                if (34 == i) {
                                        resYearLastToUse = n;
                                        printf("   If reservoir operation start year is considered (resYearOpt == 1), use this year as last status of reservoirs: %d).\n", resYearLastToUse);
                                }

                                if (35 == i) {
                                        antNatOpt = n;
                                        switch (n) {
                                        case 0:
                                                printf("   Anthropogenic run (reservoirs are considered (ResOpt==1) or treated as global lakes (ResOpt==0) (*).\n");
                                                break;
                                        case 1:
                                                printf("   Naturalized run (no consideration of reservoirs and water use): set subtract_use = 0.\n");
                                                subtract_use = 0;
                                                break;
                                        }
                                }

                                if (36 == i) {
                                        resNUsMeanYearFirst = n;
                                        printf("   First year %d of average net abstraction from surface water NUs for reservoir algorithm calculation of mean demand (e.g. G_NUs_1971_2000.UNF0, with option 23: resOpt == 1)\n", resNUsMeanYearFirst);
                                }

                                if (37 == i) {
                                        resNUsMeanYearLast = n;
                                        printf("   Last year %d of average net abstraction from surface water NUs for reservoir algorithm calculation of mean demand (e.g. G_NUs_1971_2000.UNF0, with option 23: resOpt == 1)\n", resNUsMeanYearLast);
                                }


			}
            // end if (string search for "Value")
        // end while (line reading from option file)
		fclose(file_ptr);
	} else {
		printf("\nError while opening file '%s'.\n", filename_opt);
		printf("Using default options.\n");
	}
    if (((0 == time_series) || (1 == time_series)) && (cloud < 2)) {
        cerr << "Wrong setup in OPTIONS.DAT: Sunshine percentage or cloudiness can only be used in combination with monthly climate input.\n" << endl;
        cerr << "(((0 == time_series) || (1 == time_series)) && (cloud < 2))" << endl;
        exit (-1);
    }
    if ((1 == use_kc) && (1 == calc_albedo)) {
        cerr << "Wrong setup in OPTIONS.DAT: Crop coefficients (Kc) cannot be used in combination with variable albedo.\n" << endl;
        cerr << "((1 == use_kc) && (1 == calc_albedo))" << endl;
        exit (-1);
    }

    // FP20161018N002 Reservoir operation start years
    if (resYearFirstToUse > resYearLastToUse) {
        cerr << "Wrong setup in OPTIONS.DAT: Reservoir operation last year to use: " << resYearLastToUse << " is smaller than first year to use: " << resYearFirstToUse << endl;
        cerr << "Please change accordingly: last year AFTER first year" << endl;
        cerr << "(resYearFirstToUse > resYearLastToUse)" << endl;
        exit (-1);
    }
    if ( (resYearReference > resYearLastToUse) || (resYearReference < resYearFirstToUse) ) {
        cerr << "Wrong setup in OPTIONS.DAT: Reservoir operation reference year " << resYearReference << " is outside limits of last year to use: " << resYearLastToUse << " and first year to use: " << resYearFirstToUse;
        cerr << "Please change reference year or create files for last or first year" << endl;
        cerr << "( (resYearReference > resYearLastToUse) || (resYearReference < resYearFirstToUse) )" << endl;
        exit (-1);
    }

    // FP20161018N003 Enable reading input data from climate land mask
    if ((6 == time_series) && (1 == climate_spatial_resolution)) {
        cerr << "Wrong setup in OPTIONS.DAT: Calculation of averages from climate land mask not yet implemented, exit program execution.\n" << endl;
        cerr << "((6 == time_series) && (1 == climate_spatial_resolution))" << endl;
        exit (-1);
    }

	// read TIME.DAT
	if (time_series != 5) {
		file_ptr = fopen(filename_time, "r");
		if (file_ptr != NULL) {
			printf("Reading %s:\n", filename_time);
			fgets(string, 250, file_ptr);
			fscanf(file_ptr, "%hi", &start_year);
			fgets(string, 250, file_ptr);
			fscanf(file_ptr, "%hi", &end_year);
			fgets(string, 250, file_ptr);
			fscanf(file_ptr, "%hi", &evalStartYear);
			fgets(string, 250, file_ptr);
			fscanf(file_ptr, "%hi", &init_years);
			fgets(string, 250, file_ptr);
			fscanf(file_ptr, "%hi", &time_step);
			fgets(string, 250, file_ptr);

			// ==== for variable land cover map
			if (landCoverOpt) {

				fscanf(file_ptr, "%hi", &landCoverYearsInList);
				if (landCoverYearsInList==0) {
					cerr << "no year list for  variable land cover map!" << endl;
					exit (-1);
				}
				fgets(string, 250, file_ptr);

				landCoverYears      = new short[landCoverYearsInList];
				landCoverYearsStart = new short[landCoverYearsInList];
				landCoverYearsEnd   = new short[landCoverYearsInList];

				for (int t_step=0; t_step<landCoverYearsInList; t_step++)
					fscanf(file_ptr, "%hi", &landCoverYears[t_step]);

				for (int t_step=0; t_step<landCoverYearsInList; t_step++) {
					if (t_step==0)
						if (landCoverYears[t_step]<start_year) landCoverYearsStart[t_step] = landCoverYears[t_step];
						else landCoverYearsStart[t_step] = start_year;
					else
						landCoverYearsStart[t_step] = landCoverYearsEnd[t_step-1]+1;

					if (t_step==landCoverYearsInList-1)
						landCoverYearsEnd[t_step]=end_year;
					else
						landCoverYearsEnd[t_step]=(int)((landCoverYears[t_step+1]-landCoverYears[t_step])/2+0.5)+landCoverYears[t_step]-1;
				} // for(year)

			}

			// ==== for variable land cover map

			fclose(file_ptr);
		} else {
			printf("File '%s' not found! \n", filename_time);
			start_year = 1961;
			end_year = 1990;
			evalStartYear = 1956;
			init_years = 0;
			time_step = 1;
			landCoverOpt = 0; // ohne Liste von Jahren kann man keine Karten einlesen
			cout <<"   const land cover map be used.\n";
		}
	} else {
		start_year = 1;
		end_year = 1;
		evalStartYear = 1;
		init_years = 5;
		time_step = 1;
		landCoverOpt = 0; // ohne Liste von Jahren kann man keine Karten einlesen
		cout <<"   const land cover map will be used.\n";
	}
	printf("   Start year for calculations: %d\n", start_year);
	printf("   Start year for evaluations:  %d\n", evalStartYear);
	printf("   End year   : %d\n", end_year);
	printf("   Init. years: %d\n", init_years);
	printf("   Time step  : %d\n", time_step);

	if (landCoverOpt) {
		printf("   years for land cover map  : %d\n", landCoverYearsInList);
		printf("   year\t:\tfrom\t-\tto\n");
		for (int t_step=0; t_step<landCoverYearsInList; t_step++)
			printf("   %d\t:\t%d\t-\t%d\n", landCoverYears[t_step], landCoverYearsStart[t_step], landCoverYearsEnd[t_step]);
	}

        if (6 == time_series)
		mid_year = evalStartYear + (end_year - start_year) / 2;
	else
		mid_year = -99;

	// read OUTPUT_OPTIONS.DAT
	// set default output options
	// monthly binary files
	outPrec                    = false;                        // 0
	outAET                     = false;                        // 1
	outCellAET                 = false;                        // 2
	outCellRunoff              = false;                        // 3
        outPotCellRunoffAnnual     = false;                        // 4
	outCellSurface             = false;                        // 5
	outRunoff                  = false;                        // 6
	outUrbanRunoff             = false;                        // 7
	outSurfaceRunoff           = false;                        // 8
	outGWRunoff                = false;                        // 9
	outGWRecharge              = false;                        // 10
	outSingleStorages          = false;                        // 11
	//outGWStorage               = false;                        // 12 included in 11
        outMinMaxRiverAvail	   = false;                        // 12
	outSoilWater               = false;                        // 13
	outLAI                     = false;                        // 14
	outAlbedo                  = false;                        // 15
	outInterception            = false;                        // 16
	outCanopyWater             = false;                        // 17
	outmaxCanopyWater          = false;                        // 18
        outRiverPET                = false;                        // 19
	outNetShortWaveRadiation   = false;                        // 20
	outNetLongWaveRadiation    = false;                        // 21
	outNetRadiation            = false;                        // 22
	outOpenWaterPET            = false;                        // 23
	outOpenWaterEvap           = false;                        // 24
	outPET                     = false;                        // 25
	outTotalPET                = false;                        // 26
	outRiverAvail              = false;                        // 27
	outRiverInUpstream         = false;                        // 28
	outRiverVelo               = false;                        // 29
	outSnowCover               = false;                        // 30
	outSWE                     = false;                        // 31
	outSnowFall                = false;                        // 32
	outSnowMelt                = false;                        // 33
	outSnowEvap                = false;                        // 34
	outSurfStor                = false;                        // 35
        outLocLakeStorage				   = false;                        // 36 x  included in 11
        outLocWetlStorage					   = false;                        // 37 x  included in 11
	outGloLakeStorage		   = false;                        // 38 x  included in 11
	outGloWetlStorage		   = false;                        // 39 x  included in 11
	outResStorage		       = false;                        // 40 x  included in 11
	outRiverStorage			   = false;                        // 41 x  included in 11
    outTotalWaterInStorages_km3		= false;                        // 42 x
	outTotalWaterInStorages_mm 	= false;                        // 43 x
	outActualNAS               = false;                        // 44
	outActualNAG               = false;                        // 45
    outWCa                     = false;                        // 46
    outWCaInCellUsed           = false;                        // 47
	outAllocUsein2ndCell                = false;                        // 48
    outSatisAllocUsein2ndCell = false;                         // 49
	outConsWaterUse            = false;                        // 50
	outUnsatUseSW                = false;                        // 51
	outUnsatUseSWPrev            = false;                        // 52
        outGwrSwb            = false;                        // 53
        outFswb            = false;                        // 54
        outLandAreaFrac            = false;                        // 55
        outTemperature            = false;                        // 56
        outSunshine            = false;                        // 57
        outCellAETWCa              = false;                        // 58
        outLocWetlExt              = false;                        // 59
        outGloWetlExt              = false;                        // 60
        outTotalGWR                = false;                        // 61
	// daily binary files
        outPrecDaily               = false;                        // 62
        outCellAETDaily            = false;                        // 63
        outCellRunoffDaily         = false;                        // 64
        outCellSurfaceDaily        = false;                        // 65
        outSurfaceRunoffDaily      = false;                        // 66
        outGWRunoffDaily           = false;                        // 67
        outGWRechargeDaily         = false;                        // 68
        outGWStorageDaily          = false;                        // 69
        outSoilWaterDaily          = false;                        // 70
        outLAIDaily                = false;                        // 71
        outAlbedoDaily             = false;                        // 72
        outInterceptionDaily       = false;                        // 73
        outCanopyWaterDaily        = false;                        // 74
        outmaxCanopyWaterDaily     = false;                        // 75
        outExtRadiationDaily       = false;                        // 76
        outShortDownRadiationDaily = false;                        // 77
        outShortUpRadiationDaily   = false;                        // 78
        outNetShortWaveRadiationDaily = false;                     // 79
        outLongDownRadiationDaily  = false;                        // 80
        outLongUpRadiationDaily    = false;                        // 81
        outNetLongWaveRadiationDaily = false;                      // 82
        outNetRadiationDaily       = false;                        // 83
        outPETDaily                = false;                        // 84
        outTotalPETDaily           = false;                        // 85
        outRiverAvailDaily         = false;                        // 86
        outRiverVeloDaily          = false;                        // 87
        outSnowCoverDaily          = false;                        // 88
        outSWEDaily                = false;                        // 89
        outSnowFallDaily           = false;                        // 90
        outSnowMeltDaily           = false;                        // 91
        outSurfStorDaily           = false;                        // 92
    outSingleStoragesDaily     = false;                            // 93
    outTotalWaterInStoragesDaily_km3 = false;                      // 94
    outTotalWaterInStoragesDaily_mm = false;                        // 95
    outTotalWaterInStoragesStartEndDaily_km3    = false;            //96
        outGwrSwbDaily				= false;	   // 97
        outFswbDaily				= false;	   // 98
        outLandAreaFracDaily			= false;	   // 99
        outGwrunSurfrunDaily        = false;                        // 100
        outCellAETWCaDaily          = false;                        // 101
	// additional binary file
    outGWFactor                = false;                        // 102
    outRGmax                   = false;                        // 103
    outmaxSoilWater            = false;                        // 104
    outRootingDepth            = false;                        // 105
    outClcl                    = false;                        // 106
    outLAImax                  = false;                        // 107
    outAirFrost                = false;                        // 108
    outSurfaceFrost            = false;                        // 109
    outRH                      = false;                        // 110
	// ASCII files
    outAllUpStations           = false;                        // 111
    outDirectUpStations        = false;                        // 112
    outRainDays                = false;                        // 113
    outStationDischargeAnnual  = false;                        // 114
    outStationDischargeMonthly = false;                        // 115
    outStationList             = false;                        // 116 : not implemented yet
    outSuperbasinClimate       = false;                        // 117
    outResvoirMonthly          = false;                        // 118 x
	// additional ASCII files (option: save daily values)
    outDailyValues             = false;                        // 119
    outDailyInterception       = false;                        // 120
    outStationDischargeDaily   = false;                        // 121
    outStationVelocityDaily    = false;                        // 122
    outLocLakeStorageDaily     = false;                        // 123
    outGloLakeStorageDaily     = false;                        // 124
    outLocWetlStorageDaily     = false;                        // 125
    outGloWetlStorageDaily     = false;                        // 126
    outResStorageDaily         = false;                        // 127 x

    scoutTemp                   = false;                        //sc0
    scoutExtRad                 = false;                        //sc1
    scoutShortDownRad             = false;                        //sc2
    scoutAlbedo                 = false;                        //sc3
    scoutShortUpRad            = false;                        //sc4
    scoutNetShortRad            = false;                        //sc5
    scoutLongDownRad              = false;                        //sc6
    scoutNetLongRad             = false;                        //sc7
    scoutLongUpRad             = false;                        //sc8
    scoutNetRad                 = false;                        //sc9
    scoutLAI                    = false;                        //sc10
    scoutKc                     = false;                        //sc11
    scoutLandPET                = false;                        //sc12
    scoutCellPET                = false;                        //sc13
    scoutPrecip                 = false;                        //sc14
    scoutcPrecip                = false;                        //sc15
    scoutCanopyWater            = false;                        //sc16
    scoutmaxCanopyWater         = false;                        //sc17
    scoutInterception           = false;                        //sc18
    scoutSnowfall               = false;                        //sc19
    scoutSnowWater              = false;                        //sc20
    scoutSnowCov                = false;                        //sc21
    scoutSnowmelt               = false;                        //sc22
    scoutSoilWater              = false;                        //sc23
    scoutSurfaceRunoff          = false;                        //sc24
    scoutGwRunoff               = false;                        //sc25
    scoutGwRecharge             = false;                        //sc26
    scoutCellAET                = false;                        //sc27
    scoutCellRunoffkm3          = false;                        //sc28
    scoutCellRunoffmm           = false;                        //sc29
    scoutCellSRunoff            = false;                        //sc30
    scoutQ                      = false;                        //sc31
    scoutFlowVelo               = false;                        //sc32
    scoutLocLake                = false;                        //sc33
    scoutLocWet                 = false;                        //sc34
    scoutGloLake                = false;                        //sc35
    scoutReservoir              = false;                        //sc36
    scoutGloWet                 = false;                        //sc37
    scoutRiver                  = false;                        //sc38
    scoutSurfaceStor            = false;                        //sc39
    scoutTWSkm3                 = false;                        //sc40
    scoutTWSmm                  = false;                        //sc41
    scoutGwStor                 = false;                        //sc42
    scoutCanopyStor             = false;                        //sc43
    scoutSoilStor               = false;                        //sc44
    scoutSnowStor               = false;                        //sc45
    scoutGwrSwb                 = false;                        //sc46
    scoutFswb                   = false;                        //sc47
    scoutLandAreaFrac           = false;                        //sc48

    ifstream outOptionsFile (filename_out_opt);
	if (outOptionsFile) {
		cout << "Reading " << filename_out_opt <<": ";
                char templine[150];
                bool setting[150];
                short set_no = 0;
		short notset_no = 0;
		// read file
		// ignore lines that do not start with '1' or '0'
                while (outOptionsFile.getline(templine, 150)) {
			if (templine[0] == '0') {
				setting[notset_no + set_no] = false;
				notset_no ++;
			}
			if (templine[0] == '1') {
				setting[notset_no + set_no] = true;
				set_no ++;
			}
		}
        if (set_no + notset_no != 128) {
			cerr << "Unexpected number of settings in file '" << filename_out_opt << "'!" << endl;
			exit (1);
		}

		// settings for output files
		// monthly/yearly binary files
		outPrec						= setting[0];
		outAET 						= setting[1];
		outCellAET					= setting[2];
		outCellRunoff		  		= setting[3];
        outPotCellRunoffAnnual  	= setting[4];
		outCellSurface		  		= setting[5];
		outRunoff 					= setting[6];
		outUrbanRunoff				= setting[7];
		outSurfaceRunoff			= setting[8];
		outGWRunoff 				= setting[9];
		outGWRecharge 				= setting[10];
		outSingleStorages 			= setting[11];
		//outGWStorage 				= setting[12];
		outMinMaxRiverAvail			= setting[12];
		outSoilWater 				= setting[13];
        outLAI			 			= setting[14];
        outAlbedo					= setting[15];
		outInterception 			= setting[16];
		outCanopyWater 				= setting[17];
		outmaxCanopyWater			= setting[18];
		outRiverPET                 = setting[19];
		outNetShortWaveRadiation    = setting[20];
		outNetLongWaveRadiation     = setting[21];
		outNetRadiation 			= setting[22];
		outOpenWaterPET 			= setting[23];
		outOpenWaterEvap 			= setting[24];
		outPET  					= setting[25];
		outTotalPET 				= setting[26];
		outRiverAvail 				= setting[27];
		outRiverInUpstream			= setting[28];  //
		outRiverVelo 				= setting[29];
		outSnowCover 				= setting[30];
		outSWE		 	    		= setting[31];
		outSnowFall 				= setting[32];
		outSnowMelt 				= setting[33];
		outSnowEvap 				= setting[34];
		outSurfStor					= setting[35];
        outLocLakeStorage		    = setting[36];
        outLocWetlStorage			= setting[37];
        outGloLakeStorage    		= setting[38];
        outGloWetlStorage    		= setting[39];
        outResStorage        		= setting[40];
		outRiverStorage      		= setting[41];  //
		outTotalWaterInStorages_km3 = setting[42];  //
		outTotalWaterInStorages_mm  = setting[43];  //
		outActualNAS 				= setting[44];
		outActualNAG 				= setting[45];
        outWCa                      = setting[46];
        outWCaInCellUsed            = setting[47];
		outAllocUsein2ndCell		= setting[48];
        outSatisAllocUsein2ndCell	= setting[49];
		outConsWaterUse 			= setting[50];
		outUnsatUseSW 				= setting[51];
		outUnsatUseSWPrev 			= setting[52];
        outGwrSwb 			        = setting[53];
        outFswb 	                = setting[54];
        outLandAreaFrac 			= setting[55];
        outTemperature  			= setting[56];
        outSunshine     			= setting[57];
        outCellAETWCa               = setting[58];
        outLocWetlExt               = setting[59];
        outGloWetlExt               = setting[60];
        outTotalGWR                 = setting[61];
		// daily binary files
        outPrecDaily				= setting[62];
        outCellAETDaily 			= setting[63];
        outCellRunoffDaily 			= setting[64];
        outCellSurfaceDaily			= setting[65];
        outSurfaceRunoffDaily 		= setting[66];
        outGWRunoffDaily 			= setting[67];
        outGWRechargeDaily 			= setting[68];
        outGWStorageDaily 			= setting[69];
        outSoilWaterDaily 			= setting[70];
        outLAIDaily		 			= setting[71];
        outAlbedoDaily				= setting[72];
        outInterceptionDaily 		= setting[73];
        outCanopyWaterDaily			= setting[74];
        outmaxCanopyWaterDaily		= setting[75];
        outExtRadiationDaily 		= setting[76];
        outShortDownRadiationDaily 	= setting[77];
        outShortUpRadiationDaily 	= setting[78];
        outNetShortWaveRadiationDaily = setting[79];
        outLongDownRadiationDaily 	= setting[80];
        outLongUpRadiationDaily 	= setting[81];
        outNetLongWaveRadiationDaily= setting[82];
        outNetRadiationDaily 		= setting[83];
        outPETDaily 				= setting[84];
        outTotalPETDaily 			= setting[85];
        outRiverAvailDaily 			= setting[86];
        outRiverVeloDaily			= setting[87];
        outSnowCoverDaily 			= setting[88];
        outSWEDaily 				= setting[89];
        outSnowFallDaily 			= setting[90];
        outSnowMeltDaily 			= setting[91];
        outSurfStorDaily 			= setting[92];
        outSingleStoragesDaily 		= setting[93];
        outTotalWaterInStoragesDaily_km3 = setting[94];
        outTotalWaterInStoragesDaily_mm = setting[95];
        outTotalWaterInStoragesStartEndDaily_km3 = setting[96];
        outGwrSwbDaily	            = setting[97];
        outFswbDaily		        = setting[98];
        outLandAreaFracDaily        = setting[99];
        outGwrunSurfrunDaily        = setting[100];
        outCellAETWCaDaily          = setting[101];

                // additional binary files
        outGWFactor 				= setting[102];
        outRGmax 					= setting[103];
        outmaxSoilWater				= setting[104];
        outRootingDepth				= setting[105];
        outClcl						= setting[106];
        outLAImax					= setting[107];
        outAirFrost					= setting[108];
        outSurfaceFrost				= setting[109];
        outRH						= setting[110];

        // ASCII files
        outAllUpStations 			= setting[111];
        outDirectUpStations 		= setting[112];
        outRainDays 				= setting[113];
        outStationDischargeAnnual 	= setting[114];
        outStationDischargeMonthly 	= setting[115];
        outStationList 				= setting[116]; //not yet implemented!!!
        outSuperbasinClimate 		= setting[117]; // Superbasin calculation currently obsolete, probably not working correctly in WaterGAP2.2b, do in postprocessing  // FP 2015-10
        outResvoirMonthly           = setting[118];  //
        // additional ASCII files (option: save daily values)
        outDailyValues 				= setting[119];
        outDailyInterception 		= setting[120];
        outStationDischargeDaily 	= setting[121];
        outStationVelocityDaily 	= setting[122];
        outLocLakeStorageDaily		= setting[123];
        outGloLakeStorageDaily		= setting[124];
        outLocWetlStorageDaily		= setting[125];
        outGloWetlStorageDaily		= setting[126];
        outResStorageDaily          = setting[127];  //

		cout << set_no << " files selected." <<endl;
	}
	else {
		printf("File '%s' not found! No output will be written!\n", filename_out_opt);
	}
	outOptionsFile.close();

        // do the same for SINGLECELL_OUTPUT_OPTIONS
        ifstream outscOptionsFile (filename_scout_opt);
        if (outscOptionsFile) {
                cout << "Reading " << filename_scout_opt <<": ";
                char templine[120];
                bool scsetting[100];
                short set_no = 0;
                short notset_no = 0;
                // read file
                // ignore lines that do not start with '1' or '0'
                while (outscOptionsFile.getline(templine, 120)) {
                        if (templine[0] == '0') {
                                scsetting[notset_no + set_no] = false;
                                notset_no ++;
                        }
                        if (templine[0] == '1') {
                                scsetting[notset_no + set_no] = true;
                                set_no ++;
                        }
                }
        if (set_no + notset_no != 49) {
                        cerr << "Unexpected number of settings in file '" << filename_scout_opt << "'!" << endl;
                        exit (1);
                }

        if (day_store < 2) {// set all scout options to 0 to avoid unnecessary initialization.
            for (int i = 0; i < (set_no + notset_no); i++) {
                scsetting[i] = 0;
            }
        }
        else {
                // settings for output files
                // monthly/yearly binary files
                scoutTemp                   = scsetting[0];
                scoutExtRad                 = scsetting[1];
                scoutShortDownRad                = scsetting[2];
                scoutAlbedo                 = scsetting[3];
                scoutShortUpRad                = scsetting[4];
                scoutNetShortRad            = scsetting[5];
                scoutLongDownRad              = scsetting[6];
                scoutNetLongRad             = scsetting[7];
                scoutLongUpRad             = scsetting[8];
                scoutNetRad                 = scsetting[9];
                scoutLAI                    = scsetting[10];
                scoutKc                     = scsetting[11];
                scoutLandPET                = scsetting[12];
                scoutCellPET                = scsetting[13];
                scoutPrecip                 = scsetting[14];
                scoutcPrecip                = scsetting[15];
                scoutCanopyWater            = scsetting[16];
                scoutmaxCanopyWater         = scsetting[17];
                scoutInterception           = scsetting[18];
                scoutSnowfall               = scsetting[19];
                scoutSnowWater              = scsetting[20];
                scoutSnowCov                = scsetting[21];
                scoutSnowmelt               = scsetting[22];
                scoutSoilWater              = scsetting[23];
                scoutSurfaceRunoff          = scsetting[24];
                scoutGwRunoff               = scsetting[25];
                scoutGwRecharge             = scsetting[26];
                scoutCellAET                = scsetting[27];
                scoutCellRunoffkm3          = scsetting[28];
                scoutCellRunoffmm           = scsetting[29];
                scoutCellSRunoff            = scsetting[30];
                scoutQ                      = scsetting[31];
                scoutFlowVelo               = scsetting[32];
                scoutLocLake                = scsetting[33];
                scoutLocWet                 = scsetting[34];
                scoutGloLake                = scsetting[35];
                scoutReservoir              = scsetting[36];
                scoutGloWet                 = scsetting[37];
                scoutRiver                  = scsetting[38];
                scoutSurfaceStor            = scsetting[39];
                scoutTWSkm3                 = scsetting[40];
                scoutTWSmm                  = scsetting[41];
                scoutGwStor                 = scsetting[42];
                scoutCanopyStor             = scsetting[43];
                scoutSnowStor               = scsetting[44];
                scoutSoilStor               = scsetting[45];
                scoutGwrSwb                 = scsetting[46];
                scoutFswb                   = scsetting[47];
                scoutLandAreaFrac           = scsetting[48];

                if (2 == day_store)
                    cout << set_no << " variables for SINGLECELLOUT selected." <<endl;
                else
                    cout << "no SINGLECELLOUT selected." << endl;
        }
        }
        else {
                printf("File '%s' not found! No output will be written!\n", filename_scout_opt);
        }
        outscOptionsFile.close();
}
   // MODEL_SETTINGS output start
void optionClass::createModelSettingsFile(const char *outputDir)
{

  char filename[250];

	sprintf(filename, "%s/MODEL_SETTINGS.OUT", outputDir);
	ofstream output_file(filename);

	if (!output_file) {
		cerr << "Can not open file " << filename << " for writing" << endl;
		exit(-1);
	}
        output_file << "# Model time: " << getTimeString() << filename << endl;
	output_file << "------------------------------------------------\n";

        output_file << "Values from file " << filename_data << endl;
	output_file << "Input data directory: " << input_dir << endl;
  output_file << "Output data directory: " << output_dir << endl;
	output_file << "Climate data directory: " << climate_dir << endl;
	output_file << "Routing data directory: " << routing_dir << endl;
	output_file << "Water use data directory: " << water_use_dir << endl;
        output_file << "Land cover data directory: " << land_cover_dir << endl;
	output_file << "------------------------------------------------\n";

        output_file << "Values from file " << filename_time << endl;
	output_file << "Start year of calculations: " << start_year << endl;
	output_file << "End year of simulation: " << end_year << endl;
	output_file << "Start year of evaluation: " << evalStartYear << endl;
	output_file << "Number of initialisation years: " << init_years << endl;
        output_file << "Time step: " << time_step << endl;
	output_file << "------------------------------------------------\n";

        output_file << "Values from file " << filename_opt << endl;
	output_file << "1. Indicator for endian type of binary input files    : " << fileEndianType << endl;
  output_file << "2. Indicator for waterbasins to be simulated    : " << basin << endl;
  output_file << "3. Indicator for storage of grid files    : " << grid_store << endl;
  output_file << "4. Indicator for type of storage/reservoir volumes    : " << grid_store_TypeForStorages << endl;
  output_file << "5. Indicator for storage of additional daily values    : " << day_store << endl;
  output_file << "6. Indicator for calculation of time series    : " << time_series << endl;
  output_file << "7. Indicator for the use of radiation input type    : " << cloud << endl;
  output_file << "8. Calculation of daily precipitation     : " << raindays << endl;
  output_file << "9. Indicator for correction of precipitation values    : " << prec_correct << endl;
  output_file << "10. Indicator for interception    : " << intercept << endl;
  output_file << "11. Calculate variation in albedo    : " << calc_albedo << endl;
  output_file << "12. Equation for potential evapotranspiration (PET)    : " << petOpt << endl;
  output_file << "13. Consider crop coefficients (Kc) for estimation of potential evapotransporation (PET)    : " << use_kc << endl;
  output_file << "14. Calculations with variable land cover    : " << landCoverOpt << endl;
  output_file << "15. Preparation of routing files    : " << rout_prepare << endl;
  output_file << "16. Check routing time step    : " << timeStepCheckFlag << endl;
  output_file << "17. Calculate variable flow velocity    : " << riverveloOpt << endl;
  output_file << "18. Indicator for subtraction of water use from river and lake storages    : " << subtract_use << endl;
  output_file << "19. Indicator for allocation of water use    : " << use_alloc << endl;
  output_file << "20. Allow delayed satisfaction of water use from surface water    : " << delayedUseSatisfaction << endl;
  output_file << "21. Calculate arid/humid areas    : " << clclOpt << endl;
  output_file << "22. Calculate areal extend of permafrost    : " << permaOpt << endl;
  output_file << "23. New reservoir algorithm    : " << resOpt << endl;
  output_file << "24. Station correction status    : " << statcorrOpt << endl;
  output_file << "25. Additional groundwater recharge below surface water bodies in arid / semi-arid regions    : " << aridareaOpt << endl;
  output_file << "26. Modified routing of groundwater outflow and fast runoff    : " << fractionalRoutingOpt << endl;
  output_file << "27. Allow evaporation from rivers    : " << riverEvapoOpt << endl;
  output_file << "28. Take into account NUS of riparian cells    : " << aggrNUsGloLakResOpt << endl;
  output_file << "29. Correction of G_gammaHBV in the North China Plain    : " << gammaHBV_CorrOpt << endl;
  // FP20161018N003 Enable reading input data from climate land mask
  output_file << "30. Climate land mask status    : " << climate_spatial_resolution << endl;
  // FP20161018N002 Reservoir operation start years
  output_file << "31. Reservoir operation start year status    : " << resYearOpt << endl;
  output_file << "32. Reservoir operation reference year: " << resYearReference << endl;
  output_file << "33. Reservoir operation first year to use: " << resYearFirstToUse << endl;
  output_file << "34. Reservoir operation last year to use: " << resYearLastToUse << endl;
  output_file << "35. Anthropogenic / naturalized status    : " << antNatOpt << endl;
  output_file << "------------------------------------------------\n";

  output_file << "Values from file " << filename_out_opt << endl;
  output_file << "monthly binary files "<< endl;
  output_file << "G_PRECIPITATION_[YEAR].12.UNF0   : " << outPrec<< endl;
  output_file << "G_AET_[YEAR].12.UNF0   : " << outAET << endl;
  output_file << "G_CELL_AET_[YEAR].12.UNF0/G_LAND_AET_[YEAR].12.UNF0   : " << outCellAET<< endl;
  output_file << "G_CELL_RUNOFF_[YEAR].12.UNF0   : " << outCellRunoff<< endl;
  output_file << "G_POT_CELL_RUNOFF_[YEAR].UNF0   : " << outPotCellRunoffAnnual  << endl;
  output_file << "G_CELL_SURFACE_RUNOFF_[YEAR].12.UNF0-> cell surface runoff (incl. l/w)   : " << outCellSurface<< endl;
  output_file << "G_RUNOFF_[YEAR].12.UNF0   : " << outRunoff << endl;
  output_file << "G_URBAN_RUNOFF_[YEAR].12.UNF0   : " << outUrbanRunoff<< endl;
  output_file << "G_SURFACE_RUNOFF_[YEAR].12.UNF0   : " << outSurfaceRunoff<< endl;
  output_file << "G_GW_RUNOFF_[YEAR].12.UNF0   : " << outGWRunoff << endl;
  output_file << "G_GW_RECHARGE_[YEAR].12.UNF0   : " << outGWRecharge << endl;
  output_file << "all single water storages [km3]   : " << outSingleStorages << endl;
  output_file << "G_MIN_RIVER_AVAIL_[YEAR].12.UNF0/G_MAX_RIVER_AVAIL_[YEAR].12.UNF0   : " << outMinMaxRiverAvail<< endl;
  output_file << "G_SOIL_WATER_[YEAR].12.UNF0   : " << outSoilWater << endl;
  output_file << "G_LAI_[YEAR].12.UNF0   : " << outLAI<< endl;
  output_file << "G_ALBEDO_[YEAR].12.UNF0   : " << outAlbedo<< endl;
  output_file << "G_INTERCEPTION_[YEAR].12.UNF0   : " << outInterception << endl;
  output_file << "G_CANOPY_WATER_[YEAR].12.UNF0   : " << outCanopyWater << endl;
  output_file << "G_MAX_CANOPY_WATER_[YEAR].12.UNF0   : " << outmaxCanopyWater<< endl;
  output_file << "short wave rad components   : " << outNetShortWaveRadiation << endl;
  output_file << "long wave rad components   : " << outNetLongWaveRadiation << endl;
  output_file << "G_NET_RADIATION_[YEAR].12.UNF0   : " << outNetRadiation << endl;
  output_file << "G_OPEN_WATER_PET_[YEAR].12.UNF0   : " << outOpenWaterPET << endl;
  output_file << "G_OPEN_WATER_EVAP_[YEAR].12.UNF0   : " << outOpenWaterEvap << endl;
  output_file << "G_PET_[YEAR].12.UNF0   : " << outPET << endl;
  output_file << "G_TOTAL_PET_[YEAR].12.UNF0   : " << outTotalPET<< endl;
  output_file << "G_RIVER_AVAIL_[YEAR].12.UNF0   : " << outRiverAvail << endl;
  output_file << "G_RIVER_IN_UPSTREAM_[YEAR].12.UNF0   : " << outRiverInUpstream<< endl;
  output_file << "G_RIVER_VELOCITY_[YEAR].12.UNF0   : " << outRiverVelo << endl;
  output_file << "G_SNOW_COVER_FRAC[YEAR].12.UNF0   : " << outSnowCover << endl;
  output_file << "G_SNOW_WATER_EQUIV[YEAR].12.UNF0   : " << outSWE<< endl;
  output_file << "G_SNOW_FALL_[YEAR].12.UNF0   : " << outSnowFall << endl;
  output_file << "G_SNOW_MELT_[YEAR].12.UNF0   : " << outSnowMelt << endl;
  output_file << "G_SNOW_EVAP_[YEAR].12.UNF0   : " << outSnowEvap << endl;
  output_file << "G_SURFACE_WATER_STORAGE_[YEAR].12.UNF0   : " << outSurfStor<< endl;
  output_file << "G_LOC_LAKE_STOR_km3_[YEAR].12.UNF0   : " << outLocLakeStorage    << endl;
  output_file << "G_LOC_WETL_STOR_km3_[YEAR].12.UNF0   : " << outLocWetlStorage    << endl;
  output_file << "G_GLO_LAKE_STORAGE_km3_[YEAR].12.UNF0   : " << outGloLakeStorage    << endl;
  output_file << "G_GLO_WETL_STORAGE_km3_[YEAR].12.UNF0   : " << outGloWetlStorage    << endl;
  output_file << "G_RES_STORAGE_km3_[YEAR].12.UNF0   : " << outResStorage        << endl;
  output_file << "G_RIVER_STORAGE_km3_[YEAR].12.UNF0   : " << outRiverStorage      << endl;
  output_file << "G_TOTAL_STORAGES_km3[YEAR].12.UNF0   : " << outTotalWaterInStorages_km3<< endl;
  output_file << "G_TOTAL_STORAGES_mm_[YEAR].12.UNF0   : " << outTotalWaterInStorages_mm<< endl;
  output_file << "G_ACTUAL_NAS_[YEAR].12.UNF0   : " << outActualNAS << endl;
  output_file << "G_ACTUAL_NAG_[YEAR].12.UNF0   : " << outActualNAG << endl;
  output_file << "G_ALLOC_USE_IN_2NDCELL_[YEAR].12.UNF0   : " << outAllocUsein2ndCell<< endl;
  output_file << "G_SATIS_ALLOC_USE_IN_2NDCELL_[YEAR].12.UNF0   : " << outSatisAllocUsein2ndCell << endl;
  output_file << "G_CONS_WATER_USE_[YEAR].UNF0   : " << outConsWaterUse << endl;
  output_file << "G_UNSAT_USE_SW_[YEAR].UNF0   : " << outUnsatUseSW << endl;
  output_file << "G_UNSAT_USE_SW_PREV_[YEAR].UNF0   : " << outUnsatUseSWPrev << endl;
  output_file << "G_GWR_SURFACE_WATER_BODIES_mm_[YEAR].12.UNF0   : " << outGwrSwb    << endl;
  output_file << "G_FRACTION_SURFACE_WATER_BODIES_[YEAR].12.UNF0   : " << outFswb    << endl;
  output_file << "G_LAND_AREA_FRACTION_[YEAR].12.UNF0   : " << outLandAreaFrac    << endl;
  output_file << "G_TEMPERATURE_[YEAR].12.UNF0   : " << outTemperature    << endl;
  output_file << "G_SUNSHINE_km3_[YEAR].12.UNF0   : " << outSunshine    << endl;
  output_file << "G_ACTUAL_WATER_CONSUMPTION_[YEAR].12.UNF0   : " << outWCa    << endl;
  output_file << "G_ACTUAL_WATER_CONSUMPTION_INCELLUSED_[YEAR].12.UNF0   : " << outWCaInCellUsed    << endl;
  output_file << "G_CELLAET_CONSUSE_km3_[YEAR].12.UNF0   : " << outCellAETWCa    << endl;
  output_file << "G_LOC_WETL_EXTENT_km2_[YEAR].12.UNF0   : " << outLocWetlExt    << endl;
  output_file << "G_GLO_WETL_EXTENT_km2_[YEAR].12.UNF0   : " << outGloWetlExt    << endl;
  output_file << "G_TOTAL_GW_RECHARGE_[YEAR].12.UNF0   : " << outTotalGWR    << endl;
  output_file << "daily binary files "<< endl;
  output_file << "G_PRECIPITATION_[YEAR]*.UNF0   : " << outPrecDaily<< endl;
  output_file << "G_CELL_AET_[YEAR]*.UNF0   : " << outCellAETDaily << endl;
  output_file << "G_CELL_RUNOFF_[YEAR]*.UNF0   : " << outCellRunoffDaily << endl;
  output_file << "G_CELL_SURFACE_RUNOFF_[YEAR]*.UNF0   : " << outCellSurfaceDaily<< endl;
  output_file << "G_SURFACE_RUNOFF_[YEAR]*.UNF0   : " << outSurfaceRunoffDaily << endl;
  output_file << "G_GW_RUNOFF_[YEAR]*.UNF0   : " << outGWRunoffDaily << endl;
  output_file << "G_GW_RECHARGE_[YEAR]*.UNF0   : " << outGWRechargeDaily << endl;
  output_file << "G_GW_STORAGE_[YEAR]*.UNF0   : " << outGWStorageDaily << endl;
  output_file << "G_SOIL_WATER_[YEAR]*.UNF0   : " << outSoilWaterDaily << endl;
  output_file << "G_LAI_[YEAR]*.UNF0/G_KC_[YEAR]*.UNF0   : " << outLAIDaily<< endl;
  output_file << "G_ALBEDO_[YEAR]*.UNF0   : " << outAlbedoDaily<< endl;
  output_file << "G_INTERCEPTION_[YEAR]*.UNF0   : " << outInterceptionDaily << endl;
  output_file << "G_CANOPY_WATER_[YEAR]*.UNF0   : " << outCanopyWaterDaily<< endl;
  output_file << "G_MAX_CANOPY_WATER_[YEAR]*.UNF0   : " << outmaxCanopyWaterDaily<< endl;
  output_file << "G_PET_[YEAR]*.UNF0   : " << outPETDaily << endl;
  output_file << "G_TOTAL_PET_[YEAR]*.UNF0   : " << outTotalPETDaily << endl;
  output_file << "G_RIVER_AVAIL_[YEAR]*.UNF0   : " << outRiverAvailDaily << endl;
  output_file << "G_RIVER_VELOCITY_[YEAR]*.UNF0   : " << outRiverVeloDaily<< endl;
  output_file << "G_SNOW_COVER_FRAC[YEAR]*.UNF0   : " << outSnowCoverDaily << endl;
  output_file << "G_SNOW_WATER_EQUIV[YEAR]*.UNF0   : " << outSWEDaily << endl;
  output_file << "G_SNOW_FALL_[YEAR]*.UNF0   : " << outSnowFallDaily << endl;
  output_file << "G_SNOW_MELT_[YEAR]*.UNF0   : " << outSnowMeltDaily << endl;
  output_file << "G_SURFACE_WATER_STORAGE_[YEAR]*.UNF0   : " << outSurfStorDaily << endl;
  output_file << "all single daily water storages [km3]   : " << outSingleStoragesDaily << endl;
  output_file << "G_TOTAL_STORAGE_km3_[YEAR]*.UNF0   : " << outTotalWaterInStoragesDaily_km3 << endl;
  output_file << "G_TOTAL_STORAGE_mm_[YEAR]*.UNF0   : " << outTotalWaterInStoragesDaily_mm << endl;
  output_file << "G_TOTAL_STORAGE_STARTEND_km3_[YEAR].2.UNF0   : " << outTotalWaterInStoragesStartEndDaily_km3 << endl;
  output_file << "G_EXT_RAD_[YEAR]*.UNF0   : " << outExtRadiationDaily << endl;
  output_file << "G_SHORTWAVE_DOWN_RAD_[YEAR]*.UNF0   : " << outShortDownRadiationDaily << endl;
  output_file << "G_SHORTWAVE_UP_RAD_[YEAR]*.UNF0   : " << outShortUpRadiationDaily << endl;
  output_file << "G_NET_SHORTWAVE_RAD_[YEAR]*.UNF0   : " << outNetShortWaveRadiationDaily << endl;
  output_file << "G_LONGWAVE_DOWN_RAD_[YEAR]*.UNF0   : " << outLongDownRadiationDaily << endl;
  output_file << "G_LONGWAVE_UP_RAD_[YEAR]*.UNF0   : " << outLongUpRadiationDaily << endl;
  output_file << "G_NET_LONGWAVE_RAD_[YEAR]*.UNF0   : " << outNetLongWaveRadiationDaily << endl;
  output_file << "G_NET_RADIATION_[YEAR]*.UNF0   : " << outNetRadiationDaily << endl;
  output_file << "G_GWR_SURFACE_WATER_BODIES_mm_[YEAR]*.UNF0	: " << outGwrSwbDaily << endl;
  output_file << "G_FRACTION_SURFACE_WATER_BODIES_[YEAR]*.UNF0	: " << outFswbDaily << endl;
  output_file << "G_LAND_AREA_FRACTION_[YEAR]*.UNF0	: " << outLandAreaFracDaily << endl;
  output_file << "G_GWRUN_SURFRUN_[YEAR]*.UNF0	: " << outGwrunSurfrunDaily << endl;
  output_file << "G_CELLAET_CONSUSE_[YEAR]*.UNF0	: " << outCellAETWCaDaily << endl;
  output_file << "additional binary files "<< endl;
  output_file << "G_GW_FACTOR.UNF0   : " << outGWFactor << endl;
  output_file << "G_RG_max.UNF2   : " << outRGmax << endl;
  output_file << "G_Smax.UNF0   : " << outmaxSoilWater<< endl;
  output_file << "G_RootDepth.UNF0   : " << outRootingDepth<< endl;
  output_file << "G_KOEPPEN.UNF2   : " << outClcl<< endl;
  output_file << "G_LAI_MAX.UNF0/G_LAI_MIN.UNF0   : " << outLAImax<< endl;
  output_file << "G_AIRFROST_[YEAR].UNF0   : " << outAirFrost<< endl;
  output_file << "G_SURFACEFROST_[YEAR].UNF0   : " << outSurfaceFrost<< endl;
  output_file << "G_RH_[YEAR].UNF0   : " << outRH<< endl;
  output_file << "ASCII files"<< endl;
  output_file << "ALL_UPSTREAM_STATIONS.DAT, NO_ALL_UPSTREAM_STATIONS.DAT   : " << outAllUpStations << endl;
  output_file << "DIRECT_UPSTREAM_STATIONS.DAT, NO_DIRECT_UPSTREAM_STATIONS.DAT   : " << outDirectUpStations << endl;
  output_file << "RAIN_DAYS[XX].DAT (RAIN_DAYS28.DAT, RAIN_DAYS30.DAT, RAIN_DAYS31.DAT)   : " << outRainDays << endl;
  output_file << "STATION_DISCHARGE_ANNUAL.DAT   : " << outStationDischargeAnnual << endl;
  output_file << "STATION_DISCHARGE_MONTHLY.DAT   : " << outStationDischargeMonthly << endl;
  output_file << "STATION_LIST.OUT (not yet implemented -> will be written in any case!)   : " << outStationList << endl;
  output_file << "SUPERBASIN_CLIMATE.[YEAR].OUT  (2015-10 obsolete, not written) : " << outSuperbasinClimate << endl; // Currently obsolete, probably not working correctly in WaterGAP2.2b, do in postprocessing  // FP 2015
  output_file << "RESERVOIR_INFLOW_MONTHLY.OUT, RESERVOIR_OUTFLOW_MONTHLY.OUT   : " << outResvoirMonthly << endl;
  output_file << "additional ASCII files (option: save daily values)"<< endl;
  output_file << "DAILY_VALUES_[YEAR].XXX   : " << outDailyValues << endl;
  output_file << "DAILY_INTERCEPTION.XXX   : " << outDailyInterception << endl;
  output_file << "STATION_DISCHARGE_DAILY.OUT   : " << outStationDischargeDaily << endl;
  output_file << "STATION_VELOCITY_DAILY.OUT   : " << outStationVelocityDaily << endl;
  output_file << "LOC_LAKE_STORAGE.OUT   : " << outLocLakeStorageDaily<< endl;
  output_file << "GLO_LAKE_STORAGE.OUT   : " << outGloLakeStorageDaily<< endl;
  output_file << "LOC_WETL_STORAGE.OUT   : " << outLocWetlStorageDaily<< endl;
  output_file << "GLO_WETL_STORAGE.OUT   : " << outGloWetlStorageDaily<< endl;
  output_file << "RESERVOIR_STORAGE.OUT   : " << outResStorageDaily<< endl;
  output_file << "------------------------------------------------\n";
  output_file << "End of MODEL_SETTINGS.OUT\n";
}

void optionClass::closeModelSettingsFile()
{
	modelSettingsFile.close();
}
    // HMS MODEL_SETTINGS output end

   // HMS 2013-11-21 SINGLECELL output start
short optionClass::createSingleCellFileList(const char *inputDir, const char *outputDir)
{
  extern short G_toBeCalculated[ng];
  extern short G_SingleCelltoBeCalculated[ng];
  extern geoClass geo;
  extern routingClass routing;
  extern landClass land;
  char filename[250];

	// open file which contains IMAGE-nrs of single cells
	ifstream singlecells_file(filename_singlecells);

	if (!singlecells_file) {
		cerr << "Can not open file " << filename_singlecells << endl;
		exit(-1);
	}
	// open file for output
	sprintf(filename, "%s/SINGLECELL_LIST.OUT", outputDir);
	ofstream output_file(filename);

	if (!output_file) {
		cerr << "Can not open file " << filename << " for writing" << endl;
		exit(-1);
	}
	output_file << "# " << getTimeString() << filename << endl;
	output_file << "------------------------------------------------\n";
  short cellAlreadyInList;
  float longitude, latitude;
  char string[250];
	short singlecellnr = 0;
	const short maxLength = 25;
  short length;
	int k;


	extStack < int >stackSingleCellNum;
	int scn, scn2;
  Stack < char *>stackNamePointer;
	char *SingleCellName;


	do {	singlecells_file >> string;
    if (singlecells_file) {
			singlecells_file >> longitude >> latitude;

    length = strlen(string);
			if (length >= maxLength) {
				SingleCellName = new char[maxLength + 1];

				length = maxLength;
				strncpy(SingleCellName, string, length);
			} else {
				SingleCellName = new char[length + 1];

				strncpy(SingleCellName, string, length);
				SingleCellName[length] = '\0';
			}


			output_file << "Single cell name:       " << SingleCellName << endl;
			output_file << "Longitude/Latitude: " << longitude << " / " << latitude << endl;

			// determine cell number of station
			scn = geo.cellNumByLonLat(longitude, latitude);
			output_file << "IMAGE number:        " << scn << endl;
			if (scn <= 1)
				output_file << "Not part of the landmask! Skipping! " << endl;
			else {
				// is that cell number already in the list ?
				cellAlreadyInList = 0;
				for (k = 1; k <= stackSingleCellNum.getNumberOfElements(); k++) {
					stackSingleCellNum.getElement(k, scn2);
					if (scn == scn2) {
						cellAlreadyInList = 1;
						break;
					}
				}

				if (cellAlreadyInList) {
					output_file << "Cell is already defined at Nr.: "
						<< stackSingleCellNum.getNumberOfElements() - (k - 1)
						<< endl << "Skipping that cell..." << endl;
				} else {
					// everything OK, cell should be used

          if (2 == basin){ // HMS 2011-11-21 to calculate this cells anyway.
            G_toBeCalculated[scn-1] = 1;
            G_SingleCelltoBeCalculated[scn-1] = 1;
          }
					output_file << "Cell area according to position of cell (GAREA.UNF0): " << geo.areaOfCellByArrayPos(scn) << endl;
          output_file << "Land fraction according to GFREQ: " << (float)geo.G_landfreq_const[scn-1] << endl; // FP20161018N002
          output_file << "Water fraction according to GFREQW: " << (float)geo.G_fwaterfreq_const[scn-1] << endl; // FP20161018N002
            for (int n = 0; n < ng; ++n) {
              if (n == scn - 1) {
                output_file << "Fraction of global wetlands according to G_GLOWET: " << (float)routing.G_glo_wetland[n] << endl;
                output_file << "Fraction of local wetlands according to G_LOCWET: " << (float)routing.G_loc_wetland[n] << endl;
                output_file << "Fraction of global lakes according to G_GLOLAK: " << routing.G_glo_lake[n] << endl;
                output_file << "Fraction of lokal lakes according to G_LOCLAK: " << routing.G_loc_lake[n] << endl;
                output_file << "Area of lakes according to G_LAKAREA: " << routing.G_lake_area[n] << endl;
                output_file << "Area of reservoirs and regulated lakes according to G_RESAREA: " << routing.G_reservoir_area[n] << endl;
                output_file << "Land cover type: " << (float)land.G_landCover[n] << endl;
               }
            }

          stackSingleCellNum.push(scn);
					stackNamePointer.push(SingleCellName);
			    singlecellnr++;
			output_file << "Single Cell Output Nr.: " << singlecellnr << endl;
			output_file << "------------------------------------------------\n";
	        }
	     }
      }
    } while (singlecells_file);
	singlecells_file.close();

   cellname = new char *[stackNamePointer.getNumberOfElements()];

	for (int i = stackNamePointer.getNumberOfElements(); i >= 1; i--)
		if (stackNamePointer.pop(SingleCellName))
			cellname[i - 1] = SingleCellName;

	cellnr = new int[stackSingleCellNum.getNumberOfElements()];
	for (int i = stackSingleCellNum.getNumberOfElements(); i >= 1; i--)
		if (stackSingleCellNum.pop(scn))
			cellnr[i - 1] = scn;

		//	stackNamePointer.pop(SingleCellName);

	numberOfSingleCells = singlecellnr;
return numberOfSingleCells;
}


optionClass::~optionClass()
{
	delete[] filename_data; 	filename_data 		= NULL;
	delete[] filename_time; 	filename_time 		= NULL;
	delete[] filename_opt; 		filename_opt  		= NULL;
	delete[] filename_out_opt; 	filename_out_opt	= NULL;
	delete[] filename_stations; 	filename_stations	= NULL;
	delete[] filename_singlecells; 	filename_singlecells	= NULL;
        delete[] filename_scout_opt; 	filename_scout_opt	= NULL;

	if (landCoverOpt) {
		delete[] landCoverYears;      landCoverYears      = NULL;
		delete[] landCoverYearsStart; landCoverYearsStart = NULL;
		delete[] landCoverYearsEnd;   landCoverYearsEnd   = NULL;
	}
}

