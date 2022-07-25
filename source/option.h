
/***********************************************************************
*
* see former changes at file option.h.versioninfos.txt
*
***********************************************************************/
#if !defined (_option_h_)
#define _option_h_
#include <fstream>          // MODEL_SETTINGS output
using namespace std;        // MODEL_SETTINGS output

class optionClass {
  public:
	char input_dir[250];
	char output_dir[250];
	char climate_dir[250];
	char routing_dir[250];
	char water_use_dir[250];
	char land_cover_dir[250];
	char *filename_stations;
	char *filename_singlecells; // HMS 2013-11-21 	
        char *filename_scout_opt; // HMS 2013-11-21
	char **cellname;  // HMS 2013-11-21
	
    int *cellnr;  // HMS 2013-11-21  
	short start_year;
	short end_year;
	short evalStartYear;
	short init_years;
	short time_step;
	short mid_year;
	short landCoverYearsInList; // vor variable land cover map : number of maps
	short *landCoverYears;      // vor variable land cover map : list of years with map
	short *landCoverYearsStart; // vor variable land cover map : begin of period with one map
	short *landCoverYearsEnd;   // vor variable land cover map : end of period with one map

	short fileEndianType;
	short basin;
	short grid_store;
	short grid_store_TypeForStorages;
	short day_store;
	short time_series;
	short cloud;
	short raindays;
	short prec_correct;
	short intercept;
	short calc_albedo;
	short use_kc;
	short landCoverOpt;
	short rout_prepare;
	short timeStepCheckFlag;
	short riverveloOpt;
	short subtract_use;
	short use_alloc;
	short delayedUseSatisfaction;
	short petOpt;
	short clclOpt;
	short permaOpt;
    short resOpt;
    short statcorrOpt;
    short numberOfSingleCells; // HMS 2013-11-21
    short aridareaOpt;
    short fractionalRoutingOpt;
	short riverEvapoOpt;
	short aggrNUsGloLakResOpt;
	short gammaHBV_CorrOpt;

    // FP20161018N003 Enable reading input data from climate land mask
    // Extended grids for climate land mask (union of Standard/SLM and WATCH CRU Land Mask/WLM)
    short climate_spatial_resolution;

    // FP20161018N002 Reservoir operation start years
    // Treat anthropogenic vs. naturalized runs
    short antNatOpt;
    // Settings for how to consider reservoir and regulated lakes
    short resYearOpt;
    short resYearReference;
    short resYearFirstToUse;
    short resYearLastToUse;
    // FP20161018N004 Net use / net abstraction for reservoir algorithm: new filename
    // Code for reference period G_NUs_yrfirst4digits_yrlast4digits.UNF, e.g. G_NUs_1971_2000.UNF
    // Now: name with averaging period of full year start and end (obsolete: (YearStart 2 digits + YearEnd 2 digits)
    short resNUsMeanYearFirst;
    short resNUsMeanYearLast;

    // options from ISISMIP2.2 to integrate/change
    bool outCellAETActualUse; // postprocessing
    bool outLandWaterFractions; // replace outLandAreaFract  -> outFractions
    bool outLandStorageChange; // delete
    bool outResCapacity; //new
    bool outRunoffDaily;                             //6
// = outSurfaceRunoff + outGWRunoff;
    bool outAETDaily;
    //= outAET

	// bool variables for output options
	// monthly binary files
        bool outPrec;                               //0
        bool outAET;                                //1
        bool outCellAET;                            //2
        bool outCellRunoff;                         //3
        bool outPotCellRunoffAnnual;                //4
        bool outCellSurface;                        //5
        bool outRunoff;                             //6
        bool outUrbanRunoff;                        //7
        bool outSurfaceRunoff;                      //8
        bool outGWRunoff;                           //9
        bool outGWRecharge;                         //10
        bool outSingleStorages;                     //11
        bool outMinMaxRiverAvail;                   //12
        bool outSoilWater;                          //13
        bool outLAI;                                //14
        bool outAlbedo;                             //15
        bool outInterception;                       //16
        bool outCanopyWater;                        //17
        bool outmaxCanopyWater;                     //18
        bool outRiverPET;                           //19
        bool outNetShortWaveRadiation;              //20
        bool outNetLongWaveRadiation;               //21
        bool outNetRadiation;                       //22
        bool outOpenWaterPET;                       //23
        bool outOpenWaterEvap;                      //24
        bool outPET;                                //25
        bool outTotalPET;                           //26
        bool outRiverAvail;                         //27
        bool outRiverInUpstream;                    //28
        bool outRiverVelo;                          //29
        bool outSnowCover;                          //30
        bool outSWE;                                //31
        bool outSnowFall;                           //32
        bool outSnowMelt;                           //33
        bool outSnowEvap;                           //34
        bool outSurfStor;                           //35
        bool outLocLakeStorage;                     //36
        bool outLocWetlStorage;                     //37
        bool outGloLakeStorage;                     //38
        bool outGloWetlStorage;                     //39
        bool outResStorage;                         //40
        bool outRiverStorage;                       //41
        bool outTotalWaterInStorages_km3;           //42
        bool outTotalWaterInStorages_mm;            //43
        bool outActualNAS;                          //44
    	bool outActualNAG;                          //45
		bool outWCa;                                //46
		bool outWCaInCellUsed;       				//47
        bool outAllocUsein2ndCell;                  //48
    	bool outSatisAllocUsein2ndCell;             //49
        bool outConsWaterUse;                       //50
        bool outUnsatUseSW;                         //51
        bool outUnsatUseSWPrev;                     //52
        bool outGwrSwb;                             //53
        bool outFswb;                               //54
        bool outLandAreaFrac;                       //55
        bool outTemperature;                        //56
        bool outSunshine;                           //57
        bool outCellAETWCa;                         //58
        bool outLocWetlExt;                         //59
        bool outGloWetlExt;                         //60
        bool outTotalGWR;                           //61

        // daily binary files
        bool outPrecDaily;                          //62
        bool outCellAETDaily;                       //63
        bool outCellRunoffDaily;                    //64
        bool outCellSurfaceDaily;                   //65
        bool outSurfaceRunoffDaily;                 //66
        bool outGWRunoffDaily;                      //67
        bool outGWRechargeDaily;                    //68
        bool outGWStorageDaily;                     //69
        bool outSoilWaterDaily;                     //70
        bool outLAIDaily;                           //71
        bool outAlbedoDaily;                        //72
        bool outInterceptionDaily;                  //73
        bool outCanopyWaterDaily;                   //74
        bool outmaxCanopyWaterDaily;                //75
        bool outExtRadiationDaily;                  //76
        bool outShortDownRadiationDaily;            //77
        bool outShortUpRadiationDaily;              //78
        bool outNetShortWaveRadiationDaily;         //79
        bool outLongDownRadiationDaily;             //80
        bool outLongUpRadiationDaily;               //81
        bool outNetLongWaveRadiationDaily;          //82
        bool outNetRadiationDaily;                  //83
        bool outPETDaily;                           //84
        bool outTotalPETDaily;                      //85
        bool outRiverAvailDaily;                    //86
        bool outRiverVeloDaily;                     //87
        bool outSnowCoverDaily;                     //88
        bool outSWEDaily;                           //89
        bool outSnowFallDaily;                      //90
        bool outSnowMeltDaily;                      //91
        bool outSurfStorDaily;                      //92
    bool outSingleStoragesDaily;                    //93
    bool outTotalWaterInStoragesDaily_km3;          //94
    bool outTotalWaterInStoragesDaily_mm;           //95
    bool outTotalWaterInStoragesStartEndDaily_km3;   //96
        bool outGwrSwbDaily;                        //97
        bool outFswbDaily;                          //98
        bool outLandAreaFracDaily;                  //99
        bool outGwrunSurfrunDaily;                  //100
        bool outCellAETWCaDaily;                         //101

        // additional binary files
        bool outGWFactor;                           //102
        bool outRGmax;                              //103
        bool outmaxSoilWater;                       //104
        bool outRootingDepth;                       //105
        bool outClcl;                               //106
        bool outLAImax;                             //107
        bool outAirFrost;                           //108
        bool outSurfaceFrost;                       //109
        bool outRH;                                 //110

	// ASCII files
        bool outAllUpStations;                      //111
        bool outDirectUpStations;                   //112
        bool outRainDays;                           //113
        bool outStationDischargeAnnual;             //114
        bool outStationDischargeMonthly;            //115
        bool outStationList; 		// not implemented yet!!116
        bool outSuperbasinClimate;                  //117
        bool outResvoirMonthly;                     //118

        // additional ASCII files (option: save daily values)
        bool outDailyValues;                        //119
        bool outDailyInterception;                  //120
        bool outStationDischargeDaily;              //121
        bool outStationVelocityDaily;               //122
        bool outLocLakeStorageDaily;                //123
        bool outGloLakeStorageDaily;                //124
        bool outLocWetlStorageDaily;                //125
        bool outGloWetlStorageDaily;                //126
        bool outResStorageDaily;                    //127


// single cell output options
        bool scoutTemp;             //sc0
        bool scoutExtRad;           //sc1
        bool scoutShortDownRad;     //sc2
        bool scoutAlbedo;           //sc3
        bool scoutShortUpRad;       //sc4
        bool scoutNetShortRad;      //sc5
        bool scoutLongDownRad;      //sc6
        bool scoutLongUpRad;        //sc7
        bool scoutNetLongRad;       //sc8
        bool scoutNetRad;           //sc9
        bool scoutLAI;              //sc10
        bool scoutKc;               //sc11
        bool scoutLandPET;          //sc12
        bool scoutCellPET;          //sc13
        bool scoutPrecip;           //sc14
        bool scoutcPrecip;          //sc15
        bool scoutCanopyWater;      //sc16
        bool scoutmaxCanopyWater;   //sc17
        bool scoutInterception;     //sc18
        bool scoutSnowfall;         //sc19
        bool scoutSnowWater;        //sc20
        bool scoutSnowCov;          //sc21
        bool scoutSnowmelt;         //sc22
        bool scoutSoilWater;        //sc23
        bool scoutSurfaceRunoff;    //sc24
        bool scoutGwRunoff;         //sc25
        bool scoutGwRecharge;       //sc26
        bool scoutCellAET;          //sc27
        bool scoutCellRunoffkm3;    //sc28
        bool scoutCellRunoffmm;     //sc29
        bool scoutCellSRunoff;      //sc30
        bool scoutQ;                //sc31
        bool scoutFlowVelo;         //sc32
        bool scoutLocLake;          //sc33
        bool scoutLocWet;           //sc34
        bool scoutGloLake;          //sc35
        bool scoutReservoir;        //sc36
        bool scoutGloWet;           //sc37
        bool scoutRiver;            //sc38
        bool scoutSurfaceStor;      //sc39
        bool scoutTWSkm3;           //sc40
        bool scoutTWSmm;            //sc41
        bool scoutGwStor;           //sc42
        bool scoutCanopyStor;       //sc43
        bool scoutSnowStor;         //sc44
        bool scoutSoilStor;         //sc45
        bool scoutGwrSwb;           //sc46
        bool scoutFswb;             //sc47
        bool scoutLandAreaFrac;     //sc48

	void init(int optionc, char* optionv[]);

	~optionClass();
	
    void createModelSettingsFile(const char *outputDir);     // MODEL_SETTINGS output
    void closeModelSettingsFile();                          // MODEL_SETTINGS output
    short createSingleCellFileList(const char *inputDir, const char *outputDir); // HMS 2013-11-21 single cell output
	
  private:

	char *filename_data;
	char *filename_time;
	char *filename_opt;
	char *filename_out_opt;

  ofstream modelSettingsFile;    // MODEL_SETTINGS output
};
#endif
