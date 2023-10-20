
/***********************************************************************
*
* see former changes at file option.h.versioninfos.txt
*
***********************************************************************/
#if !defined (_option_h_)
#define _option_h_
#include <fstream>          // MODEL_SETTINGS output
#include <string>
#include <vector>
#include "configFile.h"

using namespace std;        // MODEL_SETTINGS output

class optionClass {
  public:
	string input_dir;
	string output_dir;
	string climate_dir;
	string glacier_dir;
	string routing_dir;
	string water_use_dir;
	string land_cover_dir;

	string filename_stations;
	string filename_singlecells;
    string filename_scout_opt;
	char **cellname;
	
    int *cellnr;
	short start_year;
	short end_year;
	short evalStartYear;
	short init_years;
	short time_step;
	short landCoverYearsInList; // vor variable land cover map : number of maps
	vector<short> landCoverYears;      // vor variable land cover map : list of years with map
	vector<short> landCoverYearsStart; // vor variable land cover map : begin of period with one map
	vector<short> landCoverYearsEnd;   // vor variable land cover map : end of period with one map

	short fileEndianType;
	short basin;
	short grid_store;
	short grid_store_TypeForStorages;
	short day_store;
	short time_series;
	short cloud;
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
    short numberOfSingleCells;
    short aridareaOpt;
    short fractionalRoutingOpt;
	short riverEvapoOpt;
	short aggrNUsGloLakResOpt;

    // reading input data from climate land mask
    // Extended grids for climate land mask (union of Standard/SLM and WATCH CRU Land Mask/WLM)
    short climate_spatial_resolution;

    // Treat anthropogenic vs. naturalized runs
    short antNatOpt;
    // Settings for how to consider reservoir and regulated lakes
    short resYearOpt;
    short resYearReference;
    short resYearFirstToUse;
    short resYearLastToUse;

    // Settings for water temperature calculation
    short calc_wtemp;

    // Net use / net abstraction for reservoir algorithm: new filename
    // Code for reference period G_NUs_yrfirst4digits_yrlast4digits.UNF, e.g. G_NUs_1971_2000.UNF
    // Now: name with averaging period of full year start and end (obsolete: (YearStart 2 digits + YearEnd 2 digits)
    short resNUsMeanYearFirst;
    short resNUsMeanYearLast;
	short glacierOpt;

    // options from ISISMIP2.2 to integrate/change
    bool outCellAETActualUse; // postprocessing
    bool outLandWaterFractions; // replace outLandAreaFract  -> outFractions
    bool outLandStorageChange; // delete
    bool outResCapacity; //new
    bool outRunoffDaily;                             //6
    bool outAETDaily;

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
        bool outGwrunSurfrun;   					//62
		bool outGlacArea;                           //63
        bool outGlacAreaFrac;                       //64
        bool outGlacPrecip;                         //65
        bool outGlacRunoff;                         //66
        bool outGlacierStorage;                     //67

        // daily binary files
        bool outPrecDaily;                          //68
        bool outCellAETDaily;                       //69
        bool outCellRunoffDaily;                    //70
        bool outCellSurfaceDaily;                   //71
        bool outSurfaceRunoffDaily;                 //72
        bool outGWRunoffDaily;                      //73
        bool outGWRechargeDaily;                    //74
        bool outGWStorageDaily;                     //75
        bool outSoilWaterDaily;                     //76
        bool outLAIDaily;                           //77
        bool outAlbedoDaily;                        //78
        bool outInterceptionDaily;                  //79
        bool outCanopyWaterDaily;                   //80
        bool outmaxCanopyWaterDaily;                //81
        bool outExtRadiationDaily;                  //82
        bool outShortDownRadiationDaily;            //83
        bool outShortUpRadiationDaily;              //84
        bool outNetShortWaveRadiationDaily;         //85
        bool outLongDownRadiationDaily;             //86
        bool outLongUpRadiationDaily;               //87
        bool outNetLongWaveRadiationDaily;          //88
        bool outNetRadiationDaily;                  //89
        bool outPETDaily;                           //90
        bool outTotalPETDaily;                      //91
        bool outRiverAvailDaily;                    //92
        bool outRiverVeloDaily;                     //93
        bool outSnowCoverDaily;                     //94
        bool outSWEDaily;                           //95
        bool outSnowFallDaily;                      //96
        bool outSnowMeltDaily;                      //97
        bool outSurfStorDaily;                      //98
    bool outSingleStoragesDaily;                    //99
    bool outTotalWaterInStoragesDaily_km3;          //100
    bool outTotalWaterInStoragesDaily_mm;           //101
    bool outTotalWaterInStoragesStartEndDaily_km3;   //102
        bool outGwrSwbDaily;                        //103
        bool outFswbDaily;                          //104
        bool outLandAreaFracDaily;                  //105
        bool outGwrunSurfrunDaily;                  //106
        bool outCellAETWCaDaily;                    //107
        bool outGlacierStorageDaily;				//108
        bool outGlacierStorageDaily_mm;             //109
        bool outRedTemp;                            //110
        bool outTemp;                               //111

        // additional binary files
        bool outGWFactor;                           //112
        bool outRGmax;                              //113
        bool outmaxSoilWater;                       //114
        bool outRootingDepth;                       //115
        bool outClcl;                               //116
        bool outLAImax;                             //117
        bool outAirFrost;                           //118
        bool outSurfaceFrost;                       //119
        bool outRH;                                 //120

	// ASCII files
        bool outAllUpStations;                      //121
        bool outDirectUpStations;                   //122
        bool outRainDays;                           //123
        bool outStationDischargeAnnual;             //124
        bool outStationDischargeMonthly;            //125
        bool outStationList; 		// not implemented yet!!126
        bool outSuperbasinClimate;                  //127
        bool outResvoirMonthly;                     //128

        // additional ASCII files (option: save daily values)
        bool outDailyValues;                        //129
        bool outDailyInterception;                  //130
        bool outStationDischargeDaily;              //131
        bool outStationVelocityDaily;               //132
        bool outLocLakeStorageDaily;                //133
        bool outGloLakeStorageDaily;                //134
        bool outLocWetlStorageDaily;                //135
        bool outGloWetlStorageDaily;                //136
        bool outResStorageDaily;                    //137


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

    // output options for water temperature
    bool outWaterTempDaily;			 //138
    bool outWaterTempMonthly;        //139
    bool outWaterTempDailyAllSWB;    //140
    bool outWaterTempMonthlyAllSWB;	 //141

    void init(ConfigFile configFile);

	~optionClass();
	
    void createModelSettingsFile(const string outputDir);     // MODEL_SETTINGS output
    void closeModelSettingsFile();                          // MODEL_SETTINGS output
    short createSingleCellFileList(const string inputDir, const string outputDir); // HMS 2013-11-21 single cell output
	
  private:

	string filename_data;
	string filename_time;
	string filename_opt;
	string filename_out_opt;

    ofstream modelSettingsFile;    // MODEL_SETTINGS output
};
#endif
