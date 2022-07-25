#include "def.h"

class climateClass {
  public:
	climateClass(void);
	~climateClass(void);

  // store daily values
	float (*G_temperature_d)[31];   // [oC]
	float (*G_tmin_d)[31];          // [oC]
	float (*G_tmax_d)[31];          // [oC]
	float (*G_precipitation_d)[31]; // [mm]

  // if daily values are given use these parameters
	float (*G_shortwave_d)[31]; // [W/m2]
	float (*G_longwave_d)[31];  // [W/m2]
	
	// memory only allocated if grids are used
	// for daily values
	unsigned char  (*G_windspeed_d)[31];  // [1/10 m/s]
	unsigned short (*G_vappres_d)[31];    // [1/10 kPa]
	unsigned short (*G_temprange_d)[31];  // [oC] only for Hargreaves


  // store monthly values
	short (*G_temperature)[12];     // [1/100 oC]
	short (*G_tmin)[12];            // [1/100 oC]
	short (*G_tmax)[12];            // [1/100 oC]
	short (*G_precipitation)[12];   // [mm]
	unsigned char (*G_raindays)[12];// [-]

	// the following array should be dynamic, because each simulation
	// only needs one of them!
	short (*G_shortwave)[12];   // [W/m2]
	short (*G_longwave)[12];    // [W/m2]
	short (*G_sunshine)[12];    // [W/m2]
	char (*G_cloud)[12];        // [per cent]
	// for monthly values
	unsigned char (*G_windspeed)[12];     // [1/10 m/s]
	unsigned short (*G_vappres)[12];      // [1/10 kPa]
	unsigned short (*G_temprange)[12];    // [1/100 oC] only for Hargreaves
	
	// for precipitation correction
	float (*G_prec_correct)[12];
        // for new precipitation correction   //Adam Lettenmaier precipcorr start
	float Corr_factor_new[ng][12];    
        // for new precipitation correction
	float (*G_temp_mean)[12];
	float (*R_snow_mean)[12];
	float (*R_snow)[12];
        float (*G_temperature_float)[12];      //Adam Lettenmaier precipcorr end
     

	void createTempAvg(const char *climate_dir, const short start_year, const short end_year);
	void createPrecAvg(const char *climate_dir, const short start_year, const short end_year);
	void createRaindaysAvg(const char *climate_dir, const short start_year, const short end_year);
	void createRadiationAvg(const char *climate_dir, const short start_year, const short end_year);
	void createSunshineAvg(const char *climate_dir, const short start_year, const short end_year);
	void createCloudAvg(const char *climate_dir, const short start_year, const short end_year);
	// read daily or monthly climate data
	void read_climate_data_monthly(short actual_year, short data_m_d);
	void read_climate_data_daily(short month, short actual_year, short data_m_d);

        // FP20161018N003 Enable reading input data from climate land mask
        void read_climate_data_monthly_withConvert(short actual_year, short data_m_d);
        void read_climate_data_daily_withConvert(short month, short actual_year, short data_m_d);

        void read_climate_longtermAvg();
        void read_climate_longtermAvg_withConvert();


        void readPrecAdjustFactors(const char *input_dir, const char *climate_dir);   //Adam Lettenmaier precipcorr
	void adjustPrecipitation();
	void initPM();
	void initHG();
	void init();
	void cleanup();

  private:
	 bool precCorrGridIsRead;

};
