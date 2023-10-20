#if !defined (_climate_h_)
#define _climate_h_

#include "def.h"
#include "grid.h"
#include <string>

class climateClass {
  public:
		climateClass(void);
		~climateClass(void);

		// store daily values
		Daily31Grid<false, float> G_temperature_d;   // [oC]
		Daily31Grid<true,float> G_tmin_d;          // [oC]
		Daily31Grid<true,float> G_tmax_d;          // [oC]
		Daily31Grid<false, float> G_precipitation_d; // [mm]

		// if daily values are given use these parameters
		Daily31Grid<false, float> G_shortwave_d; // [W/m2]
		Daily31Grid<false, float> G_longwave_d;  // [W/m2]

		// memory only allocated if grids are used
		// for daily values
		Daily31Grid<false, unsigned char> G_windspeed_d;  // [1/10 m/s]
		Daily31Grid<false, unsigned short> G_vappres_d;    // [1/10 kPa]
		Daily31Grid<false, unsigned short> G_temprange_d;  // [oC] only for Hargreaves

		// store monthly values
		MonthlyGrid<false, short> G_temperature;     // [1/100 oC]
		MonthlyGrid<false, short> G_precipitation;   // [mm]
		MonthlyGrid<false, unsigned char> G_raindays;// [-]

		// the following array should be dynamic, because each simulation
		// only needs one of them!
		MonthlyGrid<false, short> G_shortwave;   // [W/m2]
		MonthlyGrid<false, short> G_longwave;    // [W/m2]
		MonthlyGrid<false, short> G_sunshine;    // [W/m2]
		MonthlyGrid<false, char> G_cloud;        // [per cent]

		// for monthly values
		MonthlyGrid<false, unsigned char> G_windspeed;     // [1/10 m/s]

		// for precipitation correction
		MonthlyGrid<false,float> G_prec_correct;

    // for new precipitation correction   //Adam Lettenmaier precipcorr start
		MonthlyGrid<false,float> Corr_factor_new;

    // for new precipitation correction
		MonthlyGrid<true,float> G_temp_mean;
		MonthlyGrid<true,float> R_snow_mean;
		MonthlyGrid<true,float> R_snow;
		MonthlyGrid<true,float> G_temperature_float;      //Adam Lettenmaier precipcorr end

        // Groundwater Temperature (yearly mean of AirTempC)
        Grid<double>G_GW_Temperature_y;

    // for annual Temp reduction factor (temp_diff)

        Grid<double> G_temp_diff; //[oC]

		// read daily or monthly climate data
		void read_climate_data_daily(short month, short actual_year);
		void read_climate_data_daily_withConvert(short month, short actual_year);

		void read_climate_longtermAvg();
		void read_climate_longtermAvg_withConvert();
		void read_annual_temp_reduction(short actual_year);

		void readPrecAdjustFactors(const std::string input_dir, const std::string climate_dir);   //Adam Lettenmaier precipcorr
		void adjustPrecipitation();
		void init();
		void initHG();
		void initPM();
        void initMD();
		void cleanup();

  private:
	 bool precCorrGridIsRead;

};
#endif
