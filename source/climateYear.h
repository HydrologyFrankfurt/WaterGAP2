#if !defined (_climateYear_h_)
#define _climateYear_h_

#include "def.h"

class climateYearClass {
  public:
	climateYearClass(void);
	~climateYearClass(void);

	// store daily year values
  DailyGrid<true,float> G_temperature_d365;
  DailyGrid<true,float> G_precipitation_d365;
  DailyGrid<true,float> G_shortwave_d365;
  DailyGrid<true,float> G_longwave_d365;

  void read_climate_data_daily_per_year(short actual_year);
  void split_climate_data_from_year_to_daily(short dayOfTheYear, short numberDaysInMonth);
	
	void init();
	void cleanup();
};
#endif
