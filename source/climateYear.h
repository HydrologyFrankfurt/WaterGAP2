#include "def.h"

class climateYearClass {
  public:
	climateYearClass(void);
	~climateYearClass(void);

	// store daily year values
        float (*G_temperature_d365)[365];
        float (*G_precipitation_d365)[365];
        float (*G_shortwave_d365)[365];
        float (*G_longwave_d365)[365];

  void read_climate_data_daily_per_year(short actual_year, short data_d);
  void split_climate_data_from_year_to_daily(short dayOfTheYear, short numberDaysInMonth);
	
	void init();
	void cleanup();
};
