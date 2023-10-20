#include "def.h"
#include "grid.h"
#include <string>

/*

glacierYear.cpp and .h introduced by Hannes Mueller Schmied based on climateYear.cpp and .h by
Heike Hoffmann-Dobrev (HHD) to read in .365 HYOGA2 glacier input. The original files were modified
by Denise Caceres to read in .365 and annual glacier input from the global glacier model described
in Marzeion et al. (2012): Past and future sea-level change from the surface mass balance of glaciers.

Usage of the glacierYearClass class:
This class is used whenever option glacierOpt (take glaciers into account) is set to 1.

 */

/**
 * @brief A class for reading glacier-related data.
 *
 * This class enables reading gridded (WATCH-CRU land-sea mask) time series of annual glacier area,
 * daily precipitation on glacier area and daily glacier mass change.
 *
 */

class glacierYearClass {
    public:
        glacierYearClass(void);
        ~glacierYearClass(void);

        DailyGrid<true, float> G_glacier_mass_change_d365;
        DailyGrid<true, float> G_precipitation_on_glacier_d365;
        Grid<> G_glacier_area;

        /**
        * @brief Read initial annual glacier area
        *
        * @param start_year the first year of the simulation period
        */
        void read_initial_glacier_area(short start_year);

        /**
        * @brief Read daily (.365) glacier data
        *
        * @param actual_year current simulation year
        */
        void read_glacier_data_daily_per_year(short actual_year);

        /**
        * @brief Read annual glacier data
        *
        * @param actual_year current simulation year
        */
        void read_glacier_data_yearly_per_year(short actual_year);
	  
        void init();
        void cleanup();
};
