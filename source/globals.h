// OE 01-2018
// OE 10-2020: adjust to aems_dev_2020
//#include "gridio.h"
//#include "grid_io_adapters.h"
//#include "grid.h"
#include "daily.h"
#include "rout_prepare.h"
#include "calib_basins.h"
#include "upstream_stations.h"
#include "routing.h"
#include "option.h"
#include "gw_frac.h"
#include "s_max.h"
#include "lai.h"
#include "timestring.h"
#include "calibration.h"
#include "climate.h"
#include "climateYear.h"
#include "glacierYear.h"
#include "geo.h"
#include "land.h"
#include "permafrost.h"
#include "clcl.h"
// alle headers von den hier definierten Klassen muessen eingebunden werden

#if !defined (_globals_h_)
#define _globals_h_  
// input grids

//---------ARID GROUNDWATER recharge calculations--------
extern Grid<short> G_aindex;
extern Grid<char> G_LDD;
extern Grid<short> G_toBeCalculated; // 0: no calculation
extern Grid<short> G_SingleCelltoBeCalculated;
extern Grid<signed short> G_sbasin; // superbasins
//extern float G_float[ng];
extern Grid<float> G_float;
//extern double timediff(struct timeval tv2_local, struct timeval tv1_local);


extern cbasinClass cbasin;
extern optionClass options;
extern laiClass lai;
extern routingClass routing;
extern calibGammaClass calibGamma;
extern climateClass climate;
extern climateYearClass climateYear;
extern glacierYearClass glacierYear;
extern upstreamStationClass directUpstSt;
extern upstreamStationClass allUpstSt;
extern dailyWaterBalanceClass dailyWaterBalance;
extern geoClass geo;
extern landClass land;
extern soilWatCapClass maxSoilWaterCap;
extern groundwaterFactorClass GW;
extern clclClass clcl;
//}

#endif


