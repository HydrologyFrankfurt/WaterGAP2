// OE 01-2018 (based on lines 104-134 in watergap.cpp)
#include "def.h"
#include "globals.h"


//---------ARID GROUNDWATER recharge calculations--------
Grid<short> G_aindex;
Grid<char> G_LDD;

Grid<short> G_toBeCalculated; // 0: no calculation
Grid<short> G_SingleCelltoBeCalculated; // HMS single cell output 2013-11-21

Grid<signed short> G_sbasin; // superbasins
 //float G_float[ng];
Grid<float> G_float;
// double timediff(struct timeval tv2_local, struct timeval tv1_local);

 cbasinClass cbasin;
 optionClass options;
 laiClass lai;
 routingClass routing;
 calibGammaClass calibGamma;
 climateClass climate;
 climateYearClass climateYear;
 glacierYearClass glacierYear;
 upstreamStationClass directUpstSt;
 upstreamStationClass allUpstSt;
 dailyWaterBalanceClass dailyWaterBalance;
 geoClass geo;
 landClass land;
 soilWatCapClass maxSoilWaterCap;
 groundwaterFactorClass GW;
 clclClass clcl;




