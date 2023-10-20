/**
* @brief calculation of water temperatures
*
* calculation of water temperature in the same order as WaterGAP in routing.cpp
* in contrast to WaterGAP (as of May 2020) the water temperature calculation of river considers the evaporation
*
* implementation of approach from:
* - Beek et al. (2012) doi:10.1029/2012WR011819
* - Wanders et al. (2019) doi:10.1029/2018WR023250
*
* @author Sebastian Ackermann
*/


#include <cmath>
#include <algorithm>

#include "calcWaterTemp.h"


using namespace std;

//constants
const double convK = 273.15;					// conversion [degC] <-> [K]
const double HeatCapW = 4180.;					// heat capacity water [J/(kg*K)]
const double rhoW = 1000.*1.e9;					// density water [kg/km3]
const int K_h = 20.;                            // turbulent heat exchange coefficient (water to air) [J/(s*m2*K] see Beek et al. (2012)
const int K_hIce = 8.;                          // turbulent heat exchange coefficient (water to Ice) [J/(s*m2*K] see Beek et al. (2012)
const double precipTempCorr = 1.5;              // temperature correction of precipitation
const double lamda_f = 333.4*1.e3;              // latent heat of fusion of ice [J/kg] see Beek et al. (2012)
const double albedoW = 0.08;                    // albedo water 0.08
const double albedoI = 0.6;                     // albedo ice 0.6

/**
* @brief Calculate water temperatures of surface water bodies (swb)
*
* Calculate water temperatures of local lakes, local wetlands, global lakes, reservoirs, global wetlands
* and river (cell outflow)
*
* @author Sebastian Ackermann
*
* @param const transfer &t: contains all variables needed from routing
*
* @return vector containing temperatures of surface body waters and ice thickness
* temperatures [degC]: [0] outflow/river, [1] locLake, [2] locWetland, [3] gloLake, [4] reservoir, [5] gloWetland if non existent in cell -9999. as temp
* ice thickness [km]: [6] river, [7] locLake, [8] locWetland, [9] gloLake, [10] reservoir, [11] gloWetland
*/
vector<double> calcWaterTempClass::calcWTemp(const transfer &t)
{
    //variables
    vector<double> calcReturn (12, -9999.);    // water temperature: [0] outflow, [1] locLake, [2] locWetland, [3] gloLake, [4] reservoir, [5] gloWetland; if non existent in cell -9999. as temp
                                                        // ice thickness: [6] river, [7] locLake, [8] locWetland, [9] gloLake, [10] reservoir, [11] gloWetland
    double RiverTemp = 0.;                              // river temperature [degC]
    double locLakeTemp = 0.;                            // locLake temperature [degC]
    double locWetlandTemp = 0.;                         // locWetland temperature [degC]
    double gloLakeTemp = 0.;                            // gloLake temperature [degC]
    double reservoirTemp = 0.;                          // reservoir temperature [degC]
    double gloWetlandTemp = 0.;                         // gloWetland temperature [degC]
    double latHeatVapor = 0.;                           // latent heat of vaporization [J/kg]
    double NetShortRadIce = 0.;                         // net shortwave radiation ice[W/m2]
    double NetShortRadWater = 0.;                       // net shortwave radiation water [W/m2]
    double NetLongRad = 0.;                             // net longwave radiation [W/m2]
    double NetRadIce = 0.;                              // net radiation ice [W/m2]
    double NetRadWater = 0.;                            // net radiation water  [W/m2]
    double MixTemp = -9999.;                            // mixing temperature of inflow and LocWetland [degC]
    double totCellInflow = 0.;							// total volume of inflows from other cells [km3]
    double cellEnergyInflow = 0.;						// energy coming from other cells via water inflow [J]
    double cellInflowTemp = 0.;                         // Temperature of cell inflow [degC]
    bool gw_used = false;                               // indicator if groundwater already flown into prior swb
    vector <double> returnvalues (2,0.);                // container for return values: Temp.index = 0, ice thickness index = 1

    // total river inflows from surrounding cells
    for (int i = 0; i < 9; i++) {
        totCellInflow += t.cellInflows[i];
    }

    // energy of inflows from surrounding cells
    for (int i = 0; i < 9; i++) {
        cellEnergyInflow += t.cellInflows[i] * rhoW * HeatCapW * (t.cellInflowsTemp[i] + convK);
    }

    if(totCellInflow > 0.)
        cellInflowTemp = cellEnergyInflow / (totCellInflow * rhoW * HeatCapW) - convK;

    // see daily.cpp line 364: // latent heat of evaporation
    if (t.AirTempC > 0.) {
        // latent heat of vaporization of water
        latHeatVapor = (2.501 - 0.002361 * t.AirTempC) * 1.e6;	// MJ -> J
    }
    else {
        latHeatVapor = 0.;	// after Beek et al. (2012)
    }

    // short-/longwave radiation
    NetLongRad = t.LongRad - 1. * 0.00000005670374419 * pow((t.AirTempC + convK), 4.);  // emissivity = 1; Stephan-Boltzmann-Const. = 0.00000005670374419 [W/(m2*K^4)]

    NetShortRadWater = t.ShortRad * (1. - albedoW);
    NetRadWater = (NetShortRadWater + NetLongRad);

    NetShortRadIce = t.ShortRad * (1. - albedoI);
    NetRadIce = (NetShortRadIce + NetLongRad);

    // locLake temperature
    if (t.locLakeArea > 0.) {
        returnvalues = calcLocLakeTemp(t, NetRadIce, NetRadWater, latHeatVapor);
        gw_used = true;
        locLakeTemp = returnvalues[0];      // temperature [degC]
        calcReturn[1] = locLakeTemp;
        calcReturn[0] = locLakeTemp;
        calcReturn[7] = returnvalues[1];    // ice thickness [km]
    }

    // locWetland temperature
    if (t.locWetlandArea > 0.){
        returnvalues = calcLocWetlandTemp(t, NetRadIce, NetRadWater, latHeatVapor, locLakeTemp, gw_used);
        gw_used = true;
        locWetlandTemp = returnvalues[0];       // temperature [degC]
        calcReturn[2] = locWetlandTemp;
        calcReturn[0] = locWetlandTemp;
        calcReturn[8] = returnvalues[1];        // ice thickness [km]
    }

    // gloLake temperature
    // mixing temperature of locWetland/locLake and cell inflow
    if (t.gloLakeArea > 0.){
        if ((totCellInflow >= 0.) && (t.outflowLocWetland > 0.)) {     // MixTemp River (if existent) & LocWetland (if existent)
            MixTemp = ((cellInflowTemp + convK) * totCellInflow + (locWetlandTemp + convK) * t.outflowLocWetland) / (totCellInflow + t.outflowLocWetland) - convK;
            if (MixTemp < 0.)
                MixTemp = 0.;
        }
        else if (totCellInflow > 0.){    // MixTemp River & locLake (if existent)
            MixTemp = ((cellInflowTemp + convK) * totCellInflow + (locLakeTemp + convK) * t.inflowLocWetland) / (totCellInflow + t.inflowLocWetland) - convK;
            if (MixTemp < 0.)
                MixTemp = 0.;
        }
        else if (t.inflowLocWetland > 0.)      // if only locLake
            MixTemp = locLakeTemp;

        else MixTemp = -9999.;   //if absolutely no inflow

        returnvalues = calcGloLakeTemp(t, NetRadIce, NetRadWater, latHeatVapor, MixTemp, gw_used);
        gw_used = true;
        gloLakeTemp = returnvalues[0];       // temperature [degC]
        calcReturn[3] = gloLakeTemp;
        calcReturn[0] = gloLakeTemp;
        calcReturn[9] = returnvalues[1];        // ice thickness [km]
    }

    // reservoir temperature
    if (t.reservoirArea > 0.){
        if (t.gloLakeArea <= 0.) {  // in no gloLake
            if ((totCellInflow >= 0.) && (t.outflowLocWetland > 0.)){     // MixTemp River (if existent) & LocWetland (if existent)
                MixTemp = ((cellInflowTemp + convK) * totCellInflow + (locWetlandTemp + convK) * t.outflowLocWetland) / (totCellInflow + t.outflowLocWetland) - convK;
                if (MixTemp < 0.)
                    MixTemp = 0.;
            }
            else if (totCellInflow > 0.){           // MixTemp River & locLake (if existent)
                MixTemp = ((cellInflowTemp + convK) * totCellInflow + (locLakeTemp + convK) * t.inflowLocWetland) / (totCellInflow + t.inflowLocWetland) - convK;
                if (MixTemp < 0.)
                    MixTemp = 0.;
            }
            else if (t.inflowLocWetland > 0.)      // if only locLake
                MixTemp = locLakeTemp;

            else MixTemp = -9999.;  // if absolutely no inflow
        }
        else MixTemp = gloLakeTemp;

        returnvalues = calcReservoirTemp(t, NetRadIce, NetRadWater, latHeatVapor, MixTemp, gw_used);
        gw_used = true;
        reservoirTemp = returnvalues[0];       // temperature [degC]
        calcReturn[4] = reservoirTemp;
        calcReturn[0] = reservoirTemp;
        calcReturn[10] = returnvalues[1];        // ice thickness [km]
    }

    // Global Wetland Temperature
    if (t.gloWetlandArea > 0.){
        if (t.reservoirArea > 0.)
            MixTemp = reservoirTemp;
        else if (t.gloLakeArea <= 0.) {     // if no gloLake
            if ((totCellInflow >= 0.) && (t.outflowLocWetland > 0.)){     // MixTemp River (if existent) & LocWetland (if existent)
                MixTemp = ((cellInflowTemp + convK) * totCellInflow + (locWetlandTemp + convK) * t.outflowLocWetland) / (totCellInflow + t.outflowLocWetland) - convK;
                if (MixTemp < 0.)
                    MixTemp = 0.;
            }
            else if (totCellInflow > 0.){     // MixTemp River & locLake (if existent)
                MixTemp = ((cellInflowTemp + convK) * totCellInflow + (locLakeTemp + convK) * t.inflowLocWetland) / (totCellInflow + t.inflowLocWetland) - convK;
                if (MixTemp < 0.)
                    MixTemp = 0.;
            }
            else if (t.inflowLocWetland > 0.)       // if only locLake
                MixTemp = locLakeTemp;

            else MixTemp = -9999.;   //if absolutely no inflow
        }
        else if (t.gloLakeArea > 0.)
            MixTemp = gloLakeTemp;

        returnvalues = calcGloWetlandTemp(t, NetRadIce, NetRadWater, latHeatVapor, MixTemp, gw_used);
        gw_used = true;
        gloWetlandTemp = returnvalues[0];       // temperature [degC]
        calcReturn[5] = gloWetlandTemp;
        calcReturn[0] = gloWetlandTemp;
        calcReturn[11] = returnvalues[1];       // ice thickness [km]
    }

    // River temperature
    if (t.hasOutflow || t.RiverLengthCheck != 55.) {    // cell has outflow? or RiverLength != 55 to check if river existent in cell (55 marks outlet or inland sink cells but is altered by meandering ratio see rout_prepare.cpp)
        if (t.gloWetlandArea > 0.)
            MixTemp = gloWetlandTemp;
        else if (t.reservoirArea > 0.)
            MixTemp = reservoirTemp;
        else if (t.gloLakeArea > 0.)
            MixTemp = gloLakeTemp;
        else if (t.gloLakeArea <= 0.) {
            if ((totCellInflow >= 0.) && (t.outflowLocWetland > 0.)){     // MixTemp River (if existent) & LocWetland (if existent)
                MixTemp = ((cellInflowTemp + convK) * totCellInflow + (locWetlandTemp + convK) * t.outflowLocWetland) / (totCellInflow + t.outflowLocWetland) - convK;
                if (MixTemp < 0.)
                    MixTemp = 0.;
            }
            else if (totCellInflow > 0.){     // MixTemp River & locLake (if existent)
                MixTemp = ((cellInflowTemp + convK) * totCellInflow + (locLakeTemp + convK) * t.inflowLocWetland) / (totCellInflow + t.inflowLocWetland) - convK;
                if (MixTemp < 0.)
                    MixTemp = 0.;
            }
            else if (t.inflowLocWetland > 0.)       // if only locLake
                MixTemp = locLakeTemp;

            else MixTemp = -9999.;   //if absolutely no inflow
        }

        returnvalues = calcRiverTemp(t, NetRadIce, NetRadWater, latHeatVapor, MixTemp, gw_used);

        calcReturn[0] = returnvalues[0];        // water temperature [degC]
        calcReturn[6] = returnvalues[1];        // ice thickness [km]
    }

	return calcReturn;
}

/**
* @brief calculates locLake water temperature
*
* calculation of local lake water temperature
* called by function calcWTemp
*
*
* @author Sebastian Ackermann
*
* @param const transfer &t: contains all variables needed from routing
* @param const double &NetRadIce: net radiation if ice is present
* @param const double &NetRadWater: net radiation if no ice is present
* @param const double &latHeatVapor: latent heat of evaporation
*
* @return local lake water temperature [degC], ice thickness [km]
*/
vector<double> calcWaterTempClass::calcLocLakeTemp(const transfer &t, const double &NetRadIce, const double &NetRadWater, const double &latHeatVapor){
    double locLakeArea_km2 = 0.;                 // locLake area [km2]
    double locLakeEnergyStor = 0.;               // stored energy in locLake [J]
    double locLakeEnergyRad = 0.;                // net short & longwave radiation energy [J]
    double locLakeSensHeatflux = 0.;             // sensible heat flux [J]
    double locLakeLatHeatflux = 0.;              // latent heat flux [J]
    double locLakePrecip = 0.;                   // precipitation on locLake fraction [km3]
    double locLakeEnergyPrecip = 0.;             // energy via precipitation [J]
    double locLakeSurfRunoff = 0.;               // surface runoff into locLake [km3]
    double locLakeGWRunoff = t.localGWRunoff;    // gw runoff into locLake [km3]
    double EnergyLocalSurfRunoff = 0.;           // energy via surface runoff [J]
    double EnergyLocalGWRunoff = 0.;             // energy via groundwater runoff [J]
    double locLaketotEnergy = 0.;                // total energy in locLake water [J]
    double locLakeVol = 0.;                      // locLake volume to calculate temperature [km3]
    double locLakeTemp = 0.;                     // return value temperature [degC]
    double iceThickness = t.iceThickness[1];     // ice thickness [km]
    double inflowTemp = 0.;                      // MixTemp of inflows [K]
    bool iceMax = false;                         // true if max. ice volume in lake and no inflows
    vector<double> returnvector (2,0.);          // returns values

    locLakeArea_km2 = t.CellArea * (t.locLakeArea / 100.) * t.locLakeAreaRedFactor; // [%] -> [km2]

    locLakePrecip = t.openWaterPrecip / 1000000. * locLakeArea_km2;       // [mm] -> [km3]

    // surface runoff into locLake [km3] from fractional routing in routing.cpp combines all flows from locLakes locWetlands and gloWetlands bc. of G_fswb_catchment
    locLakeSurfRunoff = t.localRunoff - t.localGWRunoff;

    // if cell has no river computation then localGWRunoffIntoRiver and surfaceRunoffIntoRiver from fractional routing must be routed into existing lakes/wetlands in cell
    if (!t.hasOutflow && t.RiverLengthCheck == 55.) {
        locLakeSurfRunoff += t.localRunoffIntoRiver - t.localGWRunoffIntoRiver;
        locLakeGWRunoff += t.localGWRunoffIntoRiver;
    }

    locLakeEnergyStor = t.locLakeStorPrevStep * rhoW * HeatCapW * (t.locLakeTempPrevStep + convK);

    // conditions for ice calculation: negative AirTemp or previous ice exists
    if(t.AirTempC <= 0. || iceThickness > 0.){
        iceThickness += ((-(K_hIce * t.locLakeTempPrevStep) + (K_h * (-t.AirTempC)) - NetRadIce) / (lamda_f * rhoW / 1.e9) * 86400.) / 1000.; // [km]
        if (iceThickness * locLakeArea_km2 >= t.locLakeStorPrevStep) {
            iceThickness = t.locLakeStorPrevStep / locLakeArea_km2;     // max ice thickness = current water depth --> completely frozen if no inflows
            if(locLakePrecip + locLakeSurfRunoff + t.localGWRunoff <= 0.)
                iceMax = true;
        }
    }

    if (!iceMax) {
        if (iceThickness > 0.) {
            locLakeSensHeatflux = K_hIce * t.locLakeTempPrevStep * 1.e6 * 86400. * locLakeArea_km2;     // [s] -> [day] and [m2] -> [km2]

            // energy storage - sensHeatflux
            locLaketotEnergy = locLakeEnergyStor - locLakeSensHeatflux;

            // volume balance
            locLakeVol = t.locLakeStorPrevStep;

            // not important if ice decreases water volume or not, same results
            locLakeTemp = locLaketotEnergy / (locLakeVol * rhoW * HeatCapW) - convK; // [degC]

            if (locLakeTemp < 0.)       // cannot be below 0 degC for MixTemp calculation
                locLakeTemp = 0.;

            //MixTemp of inflows
            if (locLakePrecip + locLakeSurfRunoff + t.localGWRunoff > 0.) {  // no inflow then no MixTemp
                inflowTemp = (locLakePrecip * (max(t.AirTempC - precipTempCorr, 0.) + convK) + locLakeSurfRunoff * (max(t.AirTempC - precipTempCorr, 0.)
                             + convK) + t.localGWRunoff * (t.GWTemp + convK)) / (locLakePrecip + locLakeSurfRunoff + t.localGWRunoff);    // [K]

                if ((locLakeVol - iceThickness * locLakeArea_km2) > 0.) {   // if no water left (completely frozen) --> locLakeTemp = inflowTemp
                    locLakeTemp = (inflowTemp * (locLakePrecip + locLakeSurfRunoff + t.localGWRunoff) + (locLakeVol - iceThickness * locLakeArea_km2) * (locLakeTemp + convK)) /
                                  ((locLakeVol - iceThickness * locLakeArea_km2) + (locLakePrecip + locLakeSurfRunoff + t.localGWRunoff)) - convK;
                }
                else
                    locLakeTemp = inflowTemp - convK;   // [degC]
            }
        }
        else {  // no ice

            EnergyLocalSurfRunoff = locLakeSurfRunoff * rhoW * HeatCapW * (max(t.AirTempC - precipTempCorr, 0.) + convK);

            locLakeEnergyPrecip = locLakePrecip * rhoW * HeatCapW * (max(t.AirTempC - precipTempCorr, 0.) + convK); // precip never below freezing temp. of 0 degC see Beek et al. (2012)

            EnergyLocalGWRunoff = locLakeGWRunoff * rhoW * HeatCapW * (t.GWTemp + convK);

            locLakeEnergyRad = NetRadWater * 1.e6 * 86400. * locLakeArea_km2;  // Radiation [W/m2] --> [s] -> [day] and [m2] -> [km2]

            locLakeSensHeatflux = K_h * (t.locLakeTempPrevStep - t.AirTempC) * 1.e6 * 86400. * locLakeArea_km2; // [s] -> [day] and [m2] -> [km2]

            locLakeLatHeatflux = latHeatVapor * rhoW * (t.dailyLocLakeEvapo / 1000000. * locLakeArea_km2);      // [mm] -> [km3]

            // energy balance
            locLaketotEnergy = locLakeEnergyStor + locLakeEnergyRad - locLakeLatHeatflux +
                               locLakeEnergyPrecip + EnergyLocalSurfRunoff + EnergyLocalGWRunoff;

            // volume balance
            locLakeVol = t.locLakeStorPrevStep + locLakePrecip + locLakeSurfRunoff + t.localGWRunoff -
                         (t.dailyLocLakeEvapo / 1000000. * locLakeArea_km2);
            // water temperature
            if (locLakeVol > 0.) {     // for calculation safety only, locLakeVol must not be below 0

                locLakeTemp = (locLaketotEnergy / (locLakeVol * rhoW * HeatCapW)) - convK;    // [K] -> [degC]

                locLaketotEnergy -= locLakeSensHeatflux;

                // sensible heatflux can only cool/heat water to air temperature; sensible heat uses water temp. from previous day -> different delta T
                // if locLakeTemp > AirTempC and positive sensible heatflux then only cooling to AirTempC is possible
                // if locLakeTemp < AirTempC and negative sensible heatflux then only warming to AirTempC is possible
                // implemented to counter instability
                // if volume very small computation can get unstable -> reset to air temp as last resort
                if (locLakeTemp >= t.AirTempC) {
                    locLakeTemp = (locLaketotEnergy / (locLakeVol * rhoW * HeatCapW)) - convK;    // [K] -> [degC]
                    if (locLakeTemp < t.AirTempC)
                        locLakeTemp = t.AirTempC;
                }
                else {  // (locLakeTemp < t.AirTempC)
                    locLakeTemp = (locLaketotEnergy / (locLakeVol * rhoW * HeatCapW)) - convK;    // [K] -> [degC]
                    if (locLakeTemp > t.AirTempC)
                        locLakeTemp = t.AirTempC;
                }
            }
            else
                locLakeTemp = -9999.;  // no volume in locLake
        }
    }
    else{
        if(iceThickness > 0.)   // really ice there or no volume to begin with (not asked at beginning)
      	    locLakeTemp = 0.;   // max ice thickness = current water depth & no inflows --> completely frozen --> Temp = 0
        else
      	    locLakeTemp = -9999.;   // no volume in locLake
    }

    if (locLakeTemp > 60.)     // highest ever measured air temp. 56 degC
        locLakeTemp = t.AirTempC;

    if(locLakeTemp < 1.e-10 && locLakeTemp != -9999.)     // to eliminate very small temperature values due to calculation errors e.g. 1.e-14 or negative temperatures
        locLakeTemp = 0.;

    if(iceThickness < 0.)
        iceThickness = 0.;      // if ice thickness < 0. then ice thickness = 0 for next step

    returnvector[0] = locLakeTemp;
    returnvector[1] = iceThickness;

    return returnvector;
}

/**
* @brief calculates locWetland water temperature
*
* calculation of local wetland water temperature
* called by function calcWTemp
*
*
* @author Sebastian Ackermann
*
* @param const transfer &t: contains all variables needed from routing
* @param const double &NetRadIce: net radiation if ice is present
* @param const double &NetRadWater: net radiation if no ice is present
* @param const double &latHeatVapor: latent heat of evaporation
* @param const double &locLakeTemp: temperature of locLake (if exists) in cell
* @param const bool &gw_used: flag indicating if groundwater already flown into prev. swb
*
* @return local wetland water temperature [degC], ice thickness [km]
*/
vector<double> calcWaterTempClass::calcLocWetlandTemp(const transfer &t, const double &NetRadIce, const double &NetRadWater, const double &latHeatVapor, const double &locLakeTemp, const bool &gw_used) {
    double locWetlandArea_km2 = 0.;              // locWetland area [km2]
    double locWetlandEnergyStor = 0.;            // stored energy in locWetland [J]
    double locWetlandEnergyRad = 0.;             // net short & longwave radiation energy [J]
    double locWetlandSensHeatflux = 0.;          // sensible heat flux [J]
    double locWetlandLatHeatflux = 0.;           // latent heat flux [J]
    double locWetlandPrecip = 0.;                // precipitation on locWetland fraction [km3]
    double locWetlandEnergyPrecip = 0.;          // energy via precipitation [J]
    double locWetlandtotEnergy = 0.;             // total energy in locWetland water [J]
    double inflowLocWetlandEnergy = 0.;          // Energy inflow from locLake [J]
    double locWetlandSurfRunoff = 0.;            // surface runoff into locWetland [km3]
    double locWetlandGWRunoff = 0.;              // gw runoff into locWetland [km3]
    double EnergyLocalSurfRunoff = 0.;           // energy via surface runoff [J]
    double EnergyLocalGWRunoff = 0.;             // energy via groundwater runoff [J]
    double locWetlandVol = 0.;                   // locWetland volume to calculate temperature [km3]
    double locWetlandTemp = 0.;                  // return value temperature [degC]
    double iceThickness = t.iceThickness[2];     // ice thickness [km]
    double inflowTemp = 0.;                      // MixTemp of inflows [K]
    bool iceMax = false;                         // true if max. ice volume in lake and no inflows
    double sumRunoff = 0.;                       // sum of gw and surface runoff; 0 if locLake present
    vector<double> returnvector(2,0.);           // returns values

    locWetlandArea_km2 = t.CellArea * (t.locWetlandArea / 100.) * t.locWetlandAreaRedFactor; // [%] -> [km2]

    locWetlandPrecip = t.openWaterPrecip / 1000000. * locWetlandArea_km2;       // [mm] -> [km3]

    locWetlandEnergyStor = t.locWetlandStorPrevStep * rhoW * HeatCapW * (t.locWetlandTempPrevStep + convK);

    // case if local runoff from fractional routing is not flowing into locLake
    if (!gw_used) {
        // if cell has no river computation then localGWRunoffIntoRiver and surfaceRunoffIntoRiver from fractional routing must be routed into existing lakes/wetlands in cell
        if (!t.hasOutflow && t.RiverLengthCheck == 55.) {
            locWetlandSurfRunoff += t.localRunoffIntoRiver - t.localGWRunoffIntoRiver;
            locWetlandGWRunoff += t.localGWRunoffIntoRiver;
        }
        locWetlandSurfRunoff += t.localRunoff - t.localGWRunoff;
        EnergyLocalSurfRunoff = locWetlandSurfRunoff * rhoW * HeatCapW * (max(t.AirTempC - precipTempCorr, 0.) + convK);
        locWetlandGWRunoff += t.localGWRunoff;
        EnergyLocalGWRunoff = locWetlandGWRunoff * rhoW * HeatCapW * (t.GWTemp + convK);
        sumRunoff = locWetlandGWRunoff + locWetlandSurfRunoff;
    }

    // conditions for ice calculation: negative AirTemp or previous ice exists
    if(t.AirTempC <= 0. || iceThickness > 0.){
        iceThickness += ((-(K_hIce * t.locWetlandTempPrevStep) + (K_h * (-t.AirTempC)) - NetRadIce) / (lamda_f * rhoW / 1.e9) * 86400.) / 1000.; // [km]
        if (iceThickness * locWetlandArea_km2 >= t.locWetlandStorPrevStep) {
            iceThickness = t.locWetlandStorPrevStep / locWetlandArea_km2;     // max ice thickness = current water depth --> completely frozen if no inflows
            if(locWetlandPrecip + sumRunoff + t.inflowLocWetland <= 0.)
                iceMax = true;
        }
    }

    if (!iceMax) {
        if (iceThickness > 0.){
            locWetlandSensHeatflux = K_hIce * t.locWetlandTempPrevStep * 1.e6 * 86400. * locWetlandArea_km2;     // [s] -> [day] and [m2] -> [km2]

            // energy storage - sensHeatflux
            locWetlandtotEnergy = locWetlandEnergyStor - locWetlandSensHeatflux;
            // volume balance
            locWetlandVol = t.locWetlandStorPrevStep;
            // not important if ice decreases water volume or not, same results
            locWetlandTemp = locWetlandtotEnergy / (locWetlandVol * rhoW * HeatCapW) - convK; // [degC]

            if (locWetlandTemp < 0.)       // cannot be below 0 degC for MixTemp calculation
                locWetlandTemp = 0.;

            //MixTemp of inflows
            if (locWetlandPrecip + sumRunoff + t.inflowLocWetland > 0.) {  // no inflow then no MixTemp
                inflowTemp = (locWetlandPrecip * (max(t.AirTempC - precipTempCorr, 0.) + convK) + locWetlandSurfRunoff * (max(t.AirTempC - precipTempCorr, 0.)
                              + convK) + locWetlandGWRunoff * (t.GWTemp + convK) + t.inflowLocWetland * (locLakeTemp + convK)) / (locWetlandPrecip + sumRunoff + t.inflowLocWetland);    // [K]

                if ((locWetlandVol - iceThickness * locWetlandArea_km2) > 0.) {   // if no water left (completely frozen) --> locWetlandTemp = inflowTemp
                    locWetlandTemp = (inflowTemp * (locWetlandPrecip + sumRunoff + t.inflowLocWetland) + (locWetlandVol - iceThickness * locWetlandArea_km2) * (locWetlandTemp + convK)) /
                                  ((locWetlandVol - iceThickness * locWetlandArea_km2) + (locWetlandPrecip + sumRunoff + t.inflowLocWetland)) - convK;  // [degC]
                }
                else
                    locWetlandTemp = inflowTemp - convK;   // [degC]
            }
        }
        else{   // no ice
            locWetlandEnergyPrecip = locWetlandPrecip * rhoW * HeatCapW * (max(t.AirTempC - precipTempCorr, 0.) + convK);

            locWetlandEnergyRad = NetRadWater * 1.e6 * 86400. * locWetlandArea_km2;  // Radiation [W/m2] --> [s] -> [day] and [m2] -> [km2]

            locWetlandSensHeatflux = K_h * ((t.locWetlandTempPrevStep + convK) - (t.AirTempC + convK)) * 1.e6 * 86400. * locWetlandArea_km2; // [m2] -> [km2] and [s] -> [day]

            locWetlandLatHeatflux = latHeatVapor * rhoW * (t.dailyLocWetlandEvapo / 1000000. * locWetlandArea_km2);      // [mm] -> [km3]

            if (t.locLakeArea > 0.)
                inflowLocWetlandEnergy = t.inflowLocWetland * rhoW * HeatCapW * (locLakeTemp + convK);

            locWetlandtotEnergy = locWetlandEnergyStor + locWetlandEnergyRad - locWetlandLatHeatflux +
                           locWetlandEnergyPrecip + inflowLocWetlandEnergy + EnergyLocalSurfRunoff + EnergyLocalGWRunoff;

            //water temperature
            locWetlandVol = t.locWetlandStorPrevStep + locWetlandPrecip - (t.dailyLocWetlandEvapo / 1000000. * locWetlandArea_km2) + t.inflowLocWetland;

            if (!gw_used)
                locWetlandVol += locWetlandSurfRunoff + locWetlandGWRunoff;

            if (locWetlandVol > 0.) {      // for calculation safety only, locWetlandVol must not be below 0
                locWetlandTemp = (locWetlandtotEnergy / (locWetlandVol * rhoW * HeatCapW)) - convK;    // [K] -> [degC]

                locWetlandtotEnergy -= locWetlandSensHeatflux;

                // sensible heatflux can only cool/heat water to air temperature; sensible heat uses water temp. from previous day -> different delta T
                // if locLakeTemp > AirTempC and positive sensible heatflux then only cooling to AirTempC is possible
                // if locLakeTemp < AirTempC and negative sensible heatflux then only warming to AirTempC is possible
                // implemented to counter instability
                // if volume very small computation can get unstable -> reset to air temp as last resort
                if (locWetlandTemp >= t.AirTempC) {
                    locWetlandTemp = (locWetlandtotEnergy / (locWetlandVol * rhoW * HeatCapW)) - convK;    // [K] -> [degC]
                    if (locWetlandTemp < t.AirTempC)
                        locWetlandTemp = t.AirTempC;
                }
                else {  // (locWetlandTemp < t.AirTempC)
                    locWetlandTemp = (locWetlandtotEnergy / (locWetlandVol * rhoW * HeatCapW)) - convK;    // [K] -> [degC]
                    if (locWetlandTemp > t.AirTempC)
                        locWetlandTemp = t.AirTempC;
                }
            }
            else
                locWetlandTemp = -9999.;    // no volume in locWetland
        }
    }
    else{
        if(iceThickness > 0.)   // really ice there or no volume to begin with (not asked at beginning)
            locWetlandTemp = 0.;   // max ice thickness = current water depth & no inflows --> completely frozen --> Temp = 0
        else
            locWetlandTemp = -9999.;    // no volume in locWetland
    }

    if (locWetlandTemp > 60.)     // highest ever measured air temp. 56 degC
        locWetlandTemp = t.AirTempC;

    if(locWetlandTemp < 1.e-10 && locWetlandTemp != -9999.)     // to eliminate very small temperature values due to calculation errors e.g. 1.e-14 or negative temperatures
        locWetlandTemp = 0.;

    if(iceThickness < 0.)
        iceThickness = 0.;      // if ice thickness < 0. then ice thickness = 0 for next step

    returnvector[0] = locWetlandTemp;
    returnvector[1] = iceThickness;

    return returnvector;
}

/**
* @brief calculates gloLake water temperature
*
* calculation of global lake water temperature
* called by function calcWTemp
*
*
* @author Sebastian Ackermann
*
* @param const transfer &t: contains all variables needed from routing
* @param const double &NetRadIce: net radiation if ice is present
* @param const double &NetRadWater: net radiation if no ice is present
* @param const double &latHeatVapor: latent heat of evaporation
* @param const double &MixTemp: mixing temperature of inflows into gloLake
* @param const bool &gw_used: flag indicating if groundwater already flown into prev. swb
*
* @return global lake water temperature [degC], ice thickness [km]
*/
vector<double> calcWaterTempClass::calcGloLakeTemp(const transfer &t, const double &NetRadIce, const double &NetRadWater, const double &latHeatVapor, const double &MixTemp, const bool &gw_used) {

    double gloLakeEnergyStor = 0.;               // stored energy in gloLake [J]
    double gloLakeEnergyRad = 0.;                // net short & longwave radiation energy [J]
    double gloLakeSensHeatflux = 0.;             // sensible heat flux [J]
    double gloLakeLatHeatflux = 0.;              // latent heat flux [J]
    double gloLakePrecip = 0.;                   // precipitation on gloLake [km3]
    double gloLakeEnergyPrecip = 0.;             // energy via precipitation [J]
    double inflowGloLakeEnergy = 0.;             // Energy of inflows [J]
    double gloLaketotEnergy = 0.;                // total energy in gloLake water [J]
    double gloLakeSurfRunoff = 0.;               // surface runoff into gloLake [km3]
    double gloLakeGWRunoff = 0.;                 // gw runoff into gloLake [km3]
    double EnergyLocalSurfRunoff = 0.;           // energy via surface runoff [J]
    double EnergyLocalGWRunoff = 0.;             // energy via groundwater runoff [J]
    double gloLakeStor = 0.;                     // volume stored in gloLake that can interact with atmosphere [km3]
    double gloLakeVol = 0.;                      // gloLake volume to calculate temperature [km3]
    double gloLakeTemp = 0.;                     // return value temperature [degC]
    double iceThickness = t.iceThickness[3];     // ice thickness [km]
    double inflowTemp = 0.;                      // MixTemp of inflows [K]
    bool iceMax = false;                         // true if max. ice volume in lake and no inflows
    double sumRunoff = 0.;                       // sum of gw and surface runoff; 0 if locLake present
    vector<double> returnvector(2,0.);           // returns values

    //calculating thermocline of gloLake(0) and ask which volume is smaller
    gloLakeStor =  min(t.gloLakeStorPrevStep, calcThermocline(t.gloLakeArea, 0));

    gloLakePrecip = t.openWaterPrecip / 1000000. * t.gloLakeArea;       // [mm] -> [km3]

    gloLakeEnergyStor = gloLakeStor * rhoW * HeatCapW * (t.gloLakeTempPrevStep + convK);

    // case if local runoff from fractional routing is not flowing into locLake or locWetland
    if (!gw_used) {
        // if cell has no river computation then localGWRunoffIntoRiver and surfaceRunoffIntoRiver from fractional routing must be routed into existing lakes/wetlands in cell
        if (!t.hasOutflow && t.RiverLengthCheck == 55.) {
            gloLakeSurfRunoff += t.localRunoffIntoRiver - t.localGWRunoffIntoRiver;
            gloLakeGWRunoff += t.localGWRunoffIntoRiver;
        }
        gloLakeSurfRunoff += t.localRunoff - t.localGWRunoff;
        EnergyLocalSurfRunoff = gloLakeSurfRunoff * rhoW * HeatCapW * (max(t.AirTempC - precipTempCorr, 0.) + convK);
        gloLakeGWRunoff += t.localGWRunoff;
        EnergyLocalGWRunoff = gloLakeGWRunoff * rhoW * HeatCapW * (t.GWTemp + convK);
        sumRunoff = gloLakeSurfRunoff + gloLakeGWRunoff;
    }

    // conditions for ice calculation: negative AirTemp or previous ice exists
    if(t.AirTempC <= 0. || iceThickness > 0.){
        iceThickness += ((-(K_hIce * t.gloLakeTempPrevStep) + (K_h * (-t.AirTempC)) - NetRadIce) / (lamda_f * rhoW / 1.e9) * 86400.) / 1000.; // [km]
        if (iceThickness * t.gloLakeArea >= gloLakeStor) {
            iceThickness = gloLakeStor / t.gloLakeArea;     // max ice thickness = current water depth --> completely frozen if no inflows
            if(gloLakePrecip + sumRunoff + t.inflowGloLake <= 0.)
                iceMax = true;
        }
    }

    if (!iceMax) {
        if (iceThickness > 0.){
            gloLakeSensHeatflux = K_hIce * t.gloLakeTempPrevStep * 1.e6 * 86400. * t.gloLakeArea;     // [s] -> [day] and [m2] -> [km2]

            // energy storage - sensHeatflux
            gloLaketotEnergy = gloLakeEnergyStor - gloLakeSensHeatflux;
            // volume balance
            gloLakeVol = gloLakeStor;
            // not important if ice decreases water volume or not, same results
            gloLakeTemp = gloLaketotEnergy / (gloLakeVol * rhoW * HeatCapW) - convK; // [degC]

            if (gloLakeTemp < 0.)       // cannot be below 0 degC for MixTemp calculation
                gloLakeTemp = 0.;

            // MixTemp of inflows
            if (gloLakePrecip + sumRunoff + t.inflowGloLake > 0.) {  // no inflow then no MixTemp
                if (MixTemp < 0.){      // <0 (-9999.) means no inflow from river/locWetland... (see function calcWTemp)
                    inflowTemp = (gloLakePrecip * (max(t.AirTempC - precipTempCorr, 0.) + convK) + gloLakeSurfRunoff * (max(t.AirTempC - precipTempCorr, 0.)
                        + convK) + gloLakeGWRunoff * (t.GWTemp + convK)) / (gloLakePrecip + sumRunoff);    // [K]
                }
                else {
                    inflowTemp = (gloLakePrecip * (max(t.AirTempC - precipTempCorr, 0.) + convK) + gloLakeSurfRunoff * (max(t.AirTempC - precipTempCorr, 0.)
                        + convK) + gloLakeGWRunoff * (t.GWTemp + convK) + t.inflowGloLake * (MixTemp + convK)) / (gloLakePrecip + sumRunoff + t.inflowGloLake);    // [K]
                }

                if ((gloLakeVol - iceThickness * t.gloLakeArea) > 0.) {   // if no water left (completely frozen) --> gloLakeTemp = inflowTemp
                    gloLakeTemp = (inflowTemp * (gloLakePrecip + sumRunoff + t.inflowGloLake) + (gloLakeVol - iceThickness * t.gloLakeArea) * (gloLakeTemp + convK)) /
                                     ((gloLakeVol - iceThickness * t.gloLakeArea) + (gloLakePrecip + sumRunoff + t.inflowGloLake)) - convK;  // [degC]
                }
                else
                    gloLakeTemp = inflowTemp - convK;   // [degC]
            }
        }
        else{   // no ice
            gloLakeEnergyPrecip = gloLakePrecip * rhoW * HeatCapW * (max(t.AirTempC - precipTempCorr, 0.) + convK);

            gloLakeEnergyRad = NetRadWater * 1.e6 * 86400. * t.gloLakeArea;  // Radiation [W/m2] --> [s] -> [day] and [m2] -> [km2]

            gloLakeSensHeatflux = K_h * ((t.gloLakeTempPrevStep + convK) - (t.AirTempC + convK)) * 1.e6 * 86400. *
                                  t.gloLakeArea; // [m2] -> [km2] and [s] -> [day]

            gloLakeLatHeatflux = latHeatVapor * rhoW * (t.dailyGloLakeEvapo / 1000000. * t.gloLakeArea);      // [mm] -> [km3]

            if (MixTemp >= 0.)     // if < 0. (-9999.) then no inflows (see function calcWTemp)
                inflowGloLakeEnergy = ((MixTemp + convK) * t.inflowGloLake * rhoW * HeatCapW);

            //energy balance without sensible Heatflux
            gloLaketotEnergy = gloLakeEnergyStor + gloLakeEnergyRad - gloLakeLatHeatflux +
                               gloLakeEnergyPrecip + inflowGloLakeEnergy + EnergyLocalSurfRunoff + EnergyLocalGWRunoff;

            //water temperature
            gloLakeVol = gloLakeStor + gloLakePrecip - (t.dailyGloLakeEvapo / 1000000. * t.gloLakeArea) + t.inflowGloLake;

            if (!gw_used)
                gloLakeVol += gloLakeSurfRunoff + gloLakeGWRunoff;

            if (gloLakeVol > 0.) {      // for calculation safety only, gloLakeVol must not be below 0
                gloLakeTemp = (gloLaketotEnergy / (gloLakeVol * rhoW * HeatCapW)) - convK;    // [K] -> [degC]

                //energy balance with sensible Heatflux
                gloLaketotEnergy -= gloLakeSensHeatflux;

                // sensible heatflux can only cool/heat water to air temperature; sensible heat uses water temp. from previous day -> different delta T
                // if gloLakeTemp > AirTempC and positive sensible heatflux then only cooling to AirTempC is possible
                // if gloLakeTemp < AirTempC and negative sensible heatflux then only warming to AirTempC is possible
                // implemented to counter instability
                // if volume very small computation can get unstable -> reset to air temp as last resort
                if(gloLakeTemp >= t.AirTempC){
                    gloLakeTemp = (gloLaketotEnergy / (gloLakeVol * rhoW * HeatCapW)) - convK;    // [K] -> [degC]
                    if(gloLakeTemp < t.AirTempC)
                        gloLakeTemp = t.AirTempC;
                }
                else {  //(gloLakeTemp < t.AirTempC)
                    gloLakeTemp = (gloLaketotEnergy / (gloLakeVol * rhoW * HeatCapW)) - convK;    // [K] -> [degC]
                    if(gloLakeTemp > t.AirTempC)
                        gloLakeTemp = t.AirTempC;
                }
            }
            else
                gloLakeTemp = -9999.;   // no volume in gloLake
        }
    }
    else{
        if(iceThickness > 0.)   // really ice there or no volume to begin with (not asked at beginning)
            gloLakeTemp = 0.;   // max ice thickness = current water depth & no inflows --> completely frozen --> Temp = 0. degC
        else
            gloLakeTemp = -9999.;   // no volume in gloLake
    }

    if (gloLakeTemp > 60.)     // highest ever measured air temp. 56 degC
        gloLakeTemp = t.AirTempC;

    if(gloLakeTemp < 1.e-10 && gloLakeTemp != -9999.)     // to eliminate very small temperature values due to calculation errors e.g. 1.e-14 or negative temperatures
        gloLakeTemp = 0.;

    if(iceThickness < 0.)
        iceThickness = 0.;      // if ice thickness < 0. then ice thickness = 0. for next step

    returnvector[0] = gloLakeTemp;
    returnvector[1] = iceThickness;

    return returnvector;
}

/**
* @brief calculates reservoir water temperature
*
* calculation of reservoir water temperature
* called by function calcWTemp
*
*
* @author Sebastian Ackermann
*
* @param const transfer &t: contains all variables needed from routing
* @param const double &NetRadIce: net radiation if ice is present
* @param const double &NetRadWater: net radiation if no ice is present
* @param const double &latHeatVapor: latent heat of evaporation
* @param const double &MixTemp: mixing temperature of inflows into reservoir
* @param const bool &gw_used: flag indicating if groundwater already flown into prev. swb
*
* @return reservoir water temperature [degC], ice thickness [km]
*/
vector<double> calcWaterTempClass::calcReservoirTemp(const transfer &t, const double &NetRadIce, const double &NetRadWater, const double &latHeatVapor, const double &MixTemp, const bool &gw_used){

    double ReservoirEnergyStor = 0.;             // stored energy in reservoir [J]
    double ReservoirEnergyRad = 0.;              // net short & longwave radiation energy [J]
    double ReservoirSensHeatflux = 0.;           // sensible heat flux [J]
    double ReservoirLatHeatflux = 0.;            // latent heat flux [J]
    double ReservoirPrecip = 0.;                 // precipitation on reservoir [km3]
    double ReservoirEnergyPrecip = 0.;           // energy via precipitation [J]
    double inflowReservoirEnergy = 0.;           // Energy of inflows [J]
    double ReservoirtotEnergy = 0.;              // total energy in reservoir water [J]
    double ReservoirSurfRunoff = 0.;             // surface runoff into reservoir [km3]
    double ReservoirGWRunoff = 0.;               // gw runoff into reservoir [km3]
    double EnergyLocalSurfRunoff = 0.;           // energy via surface runoff [J]
    double EnergyLocalGWRunoff = 0.;             // energy via groundwater runoff [J]
    double ReservoirStor = 0.;                   // volume stored in reservoir that can interact with atmosphere [km3]
    double ReservoirVol = 0.;                    // reservoir volume to calculate temperature [km3]
    double ReservoirTemp = 0.;                   // return value temperature [degC]
    double iceThickness = t.iceThickness[4];     // ice thickness [km]
    double inflowTemp = 0.;                      // MixTemp of inflows [K]
    bool iceMax = false;                         // true if max. ice volume in lake and no inflows
    double sumRunoff = 0.;                       // sum of gw and surface runoff; 0 if locLake present
    vector<double> returnvector(2,0.);           // returns values

    //calculating thermocline of reservoir(1) and ask which volume is smaller
    ReservoirStor = min(t.reservoirStorPrevStep, calcThermocline(t.reservoirArea, 1));

    ReservoirEnergyStor = ReservoirStor * rhoW * HeatCapW * (t.reservoirTempPrevStep + convK);

    ReservoirPrecip = t.openWaterPrecip / 1000000. * t.reservoirArea;       // [mm] -> [km3]

    // case if local runoff from fractional routing is not flowing into locLake, locWetland or gloLake
    if (!gw_used) {
        // if cell has no river computation then localGWRunoffIntoRiver and surfaceRunoffIntoRiver from fractional routing must be routed into existing lakes/wetlands in cell
        if (!t.hasOutflow && t.RiverLengthCheck == 55.) {
            ReservoirSurfRunoff += t.localRunoffIntoRiver - t.localGWRunoffIntoRiver;
            ReservoirGWRunoff += t.localGWRunoffIntoRiver;
        }
        ReservoirSurfRunoff += t.localRunoff - t.localGWRunoff;
        EnergyLocalSurfRunoff = ReservoirSurfRunoff * rhoW * HeatCapW * (max(t.AirTempC - precipTempCorr, 0.) + convK);
        ReservoirGWRunoff += t.localGWRunoff;
        EnergyLocalGWRunoff = ReservoirGWRunoff * rhoW * HeatCapW * (t.GWTemp + convK);
        sumRunoff = ReservoirSurfRunoff + ReservoirGWRunoff;
    }

    // conditions for ice calculation: negative AirTemp or previous ice exists
    if(t.AirTempC <= 0. || iceThickness > 0.){
        iceThickness += ((-(K_hIce * t.reservoirTempPrevStep) + (K_h * (-t.AirTempC)) - NetRadIce) / (lamda_f * rhoW / 1.e9) * 86400.) / 1000.; // [km]
        if (iceThickness * t.reservoirArea >= ReservoirStor) {
            iceThickness = ReservoirStor / t.reservoirArea;     // max ice thickness = current water depth --> completely frozen if no inflows
            if(ReservoirPrecip + sumRunoff + t.inflowReservoir <= 0.)
                iceMax = true;
        }
    }

    if (!iceMax) {
        if (iceThickness > 0.){
            ReservoirSensHeatflux = K_hIce * t.reservoirTempPrevStep * 1.e6 * 86400. * t.reservoirArea;     // [s] -> [day] and [m2] -> [km2]

            // energy storage - sensHeatflux
            ReservoirtotEnergy = ReservoirEnergyStor - ReservoirSensHeatflux;
            // volume balance
            ReservoirVol = ReservoirStor;
            // not important if ice decreases water volume or not, same results
            ReservoirTemp = ReservoirtotEnergy / (ReservoirVol * rhoW * HeatCapW) - convK; // [degC]

            if (ReservoirTemp < 0.)       // cannot be below 0 degC for MixTemp calculation
                ReservoirTemp = 0.;

            // MixTemp of inflows
            if (ReservoirPrecip + sumRunoff + t.inflowReservoir > 0.) {  // no inflow then no MixTemp
                if (MixTemp < 0.) {      // < 0. (-9999.) means no inflow from river/locWetland... (see function calcWTemp)
                    inflowTemp = (ReservoirPrecip * (max(t.AirTempC - precipTempCorr, 0.) + convK) + ReservoirSurfRunoff * (max(t.AirTempC - precipTempCorr, 0.)
                                   + convK) + ReservoirGWRunoff * (t.GWTemp + convK)) / (ReservoirPrecip + sumRunoff);    // [K]
                }
                else{
                    inflowTemp = (ReservoirPrecip * (max(t.AirTempC - precipTempCorr, 0.) + convK) + ReservoirSurfRunoff * (max(t.AirTempC - precipTempCorr, 0.)
                                   + convK) + ReservoirGWRunoff * (t.GWTemp + convK) + t.inflowReservoir * (MixTemp + convK)) / (ReservoirPrecip + sumRunoff + t.inflowReservoir);    // [K]
                }

                if ((ReservoirVol - iceThickness * t.reservoirArea) > 0.) {   // if no water left (completely frozen) --> reservoirTemp = inflowTemp
                    ReservoirTemp = (inflowTemp * (ReservoirPrecip + sumRunoff + t.inflowReservoir) + (ReservoirVol - iceThickness * t.reservoirArea) * (ReservoirTemp + convK)) /
                                    ((ReservoirVol - iceThickness * t.reservoirArea) + (ReservoirPrecip + sumRunoff + t.inflowReservoir)) - convK;  // [degC]
                }
                else
                    ReservoirTemp = inflowTemp - convK;   // [degC]
            }
        }
        else{   // no ice
            ReservoirEnergyRad = NetRadWater * 1.e6 * 86400. * t.reservoirArea;  // Radiation [W/m2] --> [s] -> [day] and [m2] -> [km2]

            ReservoirSensHeatflux = K_h * ((t.reservoirTempPrevStep + convK) - (t.AirTempC + convK)) * 1.e6 * 86400. * t.reservoirArea; // [m2] -> [km2] and [s] -> [day]

            ReservoirLatHeatflux = latHeatVapor * rhoW * (t.dailyReservoirEvapo / 1000000. * t.reservoirArea);      // [mm] -> [km3]

            ReservoirEnergyPrecip = ReservoirPrecip * rhoW * HeatCapW * (max(t.AirTempC - precipTempCorr, 0.) + convK);

            if (MixTemp >= 0.)     // if < 0. (-9999.) then no inflows (see function calcWTemp)
                inflowReservoirEnergy = ((MixTemp + convK) * t.inflowReservoir * rhoW * HeatCapW);

            //energy balance without sensible heatflux
            ReservoirtotEnergy = ReservoirEnergyStor + ReservoirEnergyRad - ReservoirLatHeatflux +
                                 ReservoirEnergyPrecip + inflowReservoirEnergy + EnergyLocalSurfRunoff + EnergyLocalGWRunoff;

            //water temperature
            ReservoirVol = ReservoirStor + ReservoirPrecip - (t.dailyReservoirEvapo / 1000000. * t.reservoirArea) + t.inflowReservoir;
            //ReservoirVol = t.reservoirStorPrevStep + ReservoirPrecip - (t.dailyReservoirEvapo / 1000000. * t.reservoirArea) + t.inflowReservoir;

            if (!gw_used)
                ReservoirVol += ReservoirSurfRunoff + ReservoirGWRunoff;

            if (ReservoirVol > 0.) {      // for calculation safety only, reservoirVol must not be below 0
                ReservoirTemp = (ReservoirtotEnergy / (ReservoirVol * rhoW * HeatCapW)) - convK;    // [K] -> [degC]
                //energy balance with sensible heatflux
                ReservoirtotEnergy -= ReservoirSensHeatflux;
                //explanation see gloLake
                if(ReservoirTemp >= t.AirTempC){
                    ReservoirTemp = (ReservoirtotEnergy / (ReservoirVol * rhoW * HeatCapW)) - convK;    // [K] -> [degC]
                    if(ReservoirTemp < t.AirTempC)
                        ReservoirTemp = t.AirTempC;
                }
                else {   //(ReservoirTemp < t.AirTempC)
                    ReservoirTemp = (ReservoirtotEnergy / (ReservoirVol * rhoW * HeatCapW)) - convK;    // [K] -> [degC]
                    if(ReservoirTemp > t.AirTempC)
                        ReservoirTemp = t.AirTempC;
                }
            }
            else
                ReservoirTemp = -9999.;     // no volume in reservoir
        }
    }
    else{
        if(iceThickness > 0.)   // really ice there or no volume to begin with (not asked at beginning)
            ReservoirTemp = 0.;   // max ice thickness = current water depth & no inflows --> completely frozen --> Temp = 0. degC
        else
            ReservoirTemp = -9999.;     // no volume in reservoir
    }

    if (ReservoirTemp > 60.)     // highest ever measured air temp. 56 degC
        ReservoirTemp = t.AirTempC;

    if(ReservoirTemp < 1.e-10 && ReservoirTemp != -9999.)     // to eliminate very small temperature values due to calculation errors e.g. 1.e-14 or negative temperatures
        ReservoirTemp = 0.;

    if(iceThickness < 0.)
        iceThickness = 0.;      // if ice thickness < 0. then ice thickness = 0. for next step

    returnvector[0] = ReservoirTemp;
    returnvector[1] = iceThickness;

    return returnvector;
}

/**
* @brief calculates gloWetland water temperature
*
* calculation of global wetland water temperature
* called by function calcWTemp
*
*
* @author Sebastian Ackermann
*
* @param const transfer &t: contains all variables needed from routing
* @param const double &NetRadIce: net radiation if ice is present
* @param const double &NetRadWater: net radiation if no ice is present
* @param const double &latHeatVapor: latent heat of evaporation
* @param const double &MixTemp: mixing temperature of inflows into gloWetland
* @param const bool &gw_used: flag indicating if groundwater already flown into prev. swb
*
* @return global wetland water temperature [degC], ice thickness [km]
*/
vector<double> calcWaterTempClass::calcGloWetlandTemp(const transfer &t, const double &NetRadIce, const double &NetRadWater, const double &latHeatVapor, const double &MixTemp, const bool &gw_used){

    double gloWetlandEnergyStor = 0.;            // stored energy in gloWetland [J]
    double gloWetlandEnergyRad = 0.;             // net short & longwave radiation energy [J]
    double gloWetlandSensHeatflux = 0.;          // sensible heat flux [J]
    double gloWetlandLatHeatflux = 0.;           // latent heat flux [J]
    double gloWetlandPrecip = 0.;                // precipitation on gloWetland [km3]
    double gloWetlandEnergyPrecip = 0.;          // energy via precipitation [J]
    double inflowGloWetlandEnergy = 0.;          // Energy of inflows [J]
    double gloWetlandtotEnergy = 0.;             // total energy in gloWetland water [J]
    double gloWetlandSurfRunoff = 0.;            // surface runoff into gloWetland [km3]
    double gloWetlandGWRunoff = 0.;              // gw runoff into gloWetland [km3]
    double EnergyLocalSurfRunoff = 0.;           // energy via surface runoff [J]
    double EnergyLocalGWRunoff = 0.;             // energy via groundwater runoff [J]
    double gloWetlandArea_km2 = 0.;              // gloWetland Area [km2]
    double gloWetlandVol = 0.;                   // gloWetland volume to calculate temperature [km3]
    double gloWetlandTemp = 0.;                  // return value temperature [degC]
    double iceThickness = t.iceThickness[5];     // ice thickness [km]
    double inflowTemp = 0.;                      // MixTemp of inflows [K]
    bool iceMax = false;                         // true if max. ice volume in lake and no inflows
    double sumRunoff = 0.;                       // sum of gw and surface runoff; 0 if locLake present
    vector<double> returnvector(2, 0.);          // returns values

    gloWetlandArea_km2 = t.CellArea * (t.gloWetlandArea / 100.) * t.gloWetlandAreaReductionFactor; // [%] -> [km2]

    gloWetlandEnergyStor = t.gloWetlandStorPrevStep * rhoW * HeatCapW * (t.gloWetlandTempPrevStep + convK);

    gloWetlandPrecip = t.openWaterPrecip / 1000000. * gloWetlandArea_km2;       // [mm] -> [km3]

    // case if local runoff from fractional routing is not flowing into locLake, locWetland, gloLake, reservoir
    if (!gw_used) {
        // if cell has no river computation then localGWRunoffIntoRiver and surfaceRunoffIntoRiver from fractional routing must be routed into existing lakes/wetlands in cell
        if (!t.hasOutflow && t.RiverLengthCheck == 55.) {
            gloWetlandSurfRunoff += t.localRunoffIntoRiver - t.localGWRunoffIntoRiver;
            gloWetlandGWRunoff += t.localGWRunoffIntoRiver;
        }
        gloWetlandSurfRunoff += t.localRunoff - t.localGWRunoff;
        EnergyLocalSurfRunoff = gloWetlandSurfRunoff * rhoW * HeatCapW * (max(t.AirTempC - precipTempCorr, 0.) + convK);
        gloWetlandGWRunoff += t.localGWRunoff;
        EnergyLocalGWRunoff = gloWetlandGWRunoff * rhoW * HeatCapW * (t.GWTemp + convK);
        sumRunoff = gloWetlandSurfRunoff + gloWetlandGWRunoff;
    }

    // conditions for ice calculation: negative AirTemp or previous ice exists
    if(t.AirTempC <= 0. || iceThickness > 0.){
        iceThickness += ((-(K_hIce * t.gloWetlandTempPrevStep) + (K_h * (-t.AirTempC)) - NetRadIce) / (lamda_f * rhoW / 1.e9) * 86400.) / 1000.; // [km]
        if (iceThickness * gloWetlandArea_km2 >= t.gloWetlandStorPrevStep) {
            iceThickness = t.gloWetlandStorPrevStep / gloWetlandArea_km2;     // max ice thickness = current water depth --> completely frozen if no inflows
            if(gloWetlandPrecip + sumRunoff + t.gloWetlandInflow <= 0.)
                iceMax = true;
        }
    }

    if (!iceMax) {
        if (iceThickness > 0.){
            gloWetlandSensHeatflux = K_hIce * t.gloWetlandTempPrevStep  * 1.e6 * 86400. * gloWetlandArea_km2;     // [s] -> [day] and [m2] -> [km2]
            // energy storage - sensHeatflux
            gloWetlandtotEnergy = gloWetlandEnergyStor - gloWetlandSensHeatflux;
            // volume balance
            gloWetlandVol = t.gloWetlandStorPrevStep;
            // not important if ice decreases water volume or not, same results
            gloWetlandTemp = gloWetlandtotEnergy / (gloWetlandVol * rhoW * HeatCapW) - convK; // [degC]

            if (gloWetlandTemp < 0.)       // cannot be below 0. degC for MixTemp calculation
                gloWetlandTemp = 0.;

            // MixTemp of inflows
            if (gloWetlandPrecip + sumRunoff + t.gloWetlandInflow > 0.) {  // no inflow then no MixTemp
                if (MixTemp < 0.) {      // < 0. (-9999.) means no inflow from river/locWetland... (see function calcWTemp)
                    inflowTemp = (gloWetlandPrecip * (max(t.AirTempC - precipTempCorr, 0.) + convK) + gloWetlandSurfRunoff * (max(t.AirTempC - precipTempCorr, 0.)
                                   + convK) + gloWetlandGWRunoff * (t.GWTemp + convK)) / (gloWetlandPrecip + sumRunoff);    // [K]
                }
                else{
                    inflowTemp = (gloWetlandPrecip * (max(t.AirTempC - precipTempCorr, 0.) + convK) + gloWetlandSurfRunoff * (max(t.AirTempC - precipTempCorr, 0.)
                                   + convK) + gloWetlandGWRunoff * (t.GWTemp + convK) + t.gloWetlandInflow * (MixTemp + convK)) / (gloWetlandPrecip + sumRunoff + t.gloWetlandInflow);    // [K]
                }

                if ((gloWetlandVol - iceThickness * gloWetlandArea_km2) > 0.) {   // if no water left (completely frozen) --> gloWetlandTemp = inflowTemp
                    gloWetlandTemp = (inflowTemp * (gloWetlandPrecip + sumRunoff + t.gloWetlandInflow) + (gloWetlandVol - iceThickness * gloWetlandArea_km2) * (gloWetlandTemp + convK)) /
                                     ((gloWetlandVol - iceThickness * gloWetlandArea_km2) + (gloWetlandPrecip + sumRunoff + t.gloWetlandInflow)) - convK;  // [degC]
                }
                else
                    gloWetlandTemp = inflowTemp - convK;   // [degC]
            }
        }
        else{   // no ice
            gloWetlandEnergyRad = NetRadWater * 1.e6 * 86400. * gloWetlandArea_km2;  // Radiation [W/m2] --> [s] -> [day] and [m2] -> [km2]

            gloWetlandSensHeatflux = K_h * ((t.gloWetlandTempPrevStep + convK) - (t.AirTempC + convK)) * 1.e6 * 86400. * gloWetlandArea_km2; // [m2] -> [km2] and [s] -> [day]

            gloWetlandLatHeatflux = latHeatVapor * rhoW * (t.dailyGloWetlandEvapo / 1000000. * gloWetlandArea_km2);      // [mm] -> [km3]

            gloWetlandEnergyPrecip = gloWetlandPrecip * rhoW * HeatCapW * (max(t.AirTempC - precipTempCorr, 0.) + convK);

            if (MixTemp >= 0.)     // if < 0. (-9999.) then no inflows (see function calcWTemp)
                inflowGloWetlandEnergy = ((MixTemp + convK) * t.gloWetlandInflow * rhoW * HeatCapW);


            //energy balance without sensible Heatflux
            gloWetlandtotEnergy = gloWetlandEnergyStor + gloWetlandEnergyRad - gloWetlandLatHeatflux +
                                  gloWetlandEnergyPrecip + inflowGloWetlandEnergy + EnergyLocalSurfRunoff + EnergyLocalGWRunoff;

            gloWetlandVol = t.gloWetlandStorPrevStep + gloWetlandPrecip - (t.dailyGloWetlandEvapo / 1000000. * gloWetlandArea_km2) + t.gloWetlandInflow;
            if (!gw_used)
                gloWetlandVol += gloWetlandSurfRunoff + gloWetlandGWRunoff;

            if (gloWetlandVol > 0.) {      // for calculation safety only, gloWetlandVol must not be below 0
                gloWetlandTemp = (gloWetlandtotEnergy / (gloWetlandVol * rhoW * HeatCapW)) - convK;    // [K] -> [degC]
                //energy balance with sensible Heatflux
                gloWetlandtotEnergy -= gloWetlandSensHeatflux;
                //explanation see gloLake
                if(gloWetlandTemp >= t.AirTempC){
                    gloWetlandTemp = (gloWetlandtotEnergy / (gloWetlandVol * rhoW * HeatCapW)) - convK;    // [K] -> [degC]
                    if(gloWetlandTemp < t.AirTempC)
                        gloWetlandTemp = t.AirTempC;
                }
                else {   //(gloWetlandTemp < t.AirTempC)
                    gloWetlandTemp = (gloWetlandtotEnergy / (gloWetlandVol * rhoW * HeatCapW)) - convK;    // [K] -> [degC]
                    if(gloWetlandTemp > t.AirTempC)
                        gloWetlandTemp = t.AirTempC;
                }
            }
            else
                gloWetlandTemp = -9999.;    // no volume in gloWetland
        }
    }
    else{
        if(iceThickness > 0.)   // really ice there or no volume to begin with (not asked at beginning)
            gloWetlandTemp = 0.;   // max ice thickness = current water depth & no inflows --> completely frozen --> Temp = 0. degC
        else
            gloWetlandTemp = -9999.;    // no volume in gloWetland
    }

    if (gloWetlandTemp > 60.)     // highest ever measured air temp. 56 degC
        gloWetlandTemp = t.AirTempC;

    if(gloWetlandTemp < 1.e-10 && gloWetlandTemp != -9999.)     // to eliminate very small temperature values due to calculation errors e.g. 1.e-14 or negative temperatures
        gloWetlandTemp = 0.;

    if(iceThickness < 0.)
        iceThickness = 0.;      // if ice thickness < 0. then ice thickness = 0. for next step

    returnvector[0] = gloWetlandTemp;
    returnvector[1] = iceThickness;

    return returnvector;
}

/**
* @brief calculates river water temperature
*
* calculation of river water temperature
* called by function calcWTemp
*
*
* @author Sebastian Ackermann
*
* @param const transfer &t: contains all variables needed from routing
* @param const double &NetRadIce: net radiation if ice is present
* @param const double &NetRadWater: net radiation if no ice is present
* @param const double &latHeatVapor: latent heat of evaporation
* @param const double &MixTemp: mixing temperature of inflows into river
* @param const bool &gw_used: flag indicating if groundwater already flown into prev. swb
*
* @return river water temperature [degC], ice thickness [km]
*/
vector<double> calcWaterTempClass::calcRiverTemp(const transfer &t, const double &NetRadIce, const double &NetRadWater, const double &latHeatVapor, const double &MixTemp, const bool &gw_used){
    //variables
    double RiverArea = 0.;								// area of river [km2]
    double RiverDepth = 0.;								// river depth [km]
    double crossSectionalArea = 0.;						// cross section area of river [km2]
    double RiverPrecip = 0.;							// precipitation on river fraction [km3]
    double RiverSurfRunoff = 0.;						// surface runoff into river [km3]
    double RiverGWRunoff = t.localRunoffIntoRiver;      // gw runoff into river [km3]
    double RiverEvapo = 0.;								// river evaporation [km3]
    double RiverEnergyPrecip = 0.;						// energy via precipitation [J]
    double RiverEnergySurfRunoff = 0.;					// energy via surface runoff [J]
    double RiverEnergyGWRunoff = 0.;					// energy via groundwater runoff [J]
    double RiverEnergyRad = 0.;							// net short & longwave radiation energy [J]
    double sensHeatflux = 0.;							// sensible heat flux [J]
    double latHeatflux = 0.;							// latent heat flux [J]
    double RiverEnergyStor = 0.;						// stored energy in river [J]
    double inflowRiverEnergy = 0.;                      // energy from inflow [J]
    double RiverTotEnergy = 0.;							// total energy in river water [J]
    double RiverVol = 0.;                               // river volume to calculate temperature [km3]
    double RiverWTemp = 0.;								// return value: water temperature river [degC]
    double iceThickness = t.iceThickness[0];            // ice thickness [km]
    bool iceMax = false;                                // true if max. ice volume in lake and no inflows
    double inflowTemp = 0.;                             // MixTemp of inflows [K]
    double withdrawal = 0.;                             // withdrawal out of river for power plant cooling [km3]
    double powerPlantReturnFlow = 0.;                   // return flow into river from power plant [km3]
    double PPReturnFlowTempIncr = 3.;                   // temperature increase of return flow from power plant [degC]
    vector<double> returnvector(2, 0.);                 // returns values


    crossSectionalArea = t.RiverStor / t.RiverLength;    // [km3/km = km2]

    if (crossSectionalArea > 0.) {
        // computation of river depth see routing.cpp
        RiverDepth = -t.RiverBottomWidth / (4. * 1000.) + sqrt(t.RiverBottomWidth / 1000. * t.RiverBottomWidth / (16. * 1000.) + 0.5 * crossSectionalArea); // riverBottomWidth: [m] -> [km]
        // computation of river area [km2]
        RiverArea = t.RiverStor / RiverDepth;
    }
    else {
        RiverDepth = 0.;
        RiverArea = 0.;
    }

    // stored energy
    RiverEnergyStor = t.RiverStor * rhoW * HeatCapW * (t.cellInflowsTemp[4] + convK);

    RiverPrecip = t.openWaterPrecip * RiverArea / 1000000.; //[mm] -> [km3]

    RiverSurfRunoff = t.localRunoffIntoRiver - t.localGWRunoffIntoRiver;

    // if no lakes, wetlands, reservoirs at all in cell
    if (!gw_used) {
        RiverSurfRunoff += (t.localRunoff - t.localGWRunoff);
        RiverGWRunoff += t.localGWRunoff;     // so both volumes (surf and gw) are accounted for in total volume
    }

    // conditions for ice calculation: negative AirTemp or previous ice exists
    if(t.AirTempC <= 0. || iceThickness >= 5.e-6){   // ice thickness greater than 5 mm to simulate ice breakup in moving water
        if(RiverArea > 0.) {
            iceThickness += ((-(K_hIce * t.RiverStor) + (K_h * (-t.AirTempC)) - NetRadIce) / (lamda_f * rhoW / 1.e9) * 86400.) / 1000.; // [km]

            if (iceThickness * RiverArea >= t.RiverStor) {
                iceThickness = t.RiverStor / RiverArea;     // max ice thickness = current water depth --> completely frozen if no inflows
                if (RiverPrecip + RiverSurfRunoff + RiverGWRunoff + t.RiverInflow <= 0.)
                    iceMax = true;
            }
        }
        else
            iceThickness = 0.;
    }

    if (!iceMax) {
        if (iceThickness >= 5.e-6) {
            sensHeatflux = K_hIce * t.cellInflowsTemp[4] * 1.e6 * 86400. * RiverArea;     // [s] -> [day] and [m2] -> [km2]

            // energy storage - sensHeatflux
            RiverTotEnergy = RiverEnergyStor - sensHeatflux;
            // volume balance
            RiverVol = t.RiverStor;
            // not important if ice decreases water volume or not, same results
            RiverWTemp = RiverTotEnergy / (RiverVol * rhoW * HeatCapW) - convK; // [degC]

            if (RiverWTemp < 0.)       // cannot be below 0. degC for MixTemp calculation
                RiverWTemp = 0.;

            // MixTemp of inflows
            if (RiverPrecip + RiverSurfRunoff + RiverGWRunoff + t.RiverInflow > 0.) {  // no inflow then no MixTemp
                if (MixTemp < 0.) {      // < 0. (-9999.) means no inflow from river/locWetland ... (see function calcWTemp)
                    inflowTemp = (RiverPrecip * (max(t.AirTempC - precipTempCorr, 0.) + convK) + RiverSurfRunoff * (max(t.AirTempC - precipTempCorr, 0.)
                                   + convK) + RiverGWRunoff * (t.GWTemp + convK)) / (RiverPrecip + RiverSurfRunoff + RiverGWRunoff);    // [K]
                }
                else {
                    inflowTemp = (RiverPrecip * (max(t.AirTempC - precipTempCorr, 0.) + convK) + RiverSurfRunoff * (max(t.AirTempC - precipTempCorr, 0.)
                                   + convK) + RiverGWRunoff * (t.GWTemp + convK) + t.RiverInflow * (MixTemp + convK)) / (RiverPrecip + RiverSurfRunoff + RiverGWRunoff + t.RiverInflow);    // [K]
                }

                if ((RiverVol - iceThickness * RiverArea) > 0.) {   // if no water left (completely frozen) --> locLakeTemp = inflowTemp
                    RiverWTemp = (inflowTemp * (RiverPrecip + RiverSurfRunoff + RiverGWRunoff + t.RiverInflow) + (RiverVol - iceThickness * RiverArea) * (RiverWTemp + convK)) /
                                 ((RiverVol - iceThickness * RiverArea) + (RiverPrecip + RiverSurfRunoff + RiverGWRunoff + t.RiverInflow)) - convK;  // [degC]
                }
                else
                    RiverWTemp = inflowTemp - convK;   // [degC]
            }
        }
        else {  // no ice
            RiverEnergySurfRunoff = RiverSurfRunoff * rhoW * HeatCapW * (max(t.AirTempC - precipTempCorr, 0.) + convK);

            RiverEnergyGWRunoff = RiverGWRunoff * rhoW * HeatCapW * (t.GWTemp + convK);

            RiverEvapo = ((1.0 - t.CellCorrFact) * t.openWaterPrecip + (t.CellCorrFact * t.openWaterPET)) * RiverArea / 1000000.;        // see routing.cpp:3121

            if (RiverEvapo >= t.RiverStor + t.RiverInflow + RiverPrecip + RiverSurfRunoff + RiverGWRunoff) {   // only present water and inflows can evaporate not more
                returnvector[0] = -9999.;   // Evapo greater or equal than total water available in river, then no water in river for computation
                returnvector[1] = 0.;       // needed bc. evaporation not included in river volume computation in WaterGAP (as of May 2020)
                return returnvector;
            }

            RiverEnergyPrecip = RiverPrecip * rhoW * HeatCapW * (max(t.AirTempC - precipTempCorr, 0.) + convK); // precip never below 0. degC

            // energy from radiation
            RiverEnergyRad = NetRadWater * 1.e6 * 86400. * RiverArea;   // Radiation [W/m2] --> [s] -> [day] and [m2] -> [km2]

            // heat fluxes
            sensHeatflux = K_h * ((t.cellInflowsTemp[4] + convK) - (t.AirTempC + convK)) * 1.e6 * 86400. * RiverArea;        // [m2] -> [km2] and [s] -> [day]

            latHeatflux = latHeatVapor * rhoW * RiverEvapo;

            // energy from inflow
            if (MixTemp >= 0.)     // if < 0 (-9999.) then no inflows
                inflowRiverEnergy = t.RiverInflow * rhoW * HeatCapW * (MixTemp + convK);

            //energy balance without sensible heatflux
            RiverTotEnergy = RiverEnergyStor + RiverEnergyRad - latHeatflux + inflowRiverEnergy +
                             RiverEnergyPrecip + RiverEnergySurfRunoff + RiverEnergyGWRunoff;

            // water temperature river
            RiverVol = t.RiverStor + t.RiverInflow + RiverPrecip + RiverSurfRunoff + RiverGWRunoff -
                       RiverEvapo; // Evapo not accounted in WaterGAP currently (as of May 2020) -> different Vol.-balance than WaterGAP

            if (RiverVol > 0.) {      // for calculation safety only, RiverVol must not be below 0
                RiverWTemp = (RiverTotEnergy / (RiverVol * rhoW * HeatCapW)) - convK;        // [K] -> [degC]
                //energy balance with sensible heatflux
                RiverTotEnergy -= sensHeatflux;

                //explanation see gloLake
                if(RiverWTemp >= t.AirTempC){
                    RiverWTemp = (RiverTotEnergy / (RiverVol * rhoW * HeatCapW)) - convK;    // [K] -> [degC]
                    if(RiverWTemp < t.AirTempC)
                        RiverWTemp = t.AirTempC;
                }
                else {   //(RiverWTemp < t.AirTempC)
                    RiverWTemp = (RiverTotEnergy / (RiverVol * rhoW * HeatCapW)) - convK;    // [K] -> [degC]
                    if(RiverWTemp > t.AirTempC)
                        RiverWTemp = t.AirTempC;
                }
            }
            else
                RiverWTemp = -9999.;    // no volume in river
        }   //end no ice

        // water temperature calculation with power plants
        if(RiverWTemp >= 0. && t.powerPlantWW > 0. && t.calcTempPP) {
            if(iceThickness >= 5.e-6)
                RiverVol -= iceThickness * RiverArea;       // if ice adapt RiverVol

            if(RiverVol > 0.) {
                if (RiverVol < t.powerPlantWW)        // if RiverVol < withdrawal from power plant
                    withdrawal = RiverVol;      // set max. possible withdrawal
                else
                    withdrawal = t.powerPlantWW;

                if (withdrawal > t.powerPlantWC)                  // if adapted withdrawal < consumption of power plant ...
                    powerPlantReturnFlow = withdrawal - t.powerPlantWC;
                else
                    powerPlantReturnFlow = 0.;                    // ... returnFlow from power plant into river = 0.

                // if there is a returnFlow calculate mixing temperature btw. River and returnFlow (+ PPReturnFlowTempIncr [degC])
                if (powerPlantReturnFlow > 0.) {
                    RiverWTemp = (((RiverVol - withdrawal) * (RiverWTemp + convK) + powerPlantReturnFlow * (RiverWTemp + PPReturnFlowTempIncr + convK)) / (RiverVol - withdrawal + powerPlantReturnFlow)) - convK; // [degC]
                }
            }
        }
    }
    else {
        if(iceThickness > 0.)   // really ice there or no volume to begin with (not asked at beginning)
            RiverWTemp = 0.;    // max ice thickness = current water depth & no inflows --> completely frozen --> Temp = 0. degC
        else
            RiverWTemp = -9999.;    // no volume in river
    }

    if (RiverWTemp > 60.)     // highest ever measured air temp. 56 degC
        RiverWTemp = t.AirTempC;

    if(RiverWTemp < 1.e-10 && RiverWTemp != -9999.)     // to eliminate very small temperature values due to calculation errors e.g. 1.e-14 or negative temperatures
        RiverWTemp = 0.;

    if(iceThickness < 5.e-6)
        iceThickness = 0.;      // if ice thickness < 5 mm then ice thickness = 0. for next step

    returnvector[0] = RiverWTemp;
    returnvector[1] = iceThickness;

    return returnvector;
}

/**
* @brief calculates volume of SWB above thermocline
*
* implemented according to Wanders et al. (2019) doi:10.1029/2018WR023250
*
* @author Sebastian Ackermann
*
* @param swb surface area
* @param type of SWB: gloLake(0), reservoir(1)
*
* @return volume of SWB above thermocline [km3]
*/
double calcWaterTempClass::calcThermocline(const double &area, const short &swbType){
    double f = 0.;                  // fetch length

    //gloLake (assumed as rectangle -> longest length equals diagonal)
    if (swbType == 0)
        f = sqrt(2 * area);

    // reservoir (assumed as equilateral triangle -> longest length equals height of triangle)
    else if (swbType == 1)
        f = sqrt(3 * area / sqrt(3));

    return area * (9.52 * pow(f, 0.425)); // see Wanders et al. (2019) depth of thermocline: D_t = 9.52 * f^0.425
}