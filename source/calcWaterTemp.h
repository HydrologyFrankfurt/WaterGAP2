#if !defined (_calcWaterTemp_h_)
#define _calcWaterTemp_h_
#pragma once

#include <vector>

using namespace std;

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

class calcWaterTempClass {
public:
    struct transfer {           // structure to transfer all values needed from routing.cpp to calcWaterTemp.cpp ( easily expandable)
        // power plants
        bool calcTempPP;        // true if power plants must be considered in water temperature calculation
        double powerPlantWW;    // water withdrawal from river to cool power plant
        double powerPlantWC;    // water consumption from power plant due to cooling

        //misc
        double ShortRad;        // shortwave radiation
        double LongRad;         // longwave radiation
        double AirTempC;        // air temperature
        double CellCorrFact;    // cell correction factor
        double openWaterPET;    // open water potential evaporation
        double openWaterPrecip; // open water precipitation
        double CellArea;        // area of cell
        double GWTemp;          // groundwater temperature
        vector<double> cellInflows;         // water volume flowing into cell from surrounding cells
        vector<double> cellInflowsTemp;     // water temperature of water flowing into cell from surrounding cells
        vector<double> iceThickness{0.,0.,0.,0.,0.,0.};     // ice thickness of swb in cell [0] river, [1] locLake, [2] locWetland, [3] gloLake, [4] reservoir, [5] gloWetland

        //river
        double RiverStor;               // volume stored in river
        double localGWRunoffIntoRiver;  // groundwater runoff into river
        double localRunoffIntoRiver;    // surface water runoff into river
        double RiverBottomWidth;        // river width at bottom of river
        double RiverLength;             // length of river
        double RiverLengthCheck;        // length of river to check if river exists (!=55)
        double RiverInflow;             // inflow into river
        bool hasOutflow;                // true if cell flows into other cell

        //locLake
        double locLakeArea;             // surface area of locLake
        double locLakeStorPrevStep;     // stored volume in locLake at prev. step
        double locLakeTempPrevStep;     // water temperature of locLake at prev. step
        double dailyLocLakeEvapo;       // locLake evaporation
        double locLakeAreaRedFactor;    // area reduction factor for locLake
        double localRunoff;             // surface runoff into locLake
        double localGWRunoff;           // groundwater runoff into locLake

        //locWetland
        double locWetlandArea;          // surface area of locWetland
        double locWetlandStorPrevStep;  // stored volumne in locWetland at prev. step
        double locWetlandTempPrevStep;  // water temperature of locWetland at prev. step
        double dailyLocWetlandEvapo;    // locWetland evaporation
        double locWetlandAreaRedFactor; // area reduction factor of locWetland
        double inflowLocWetland;        // inflow into locWetland from locLake
        double outflowLocWetland;       // outflow of locWetland

        //gloLake
        double gloLakeArea;             // surface area of gloLake
        double gloLakeStorPrevStep;     // stored volume in gloLake at prev. step
        double gloLakeTempPrevStep;     // water temperature of gloLake at prev. step
        double dailyGloLakeEvapo;       // gloLake evaporation
        double inflowGloLake;           // inflows into glolake from other cells and swb

        //reservoir
        double inflowReservoir;         // inflow into reservoir from other swb
        double reservoirArea;           // surface area of reservoir
        double reservoirStorPrevStep;   // stored volume in reservoir at prev. step
        double reservoirTempPrevStep;   // water temperature of reservoir at prev. step
        double dailyReservoirEvapo;     // reservoir evaporation

        //gloWetland
        double gloWetlandInflow;                // inflow into gloWetland from other swb
        double gloWetlandArea;                  // surface area of gloWetland
        double gloWetlandStorPrevStep;          // stored volume in gloWetland at prev. step
        double gloWetlandTempPrevStep;          // water temperature of gloWetland at prev. step
        double dailyGloWetlandEvapo;            // gloWetland evaporation
        double gloWetlandAreaReductionFactor;   // area reduction factor of gloWetland

    };

    transfer tsf;       // instance of structure transfer

    vector<double> calcWTemp(const transfer &t);

    vector<double> calcLocLakeTemp(const transfer &t, const double &NetRadIce, const double &NetRadWater, const double &latHeatVapor);

    vector<double> calcLocWetlandTemp(const transfer &t, const double &NetRadIce, const double &NetRadWater, const double &latHeatVapor, const double &locLakeTemp, const bool &gw_used);

    vector<double> calcGloLakeTemp(const transfer &t, const double &NetRadIce, const double &NetRadWater, const double &latHeatVapor, const double &MixTemp, const bool &gw_used);

    vector<double> calcReservoirTemp(const transfer &t, const double &NetRadIce, const double &NetRadWater, const double &latHeatVapor, const double &MixTemp, const bool &gw_used);

    vector<double> calcGloWetlandTemp(const transfer &t, const double &NetRadIce, const double &NetRadWater, const double &latHeatVapor, const double &MixTemp, const bool &gw_used);

    vector<double> calcRiverTemp(const transfer &t, const double &NetRadIce, const double &NetRadWater, const double &latHeatVapor, const double &MixTemp, const bool &gw_used);

    double calcThermocline (const double &area, const short &swbType);
};
#endif
