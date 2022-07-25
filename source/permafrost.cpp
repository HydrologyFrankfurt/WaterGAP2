/**********************************************************************
*
*permafrost calculations developed by Tim aus der Beek: Aus der Beek, T.;
*Teichert, E. (2008) Global Permafrost Distribution in the Past, Present and
*Future Proc., 9th International Conference on Permafrost, Fairbanks (Alaska),
*29.6.-3.7., adapted to work with daily climate input by Stephanie Eisner
*
***********************************************************************/
#include <cstdio>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <cstring>

#include "gridio.h"
#include "permafrost.h"
#include "daily.h"
#include "timestring.h"
#include "option.h"
#include "climate.h"
#include "geo.h"
#include "math.h"
#include "gw_frac.h"

using namespace std;


extern gridioClass gridIO;
extern optionClass options;
extern climateClass climate;
extern groundwaterFactorClass GW;

permaClass::permaClass() {
}

void permaClass::calcFrostNumber(const short year, const int n)
{
	const double pi = 3.141592653589793;
	// parameters for permafrost  (Tim: sind wirklich alle Parameter notwendig hier??? -> lokale Variablen??
	//short temp_amplitude_minmax[ng][12]; // monthly temperature min/max amplitude [°C]
	float temp_annualmean[ng]; // annual mean air temperature [°C]
	float temp_annualmean_surface; // annual mean surface temperature [°C]
	float temp_annualhot; // air temperature of the years hottest month [°C]
	float temp_annualcold; // air temperature of the years coldest month [°C]
	float temp_annualamplitude; // intra annual air temperature amplitude [°C]
	float temp_annualamplitude_surface; // intra annual surface temperature amplitude [°C]
	float temp_summermean; // mean summer air temperature [°C]
	float temp_wintermean[ng]; // mean winter air temperature [°C]
	float temp_wintermean_surface; // mean winter surface temperature [°C]
	float frostangle; // point @ x-axis where temp curve crosses 0 [-]
	float length_summer; // length of summer [days]
	float length_winter; // length of winter [days]
	float DDFreeze; // air freezing index [-]
	float DDFreeze_surface; // surface freezing index [-]
	float DDThaw; // air thawing index [-]
	float snow_depth[ng]; // mean snow depth [mm]
	float snow_dampingdepth[ng]; // damping depth of snow [m]
	float snow_density[ng]; // snow density [kg/m³]
	float snow_conductivity[ng]; // thermal conductivity of snow [w/m°C]
	float snow_capacity[ng]; // specific heat capacity of snow [J/kg°C]
	float snow_diffusivity[ng]; // thermal diffusivity of snow [cm²/s]
	int temp_length; // length of annual temperature cycle [days]
	float latitude; 
	short k;

	// initialisation	
	//for (short month = 0; month < 12; month++) temp_amplitude_minmax[n][month] = 0;
	temp_annualmean[n] = 0.0;
	temp_wintermean[n] = 0.0;
	frost_number_air[n] = -9999.9;
	snow_depth[n] = 0.0;	// snow depth immer wieder auf '0'??? Kein Übertrag vom vorhergehenden Jahr???
	snow_density[n] = 0.0;
	snow_conductivity[n] = 0.0;
	snow_capacity[n] = 0.0;
	snow_diffusivity[n] = 0.0;
	snow_dampingdepth[n] = 0.0;
	frost_number_surface[n] = -9999.9;
	temp_annualhot 		= -100000.0;
	temp_annualcold 	= 100000.0;
	temp_annualamplitude= 0.0;
	frostangle 			= 0.0;
	temp_summermean 	= 0.0;
	length_summer 		= 0.0;
	length_winter 		= 0.0;
	DDThaw 				= 0.0;
	DDFreeze 			= 0.0;	
	latitude 			= 0;
	k 					= 0;

	// calculate monthly values 
	for (short month = 0; month < 12; month++) {

		// check for years hottest month
		//temp_amplitude_minmax[n][month] = climate.G_tmax[n][month] - climate.G_tmin[n][month];
		if (climate.G_temperature[n][month]/100.0 > temp_annualhot)
			temp_annualhot = climate.G_temperature[n][month]/100.0;
		// check for years coldest month
		//temp_amplitude_minmax[n][month] = abs(temp_amplitude_minmax[n][month]);
		if (climate.G_temperature[n][month]/100.0 < temp_annualcold)
			temp_annualcold = climate.G_temperature[n][month]/100.0;
		
	}

	temp_annualmean[n] = (temp_annualhot + temp_annualcold)/2.0;
	temp_annualamplitude = (temp_annualhot - temp_annualcold)/2.0;
	frostangle = -1.0 * temp_annualmean[n]/ temp_annualamplitude;
		
	if (frostangle < -0.99) {frostangle = -0.99;}			// frostangle < -1 causes NANs!!!
	if (frostangle > 0.99) {frostangle = 0.99;} 	   // frostangle > 1 causes NANs!!!
	frostangle = acos(frostangle);
	
	temp_summermean = sin(frostangle)/frostangle;
	temp_summermean = temp_annualmean[n] + (temp_annualamplitude * temp_summermean);
	temp_wintermean[n] = sin(frostangle)/(pi - frostangle);
	temp_wintermean[n] = temp_annualmean[n] - (temp_annualamplitude * temp_wintermean[n]);
	length_summer = 365.0 * (frostangle/pi);
	length_winter = 365.0 - length_summer;
	
	DDThaw = abs(temp_summermean) * length_summer;
	DDFreeze = abs(temp_wintermean[n]) * length_winter;
	// check this Tim: sqrt of 0 is possible !!
	if (DDFreeze == 0.0)
		{DDFreeze = 0.01;} 	// sqrt of 0 not possible
	if (temp_annualcold > 0.0) 
		{DDFreeze = 100.0; DDThaw = 9000.0;} // no permafrost anyway
	if (temp_annualhot < 0.0) 
		{DDFreeze = 9000.0; DDThaw = 100.0;} // no thawing anyway
	
	frost_number_air[n]	= sqrt(DDFreeze) / (sqrt(DDFreeze) + sqrt(DDThaw)); 

	//	SURFACE Frost Number
	extern geoClass geo;
		
	// use 'G_row' here as value not as index in an array
	//latitude = -((geo.G_row[n]-1)/2-90);	
	latitude = -(geo.G_row[n] / 2.0 - 90.25);	
	// snow_density[n] =  100.0;
	if (temp_annualmean[n] <= 4.142)
		snow_density[n] =  545.0 * std::pow((5.0 - temp_annualmean[n]), -1.15) + 50.0; // ref:Meister 1985
	else 
		snow_density[n] = 700.0;
	
	for (int month = 0; month < 12; month++){ 
		if (climate.G_temperature[n][month]/100.0 <= 0.0) 
			k++;	// count month with temp < 0
	}

	int j=0;

	for (int month = 0; month < 12; month++){	
		if (climate.G_temperature[n][month] <= 0.0)  {
			j++;
			// if (n==7651) cout <<snow_depth[n]<< "\t"<<"month:\t" << month<< "\t" << (climate.G_precipitation[n][month]/(snow_density[n]/1000.0)) *(k - (j - 1))<< endl;
			snow_depth[n] += (climate.G_precipitation[n][month]/(snow_density[n]/1000.0)) *(k - (j - 1));	// precip as snow in cold months
		}
		//else {snow_depth[n] += 0.01;}
	}  // for (month)
	
	if (k == 0) 
		snow_depth[n]= 0.0;
	else 
		snow_depth[n] = (0.5*(1-cos(2.0*latitude*pi/180.0)))*(snow_depth[n]/k);

	// if (n==7651) cout << "snow_depth[n]: " << snow_depth[n] << endl;

	// snow_conductivity[n] = 2.1 * 0.01 + 4.2 * 0.0001 * snow_density[n] + 2.2 * 0.000000001 * pow(snow_density[n],3);	// Ref: Paterson 1981
	if (snow_density[n] > 156.0)
		snow_conductivity[n] = 0.138 - (1.01 * snow_density[n]/1000) + (3.233 * (snow_density[n]/1000)*(snow_density[n]/1000));	
	else 
		snow_conductivity[n]= 0.023 + 0.234*(snow_density[n]/1000);	 //Ref: Sturm 1997
	snow_capacity[n] = 7.79 * temp_wintermean[n] + 2115.0;	
	snow_diffusivity[n] = 1000000 * snow_conductivity[n] / (snow_capacity[n] * snow_density[n]); // small values 
	snow_dampingdepth[n] = sqrt(snow_diffusivity[n]/1000000 * 365.0 *86400.0/ pi);
	temp_annualamplitude_surface = temp_annualamplitude * exp(-snow_depth[n]/snow_dampingdepth[n]/1000.0);
	temp_wintermean_surface = temp_annualmean[n] - temp_annualamplitude_surface * (sin(frostangle)/(pi - frostangle));
	DDFreeze_surface = temp_wintermean_surface * length_winter;
	temp_annualmean_surface = (DDThaw - DDFreeze_surface) / 365.0;
	
	if (DDFreeze_surface == 0.0)
			{DDFreeze_surface = 0.01;} 	// sqrt of 0 not possible
	if (temp_annualcold > 0.0) 
			{DDFreeze_surface = 100.0; DDThaw = 9000.0;} // no permafrost anyway
	if (temp_annualhot < 0.0) 
			{DDFreeze_surface = 9000.0; DDThaw = 100.0;} // no thawing anyway
	frost_number_surface[n] = sqrt(abs(DDFreeze_surface)) / (sqrt(abs(DDFreeze_surface)) + sqrt(abs(DDThaw)));
		
	//if (n==5) cout << "frost_number_surface[n]: " << frost_number_surface[n]<<endl 
	//		<< "Density:meister Conductivity:sturm Angle: acos " << endl;
	
} // end of calcFrostNumber()

float permaClass::calcgwFactor(const int n, float gwFactor_old, float gw_permafactor_old)
{
	float gwFactor_new, gw_permafactor_new;
	
	if(gw_permafactor_old==0.)
		gwFactor_new = gwFactor_old;
	else
		gwFactor_new = gwFactor_old / gw_permafactor_old;
	
	// continuous permafrost
	if ( frost_number_surface[n] >= 0.67 )
		gw_permafactor_new = 0.05;
	// extensive discontinuous permafrost
	else if ( (frost_number_surface[n] < 0.67) && (frost_number_surface[n] >= 0.6) )
		gw_permafactor_new = 0.3;
	// sporadic discontinuous permafrost
	else if ( (frost_number_surface[n] < 0.6) && (frost_number_surface[n] >= 0.5) )
		gw_permafactor_new = 0.7;
	// ~ isolated patches of permafrost
	else if ( (frost_number_surface[n] < 0.5) && (frost_number_surface[n] >= 0.45) )
		gw_permafactor_new = 0.95;
	// no permafrost at all
	else
		gw_permafactor_new = 1.;
	
	GW.setpermaFactor(n,gw_permafactor_new);
	
	return ( gwFactor_new * gw_permafactor_new );	
		
} // end of calcgwFactor()
 
void permaClass::writeYearlyGrids(char *outputDir, const int year)
{
	if (2 != options.grid_store && 4 != options.grid_store && 6 != options.grid_store) return;

	char filename[250];
        //double *temp = new double[ng_land];
        double *temp = new double[ng];

	if (options.outAirFrost) {
		sprintf(filename, "%s/G_AIRFROST_%d.UNF0", outputDir, year);
		gridIO.writeUnfFile(filename, ng, frost_number_air);
	}
	if (options.outSurfaceFrost) {	
		sprintf(filename, "%s/G_SURFACEFROST_%d.UNF0", outputDir, year);
		gridIO.writeUnfFile(filename, ng, frost_number_surface);
	}


	if (options.outGWFactor) { // new output options 3.1
                //for (int i=0; i<ng_land; ++i) temp[i]=GW.getgwFactor(i);
            for (int i=0; i<ng; ++i) temp[i]=GW.getgwFactor(i);
		sprintf(filename, "%s/G_GW_FACTOR_%d.UNF0", outputDir, year);
                //gridIO.writeUnfFile(filename, ng_land, temp);
                gridIO.writeUnfFile(filename, ng, temp);
	}

	delete[] temp; temp = NULL;
}



