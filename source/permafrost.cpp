/**********************************************************************
*
*permafrost calculations developed by Tim aus der Beek: Aus der Beek, T.;
*Teichert, E. (2008) Global Permafrost Distribution in the Past, Present and
*Future Proc., 9th International Conference on Permafrost, Fairbanks (Alaska),
*29.6.-3.7., adapted to work with daily climate input by Stephanie Eisner
*
***********************************************************************/
//#include <cstdio>
#include <cmath>
//#include <iostream>
//#include <iomanip>
//#include <fstream>
//#include <vector>
#include <cstdlib>
#include <cstring>
#include "common.h"
#include "globals.h"
#include "math.h"
#include "grid.h"

using namespace std;

//extern optionClass options;
//extern climateClass climate;
//extern groundwaterFactorClass GW;

permaClass::permaClass() {
}

void permaClass::calcFrostNumber(const short year, const int n)
{
	const double pi = 3.141592653589793;
	// parameters for permafrost  (Tim: sind wirklich alle Parameter notwendig hier??? -> lokale Variablen??

	float temp_annualmean = 0.0; // annual mean air temperature [°C]
	float temp_annualmean_surface; // annual mean surface temperature [°C]
	float temp_annualhot = -100000.0; // air temperature of the years hottest month [°C]
	float temp_annualcold = 100000.0; // air temperature of the years coldest month [°C]
	float temp_annualamplitude = 0.0; // intra annual air temperature amplitude [°C]
	float temp_annualamplitude_surface; // intra annual surface temperature amplitude [°C]
	float temp_summermean = 0.0; // mean summer air temperature [°C]
	float temp_wintermean = 0.0; // mean winter air temperature [°C]
	float temp_wintermean_surface; // mean winter surface temperature [°C]
	float frostangle = 0.0; // point @ x-axis where temp curve crosses 0 [-]
	float length_summer = 0.0; // length of summer [days]
	float length_winter = 0.0; // length of winter [days]
	float DDFreeze = 0.0; // air freezing index [-]
	float DDFreeze_surface; // surface freezing index [-]
	float DDThaw = 0.0; // air thawing index [-]
	float snow_depth = 0.0; // mean snow depth [mm] snow depth immer wieder auf '0'??? Kein Übertrag vom vorhergehenden Jahr???
	float snow_dampingdepth = 0.0; // damping depth of snow [m]
	float snow_density = 0.0; // snow density [kg/m³]
	float snow_conductivity = 0.0; // thermal conductivity of snow [w/m°C]
	float snow_capacity = 0.0; // specific heat capacity of snow [J/kg°C]
	float snow_diffusivity = 0.0; // thermal diffusivity of snow [cm²/s]
	int temp_length; // length of annual temperature cycle [days]
	float latitude = 0;
	short k = 0;

	// initialisation of the current cell values
	frost_number_air[n] = -9999.9;
	frost_number_surface[n] = -9999.9;

	// calculate monthly values 
	for (short month = 0; month < 12; month++) {

		// check for years hottest month
		if (climate.G_temperature(n,month)/100.0 > temp_annualhot)
			temp_annualhot = climate.G_temperature(n,month)/100.0;
		// check for years coldest month
		if (climate.G_temperature(n,month)/100.0 < temp_annualcold)
			temp_annualcold = climate.G_temperature(n,month)/100.0;
		
	}

	temp_annualmean = (temp_annualhot + temp_annualcold)/2.0;
	temp_annualamplitude = (temp_annualhot - temp_annualcold)/2.0;
	frostangle = -1.0 * temp_annualmean/ temp_annualamplitude;
		
	if (frostangle < -0.99) {frostangle = -0.99;}			// frostangle < -1 causes NANs!!!
	if (frostangle > 0.99) {frostangle = 0.99;} 	   // frostangle > 1 causes NANs!!!
	frostangle = acos(frostangle);
	
	temp_summermean = sin(frostangle)/frostangle;
	temp_summermean = temp_annualmean + (temp_annualamplitude * temp_summermean);
	temp_wintermean = sin(frostangle)/(pi - frostangle);
	temp_wintermean = temp_annualmean - (temp_annualamplitude * temp_wintermean);
	length_summer = 365.0 * (frostangle/pi);
	length_winter = 365.0 - length_summer;
	
	DDThaw = abs(temp_summermean) * length_summer;
	DDFreeze = abs(temp_wintermean) * length_winter;

	// check this Tim: sqrt of 0 is possible !!
	if (DDFreeze == 0.0)
	{
		// sqrt of 0 not possible
		DDFreeze = 0.01;
	}
	if (temp_annualcold > 0.0) 
	{
		// no permafrost anyway
		DDFreeze = 100.0;
		DDThaw = 9000.0;
	}
	if (temp_annualhot < 0.0) 
	{
		// no thawing anyway
		DDFreeze = 9000.0;
		DDThaw = 100.0;
	}
	
	frost_number_air[n]	= sqrt(DDFreeze) / (sqrt(DDFreeze) + sqrt(DDThaw)); 

	//	SURFACE Frost Number
	//extern geoClass geo;
		
	// use 'G_row' here as value not as index in an array
	latitude = -(geo.G_row[n] / 2.0 - 90.25);	
	if (temp_annualmean <= 4.142)
		snow_density =  545.0 * std::pow((5.0 - temp_annualmean), -1.15) + 50.0; // ref:Meister 1985
	else 
		snow_density = 700.0;
	
	for (int month = 0; month < 12; month++){ 
		if (climate.G_temperature(n,month)/100.0 <= 0.0)
			k++;	// count month with temp < 0
	}

	int j=0;

	for (int month = 0; month < 12; month++){	
		if (climate.G_temperature(n,month) <= 0.0)  {
			j++;
			snow_depth += (climate.G_precipitation(n,month)/(snow_density/1000.0)) *(k - (j - 1));	// precip as snow in cold months
		}
	}
	
	if (k == 0) 
		snow_depth= 0.0;
	else 
		snow_depth = (0.5*(1-cos(2.0*latitude*pi/180.0)))*(snow_depth/k);

	if (snow_density > 156.0)
		snow_conductivity = 0.138 - (1.01 * snow_density/1000) + (3.233 * (snow_density/1000)*(snow_density/1000));
	else 
		snow_conductivity = 0.023 + 0.234*(snow_density/1000);	 //Ref: Sturm 1997

	snow_capacity = 7.79 * temp_wintermean + 2115.0;
	snow_diffusivity = 1000000 * snow_conductivity / (snow_capacity * snow_density); // small values
	snow_dampingdepth = sqrt(snow_diffusivity/1000000 * 365.0 *86400.0/ pi);
	temp_annualamplitude_surface = temp_annualamplitude * exp(-snow_depth/snow_dampingdepth/1000.0);
	temp_wintermean_surface = temp_annualmean - temp_annualamplitude_surface * (sin(frostangle)/(pi - frostangle));
	DDFreeze_surface = temp_wintermean_surface * length_winter;
	temp_annualmean_surface = (DDThaw - DDFreeze_surface) / 365.0;
	
	if (DDFreeze_surface == 0.0)
	{
		// sqrt of 0 not possible
		DDFreeze_surface = 0.01;
	}
	if (temp_annualcold > 0.0) 
	{
		// no permafrost anyway
		DDFreeze_surface = 100.0;
		DDThaw = 9000.0;
	}
	if (temp_annualhot < 0.0) 
	{
		// no thawing anyway
		DDFreeze_surface = 9000.0;
		DDThaw = 100.0;
	}

	frost_number_surface[n] = sqrt(abs(DDFreeze_surface)) / (sqrt(abs(DDFreeze_surface)) + sqrt(abs(DDThaw)));
}

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
		
}
 
void permaClass::writeYearlyGrids(const std::string outputDirectory, const int year)
{
	if (2 != options.grid_store && 4 != options.grid_store && 6 != options.grid_store) return;

	Grid<> GW_factor;

	if (options.outAirFrost) {
		frost_number_air.write(outputDirectory + "/G_AIRFROST_" + to_string(year) + ".UNF0");
	}
	if (options.outSurfaceFrost) {	
		frost_number_surface.write(outputDirectory + "/G_SURFACEFROST_" + to_string(year) + ".UNF0");
	}

	if (options.outGWFactor) {
        for (int i=0; i<ng; ++i)
			GW_factor[i]=GW.getgwFactor(i);
		GW_factor.write(outputDirectory + "/G_GW_FACTOR_" + to_string(year) + ".UNF0");
	}
}



