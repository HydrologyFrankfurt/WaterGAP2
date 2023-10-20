/***********************************************************************
*
*clcl is introduced by Kerstin Verzano (bevore marriage: Kerstin Schulze) in order to use climate
*data dependend classification of Köppen-climate zones (e.g. for defining arid-
*humid regions as well as for Köppen-based alpha coefficient). It is
*described in the Thesis of Kerstin: Verzano, K. (2009), Climate Change Impacts
*on Flood Related Hydrological Processes: Further Development and Application
*of a Global Scale Hydrological Model, 166 pp., University of Kassel.
*
***********************************************************************/
#include <cstdio>
#include <cmath>
//#include <iostream>
//#include <iomanip>
//#include <fstream>
//#include <vector>
#include <cstdlib>
#include <cstring>

//#include "clcl.h"
//#include "option.h"
#include "common.h"
#include "math.h"
//#include "climate.h"
//#include "geo.h"
#include "def.h"
#include "globals.h"

using namespace std;

//extern climateClass climate;
//extern geoClass geo;

clclClass::clclClass() {
}

void clclClass::init()
{
	for (int n = 0; n < ng; n++) {
		clcl_alpha[n] = 0.;
		cls[n] = 0;
	}
}


void clclClass::cleanup()
{
}

void clclClass::calcKoep()
{
	// koeppen climate classification
	// the following 28 classes are calculated, 30 can be calculated:
	// *A*
	// 11 ... Af
	// 12 ... Am
	// 13 ... Aw
	// *B*
	// 20 ... BW			// 21 ... BWk
							// 22 ... BWh
	// 25 ... BS			// 26 ... BSk
							// 27 ... BSh
	// *C*
	// 31 ... Cfa
	// 32 ... Cfb
	// 33 ... Cfc
	// 34 ... Csa
	// 35 ... Csb
	// 36 ... Csc
	// 37 ... Cwa
	// 38 ... Cwb
	// 39 ... Cwc
	// *D*
	// 41 ... Dfa
	// 42 ... Dfb
	// 43 ... Dfc
	// 44 ... Dfd
	// 45 ... Dsa
	// 46 ... Dsb
	// 47 ... Dsc
	// 48 ... Dsd
	// 49 ... Dwa
	// 50 ... Dwb
	// 51 ... Dwc
	// 52 ... Dwd
	// *E*
	// 61 ... EF
	// 62 ... ET
	// -------------------------------------------------------------------//

	// koeppen parameters

	Grid<float> t_min;     // coldest monthly temperature [°C]
	Grid<float> t_max;     // warmerst monthly temperature [°C]
	Grid<float> t_ann;     // mean annual temperature [°C]
	Grid<short>t_cnt;	 // number of month with temp >= 10°C
	Grid<float> p_min;	 // prec min value of all months [mm]
	Grid<float> p_max;	 // prec max value of all months  [mm]
	Grid<short>p_min_loc; // month of prec min value
	Grid<float> p_sum;     // sum yearly prec
	Grid<float> p_som;     // mean summer prec [mm]
	Grid<float> p_win;     // mean winter prec [mm]
	Grid<float> p_som_max; //	 summer max
	Grid<float> p_win_max; //	 summer min
	Grid<short> p_type; 	  // winter summer or constant rain
	
	// initialisation
	for (int n=0; n<ng; n++) {
		t_min[n]= 100.0;
		t_max[n]= -100.0;
		t_ann[n]= 0.0;
		t_cnt[n]= 0;
		p_min[n]= 10000.0;
		p_max[n]= 0.0;
		p_min_loc[n]= 0;
		p_sum[n]= 0.0;
		p_som[n]= 0.0;
		p_win[n]= 0.0;
		p_som_max[n]= 0.0;
		p_som_max[n]= 0.0;
		p_type[n]= -9999;
		cls[n]=0;
	}

	for (int n = 0; n < ng; ++n) 
	{
		for (int month = 0; month < 12; month++) {

			// temperature count: number of months with t >= 10°C
			if ((climate.G_temperature(n,month)/100.0) >= 10.0){
				t_cnt[n] += 1;}
			// determining temperature and precipitation minima and maxima
			if ((climate.G_temperature(n,month)/100.0) < t_min[n]){
				t_min[n] = (climate.G_temperature(n,month)/100.0);}
			if ((climate.G_temperature(n,month)/100.0) > t_max[n]){
				t_max[n] = climate.G_temperature(n,month)/100.0;}

			if (climate.G_precipitation(n,month) < p_min[n]){
				p_min[n] = climate.G_precipitation(n,month);
				p_min_loc[n]= month;
			}
			if (climate.G_precipitation(n,month) > p_max[n]){
				p_max[n] = climate.G_precipitation(n,month);}

			//calculating mean annual temp
			t_ann[n]+= (climate.G_temperature(n,month)/1200.0);  // temp in 0.01 °C, 12 monate

			// calculating annual and seasonal prec
			p_sum[n]+= climate.G_precipitation(n,month);

			// p som enthält jeweils sommer-summen,dh. noerdl hemisphere monate 4-9
			// suedl hemisphere monate 1-3+10-12; vice versa for winter
			// p_som_max enthält max wert des sommers (noerdl hem aus 4-9, suedl aus 1-3,10-12
			// p_win_max enthaelt max wert des winters (noerdl 1-3,10-12; suedl. 4-9)
			// northern hemisphere
			if (-(geo.G_row[n] / 2.0 - 90.25) > 0.0) {
				//wg2
				if (3<month&&month<10){
					p_som[n]+=climate.G_precipitation(n,month);
					if (climate.G_precipitation(n,month) > p_som_max[n]){
						p_som_max[n] = climate.G_precipitation(n,month);}
				}
				if ((1<=month&&month<4)||(9<month&&month<=12)){
					p_win[n]+=climate.G_precipitation(n,month);
					if (climate.G_precipitation(n,month) > p_win_max[n]){
						p_win_max[n] = climate.G_precipitation(n,month);}
				}
			}
			// southern hemisphere
			else {
				if (3<month&&month<10){
					p_win[n]+=climate.G_precipitation(n,month);
					if (climate.G_precipitation(n,month) > p_win_max[n]){
						p_win_max[n] = climate.G_precipitation(n,month);}
				}
				if ((1<=month&&month<4)||(9<month&&month<=12)){
					p_som[n]+=climate.G_precipitation(n,month);
					if (climate.G_precipitation(n,month) > p_som_max[n]){
						p_som_max[n] = climate.G_precipitation(n,month);}
				}
			}
		}

		//determining hemispherical summer-, winter- or constant rain
		//summer-rain: 2/3 of p_sum in summer -> 1
		//winter-rain: 2/3 of p_sum in winter -> 2
		//constant: neither summer nor winter -> 0
		if (p_win[n] >= (2./3.)*p_sum[n]){
			p_type[n] = 2;}   // winter rain
		else if (p_som[n] >= (2./3.)*p_sum[n]){
			p_type[n] = 1;}   // summer rain
		else{
			p_type[n] = 0;    // constant rain
		}
	}

	for (int n = 0; n < ng; ++n) 
	{

		// determinig the 5 basic climate classifications
		// A = 10 ... tropical
		// B = 20 ... dry
		// C = 30 ... warm
		// D = 40 ... snow
		// E = 60 ... ice

		if (t_max[n] < 10){ // hottest month lower than 10°C (1)
			cls[n] = 60;} // E = 60 ... polar
		else {
			if ((p_type[n] == 1) && ((p_sum[n]/10.0) < (2.0*t_ann[n]+28.0))){ // summer-rain and p_sum[n] in cm lower than arid limit
				cls[n] = 20; }// B = 20 ... dry
			else if ((p_type[n] == 2) && ((p_sum[n]/10.0) < (2.0*t_ann[n]))){ //winter-rain and p_sum[n] in cm lower than arid limit
				cls[n] = 20;} // B = 20 ... dry
			else if ((p_type[n] == 0) && ((p_sum[n]/10.0 ) < (2.*t_ann[n]+14.0))){ // constant-rain and p_sum[n] in cm lower than arid limit
				cls[n] = 20;} // B = 20 ... dry
			else if (t_min[n] >= 18.0){ // coldest month higher equals tropical limit
				cls[n] = 10;} // A = 10 ... tropical
			else if (t_min[n] <= -3.0){ // coldest month lower equals snow limit
				cls[n] = 40;} // D = 40 ... snow
			else if ((t_min[n] > -3.0) && (t_min[n] < 18.0)){ // -3°C < coldest month < 18°C
				cls[n] = 30;} // C = 30 ... warm
			else {
				cls[n] = 99;
				cerr << "climate characteristic of grid cell " << n << " does not fit in any of the Koeppen classes" << endl;
			}
		}
	}

	// second letters and third koeppen letters:
	for (int n = 0; n < ng; ++n) 
	{

		//  splitting polar climate E into frost EF and tundra ET
		//  EF = 61 ... polar frost, no vegetation
		//  ET = 62 ... polar tundra, dwarf tree species and mosses

		if (cls[n] == 60 && t_max[n] < 0.0){
			cls[n] = 61;}  //EF ... polar frost -> finished
		else if (cls[n] == 60 && t_max[n] >= 0.0){
			cls[n] = 62;}  //ET ... polar tundra -> finished

		// splitting dry climate B into desert BW and steppe BS
		// BW = 20 ... desert
		// BS = 25 ... steppe, vegetation bushes to grassland

		if (cls[n] == 20){
			if (p_type[n] == 1 && p_sum[n]/10.0 > t_ann[n]+14.0){
				cls[n] = 25;} // BS ... steppe
			else if (p_type[n] == 2 && p_sum[n]/10.0 > t_ann[n]){
				cls[n] = 25;} // BS ... steppe
			else if (p_type[n] == 0 && p_sum[n]/10.0 > t_ann[n]+7.0){
				cls[n] = 25;} // BS ... steppe
			else{
				cls[n] = 20;} // BW ... desert
		}

		// splitting the tropical climate A based on precipitation amounts
		// Af = 11 ... tropical wet
		// Am = 12 ... tropical monsoonal, Af and Am both tropical evergreen rainforest
		// Aw = 13 ... tropical savanna (summer-dry/winter-dry not distinguished)
		
		//(5)
		if (cls[n] == 10 && p_min[n] >= 60.0){ // no dry season 
			cls[n] = 11;} // Af ... tropical wet -> finished
		//(6)
		else if (cls[n] == 10 && p_min[n] < 60.0 && p_sum[n] >= 25.0*(100.0-p_min[n])){
			cls[n] = 12;} // Am ... tropical monsoonal; dry season compensted by annaul prec -> finished
		//(7)
		else if (cls[n]==10){
			cls[n]= 13;} // Aw ... tropical savanna  -> finished
		

		 // splitting temperate climate C and boreal forest/snow climate D based on precipitation
		 // Cf/Df = 31/41 ... subtropical/continental without dry period
		 // Cs/Ds = 34/44 ... mediterranean/continental with dry summer
		 // Cw/Dw = 37/47 ... subtropical/continental with dry winter

		 if (cls[n] == 30 || cls[n] == 40){
			if (-(geo.G_row[n] / 2.0 - 90.25) > 0.0) {
				//wg2
				// northern hemisphere summer-dry
				if (p_min_loc[n] >= 4 && p_min_loc[n] <= 9){
					if (p_min[n] < p_win_max[n]/3.0 && p_min[n] < 40.0){ // summer minimum lower than 1/3 of winter maximum and lower than 40 mm/d
						cls[n] = cls[n] + 4;} // Cs/Ds ... mediterranean/continental with dry summer
					else{
						cls[n] = cls[n]+1;}// Cf/Df ... subtropical/continental without dry period
				}
				// northern hemisphere winter-dry
				else {
					if (p_min[n] < p_som_max[n]/10.0 ){ // winter minimum lower than 1/10 of northern summer maximum
						cls[n] = cls[n]+7;} // Cw/Dw ... subtropical/continental with dry winter
					else {
						cls[n] = cls[n]+1;} // Cf/Df ... subtropical/continental without dry period
				}
			}

			else {
				// southern hemisphere
				// southern hemisphere winter-dry
				if (p_min_loc[n] >= 4 && p_min_loc[n] <= 9){
					if (p_min[n] < p_som_max[n]/10.0 ){ // winter minimum lower than 1/10 of southern summer maximum
						cls[n] = cls[n]+7;} // Cw/Dw ... subtropical/continental with dry winter
					else{
						cls[n] = cls[n]+1;} // Cf/Df ... subtropical/continental without dry period
				}
				// southern hemisphere summer-dry
				else {
					if (p_min[n] < p_win_max[n]/3.0 && p_min[n] < 40.0){ // summer minimum lower than 1/3 of southern winter maximum and lower than 40 mm/d
						cls[n] = cls[n]+4;} // Cs/Ds ... mediterranean/continental with dry summer
					else {
						cls[n] = cls[n]+1;} // Cf/Df ... subtropical/continental without dry period

				}
			}

		 }

		 // third letter: splitting warm climate C and snow climate D based on temperature
		 // Cfa/Dfa = 31/41 ... hot summer
		 // Csa/Dsa = 34/44 ... hot summer
		 // Cwa/Dwa = 37/47 ... hot summer
		 // Cfb/Dfb = 32/42 ... warm summer
		 // Csb/Dsb = 35/45 ... warm summer
		 // Cwb/Dwb = 38/48 ... warm summer
		 // Cfc/Dfc = 33/43 ... cool summer
		 // Csc/Dsc = 36/46 ... cool summer
		 // Cwc/Dwc = 39/49 ... cool summer
		 // Dfd,Dsd,Dwd = 51,54,57 ... extrem continental

		 if (cls[n] > 30 && cls[n] < 50){
			 if (t_max[n] < 22.0){
				 if (t_cnt[n] >= 4){
					 cls[n] = cls[n]+1;} // C(fsw)b/D(fsw)b -> finished
				 else {
					 if (t_min[n] > -38.0){
						 cls[n] = cls[n]+2;} // C(fsw)c/D(fsw)c -> finished
					 else{
						cls[n] = cls[n]+10;} // D(fsw)d -> finished
				}
			 }
		 }

		  // sort snow climate D by descending letters
		  // Dfa = 41
		  // Dfb = 42
		  // Dfc = 43
		  // Dfd = 44
		  // Dsa = 45
		  // Dsb = 46
		  // Dsc = 47
		  // Dsd = 48
		  // Dwa = 49
		  // Dwb = 50
		  // Dwc = 51
		  // Dwd = 52

		  if (cls[n] > 40 && cls[n] < 60)
		  	cls[n] = cls[n]+100;

		  if (cls[n] == 141) cls[n] = 41;
		  if (cls[n] == 142) cls[n] = 42;
		  if (cls[n] == 143) cls[n] = 43;
		  if (cls[n] == 151) cls[n] = 44;
		  if (cls[n] == 144) cls[n] = 45;
		  if (cls[n] == 145) cls[n] = 46;
		  if (cls[n] == 146) cls[n] = 47;
		  if (cls[n] == 154) cls[n] = 48;
		  if (cls[n] == 147) cls[n] = 49;
		  if (cls[n] == 148) cls[n] = 50;
		  if (cls[n] == 149) cls[n] = 51;
		  if (cls[n] == 157) cls[n] = 52;
	}

	cout << "***** climate classification acc. to Koeppen has been carried out *****" << endl;
}

// assign alpha factor for PT method in PET calculations
void clclClass::alphaKoep()
{
	// assign alpha factor for different cls-classes
	for (int n = 0; n < ng; ++n) {
		switch (cls[n]){
			//based on Class A Pan evaporation measurements
			case 11:
			case 20:
				clcl_alpha[n] = 2.0;
				break;
			case 12:
				clcl_alpha[n] = 1.65;
				break;
			case 13:
				clcl_alpha[n] = 1.4;
				break;
			case 25:
				clcl_alpha[n] = 1.85;
				break;
			case 31:
			case 32:
				clcl_alpha[n] = 1.45;
				break;
			case 34:
			case 35:
				clcl_alpha[n] = 1.65;
				break;
			case 37:
				clcl_alpha[n] = 1.37;
				break;
			case 36:
			case 38:			//Cwb
			case 39:			//[Cwc]
			case 49:			//Dwa
				clcl_alpha[n] = 1.26;
				break;
			case 41: 			//Dfa
			case 42: 			//Dfb
			case 50:			//Dwb
				clcl_alpha[n] = 1.3;
				break;
			case 45:
			case 33: 			//[Csc]
				clcl_alpha[n] = 1.2;
				break;
			case 43:
			case 44:
			case 46:
			case 47:
			case 51:
			case 62:
				clcl_alpha[n] = 1.1;
				break;
			case 48:
			case 52:
				clcl_alpha[n] = 1.15;
				break;
			case 61:
				clcl_alpha[n] = 1.0;
				break;
		}
	}
	
}



