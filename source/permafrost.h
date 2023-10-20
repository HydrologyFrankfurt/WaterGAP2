#if !defined (_permafrost_h_)
#define _permafrost_h_

#pragma once
#include <vector>

#include "def.h"
#include "grid.h"

using namespace std;

class permaClass {
	public:
		permaClass();
		void calcFrostNumber(const short year, const int n);
		float calcgwFactor(const int n, float gwFactor_old, float gw_permafactor_old);
		void writeYearlyGrids(const std::string outputDir, const int year);

	private:
		Grid<float> frost_number_air;       // air frost number [-]
		Grid<float> frost_number_surface;   // surface frost number [-]
 };
#endif
