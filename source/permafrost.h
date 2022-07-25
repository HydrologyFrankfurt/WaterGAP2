#if !defined (_permafrost_h_)
#define _permafrost_h_


#include <vector>

#include "def.h"
#include "gridio.h"
using namespace std;

class permaClass {
	public:
	permaClass();
	void calcFrostNumber(const short year, const int n);
	float calcgwFactor(const int n, float gwFactor_old, float gw_permafactor_old);
		
	//float (*G_monthlyTempPerma)[12];
	
	void writeYearlyGrids(char *outputDir, const int year);

	private:
	float frost_number_air[ng];       // air frost number [-]
	float frost_number_surface[ng];   // surface frost number [-]

};


#endif

