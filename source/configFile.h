#if !defined (_configfile_h_)
#define _configfile_h_
#include "common.h"
class ConfigFile{
   
  public:
    ConfigFile(std::string, int year, int month, std::string); // OE 01-2018: additional input int year and  int month , std::string
    ~ConfigFile();
		
		std::string configfile;
		std::string progName; //OE
		std::string startvaluefile; 
		std::string calibrationfile;
		std::string parameterfile;
		std::string snowInElevationfile; // snow in elevation
		std::string additionalfile; // additional input: G_days_since_start, G_growingStatus and G_PrecSum
		
		std::string outputmeanfile;
		std::string outputlastdayfile;
		std::string outputdailyfile;
 		std::string outputsnowlastdayfile; // snow in elevation
 		std::string outputadditionalfile; // additional output: G_days_since_start, G_growingStatus ,G_PrecSum
		std::string outputparameter; //OE
    
		std::string runtimeoptionsfile;
		std::string outputoptionsfile;
		std::string routingoptionsfile;
		std::string stationsfile;		
		
		std::string inputDir;
		std::string outputDir;
        std::string dischargeDir;
        std::string snowCoverFracDir;
	std::string reanPrecipDir;
		std::string climateDir;
		std::string glacierDir;
		std::string climateSLRDDir;
		std::string routingDir;
		std::string wateruseDir;
		std::string landcoverDir;
		
		int startMonth;
		int startYear;
		int endMonth;
		int endYear;
		int timeStep;
		int numInitYears;
};

#endif
