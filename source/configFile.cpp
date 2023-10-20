#include "configFile.h"

//ConfigFile("","");
//ConfigFile("","",2008,7);
// ordne das folgendermaßen: zuerst variablen, dann defaults!
//ConfigFile::ConfigFile(std::string cf, std::string pn, int year=2003, int month=5)
ConfigFile::ConfigFile(std::string cf, int year, int month, std::string pn) // OE 01-2018: additional input int year and  int month,std::string)
{
try{    
  configfile = cf;
	progName = pn;
	startMonth = 0;
	startYear = 0;
	endMonth = 0;
	endYear = 0;
	timeStep = 0;
	numInitYears=0;
	
  std::cout <<"load config file <"<<configfile<<">"<<std::endl;
	// Header einlesen
	
  std::ifstream stream;
  stream.open(configfile.c_str(), std::ios::binary);
  if(!stream.good())
    throw(Exception("error by opening file"));
  stream.exceptions(std::ios::badbit|std::ios::failbit);

	std::string line;
	for(;;)
	{
		getline(stream, line);
		std::stringstream ss(line);
		std::string tag;

    char c;
    ss>>c;
    ss.putback(c);

		if(c!='#') // comment sign
		{
			ss>>tag;

			if(tag == "wghm_state")
				ss>>startvaluefile;

			if(tag == "calibration_parameters")
				ss>>calibrationfile;

            if(tag == "param_json") // HMS
                ss>>parameterfile;
			
			if(tag == "snowInElevation_startvalues") // AEMS
			       ss>>snowInElevationfile;
			
			if(tag == "additionalOutIn_startvalues") // AEMS
			       ss>>additionalfile;

			if(tag == "output_state_mean")
				ss>>outputmeanfile;

			if(tag == "output_state_lastday")
				ss>>outputlastdayfile;
			
			if(tag == "output_state_day")
				ss>>outputdailyfile;
			
			if(tag == "output_snowInElevation_lastday") // AEMS
			      ss>>outputsnowlastdayfile;
			
			if(tag == "additionalOutIn_lastday") // AEMS
			      ss>>outputadditionalfile;
            
            if(tag == "output_calibration_parameters") // OE
            ss>>outputparameter;

			if(tag == "start_month")
			{
				std::string no;
				ss>>no;
				if (progName == "OL") // OE
					startMonth = atoi(no.c_str());
				else
					startMonth = month; // OE
			}

			if(tag == "start_year")
			{
				std::string no;
				ss>>no;
				if (progName == "OL") // OE
					startYear = atoi(no.c_str());
				else
					startYear = year; // OE
			}
			
			if(tag == "end_month")
			{
				std::string no;
				ss>>no;
				if (progName == "OL") // OE
					endMonth = atoi(no.c_str());
				else
					endMonth = month;// OE
			}

			if(tag == "end_year")
			{
				std::string no;
				ss>>no;
				if (progName == "OL") // OE
					endYear = atoi(no.c_str());
				else
					endYear = year;// OE
			}

			if(tag == "time_step")
			{
				std::string no;
				ss>>no;
				timeStep = atoi(no.c_str());
			}

			if(tag == "num_init_years")
			{
				std::string no;
				ss>>no;
				numInitYears = atoi(no.c_str());
			}

			if(tag == "runtime_options")
				ss>>runtimeoptionsfile;

			if(tag == "output_options")
				ss>>outputoptionsfile;

			if(tag == "routing")
				ss>>routingoptionsfile;

			if(tag == "stations")
				ss>>stationsfile;

			if(tag == "input_dir")
				ss>>inputDir;

			if(tag == "output_dir")
				ss>>outputDir;

            if(tag == "discharge_dir")
            ss>>dischargeDir;

	        if(tag == "snowCoverFrac_dir")
	            ss>>snowCoverFracDir;        
            
            if(tag == "snowCoverFrac_dir")
            ss>>snowCoverFracDir;

if(tag == "reanPrecip_dir")
            ss>>reanPrecipDir;

            
			if(tag == "climate_dir")
				ss>>climateDir;

			if(tag == "glacier_dir")
			    ss>>glacierDir;
            
            if(tag == "climateSLRD_dir")
                ss>>climateSLRDDir;

			if(tag == "routing_dir")
				ss>>routingDir;

			if(tag == "water_use_dir")
				ss>>wateruseDir;
			
			if(tag == "land_cover_dir")
				ss>>landcoverDir;
			
			// Header fertig ?
			if(tag == "end_of_head")
				break;
		}
	}
	stream.close();
	
}
catch(std::exception &e)
{
  throw(Exception(std::string("In ConfigFile::ConfigFile():\n")+e.what()));
}
}


ConfigFile::~ConfigFile()
{
}

