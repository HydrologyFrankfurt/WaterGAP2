#if !defined (_calib_param_h_)
#define _calib_param_h_

//############################################################################################################
// Tool:	"calib_param.cpp/h"
// Purpose:	Treat WaterGAP (external) Calibration Parameters
//
// Author:		Felix Portmann
// Affiliation:	Institute of Physical Geography, Goethe-University Frankfurt am Main, Germany
//
// Status and details:
// see calib_param.cpp
//############################################################################################################

//
// Includes for json11 library and namespace
// https://github.com/dropbox/json11.git 
//

#include <iostream>     // e.g. std::cout, std::endl
#include "json11.hpp"
// includes in json11.hpp
// #include <string>	// std::string, std::to_string
// #include <vector>
// #include <map>
// #include <memory>
// #include <initializer_list>

// Potential further includes for testing of json11 library 
//#include <cstdio>
//#include <cstring>
//#include <sstream>
//#include <cassert>
//#include <list>
//#include <set>
//#include <unordered_map>


//
// Includes for WaterGAP
// Felix Portmann, Institute of Physical Geography, Geothe-University Frankfurt am Main
//	using previous WaterGAP code in files watergap.cpp, gridio.cpp/h and timestring.cpp/h
//	from other developers
//
//#include <string> // already in json11.hpp
//#include <vector> // already in json11.hpp
// Treating file streams
#include <fstream> // e.g. std::ifstream, std::ofstream, std::ios
using std::ifstream;
using std::ofstream;
using std::ios;
// from <iostream>
using std::cout;
using std::cerr;
using std::endl;
// from <string>
using std::to_string;
#include <iomanip>      // std::setprecision (for testing)
using std::setprecision;

// #include <cstdbool>      // std::bool (for testing)

// basic number of grid cells of WaterGAP
// ng = 66896 for SLM
// ng = 67420 for WLM
#include "def.h"

// for time
#include <sys/time.h> //for gettimeofday()
#include "timestring.h"

// PERHAPS OF WaterGAP Version AEMS (University of Bonn)
// #include <vector>
// #include <string>
// #include <iostream>
// #include <fstream> 
// #include <sstream> 
// #include <vector>
// #include <stdlib.h>

// Include simple exception message handling of Torsten Mayer-Guerr (University of Bonn, Germany)
#include "exception.h"

#include "json11.hpp"

//size_t n_descriptors = 13;
//size_t n_ordinators = 2;
//size_t n_parameters = 26;

#define n_descriptors 13
#define n_ordinators 2
#define n_parameters 26

// Set precision (number of digits) for writing parameter data with cout
// short n_digits_coutprecision = 20;
#define n_digits_coutprecision 20

// test writing
// number of elements to write out during test reading
// long int n_elements_test = 10; // number of elements to write out during test reading
#define n_elements_test 10


// Calibration parameters (for external calibration)
// Included in standard WaterGAP calibration: P_GAMRUN_C, P_CFA, P_CFS
// Enumeration of short names for use in functions/calls
// FP 2015-10
enum eCalibParam {
	P_GAMRUN_C,
	P_CFA,
	P_CFS,
	M_ROOT_D,
	M_RIVRGH_C,
	P_LAK_D,
	P_WET_D,
	P_SWOUTF_C,
	M_EVAREDEX,
	M_NETRAD,
	P_PTC_HUM,
	P_PTC_ARI,
	P_PET_MXDY,
	P_MCWH,
	M_LAI,
	P_T_SNOWFZ,
	P_T_SNOWMT,
	M_DEGDAY_F,
	P_T_GRADNT,
	M_GW_F,
	M_RG_MAX,
	P_PCRITGWA,
	P_GWOUTF_C,
	M_NETABSSW,
	M_NETABSGW,
	M_PREC,
	First = P_GAMRUN_C,
	Last = M_PREC
};
// end eCalibParam

// Operators necessary to iterate over enum eCalibParam
// http://stackoverflow.com/questions/8498300/allow-for-range-based-for-with-enum-classes
// 2015-10-21
eCalibParam operator++( eCalibParam& x );
eCalibParam operator*(eCalibParam c);
eCalibParam begin(eCalibParam r);
 // end iterator needs to return one past the end!
eCalibParam end(eCalibParam r);


// Class of Calibration parameters
// (possibly from external calibration, e.g. with multi-criteria decision Borg algorithm)
class calibParamClass
{

// Methods available at other sources, e.g. watergap.cpp
public:

	// Fill the vectors of names
	void defineVNamesCalibParamJson();

	// Read cell-specific calibration parameters from text file (structure JSON)
	// Parsing according to json11.hpp
	void readJson(std::string str_pathfilename);

	// Call of detailed checks of parsed JSON object
	void callCheckJson();

	// Get cell-specific calibration parameter value for a specific parameter
	inline double getValue(eCalibParam eCalibParam_local, int n_local) {
		// Requisite: Check after parsing, within readjson()
		// that calibParamJson.is_object() == true

	//		cout << "CHECK START calibParam.getValue()" << endl;
	//		cout << "CHECK eCalibParam_local = " << eCalibParam_local << endl;
	//		cout << "CHECK n_local = " << n_local << endl;

			double dbl_value = 0.0;

			switch(eCalibParam_local){

				case P_GAMRUN_C:
	//				cout << "P_GAMRUN_C" << endl;
					dbl_value = calibParamJson["gammaHBV_runoff_coeff"][n_local].number_value();
					break;
				case P_CFA:
	//				cout << "P_CFA" << endl;
					dbl_value = calibParamJson["CFA_cellCorrFactor"][n_local].number_value();
					break;
				case P_CFS:
	//				cout << "P_CFS" << endl;
					dbl_value = calibParamJson["CFS_statCorrFactor"][n_local].number_value();
					break;
				case M_ROOT_D:
	//				cout << "M_ROOT_D" << endl;
					dbl_value = calibParamJson["root_depth_multiplier"][n_local].number_value();
					break;
				case M_RIVRGH_C:
	//				cout << "M_RIVRGH_C" << endl;
					dbl_value = calibParamJson["river_roughness_coeff_mult"][n_local].number_value();
					break;
				case P_LAK_D:
	//				cout << "P_LAK_D" << endl;
					dbl_value = calibParamJson["lake_depth"][n_local].number_value();
					break;
				case P_WET_D:
	//				cout << "P_WET_D" << endl;
					dbl_value = calibParamJson["wetland_depth"][n_local].number_value();
					break;
				case P_SWOUTF_C:
	//				cout << "P_SWOUTF_C" << endl;
					dbl_value = calibParamJson["surfacewater_outflow_coefficient"][n_local].number_value();
					break;
				case M_EVAREDEX:
	//				cout << "M_EVAREDEX" << endl;
					dbl_value = calibParamJson["evapo_red_fact_exp_mult"][n_local].number_value();
					break;
				case M_NETRAD:
	//				cout << "M_NETRAD" << endl;
					dbl_value = calibParamJson["net_radiation_mult"][n_local].number_value();
					break;
				case P_PTC_HUM:
	//				cout << "P_PTC_HUM" << endl;
					dbl_value = calibParamJson["PT_coeff_humid"][n_local].number_value();
					break;
				case P_PTC_ARI:
	//				cout << "P_PTC_ARI" << endl;
					dbl_value = calibParamJson["PT_coeff_arid"][n_local].number_value();
					break;
				case P_PET_MXDY:
	//				cout << "P_PET_MXDY" << endl;
					dbl_value = calibParamJson["max_daily_PET"][n_local].number_value();
					break;
				case P_MCWH:
	//				cout << "P_MCWH" << endl;
					dbl_value = calibParamJson["mcwh"][n_local].number_value();
					break;
				case M_LAI:
	//				cout << "M_LAI" << endl;
					dbl_value = calibParamJson["LAI_mult"][n_local].number_value();
					break;
				case P_T_SNOWFZ:
	//				cout << "P_T_SNOWFZ" << endl;
					dbl_value = calibParamJson["snow_freeze_temp"][n_local].number_value();
					break;
				case P_T_SNOWMT:
	//				cout << "P_T_SNOWMT" << endl;
					dbl_value = calibParamJson["snow_melt_temp"][n_local].number_value();
					break;
				case M_DEGDAY_F:
	//				cout << "M_DEGDAY_F" << endl;
					dbl_value = calibParamJson["degree_day_factor_mult"][n_local].number_value();
					break;
				case P_T_GRADNT:
	//				cout << "P_T_GRADNT" << endl;
					dbl_value = calibParamJson["temperature_gradient"][n_local].number_value();
					break;
				case M_GW_F:
	//				cout << "M_GW_F" << endl;
					dbl_value = calibParamJson["gw_factor_mult"][n_local].number_value();
					break;
				case M_RG_MAX:
	//				cout << "M_RG_MAX" << endl;
					dbl_value = calibParamJson["rg_max_mult"][n_local].number_value();
					break;
				case P_PCRITGWA:
	//				cout << "P_PCRITGWA" << endl;
					dbl_value = calibParamJson["pcrit_aridgw"][n_local].number_value();
					break;
				case P_GWOUTF_C:
	//				cout << "P_GWOUTF_C" << endl;
					dbl_value = calibParamJson["groundwater_outflow_coeff"][n_local].number_value();
					break;
				case M_NETABSSW:
	//				cout << "M_NETABSSW" << endl;
					dbl_value = calibParamJson["net_abstraction_surfacewater_mult"][n_local].number_value();
					break;
				case M_NETABSGW:
	//				cout << "M_NETABSGW" << endl;
					dbl_value = calibParamJson["net_abstraction_groundwater_mult"][n_local].number_value();
					break;
				case M_PREC:
	//				cout << "M_PREC" << endl;
					dbl_value = calibParamJson["precip_mult"][n_local].number_value();
					break;
				default:
					cerr << "calibParamClass::getValue() - unknown parameter name/number '" << eCalibParam_local << "' set value to zero" << endl;
					// dbl_value = 0.0; // as initialized
					break;
			}
			// end switch

	//		cout << "CHECK dbl_value = " << dbl_value << endl;
			return dbl_value;

	}
	// end calibParamClass::getValue()


	// Print cell-specific calibration parameter value for a specific parameter
	// With precision as defined in header
	void printValue(eCalibParam eCalibParam_local, int n_local);

	// Print cell-specific calibration parameter value of all parameters
	void printValueAllParam(int n_local);

	// Check the speed of iteration over all values of all calibration parameters
	void checkJsonSpeed();

// Standard: use methods only within class and subclasses = "protected"
protected:

	// Detailed check of parsed JSON object
	void checkJson(json11::Json calibParamJson);

	// Check ascending GCRC sequence in parsed JSON object
	void checkJsonGCRCseq(json11::Json calibParamJson);

	// JSON object for storing the cell-specific calibration parameters
	// with other auxiliary variables (descriptors & ordinators)
	json11::Json calibParamJson;

	// Iterate over all cell-specific calibration parameter value of all parameters
	// For checking the performance
	void iterateGetAllValueAllParam();

	// Conversion of directory as string to directory as string
//	void calibParamClass::convert_to_str(std::string str_string, char char250);


//  // CalibrationParameterFile();


//  void read(std::string);
  
//  inline double rootDepthMult(){return _rootDepthMult;}
//  inline float runoffCoefficient(){return _runoffCoefficient;}
//  inline double riverRoughnessCoeffMult(){return _riverRoughnessCoeffMult;}
//  inline double lakeDepth(){return _lakeDepth;}
//  inline double wetlandDepth(){return _wetlandDepth;}
//  inline double swOutflowCoeff(){return _swOutflowCoeff;}
//  inline double netRadiationMult(){return _netRadiationMult;}
//  inline double PTcoeffHumid(){return _PTcoeffHumid;}
//  inline double PTcoeffArid(){return _PTcoeffArid;}
//  inline double maxDailyPET(){return _maxDailyPET;}
//  inline double mcwh(){return _mcwh;}
//  inline float laiMult(){return _laiMult;}
//  inline double snowFreezeTemp(){return _snowFreezeTemp;}
//  inline double snowMeltTemp(){return _snowMeltTemp;}
//  inline double degreeDayFactorMult(){return _degreeDayFactorMult;}
//  inline double temperatureGradient(){return _temperatureGradient;}
//  inline float gwFactorMult(){return _gwFactorMult;}
//  inline float rgMaxMult(){return _rgMaxMult;}
//  inline double pCritAridGW(){return _pCritAridGW;}
//  inline double gwOutflowCoeff(){return _gwOutflowCoeff;}
//  inline float netAbstractionSWMult(){return _netAbstractionSWMult;}
//  inline float netAbstractionGWMult(){return _netAbstractionGWMult;}
//  inline double precipMult(){return _precipMult;}

	std::string str_string;
	char char250[250];

//double _rootDepthMult;
//float _runoffCoefficient;
//double _riverRoughnessCoeffMult;
//double _lakeDepth;
//double _wetlandDepth;
//double _swOutflowCoeff;
//double _netRadiationMult;
//double _PTcoeffHumid;
//double _PTcoeffArid;
//double _maxDailyPET;
//double _mcwh;
//float _laiMult;
//double _snowFreezeTemp;
//double _snowMeltTemp;
//double _degreeDayFactorMult;
//double _temperatureGradient;
//float _gwFactorMult;
//float _rgMaxMult;
//double _pCritAridGW;
//double _gwOutflowCoeff;
//float _netAbstractionSWMult;
//float _netAbstractionGWMult;
//double _precipMult;
};
// end class calibParam

#endif
