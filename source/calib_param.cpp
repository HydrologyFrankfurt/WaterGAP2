//############################################################################################################
// Tool:	"calib_param.cpp/h"
// Purpose:	Treat WaterGAP (external) Calibration Parameters
//			for multicriteria-optimization of cell-specific parameter values
//			for WaterGAP2.2b
//			e.g.:
//			  Read JSON text file
//			  Parse to JSON object
//			  Check integrity and consistency of the data
//			  Get values of individual parameters
//
// Using free JSON library for C++11
// from https://github.com/dropbox/json11.git
//
// Author:		Felix Portmann
// Affiliation:	Institute of Physical Geography, Goethe-University Frankfurt am Main, Germany
//
// Status:
// 2015-10-21	Working basic methods to read, test and access JSON object
// 2015-10-14	Establishment from wg_param2jasonfile.cpp
//############################################################################################################

// Include headers, eg. jason11.hpp
#include "calib_param.h"
//#include <cstdio>


// original namespace of json11 library
using namespace json11;
using std::string;  // 

// for time accounting
extern double timediff(struct timeval tv2_local, struct timeval tv1_local);
// Difference in time
struct timeval tv1_local, tv2_local; // for gettimeofday (sys/time.h)


// to be read as parameters when calling the program
// directories
//std::string str_inputdir; // e.g. "/home/temp90/nobackup/portmann/WaterGAP2_2/JSON/JSON_WG_param2jsonfile"
//std::string str_outputdir; // e.g. "/home/temp90/nobackup/portmann/WaterGAP2_2/JSON/JSON_WG_param2jsonfile"
// // filenames
//std::string str_filename_in_paramglobal; // e.g. "calib_parameters_initial_global.json"
//std::string str_filename_in_paramcell; // e.g. "calib_parameters_initial_cell.json"

// // to be composed later
//std::string str_filename; // unspecific string for filenames (composed later)
//std::string str_pathfilename;

// String for error messages, e.g. during parsing
std::string str_err;
	

// definitions specific to WaterGAP (when possibly testing file opening procedure)
// char filename[250];
// char inputdir[250];


//
// Methods
//

// Operators necessary to iterate over enum eCalibParam
// http://stackoverflow.com/questions/8498300/allow-for-range-based-for-with-enum-classes
// 2015-10-21
eCalibParam operator++( eCalibParam& x ){ return x = (eCalibParam)(((int)(x) + 1)); }
eCalibParam operator*(eCalibParam c) {return c;}
eCalibParam begin(eCalibParam r) {return eCalibParam::First;}
// end iterator needs to return one past the end!
eCalibParam end(eCalibParam r)   {return eCalibParam(int(eCalibParam::Last) + 1);}


// Conversion of directory as character to directory as string
//void calibParamClass::convert_to_str(std::string str_string, char char250)
//{
//  str_string = string(char250);
//}
// end calibParam::convert_to_str()

// Define vectors for names
std::vector<std::string> vct_descriptor_name(n_descriptors);
std::vector<std::string> vct_descriptor_type(n_descriptors);
std::vector<std::string> vct_ordinator(n_ordinators);
std::vector<std::string> vct_parameter(n_parameters);


// Fill the vectors of names
void calibParamClass::defineVNamesCalibParamJson() {

	try {
//		cout << "START calibParamClass::defineVNamesCalibParamJson()" << endl; // CHECK

		//
		// Descriptors
		//

		// Fill vector of descriptor names

		for (int i=0; i < n_descriptors; i++) {
			string tmp_string;
			switch (i) {
				case 0: tmp_string=string("file_encoding"); break;
				case 1: tmp_string=string("n_descriptors"); break;
				case 2: tmp_string=string("n_ordinators"); break;
				case 3: tmp_string=string("n_parameters"); break;
				case 4: tmp_string=string("ng_param"); break;
				case 5: tmp_string=string("watergap_landmask"); break;
				case 6: tmp_string=string("watergap_version"); break;
				case 7: tmp_string=string("reference_year"); break;
				case 8: tmp_string=string("reference_month"); break;
				case 9: tmp_string=string("creation_datetime"); break;
				case 10: tmp_string=string("creation_institution"); break;
				case 11: tmp_string=string("creation_staff"); break;
				case 12: tmp_string=string("comments"); break;
			}
			vct_descriptor_name[i]=tmp_string;
			// check
//			cout << "vct_descriptor_name[" << i << "]: " << vct_descriptor_name[i] << endl;
		}



		// Fill vector of descriptor types

		for (int i=0; i < n_descriptors; i++) {
			string tmp_string;
			switch (i) {
				case 0: tmp_string=string("string"); break;
				case 1: tmp_string=string("number"); break;
				case 2: tmp_string=string("number"); break;
				case 3: tmp_string=string("number"); break;
				case 4: tmp_string=string("number"); break;
				case 5: tmp_string=string("string"); break;
				case 6: tmp_string=string("string"); break;
				case 7: tmp_string=string("number"); break;
				case 8: tmp_string=string("number"); break;
				case 9: tmp_string=string("string"); break;
				case 10: tmp_string=string("string"); break;
				case 11: tmp_string=string("string"); break;
				case 12: tmp_string=string("string"); break;
			}
			vct_descriptor_type[i]=tmp_string;
			// check
//			cout << "vct_descriptor_type[" << i << "]: " << vct_descriptor_type[i] << endl;
		}


		//
		// Ordinators
		//

		//
		// Fill vector of ordinator names
		//
		cout << "Define vector of ordinator names (n_ordinators: " << n_ordinators << ")" << endl;
		for (int i=0; i < n_ordinators; i++) {
			string tmp_string;
			switch (i) {
				case 0: tmp_string=string("gcrc_cellnumber"); break; // gcrc = n + 1
				case 1: tmp_string=string("arc_id"); break; // ArcID from file
				}
			vct_ordinator[i]=tmp_string;
//			cout << "vct_ordinator[" << i << "]: " << vct_ordinator[i] << endl;
		}


		//
		// Calibration parameters
		//

		//
		// Fill vector of calibration parameter names (as in file)
		//
		cout << "Define vector of calibration parameter names in file (n_parameters: " << n_parameters << ")" << endl;
		for (int i=0; i < n_parameters; i++) {
			string tmp_string;
			switch (i) {
				case 0: tmp_string=string("gammaHBV_runoff_coeff"); break; // daily.G_gammaHBV[ng]
				case 1: tmp_string=string("CFA_cellCorrFactor"); break; // daily.G_cellCorrFactor[ng]
				case 2: tmp_string=string("CFS_statCorrFactor"); break; // routing.G_statCorrFactor[ng]
				case 3: tmp_string=string("root_depth_multiplier"); break;
				case 4: tmp_string=string("river_roughness_coeff_mult"); break;
				case 5: tmp_string=string("lake_depth"); break;
				case 6: tmp_string=string("wetland_depth"); break;
				case 7: tmp_string=string("surfacewater_outflow_coefficient"); break;
				case 8: tmp_string=string("evapo_red_fact_exp_mult"); break;
				case 9: tmp_string=string("net_radiation_mult"); break;
				case 10: tmp_string=string("PT_coeff_humid"); break;
				case 11: tmp_string=string("PT_coeff_arid"); break;
				case 12: tmp_string=string("max_daily_PET"); break;
				case 13: tmp_string=string("mcwh"); break;
				case 14: tmp_string=string("LAI_mult"); break;
				case 15: tmp_string=string("snow_freeze_temp"); break;
				case 16: tmp_string=string("snow_melt_temp"); break;
				case 17: tmp_string=string("degree_day_factor_mult"); break;
				case 18: tmp_string=string("temperature_gradient"); break;
				case 19: tmp_string=string("gw_factor_mult"); break;
				case 20: tmp_string=string("rg_max_mult"); break;
				case 21: tmp_string=string("pcrit_aridgw"); break;
				case 22: tmp_string=string("groundwater_outflow_coeff"); break;
				case 23: tmp_string=string("net_abstraction_surfacewater_mult"); break;
				case 24: tmp_string=string("net_abstraction_groundwater_mult"); break;
				case 25: tmp_string=string("precip_mult"); break;
			}
			vct_parameter[i]=tmp_string;
//			cout << "vct_parameter[" << i << "]: " << vct_parameter[i] << endl;
		}

	}
	// end try

	// Catch exception if necessary
	catch(std::exception &e)
	{
		throw(Exception(std::string("Exception in calibParamClass::defineVNamesCalibParamJson() :\n")+e.what()));
	}
	// end catch

}
// end calibParam::defineVNamesCalibParamJson()


// Read cell-specific calibration parameters from text file (structure JSON)
// Parsing according to json11.cpp & json11.hpp
void calibParamClass::readJson(std::string str_pathfilename_local) {

	try {
	//	cout << "START calibParamClass::readjson(" << str_pathfilename_local << ")" << endl; // CHECK

		// Reading JSON file with parameters:
		//	(1) Reading the parameter textfile into a string
		//	(2) Parsing from string to Json object
		//	If desired:
		//	(3) Checking the Json object elements values
		//

		// string of file content for easy transfer to Json objects
		string str_filecontent;
		// string of current line read
		string str_newline;

		// to check ng in file
		int ng_fromfile=0;

		//
		// (1) Reading the calbration parameter textfile output file into a string
		//
		cout << "Reading calibration parameter textfile with JSON format" << endl;
		// Use strings
		// str_pathfilename = str_inputdir+"/"+str_filename_in_paramcell;
		cout << "Reading from file '" << str_pathfilename_local << "'" << endl;

		std::ifstream infile_calibParam_Json_text (str_pathfilename_local, ios::in|ios::binary);

		// clear string of possible previous content
		str_filecontent.clear();

		// measure time
		gettimeofday(&tv1_local, NULL);

		if (infile_calibParam_Json_text.is_open()){

	//		cout << "Opened file '" << str_pathfilename_local << "'" << endl;

			// read all data from streamfile object to character
			while ( getline(infile_calibParam_Json_text, str_newline, '\n') ) {
				// cout << "newline : " << str_newline << endl;
				str_filecontent += str_newline;
			}

			// cout << "str_filecontent : " << str_filecontent << endl;
	//		cout << "Closing read-in textfile '" << str_pathfilename_local << "'" << endl;
			infile_calibParam_Json_text.close();

		} else std::cout << "Unable to open file '" << str_pathfilename_local << "'" << endl;

		// measure time
		gettimeofday(&tv2_local, NULL);
		cout << "Time difference in Wallclock Time (seconds): = " << timediff(tv2_local, tv1_local) << endl;


		//
		// (2) Parsing from string to Json object
		//
		cout << "Parsing from string to JSON object" << endl;
		gettimeofday(&tv1_local, NULL);

		// automatically define the type of calibParamJson (e.g. ARRAY only)
		calibParamJson = json11::Json::parse(str_filecontent, str_err);
		if (str_err.length()>0){
			cout << "ERROR in parsing: '" << str_err << "'" << endl;
		}

		gettimeofday(&tv2_local, NULL);
		cout << "Time difference in Wallclock Time (seconds): = " << timediff(tv2_local, tv1_local) << endl;


		//
		// (3) Checking the Json object elements values
		//
		cout << "Overall checking of parsed JSON object" << endl;

		if (calibParamJson.is_object()) {
			cout << "OK 1 - JSON Object generated" << endl;
	//		cout << "dump calibParamJson: " << calibParamJson.dump() << "\n";
			ng_fromfile = calibParamJson["ng_param"].number_value();
			if (ng != ng_fromfile) {
				cerr << "ERROR readJson(): Predefined number of grid cells " << ng << " does not match number found in JSON object header " << ng_fromfile << endl;
				cerr << "Something went wrong with the file. Please check the file!" << endl;
				exit (-1);
			}
			else {
				cout << "OK 2 - Predefined number of grid cells also found in JSON object header: " << ng_fromfile << endl;
			}


		}
		// end if object (expected type)
		else {
			cerr << "ERROR readJson(): calibParamJson.is_object() == false" << endl;
			cerr << "Something went wrong with the file. Please check the file!" << endl;
			exit (-1);
		}

	}
	// end try

	// Catch exception if necessary
	catch(std::exception &e)
	{
		throw(Exception(std::string("Exception in calibParamClass::readJson() :\n")+e.what()));
	}
	// end catch

}
// end calibParamClass::readJson() {


// Call of detailed check of parsed JSON object
void calibParamClass::callCheckJson() {

	try {

		//	cout << "START calibParam.callCheckJson()" << endl;

		// Check format of JSON sub-objects
		checkJson(calibParamJson);

		// Check correct sequence of GCRC cell numbers
		checkJsonGCRCseq(calibParamJson);

	}
	// end try

	// Catch exception if necessary
	catch(std::exception &e)
	{
		throw(Exception(std::string("Exception in calibParamClass::callCheckJson() :\n")+e.what()));
	}
	// end catch

}
// end calibParamClass::callCheckJson();


// Detailed check of parsed JSON object
void calibParamClass::checkJson(json11::Json calibParamJson) {

	try {
	//	cout << "START calibParam.checkJson()" << endl;

		// (3.1) Loop over descriptors
		for (int i=0; i<n_descriptors; i++) {

			// cout << "i: " << i << endl;
			cout << "vct_descriptor_name[" << i << "]: " << vct_descriptor_name[i] << endl;
			// cout << "dump value: " << calibParamJson[vct_descriptor_name[i]].dump() << endl;

			if (calibParamJson[vct_descriptor_name[i]].is_string()) {
				cout << "calibParamJson[\"" << vct_descriptor_name[i] << "\"].is_string() == true" << endl;
				cout << "dump: " << calibParamJson[vct_descriptor_name[i]].dump() << endl;
			}
			// end if string (expected type)
			if (calibParamJson[vct_descriptor_name[i]].is_number()) {
				cout << "calibParamJson[\"" << vct_descriptor_name[i] << "\"].is_number() == true" << endl;
				cout << "dump: " << calibParamJson[vct_descriptor_name[i]].dump() << endl;
			}
			// end if number (expected type)

			// Erroneous array
			if (calibParamJson[vct_descriptor_name[i]].is_array()) {
				cout << "ATTENTION/ERROR: calibParamJson[\"" << vct_descriptor_name[i] << "\"].is_array() == true" << endl;
				cout << "dump: " << calibParamJson[vct_descriptor_name[i]].dump() << endl;
			}

			// Erroneous object
			if (calibParamJson[vct_descriptor_name[i]].is_object()) {
				cout << "ATTENTION/ERROR: calibParamJson[\"" << vct_descriptor_name[i] << "\"].is_object() == true" << endl;
				cout << "dump: " << calibParamJson[vct_descriptor_name[i]].dump() << endl;
			}

		}
		// end for loop over descriptors


		// (3.2) Loop over ordinators
		for (int i=0; i<n_ordinators; i++) {

			// cout << "i: " << i << endl;
			cout << "vct_ordinator[" << i << "]: " << vct_ordinator[i] << endl;
			// cout << "dump value (possibly array): " << calibParamJson[vct_ordinator[i]].dump() << endl;

			if (calibParamJson[vct_ordinator[i]].is_array()) {
				cout << "calibParamJson[\"" << vct_ordinator[i] << "\"].is_array() == true" << endl;

					// check
					for (int n=0; n < n_elements_test; n++) {
						cout << "pos: " << n  << setprecision(n_digits_coutprecision) << " number_value(): " << calibParamJson[vct_ordinator[i]][n].number_value() << endl;
					}

			}
			// end if array (expected type)

			if (calibParamJson[vct_ordinator[i]].is_object()) {
				cout << "ATTENTION/ERROR: calibParamJson[\"" << vct_ordinator[i] << "\"].is_object() == true" << endl;
				cout << "dump: " << calibParamJson[vct_ordinator[i]].dump() << endl;
			}

		}
		// end for loop over ordinators


		// (3.3) Loop over parameters
		for (int i=0; i<n_parameters; i++) {

			cout << "i: " << i << endl;
			cout << "vct_parameter[" << i << "]: " << vct_parameter[i] << endl;
			// cout << "dump value (possibly array): " << calibParamJson[vct_parameter[i]].dump() << endl;

		   if (calibParamJson[vct_parameter[i]].is_array()) {
				cout << "calibParamJson[\"" << vct_parameter[i] << "\"].is_array() == true" << endl;

					// check
					for (int n=0; n < n_elements_test; n++) {
						cout << "pos: " << n << setprecision(n_digits_coutprecision) << " number_value(): " << calibParamJson[vct_parameter[i]][n].number_value() << endl;
					}
	// full array
	//		 std::cout << "dump calibParamJson[" << vct_parameter[i] << "]: " << calibParamJson[vct_parameter[i]].dump() << "\n";
	//		for (auto &k : calibParamJson[vct_parameter[i]].array_items()) {
	//			std::cout << "    - " << k.dump() << "\n";
	//		}

			}
			// end if array (expected type)

			// Erroneous object
			if (calibParamJson[vct_parameter[i]].is_object()) {
				cout << "ATTENTION/ERROR: calibParamJson[\"" << vct_parameter[i] << "\"].is_object() == true" << endl;
				cout << "dump: " << calibParamJson[vct_parameter[i]].dump() << endl;
			}


		}
		// end for loop over parameters

	}
	// end try

	// Catch exception if necessary
	catch(std::exception &e)
	{
		throw(Exception(std::string("Exception in calibParamClass::checkJson() :\n")+e.what()));
	}
	// end catch

}
// end calibParamClass::checkJson()


// Check ascending GCRC sequence in parsed JSON object
void calibParamClass::checkJsonGCRCseq(json11::Json calibParamJson) {

	try {
	//	cout << "START calibParam.checkJsonGCRCseq()" << endl;

		double diff = 0;
		short indicator_diff_wrong = 0;
		int gcrc_firstwrong = 0;
		int gcrc_lastwrong = 0;
		int n_local = 1;

		// Use the number of grid cells in JSON object
		int ng_fromfile = 0;
		ng_fromfile = calibParamJson["ng_param"].number_value();

		// Correct limits (GCRC is n+1)
		if ( (1 == calibParamJson["gcrc_cellnumber"][0].number_value() )
			 || (ng_fromfile == calibParamJson["gcrc_cellnumber"][ng].number_value() )  )  {

			// Loop over gcrc-numbers (1 ... ng) in JSON object (from calibration file)
			for (n_local=1; n_local<ng; n_local++) {

				// Calculate difference (increment) between consecutive cells
				// current minus previous
				diff = calibParamJson["gcrc_cellnumber"][n_local].number_value()
					- calibParamJson["gcrc_cellnumber"][n_local-1].number_value();

				// correct difference is 1
				// if not, memorize the case
				if (1.0 != diff) {
					indicator_diff_wrong = 1;
					if (0 == indicator_diff_wrong) {
						gcrc_firstwrong = n_local+1; // GCRC is n+1
					}
					gcrc_lastwrong = n_local+1; // GCRC is n+1
				}

			}
			// end for loop over grid cells

			if (indicator_diff_wrong) {
				cerr << "ERROR: GCRC cell numbers increment unequal to 1, first gcrc:" << gcrc_firstwrong << " last gcrc:" << gcrc_lastwrong << endl;
				exit (-1);
			}
			else {
				cout << "GCRC cell numbers increment is to 1 and within limits (lower 1, upper " << ng_fromfile << ")" << endl;
			}

		}
		else {
			cerr << "ERROR: invalid limits of GCRC cell numbers, outside lower 1 to upper " << ng_fromfile << endl;
			exit (-1);
		}

	}
	// end try

	// Catch exception if necessary
	catch(std::exception &e)
	{
		throw(Exception(std::string("Exception in calibParamClass::checkJsonGCRCseq() :\n")+e.what()));
	}
	// end catch

}
// end calibParamClass::checkJsonGCRC()


// Currently in header as inline function withiut try & catch exception
//// Get cell-specific calibration parameter value for a specific parameter
// double calibParamClass::getValue(eCalibParam eCalibParam_local, int n_local) {
//	// Requisite: Check after parsing, within readjson()
//	// that calibParamJson.is_object() == true

//	// Possible implicit checks:
//	//	calibParamJson is object
//	//	n_local is within valid range 0 ... ng-1
//	try {

////		cout << "CHECK START calibParam.getValue()" << endl;
////		cout << "CHECK eCalibParam_local = " << eCalibParam_local << endl;
////		cout << "CHECK n_local = " << n_local << endl;

//		double dbl_value = 0.0;

//		switch(eCalibParam_local){

//			case P_GAMRUN_C:
////				cout << "P_GAMRUN_C" << endl;
//				dbl_value = calibParamJson["gammaHBV_runoff_coeff"][n_local].number_value();
//				break;
//			case P_CFA:
////				cout << "P_CFA" << endl;
//				dbl_value = calibParamJson["CFA_cellCorrFactor"][n_local].number_value();
//				break;
//			case P_CFS:
////				cout << "P_CFS" << endl;
//				dbl_value = calibParamJson["CFS_statCorrFactor"][n_local].number_value();
//				break;
//			case M_ROOT_D:
////				cout << "M_ROOT_D" << endl;
//				dbl_value = calibParamJson["root_depth_multiplier"][n_local].number_value();
//				break;
//			case M_RIVRGH_C:
////				cout << "M_RIVRGH_C" << endl;
//				dbl_value = calibParamJson["river_roughness_coeff_mult"][n_local].number_value();
//				break;
//			case P_LAK_D:
////				cout << "P_LAK_D" << endl;
//				dbl_value = calibParamJson["lake_depth"][n_local].number_value();
//				break;
//			case P_WET_D:
////				cout << "P_WET_D" << endl;
//				dbl_value = calibParamJson["wetland_depth"][n_local].number_value();
//				break;
//			case P_SWOUTF_C:
////				cout << "P_SWOUTF_C" << endl;
//				dbl_value = calibParamJson["surfacewater_outflow_coefficient"][n_local].number_value();
//				break;
//			case M_EVAREDEX:
////				cout << "M_EVAREDEX" << endl;
//				dbl_value = calibParamJson["evapo_red_fact_exp_mult"][n_local].number_value();
//				break;
//			case M_NETRAD:
////				cout << "M_NETRAD" << endl;
//				dbl_value = calibParamJson["net_radiation_mult"][n_local].number_value();
//				break;
//			case P_PTC_HUM:
////				cout << "P_PTC_HUM" << endl;
//				dbl_value = calibParamJson["PT_coeff_humid"][n_local].number_value();
//				break;
//			case P_PTC_ARI:
////				cout << "P_PTC_ARI" << endl;
//				dbl_value = calibParamJson["PT_coeff_arid"][n_local].number_value();
//				break;
//			case P_PET_MXDY:
////				cout << "P_PET_MXDY" << endl;
//				dbl_value = calibParamJson["max_daily_PET"][n_local].number_value();
//				break;
//			case P_MCWH:
////				cout << "P_MCWH" << endl;
//				dbl_value = calibParamJson["mcwh"][n_local].number_value();
//				break;
//			case M_LAI:
////				cout << "M_LAI" << endl;
//				dbl_value = calibParamJson["LAI_mult"][n_local].number_value();
//				break;
//			case P_T_SNOWFZ:
////				cout << "P_T_SNOWFZ" << endl;
//				dbl_value = calibParamJson["snow_freeze_temp"][n_local].number_value();
//				break;
//			case P_T_SNOWMT:
////				cout << "P_T_SNOWMT" << endl;
//				dbl_value = calibParamJson["snow_melt_temp"][n_local].number_value();
//				break;
//			case M_DEGDAY_F:
////				cout << "M_DEGDAY_F" << endl;
//				dbl_value = calibParamJson["degree_day_factor_mult"][n_local].number_value();
//				break;
//			case P_T_GRADNT:
////				cout << "P_T_GRADNT" << endl;
//				dbl_value = calibParamJson["temperature_gradient"][n_local].number_value();
//				break;
//			case M_GW_F:
////				cout << "M_GW_F" << endl;
//				dbl_value = calibParamJson["gw_factor_mult"][n_local].number_value();
//				break;
//			case M_RG_MAX:
////				cout << "M_RG_MAX" << endl;
//				dbl_value = calibParamJson["rg_max_mult"][n_local].number_value();
//				break;
//			case P_PCRITGWA:
////				cout << "P_PCRITGWA" << endl;
//				dbl_value = calibParamJson["pcrit_aridgw"][n_local].number_value();
//				break;
//			case P_GWOUTF_C:
////				cout << "P_GWOUTF_C" << endl;
//				dbl_value = calibParamJson["groundwater_outflow_coeff"][n_local].number_value();
//				break;
//			case M_NETABSSW:
////				cout << "M_NETABSSW" << endl;
//				dbl_value = calibParamJson["net_abstraction_surfacewater_mult"][n_local].number_value();
//				break;
//			case M_NETABSGW:
////				cout << "M_NETABSGW" << endl;
//				dbl_value = calibParamJson["net_abstraction_groundwater_mult"][n_local].number_value();
//				break;
//			case M_PREC:
////				cout << "M_PREC" << endl;
//				dbl_value = calibParamJson["precip_mult"][n_local].number_value();
//				break;
//			default:
//				cerr << "calibParamClass::getValue() - unknown parameter name/number '" << eCalibParam_local << "' set value to zero" << endl;
//				// dbl_value = 0.0; // as initialized
//				break;
//		}
//		// end switch

////		cout << "CHECK dbl_value = " << dbl_value << endl;
//		return dbl_value;

//	}
//	// end try

//	// Catch exception if necessary
//	catch(std::exception &e)
//	{
//		throw(Exception(std::string("Exception in calibParamClass::getValue() :\n")+e.what()));
//	}
//	// end catch

//}
//// end calibParamClass::getValue()


// Print cell-specific calibration parameter value for a specific parameter
// With precision as defined in header
void calibParamClass::printValue(eCalibParam eCalibParam_local, int n_local) {
	// Requisite: Check after parsing, within readjson()
	// that calibParamJson.is_object() == true

	// Possible implicit checks:
	//	calibParamJson is object
	//	n_local is within valid range 0 ... ng-1
	try {

//		cout << "CHECK calibParam.printValue(cell_n:" << n_local << ")" << endl;

		switch(eCalibParam_local){

			case P_GAMRUN_C:
				cout << "P_GAMRUN_C " << setprecision(n_digits_coutprecision) <<
						calibParamJson["gammaHBV_runoff_coeff"][n_local].number_value() << endl;
				break;
			case P_CFA:
				cout << "P_CFA " << setprecision(n_digits_coutprecision) <<
						calibParamJson["CFA_cellCorrFactor"][n_local].number_value() << endl;
				break;
			case P_CFS:
				cout << "P_CFS " << setprecision(n_digits_coutprecision) <<
						calibParamJson["CFS_statCorrFactor"][n_local].number_value() << endl;
				break;
			case M_ROOT_D:
				cout << "M_ROOT_D " << setprecision(n_digits_coutprecision) <<
						calibParamJson["root_depth_multiplier"][n_local].number_value() << endl;
				break;
			case M_RIVRGH_C:
				cout << "M_RIVRGH_C " << setprecision(n_digits_coutprecision) <<
						calibParamJson["river_roughness_coeff_mult"][n_local].number_value() << endl;
				break;
			case P_LAK_D:
				cout << "P_LAK_D " << setprecision(n_digits_coutprecision) <<
						calibParamJson["lake_depth"][n_local].number_value() << endl;
				break;
			case P_WET_D:
				cout << "P_WET_D " << setprecision(n_digits_coutprecision) <<
						calibParamJson["wetland_depth"][n_local].number_value() << endl;
				break;
			case P_SWOUTF_C:
				cout << "P_SWOUTF_C " << setprecision(n_digits_coutprecision) <<
						calibParamJson["surfacewater_outflow_coefficient"][n_local].number_value() << endl;
				break;
			case M_EVAREDEX:
				cout << "M_EVAREDEX " << setprecision(n_digits_coutprecision) <<
						calibParamJson["evapo_red_fact_exp_mult"][n_local].number_value() << endl;
				break;
			case M_NETRAD:
				cout << "M_NETRAD " << setprecision(n_digits_coutprecision) <<
						calibParamJson["net_radiation_mult"][n_local].number_value() << endl;
				break;
			case P_PTC_HUM:
				cout << "P_PTC_HUM " << setprecision(n_digits_coutprecision) <<
						calibParamJson["PT_coeff_humid"][n_local].number_value() << endl;
				break;
			case P_PTC_ARI:
				cout << "P_PTC_ARI " << setprecision(n_digits_coutprecision) <<
						calibParamJson["PT_coeff_arid"][n_local].number_value() << endl;
				break;
			case P_PET_MXDY:
				cout << "P_PET_MXDY " << setprecision(n_digits_coutprecision) <<
						calibParamJson["max_daily_PET"][n_local].number_value() << endl;
				break;
			case P_MCWH:
				cout << "P_MCWH " << setprecision(n_digits_coutprecision) <<
						calibParamJson["mcwh"][n_local].number_value() << endl;
				break;
			case M_LAI:
				cout << "M_LAI " << setprecision(n_digits_coutprecision) <<
						calibParamJson["LAI_mult"][n_local].number_value() << endl;
				break;
			case P_T_SNOWFZ:
				cout << "P_T_SNOWFZ " << setprecision(n_digits_coutprecision) <<
						calibParamJson["snow_freeze_temp"][n_local].number_value() << endl;
				break;
			case P_T_SNOWMT:
				cout << "P_T_SNOWMT " << setprecision(n_digits_coutprecision) <<
						calibParamJson["snow_melt_temp"][n_local].number_value() << endl;
				break;
			case M_DEGDAY_F:
				cout << "M_DEGDAY_F " << setprecision(n_digits_coutprecision) <<
						calibParamJson["degree_day_factor_mult"][n_local].number_value() << endl;
				break;
			case P_T_GRADNT:
				cout << "P_T_GRADNT " << setprecision(n_digits_coutprecision) <<
						calibParamJson["temperature_gradient"][n_local].number_value() << endl;
				break;
			case M_GW_F:
				cout << "M_GW_F " << setprecision(n_digits_coutprecision) <<
						calibParamJson["gw_factor_mult"][n_local].number_value() << endl;
				break;
			case M_RG_MAX:
				cout << "M_RG_MAX " << setprecision(n_digits_coutprecision) <<
						calibParamJson["rg_max_mult"][n_local].number_value() << endl;
				break;
			case P_PCRITGWA:
				cout << "P_PCRITGWA " << setprecision(n_digits_coutprecision) <<
						calibParamJson["pcrit_aridgw"][n_local].number_value() << endl;
				break;
			case P_GWOUTF_C:
				cout << "P_GWOUTF_C " << setprecision(n_digits_coutprecision) <<
						calibParamJson["groundwater_outflow_coeff"][n_local].number_value() << endl;
				break;
			case M_NETABSSW:
				cout << "M_NETABSSW " << setprecision(n_digits_coutprecision) <<
						calibParamJson["net_abstraction_surfacewater_mult"][n_local].number_value() << endl;
				break;
			case M_NETABSGW:
				cout << "M_NETABSGW " << setprecision(n_digits_coutprecision) <<
						calibParamJson["net_abstraction_groundwater_mult"][n_local].number_value() << endl;
				break;
			case M_PREC:
				cout << "M_PREC " << setprecision(n_digits_coutprecision) <<
						calibParamJson["precip_mult"][n_local].number_value() << endl;
				break;
			default:
				cerr << "calibParamClass::printValue() - unknown parameter name/number '" << eCalibParam_local << endl;
				break;
		}
		// end switch

	}
	// end try

	// Catch exception if necessary
	catch(std::exception &e)
	{
		throw(Exception(std::string("Exception in calibParamClass::printValue() :\n")+e.what()));
	}
	// end catch

}
// end calibParamClass::printValue()


// Print cell-specific calibration parameter value of all parameters
void calibParamClass::printValueAllParam(int n_local) {

	try {
		cout << "calibParam.printValueAllParam(n_cell:" << n_local << ")" << endl;
		// Loop over enum of calibration parameters
		// Syntax: http://stackoverflow.com/questions/8498300/allow-for-range-based-for-with-enum-classes
		for (const auto &eCalibParam_local : eCalibParam()) {
			printValue(eCalibParam_local, n_local);
		}
	}
	// end try

	// Catch exception if necessary
	catch(std::exception &e)
	{
		throw(Exception(std::string("Exception in calibParamClass::printValueAllParam() :\n")+e.what()));
	}
	// end catch

}
// end calibParamClass::printValueAllParam()


// Check the speed of iteration over all values of all calibration parameters
void calibParamClass::checkJsonSpeed() {

	try {

		cout << "START calibParam.checkJsonSpeed()" << endl;

		// start time
		gettimeofday(&tv1_local, NULL);

		iterateGetAllValueAllParam();

		// end time
		gettimeofday(&tv2_local, NULL);
		cout << "Time difference in Wallclock Time (seconds): = " << timediff(tv2_local, tv1_local) << endl;

	}
	// end try

	// Catch exception if necessary
	catch(std::exception &e)
	{
		throw(Exception(std::string("Exception in calibParamClass::checkJsonSpeed() :\n")+e.what()));
	}
	// end catch

}
// end calibParamClass::checkJsonSpeed()


// Iterate over all cell-specific calibration parameter value of all parameters
// For checking the performance
void calibParamClass::iterateGetAllValueAllParam() {

	try {
		cout << "START calibParam.iterateAllValueAllParam()" << endl;
		double dbl_local = 0;
		int n_local = 0;
		// Loop over enum of calibration parameters
		// Syntax: http://stackoverflow.com/questions/8498300/allow-for-range-based-for-with-enum-classes
		for (const auto &eCalibParam_local : eCalibParam()) {

			// Loop over grid cells
			for (n_local=0; n_local<calibParamJson["ng_param"].number_value(); n_local++) {
				dbl_local = getValue(eCalibParam_local, n_local);
			}
			// end loop over grid cells

		}
		// end loop over parameters
	}
	// end try

	// Catch exception if necessary
	catch(std::exception &e)
	{
		throw(Exception(std::string("Exception in calibParamClass::iterateGetAllValueAllParam() :\n")+e.what()));
	}
	// end catch

}
// end calibParamClass::iterateGetAllValueAllParam()

