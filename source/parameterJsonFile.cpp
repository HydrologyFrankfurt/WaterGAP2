#include "parameterJsonFile.h"
// based on Kerstin Schulze (02/2019) - schulze@geod.uni-bonn.de
// Olga Engels,Kerstin Schulze 12.03.2020: allow ALL calPar being spatially variable

parameterJsonFile::parameterJsonFile()
{ 
    // Standard Values of GlobalCDA Project (02/2019)
    std::vector<double> x(ng,1);
    _gammaHBV=x;
    _CFA=x;
    _CFS=x;
    _rootDepthMult=x;
    _riverRoughnessCoeffMult=x;
    _lakeDepth=x;
    _wetlandDepth=x;
    _swOutflowCoeff=x;
    _evapoRedFactExpMult=x;
    _netRadiationMult=x;
    _PTcoeffHumid=x;
    _PTcoeffArid=x;
    _maxDailyPET=x;
    _mcwh=x;
    _laiMult=x;
    _snowFreezeTemp=x;
    _snowMeltTemp=x;
    _degreeDayFactorMult=x;
    _temperatureGradient=x;
    _gwFactorMult=x;
    _rgMaxMult=x;
    _pCritAridGW=x;
    _gwOutflowCoeff=x;
    _netAbstractionSWMult=x;
    _netAbstractionGWMult=x;
    _precipMult=x;
    for(int i=0; i<ng; i++){
        _gammaHBV.at(i)=1;
        _CFA.at(i)=1.0;
        _CFS.at(i)=1.0;
        _rootDepthMult.at(i)=1;
        _riverRoughnessCoeffMult.at(i)=3;
        _lakeDepth.at(i)=5;
        _wetlandDepth.at(i)=2;
        _swOutflowCoeff.at(i)=0.01;
        _evapoRedFactExpMult.at(i)=1;
        _netRadiationMult.at(i)=1;
        _PTcoeffHumid.at(i)=1.26;
        _PTcoeffArid.at(i)=1.74;
        _maxDailyPET.at(i)=15;
        _mcwh.at(i)=0.3;
        _laiMult.at(i)=1;
        _snowFreezeTemp.at(i)=2;
        _snowMeltTemp.at(i)=0;
        _degreeDayFactorMult.at(i)=1;
        _temperatureGradient.at(i)=0.006;
        _gwFactorMult.at(i)=1;
        _rgMaxMult.at(i)=1;
        _pCritAridGW.at(i)=12.5;
        _gwOutflowCoeff.at(i)=0.01;
        _netAbstractionSWMult.at(i)=1;
        _netAbstractionGWMult.at(i)=1;
        _precipMult.at(i)=1;
    }
}
parameterJsonFile::parameterJsonFile(std::string filenameJson)
{
    calibParamClass calParam;
    calParam.defineVNamesCalibParamJson();
    calParam.readJson(filenameJson);
    
    std::vector<double> x(ng,1);
    _gammaHBV=x;
    _CFA=x;
    _CFS=x;
    _rootDepthMult=x;
    _riverRoughnessCoeffMult=x;
    _lakeDepth=x;
    _wetlandDepth=x;
    _swOutflowCoeff=x;
    _evapoRedFactExpMult=x;
    _netRadiationMult=x;
    _PTcoeffHumid=x;
    _PTcoeffArid=x;
    _maxDailyPET=x;
    _mcwh=x;
    _laiMult=x;
    _snowFreezeTemp=x;
    _snowMeltTemp=x;
    _degreeDayFactorMult=x;
    _temperatureGradient=x;
    _gwFactorMult=x;
    _rgMaxMult=x;
    _pCritAridGW=x;
    _gwOutflowCoeff=x;
    _netAbstractionSWMult=x;
    _netAbstractionGWMult=x;
    _precipMult=x;
    for(int n=0; n<ng; n++){
        _gammaHBV.at(n)=calParam.getValue(P_GAMRUN_C,n);
        _CFA.at(n)=calParam.getValue(P_CFA,n);
        _CFS.at(n)=calParam.getValue(P_CFS,n);
        _rootDepthMult.at(n)            =calParam.getValue(M_ROOT_D,n);
        _riverRoughnessCoeffMult.at(n)=calParam.getValue(M_RIVRGH_C,n);
        _lakeDepth.at(n)                =calParam.getValue(P_LAK_D,n);
        _wetlandDepth.at(n)            =calParam.getValue(P_WET_D,n);
        _swOutflowCoeff.at(n)            =calParam.getValue(P_SWOUTF_C,n);
        _evapoRedFactExpMult.at(n)    =calParam.getValue(M_EVAREDEX,n);
        _netRadiationMult.at(n)        =calParam.getValue(M_NETRAD,n);
        _PTcoeffHumid.at(n)            =calParam.getValue(P_PTC_HUM,n);
        _PTcoeffArid.at(n)            =calParam.getValue(P_PTC_ARI,n);
        _maxDailyPET.at(n)            =calParam.getValue(P_PET_MXDY,n);
        _mcwh.at(n)                    =calParam.getValue(P_MCWH,n);
        _laiMult.at(n)                =calParam.getValue(M_LAI,n);
        _snowFreezeTemp.at(n)            =calParam.getValue(P_T_SNOWFZ,n);
        _snowMeltTemp.at(n)            =calParam.getValue(P_T_SNOWMT,n);
        _degreeDayFactorMult.at(n)    =calParam.getValue(M_DEGDAY_F,n);
        _temperatureGradient.at(n)    =calParam.getValue(P_T_GRADNT,n);
        _gwFactorMult.at(n)            =calParam.getValue(M_GW_F,n);
        _rgMaxMult.at(n)                =calParam.getValue(M_RG_MAX,n);
        _pCritAridGW.at(n)            =calParam.getValue(P_PCRITGWA,n);
        _gwOutflowCoeff.at(n)            =calParam.getValue(P_GWOUTF_C,n);
        _netAbstractionSWMult.at(n)    =calParam.getValue(M_NETABSSW,n);
        _netAbstractionGWMult.at(n)    =calParam.getValue(M_NETABSGW,n);
        _precipMult.at(n)                =calParam.getValue(M_PREC,n);
    }
}

void parameterJsonFile::parameterJsonFile_cda(Matrix &mat,int *groupmatrixindex,int nr_cda_unit, int ids)
{
    //paramJson.get_precipMult.at(0)
    //in den vorhandenen json-file updated calPar per basin reinschreiben
    for (int i=0;i<nr_cda_unit;i++){
        for(int n=0;n<ids;n++){
                if (groupmatrixindex[n+i*ids] > 0){
                    _gammaHBV.at(groupmatrixindex[n+i*ids]-1)=mat(0,i);
                    _CFA.at(groupmatrixindex[n+i*ids]-1)=mat(1,i);
                    _CFS.at(groupmatrixindex[n+i*ids]-1)=mat(2,i);
                    _rootDepthMult.at(groupmatrixindex[n+i*ids]-1)=mat(3,i);
                    _riverRoughnessCoeffMult.at(groupmatrixindex[n+i*ids]-1)=mat(4,i);
                    _lakeDepth.at(groupmatrixindex[n+i*ids]-1)=mat(5,i);
                    _wetlandDepth.at(groupmatrixindex[n+i*ids]-1)=mat(6,i);
                    _swOutflowCoeff.at(groupmatrixindex[n+i*ids]-1)=mat(7,i);
                    _evapoRedFactExpMult.at(groupmatrixindex[n+i*ids]-1)=mat(8,i);
                    _netRadiationMult.at(groupmatrixindex[n+i*ids]-1)=mat(9,i);
                    _PTcoeffHumid.at(groupmatrixindex[n+i*ids]-1)=mat(10,i);
                    _PTcoeffArid.at(groupmatrixindex[n+i*ids]-1)=mat(11,i);
                    _maxDailyPET.at(groupmatrixindex[n+i*ids]-1)=mat(12,i);
                    _mcwh.at(groupmatrixindex[n+i*ids]-1)=mat(13,i);
                    _laiMult.at(groupmatrixindex[n+i*ids]-1)=mat(14,i);
                    _snowFreezeTemp.at(groupmatrixindex[n+i*ids]-1)=mat(15,i);
                    _snowMeltTemp.at(groupmatrixindex[n+i*ids]-1)=mat(16,i);
                    _degreeDayFactorMult.at(groupmatrixindex[n+i*ids]-1)=mat(17,i);
                    _temperatureGradient.at(groupmatrixindex[n+i*ids]-1)=mat(18,i);
                    _gwFactorMult.at(groupmatrixindex[n+i*ids]-1)=mat(19,i);
                    _rgMaxMult.at(groupmatrixindex[n+i*ids]-1)=mat(20,i);
                    _pCritAridGW.at(groupmatrixindex[n+i*ids]-1)=mat(21,i);
                    _gwOutflowCoeff.at(groupmatrixindex[n+i*ids]-1)=mat(22,i);
                    _netAbstractionSWMult.at(groupmatrixindex[n+i*ids]-1)=mat(23,i);
                    _netAbstractionGWMult.at(groupmatrixindex[n+i*ids]-1)=mat(24,i);
                    _precipMult.at(groupmatrixindex[n+i*ids]-1)=mat(25,i);
                }
            }
    }
}

void parameterJsonFile::setCFAandCFSfromJson(std::string filenameJson)
{ 
    calibParamClass calParam;
    calParam.defineVNamesCalibParamJson();
    calParam.readJson(filenameJson);
    
    std::vector<double> x(ng,1);
    _CFA=x;
    _CFS=x;
    for(int n=0; n<ng; n++){
        _CFA.at(n)=calParam.getValue(P_CFA,n);
        _CFS.at(n)=calParam.getValue(P_CFS,n);
    }
}
void parameterJsonFile::setCFAandCFStoConstantValue(double value)
{ 
    calibParamClass calParam;
    calParam.defineVNamesCalibParamJson();
    
    std::vector<double> x(ng,1);
    _CFA=x;
    _CFS=x;
    for(int n=0; n<ng; n++){
        _CFA.at(n)=value;
        _CFS.at(n)=value;
    }
}

void parameterJsonFile::setGammaHBVfromJson(std::string filenameJson)
{ 
    calibParamClass calParam;
    calParam.defineVNamesCalibParamJson();
    calParam.readJson(filenameJson);
    
    std::vector<double> x(ng,1);
    _gammaHBV=x;
    for(int n=0; n<ng; n++){
        _gammaHBV.at(n)=calParam.getValue(P_GAMRUN_C,n);
    }
}

// new parameter version used within the GlobalCDA project (for WaterGAP 2.2d)
void parameterJsonFile::save(std::string filenameOutput, std::string filenameInputArcID)
{
    //gridioClass gridIO; //OE uncommented
    // OE: extract ng from def.h
    long int n_elements = ng; //ng_standardMask;
    
    // DESCRIPTOR
    char chrarr_time[250]; // for getISOdateTime (timestring.h/cpp)
    getISOdateTime(chrarr_time, sizeof(chrarr_time));
    // define and fill vector of descriptor names
    std::vector<string> vct_descriptor_name(n_descriptors);
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
        // cout << "vct_descriptor_name[" << i << "]: " << vct_descriptor_name[i] << endl;
    }
    // define and fill vector of descriptor values
    // ATTENTION: For numbers (unless they are written as string in quotes):
    // (1) Do not write leading zeroes (e.g. 0123)! (correct: 123)
    // (2) Do not write singular decimal points (e.g. .1, 1.) (correct: 0.1, 1.0)
    // These lead to an error when parsing (no Json object is created)
    std::vector<string> vct_descriptor_value(n_descriptors);
    for (int i=0; i < n_descriptors; i++) {
        string tmp_string;
        switch (i) {
            case 0: tmp_string=string("UTF-8 without BOM"); break;
            case 1: tmp_string=to_string(n_descriptors); break;
            case 2: tmp_string=to_string(n_ordinators); break;
            case 3: tmp_string=to_string(n_parameters); break;
            case 4: tmp_string=to_string(ng); break; // number
            case 5: tmp_string=string("WLM"); break;
            case 6: tmp_string=string("WaterGAP2.2b"); break;
            case 7: tmp_string=string("1901"); break; // number
            case 8: tmp_string=string("1"); break; // number, no leading zeroes, nor singular decimal points!
            case 9: tmp_string=string(chrarr_time); break; // creation time from system time
            case 10: tmp_string=string("BonnUniversity-IGG"); break;
            case 11: tmp_string=string("APMG"); break;
            case 12: tmp_string=string("Cells filled with read-in data: (1) G_GAMMA_HBV.UNF0, G_CORR_FACTOR.UNF0, G_STAT_CORR.UNF0 (2) global parameters"); break;
        }
        vct_descriptor_value[i]=tmp_string;
        // check
        // cout << "vct_descriptor_value[" << i << "]: " << vct_descriptor_value[i] << endl;
    }
    // define and fill vector of descriptor types
    std::vector<string> vct_descriptor_type(n_descriptors);
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
        // cout << "vct_descriptor_type[" << i << "]: " << vct_descriptor_type[i] << endl;
    }
    
    // ORDINATOR (gcrc_cellnumber, arcID)
    // Define and fill vector of ordinator names
    // cout << "Define vector of ordinator names (n_ordinators: " << n_ordinators << ")" << endl;
    std::vector<string> vct_ordinator(n_ordinators);
    for (int i=0; i < n_ordinators; i++) {
        string tmp_string;
        switch (i) {
            case 0: tmp_string=string("gcrc_cellnumber"); break; // gcrc = n + 1
            case 1: tmp_string=string("arc_id"); break; // ArcID from file
        }
        vct_ordinator[i]=tmp_string;
        // cout << "vct_ordinator[" << i << "]: " << vct_ordinator[i] << endl;
    }
    
    // PARAMETERS
    // Define and fill vector of parameter names
    // cout << "Define vector of parameter names (n_parameters: " << n_parameters << ")" << endl;
    std::vector<string> vct_parameter(n_parameters);
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
        // cout << "vct_parameter[" << i << "]: " << vct_parameter[i] << endl;
    }
    
    std::string str_pathfilename = filenameInputArcID;
    // std::cout << "filename: " << str_pathfilename << endl;
    std::ifstream infile_arcid_gcrc_text(str_pathfilename, ios::in|ios::binary);
    string str_newline;
    long int arcid_in_tmp, gcrc_in_tmp;
    long int G_arcid_in[ng], G_gcrc_in[ng];
    if (infile_arcid_gcrc_text.is_open()){
        // cout << "File opened: " << str_pathfilename << endl;
        long int n_lines=0;
        // read 1st comment line
        if (!infile_arcid_gcrc_text.eof()) {
            getline(infile_arcid_gcrc_text, str_newline, '\n');
            // cout << "Header line: " << str_newline << endl;
        }
        // read all data from streamfile object to string and transfer to arrays
        while ( !infile_arcid_gcrc_text.eof() && (infile_arcid_gcrc_text >> arcid_in_tmp >> gcrc_in_tmp) ) {
            // Restore
            G_arcid_in[n_lines] = arcid_in_tmp;
            G_gcrc_in[n_lines] = gcrc_in_tmp;
            // next line number / index
            n_lines +=1;
        }
        // end whileloop over new lines
        //cout << "n_lines: " << n_lines << endl;
        // Close file
        infile_arcid_gcrc_text.close();
    } else std::cout << "Unable to open file " << str_pathfilename << endl;
    
    // Writing parameter textfile
    // Define separators for JSON object text file
    string str_begObj_Json = string("{");
    string str_endObj_Json = string("}");
    string str_begendKeyVal_Json = string("\"");
    string str_sepKeyVal_Json = string(": ");
    string str_sepData_Json = string(",");
    string str_begArray_Json = string("[");
    string str_endArray_Json = string("]");
    // create outfile stream pointer
    // Use strings
    str_pathfilename = filenameOutput;
    // cout << "      Write to file '" << str_pathfilename << "'" << endl;
    std::ofstream outfile_parameters_Json_text(str_pathfilename, ios::out|ios::binary);
    // fill outfile stream pointer
    if (outfile_parameters_Json_text.is_open()){
        // cout << "      Opened output file '" << str_pathfilename << "'" << endl;
        // Start object
        outfile_parameters_Json_text << str_begObj_Json << endl;
        // Write descriptors (one per line)
        // ATTENTION: Do not write leading zeroes with numbers (unless they are written as string in quotes)!
        // These lead to an error when parsing (no Json object is created)
        for (int i=0; i < n_descriptors; i++) {
            // Start each descriptor with the key = name of the descriptor
            outfile_parameters_Json_text << str_begendKeyVal_Json << vct_descriptor_name[i] << str_begendKeyVal_Json;
            // Separator key from value
            outfile_parameters_Json_text << str_sepKeyVal_Json;
            // Write descriptor value
            // strings with separators
            if (vct_descriptor_type[i] == string("string")) {
                outfile_parameters_Json_text << str_begendKeyVal_Json << vct_descriptor_value[i] << str_begendKeyVal_Json;
            }
            // numbers without separators
            else if (vct_descriptor_type[i] == string("number")) {
                outfile_parameters_Json_text << vct_descriptor_value[i];
            }
            // Separating descriptors as necessary
            if (i < n_descriptors-1) {
                // with separator and endofline/LF
                outfile_parameters_Json_text << str_sepData_Json << endl;
            }
        }// end for loop over descriptors
        // End list of descriptors possibly without separator (if written alone)
        // With separator when ordinators or parameters follow:
        if ((n_ordinators > 0) || (n_parameters > 0)) {
            outfile_parameters_Json_text << str_sepData_Json << endl;
        }
        
        //
        // Write ordinators (one per line)
        //
        for (int i=0; i < n_ordinators; i++) {
            // Start each ordinator array with the key = name of the array
            outfile_parameters_Json_text << str_begendKeyVal_Json << vct_ordinator[i] << str_begendKeyVal_Json;
            // Separator key from value
            outfile_parameters_Json_text << str_sepKeyVal_Json;
            // Writing array - individually
            // (1.1) start array
            outfile_parameters_Json_text << str_begArray_Json;
            // (1.2) write array
            for (int n=0; n < n_elements; n++) {
                // grcrc number of grid cell is: cell index n + 1
                long int gcrc_out = n+1;
                // if read-in values of gcrc correspond to thos automatically generated,
                // then write the gcrc numbers and the corresponding ArcID
                if (gcrc_out == G_gcrc_in[n]) {
                    switch (i) {
                        case 0: outfile_parameters_Json_text << gcrc_out; break; // gcrc = n + 1
                        case 1: outfile_parameters_Json_text << G_arcid_in[n]; break; // arcid
                    }
                } else {
                    outfile_parameters_Json_text << "inconsistent_gcrc_in_at_n" << n;
                } // end ifcorrect gcrc
                // Separator (except after last item)
                if (n < n_elements-1) {
                    outfile_parameters_Json_text << str_sepData_Json;
                }
            } // end for loop
            // (1.3) end array
            outfile_parameters_Json_text << str_endArray_Json;
            // (1.4) Separating arrays as necessary
            if (i < n_ordinators-1) {
                // with separator and endofline/LF
                outfile_parameters_Json_text << str_sepData_Json << endl;
            }
        } // end for loop over ordinators
        // With separator when parameters follow:
        if (n_parameters > 0) {
            outfile_parameters_Json_text << str_sepData_Json << endl;
        }
        //
        // Write calibration parameters (one per line)
        //
        for (int i=0; i < n_parameters; i++) {
            // Start each parameter array with the key = name of the array
            outfile_parameters_Json_text << str_begendKeyVal_Json << vct_parameter[i] << str_begendKeyVal_Json;
            // Separator key from value
            outfile_parameters_Json_text << str_sepKeyVal_Json;
            // (1.1) start array
            outfile_parameters_Json_text << str_begArray_Json;
            // Write depending on type of data to write
            // Matrix parameterMatrix = _calParamJson.getParameters();
            double value;
            for (int n=0; n < n_elements; n++) {
                switch (i) {
                    case 0:  value=_gammaHBV.at(n); break;
                    case 1:  value=_CFA.at(n); break;
                    case 2:  value=_CFS.at(n); break;
                    case 3:  value=_rootDepthMult.at(n); break;
                    case 4:  value=_riverRoughnessCoeffMult.at(n); break;
                    case 5:  value=_lakeDepth.at(n); break;
                    case 6:  value=_wetlandDepth.at(n); break;
                    case 7:  value=_swOutflowCoeff.at(n); break;
                    case 8:  value=_evapoRedFactExpMult.at(n); break;
                    case 9:  value=_netRadiationMult.at(n); break;
                    case 10: value=_PTcoeffHumid.at(n); break;
                    case 11: value=_PTcoeffArid.at(n); break;
                    case 12: value=_maxDailyPET.at(n); break;
                    case 13: value=_mcwh.at(n); break;
                    case 14: value=_laiMult.at(n); break;
                    case 15: value=_snowFreezeTemp.at(n); break;
                    case 16: value=_snowMeltTemp.at(n); break;
                    case 17: value=_degreeDayFactorMult.at(n); break;
                    case 18: value=_temperatureGradient.at(n); break;
                    case 19: value=_gwFactorMult.at(n); break;
                    case 20: value=_rgMaxMult.at(n); break;
                    case 21: value=_pCritAridGW.at(n); break;
                    case 22: value=_gwOutflowCoeff.at(n); break;
                    case 23: value=_netAbstractionSWMult.at(n); break;
                    case 24: value=_netAbstractionGWMult.at(n); break;
                    case 25: value=_precipMult.at(n); break;
                }
                outfile_parameters_Json_text << value;
                if (n < n_elements-1) {
                    outfile_parameters_Json_text << str_sepData_Json;
                }
            }
            // (1.3) end array
            outfile_parameters_Json_text << str_endArray_Json;
            // (1.4) Separating arrays as necessary
            if (i < n_parameters-1) {
                // with separator and endofline/LF
                outfile_parameters_Json_text << str_sepData_Json << endl;
            }
        } // end for loop over parameters
        // write endofline/LF after last line
        outfile_parameters_Json_text << endl;
        // End object
        outfile_parameters_Json_text << str_endObj_Json;
        // Close file
        // cout << "      Closing cell-specific parameter textfile '" << str_pathfilename << "'" << endl;
        outfile_parameters_Json_text.close();
    } else std::cout << "Unable to open file '" << str_pathfilename << "'" << endl;
    // End writing text file of JSON object with descriptors and parameters
} // end save function

std::vector<string> parameterJsonFile::getNamesCalPar()
{
    try
    {
        std::vector<std::string> vct_parameter(n_parameters);
        vct_parameter[0]=string("gammaHBV_runoff_coeff");
        vct_parameter[1]=string("CFA_cellCorrFactor");
        vct_parameter[2]=string("CFS_statCorrFactor");
        vct_parameter[3]=string("root_depth_multiplier");
        vct_parameter[4]=string("river_roughness_coeff_mult");
        vct_parameter[5]=string("lake_depth");
        vct_parameter[6]=string("wetland_depth");
        vct_parameter[7]=string("surfacewater_outflow_coefficient");
        vct_parameter[8]=string("evapo_red_fact_exp_mult");
        vct_parameter[9]=string("net_radiation_mult");
        vct_parameter[10]=string("PT_coeff_humid");
        vct_parameter[11]=string("PT_coeff_arid");
        vct_parameter[12]=string("max_daily_PET");
        vct_parameter[13]=string("mcwh");
        vct_parameter[14]=string("LAI_mult");
        vct_parameter[15]=string("snow_freeze_temp");
        vct_parameter[16]=string("snow_melt_temp");
        vct_parameter[17]=string("degree_day_factor_mult");
        vct_parameter[18]=string("temperature_gradient");
        vct_parameter[19]=string("gw_factor_mult");
        vct_parameter[20]=string("rg_max_mult");
        vct_parameter[21]=string("pcrit_aridgw");
        vct_parameter[22]=string("groundwater_outflow_coeff");
        vct_parameter[23]=string("net_abstraction_surfacewater_mult");
        vct_parameter[24]=string("net_abstraction_groundwater_mult");
        vct_parameter[25]=string("precip_mult");
        return vct_parameter;
    }
    catch(std::exception &e)
    {
        throw(Exception(std::string("In CalibrationParameterFile::getParameters():\n")+e.what()));
    }
}

// Save-function for CDA-analysis in Bonn
void parameterJsonFile::save_cda_txt(std::string filename, double nr_subbasins, Matrix &mat)
{
    
    std::vector<string> vct_parameter(n_parameters);
    parameterJsonFile paramJson;
    vct_parameter = paramJson.getNamesCalPar();
    
    std::ofstream stream;
    stream.open(filename.c_str(), std::ios::binary);
    if(!stream.good())
        throw(Exception("error by opening file"));
    stream.exceptions(std::ios::badbit|std::ios::failbit);
    
    int width = 33;
    //stream.precision(16);
      //cout << "CHECK: nr_subbasins = " << nr_subbasins << endl;
    // CH_laf:
    for(int i=0;i<n_parameters;i++)
    {
        stream << std::setw(width) << std::setfill(' ') << vct_parameter[i]; //.at(i)
        for(int j=0;j<nr_subbasins;j++) {
            //cout << "CHECK: nr_subbasins = " << j << endl;
            stream << std::setw(width) << std::setfill(' ') << std::scientific << mat(i, j); //std::scientific
        }
        stream << std::endl;
    }
    
    
    stream.close();
    
}
