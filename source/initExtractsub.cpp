// OE 01-2018 (based among others on  MS 06/2014 (enKF2wghmState_dailyDA.cpp))
//#include <iostream>
//#include <iomanip>
//#include <fstream>
//#include <vector>
#include<iterator>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include "def.h"


#include "common.h"
#include "globals.h"
#include "initExtractsub.h"
#include "parameterJsonFile.h"

void init_extractsub_ (const char* string_ids,const char* string_config, WghmStateFile *& wghmState, double *& output, const char* str_calPar_file,long *oy,calibParamClass *& calParam, long *total_nr_calPar,long *calpar_size,int *calPar_index,const char* calPar_filename,int *groupmatrixindex,const char* pathTowghmMean,WghmStateFile *& wghmMean)

{
    try{
        /*
         ********************************************************************************
         read cell information of specified basin
         ********************************************************************************
         */
        std::string fileName(string_ids);
        std::ifstream stream;
        //create an input file stream
        stream.open(fileName.c_str(), std::ios::binary);
        if(!stream.good())  // check if opened correctly
        throw(Exception("error by opening file"));;
        std::vector<uint> ids;
        std::string line;  // skip header
        getline(stream, line);
        // read linewise data: ID, longitude and latitude
        while(!stream.eof())  // read data until end of file
        {
            int ID;
            double lambda;
            double phi;
            getline(stream,line); // get current line
            if(!line.empty())  // if line is not empty
            {
                std::stringstream  ss(line);
                ss>> ID >> lambda >> phi; // read 1 double, 1 double and another double
                ids.push_back( ID );
                //cerr << ID << ", " << lambda << ", " << phi << endl ;
            }
        }
        stream.close();
        std::string fileNameConfig(string_config);
        wghmState=new WghmStateFile(ng,1);
        wghmState->load(fileNameConfig);
        
        std::string fileNameWghmMean(pathTowghmMean);
        wghmMean=new WghmStateFile(ng,1);
        wghmMean->load(fileNameWghmMean);
        
        calParam=new calibParamClass;
        calParam->defineVNamesCalibParamJson();
        std::string fileNameCalParFile(str_calPar_file);
        calParam->readJson(fileNameCalParFile);
        /*
         ********************************************************************************
         write states and parameter values to output array
         ********************************************************************************
         */
        int nr_par;
        nr_par=*calpar_size;
        output=new double [ids.size()*10+nr_par];
        std::cout <<"size is "<< ids.size()*10+nr_par << std::endl;
        int j = 0;
        for(int i=0;i<ids.size();i++)
        {
            output[j++] = wghmState->cell(ids.at(i)-1).mean().canopy(0)- wghmMean->cell(ids.at(i)-1).canopy(0);
            output[j++] = wghmState->cell(ids.at(i)-1).mean().snow(0)- wghmMean->cell(ids.at(i)-1).snow(0);
            output[j++] = wghmState->cell(ids.at(i)-1).mean().soil(0) - wghmMean->cell(ids.at(i)-1).soil(0);
            output[j++] = wghmState->cell(ids.at(i)-1).mean().locallake(0)- wghmMean->cell(ids.at(i)-1).locallake(0);
            output[j++] = wghmState->cell(ids.at(i)-1).mean().localwetland(0)- wghmMean->cell(ids.at(i)-1).localwetland(0);
            output[j++] = wghmState->cell(ids.at(i)-1).mean().globallake(0)- wghmMean->cell(ids.at(i)-1).globallake(0);
            output[j++] = wghmState->cell(ids.at(i)-1).mean().globalwetland(0)- wghmMean->cell(ids.at(i)-1).globalwetland(0);
            output[j++] = wghmState->cell(ids.at(i)-1).mean().reservoir(0) - wghmMean->cell(ids.at(i)-1).reservoir(0);
            output[j++] = wghmState->cell(ids.at(i)-1).mean().river(0)- wghmMean->cell(ids.at(i)-1).river(0);
            output[j++] = wghmState->cell(ids.at(i)-1).mean().groundwater(0)- wghmMean->cell(ids.at(i)-1).groundwater(0);
        }
        
        if (nr_par > 0){
            {
                int total_nr_par;
                total_nr_par=*total_nr_calPar;
                parameterJsonFile paramJson;
                std::vector<string> vct_parameter(total_nr_par);
                vct_parameter = paramJson.getNamesCalPar();
            
                //            for(int i=0;i<5;i++){
                //                std::cout <<"CHECK: calPar_index in cpp = "<< calPar_index[i] << std::endl;
                //            }
            
                int nr_cda_unit;
                nr_cda_unit=*oy;
                // determine how many cells per subbasin
                int ind_size[nr_cda_unit];
                for (int i=0;i<nr_cda_unit;i++){
                    int count=0;
                    for(int j=0;j<ids.size();j++){
                    if (groupmatrixindex[j+i*ids.size()] != 0){
                        count++;
                    }
                    }
                    ind_size[i]= count;
                    cout << "ind_size in initExtractsub.cpp[" << i << "]: " << ind_size[i] << endl;
                }
                /*
                 ********************************************************************************
                 read calPar information
                 ********************************************************************************
                 */
                int count=0;
                // loop over all subbasins
                for (int i=0;i<nr_cda_unit;i++){
                    // loop over all calibration Parameters (=26)
                    for(int j=0;j<total_nr_par;j++){
                        double summe=0;
                        // use only those calPar which have the corresponding calPar_index==1
                        if (calPar_index[j+i*total_nr_par] == 1){
                            // check the calPar name to compute subbasin-mean
                            if (vct_parameter[j] == "gammaHBV_runoff_coeff") {
                                for(int n=0;n<ids.size();n++){
                                    // compute the mean only over the cells within the subbasin i
                                    if (groupmatrixindex[n+i*ids.size()] != 0){
                                        summe += calParam->getValue(P_GAMRUN_C,ids.at(n)-1);
                                    }
                                }
                            }
                            else if (vct_parameter[j] == "CFA_cellCorrFactor") {
                                for(int n=0;n<ids.size();n++){
                                    // compute the mean only over the cells within the subbasin i
                                    if (groupmatrixindex[n+i*ids.size()] != 0){
                                        summe += calParam->getValue(P_CFA,ids.at(n)-1);
                                    }
                                }
                            }
                            else if (vct_parameter[j] == "CFS_statCorrFactor") {
                                for(int n=0;n<ids.size();n++){
                                    // compute the mean only over the cells within the subbasin i
                                    if (groupmatrixindex[n+i*ids.size()] != 0){
                                        summe += calParam->getValue(P_CFS,ids.at(n)-1);
                                    }
                                }
                            }
                            else if (vct_parameter[j] == "root_depth_multiplier") {
                                for(int n=0;n<ids.size();n++){
                                    // compute the mean only over the cells within the subbasin i
                                    if (groupmatrixindex[n+i*ids.size()] != 0){
                                        summe += calParam->getValue(M_ROOT_D,ids.at(n)-1);
                                    }
                                }
                            }
                            else if (vct_parameter[j] == "river_roughness_coeff_mult") {
                                for(int n=0;n<ids.size();n++){
                                    // compute the mean only over the cells within the subbasin i
                                    if (groupmatrixindex[n+i*ids.size()] != 0){
                                        summe += calParam->getValue(M_RIVRGH_C,ids.at(n)-1);
                                    }
                                }
                            }
                            else if (vct_parameter[j] == "lake_depth") {
                                for(int n=0;n<ids.size();n++){
                                    // compute the mean only over the cells within the subbasin i
                                    if (groupmatrixindex[n+i*ids.size()] != 0){
                                        summe += calParam->getValue(P_LAK_D,ids.at(n)-1);
                                    }
                                }
                            }
                            else if (vct_parameter[j] == "wetland_depth") {
                                for(int n=0;n<ids.size();n++){
                                    // compute the mean only over the cells within the subbasin i
                                    if (groupmatrixindex[n+i*ids.size()] != 0){
                                        summe += calParam->getValue(P_WET_D,ids.at(n)-1);
                                    }
                                }
                            }
                            else if (vct_parameter[j] == "surfacewater_outflow_coefficient") {
                                for(int n=0;n<ids.size();n++){
                                    // compute the mean only over the cells within the subbasin i
                                    if (groupmatrixindex[n+i*ids.size()] != 0){
                                        summe += calParam->getValue(P_SWOUTF_C,ids.at(n)-1);
                                    }
                                }
                            }
                            else if (vct_parameter[j] == "evapo_red_fact_exp_mult") {
                                for(int n=0;n<ids.size();n++){
                                    // compute the mean only over the cells within the subbasin i
                                    if (groupmatrixindex[n+i*ids.size()] != 0){
                                        summe += calParam->getValue(M_EVAREDEX,ids.at(n)-1);
                                    }
                                }
                            }
                            else if (vct_parameter[j] == "net_radiation_mult") {
                                for(int n=0;n<ids.size();n++){
                                    // compute the mean only over the cells within the subbasin i
                                    if (groupmatrixindex[n+i*ids.size()] != 0){
                                        summe += calParam->getValue(M_NETRAD,ids.at(n)-1);
                                    }
                                }
                            }
                            else if (vct_parameter[j] == "PT_coeff_humid") {
                                for(int n=0;n<ids.size();n++){
                                    // compute the mean only over the cells within the subbasin i
                                    if (groupmatrixindex[n+i*ids.size()] != 0){
                                        summe += calParam->getValue(P_PTC_HUM,ids.at(n)-1);
                                    }
                                }
                            }
                            else if (vct_parameter[j] == "PT_coeff_arid") {
                                for(int n=0;n<ids.size();n++){
                                    // compute the mean only over the cells within the subbasin i
                                    if (groupmatrixindex[n+i*ids.size()] != 0){
                                        summe += calParam->getValue(P_PTC_ARI,ids.at(n)-1);
                                    }
                                }
                            }
                            else if (vct_parameter[j] == "max_daily_PET") {
                                for(int n=0;n<ids.size();n++){
                                    // compute the mean only over the cells within the subbasin i
                                    if (groupmatrixindex[n+i*ids.size()] != 0){
                                        summe += calParam->getValue(P_PET_MXDY,ids.at(n)-1);
                                    }
                                }
                            }
                            else if (vct_parameter[j] == "mcwh") {
                                for(int n=0;n<ids.size();n++){
                                    // compute the mean only over the cells within the subbasin i
                                    if (groupmatrixindex[n+i*ids.size()] != 0){
                                        summe += calParam->getValue(P_MCWH,ids.at(n)-1);
                                    }
                                }
                            }
                            else if (vct_parameter[j] == "LAI_mult") {
                                for(int n=0;n<ids.size();n++){
                                    // compute the mean only over the cells within the subbasin i
                                    if (groupmatrixindex[n+i*ids.size()] != 0){
                                        summe += calParam->getValue(M_LAI,ids.at(n)-1);
                                    }
                                }
                            }
                            else if (vct_parameter[j] == "snow_freeze_temp") {
                                for(int n=0;n<ids.size();n++){
                                    // compute the mean only over the cells within the subbasin i
                                    if (groupmatrixindex[n+i*ids.size()] != 0){
                                        summe += calParam->getValue(P_T_SNOWFZ,ids.at(n)-1);
                                    }
                                }
                            }
                            else if (vct_parameter[j] == "snow_melt_temp") {
                                for(int n=0;n<ids.size();n++){
                                    // compute the mean only over the cells within the subbasin i
                                    if (groupmatrixindex[n+i*ids.size()] != 0){
                                        summe += calParam->getValue(P_T_SNOWMT,ids.at(n)-1);
                                    }
                                }
                            }
                            else if (vct_parameter[j] == "degree_day_factor_mult") {
                                for(int n=0;n<ids.size();n++){
                                    // compute the mean only over the cells within the subbasin i
                                    if (groupmatrixindex[n+i*ids.size()] != 0){
                                        summe += calParam->getValue(M_DEGDAY_F,ids.at(n)-1);
                                    }
                                }
                            }
                            else if (vct_parameter[j] == "temperature_gradient") {
                                for(int n=0;n<ids.size();n++){
                                    // compute the mean only over the cells within the subbasin i
                                    if (groupmatrixindex[n+i*ids.size()] != 0){
                                        summe += calParam->getValue(P_T_GRADNT,ids.at(n)-1);
                                    }
                                }
                            }
                            else if (vct_parameter[j] == "gw_factor_mult") {
                                for(int n=0;n<ids.size();n++){
                                    // compute the mean only over the cells within the subbasin i
                                    if (groupmatrixindex[n+i*ids.size()] != 0){
                                        summe += calParam->getValue(M_GW_F,ids.at(n)-1);
                                    }
                                }
                            }
                            else if (vct_parameter[j] == "rg_max_mult") {
                                for(int n=0;n<ids.size();n++){
                                    // compute the mean only over the cells within the subbasin i
                                    if (groupmatrixindex[n+i*ids.size()] != 0){
                                        summe +=  calParam->getValue(M_RG_MAX,ids.at(n)-1);
                                    }
                                }
                            }
                            else if (vct_parameter[j] == "pcrit_aridgw") {
                                for(int n=0;n<ids.size();n++){
                                    // compute the mean only over the cells within the subbasin i
                                    if (groupmatrixindex[n+i*ids.size()] != 0){
                                        summe +=  calParam->getValue(P_PCRITGWA,ids.at(n)-1);
                                    }
                                }
                            }
                            else if (vct_parameter[j] == "groundwater_outflow_coeff") {
                                for(int n=0;n<ids.size();n++){
                                    // compute the mean only over the cells within the subbasin i
                                    if (groupmatrixindex[n+i*ids.size()] != 0){
                                    summe +=  calParam->getValue(P_GWOUTF_C,ids.at(n)-1);
                                    }
                                }
                            }
                            else if (vct_parameter[j] == "net_abstraction_surfacewater_mult") {
                                for(int n=0;n<ids.size();n++){
                                    // compute the mean only over the cells within the subbasin i
                                    if (groupmatrixindex[n+i*ids.size()] != 0){
                                        summe +=  calParam->getValue(M_NETABSSW,ids.at(n)-1);
                                    }
                                }
                            }
                            else if (vct_parameter[j] == "net_abstraction_groundwater_mult") {
                                for(int n=0;n<ids.size();n++){
                                    // compute the mean only over the cells within the subbasin i
                                    if (groupmatrixindex[n+i*ids.size()] != 0){
                                        summe +=  calParam->getValue(M_NETABSGW,ids.at(n)-1);
                                    }
                                }
                            }
                            else if (vct_parameter[j] == "precip_mult") {
                                for(int n=0;n<ids.size();n++){
                                    // compute the mean only over the cells within the subbasin i
                                    if (groupmatrixindex[n+i*ids.size()] != 0){
                                        summe +=  calParam->getValue(M_PREC,ids.at(n)-1);
                                    }
                                }
                            }
                            // copy calibration parameters to output vector:
                            std::cout <<"CHECK: vct_parameter "<< vct_parameter[j] << std::endl;
                            std::cout <<"CHECK: count " << count << std::endl;
                            output[ids.size()*10 + count] = summe/ind_size[i];
                            count++;
                        }
                    }
                }
            }
            delete wghmState;
            wghmState=NULL;
            delete wghmMean;
            wghmMean=NULL;
            delete calParam;
            calParam=NULL;
        }
    }
        catch(std::exception &e)
        {
            throw(Exception(std::string("In watergap::main():\n")+e.what()));
        }
        
        return;
    }
    
    
    


