// OE 01-2018 (based on  MS 06/2014: enKF2wghmState_dailyDA.cpp and AEMS: watergap.cpp)
// OE 08-2019: modifications related to the calPar.json
// OE: 02-2020: add wghm temporal meanfield that was subtracted from absolute storages (reservoirs, snow, soil). Comment the testing for the negative values of anomaly storages (all, except reservoirs, snow, soil)
#include<iterator>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include "def.h"

#include "common.h"
//#include "globals.h"
#include "enKF2wghmState.h"
#include "matrix.h"

// all classes contain objects before assimilation, only double *field is an output after the DA!
void enkf_wghmstate_ (const char* s,double *field,double *prediction, ConfigFile *& configFile, WghmStateFile *& wghmState, AdditionalOutputInputFile *& additionalOutIn, SnowInElevationFile *& snow_in_elevation,long *step, long *total_steps,long *year, long *month,long *ny,double *& output,WghmStateFile *& wghmStateMean, const char* s2,calibParamClass *& calParam, long *calpar_size,const char* calparsample, const char* path_calpar_IDs,const char* calPar_filename, WghmStateFile *& wghmMean,double *calpar_range,long *oy, long *total_nr_calPar,int *calPar_index,int *groupmatrixindex)
{
    try{
        wghmStateMean=new WghmStateFile(ng,1);
//        std::cout << std::fixed; // to set the setprecision for double
//        std::cout << std::setprecision(16); // to set the setprecision for double
        *wghmStateMean = *wghmState;
        wghmStateMean->mean_inplace(); //monthly global mean
        wghmState->lastDay_inplace();// last global day
        
        int yr;int mo;int stepping;int lastStep;int size_doublearray; int nr_par;
        yr=*year;mo=*month;lastStep=*total_steps;stepping=*step;size_doublearray=*ny;nr_par=*calpar_size;
        
        std::cout << " stepping " << stepping << " from " << lastStep << std::endl;
        std::stringstream ss;
        ss << "_" << yr <<"-"<<std::setw(2)<<std::setfill('0')<<mo;
        std::string date = ss.str();
        
        // AEMS: Save monthly mean before assimilation
        if (!configFile->outputmeanfile.empty())
        {
            std::string fn = configFile->outputmeanfile;
            std::string::size_type pos = fn.rfind('.');
            if(pos == std::string::npos)
                fn = fn + date;
            else
                fn.insert(pos,date);
            std::cout << "Save monthly mean before assimilation to <"<<fn<<"> ... "<< std::flush;
            wghmStateMean->saveMean(fn);//,wghmMean
            std::cout <<"done."<< std::endl << std::flush;
        }
        
        
        /*
         ********************************************************************************
         read cell information of specified basin
         ********************************************************************************
         */
        std::string fileName(s);
        std::ifstream stream;
        //create an input file stream
        stream.open(fileName.c_str(), std::ios::binary);
        if(!stream.good())  // check if opend correct
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
        int numStates = (size_doublearray-nr_par)/10;
        std::cout <<"NUMSTATES:"<< numStates << std::endl << std::flush;
        /*
         ********************************************************************************
         replace wghm states of the last day of a month with DA wghm state. Only absolute storages (reservoirs, snow, soil) are not allowed to be negative--> here set them to zero, if negative
         ********************************************************************************
         */
        // update last day (absolute values) = last day + modification of model prediction (=update vector-prediction vector)
        int count = 0;
        // HG (2020/07) start: set lower and upper limits of storages
        for(int i=0;i<numStates;i++){
            wghmState->cell(ids.at(i)-1).canopy(0) += (field[count]-prediction[count++]);
            if(wghmState->cell(ids.at(i)-1).canopy(0)<0.) // check if storage is negative
                wghmState->cell(ids.at(i)-1).canopy(0) = 0.; // set storage to zero
            wghmState->cell(ids.at(i)-1).snow(0) += (field[count]-prediction[count++]);
            if(wghmState->cell(ids.at(i)-1).snow(0)<0.)
                wghmState->cell(ids.at(i)-1).snow(0) = 0.;
            if(wghmState->cell(ids.at(i)-1).snow(0)>1000.)
                wghmState->cell(ids.at(i)-1).snow(0) = 1000.;
            wghmState->cell(ids.at(i)-1).soil(0) += (field[count]-prediction[count++]);
            if(wghmState->cell(ids.at(i)-1).soil(0)<0.)
                wghmState->cell(ids.at(i)-1).soil(0) = 0.;
            wghmState->cell(ids.at(i)-1).locallake(0) += (field[count]-prediction[count++]);
            //HG: if(wghmState->cell(ids.at(i)-1).locallake(0)<0)
                //HG: wghmState->cell(ids.at(i)-1).locallake(0) = 0.;
            wghmState->cell(ids.at(i)-1).localwetland(0) += (field[count]-prediction[count++]);
            if( wghmState->cell(ids.at(i)-1).localwetland(0)<0.)
                wghmState->cell(ids.at(i)-1).localwetland(0) = 0.;
            wghmState->cell(ids.at(i)-1).globallake(0) += (field[count]-prediction[count++]);
            //HG: if( wghmState->cell(ids.at(i)-1).globallake(0)<0)
                //HG: wghmState->cell(ids.at(i)-1).globallake(0) = 0.;
            wghmState->cell(ids.at(i)-1).globalwetland(0) += (field[count]-prediction[count++]);
            if(wghmState->cell(ids.at(i)-1).globalwetland(0)<0.)
                wghmState->cell(ids.at(i)-1).globalwetland(0) = 0.;
            wghmState->cell(ids.at(i)-1).reservoir(0) += (field[count]-prediction[count++]);
            if(wghmState->cell(ids.at(i)-1).reservoir(0)<0.)
                wghmState->cell(ids.at(i)-1).reservoir(0) = 0.;
            wghmState->cell(ids.at(i)-1).river(0) += (field[count]-prediction[count++]);
            if(wghmState->cell(ids.at(i)-1).river(0)<0.)
                wghmState->cell(ids.at(i)-1).river(0) = 0.;
            wghmState->cell(ids.at(i)-1).groundwater(0) += (field[count]-prediction[count++]);
            // no check if storage is negative because of groundwater depletion (negative storages are allowed)
        }
        // HG (2020/07) end: set lower and upper limits of storage
        if (!configFile->outputlastdayfile.empty())
        {
            std::cout <<"done with replacing the last whgmState day "<< count  << std::endl;
            std::string fn = configFile->outputlastdayfile;
            std::cout << "Save last day to <"<<fn<<"> ... "<< std::flush;
            wghmState->saveDay(fn,0);//,wghmMean
            std::cout <<"done."<< std::endl << std::flush;
        }
        
        /*
         ********************************************************************************
         calculate update of calibration parameters
         ********************************************************************************
         */
        if(nr_par > 0) // define parameter values new only if calibration approach is chosen
        {
            int total_nr_par=*total_nr_calPar;
            std::vector<std::string> vct_parameter(total_nr_par);
            parameterJsonFile paramJson(configFile->calibrationfile);
            //cout << "CHECK: paramJson.get_gammaHBV.at(0) = " << paramJson.get_gammaHBV().at(0)<< endl;
            vct_parameter = paramJson.getNamesCalPar();
            int nr_cda_unit = *oy;
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
                cout << "ind_size[" << i << "]: " << ind_size[i] << endl;
            }
            
            Matrix range(total_nr_par,2);
            for (int i=0;i<2;i++){
                for(int j=0;j<total_nr_par;j++){
                    range(j,i)=calpar_range[j+i*total_nr_par];
                }
            }
            //             ********************************************************************************
            //             set calPar based on calibrated values and init-values otherwise
            //             ********************************************************************************
            Matrix mat(total_nr_par,nr_cda_unit);
            // over subbasins
            int index = 0;
            for (int i=0;i<nr_cda_unit;i++){
                // over 26 calPar
                for(int j=0;j<total_nr_par;j++){
                    double summe=0;
                    // use only those calPar which have the corresponding calPar_index==1
                    if (calPar_index[j+i*total_nr_par] == 1){
                        mat(j,i)= field[size_doublearray-nr_par+index];
                        cout << "field[" << size_doublearray-nr_par+index << "]: " << field[size_doublearray-nr_par+index] << endl;
                        //is calibration parameter out of range?
                        if(mat(j,i)<range(j,0))
                            mat(j,i) = range(j,0);
                        if(mat(j,i)>range(j,1))
                            mat(j,i) = range(j,1);
                        index++;
                    }
                    //compute mean over subbasin from json-File (original input) if this parameter is not calibrated within CDA
                    else if (calPar_index[j+i*total_nr_par] == 0){
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
                        mat(j,i) = summe/ind_size[i];
                    }
                    
                }
            }
//            for (int i=0;i<nr_cda_unit;i++){
//                // over 26 calPar
//                for(int j=0;j<total_nr_par;j++){
//                    cout << "CHECK: mat[" << j <<  i << "]: " << mat(j,i) << endl;
//                }
//            }
            
            /*
             ********************************************************************************
             set calPar based on calibrated values and init-values otherwise
             ********************************************************************************
             */

            
            //     save calibration parameters
            //    ---------------------------------
            std::string outNamecalibParametersTimeEvolution;
            std::string outNamecalibParametersTimeEvolution1 = calparsample; outNamecalibParametersTimeEvolution.append(outNamecalibParametersTimeEvolution1);outNamecalibParametersTimeEvolution.append(date);outNamecalibParametersTimeEvolution.append(".txt");
            std::cout << "save time evolution of calibration parameters to file <" << outNamecalibParametersTimeEvolution << ">" << std::flush;
            paramJson.save_cda_txt(outNamecalibParametersTimeEvolution,nr_cda_unit,mat);
            std::cout <<"done."<< std::endl << std::flush;
            
            //std::cout <<"path_calpar_IDs: "<< path_calpar_IDs  << std::endl;
            paramJson.parameterJsonFile_cda(mat,groupmatrixindex,nr_cda_unit,ids.size());
            std::string fn = configFile->outputparameter;
            std::cout << "save calibration parameters to file <"<<fn<<"> ... "<< std::flush;
            std::string filenameInputArcID = path_calpar_IDs;
            paramJson.save(fn,filenameInputArcID);
            std::cout <<"done."<< std::endl << std::flush;
        }
        
        /*
         ********************************************************************************
         calculate update of snow in elevation + add wghm temporal field
         ********************************************************************************
         */
        // MS 06/2014 (enKF2wghmState_dailyDA.cpp) program applied to this context
        // sort updated value into row ids.at(i)-1 from snow_in_elevation.snowInElevation !!!
        // update of snow in elevation
        output=new double [numStates];
        count = 1;
        for(int i=0;i<numStates;i++)
        {
            // total snow storage before assimilation equal to zero?
            if(wghmStateMean->cell(ids.at(i)-1).snow(0) == 0)
            {
                for(int j=1;j<101;j++)
                {
                    // new snow in elevaltion value        =  1/100 of assimilated total snow value
                    snow_in_elevation->snowInElevation(ids.at(i)-1,j) =  (field[count] + wghmMean->cell(ids.at(i)-1).snow(0))/100;
                    if(snow_in_elevation->snowInElevation(ids.at(i)-1,j)<0) // check if storage is negative
                        snow_in_elevation->snowInElevation(ids.at(i)-1,j) = 0.; // set storage to zero
                    if(snow_in_elevation->snowInElevation(ids.at(i)-1,j)>1000.) // // HG (2020/07) check if storage is >1000
                        snow_in_elevation->snowInElevation(ids.at(i)-1,j) = 1000.; // // HG (2020/07) set storage to 1000
                }
            }
            else
            {
                // factor for snow in elevation: absolute snow after assimilation devided by absolute snow before assimilation
                output[i] = (field[count] + wghmMean->cell(ids.at(i)-1).snow(0)) / wghmStateMean->cell(ids.at(i)-1).snow(0);
                for(int j=1;j<101;j++)
                {
                    snow_in_elevation->snowInElevation(ids.at(i)-1,j) *=  output[i];
                    if(snow_in_elevation->snowInElevation(ids.at(i)-1,j)<0.) // check if storage is negative
                        snow_in_elevation->snowInElevation(ids.at(i)-1,j) = 0.; // set storage to zero
                    if(snow_in_elevation->snowInElevation(ids.at(i)-1,j)>1000.) // // HG (2020/07) check if storage is >1000
                        snow_in_elevation->snowInElevation(ids.at(i)-1,j) = 1000.; // // HG (2020/07) set storage to 1000
                }
            }
            count += 10;
        }
        
        if (!configFile->outputsnowlastdayfile.empty()) // AEMS
        {
            std::string fn = configFile->outputsnowlastdayfile;
            std::cout << "Save snow in elevation of last day to <"<<fn<<"> ... "<< std::flush;
            snow_in_elevation->save(fn);
            std::cout <<"done."<< std::endl << std::flush;
        }
        /*
         ********************************************************************************
         replace the monthly wghm states with DA wghm state + add back the temporal wghm meanfield to reservoirs, snow, soil, which is subtracted in extract_sub.cpp. We reduced the temporal men field from all compartmets to calibrate them to the same time span.  Only absolute storages (reservoirs, snow, soil) are not allowed to be negative--> here set them to zero, if negative
         ********************************************************************************
         */
        count = 0;
        // HG (2020/07) start: set lower and upper limits of storage
        for(int i=0;i<numStates;i++){ //size_doublearray-1
            wghmStateMean->cell(ids.at(i)-1).canopy(0) = field[count++]+ wghmMean->cell(ids.at(i)-1).canopy(0);
            if(wghmStateMean->cell(ids.at(i)-1).canopy(0)<0.) // check if storage is negative
                wghmStateMean->cell(ids.at(i)-1).canopy(0) = 0.; // set storage to zero
            
            wghmStateMean->cell(ids.at(i)-1).snow(0)= field[count++] + wghmMean->cell(ids.at(i)-1).snow(0);
            if(wghmStateMean->cell(ids.at(i)-1).snow(0)<0.)
                wghmStateMean->cell(ids.at(i)-1).snow(0) = 0.;
            if(wghmStateMean->cell(ids.at(i)-1).snow(0)>1000.)
                wghmStateMean->cell(ids.at(i)-1).snow(0) = 1000.;
            
            wghmStateMean->cell(ids.at(i)-1).soil(0)= field[count++] + wghmMean->cell(ids.at(i)-1).soil(0);
            if(wghmStateMean->cell(ids.at(i)-1).soil(0)<0.)
                wghmStateMean->cell(ids.at(i)-1).soil(0) = 0.;
            
            wghmStateMean->cell(ids.at(i)-1).locallake(0)= field[count++]+ wghmMean->cell(ids.at(i)-1).locallake(0);
            //    if(wghmStateMean->cell(ids.at(i)-1).locallake(0)<0)
            //      wghmStateMean->cell(ids.at(i)-1).locallake(0) = 0.;
            
            wghmStateMean->cell(ids.at(i)-1).localwetland(0)= field[count++]+ wghmMean->cell(ids.at(i)-1).localwetland(0);
            if( wghmStateMean->cell(ids.at(i)-1).localwetland(0)<0.)
                wghmStateMean->cell(ids.at(i)-1).localwetland(0) = 0.;
            
            wghmStateMean->cell(ids.at(i)-1).globallake(0)= field[count++]+ wghmMean->cell(ids.at(i)-1).globallake(0);
            //    if( wghmStateMean->cell(ids.at(i)-1).globallake(0)<0)
            //      wghmStateMean->cell(ids.at(i)-1).globallake(0) = 0.;
            
            wghmStateMean->cell(ids.at(i)-1).globalwetland(0)= field[count++]+ wghmMean->cell(ids.at(i)-1).globalwetland(0);
            if(wghmStateMean->cell(ids.at(i)-1).globalwetland(0)<0.)
                wghmStateMean->cell(ids.at(i)-1).globalwetland(0) = 0.;
            
            wghmStateMean->cell(ids.at(i)-1).reservoir(0)= field[count++] + wghmMean->cell(ids.at(i)-1).reservoir(0);
            if(wghmStateMean->cell(ids.at(i)-1).reservoir(0)<0.)
                wghmStateMean->cell(ids.at(i)-1).reservoir(0) = 0.;
            
            wghmStateMean->cell(ids.at(i)-1).river(0)= field[count++] + wghmMean->cell(ids.at(i)-1).river(0);

            if(wghmStateMean->cell(ids.at(i)-1).river(0)<0.)
                wghmStateMean->cell(ids.at(i)-1).river(0) = 0.;
            
            wghmStateMean->cell(ids.at(i)-1).groundwater(0)= field[count++] + wghmMean->cell(ids.at(i)-1).groundwater(0);
            // no check if storage is negative because of groundwater depletion (negative storages are allowed)
        }
        // HG (2020/07) end: set lower and upper limits of storages
        std::cout <<"done with replacing the monthly mean whgmState "<< count  << std::endl;
        
        std::string str;
        std::string str11 = s2; // path_to_mean_update
        std::string str12 = "states_mean_update_";
        std::string str2 = configFile->outputmeanfile;
        str.append(str11);str.append(str12);
        // extract the ensemble number
        if (str2.substr(str2.find_last_of(".") + 1) == "nc")
        {
            str.append(str2.end()-6,str2.end()-3);
            str.append(date);str.append(".nc");
        }
        else if (str2.substr(str2.find_last_of(".") + 1) == "txt")
        {
            str.append(str2.end()-7,str2.end()-4);
            str.append(date);str.append(".txt");
        }
        std::cout << "Save monthly mean after assimilation to <"<<str<<"> ... "<< std::flush;
        wghmStateMean->saveMean(str);//, wghmMean
        std::cout <<"done."<< std::endl << std::flush;
        
        if (!configFile->outputadditionalfile.empty()) // AEMS
        {
            std::string fn = configFile->outputadditionalfile;
            std::cout << "Save G_days_since_start, G_growingStatus and G_PrecSum of last day to <"<<fn<<"> ... "<< std::flush;
            additionalOutIn->save(fn);
            std::cout <<"done."<< std::endl << std::flush;
        }
        
        delete output;
        output=NULL;
        std::cout << "die Adresse vom freigegebenen Speicher fuer Faktor ist " << output << ">" << std::endl;
        delete wghmStateMean;
        wghmStateMean=NULL;
        std::cout << "die Adresse vom freigegebenen Speicher fuer wghmStateMean ist " << wghmStateMean << ">" << std::endl;
        if (stepping == lastStep-1) {
            delete wghmState;
            wghmState=NULL;
            std::cout << "die Adresse vom freigegebenen Speicher fuer wghmState ist " << wghmState << ">" << std::endl;
            delete calParam;
            calParam=NULL;
            std::cout << "die Adresse vom freigegebenen Speicher fuer calParam ist " << calParam << ">" << std::endl;
            delete additionalOutIn;
            additionalOutIn=NULL;
            std::cout << "die Adresse vom freigegebenen Speicher fuer additionalOutIn ist " << additionalOutIn << ">" << std::endl;
            delete snow_in_elevation;
            snow_in_elevation=NULL;
            std::cout << "die Adresse vom freigegebenen Speicher fuer snow_in_elevation ist " << snow_in_elevation << ">" << std::endl;
            delete configFile;
            configFile=NULL;
            std::cout << "die Adresse vom freigegebenen Speicher fuer configFile ist " << configFile << ">" << std::endl;
        }
        
    }
    
    catch(std::exception &e)
    {
        throw(Exception(std::string("In watergap::main():\n")+e.what()));
    }
    
    return;
}



// 2 remove
void enkf_wghmstate(const char* s,double *field,double *prediction, ConfigFile *& configFile, WghmStateFile *& wghmState, AdditionalOutputInputFile *& additionalOutIn, SnowInElevationFile *& snow_in_elevation,long *step, long *total_steps,long *year, long *month,long *ny,double *& output,WghmStateFile *& wghmStateMean, const char* s2,calibParamClass *& calParam, long *calpar_size,const char* calparsample, const char* path_calpar_IDs,const char* calPar_filename, WghmStateFile *& wghmMean,double *calpar_range,long *oy, long *total_nr_calPar,int *calPar_index,int *groupmatrixindex)
{
    enkf_wghmstate_(s,field,prediction,configFile,wghmState, additionalOutIn, snow_in_elevation, step,total_steps,year,month,ny,output,wghmStateMean,s2,calParam,calpar_size,calparsample,path_calpar_IDs,calPar_filename,wghmMean,calpar_range,oy,total_nr_calPar,calPar_index,groupmatrixindex); //keine funkition wird hier definiert, sondern nur aufgerufen
}


