#include "additionalOutputInputFile.h"
#include "common.h"

// AY, include ------->
#include <netcdf.h>
#include <stdio.h>
#include <cstring>
#define ERRCODE 2
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(ERRCODE);}
// <-------- AY

AdditionalOutputInputFile::AdditionalOutputInputFile()
{
    _additionalOutputInput.resize(ng); //ng=66896 => all cells

    // CH_laf:
    for(int i=0; i<ng;i++)
    {
        // [i][0]:G_days_since_start
        // [i][1]:G_GrowingStatus
        // [i][2]:G_PrecSum
        // [i][3]:G_totalUnsatisfiedUse
        // [i][4]:G_UnsatisfiedUsePrevYear
        // [i][5]:K_release
        // [i][6]:G_landAreaFrac
        // [i][7]: G_landAreaFracPrevTimestep
        // [i][8]:G_locWetlAreaReductionFactor          !!! obsolete
        // [i][9]:G_gloLakeEvapoReductionFactor         !!! obsolete
        // [i][10]:G_groundwaterStorage                 !!! obsolete
        // [i][11]:G_locLakeAreaReductionFactor         !!! obsolete
        // [i][12]:G_gloWetlAreaReductionFactor         !!! obsolete
        // [i][13]:G_gloResEvapoReductionFactor         !!! obsolete
        // [i][14]:G_fwsbInit[n]
        // [i][15]:G_gloWetlStorage[n]                  !!! obsolete
        // [i][16]:G_fswbLandAreaFrac                   !!! obsolete
        // [i][17]:G_locWetlStorage                     !!! obsolete
        // [i][18]:G_locLakeStorage                     !!! obsolete
        // [i][19]:G_riverStorage                       !!! obsolete
        // [i][20]:G_soilwatercontent[n]                !!! obsolete
        // [i][21]:G_gloLakeStorage[n]                  !!! obsolete
        // [i][22]: G_canopywatercontent[n]             !!! obsolete
        // [i][23] : G_gloResStorage[n]                 !!! obsolete
        //[i][24]: G_snow[n]                            !!! obsolete
        //[i][25]: G_withdrawalIrrigFromSwb[n]
        //[i][26]: G_consumptiveUseIrrigFromSwb[n]
        //[i](27]: G_unsatallocuse
        //[i][28]: G_AllocatedUse
        //[i][29]: G_Secondcell[n]
        //[i][30]; G_daily_UnsatAllocUseNextDay[n]
        //[i][31]: G_daily_allocatedUseNextDay[n]
        //[i][32]: G_AllocUseToNeigborCell[n]
        //[i][33]: G_fswbLandAreaFracNextTimestep[n]
        //[i][34]: G_PrevUnsatAllocUse[n]
        //[i][35]: G_dailyRemainigUse[n]
        //[i][36]: G_dailyAllocatedUse[n]
        //[i][37]: G_PrevTotalUnstatisfiedUse[n]
        //[i][38]: G_reducedReturnFlow[n]
        //[i][39]: G_unsatisfiedNAsFromIrrig[n]
        //[i][40]: G_dailySatisAllocatedUseInSecondCell[n]
        //[i][41]: G_dailyAllocatedUse[n]
        //[i][42]: G_unsatUseRiparian[n]
        //[i][43]: G_AllocatedUse[secondCell]
        //[i][44]: G_glores_prevyear[n]
        //[i][45]: G_fLocLake[n]
        //[i][46]: G_fGloWet[n]
        //[i][47]: G_fLocWet[n]
        //[i][48]: G_unsatUseRiparian[i]                ???
        //[i][49]: G_unsatisfiedNAsFromOtherSectors[n]
        //[i][50]: G_reducedReturnFlowPrevYear [n]
        //[i][51] G_unsatisfiedNAsFromIrrigPrevYear[n]
        //[i][52] G_unsatisfiedNAsFromOtherSectorsPrevYear[n]

        _additionalOutputInput.at(i).resize(53);
        for(int j=0; j<5;j++)
        {
            _additionalOutputInput[i][j] = 0.;
        }
        _additionalOutputInput[i][5] = 0.1;
        _additionalOutputInput[i][6] = 0.;
        _additionalOutputInput[i][7] = 0.;
        _additionalOutputInput[i][8] = 1.;
        _additionalOutputInput[i][9] = 1.;
        _additionalOutputInput[i][10] = 0.;
        _additionalOutputInput[i][11] = 1.;
        _additionalOutputInput[i][12] = 1.;
        _additionalOutputInput[i][13] = 1.;
        _additionalOutputInput[i][14] = 0.;
        _additionalOutputInput[i][15] = 0.;
        _additionalOutputInput[i][16] = 0.;
        _additionalOutputInput[i][17] = 0.;
        _additionalOutputInput[i][18] = 0.;
        _additionalOutputInput[i][19] = 0.;
        _additionalOutputInput[i][20] = 0.;
        _additionalOutputInput[i][21] = 0.;
        _additionalOutputInput[i][22] = 0.;
        _additionalOutputInput[i][23] = 0.;
        _additionalOutputInput[i][24] = 0.;
        _additionalOutputInput[i][25] = 0.;
        _additionalOutputInput[i][26] = 0.;
        _additionalOutputInput[i][27] = 0.;
        _additionalOutputInput[i][28] = 0.;
        _additionalOutputInput[i][29] = 0.;
        _additionalOutputInput[i][30] = 0.;
        _additionalOutputInput[i][31] = 0.;
        _additionalOutputInput[i][32] = 0.;
        _additionalOutputInput[i][33] = 0.;
        _additionalOutputInput[i][34] = 0.;
        _additionalOutputInput[i][35] = 0.;
        _additionalOutputInput[i][36] = 0.;
        _additionalOutputInput[i][37] = 0.;
        _additionalOutputInput[i][38] = 0.;
        _additionalOutputInput[i][39] = 0.;
        _additionalOutputInput[i][40] = 0.;
        _additionalOutputInput[i][41] = 0.;
        _additionalOutputInput[i][42] = 0.;
        _additionalOutputInput[i][43] = 0.;
        _additionalOutputInput[i][44] = 0.;
        _additionalOutputInput[i][45] = 0.;
        _additionalOutputInput[i][46] = 0.;
        _additionalOutputInput[i][47] = 0.;
        _additionalOutputInput[i][48] = 0.;
        _additionalOutputInput[i][49] = 0.;
        _additionalOutputInput[i][50] = 0.;
        _additionalOutputInput[i][51] = 0.;
        _additionalOutputInput[i][52] = 0.;

    }
}


void AdditionalOutputInputFile::load(std::string filename)
{
    // AY, add load from netcdf file -------->
    
    // determine extension of filename-string
    if (filename.substr(filename.find_last_of(".") + 1) == "nc")
    {
        int ncid, varid;
        int retval;
        int dimid[2];
        size_t dimensions[2];
        int linesIN, colsIN;
        
        // open file
        if ((retval = nc_open(filename.c_str(), NC_NOWRITE, &ncid)))
            ERR(retval);
        
        // get id of the variable
        if ((retval = nc_inq_varid(ncid, "additionalOutputInput", &varid)))
            ERR(retval);
        
        // get size of the variable
        if ((retval = nc_inq_vardimid(ncid, varid, dimid)))
            ERR(retval);
        
        if ((retval = nc_inq_dimlen(ncid, dimid[0], &dimensions[0])))
            ERR(retval);
        
        if ((retval = nc_inq_dimlen(ncid, dimid[1], &dimensions[1])))
            ERR(retval);
        
        linesIN = dimensions[0];
        colsIN = dimensions[1];
        
        // read data ELEMENT WISE - 103s
        for (unsigned int j=0;j<colsIN;j++)
        {
            for (unsigned int i=0;i<linesIN;i++)
            {
                const size_t add_idx[] = {i,j};
                if ((retval = nc_get_var1_double(ncid, varid, add_idx, &additionalOutputInput(i,j))))
                    ERR(retval);
            }
        }
        
        // close file
        if ((retval = nc_close(ncid)))
            ERR(retval);
    }
    else if (filename.substr(filename.find_last_of(".") + 1) == "txt")
    {
        // <---------- AY
    try{
        std::ifstream stream;
        stream.open(filename.c_str(), std::ios::binary);
        if(!stream.good())
            throw(Exception("error by opening file"));
        stream.exceptions(std::ios::badbit|std::ios::failbit);

        int l=0;
        while(!stream.eof())
        {
            std::string line;

            try
            {
                getline(stream, line);
            }
            catch(std::exception &e)
            {
                if(stream.eof())
                    break;
                throw(Exception(e.what()));
            }

            if(line.empty())
                break;

            std::stringstream ss(line);
            char c;
            ss>>c;
            ss.putback(c);

            if(!isalpha(c))
            {
                int id;
                ss>>id;

                for(int k=0;k<53;k++)
                {
                    ss>>additionalOutputInput(l,k);
                }
                l++;
            }

        }

    }
    catch(std::exception &e)
    {
        throw(Exception(std::string("In AdditionalOutputInputFile::load():\n")+e.what()));
    }
        } // AY, don't forget to close the if-loop
}


void AdditionalOutputInputFile::save(std::string filename)
{
    // AY, add load from netcdf file -------->
    
    // determine extension of filename-string
    if (filename.substr(filename.find_last_of(".") + 1) == "nc")
    {
        try
        {
            // define additional variables for netcdf files
            int NDIMS=2;
            int linesIN, colsIN;
            int ncid, col_dimid, row_dimid, varid_lon, varid_lat, varid_ewh;
            int dimids[NDIMS];
            size_t chunks[NDIMS];
            int shuffle, deflate, deflate_level;
            int x, y, retval;
            int varid[1];
            
            // set additional variables
            linesIN=ng;
            colsIN =53;
            shuffle = NC_SHUFFLE;
            deflate = 1;
            deflate_level = 1;
            
            // write into a netcdf file
            
            // initialize the nc file
            if ((retval = nc_create(filename.c_str(), NC_NETCDF4, &ncid)))
                ERR(retval);
            
            // define dimensions
            if ((retval = nc_def_dim(ncid, "rows", linesIN, &row_dimid)))
                ERR(retval);
            if ((retval = nc_def_dim(ncid, "columns", colsIN, &col_dimid)))
                ERR(retval);
            
            dimids[0] = row_dimid;
            dimids[1] = col_dimid;
            chunks[0] = 180;
            chunks[1] = 1;
            
            // define variables
            if ((retval = nc_def_var(ncid, "additionalOutputInput", NC_DOUBLE, NDIMS,
                                     dimids, &varid[0])))
                ERR(retval);
            
            // define compression
            if ((retval = nc_def_var_chunking(ncid, varid[0], 0, &chunks[0])))
                ERR(retval);
            
            if ((retval = nc_def_var_deflate(ncid, varid[0], shuffle, deflate, deflate_level)))
                ERR(retval);
            
            // assign values to the variable ELEMENT WISE
            for (unsigned int i=0;i<linesIN;i++)
            {
                for (unsigned int j=0;j<colsIN;j++)
                {
                    const size_t add_idx[] = {i,j};
                    if ((retval = nc_put_var1_double(ncid, varid[0], add_idx, &additionalOutputInput(i,j))))
                        ERR(retval);
                }
            }
            
            // close the nc file
            if ((retval = nc_close(ncid)))
                ERR(retval);
            
        }
        catch(std::exception &e)
        {
            throw(Exception(std::string("In AdditionalOutputInputFile::save():\n")+e.what()));
        }
    }
    else if (filename.substr(filename.find_last_of(".") + 1) == "txt")
    {
        // <---------- AY
    try{
        std::ofstream stream;
        stream.open(filename.c_str(), std::ios::binary);
        if(!stream.good())
            throw(Exception("error by opening file"));
        stream.exceptions(std::ios::badbit|std::ios::failbit);

        int width = 33;

        // CH_laf:
        //Header schreiben:
        stream << "additional out- and input" <<std::endl;
        stream << std::setw(6)  << std::setfill(' ') << "ID";
        stream << std::setw(width) << std::setfill(' ') << "G_days_since_start";
        stream << std::setw(width) << std::setfill(' ') << "G_GrowingStatus";
        stream << std::setw(width) << std::setfill(' ') << "G_PrecSum";
        stream << std::setw(width) << std::setfill(' ') << "G_totalUnsatisfiedUse";
        stream << std::setw(width) << std::setfill(' ') << "G_UnsatisfiedUsePrevYear";
        stream << std::setw(width) << std::setfill(' ') << "K_release";
        stream << std::setw(width) << std::setfill(' ') << "G_landAreaFrac";
        stream << std::setw(width) << std::setfill(' ') << "G_landAreaFracPrevTimestep";
        stream << std::setw(width) << std::setfill(' ') << "G_locWetlAreaReductionFactor";
        stream << std::setw(width) << std::setfill(' ') << "G_gloLakeEvapoReductionFactor";
        stream << std::setw(width) << std::setfill(' ') << "G_groundwaterStorage";
        stream << std::setw(width) << std::setfill(' ') << "G_locLakeAreaReductionFactor";
        stream << std::setw(width) << std::setfill(' ') << "G_gloWetlAreaReductionFactor";
        stream << std::setw(width) << std::setfill(' ') << "G_gloResEvapoReductionFactor";
        stream << std::setw(width) << std::setfill(' ') << "G_fswbInit[n]";
        stream << std::setw(width) << std::setfill(' ') << "G_gloWetlStorage";
        stream << std::setw(width) << std::setfill(' ') << "G_fswbLandAreaFrac";
        stream << std::setw(width) << std::setfill(' ') << "G_locWetlStorage[n]";
        stream << std::setw(width) << std::setfill(' ') << "G_locLakeStorage[n]";
        stream << std::setw(width) << std::setfill(' ') << "G_riverStorage[n]";
        stream << std::setw(width) << std::setfill(' ') << "G_soilwatercontent[n]";
        stream << std::setw(width) << std::setfill(' ') << "G_gloLakeStorage[n]";
        stream << std::setw(width) << std::setfill(' ') << "G_canopywatercontent[n]";
        stream << std::setw(width) << std::setfill(' ') << "G_gloResStorage[n]";
        stream << std::setw(width) << std::setfill(' ') << "G_snow[n]";
        stream << std::setw(width) << std::setfill(' ') << "G_withdrawalIrrigFromSwb[n]";
        stream << std::setw(width) << std::setfill(' ') << "G_consumptiveUseIrrigFromSwb[n]";
        stream << std::setw(width) << std::setfill(' ') << "G_unsatAllocUSE";
        stream << std::setw(width) << std::setfill(' ') << "G_AllocUSE";
        stream << std::setw(width) << std::setfill(' ') << "G_secondcell";
        stream << std::setw(width) << std::setfill(' ') << "G_daily_UnsatAllocUseNextDay[n]";
        stream << std::setw(width) << std::setfill(' ') << " G_daily_allocatedUseNextDay[n]";
        stream << std::setw(width) << std::setfill(' ') << " G_AllocUSETOneigborcell[n]";
        stream << std::setw(width) << std::setfill(' ') << " G_fswbLandAreaFracNextTimestep[n]";
        stream << std::setw(width) << std::setfill(' ') << " G_PrevUnsatAllocUse[n]";
        stream << std::setw(width) << std::setfill(' ') << " G_dailyRemainingUse";
        stream << std::setw(width) << std::setfill(' ') << " G_dailyAllocatedUse[n]";
        stream << std::setw(width) << std::setfill(' ') << " G_PrevTotalUnstatisfieduse[n]";
        stream << std::setw(width) << std::setfill(' ') << " G_reducedReturnFlow[n]";
        stream << std::setw(width) << std::setfill(' ') << " G_unsatisfiedNAsFromIrrig[n]";
        stream << std::setw(width) << std::setfill(' ') << " G_dailySatisAllocatedUseInSecondCell";
        stream << std::setw(width) << std::setfill(' ') << " G_dailyAllocatedUse";
        stream << std::setw(width) << std::setfill(' ') << " G_unsatUseRiparian[n]";
        stream << std::setw(width) << std::setfill(' ') << " G_ActualUse[secondCell]";
        stream << std::setw(width) << std::setfill(' ') << " G_glores_prevyear[n]";
        stream << std::setw(width) << std::setfill(' ') << " G_fLocLake[n]";
        stream << std::setw(width) << std::setfill(' ') << " G_fGloWet[n]";
        stream << std::setw(width) << std::setfill(' ') << " G_fLocWet[n]";
        stream << std::setw(width) << std::setfill(' ') << " G_unsatUseRiparian[i]";
        stream << std::setw(width) << std::setfill(' ') << " G_unsatisfiedNAsFromOtherSectors[n]";
        stream << std::setw(width) << std::setfill(' ') << " G_reducedReturnFlowPrevYear[n]";
        stream << std::setw(width) << std::setfill(' ') << " G_unsatisfiedNAsFromIrrigPrevYear[n]";
        stream << std::setw(width) << std::setfill(' ') << " G_unsatisfiedNAsFromOtherSectorsPrevYear[n]";
        stream << std::endl;
        stream.precision(16);

        // CH_laf:
        for(int i=0;i<ng;i++)
        {
            stream << std::setw(6)     << std::setfill(' ') << i+1;
            for(int j=0;j<53;j++) {
                stream << std::setw(width) << std::setfill(' ') << std::scientific << additionalOutputInput(i, j); //std::scientific
            }

            if (i<ng-1){
                stream << std::endl;
            }
        }


        stream.close();
    }
    catch(std::exception &e)
    {
        throw(Exception(std::string("In AdditionalOutputInputFile::save():\n")+e.what()));
    }
        } // AY, don't forget to close the if-loop
}

