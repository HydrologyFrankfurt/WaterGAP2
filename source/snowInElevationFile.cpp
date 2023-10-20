#include "snowInElevationFile.h"
#include "common.h"
  
// AY, include ------->
#include <netcdf.h>
#include <stdio.h>
#include <cstring>
#define ERRCODE 2
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(ERRCODE);}
// <-------- AY

SnowInElevationFile::SnowInElevationFile()
{
   _numLandCells = ng;
   _numLevels = 101;  // first column contains mean elevation (?), 100 elevation levels
   _snowInElevation.resize(ng);

   for(int i=0; i<ng;i++)
   {
     _snowInElevation.at(i).resize(_numLevels);
     for(int j=0; j<101;j++)
     {
       _snowInElevation[i][j] = 0.;
     }
   }
}
  
  
void SnowInElevationFile::load(std::string filename)
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
        if ((retval = nc_inq_varid(ncid, "snowInElevation", &varid)))
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
                const size_t snow_idx[] = {i,j};
                if ((retval = nc_get_var1_double(ncid, varid, snow_idx, &snowInElevation(i,j))))
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
      for(int k=0;k<101;k++)
      { 
        ss>>snowInElevation(l,k);
      }
      l++;
    }

  }
  
}
catch(std::exception &e)
{
  throw(Exception(std::string("In SnowInElevationFile::load():\n")+e.what()));
}
        } // AY, don't forget to close the if-loop
}


void SnowInElevationFile::save(std::string filename)
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
            colsIN =101;
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
            if ((retval = nc_def_var(ncid, "snowInElevation", NC_DOUBLE, NDIMS,
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
                    const size_t snow_idx[] = {i,j};
                    if ((retval = nc_put_var1_double(ncid, varid[0], snow_idx, &snowInElevation(i,j))))
                        ERR(retval);
                }
            }
            
            // close the nc file
            if ((retval = nc_close(ncid)))
                ERR(retval);
            
        }
        catch(std::exception &e)
        {
            throw(Exception(std::string("In SnowInElevationFile::save():\n")+e.what()));
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

  //Header schreiben:
  stream << "Snow in elevation for land cells" <<std::endl;
  stream.precision(16);
  
  //ZEILENWEISE IN DATEI SCHREIBEN MIT KOPIE!!!!
  
  for(int i=0;i<ng;i++)
  { 
    stream << i+1 <<"\t" ; // save ID
    for(int j=0;j<100;j++)
    {  
      stream << std::scientific << snowInElevation(i,j) <<"\t" ;
    }
    if(i<ng-1)
    {
      stream << std::scientific << snowInElevation(i,100) << std::endl;
    }
    else
      stream << std::scientific << snowInElevation(i,100);
  }
  
  stream.close();
}
catch(std::exception &e)
{
  throw(Exception(std::string("In SnowInElevationFile::save():\n")+e.what()));
}
        } // AY, don't forget to close the if-loop
}

