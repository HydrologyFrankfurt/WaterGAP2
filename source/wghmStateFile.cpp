#include "wghmStateFile.h"
#include "common.h"
#include "globals.h"

// AY, include ------->
#include <netcdf.h>
#include <stdio.h>
#include <cstring>
#define ERRCODE 2
#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); exit(ERRCODE);}
// <-------- AY

WghmStateFile::WghmStateFile(size_t numCells, size_t lengthOfTimeSeries)
{  
  _cells.resize(numCells);
  
  for(uint i=0;i<_cells.size();i++)
    _cells.at(i) = Cell(i+1,lengthOfTimeSeries);
  
}

WghmStateFile::WghmStateFile(const WghmStateFile &f)
{
  _cells = f._cells;  
}

WghmStateFile &WghmStateFile::operator= (const WghmStateFile &f)
{
  
  _cells = f._cells;  
  
  return *this;
  
}
  
  
void WghmStateFile::load(std::string filename)
{
    // AY, add load from netcdf file -------->
    
    // determine extension of filename-string
    if (filename.substr(filename.find_last_of(".") + 1) == "nc")
    {
        int ncid, varid[12];
        int retval;
        int dimid[2];
        size_t dimensions[2];
        
        int linesIN;
        
        // open file
        if ((retval = nc_open(filename.c_str(), NC_NOWRITE, &ncid)))
            ERR(retval);
        
        // get id of the variables
        if ((retval = nc_inq_varid(ncid, "ID", &varid[0])))
            ERR(retval);
        if ((retval = nc_inq_varid(ncid, "TWS", &varid[1])))
            ERR(retval);
        if ((retval = nc_inq_varid(ncid, "CANOPY", &varid[2])))
            ERR(retval);
        if ((retval = nc_inq_varid(ncid, "SNOW", &varid[3])))
            ERR(retval);
        if ((retval = nc_inq_varid(ncid, "SOIL", &varid[4])))
            ERR(retval);
        if ((retval = nc_inq_varid(ncid, "LOCALLAKE", &varid[5])))
            ERR(retval);
        if ((retval = nc_inq_varid(ncid, "LOCALWETLAND", &varid[6])))
            ERR(retval);
        if ((retval = nc_inq_varid(ncid, "GLOBALLAKE", &varid[7])))
            ERR(retval);
        if ((retval = nc_inq_varid(ncid, "GLOBALWETLAND", &varid[8])))
            ERR(retval);
        if ((retval = nc_inq_varid(ncid, "RESERVOIR", &varid[9])))
            ERR(retval);
        if ((retval = nc_inq_varid(ncid, "RIVER", &varid[10])))
            ERR(retval);
        if ((retval = nc_inq_varid(ncid, "GROUNDWATER", &varid[11])))
            ERR(retval);
        
        // get size of the variable and create corresponding data containers
        if ((retval = nc_inq_vardimid(ncid, varid[0], dimid)))
            ERR(retval);
        
        if ((retval = nc_inq_dimlen(ncid, dimid[0], &dimensions[0])))
            ERR(retval);
        
        if ((retval = nc_inq_dimlen(ncid, dimid[1], &dimensions[1])))
            ERR(retval);
        
        linesIN = dimensions[0];
        int data_id[linesIN];
        double data_tws[linesIN], data_canopy[linesIN],
        data_snow[linesIN], data_soil[linesIN], data_locallake[linesIN],
        data_localwetland[linesIN], data_globallake[linesIN], data_globalwetland[linesIN], data_reservoir[linesIN], data_river[linesIN], data_groundwater[linesIN];

        // read data
        if ((retval = nc_get_var_int(ncid, varid[0], &data_id[0])))
            ERR(retval);
        if ((retval = nc_get_var_double(ncid, varid[1], &data_tws[0])))
            ERR(retval);
        if ((retval = nc_get_var_double(ncid, varid[2], &data_canopy[0])))
            ERR(retval);
        if ((retval = nc_get_var_double(ncid, varid[3], &data_snow[0])))
            ERR(retval);
        if ((retval = nc_get_var_double(ncid, varid[4], &data_soil[0])))
            ERR(retval);
        if ((retval = nc_get_var_double(ncid, varid[5], &data_locallake[0]))) 
            ERR(retval);
        if ((retval = nc_get_var_double(ncid, varid[6], &data_localwetland[0]))) //AY: data_globalwetland
            ERR(retval);
        if ((retval = nc_get_var_double(ncid, varid[7], &data_globallake[0]))) //AY: data_localwetland
            ERR(retval);
        if ((retval = nc_get_var_double(ncid, varid[8], &data_globalwetland[0]))) //AY: data_globallake
            ERR(retval);
        if ((retval = nc_get_var_double(ncid, varid[9], &data_reservoir[0])))
            ERR(retval);
        if ((retval = nc_get_var_double(ncid, varid[10], &data_river[0]))) 
            ERR(retval);
        if ((retval = nc_get_var_double(ncid, varid[11], &data_groundwater[0])))
            ERR(retval);

        // close file
        if ((retval = nc_close(ncid)))
            ERR(retval);
        
        // assign to cells
        std::vector<Cell> cellsFromFile;
        
        for(int i=0;i<linesIN;i++)
        {
            Cell cell(i+1,1);
            cell.canopy(0) = data_canopy[i];
            cell.snow(0) = data_snow[i];
            cell.soil(0) = data_soil[i];
            cell.locallake(0) = data_locallake[i];
            cell.localwetland(0) = data_localwetland[i];
            cell.globallake(0) = data_globallake[i];
            cell.globalwetland(0) = data_globalwetland[i];
            cell.reservoir(0) = data_reservoir[i];
            cell.river(0) = data_river[i];
            cell.groundwater(0) = data_groundwater[i];
            
            cellsFromFile.push_back(cell);
        }
        
        if(_cells.size()==0)
            _cells = cellsFromFile;
        else
        {
            for(uint i=0;i<cellsFromFile.size();i++)
                _cells.at(cellsFromFile.at(i).id()-1) = cellsFromFile.at(i);
        }
        
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
  
  std::vector<Cell> cellsFromFile;
	
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

      Cell cell(id,1);
      double tws;
      ss>>tws;
      ss>>cell.canopy(0);
      ss>>cell.snow(0);
      ss>>cell.soil(0);
      ss>>cell.locallake(0);
      ss>>cell.localwetland(0);
      ss>>cell.globallake(0);
      ss>>cell.globalwetland(0);
      ss>>cell.reservoir(0);
      ss>>cell.river(0);
      ss>>cell.groundwater(0);
           
      cellsFromFile.push_back(cell);

    }
  }
  
	if(_cells.size()==0)
    _cells = cellsFromFile;
  else
  {
    for(uint i=0;i<cellsFromFile.size();i++)
      _cells.at(cellsFromFile.at(i).id()-1) = cellsFromFile.at(i);
  }
  
}
catch(std::exception &e)
{
  throw(Exception(std::string("In WghmStateFile::load():\n")+e.what()));
}
        } // AY, don't forget to close the if-loop
}

void WghmStateFile::saveMean(std::string filename) //,WghmStateFile wghmMean
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
            
            // set size of the nc file
            linesIN = (int) _cells.size();
            colsIN=12;
            
            char varname[12][15] = {"ID","TWS","CANOPY","SNOW","SOIL",
                "LOCALLAKE","LOCALWETLAND","GLOBALLAKE","GLOBALWETLAND","RESERVOIR",
                "RIVER","GROUNDWATER"};
            
            int varid[12];
            
            int data_id[linesIN];
            double data_tws[linesIN], data_canopy[linesIN],
            data_snow[linesIN], data_soil[linesIN], data_locallake[linesIN],
            data_localwetland[linesIN], data_globallake[linesIN], data_globalwetland[linesIN], data_reservoir[linesIN], data_river[linesIN], data_groundwater[linesIN];
            
            shuffle = NC_SHUFFLE;
            deflate = 1;
            deflate_level = 1;
            
            for(int i=0;i<(int) _cells.size();i++)
            {
                Cell cell =  _cells.at(i).mean();
                
                data_id[i] = cell.id();
                data_tws[i] = cell.tws(0);
                data_canopy[i] = cell.canopy(0);
                data_snow[i] = cell.snow(0);
                data_soil[i] = cell.soil(0);
                data_locallake[i] = cell.locallake(0);
                data_localwetland[i] = cell.localwetland(0);
                data_globallake[i] = cell.globallake(0);
                data_globalwetland[i] = cell.globalwetland(0);
                data_reservoir[i] = cell.reservoir(0);
                data_river[i] = cell.river(0);
                data_groundwater[i] = cell.groundwater(0) ;
            }
            
            // initialize the nc file
            if ((retval = nc_create(filename.c_str(), NC_NETCDF4, &ncid)))
                ERR(retval);
            
            // define dimensions
            if ((retval = nc_def_dim(ncid, "rows", linesIN, &row_dimid)))
                ERR(retval);
            if ((retval = nc_def_dim(ncid, "columns", 1, &col_dimid)))
                ERR(retval);
            
            dimids[0] = row_dimid;
            dimids[1] = col_dimid;
            chunks[0] = 180;
            chunks[1] = 1;
            
            // define variables
            for (int i=0; i<colsIN; i++)
            {
                if ((retval = nc_def_var(ncid, varname[i], NC_DOUBLE, NDIMS,
                                         dimids, &varid[i])))
                    ERR(retval);
            }
            // define compression
            for (int i=0; i<colsIN; i++)
            {
                if ((retval = nc_def_var_chunking(ncid, varid[i], 0, &chunks[0])))
                    ERR(retval);
            }
            
            for (int i=0; i<colsIN; i++)
            {
                if ((retval = nc_def_var_deflate(ncid, varid[i], shuffle, deflate, deflate_level)))
                    ERR(retval);
            }
            
            // assign values to the variables
            if ((retval = nc_put_var_int(ncid, varid[0], &data_id[0])))
                ERR(retval);
            if ((retval = nc_put_var_double(ncid, varid[1], &data_tws[0])))
                ERR(retval);
            if ((retval = nc_put_var_double(ncid, varid[2], &data_canopy[0])))
                ERR(retval);
            if ((retval = nc_put_var_double(ncid, varid[3], &data_snow[0])))
                ERR(retval);
            if ((retval = nc_put_var_double(ncid, varid[4], &data_soil[0])))
                ERR(retval);
            if ((retval = nc_put_var_double(ncid, varid[5], &data_locallake[0])))
                ERR(retval);
            if ((retval = nc_put_var_double(ncid, varid[6], &data_localwetland[0])))
                ERR(retval);
            if ((retval = nc_put_var_double(ncid, varid[7], &data_globallake[0])))
                ERR(retval);
            if ((retval = nc_put_var_double(ncid, varid[8], &data_globalwetland[0])))
                ERR(retval);
            if ((retval = nc_put_var_double(ncid, varid[9], &data_reservoir[0])))
                ERR(retval);
            if ((retval = nc_put_var_double(ncid, varid[10], &data_river[0])))
                ERR(retval);
            if ((retval = nc_put_var_double(ncid, varid[11], &data_groundwater[0])))
                ERR(retval);
            
            // close the nc file
            if ((retval = nc_close(ncid)))
                ERR(retval);
            
        }
        catch(std::exception &e)
        {
            throw(Exception(std::string("In WghmStateFile::saveDay():\n")+e.what()));
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

  uint width = 33;
  
  // Header schreiben:
  stream << "WGHM water storage states" <<std::endl;
  stream << std::setw(6) << std::setfill(' ')     << "ID";
  stream << std::setw(width) << std::setfill(' ') <<"TWS";
  stream << std::setw(width) << std::setfill(' ') <<"CANOPY";
  stream << std::setw(width) << std::setfill(' ') <<"SNOW";
  stream << std::setw(width) << std::setfill(' ') <<"SOIL";
  stream << std::setw(width) << std::setfill(' ') <<"LOCALLAKE";
  stream << std::setw(width) << std::setfill(' ') <<"LOCALWETLAND";
  stream << std::setw(width) << std::setfill(' ') <<"GLOBALLAKE";
  stream << std::setw(width) << std::setfill(' ') <<"GLOBALWETLAND";
  stream << std::setw(width) << std::setfill(' ') <<"RESERVOIR";
  stream << std::setw(width) << std::setfill(' ') <<"RIVER";
  stream << std::setw(width) << std::setfill(' ') <<"GROUNDWATER";
  stream << std::endl;

  
  for(int i=0;i<(int) _cells.size();i++)
  { 
    Cell cell =  _cells.at(i).mean();
     
    stream.precision(16);
    stream << std::setw(6)     << std::setfill(' ') << cell.id();
    stream << std::setw(width) << std::scientific << std::setfill(' ') << cell.tws(0);
    stream << std::setw(width) << std::scientific << std::setfill(' ') << cell.canopy(0);
    stream << std::setw(width) << std::scientific << std::setfill(' ') << cell.snow(0);
    stream << std::setw(width) << std::scientific << std::setfill(' ') << cell.soil(0);
    stream << std::setw(width) << std::scientific << std::setfill(' ') << cell.locallake(0);
    stream << std::setw(width) << std::scientific << std::setfill(' ') << cell.localwetland(0);
    stream << std::setw(width) << std::scientific << std::setfill(' ') << cell.globallake(0);
    stream << std::setw(width) << std::scientific << std::setfill(' ') << cell.globalwetland(0);
    stream << std::setw(width) << std::scientific << std::setfill(' ') << cell.reservoir(0);
    stream << std::setw(width) << std::scientific << std::setfill(' ') << cell.river(0);
    stream << std::setw(width) << std::scientific << std::setfill(' ') << cell.groundwater(0) << std::endl;

  }
  
  stream.close();
}
catch(std::exception &e)
{
  throw(Exception(std::string("In WghmStateFile::saveMean():\n")+e.what()));
}
        } // AY, don't forget to close the if-loop
}

void WghmStateFile::saveDay(std::string filename, int day)//,WghmStateFile wghmMean
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
            
            // set size of the nc file
            linesIN = (int) _cells.size();
            colsIN=12;
            
            char varname[12][15] = {"ID","TWS","CANOPY","SNOW","SOIL",
                "LOCALLAKE","LOCALWETLAND","GLOBALLAKE","GLOBALWETLAND","RESERVOIR",
                "RIVER","GROUNDWATER"};
            
            int varid[12];
            
            int data_id[linesIN];
            double data_tws[linesIN], data_canopy[linesIN],
            data_snow[linesIN], data_soil[linesIN], data_locallake[linesIN],
            data_localwetland[linesIN], data_globallake[linesIN], data_globalwetland[linesIN], data_reservoir[linesIN], data_river[linesIN], data_groundwater[linesIN];
            
            shuffle = NC_SHUFFLE;
            deflate = 1;
            deflate_level = 1;
            
            
            for(int i=0;i<(int) _cells.size();i++)
            {
                Cell cell = _cells.at(i);
                
                data_id[i] = cell.id();
                data_tws[i] = cell.tws(day);
                data_canopy[i] = cell.canopy(day);
                data_snow[i] = cell.snow(day);
                data_soil[i] = cell.soil(day);
                data_locallake[i] = cell.locallake(day);
                data_localwetland[i] = cell.localwetland(day);
                data_globallake[i] = cell.globallake(day);
                data_globalwetland[i] = cell.globalwetland(day);
                data_reservoir[i] = cell.reservoir(day);
                data_river[i] = cell.river(day);
                data_groundwater[i] = cell.groundwater(day) ;
            }
            
            // initialize the nc file
            // note: c-style string for netcdf
            if ((retval = nc_create(filename.c_str(), NC_NETCDF4, &ncid)))
                ERR(retval);
            
            // define dimensions
            if ((retval = nc_def_dim(ncid, "rows", linesIN, &row_dimid)))
                ERR(retval);
            if ((retval = nc_def_dim(ncid, "columns", 1, &col_dimid)))
                ERR(retval);
            
            dimids[0] = row_dimid;
            dimids[1] = col_dimid;
            chunks[0] = 180;
            chunks[1] = 1;
            
            // define variables
            for (int i=0; i<colsIN; i++)
            {
                if ((retval = nc_def_var(ncid, varname[i], NC_DOUBLE, NDIMS,
                                         dimids, &varid[i])))
                    ERR(retval);
            }
            // define compression
            for (int i=0; i<colsIN; i++)
            {
                if ((retval = nc_def_var_chunking(ncid, varid[i], 0, &chunks[0])))
                    ERR(retval);
            }
            
            for (int i=0; i<colsIN; i++)
            {
                if ((retval = nc_def_var_deflate(ncid, varid[i], shuffle, deflate, deflate_level)))
                    ERR(retval);
            }
            
            // assign values to the variables
            if ((retval = nc_put_var_int(ncid, varid[0], &data_id[0])))
                ERR(retval);
            if ((retval = nc_put_var_double(ncid, varid[1], &data_tws[0])))
                ERR(retval);
            if ((retval = nc_put_var_double(ncid, varid[2], &data_canopy[0])))
                ERR(retval);
            if ((retval = nc_put_var_double(ncid, varid[3], &data_snow[0])))
                ERR(retval);
            if ((retval = nc_put_var_double(ncid, varid[4], &data_soil[0])))
                ERR(retval);
            if ((retval = nc_put_var_double(ncid, varid[5], &data_locallake[0])))
                ERR(retval);
            if ((retval = nc_put_var_double(ncid, varid[6], &data_localwetland[0])))
                ERR(retval);
            if ((retval = nc_put_var_double(ncid, varid[7], &data_globallake[0])))
                ERR(retval);
            if ((retval = nc_put_var_double(ncid, varid[8], &data_globalwetland[0])))
                ERR(retval);
            if ((retval = nc_put_var_double(ncid, varid[9], &data_reservoir[0])))
                ERR(retval);
            if ((retval = nc_put_var_double(ncid, varid[10], &data_river[0])))
                ERR(retval);
            if ((retval = nc_put_var_double(ncid, varid[11], &data_groundwater[0])))
                ERR(retval);
            
            // close the nc file
            if ((retval = nc_close(ncid)))
                ERR(retval);
            
            printf("*** writing complete!\n");
        }
        catch(std::exception &e)
        {
            throw(Exception(std::string("In WghmStateFile::saveDay():\n")+e.what()));
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

  uint width = 33;
  
  // Header schreiben:
  stream << "WGHM water storage states" <<std::endl;
  stream << std::setw(6) << std::setfill(' ')     << "ID";
  stream << std::setw(width) << std::setfill(' ') <<"TWS";
  stream << std::setw(width) << std::setfill(' ') <<"CANOPY";
  stream << std::setw(width) << std::setfill(' ') <<"SNOW";
  stream << std::setw(width) << std::setfill(' ') <<"SOIL";
  stream << std::setw(width) << std::setfill(' ') <<"LOCALLAKE";
  stream << std::setw(width) << std::setfill(' ') <<"LOCALWETLAND";
  stream << std::setw(width) << std::setfill(' ') <<"GLOBALLAKE";
  stream << std::setw(width) << std::setfill(' ') <<"GLOBALWETLAND";
  stream << std::setw(width) << std::setfill(' ') <<"RESERVOIR";
  stream << std::setw(width) << std::setfill(' ') <<"RIVER";
  stream << std::setw(width) << std::setfill(' ') <<"GROUNDWATER";
  stream << std::endl;
  
  for(int i=0;i<(int) _cells.size();i++)
  { 
    Cell cell = _cells.at(i);
      
    stream.precision(16);
    stream << std::setw(6)     << std::setfill(' ') << cell.id();
    stream << std::setw(width) << std::scientific << std::setfill(' ') << cell.tws(day);
    stream << std::setw(width) << std::scientific << std::setfill(' ') << cell.canopy(day);
    stream << std::setw(width) << std::scientific << std::setfill(' ') << cell.snow(day);
    stream << std::setw(width) << std::scientific << std::setfill(' ') << cell.soil(day);
    stream << std::setw(width) << std::scientific << std::setfill(' ') << cell.locallake(day);
    stream << std::setw(width) << std::scientific << std::setfill(' ') << cell.localwetland(day);
    stream << std::setw(width) << std::scientific << std::setfill(' ') << cell.globallake(day);
    stream << std::setw(width) << std::scientific << std::setfill(' ') << cell.globalwetland(day);
    stream << std::setw(width) << std::scientific << std::setfill(' ') << cell.reservoir(day);
    stream << std::setw(width) << std::scientific << std::setfill(' ') << cell.river(day);
    stream << std::setw(width) << std::scientific << std::setfill(' ') << cell.groundwater(day) << std::endl;
    
  }
  
  stream.close();
}
catch(std::exception &e)
{
  throw(Exception(std::string("In WghmStateFile::saveDay():\n")+e.what()));
}
        } // AY, don't forget to close the if-loop
}

void WghmStateFile::resizeCells(size_t l)
{
  for(int i=0;i<(int) _cells.size();i++)
    _cells.at(i).resize(l);  
}

void WghmStateFile::resetCells(size_t l)
{
  for(int i=0;i<(int) _cells.size();i++)
  {
    _cells.at(i).resize(0);  
    _cells.at(i).resize(l);  
  }
}
void WghmStateFile::lastDay_inplace()
{
    for(int i=0;i<(int) _cells.size();i++)
    {
      _cells.at(i).lastDay_remains();
    }
    
}

void WghmStateFile::mean_inplace()
{
    for(int i=0;i<(int) _cells.size();i++)
    {
        _cells.at(i)=_cells.at(i).mean();
        
    }
    
}

Cell::Cell()
{
  _cellID = 0;
  resize(0);
}

Cell::Cell(size_t cellID, size_t l)
{
  _cellID = cellID;
  resize(l);
}

Cell::Cell(const Cell &c)
{
  _cellID = c._cellID;
  
  _canopy = c._canopy;
  _snow = c._snow;
  _soil = c._soil;
  _locallake = c._locallake;
  _localwetland = c._localwetland;
  _globallake = c._globallake;
  _globalwetland = c._globalwetland;
  _reservoir = c._reservoir;
  _river = c._river;
  _groundwater = c._groundwater;
}

Cell &Cell::operator= (const Cell &c)
{
  _cellID = c._cellID;
  
  _canopy 	= c._canopy;
  _snow 	= c._snow;
  _soil 	= c._soil;
  _locallake 	= c._locallake;
  _localwetland = c._localwetland;
  _globallake 	= c._globallake;
  _globalwetland= c._globalwetland;
  _reservoir 	= c._reservoir;
  _river 	= c._river;
  _groundwater 	= c._groundwater;

  return *this;
}

void Cell::resize(size_t l)
{
  _canopy.resize(l);
  _snow.resize(l);
  _soil.resize(l);
  _locallake.resize(l);
  _localwetland.resize(l);
  _globallake.resize(l);
  _globalwetland.resize(l);
  _reservoir.resize(l);
  _river.resize(l);
  _groundwater.resize(l);

}

void Cell::lastDay_remains()
{
  _canopy.erase(_canopy.begin(),_canopy.begin()+_canopy.size()-1);
  _snow.erase(_snow.begin(),_snow.begin()+_snow.size()-1);
  _soil.erase(_soil.begin(),_soil.begin()+_soil.size()-1);
  _locallake.erase(_locallake.begin(),_locallake.begin()+_locallake.size()-1);
  _localwetland.erase(_localwetland.begin(),_localwetland.begin()+_localwetland.size()-1);
  _globallake.erase(_globallake.begin(),_globallake.begin()+_globallake.size()-1);
  _globalwetland.erase(_globalwetland.begin(),_globalwetland.begin()+_globalwetland.size()-1);
  _reservoir.erase(_reservoir.begin(),_reservoir.begin()+_reservoir.size()-1);
  _river.erase(_river.begin(),_river.begin()+_river.size()-1);
  _groundwater.erase(_groundwater.begin(),_groundwater.begin()+_groundwater.size()-1);

}


Cell Cell::mean()
{
  Cell m(_cellID, 1);

  m.canopy(0) = std::accumulate(_canopy.begin(),_canopy.end(),0.0) / (double) _canopy.size();
  m.snow(0) = std::accumulate(_snow.begin(),_snow.end(),0.0) / (double) _snow.size();
  m.soil(0) = std::accumulate(_soil.begin(),_soil.end(),0.0) / (double) _soil.size();
  m.locallake(0) = std::accumulate(_locallake.begin(),_locallake.end(),0.0) / (double) _locallake.size();
  m.localwetland(0) = std::accumulate(_localwetland.begin(),_localwetland.end(),0.0) / (double) _localwetland.size();
  m.globallake(0) = std::accumulate(_globallake.begin(),_globallake.end(),0.0) / (double) _globallake.size();
  m.globalwetland(0) = std::accumulate(_globalwetland.begin(),_globalwetland.end(),0.0) / (double) _globalwetland.size();
  m.reservoir(0) = std::accumulate(_reservoir.begin(),_reservoir.end(),0.0) / (double) _reservoir.size();
  m.river(0) = std::accumulate(_river.begin(),_river.end(),0.0) / (double) _river.size();
  m.groundwater(0) = std::accumulate(_groundwater.begin(),_groundwater.end(),0.0) / (double) _groundwater.size();
  
  return m;
}


double Cell::tws(int i) //, WghmStateFile wghmMean
{
    double tws = canopy(i) + snow(i) + soil(i) + locallake(i) + localwetland(i) + globallake(i) + globalwetland(i) + reservoir(i) + river(i) + groundwater(i);//- wghmMean.canopy(i)-wghmMean.snow(i)-wghmMean.soil(i)-wghmMean.locallake(i)-wghmMean.localwetland(i)-wghmMean.globallake(i)-wghmMean.globalwetland(i)-wghmMean.reservoir(i)-wghmMean.river(i)-wghmMean.groundwater(i);
  
  return tws;
}
