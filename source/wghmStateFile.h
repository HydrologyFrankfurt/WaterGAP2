#if !defined (_wghmstatefile_h_)
#define _wghmstatefile_h_

#include "common.h"
#include <numeric>

class Cell{
  public:
    
  Cell();
  Cell(size_t cellID, size_t l);
  Cell(const Cell &c);
  Cell &operator= (const Cell &c);
  
  void resize(size_t);
  double &canopy(int i){return (_canopy.at(i));}
  double &snow(int i){return (_snow.at(i));}
  double &soil(int i){return (_soil.at(i));}
  double &locallake(int i){return (_locallake.at(i));}
  double &localwetland(int i){return (_localwetland.at(i));}
  double &globallake(int i){return (_globallake.at(i));}
  double &globalwetland(int i){return (_globalwetland.at(i));}
  double &reservoir(int i){return (_reservoir.at(i));}
  double &river(int i){return (_river.at(i));}
  double &groundwater(int i){return (_groundwater.at(i));}
  
  void lastDay_remains();
  Cell mean();
  double tws(int); //, WghmStateFile wghmMean
  int id() {return _cellID;}
  
  size_t size(){return _canopy.size();}
  
  private:
    
    size_t _cellID;
    
    std::vector<double> _canopy;
    std::vector<double> _snow;
    std::vector<double> _soil;
    std::vector<double> _locallake;
    std::vector<double> _localwetland;
    std::vector<double> _globallake;
    std::vector<double> _globalwetland;
    std::vector<double> _reservoir;
    std::vector<double> _river;
    std::vector<double> _groundwater;
};

class WghmStateFile{ // total water storage of the grid cells
public:
	
  WghmStateFile(size_t numCells=0, size_t lengthOfTimeSeries=0);
  WghmStateFile(const WghmStateFile &f);
  WghmStateFile &operator= (const WghmStateFile &f);


  void load(std::string);
  void saveMean(std::string);//, WghmStateFile wghmMean
  void saveDay(std::string, int day);//, WghmStateFile wghmMean
  void lastDay_inplace();
  void mean_inplace();
  
  size_t size() {return _cells.size();}
  Cell &cell(int i) {return (_cells.at(i));}
  void resizeCells(size_t l);
  void resetCells(size_t l=0);
  

  private:

  std::vector<Cell> 	_cells;	// wghm grid cells
};

#endif
