#if !defined (_snowInElevationFile_h_)
#define _snowInElevationFile_h_

#include "def.h"
#include <numeric>
#include "common.h"


class SnowInElevationFile{
public:

  // constructor
  SnowInElevationFile();
  
  // Gleichoperator
  
  double &snowInElevation(int i, int j){return (_snowInElevation[i][j]);}
  
  void load(std::string);
  void save(std::string);

  private:

    int _numLandCells;
    int _numLevels;
    //std::vector<float> _snowInElevation;	// [mm] array to store snow cover for 101 (?) seperate elevations within G_snow cell
    std::vector<std::vector<double> > _snowInElevation;
};

#endif
