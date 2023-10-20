#if !defined (_additionalOutputInputFile_h_)
#define _additionalOutputInputFile_h_

#include "def.h"
#include <numeric>
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>
#include "exception.h"


class AdditionalOutputInputFile{
public:
    int additionalfilestatus;
  // constructor
  AdditionalOutputInputFile();
  
  // Gleichoperator
  
  double &additionalOutputInput(int i, int j){return (_additionalOutputInput[i][j]);}
  
  void load(std::string);
  void save(std::string);

  private:

    std::vector<std::vector<double> > _additionalOutputInput;
};

#endif
