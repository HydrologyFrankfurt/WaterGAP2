#if !defined (_geo_h_)
#define _geo_h_

#include "def.h"
#include "grid.h"
#include <string>

// Introducing G_contcell and limit_contcell_pct

class geoClass {
  public:
    void init(const std::string inputDir, const short fileEndianType);

    Grid<short> G_row;
    Grid<short> G_col;
    VariableChannelGrid<360,int,true,720> GCRC;

    // _const[n] to distinguish from yearly changed routing.G_landfreq[n]/G_fwaterfreq[n]
    Grid<> G_landfreq_const; // changed to [%]
    Grid<> G_fwaterfreq_const; // changed to [%]
    Grid<> G_contfreq; // instead of sum of G_landfreq and G_fwaterfreq, in %

    // Enable reading input data from climate land mask
    Grid<int> G_cellnum_clm;
    Grid<double,false,360> area;	// area of grid cells depending on row

    // continental cells indicator mask = Logical indicator whether cells
    // are continent (1) or ocean (0) cells (Caspian Sea in WaterGAP standard land mask)
    Grid<short> G_contcell;

    // minimum percentage of cell area for continental cells
    double limit_contcell_pct = 0.000001;
    short calibrun = 0;

    inline float areaOfCell(const int cellNumber) {
      return area[G_row[cellNumber - 1] - 1];
    }

    inline double areaOfCellByArrayPos(const int n) {
      return area[G_row[n] - 1];
    }

    //used to calculate (constant) calib basin land area
    template < class T >
    float landAreaOfObjectInGrid(T objectID, Grid<T> grid, int gridSize) {

      float objectArea = 0;

      for (int n = 0; n <= gridSize - 1; n++)
        if (objectID == grid[n]){
            objectArea += area[G_row[n] - 1] * ( (G_contfreq[n] - G_fwaterfreq_const[n])/ 100.0);
        }

      return objectArea;
    }

    //used to calculate (constant) calib basin land area
    template < class T >
    float nonOceanAreaOfObjectInGrid(T objectID, Grid<T> grid, int gridSize) {

      float objectArea = 0;

      for (int n = 0; n <= gridSize - 1; n++)
        if (objectID == grid[n])
                  objectArea += area[G_row[n] - 1] * (G_contfreq[n] / 100.0);
      return objectArea;
    }

	  int cellNumByLonLat(const float longitude, const float latitude);

    // Methods to write (1) continental cell mask (2) grid cell area (3) percentage of continental area
    void writeContCellMask(const std::string output_dir);
    void writeGridAreas(const std::string output_dir);
    void writeContFreq(const std::string output_dir);

};
#endif
