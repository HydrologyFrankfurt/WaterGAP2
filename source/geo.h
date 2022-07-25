#if !defined (_geo_h_)
#define _geo_h_

#include "def.h"

// Felix Portmann 2015
// Introducing G_contcell and limit_contcell_pct
//void writeContCellMask(char *output_dir);
//void writeBasicAreas(char *output_dir);

class geoClass {
  public:
	void init(const char *inputDir, const short fileEndianType);
	short G_row[ng];
	short G_col[ng];
	int GCRC[720][360];
        // FP20161018N002 Reservoir operation start years: _const[n] to distinguish from yearly changed routing.G_landfreq[n]/G_fwaterfreq[n]
        double G_landfreq_const[ng]; // changed to [%]
        double G_fwaterfreq_const[ng]; // changed to [%]
        double G_contfreq[ng]; // instead of sum of G_landfreq and G_fwaterfreq, in % // FP20161010N001
        // FP20161018N003 Enable reading input data from climate land mask
        int *G_cellnum_clm;
        double area[360];	// area of grid cells depending on row
        short G_contcell[ng]; // continental cells indicator mask = Logical indicator whether cells are continent (1) or ocean (0) cells (Caspian Sea in WaterGAP standard land mask)
        // minimum percentage of cell area for continental cells
        double limit_contcell_pct = 0.000001;
        short calibrun = 0;

    inline float areaOfCell(const int cellNumber) {
		return area[G_row[cellNumber - 1] - 1];
	}

        // obviously not used anymore
        /*inline float landAreaOfCell(const int cellNumber) {
		return area[G_row[cellNumber - 1] - 1] * (G_landfreq[cellNumber - 1] / 100.0);
        }*/

        inline double areaOfCellByArrayPos(const int n) {
		return area[G_row[n] - 1];
	}

        // obviously not used anymore
        /*inline float landAreaOfCellByArrayPos(const int n) {
		return area[G_row[n] - 1] * (G_landfreq[n] / 100.0);
        }*/

        //used to calculate (constant) calib basin land area
    template < class T > float landAreaOfObjectInGrid(T objectID, T * grid, int gridSize) {

		float objectArea = 0;
        // FP20161018N002 Reservoir operation start years: G_landfreq_const[n] instead G_landfreq[n]
		for (int n = 0; n <= gridSize - 1; n++)
			if (objectID == grid[n])
                objectArea += area[G_row[n] - 1] * (G_landfreq_const[n] / 100.0);
		return objectArea;
	}

        //used to calculate (constant) calib basin land area
    template < class T > float nonOceanAreaOfObjectInGrid(T objectID, T * grid, int gridSize) {

		float objectArea = 0;

		for (int n = 0; n <= gridSize - 1; n++)
			if (objectID == grid[n])
                objectArea += area[G_row[n] - 1] * (G_contfreq[n] / 100.0); // FP20161010N001
		return objectArea;
	}

	int cellNumByLonLat(const float longitude, const float latitude);

    // Felix Portmann 2015 & 2016
    // Methods to write (1) continental cell mask (2) grid cell area (3) percentage of continental area
    void writeContCellMask(char *output_dir);
    void writeGridAreas(char *output_dir);
    void writeContFreq(char *output_dir); // FP20161010N001

};
#endif
