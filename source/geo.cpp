#include <cstdio>
#include "gridio.h"
#include "option.h"
#include "geo.h"
#include "calibration.h"

using namespace std;

// FP20161018N003 Enable reading input data from climate land mask
extern optionClass options;

void geoClass::init(const char *inputDir, const short fileEndianType)
{
	char filename[250];
	gridioClass gridIO;

	gridIO.setFileEndianType(fileEndianType);

	// read file with area of cells in a certain row  [km2]
	// file contains 360 lines, two for each degree of latitude
	// according to the cell size of about 0.5 degree
	// the value for the grid cell area are given in km2
	sprintf(filename, "%s/GAREA.UNF0", inputDir);
	gridIO.readUnfFile(filename, 360, area);

	// read file which gives number of row for each gridcell
	sprintf(filename, "%s/GR.UNF2", inputDir);
	gridIO.readUnfFile(filename, ng, G_row);
	sprintf(filename, "%s/GC.UNF2", inputDir);
	gridIO.readUnfFile(filename, ng, G_col);
	sprintf(filename, "%s/GCRC.UNF4", inputDir);
	gridIO.readUnfFile(filename, 360 * 720, &GCRC[0][0]);
    cout << "Reading land and (continental) water fractions and generating grid with continental fractions" << endl; // FP20161010N001
	sprintf(filename, "%s/GFREQ.UNF0", inputDir);
    gridIO.readUnfFile(filename, ng, G_landfreq_const); // FP20161018N002
	sprintf(filename, "%s/GFREQW.UNF0", inputDir);
    gridIO.readUnfFile(filename, ng, G_fwaterfreq_const); // FP20161018N002

    // Define indicator mask of continental (1) vs. ocean cells (0)
    long n_cells_zero = 0;
    cout << "Calculate continental percentages & Evaluate number of grid cells with zero continental area (Caspian Sea in Standard Land Mask)" << endl;
    for (int n = 0; n < ng; n++) {
        // Continental area as percentage of grid cell area from land and freshwater percentages
        G_contfreq[n] = G_landfreq_const[n] + G_fwaterfreq_const[n]; // FP20161010N001 // FP20161018N002
        // ocean cells (Caspian Sea in Standard Land Mask): fraction is lower than limit_contcell_pct
        if (G_contfreq[n] < limit_contcell_pct) {
            G_contcell[n] = 0;
            //printf("(this) G_contcell[n] zero %d at grid cell n: %d, geo.G_landfreq[n]: %e, geo.G_fwaterfreq[n], geo.G_contfreq[n]: %e \n", this->G_contcell[n], n, this->G_landfreq[n], this->G_fwaterfreq[n], this->G_contfreq[n]);
            //printf("(...) G_contcell[n] zero %d at grid cell n: %d, geo.G_landfreq_const[n]: %e, geo.G_fwaterfreq_const[n]: %e, geo.G_contfreq[n]: %e \n", G_contcell[n], n, G_landfreq_const[n], G_fwaterfreq_const[n], G_contfreq[n]); // FP20161018N002
            n_cells_zero += 1;
        }
        // continental cells
        else {
            G_contcell[n] = 1;
        }
    }
    printf("Number of cells with zero percentage of continental area: %d, limit %e\n", n_cells_zero, limit_contcell_pct);

    // FP20161018N003 Enable reading input data from climate land mask
    // read information for climate conversion
    if (1 == options.climate_spatial_resolution) {
        G_cellnum_clm = new int[ng];
        sprintf(filename,"%s/G_CELLNUM_CLM.UNF4",inputDir);
        gridIO.readUnfFile(filename, ng, G_cellnum_clm);
    }

#ifdef _WATERGAP_CHECKS_GLOBAL_H
        // Check Felix Portmann 2015 - for invalid data
        // Check iteration sequence
        cout << "geo.cpp - Check possible iteration sequence n++ vs. ++n" << endl;
        for (int n = 0; n < 10; n++) {
            printf("n++ iteration for n: %d\n", n);
        }
        for (int n = 0; n < 10; ++n) {
            printf("++n iteration for n: %d\n", n);
        }
#endif

}
// end init


// Felix Portmann 2015 & 2016
// Class to write grid of mask of continental cells
// FP20161010N002: remove internal selection "if (calibrun == 0)"; now selection before calling "if (geo.calibrun == 0)"
void geoClass::writeContCellMask(char *output_dir)
{
	// external
	extern gridioClass gridIO;
	 // local
	char filename[250];
		sprintf(filename, "%s/G_CONT_CELL_MASK.UNF2", output_dir);
		cout << "writing mask of continental cells (1) vs. ocean cells (0) to file " << filename << endl;
		gridIO.writeUnfFile(filename, ng, &G_contcell[0]);
}
// end writeContCellMask

// FP20161010N001
// Class to write grid of percentage of continental area
void geoClass::writeContFreq(char *output_dir)
{
    // external
    extern gridioClass gridIO;
    // local
    char filename[250];
        sprintf(filename, "%s/GCONTFREQ.UNF0", output_dir);
        cout << "writing continental share/area as percentage of grid cell area to file " << filename << endl;
        gridIO.writeUnfFile(filename, ng, &G_contfreq[0]);
}
// end writeContFreq

// Felix Portmann 2015 & 2016
// Class to write grids of characteristic basic areas of grid cells (data type: float)
// FP20161010N002: remove internal selection "if (calibrun == 0)"; now selection before calling "if (geo.calibrun == 0)"
void geoClass::writeGridAreas(char *output_dir)
{
    // external
    extern float G_float[ng];
    extern gridioClass gridIO;

     // local
    int n;
    char filename[250];

    // Initialize output grid
    for (n = 0; n < ng; n++) {
        G_float[n] = 0.0;
    }

    // Write: Cell area [km2]
    for (n = 0; n < ng; n++) {
        G_float[n] = areaOfCellByArrayPos(n);
    }
		sprintf(filename, "%s/G_CELL_AREA_OUT_km2.UNF0", output_dir);
		cout << "writing cell area in km2 (no transformation) to file " << filename << endl;
		cout << "= geo.areaOfCellByArrayPos(n)" << endl;
		gridIO.writeUnfFile(filename, ng, &G_float[0]);

		// Write: Continental area (including surface water and land, but no ocean)  [km2]
        for (n = 0; n < ng; n++) {
            G_float[n] = areaOfCellByArrayPos(n) * ( G_contfreq[n] / 100.0 ); // FP20161010N001
        }
		sprintf(filename, "%s/G_CONT_AREA_OUT_km2.UNF0", output_dir);
		cout << "writing continental area in km2 from percent fractions to file " << filename << endl;
        cout << "= geo.areaOfCellByArrayPos(n) * ( G_contfreq[0] / 100.0 ); with G_contfreq[0] = (geo.G_landfreq[n] + geo.G_fwaterfreq[n]) " << endl; // FP20161010N001
		gridIO.writeUnfFile(filename, ng, &G_float[0]);
}
// end writeGridAreas


int geoClass::cellNumByLonLat(const float longitude, const float latitude)
{
	short row, column;
	int cellNumber;

	if ((longitude < -180) || (longitude > 180)) {
		cerr << "Longitude out of range" << endl;
		return -99;
	}
	if ((latitude < -90) || (latitude > 90)) {
		cerr << "Latitude out of range" << endl;
		return -99;
	}

	column = (short) ((longitude + 180) * 2 + 1);
	row = (short) ((-latitude + 90) * 2 + 1);
	// starting from value 1 
	// according to definition in IMAGE (FORTRAN-arrays)
	// (e.g. in file GCRC.UNF4).
	// for C/C++ first value of arrays is zero!!

	cellNumber = GCRC[column - 1][row - 1];
	if (0 == cellNumber) {
		cerr << "Coordinate: "
			<< longitude << "/" << latitude << " is not part of the landmask" << endl;
		return -99;
	} else
		return cellNumber;
}
// end cellNumByLonLat
