#include <cstdio>
#include "common.h"
#include "globals.h"

using namespace std;

void geoClass::init(const string inputDir, const short fileEndianType)
{
	char filename[250];

	// read file with area of cells in a certain row  [km2]
	// file contains 360 lines, two for each degree of latitude
	// according to the cell size of about 0.5 degree
	// the value for the grid cell area are given in km2
    area.read(inputDir + "/GAREA.UNF0");

	// read file which gives number of row for each gridcell
    G_row.read(inputDir + "/GR.UNF2");

    G_col.read(inputDir + "/GC.UNF2");

    GCRC.read(inputDir + "/GCRC.UNF4");

    cout << "Reading contfreq and generating grid with water fractions" << endl;

    G_contfreq.read(inputDir + "/GCONTFREQ.UNF0");
    // generating water fractions (G_fwaterfreq_const) and making sure that G_landfreq (G_contfreq - G_fwaterfreq_const)
    // will be at least 0
    Grid<> loclakconst;
    loclakconst.read(inputDir + "/G_LOCLAK.UNF0");
    G_fwaterfreq_const += loclakconst;
    Grid<> glolakconst;
    glolakconst.read(inputDir + "/G_GLOLAK.UNF0");
    G_fwaterfreq_const += glolakconst;
    if (options.resOpt == 1) {
        Grid<> locresconst;
        locresconst.read(inputDir + "/G_LOCRES.UNF0");
        G_fwaterfreq_const += locresconst;
        Grid<> resconst;
        resconst.read(inputDir + "/G_RES/G_RES_FRAC.UNF0");
        G_fwaterfreq_const += resconst;
    }
    // Guarantee that G_fwaterfreq_const is >= G_contfreq
    for (int n = 0; n < ng; n++) {
        if (G_fwaterfreq_const[n] > G_contfreq[n]) G_contfreq[n] = G_fwaterfreq_const[n];
    }
    // Define indicator mask of continental (1) vs. ocean cells (0)
    long n_cells_zero = 0;
    cout << "Calculate continental percentages & Evaluate number of grid cells with zero continental area (Caspian Sea in Standard Land Mask)" << endl;
    for (int n = 0; n < ng; n++) {

        // ocean cells (Caspian Sea in Standard Land Mask): fraction is lower than limit_contcell_pct
        if (G_contfreq[n] < limit_contcell_pct) {
            G_contcell[n] = 0;
            n_cells_zero += 1;
        }
        // continental cells
        else {
            G_contcell[n] = 1;
        }
    }
    printf("Number of cells with zero percentage of continental area: %ld, limit %e\n", n_cells_zero, limit_contcell_pct);

    // read information for climate conversion
    if (1 == options.climate_spatial_resolution) {
        G_cellnum_clm.read(inputDir + "/G_CELLNUM_CLM.UNF4");
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


// Function to write grid of mask of continental cells
void geoClass::writeContCellMask(const string output_dir)
{
    cout << "writing mask of continental cells (1) vs. ocean cells (0) to file " << output_dir << "/G_CONT_CELL_MASK.UNF2" << endl;
    G_contcell.write(output_dir + "/G_CONT_CELL_MASK.UNF2");
}

// Function to write grid of percentage of continental area
void geoClass::writeContFreq(const string output_dir)
{
    cout << "writing continental share/area as percentage of grid cell area to file " << output_dir << "/GCONTFREQ.UNF0" << endl;
    G_contfreq.write(output_dir + "/GCONTFREQ.UNF0");
}

// Function to write grids of characteristic basic areas of grid cells (data type: float)
void geoClass::writeGridAreas(const string output_dir)
{
    extern Grid<float> G_float;

     // local
    int n;

    // Write: Cell area [km2]
    for (n = 0; n < ng; n++) {
        G_float[n] = areaOfCellByArrayPos(n);
    }

    cout << "writing cell area in km2 (no transformation) to file " << output_dir << "/G_CELL_AREA_OUT_km2.UNF0" << endl;
    cout << "= geo.areaOfCellByArrayPos(n)" << endl;
    G_float.write(output_dir + "/G_CELL_AREA_OUT_km2.UNF0");

    // Write: Continental area (including surface water and land, but no ocean)  [km2]
    for (n = 0; n < ng; n++) {
        G_float[n] = areaOfCellByArrayPos(n) * ( G_contfreq[n] / 100.0 );
    }

    cout << "writing continental area in km2 from percent fractions to file " << output_dir << "/G_CONT_AREA_OUT_km2.UNF0" << endl;
    cout << "= geo.areaOfCellByArrayPos(n) * ( G_contfreq[0] / 100.0 ); with G_contfreq[0] = (geo.G_landfreq[n] + geo.G_fwaterfreq[n]) " << endl;
    G_float.write(output_dir + "/G_CONT_AREA_OUT_km2.UNF0");
}


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

	cellNumber = GCRC(column - 1,row - 1);
	if (0 == cellNumber) {
		cerr << "Coordinate: "
			<< longitude << "/" << latitude << " is not part of the landmask" << endl;
		return -99;
	} else
		return cellNumber;
}
