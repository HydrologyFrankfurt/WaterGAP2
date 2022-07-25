
/***********************************************************************
*
* changed G_freq, G_w_freq to float by Stephanie Eisner, read in GFREQ and
* GFREQW as UNF0 (input data created by Linda Adam (IPG), Stephanie Eisner
* (CESR))
*
* because G_mean_inflow (reservoir-algorithm) is in fact mean inflow
* plus P minus PET, this was renamed to G_mean_outflow (by Hannes Müller
* Schmied) in all cases.
*
* see former changes at file rout_prepare.cpp.versioninfos.txt
*
***********************************************************************/
#include <iostream>
#include <cstdio>
#include <cmath>
#include "gridio.h"
#include "timestring.h"

#include "def.h"

using namespace std;

extern gridioClass gridIO;

#define nb_max 60000	// max. basin number for declaration of array
#define nb2_max 10000	// maximum number of basins with more than one cell

// #define debug

void replace_flow_to_ocean(int GCRC_2[722][362], short *G_row, short *G_column, char *G_LDD);

void inflow_cells(char *G_LDD,
				  int GCRC_2[722][362],
				  short *G_row, short *G_column, int G_inflow_cells[ng][9], char *output_dir);


void flow_accumulation(int G_inflow_cells[ng][9], short *G_flow_acc);

unsigned short derive_basins(char *G_LDD,
							 int G_inflow_cells[ng][9],
							 unsigned short *G_watersheds, char *output_dir);

unsigned short reindex_waterbasins(unsigned short nb,
								   short *basin_cells, unsigned short *G_watersheds);

void outflow_cell(char *G_LDD,
				  int GCRC_2[722][362],
				  short *G_row, short *G_column, int *G_outflow_cell, char *output_dir);

void neighbouring_cells(int GCRC_2[722][362], short *G_row, short *G_column, int G_neighbourCell[ng][8], char *output_dir);

void calculate_distances(float cell_distance[9][360], char *output_dir);

void calculate_river_slope(float *G_meandering_ratio, float cell_distance[9][360], //LA (2014-02-12): meandering ratio added
						   float *G_altitude,
						   int *G_outflow_cell, char *G_LDD, short *G_row, char *output_dir);

void calculate_river_length(float *G_meandering_ratio, short *G_row, short *G_column, float cell_distance[9][360], 
                            int *G_outflow_cell, char *output_dir);

void rout_order(int *G_routOrder, int G_inflow_cells[ng][9], int *G_outflow_cell, char *output_dir);

void reservoir_prepare(int *G_outflow_cell, char *input_dir, char *output_dir, short *G_row);


void prepare_routing_files(char *input_dir, char *output_dir, short arcNumbersForFlowDir, short resOpt)
{
	// read binary input files
	char filename[250];

	char G_LDD[ng];

	if (1 == arcNumbersForFlowDir) {
		// if file with flow directon contains number 
		// according to the Arc/Info conventions 
		// then convert them to the internal representation 
		// of WaterGAP

		// char *G_flowdir;
		// G_flowdir = new char[ng];
		short *G_flowdir;
		G_flowdir = new short[ng];

//		sprintf(filename, "%s/G_FLOWDIR.UNF1", input_dir);
		sprintf(filename, "%s/G_FLOWDIR.UNF2", input_dir);
		gridIO.readUnfFile(filename, ng, G_flowdir);
		for (int n = 0; n < ng; ++n) {
			switch (G_flowdir[n]) {
			case -1:// inland sinks
				G_LDD[n] = -1;
				break;
			case 8:	// SW
				G_LDD[n] = 1;
				break;
			case 4:	// S
				G_LDD[n] = 2;
				break;
			case 2:	// SE
				G_LDD[n] = 3;
				break;
			case 16:	// W
				G_LDD[n] = 4;
				break;
			case 0:	// sink
				G_LDD[n] = 5;
				break;
			case 1:	// E
				G_LDD[n] = 6;
				break;
			case 32:	// NW
				G_LDD[n] = 7;
				break;
			case 64:	// N
				G_LDD[n] = 8;
				break;
			case 128:	// NE
				G_LDD[n] = 9;
				break;
			default:	// all other cases and nodata cells 
				G_LDD[n] = 5;	// are treated as sinks
				break;
			}
		}
	} else {
		sprintf(filename, "%s/G_LDD.UNF1", input_dir);
		gridIO.readUnfFile(filename, ng, &G_LDD[0]);
	}

	short G_row[ng], G_column[ng];

	sprintf(filename, "%s/GR.UNF2", input_dir);
	gridIO.readUnfFile(filename, ng, &G_row[0]);
	sprintf(filename, "%s/GC.UNF2", input_dir);
	gridIO.readUnfFile(filename, ng, &G_column[0]);
	
	int GCRC[720][360];

	sprintf(filename, "%s/GCRC.UNF4", input_dir);
	gridIO.readUnfFile(filename, 360 * 720, &GCRC[0][0]);

	float G_freq[ng], G_w_freq[ng];

	sprintf(filename, "%s/GFREQ.UNF0", input_dir);
	gridIO.readUnfFile(filename, ng, &G_freq[0]);
	sprintf(filename, "%s/GFREQW.UNF0", input_dir);
	gridIO.readUnfFile(filename, ng, &G_w_freq[0]);
	float area[360];

	sprintf(filename, "%s/GAREA.UNF0", input_dir);
	gridIO.readUnfFile(filename, 360, &area[0]);

	// copy array, and add a row/column at the borders of the map 
	// so that easy access to the neighbouring cells is possible  
	int GCRC_2[722][362];

	for (short j = 0; j < 722; ++j) {
		GCRC_2[j][0] = 0;
		GCRC_2[j][361] = 0;
	}
	for (short i = 0; i <= 359; i++) {
		GCRC_2[0][i + 1] = GCRC[719][i];
		GCRC_2[721][i + 1] = GCRC[0][i];
		for (short j = 0; j <= 719; j++)
			GCRC_2[j + 1][i + 1] = GCRC[j][i];
	}

	// give value of 5 to those cells which contain ocean water 
	//for (n=0; n<=(ng-1); n++)
	//  if (G_freq[n] + G_w_freq[n] < 36) G_LDD[n] = 5;

	//replace_flow_to_ocean(GCRC_2, G_row, G_column, G_LDD);

	sprintf(filename, "%s/G_LDD_2.UNF1", output_dir);
	gridIO.writeUnfFile(filename, ng, &G_LDD[0]);

	// cells from which flow into the cell of interest occurs 
	int G_inflow_cells[ng][9];

	// new information saved in G_inflow_cells
	inflow_cells(G_LDD, GCRC_2, G_row, G_column, G_inflow_cells, output_dir);

	// number of upstream cells (plus the cell itself) 
	short G_flow_acc[ng];

	flow_accumulation(G_inflow_cells, G_flow_acc);

	//sprintf(filename,"%s/G_FLOW_ACC_UNCORRECTED.UNF2",output_dir);
	//gridIO.writeUnfFile(filename, ng, G_flow_acc);

	// cells with G_flow_acc == 0 did not get an flow accumulation value.
	// the reason for this is that they are part of a loop.
	// for simplicity it is assumed that these cells are pits.
	// therefore they get a LDD-value of 5.
	short int correction = 0;

	for (int n = 0; n < ng; ++n)
		if (0 == G_flow_acc[n]) {
			G_LDD[n] = 5;
			G_flow_acc[n] = 1;
			correction = 1;
		}
	sprintf(filename, "%s/G_FLOW_ACC.UNF2", output_dir);
	gridIO.writeUnfFile(filename, ng, G_flow_acc);

	if (1 == correction) {
		// Calculation of 'inflow cells' and 'flow accumulation' has to be 
		// done again, if cells have been corrected
		inflow_cells(G_LDD, GCRC_2, G_row, G_column, G_inflow_cells, output_dir);
		flow_accumulation(G_inflow_cells, G_flow_acc);
	}


	unsigned short G_watersheds[ng];
	unsigned short nb;	// number of waterbasins 

	nb = derive_basins(G_LDD, G_inflow_cells, G_watersheds, output_dir);

	// calculate number of cells per waterbasins  
	short basin_cells[nb_max];

	for (int i = 0; i < nb_max; ++i)
		basin_cells[i] = 0;
	for (int n = 0; n <ng; ++n)
		basin_cells[G_watersheds[n]]++;
	//for (i=0; i<=nb; i++)
	// printf("%d %d\n",i,basin_cells[i]);

	// give new number to watersheds:  
	// only those are considered which 
	// have more than one cell.        
	// for the other no routing has to 
	// be performed.                   
	nb = reindex_waterbasins(nb, basin_cells, G_watersheds);

	sprintf(filename, "%s/G_BASINS_2.UNF2", output_dir);
	gridIO.writeUnfFile(filename, ng, &G_watersheds[0]);


	cout << "Storing information about waterbasins ..." << endl;
	float basin_area[nb2_max];
	float coordinates[nb2_max][2];

	for (int n = 0; n <= nb; ++n) {
		basin_area[n] = 0;
		basin_cells[n] = 0;
		coordinates[n][0] = -999;
	}
	for (int n = 0; n < ng; ++n) {
		basin_area[G_watersheds[n]] += area[G_row[n] - 1] * ((G_w_freq[n] + G_freq[n]) / 100.0);	// [km2] 
		basin_cells[G_watersheds[n]]++;
	}
	for (int n = 0; n < ng; ++n) {
//		if ((5 == G_LDD[n]) && (G_watersheds[n] >= 1))
		if (( (5 == G_LDD[n]) || (-1 == G_LDD[n]) ) && (G_watersheds[n] >= 1)) 
		{
			if (coordinates[G_watersheds[n]][0] > -998)
				printf("Something is wrong with the watershed map: %d %d %f\n",
					   n, G_watersheds[n], coordinates[G_watersheds[n]][0]);
			coordinates[G_watersheds[n]][0] = ((G_column[n] - 1) / 2.0) - 179.75;
			coordinates[G_watersheds[n]][1] = -((G_row[n] - 1) / 2.0) + 89.75;
		}
	}
	FILE *file_ptr;

	sprintf(filename, "%s/BASINS.OUT", output_dir);
	file_ptr = fopen(filename, "w");
	fprintf(file_ptr, "# %s", getTimeString());
	fprintf(file_ptr, "#\n");
	fprintf(file_ptr, "# 1. column: number of waterbasin\n");
	fprintf(file_ptr, "# 2. column: number of cells of waterbasin [km2]\n");
	fprintf(file_ptr, "# 3. column: area of waterbasin [km2]\n");
	fprintf(file_ptr, "# 4. column: longitude of outflow point\n");
	fprintf(file_ptr, "# 5. column: latitude of outflow point\n");
	fprintf(file_ptr, "#\n");
	for (int n = 1; n <= nb; ++n)
		fprintf(file_ptr,
				"%4d %4d %15.5f %7.2f %7.2f\n",
				n, basin_cells[n], basin_area[n], coordinates[n][0], coordinates[n][1]);
	fclose(file_ptr);


	// calculate some additional information which are required for the
	// routing in watergap
	// these routines are independed from those above

	int G_outflow_cell[ng];

	// cell to which flow from the cell goes
	outflow_cell(G_LDD, GCRC_2, G_row, G_column, G_outflow_cell, output_dir);

	// numbers of neighbouring cells of each cell
	int G_neighbourCell[ng][8];
	neighbouring_cells(GCRC_2, G_row, G_column, G_neighbourCell, output_dir);

	// calculate distances between cells and create file which is consistent 
	// with 'inflow cells'
	float cell_distance[9][360];

	calculate_distances(cell_distance, output_dir);

	// calculate slope of flow direction
	float G_altitude[ng];

	sprintf(filename, "%s/GALTMOD.UNF0", input_dir);
	gridIO.readUnfFile(filename, ng, G_altitude);

	// calculate river length
	float G_meandering_ratio[ng];
	
	sprintf(filename, "%s/G_MEANDERING_RATIO.UNF0", input_dir);
	gridIO.readUnfFile(filename, ng, G_meandering_ratio);

    calculate_river_slope(G_meandering_ratio, cell_distance, G_altitude, G_outflow_cell, G_LDD, G_row, output_dir); //LA (2014-02-12): meandering ratio added
    calculate_river_length(G_meandering_ratio, G_row, G_column, cell_distance, G_outflow_cell, output_dir);

	int G_routOrder[ng]; 	
	rout_order(G_routOrder, G_inflow_cells, G_outflow_cell, output_dir);
	
	if (resOpt == 1) reservoir_prepare(G_outflow_cell, input_dir, output_dir, G_row);
	
}


void replace_flow_to_ocean(int GCRC_2[722][362], short *G_row, short *G_column, char *G_LDD)
{
	int n;
	short row, col;

	cout << "Assigning pits to cells with flow into ocean cells ..." << endl;

	// give value of 5 also to those cells which have 100 percent  
	// land but which have a flow direction into a cell which is   
	// 100 percent ocean water (and therefore not part of the land
	// mask)
	// CAUTION !!! If this is used, further adjustments might have to be done 
	// for inland sinks !!!
	for (n = 0; n <= (ng - 1); n++)
		//if (G_LDD[n] != 5) 
		if ( (G_LDD[n] != 5) && (G_LDD[n] != -1) ) 
		{
			row = G_row[n];
			col = G_column[n];

			// flow from NE to SW 
			if ((GCRC_2[col - 1][row + 1] == 0)
				&& (G_LDD[n] == 1)) {
				G_LDD[n] = 5;
			}
			// flow from N to S 
			if ((GCRC_2[col][row + 1] == 0)
				&& (G_LDD[n] == 2)) {
				G_LDD[n] = 5;
			}
			// flow from NW to SE
			if ((GCRC_2[col + 1][row + 1] == 0)
				&& (G_LDD[n] == 3)) {
				G_LDD[n] = 5;
			}
			// flow from E to W  
			if ((GCRC_2[col - 1][row] == 0)
				&& (G_LDD[n] == 4)) {
				G_LDD[n] = 5;
			}
			// flow from W to E    
			if ((GCRC_2[col + 1][row] == 0)
				&& (G_LDD[n] == 6)) {
				G_LDD[n] = 5;
			}
			// flow from SE to NW    
			if ((GCRC_2[col - 1][row - 1] == 0)
				&& (G_LDD[n] == 7)) {
				G_LDD[n] = 5;
			}
			// flow from S to N    
			if ((GCRC_2[col][row - 1] == 0)
				&& (G_LDD[n] == 8)) {
				G_LDD[n] = 5;
			}
			// flow from SW to NE 
			if ((GCRC_2[col + 1][row - 1] == 0)
				&& (G_LDD[n] == 9)) {
				G_LDD[n] = 5;
			}
		}
}

void inflow_cells(char *G_LDD,
				  int GCRC_2[722][362],
				  short *G_row, short *G_column, int G_inflow_cells[ng][9], char *output_dir)
{
	int a1, a2, a3, a4, a6, a7, a8, a9;
	short row, col;

	char filename[250];

	cout << "Determine cell numbers of inflowing cells ..." << endl;

	// determine the cell numbers from which inflow
	// occurs into the cell of interest             
	a1 = 0;
	a2 = 0;
	a3 = 0;
	a4 = 0;
	a6 = 0;
	a7 = 0;
	a8 = 0;
	a9 = 0;
	for (int n = 0; n < ng; ++n)
		for (short i = 0; i < 9; ++i)
			G_inflow_cells[n][i] = 0;

	for (int n = 0; n < ng; ++n) {

		row = G_row[n];
		col = G_column[n];

		// flow from NE to SW 
		if ((GCRC_2[col + 1][row - 1] != 0)
			&& (G_LDD[GCRC_2[col + 1][row - 1] - 1] == 1)) {
			G_inflow_cells[n][9 - 1] = GCRC_2[col + 1][row - 1];
			a1++;
		}
		// flow from N to S 
		if ((GCRC_2[col][row - 1] != 0)
			&& (G_LDD[GCRC_2[col][row - 1] - 1] == 2)) {
			G_inflow_cells[n][8 - 1] = GCRC_2[col][row - 1];
			a2++;
		}
		// flow from NW to SE
		if ((GCRC_2[col - 1][row - 1] != 0)
			&& (G_LDD[GCRC_2[col - 1][row - 1] - 1] == 3)) {
			G_inflow_cells[n][7 - 1] = GCRC_2[col - 1][row - 1];
			a3++;
		}
		// flow from E to W  
		if ((GCRC_2[col + 1][row] != 0)
			&& (G_LDD[GCRC_2[col + 1][row] - 1] == 4)) {
			G_inflow_cells[n][6 - 1] = GCRC_2[col + 1][row];
			a4++;
		}
		// flow from W to E    
		if ((GCRC_2[col - 1][row] != 0)
			&& (G_LDD[GCRC_2[col - 1][row] - 1] == 6)) {
			G_inflow_cells[n][4 - 1] = GCRC_2[col - 1][row];
			a6++;
		}
		// flow from SE to NW    
		if ((GCRC_2[col + 1][row + 1] != 0)
			&& (G_LDD[GCRC_2[col + 1][row + 1] - 1] == 7)) {
			G_inflow_cells[n][3 - 1] = GCRC_2[col + 1][row + 1];
			a7++;
		}
		// flow from S to N    
		if ((GCRC_2[col][row + 1] != 0)
			&& (G_LDD[GCRC_2[col][row + 1] - 1] == 8)) {
			G_inflow_cells[n][2 - 1] = GCRC_2[col][row + 1];
			a8++;
		}
		// flow from SW to NE 
		if ((GCRC_2[col - 1][row + 1] != 0)
			&& (G_LDD[GCRC_2[col - 1][row + 1] - 1] == 9)) {
			G_inflow_cells[n][1 - 1] = GCRC_2[col - 1][row + 1];
			a9++;
		}
	}
	cout << "  " << a1 + a2 + a3 + a4 + a6 + a7 + a8 + a9
		<< " cells, which have flow to other cells." << endl;
	cout << "  SW: " << a1 << endl;
	cout << "  S : " << a2 << endl;
	cout << "  SE: " << a3 << endl;
	cout << "  W : " << a4 << endl;
	cout << "  E : " << a6 << endl;
	cout << "  NW: " << a7 << endl;
	cout << "  N : " << a8 << endl;
	cout << "  NE: " << a9 << endl;

#ifdef debug
	for (n = 0; n <= ng - 1; n++) {
		printf("%ld : ", n + 1);
		for (short j = 0; j <= 8; j++)
			printf("%ld ", G_inflow_cells[n][j]);
		printf("\n");
	}
#endif

	sprintf(filename, "%s/G_INFLC.9.UNF4", output_dir);
	gridIO.writeUnfFile(filename, 9 * ng, &G_inflow_cells[0][0]);
}


void flow_accumulation(int G_inflow_cells[ng][9], short *G_flow_acc)
{
	int n, i, cells_left, cells_left_prev_step;
	short to_be_done_later, flow_acc;

	cout << "Calculating flow accumulation ..." << endl;

	// find cells which have no inflow   
	// and set accumulation value to 1 (the cell itself)
	for (n = 0; n <= (ng - 1); n++)
		G_flow_acc[n] = 1;
	// 1: cell has no inflow from upstream
	// 0: inflow cell has been found
	for (n = 0; n <= (ng - 1); n++)
		for (i = 0; i <= 8; i++)
			if (G_inflow_cells[n][i] != 0)
				G_flow_acc[n] = 0;

	cells_left = 99999;

	// loop over all gridcells: 
	// find those cells for which flow accumulation 
	// is already calculated for all inflow cells   
	do {
		cells_left_prev_step = cells_left;
		cells_left = 0;
		for (n = 0; n <= (ng - 1); n++)
			if (G_flow_acc[n] == 0)	// cell has not yet an accumulation value 
			{
				to_be_done_later = 0;
				for (i = 0; i <= 8; i++)
					if ((G_inflow_cells[n][i] != 0)
						&& (G_flow_acc[G_inflow_cells[n][i] - 1] == 0))
						// has flow acc. been calculated for all inflow cells 'i' ? 
						to_be_done_later = 1;
				if (0 == to_be_done_later) {
					flow_acc = 1;	// the cell itself 
					// sum up accumulation values of all inflow cells 
					for (i = 0; i <= 8; i++)
						if (G_inflow_cells[n][i] != 0)
							flow_acc += G_flow_acc[G_inflow_cells[n][i] - 1];
					G_flow_acc[n] = flow_acc;
				} else
					cells_left++;
			}
		cout << "  Flow accumulation: still to do: " << cells_left << " cells" << endl;
	}
	while ((cells_left > 0) && (cells_left_prev_step - cells_left != 0));
}


unsigned short derive_basins(char *G_LDD,
							 int G_inflow_cells[ng][9],
							 unsigned short *G_watersheds, char *output_dir)
{
	char filename[250];
	char G_done[ng];	// array to mark cells where calculations are finished 

	// distance to outlet cell
	// measured as number of cells
	unsigned short G_dist[ng];

	unsigned short nb;	// number of waterbasins 

	cout << "Deriving basins from flow direction map ..." << endl;

	// assign numbers to watersheds                       
	// cells which are not part of the landmask: 0  
	// all the others are numbered in order of appearance 

	int i = 0, j = 0, n;

	for (n = 0; n < ng; n++) {
		G_done[n] = -99;
		G_watersheds[n] = 0;
	}
	for (n = 0; n < ng; n++) {
		if ((0 == G_LDD[n]) || (-99 == G_LDD[n])) {
			G_done[n] = 2;
			G_watersheds[n] = 0;
			j++;

			// this is just a temporary solution                
			// this case should not occur in the final LDD file 
			G_LDD[n] = 5;
		}
		//if (5 == G_LDD[n]) 
		if ( (5 == G_LDD[n]) || (-1 == G_LDD[n]) ) 
		{
			G_done[n] = 1;
			i++;
			G_watersheds[n] = i;
		}
	}
	nb = i;
	cout << "  " << j << " cells are not part of the LDD-map." << endl;
	cout << "  " << i << " endpoints have been found." << endl;
	if (i > (nb_max - 1)) {
		cerr << "Size of array is to small in 'derive_basins'.\n";
		exit(-1);
	}
	int k = 1;

	i = 1;
	do {
		k = 0;
		for (n = 0; n <= (ng - 1); n++) {
			if (1 == G_done[n]) {
				G_dist[n] = i - 1;
				k++;
				G_done[n] = 2;
				// done == 1:  
				//     cell is already part of a watershed but flow directions 
				//     have not yet been followed.                
				// done == 2: 
				//     everything is finished for that cell                    
				for (j = 0; j <= 8; j++)
					if (G_inflow_cells[n][j] != 0) {
						G_watersheds[G_inflow_cells[n][j] - 1] = G_watersheds[n];
						G_done[G_inflow_cells[n][j] - 1] = -1;
						// a faster solution would be to set G_done = 1 here
						// but for that case G_dist is not calculated well
					}
			}
		}
		for (n = 0; n <= (ng - 1); n++)
			if (-1 == G_done[n])
				G_done[n] = 1;
		i++;
		//cout << i << ' ' << k << endl;
	}
	while (k != 0);
	i--;
	cout << "  Number of cells within the intest chain: " << i << endl;

	sprintf(filename, "%s/G_CELLS_TO_OUTLET.UNF2", output_dir);
	gridIO.writeUnfFile(filename, ng, G_dist);

	sprintf(filename, "%s/G_BASINS.UNF2", output_dir);
	gridIO.writeUnfFile(filename, ng, G_watersheds);

	return nb;	// number of waterbasins
}

unsigned short reindex_waterbasins(unsigned short nb,
								   short *basin_cells, unsigned short *G_watersheds)
{
	cout << "Reindexing waterbasins ..." << endl;

	// give new number to watersheds:  
	// only those are considered which 
	// have more than one cell.       
	// for the other no routing has to 
	// be performed.                   
	int j = 0;

	for (int i = 1; i <= nb; i++) {
		if (1 == basin_cells[i])
			for (int n = 0; n <= (ng - 1); n++) {
				if (i == G_watersheds[n]) {
					G_watersheds[n] = 0;
					break;
				}
		} else {
			j++;
			for (int n = 0; n <= (ng - 1); n++)
				if (i == G_watersheds[n])
					G_watersheds[n] = j;
		}
	}
	cout << "  " << j << " watersheds which have more than one cell." << endl;
	if (j > (nb2_max - 1)) {
		cerr << "Size of array is to small in 'reindex_waterbasins'\n";
		exit(-1);
	}

	return j;	// new number of waterbasins
}


void outflow_cell(char *G_LDD,
				  int GCRC_2[722][362],
				  short *G_row, short *G_column, int *G_outflow_cell, char *output_dir)
{
	short row, col;
	int n;

	char filename[250];

	cout << "Determine cell number of the cell to which flow goes ..." << endl;

	// determine the cell number to which flow
	// occurs 
	for (n = 0; n <= ng - 1; n++) {

		row = G_row[n];
		col = G_column[n];

		//if (5 == G_LDD[n])
		if ( (5 == G_LDD[n]) || (-1 == G_LDD[n]) )
			G_outflow_cell[n] = 0;

		// flow from NE to SW 
		if (1 == G_LDD[n])
			G_outflow_cell[n] = GCRC_2[col - 1][row + 1];

		// flow from N to S 
		if (2 == G_LDD[n])
			G_outflow_cell[n] = GCRC_2[col][row + 1];

		// flow from NW to SE
		if (3 == G_LDD[n])
			G_outflow_cell[n] = GCRC_2[col + 1][row + 1];

		// flow from E to W  
		if (4 == G_LDD[n])
			G_outflow_cell[n] = GCRC_2[col - 1][row];

		// flow from W to E  
		if (6 == G_LDD[n])
			G_outflow_cell[n] = GCRC_2[col + 1][row];

		// flow from SE to NW    
		if (7 == G_LDD[n])
			G_outflow_cell[n] = GCRC_2[col - 1][row - 1];

		// flow from S to N  
		if (8 == G_LDD[n])
			G_outflow_cell[n] = GCRC_2[col][row - 1];

		// flow from SW to NE 
		if (9 == G_LDD[n])
			G_outflow_cell[n] = GCRC_2[col + 1][row - 1];
	}
	sprintf(filename, "%s/G_OUTFLC.UNF4", output_dir);
	gridIO.writeUnfFile(filename, ng, G_outflow_cell);
}

void neighbouring_cells(int GCRC_2[722][362], short *G_row, short *G_column, int G_neighbourCell[ng][8], char *output_dir)
{
	short row, col;
	long n;

	char filename[250];
	//Ellen
	cout << "find neighbouring cells of each cell..." << endl;

	// determine the cell number to which flow
	// occurs
	for (n = 0; n <= ng - 1; n++) {
		row = G_row[n];
		col = G_column[n];

       	G_neighbourCell[n][0] = GCRC_2[col + 1][row];
       	G_neighbourCell[n][1] = GCRC_2[col + 1][row - 1];
       	G_neighbourCell[n][2] = GCRC_2[col][row - 1];
       	G_neighbourCell[n][3] = GCRC_2[col - 1][row - 1];
       	G_neighbourCell[n][4] = GCRC_2[col - 1][row];
       	G_neighbourCell[n][5] = GCRC_2[col - 1][row + 1];
       	G_neighbourCell[n][6] = GCRC_2[col][row + 1];
       	G_neighbourCell[n][7] = GCRC_2[col + 1][row + 1];

	}
	sprintf(filename, "%s/G_NEIGHBOUR_CELLS.8.UNF4", output_dir);
	gridIO.writeUnfFile(filename, ng * 8, &G_neighbourCell[0][0]);
}

void calculate_distances(float cell_distance[9][360], char *output_dir)
  // calculate vertical and horizontal distance to neighbour cell   
{
	char filename[250];

	//float cell_distance[9][360];
	// organized in the following structure: 
	// 6 7 8 
	// 3 4 5 
	// 0 1 2 
	const float pi = 3.141592653589793;
	const float earth_radius = 6371.211;	// [km] 

	// radius of a sphere with same volume 
	float vert_dist;
	float horiz_dist;
	int i;
	float l;

	cout << "Calculating distances between center of cells ..." << endl;

	vert_dist = 2 * pi * earth_radius * 0.5 / 360.0;
	for (i = 0; i <= 359; i++) {
		l = 0.25 + i * 0.5;
		horiz_dist = 2 * pi * earth_radius * sin(l * pi / 180.0) * 0.5 / 360.0;
		cell_distance[1][i] = vert_dist;
		cell_distance[3][i] = horiz_dist;
		cell_distance[4][i] = 0;
		cell_distance[5][i] = horiz_dist;
		cell_distance[7][i] = vert_dist;
	}

	//                             
	// calculate diagonal distance 
	//                             

	for (i = 0; i <= 358; i++)
		cell_distance[0][i] =
			sqrt(cell_distance[3][i] * cell_distance[3][i + 1] +
				 cell_distance[1][i] * cell_distance[1][i]);
	cell_distance[0][359] = -99;
	for (i = 1; i <= 359; i++)
		cell_distance[6][i] =
			sqrt(cell_distance[3][i] * cell_distance[3][i - 1] +
				 cell_distance[1][i] * cell_distance[1][i]);
	cell_distance[6][0] = -99;
	for (i = 0; i <= 359; i++) {
		cell_distance[2][i] = cell_distance[0][i];
		cell_distance[8][i] = cell_distance[6][i];
	}

#ifdef debug
	for (i = 0; i <= 359; i++) {
		printf("%3d ", i);
		for (short j = 0; j <= 8; j++)
			printf("%10.6f ", cell_distance[j][i]);
		printf("\n");
	}
#endif

	// store distances to file 
	sprintf(filename, "%s/GCELLDIST.9.UNF0", output_dir);
	gridIO.writeUnfFile(filename, 9 * 360, &cell_distance[0][0]);
}

void calculate_river_slope(float *G_meandering_ratio,
                           float cell_distance[9][360],
                           float *G_altitude, //G_altitude = GALTMOD.UNF0
						   int *G_outflow_cell, char *G_LDD, short *G_row, char *output_dir)
{
	cout << "Calculating river slope ...\n";
	float G_slope[ng];
	int n, n_negative = 0, n_zero = 0;

	for (n = 0; n <= ng - 1; n++) {
		//if ((G_outflow_cell[n] != 0) && (G_LDD[n] != 5)) 
		if ( (G_outflow_cell[n] != 0) && ((G_LDD[n] != 5) && (G_LDD[n] != -1)) ) 
		{
            G_slope[n] = (G_altitude[n] - G_altitude[G_outflow_cell[n] - 1])
                / (1000.0 * cell_distance[G_LDD[n] - 1][G_row[n] - 1]*((G_meandering_ratio[n]+G_meandering_ratio[G_outflow_cell[n]-1])/2)); //(2014-02-12) LA added: *((G_meandering_ratio[n]+G_meandering_ratio[G_outflow_cell[n]-1])/2))
			if (G_slope[n] < 0)
				n_negative++;
			if (G_slope[n] == 0)
				n_zero++;
		} else {
			G_slope[n] = 0;
			n_zero++;
		}
	}

	float minslope = 0.00001;

	cout << "  " << n_zero << " cells have river slope of 0 m/m.\n";
	cout << "  " << n_negative << " cells have negative river slope.\n";
	cout << "  Setting slope to a minimum of " << minslope << " m/m.\n";
	for (n = 0; n < ng; n++)
		if (G_slope[n] < minslope)
			G_slope[n] = minslope;

	char filename[250];

	sprintf(filename, "%s/G_RIVERSLOPE.UNF0", output_dir);
	gridIO.writeUnfFile(filename, ng, G_slope);
}

void calculate_river_length(float *G_meandering_ratio, short *G_row, short *G_column, float cell_distance[9][360], 
  int *G_outflow_cell, char *output_dir) {
	
	cout << "Calculating river length..."<<endl;
	float G_river_length[ng];
	
	int n_outflow;

	for (int n = 0; n < ng; n++) {
		n_outflow=G_outflow_cell[n]-1;
		if (n_outflow==-1) 
			G_river_length[n]=55.;
		else if (G_row[n]==G_row[n_outflow])
			G_river_length[n]=cell_distance[3][G_row[n]-1];  //west
		else {
			if (G_row[n]>G_row[n_outflow]) {
			 	if (G_column[n]==G_column[n_outflow])
			 		G_river_length[n]=cell_distance[7][G_row[n]-1];  //north
			 	else G_river_length[n]=cell_distance[8][G_row[n]-1];                       //north east
			}
			else {
				if (G_column[n]==G_column[n_outflow])
					G_river_length[n]=cell_distance[1][G_row[n]-1];  //south
			 	else G_river_length[n]=cell_distance[2][G_row[n]-1];                       //south east
			}
		}
		
		if (n_outflow == -1) //"-1" are outlet cells or internal sinks 
			G_river_length[n] *= G_meandering_ratio[n];
		else
			G_river_length[n] *= ((G_meandering_ratio[n]+G_meandering_ratio[n_outflow])/2);

	}
	
	char filename[250];

	sprintf(filename, "%s/G_RIVER_LENGTH.UNF0", output_dir);
	gridIO.writeUnfFile(filename, ng, G_river_length);

}

void rout_order(int *G_routOrder, int G_upstreamCells[ng][9], int *G_outflow_cell, char *output_dir) {
	
	int routing[ng];  	// count inflow cells
	int list[ng]; 	 	// count inflow cells
	
	cout << "Calculating routing order" << endl;

	for (int n = 0; n < ng; n++){ 
		G_routOrder[n]=-99;
		routing[n]=0;
	}
	int routOrder=1; // next free routing number in basin
	
	for (int n = 0; n < ng; ++n)
		for (short i = 0; i < 9; ++i) 
			if (G_upstreamCells[n][i]>0) routing[n] ++;
	
	int counter=0;
	int counterStep;
	
	do {
		counterStep=0;
		for (int n = 0; n < ng; ++n) {
			if (routing[n]==0) {G_routOrder[n]=routOrder++; counter++; counterStep++; list[n]=1;}
			else list[n]=0;
		}
		
		for (int n = 0; n < ng; ++n) {
			if (list[n]) {
				routing[n]--;
				if (G_outflow_cell[n]>0) routing[G_outflow_cell[n]-1]--;
			} 
			
		}
		
	} while (counterStep>0);

	if (counter == ng)
		cout << "Routing order is  finished. " << counter << " cells are in right routing order\n";

	for (int n = 0; n < ng; ++n)
		if (routing[n]>=0 || G_routOrder[n]==-99) {
			cerr << "cell " << n+1 << '\t' << routing[n] << '\t'<< G_routOrder[n] 
			  << ": Routing order is not finished. Something is wrong in G_upstreamCells or/and G_outflow_cell\n";
			exit(-1);
	}
	
	char filename[250];

	sprintf(filename, "%s/G_ROUT_ORDER.UNF4", output_dir);
	gridIO.writeUnfFile(filename, ng, G_routOrder);

} //rout_order()

void reservoir_prepare(int *G_outflow_cell, char *input_dir, char *output_dir, short *G_row){

	char filename[250];

	int (*G_reservoir_area);	// km2
	G_reservoir_area = new int [ng];
	if (!G_reservoir_area) cerr << "ERROR: not enough memory for G_reservoir_area\n" << endl;

        sprintf(filename, "%s/G_RESAREA.UNF0", input_dir);
	gridIO.readUnfFile(filename, ng, &G_reservoir_area[0]);

	//-----------------------------------------------

    double (*G_mean_outflow);   //mean annual inflow to global lakes plus P-PET of global lake (1961-1990)
    G_mean_outflow  	= new double [ng];
    if (!G_mean_outflow) cerr << "ERROR: not enough memory for G_mean_outflow\n" << endl;

    //read mean annual outflow of reservoirs (was called G_MEAN_INFLOW but is mean inflow plus P-PET of reservoir in fact (long term mean)
    sprintf(filename, "%s/G_MEAN_OUTFLOW.UNF0", input_dir);
    gridIO.readUnfFile(filename, ng, &G_mean_outflow[0]);

	//-----------------------------------------------
	// store inflow of relevant reservoirs in each cell and the total sum 
	double G_inflowReservoirs[ng];
	double G_alloc_coeff[ng][reservoir_dsc];//allocation coefficient (< 1 if more than one reservoir is upstream)

	// initialize
	for (int n = 0; n < ng; n++) {
		G_inflowReservoirs[n] = 0.;

		for (int i = 0; i < reservoir_dsc; i++)
			G_alloc_coeff[n][i] = 1.0;
	}
	
	// store mean inflow for each cell and each reservoir
	for (int n = 0; n < ng; n++){
		if (G_reservoir_area[n] > 0){ //cell is a reservoir
			short i=0; int downstreamCell=G_outflow_cell[n];

			while (i<reservoir_dsc && downstreamCell>0 && G_reservoir_area[downstreamCell-1]<=0) {
                G_inflowReservoirs[downstreamCell-1]+=G_mean_outflow[n];

				// next downstream cell
				downstreamCell = G_outflow_cell[downstreamCell-1]; i++;
			}

		}
	} //for(n)


	// calculate allocation coefficient
	for (int n = 0; n < ng; n++){
		if (G_reservoir_area[n] > 0) {
			short i=0; int downstreamCell=G_outflow_cell[n];
			
			while (i<reservoir_dsc && downstreamCell>0 && G_reservoir_area[downstreamCell-1]<=0) {
                G_alloc_coeff[n][i]=G_mean_outflow[n]/G_inflowReservoirs[downstreamCell-1];
				
				// next downstream cell
				downstreamCell = G_outflow_cell[downstreamCell-1]; i++;
				
			}
			
		}
		
	} // for(n)
	
	sprintf(filename, "%s/G_ALLOC_COEFF.%d.UNF0", output_dir, reservoir_dsc);
	gridIO.writeUnfFile(filename, ng*reservoir_dsc, &G_alloc_coeff[0][0]);
	
	cout << "calculation of allocation coefficient for downstream cells finished\n";
	
	//-----------------------------------------------

    double (*G_mean_outflow_monthly)[12];   //mean annual outflow to global lakes (1971-1900)
    G_mean_outflow_monthly  	= new double [ng][12];
    if (!G_mean_outflow_monthly) cerr << "ERROR: not enough memory for G_mean_outflow_monthly\n" << endl;

    //read mean annual inflow and P-PET of reservoirs (long term mean)
    sprintf(filename, "%s/G_MEAN_OUTFLOW.12.UNF0", input_dir);
    gridIO.readUnfFile(filename, ng*12, &G_mean_outflow_monthly[0][0]);
    cout << "G_MEAN_OUTFLOW read\n";
	
	char G_start_month[ng];	// start of operational year
	
	for (int n = 0; n < ng; n++){
		G_start_month[n]=1;
		
		//------------------
		/*
        // variante 1: first month (from january) when mean monthly inflow plus P - PET < mean annual inflow plus P - PET
		for (int month=0; month<12; month++) {
            if (n==0) cout <<"----\t"<< n<<'\t'<<month<< '\t'<<G_mean_outflow_monthly[n][month]<<'\t'<<G_mean_outflow[n]<<endl;
            if (G_mean_outflow_monthly[n][month]<G_mean_outflow[n]) {
				G_start_month[n]=month+1; 
				if (n==0) cout <<"++++G_start_month["<<n<<"]:"<<G_start_month[n]<<endl;
				break; 
			}
		}
		*/

		//------------------
		/*
        // variante 2: in northern hemisphere: from january to december first month when mean monthly inflow plus P - PET < mean annual inflow plus P - PET
        //           : in southern hemisphere: from december backwards first month when mean monthly inflow plus P - PET > mean annual inflow plus P - PET
		if (-(G_row[n] / 2.0 - 90.25) > 0.0) {// northern hemisphere
			for (int month=0; month<12; month++) {
                if (G_mean_outflow_monthly[n][month]<G_mean_outflow[n]) {G_start_month[n]=month+1; break;}
			} // for(month)
		}
		else { // southern hemisphere
			for (int month=11; month>=0; month--) {
                if (G_mean_outflow_monthly[n][month]>=G_mean_outflow[n]) {G_start_month[n]=month+1; break;}
			} // for(month)
		}
		*/
		//------------------
		// variante 3: begin of longest dry period
		//           : 
		int start, dry_length=0, start_=0, counter_length=0; bool last_dry=true;
		
		// finden erster Monat wo Monatswert >= Jahreswert. (beg_month)
		// dann finden laengste Trockenperiode (wo Monatswert<Jahreswert). Suchen nicht ab Januar sondern ab beg_month
                // damit erreichen wir, dass alle trockene Periode ganz sind (auch wenn die ueber Jahresende dauern)
		int beg_month=0; int month=0;
        while (month<12 && G_mean_outflow_monthly[n][month]<G_mean_outflow[n]) month++;
		beg_month=month;
		
		for (int m=0; m<12; m++) {
			month=m+beg_month; if (month>11) month-=12;
            if (G_mean_outflow_monthly[n][month]<G_mean_outflow[n] ) {
				counter_length++; 
				if (!last_dry) {start_=month; last_dry=true;}
			}
			else if (last_dry) { 
				last_dry=false; 
				if(counter_length>dry_length) {
					start=start_; 
					dry_length=counter_length;
					counter_length=0;
				}
			}
		}
		if (counter_length>dry_length) {start=start_; dry_length=counter_length;}
		
		G_start_month[n]=start+1;
		//------------------
		
	} // for(n)
	{
		int count=0;
		for (int n = 0; n < ng; n++) if (G_start_month[0]==0) count++;
		cout << "*** is 0: " <<count<<endl;
	}
	sprintf(filename, "%s/G_START_MONTH.UNF1", output_dir);
	gridIO.writeUnfFile(filename, ng, &G_start_month[0]);
	
	
	cout << "calculation of first month of release period finished\n";

    delete[] G_mean_outflow;         G_mean_outflow = NULL;
    delete[] G_mean_outflow_monthly; G_mean_outflow_monthly = NULL;
	delete[] G_reservoir_area;      G_reservoir_area = NULL;
	
	
} //reservoir_prepare()
