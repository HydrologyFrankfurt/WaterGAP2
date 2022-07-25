
/***********************************************************************
*
* see former changes at file calib_basins.cpp.versioninfos.txt
***********************************************************************/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdio>
#include <cstring>
#include <cstdlib>

#include "gridio.h"
#include "geo.h"
#include "timestring.h"
#include "def.h"
#include "stack.h"
#include "calib_basins.h"
#include "option.h"

using namespace std;

extern gridioClass gridIO;
extern geoClass geo;


short upstreamBasin(const int CellNumber,
					const signed short basinNumber,
					const int G_inflow_cells[ng][9],
					signed short *G_calib_basin, char *overlapBasinNumString);

short cbasinClass::prepare(char *input_dir, char *output_dir, char *routing_dir)
{
	extern optionClass options;
	char filename[250];

	// if streams are used instead of 'sprintf'
	// this should look like this:
	//
	//#include <strstream.h>
	// ...
	//ostrstream fname;
	//fname << input_dir << "/GCRC.UNF4" << ends;;
	//filename = fname.str();

	int G_inflow_cells[ng][9];

	sprintf(filename, "%s/G_INFLC.9.UNF4", routing_dir);
	gridIO.readUnfFile(filename, 9 * ng, &G_inflow_cells[0][0]);

	signed short G_calib_basin[ng];
	int n;

	for (n = 0; n <= (ng - 1); n++)
		G_calib_basin[n] = 0;


	// open file which contains location of stations
	ifstream station_file(options.filename_stations);

	if (!station_file) {
		cerr << "Can not open file " << options.filename_stations << endl;
		exit(-1);
	}
	// open file for output
	sprintf(filename, "%s/STATION_LIST.OUT", output_dir);
	ofstream output_file(filename);

	if (!output_file) {
		cerr << "Can not open file " << filename << " for writing" << endl;
		exit(-1);
	}
	output_file << "# " << getTimeString();

	signed short G_short[ng];
	short overlap, cellAlreadyInList;
	char overlapBasinNumString[250];
	float longitude, latitude;
	char string[250];
	const short maxLength = 15;
	short length;
	int k;

	extStack < int >stackCellNum;
	int cn, cn2;
	Stack < float >stackUpArea;
	float upArea;
	Stack < float >stackGamma;
	float gamma2;
	Stack < float >stackCellCorrFact;
	float cellCorrFact2;
	Stack < float >stackStatCorrFact;
	float statCorrFact2;
	Stack < char *>stackNamePointer;
	char *basinName;

	short basin = 0;

	do {
		station_file >> string;
		if (station_file) {
			station_file >> longitude >> latitude;
			station_file >> gamma2;
			station_file >> cellCorrFact2;
			station_file >> statCorrFact2;

			// allocate memory and copy string with the desired length into it
			length = strlen(string);
			if (length >= maxLength) {
				basinName = new char[maxLength + 1];

				length = maxLength;
				strncpy(basinName, string, length);
			} else {
				basinName = new char[length + 1];

				strncpy(basinName, string, length);
				basinName[length] = '\0';
			}

			output_file << "Station name:       " << basinName << endl;
			output_file << "Longitude/Latitude: " << longitude << " / " << latitude << endl;

			// determine cell number of station
			cn = geo.cellNumByLonLat(longitude, latitude);
			output_file << "Cell number:        " << cn << endl;
			if (cn <= 1)
				output_file << "Not part of the landmask! Skipping! " << endl;
			else {
				// is that cell number already in the list ?
				cellAlreadyInList = 0;
				for (k = 1; k <= stackCellNum.getNumberOfElements(); k++) {
					stackCellNum.getElement(k, cn2);
					if (cn == cn2) {
						cellAlreadyInList = 1;
						break;
					}
				}

				if (cellAlreadyInList) {
					output_file << "Cell is already lowest cell of basin: "
						<< stackCellNum.getNumberOfElements() - (k - 1)
						<< endl << "Skipping!" << endl;
				} else {
					// everything OK, cell should be used

					stackCellNum.push(cn);
					stackGamma.push(gamma2);
					stackCellCorrFact.push(cellCorrFact2);
					stackStatCorrFact.push(statCorrFact2);
					stackNamePointer.push(basinName);
					basin++;
					output_file << "Basin number: " << basin << endl;

					overlap = upstreamBasin(cn,
											basin,
											G_inflow_cells, G_calib_basin, overlapBasinNumString);
					switch (overlap) {
					case 0:
						break;
					case 1:
						output_file
							<< "Basin is upstream to another one defined before: \n"
							<< overlapBasinNumString << endl;
						break;
					case 2:
						output_file
							<< "Basin is downstream to the following "
							<< "basins defined before:\n" << overlapBasinNumString << endl;
						break;
					case 3:
						output_file
							<< "Basin is downstream and upstream to others defined before: \n"
							<< overlapBasinNumString << endl;
						break;
					}

					// calculate area for the different cases
					if ((0 == overlap) || (1 == overlap)) {
						upArea = geo.nonOceanAreaOfObjectInGrid(basin, G_calib_basin, ng);
						stackUpArea.push(upArea);
					}
					if ((2 == overlap) || (3 == overlap)) {
						// create a temporary grid with this basin alone.
						// this map is used to calculate the upstream area 
						// of the basin.
						for (n = 0; n <= (ng - 1); n++)
							G_short[n] = 0;
						overlap = upstreamBasin(cn,
												basin,
												G_inflow_cells, G_short, overlapBasinNumString);
						upArea = geo.nonOceanAreaOfObjectInGrid(basin, G_calib_basin, ng);
						stackUpArea.push(upArea);
					}
				}
			}
			output_file << "------------------------------------------------\n";
		}
	} while (station_file);
	station_file.close();

	numberOfCalibBasins = basin;

	// everything OK?
	if ((numberOfCalibBasins != stackCellNum.getNumberOfElements()) ||
		(numberOfCalibBasins != stackUpArea.getNumberOfElements()) ||
		(numberOfCalibBasins != stackGamma.getNumberOfElements()) ||
		(numberOfCalibBasins != stackNamePointer.getNumberOfElements()) ||
		(numberOfCalibBasins != stackCellCorrFact.getNumberOfElements()) ||
		(numberOfCalibBasins != stackStatCorrFact.getNumberOfElements())) {
		cerr << "Internal problem in cbasinClass\n";
		exit(-1);
	}
	// allocate memory for arrays and
	// copy contents of stacks into them
	int i;
	cellNum = new int[stackCellNum.getNumberOfElements()];

	for (i = stackCellNum.getNumberOfElements(); i >= 1; i--)
		if (stackCellNum.pop(cn))
			cellNum[i - 1] = cn;

	gamma = new float[stackGamma.getNumberOfElements()];

	for (i = stackGamma.getNumberOfElements(); i >= 1; i--)
		if (stackGamma.pop(gamma2))
			gamma[i - 1] = gamma2;

	cellCorrFactor = new float[stackCellCorrFact.getNumberOfElements()];

	for (i = stackCellCorrFact.getNumberOfElements(); i >= 1; i--)
		if (stackCellCorrFact.pop(cellCorrFact2))
			cellCorrFactor[i - 1] = cellCorrFact2;

	statCorrFactor = new float[stackStatCorrFact.getNumberOfElements()];

	for (i = stackStatCorrFact.getNumberOfElements(); i >= 1; i--)
		if (stackStatCorrFact.pop(statCorrFact2))
			statCorrFactor[i - 1] = statCorrFact2;

	float *upstreamArea;
	upstreamArea = new float[stackUpArea.getNumberOfElements()];

	for (i = stackUpArea.getNumberOfElements(); i >= 1; i--)
		if (stackUpArea.pop(upArea))
			upstreamArea[i - 1] = upArea;

	name = new char *[stackNamePointer.getNumberOfElements()];

	for (i = stackNamePointer.getNumberOfElements(); i >= 1; i--)
		if (stackNamePointer.pop(basinName))
			name[i - 1] = basinName;


	// calculate area of calibration areas
	calibArea = new float[numberOfCalibBasins];
	calibLandArea = new float[numberOfCalibBasins];

	for (n = 1; n <= numberOfCalibBasins; n++) {
		calibArea[n - 1]
			= geo.nonOceanAreaOfObjectInGrid((short) n, G_calib_basin, ng);
		calibLandArea[n - 1]
			= geo.landAreaOfObjectInGrid((short) n, G_calib_basin, ng);
	}

	// write the results into a file
	output_file << "# Area of basins:\n";
	output_file << "# 1. column: Basin number\n";
	output_file << "# 2. column: Total upstream area [km2]\n";
	output_file << "# 3. column: Upstream area up to next station(s).\n";
	output_file << "# 4. column: Upstream land area up to next station "
		<< "(IMAGE2.2 landmask).\n";
	for (n = 1; n <= numberOfCalibBasins; n++)
		output_file << setw(3) << n << " \t" << setprecision(2)
			<< setiosflags(ios::fixed | ios::showpoint)
			<< setw(12) << upstreamArea[n - 1] << " \t"
			<< setw(12) << calibArea[n - 1] << " \t" << setw(12) << calibLandArea[n - 1] << endl;
	output_file.close();

	// mark cell which is downstream to lowest cell of each basin.
	// this is required because transport into this cell should
	// also be calculated.
	// all of them get a value of -1 if they do not already have a value 
	// different from zero.
	int G_outflow_cell[ng];
	int outflow_cell;

	sprintf(filename, "%s/G_OUTFLC.UNF4", routing_dir);
	gridIO.readUnfFile(filename, ng, G_outflow_cell);
	for (n = 0; n <= numberOfCalibBasins - 1; n++) {
		outflow_cell = G_outflow_cell[cellNum[n] - 1];
		if (outflow_cell != 0)
			if (0 == G_calib_basin[outflow_cell - 1])
				G_calib_basin[outflow_cell - 1] = -1;
	}

	sprintf(filename, "%s/G_CALIB_BASIN.UNF2", output_dir);
	gridIO.writeUnfFile(filename, ng, G_calib_basin);

	delete upstreamArea;

	return numberOfCalibBasins;
}


short upstreamBasin(const int CellNumber,
					const signed short basinNumber,
					const int G_inflow_cells[ng][9],
					signed short *G_calib_basin, char *overlapBasinNumString)
{
	short not_yet_ready, base_value = 0;
	short upstream = 0, downstream = 0;
	char string[250];

	if (0 != G_calib_basin[CellNumber - 1]) {
#ifdef debug
		cout << "Previously defined basin " << ((int) G_calib_basin[CellNumber - 1])
			<< " is downstream to basin " << ((int) basinNumber) << endl;
#endif
		sprintf(overlapBasinNumString, "%d ", G_calib_basin[CellNumber - 1]);
		upstream = 1;	// upstream basin exits
		base_value = G_calib_basin[CellNumber - 1];
	}

	G_calib_basin[CellNumber - 1] = -1;
	do {
		not_yet_ready = 0;
		// if this variable is not set to 1 
		// within the following loop, then everything is ready

		for (int n = 0; n < ng; ++n) {
			if (-1 == G_calib_basin[n]) {
				for (short j = 0; j <= 8; ++j) {
					G_calib_basin[n] = -2;
					// cell does not need to be considered any more
					if (G_inflow_cells[n][j] != 0) {
						if (base_value == G_calib_basin[G_inflow_cells[n][j] - 1]) {
							// cell is not already part of another calibration basin
							// or the calibration basin is already recognized
							// because it is 'downstream' (see above)
							not_yet_ready = 1;
							G_calib_basin[G_inflow_cells[n][j] - 1] = -1;
						} else {
							// cell already has a value 
							//if (1 == upstream) 
							{
								// message should only be printed 
								// if this has not already been done
								// for the 'downstream' case (see above)    
#ifdef debug
								cout << "Previously defined basin "
									<< (int) G_calib_basin[G_inflow_cells[n][j] - 1]
									<< " is upstream to basin " << (int) basinNumber << endl;
#endif
								sprintf(string, "%d ", G_calib_basin[G_inflow_cells[n][j] - 1]);
								strcat(overlapBasinNumString, string);
								downstream = 2;	// downstream basin exists
							}
						}
					}
				}
			}
		}
	} while (not_yet_ready != 0);

	for (int n = 0; n <= (ng - 1); n++) {
		if (-2 == G_calib_basin[n])
			G_calib_basin[n] = basinNumber;
	}

	short overlapping=-1;

	if ((!upstream) && (!downstream))
		overlapping = 0;
	if ((upstream) && (!downstream))
		overlapping = 1;
	if ((!upstream) && (downstream))
		overlapping = 2;
	if ((upstream) && (downstream))
		overlapping = 3;

	return overlapping;
	// 0: not overlapping
	// 1: downstream basin exists
	// 2: upstream basin exists
	// 3: upstream and downstream basin exists
}

short cbasinClass::getStationNumber(int cellNumber)
{
	// return value is 0 if cellNumber does not beint to 
	// an existing station; otherwise return value is number
	// of the station

	short actualCalibStation = 0;

	for (short n = 1; n <= numberOfCalibBasins; n++) {
		if (cellNumber == cellNum[n - 1]) {
			actualCalibStation = n;
			break;
		}
	}
	return actualCalibStation;
}
