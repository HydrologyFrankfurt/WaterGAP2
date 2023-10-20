#include <cstdio>
#include <cstdlib>
#include "timestring.h"
#include "timestring.h"
#include "def.h"
#include "globals.h"
#include "common.h"


using namespace std;


upstreamStationClass::upstreamStationClass()
{
	upstreamStationList = NULL;
	numberOfUpstreamStations = NULL;
	numberOfBasins = 0;
}

void upstreamStationClass::deleteAllListElements()
{
	if (upstreamStationList != NULL) {
		listElement *temp;
		listElement *temp2;

		// go through the list, and delete all the elements
		// until an element is reached which does
		// not have a successor
		for (int n = 1; n <= numberOfBasins; n++) {
			if (upstreamStationList[n - 1] != NULL) {
				temp = upstreamStationList[n - 1]->next;
				delete upstreamStationList[n - 1];

				while (temp != NULL) {
					temp2 = temp->next;
					delete temp;

					temp = temp2;
				}
			}
		}
	}
}

upstreamStationClass::~upstreamStationClass()
{
	deleteAllListElements();
	delete[]upstreamStationList;
	delete[]numberOfUpstreamStations;
}

void upstreamStationClass::findAllStations(const std::string routing_dir)
{
	findStations(routing_dir, 1);
}

void upstreamStationClass::findDirectStations(const std::string routing_dir)
{
	// direct stations are those, which can be reached through the river
	// network without crossing a cell which has another station
	findStations(routing_dir, 0);
}

void upstreamStationClass::findStations(const std::string routing_dir, short allStationOption)
{
	char filename[250];

	VariableChannelGrid<9,int> G_inflow_cells;

	G_inflow_cells.read(routing_dir + "/G_INFLC.9.UNF4");

	int n, cbasinCellNum, cb2, inflowCell, cellsToBeDone;
	short i, cellInList;

	// if the subroutine is called more than once, 
	// it will start with an empty list each time
	deleteAllListElements();
	delete[]numberOfUpstreamStations;
	delete[]upstreamStationList;

	numberOfBasins = cbasin.numberOfCalibBasins;

	upstreamStationList = new listElement *[numberOfBasins];
	for (n = 0; n <= numberOfBasins - 1; n++)
		upstreamStationList[n] = NULL;

	Grid<char> G_done;

	// 0: nothing has been done with the cell up to now
	// 1: cell has to be check within next step
	// 2: cell is done
	for (n = 0; n <= ng - 1; n++)
		G_done[n] = 0;

	listElement *temp2;
	int cb;

	for (cb = 1; cb <= numberOfBasins; cb++) {
		cbasinCellNum = cbasin.cellNum[cb - 1];
		G_done[cbasinCellNum - 1] = 1;
		do {
			for (n = 0; n <= ng - 1; n++) {
				if (1 == G_done[n]) {
					// cell is marked to be done within this step
					G_done[n] = 2;
					for (i = 0; i <= 8; i++) {
						// go to all direct upstream cells
						inflowCell = G_inflow_cells(n,i);
						if (inflowCell != 0) {
							// is the cell in the list of cells with stations?
							cellInList = 0;
							for (cb2 = 1; cb2 <= numberOfBasins; cb2++) {
								if (inflowCell == cbasin.cellNum[cb2 - 1]) {
									cellInList = 1;
									break;
								}
							}
							if (1 == cellInList) {
								// if cell contains a station, 
								// append it to the list (as a new element)
								//cout << cb << ' ' << cb2 << endl;
								G_done[inflowCell - 1] = allStationOption;
								// if =1 is used all upstream stations are found
								// if =0 is used only the direct upstream stations are found
								listElement *temp = new listElement;

								temp->station = cb2;
								temp->next = NULL;
								if (NULL == upstreamStationList[cb - 1])
									upstreamStationList[cb - 1] = temp;
								else {
									temp2 = upstreamStationList[cb - 1];
									upstreamStationList[cb - 1] = temp;
									upstreamStationList[cb - 1]->next = temp2;
								}
							} else
								G_done[inflowCell - 1] = 1;
						}
					}
				}
			}
			// check if cells are left 
			// which have to be done in the next step
			cellsToBeDone = 0;
			for (n = 0; n <= ng - 1; n++)
				if (1 == G_done[n])
					cellsToBeDone++;
		}
		while (cellsToBeDone > 0);
	}

	// create array with number of upstream stations 
	short counter;
	numberOfUpstreamStations = new short[numberOfBasins];

	for (cb = 1; cb <= numberOfBasins; cb++) {
		counter = 0;
		if (upstreamStationList[cb - 1] != NULL) {
			counter++;
			//cout << cb << ' ' << upstreamStationList[cb-1]->station << endl;
			temp2 = upstreamStationList[cb - 1]->next;
			while (temp2 != NULL) {
				counter++;
				//cout << cb << ' ' << temp2->station << endl; 
				temp2 = temp2->next;
			}
		}
		numberOfUpstreamStations[cb - 1] = counter;
	}
}

short upstreamStationClass::getNumberOfUpstreamStations(int stn)
{
	if ((stn < 1) || (stn > numberOfBasins)) {
		cerr << "Station number out of range in upstreamStationClass: " << stn << endl;
		exit(-1);
	}
	return numberOfUpstreamStations[stn - 1];
}

short upstreamStationClass::getUpstreamStation(int stn, short i)
{

	listElement *temp;
	short counter;

	if ((stn < 1) || (stn > numberOfBasins)) {
		cerr << "Station number out of range in upstreamStationClass: " << stn << endl;
		exit(-1);
	}
	if ((i < 1) || (i > getNumberOfUpstreamStations(stn))) {
		cerr << "Out of range in upstreamStationClass: " << i << endl;
		exit(-1);
	}
	if (1 == i)
		return upstreamStationList[stn - 1]->station;
	else {
		temp = upstreamStationList[stn - 1];
		for (counter = 2; counter <= i; counter++) {
			temp = temp->next;
		}
		return temp->station;
	}
}

void upstreamStationClass::writeListToFile(const std::string filename, const std::string filename_no)
{
	ofstream outputFile(filename);
	ofstream outputFile_no(filename_no);

	if (!outputFile) {
		cerr << "Can not open file " << filename << " for writing\n";
		exit(-1);
	}
	if (!outputFile_no) {
		cerr << "Can not open file " << filename_no << " for writing\n";
		exit(-1);
	}
	outputFile << "# " << getTimeString();
	outputFile << endl;

	outputFile_no << "# " << getTimeString();
	outputFile_no <<"calibBasinNumber\tgcrcNumber\tNumberOfUpstreamStations"<< endl;
	short i;

	for (int n = 1; n <= numberOfBasins; n++) {
		outputFile << n << ":";
		outputFile_no << n <<'\t'<<cbasin.cellNum[n - 1]<<'\t'<<getNumberOfUpstreamStations(n)<<endl;
		for (i = 1; i <= getNumberOfUpstreamStations(n); i++)
			outputFile << ' ' << getUpstreamStation(n, i);
		outputFile << endl;
	}
	outputFile.close();
	outputFile_no.close();
}

void upstreamStationClass::replaceGridValues(int st, float newValue, Grid<> & Gridp, std::vector<float> upValueList)
{
	for (int m = 0; m <= ng - 1; m++)
		if (G_sbasin[m] == st) {
			Gridp[m] = newValue;
		}
	// go through the list of upstream stations.
	// if no value is given use the same, otherwise ignore.
	// negative values are considered to represent 'nodata'.
	int upSt;

	for (int n = 1; n <= getNumberOfUpstreamStations(st); n++) {
		upSt = getUpstreamStation(st, n);
		if (upValueList[upSt - 1] < 0)
			upstreamStationClass::replaceGridValues(upSt, newValue, Gridp, upValueList);
	}
}
