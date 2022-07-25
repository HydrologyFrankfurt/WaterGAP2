
/***********************************************************************
* Gamma was changed from 0.3 - 3 (2.1f, 2.1g, 2.1h) to 0.1 - 5 since 2.2 WATCH
*
* A description of the effect of additional discharge measurement stations is
* from Hunger, M., and P. Döll (2008), Value of river discharge data for
* global-scale hydrological modeling, Hydrology and Earth System Sciences,
* 12(3), 841-861, doi:10.5194/hess-12-841-2008. [online] Available from:
* http://www.hydrol-earth-syst-sci.net/12/841/2008/
*
* There is also a Word document of the calibration procedure (of 2.1f) written
* by Martin Hunger ("WaterGAP 2.1f Neuerungen, Kalibrierung, Regionalisierung")
* in 2006. A documetation of calibration process of 2.2 refvers is in
* preparation
*
* see former changes at file calibration.cpp.versioninfos.txt
*
***********************************************************************/

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <fstream>

#include "timestring.h"
#include "calib_basins.h"
#include "upstream_stations.h"
#include "geo.h"
#include "option.h"
#include "gridio.h"	// for cell based correction factors (since 2.1b)
#include "daily.h"

#include "calibration.h"
using namespace std;

extern optionClass options;

template < class T > int sign(T value)
{
	// returns sign of value 
	if (value > 0)
		return 1;
	else {
		if (value == 0)
			return 0;
		else
			return -1;
	}
}

calibGammaClass::calibGammaClass()
{
	//const float gammaUpperLimit = 31.9;
	gammaUpperLimit = 5.;
	//const float gammaLowerLimit = 0.00002;
	gammaLowerLimit = 0.1;
}

void calibGammaClass::prepareFiles()
{
	FILE *file_ptr;

	sprintf(resultFileName, "CALIBRATION.OUT");
	file_ptr = fopen(resultFileName, "w");
	if (file_ptr != NULL) {
		fprintf(file_ptr, "# %s", getTimeString());
	} else {
		cerr << "Can not open " << resultFileName << " for writing.\n";
		exit(-1);
	}
	fclose(file_ptr);

	sprintf(logFileName, "CALIBRATION.LOG");
	file_ptr = fopen(logFileName, "w");
	if (file_ptr != NULL) {
		fprintf(file_ptr, "# %s", getTimeString());
	} else {
		cerr << "Can not open " << logFileName << " for writing.\n";
		exit(-1);
	}
	fclose(file_ptr);

        sprintf(statusFileName, "CALIBSTATUS.OUT");
        file_ptr = fopen(statusFileName, "w");
        if (file_ptr != NULL) {

        } else {
                cerr << "Can not open " << statusFileName << " for writing.\n";
                exit(-1);
        }
        fclose(file_ptr);
}

short calibGammaClass::checkStations(float *gamma)
{
	// return value: 
	// 0: no calibration run 
	// otherwise: number of the station to calibrate

	extern geoClass geo;
	extern cbasinClass cbasin;
	extern upstreamStationClass allUpstSt;

	// open file which contains location of calibration station.
	// if file 'C_STATION.DAT' does not exist 
	// the model assumes that this is not a calibration run.
	// in that case normal simulation will be done.
	float gammaFromFile;

	char filename[250];

	sprintf(filename, "C_STATION.DAT");
	ifstream station_file(filename);

	if (!station_file)
		return 0;	// file not found: no calibration run
	else {
		// read data from file:
		// 1. name of the station
		// 2. longitude 
		// 3. latitude
		// (4. start value of gamma for calibration (optional)) 
		float longitude, latitude;
		int cellNumber;
		char string[250];

		station_file >> string;
		if (station_file) {
			// values for int./lat. have to be given in the file
			station_file >> longitude >> latitude;

			// start value for gamma if optional 
			station_file >> gammaFromFile;
			if (station_file.fail())
				gammaFromFile = -99;
		} else {
			cerr << "Problem while reading from file " << filename << endl;
			exit(-1);
		}

		// find station
		cellNumber = geo.cellNumByLonLat(longitude, latitude);
		if (cellNumber <= 1) {
			cerr << "Calibration station not part of the landmask!\n";
			exit(-1);
		}
		calibStationNumber = cbasin.getStationNumber(cellNumber);
		// 0: no station beints to the cell.
		if (0 == calibStationNumber) {
			cerr << "Calibration station is not available!" << endl;
			exit(-1);
		}
	}

	// this point is will only be reached if a claibration station
	// has been found.
	prepareFiles();

	ofstream logFile(logFileName, ios::app);
	ofstream resultFile(resultFileName, ios::app);

	logFile << "Station No. " << calibStationNumber << " selected for calibration." << endl;


	// are all the upstream basins calibrated ?

	short allCalibrated = 1;
	int upSt;

	for (int n = 1; n <= allUpstSt.getNumberOfUpstreamStations(calibStationNumber); n++) {
		upSt = allUpstSt.getUpstreamStation(calibStationNumber, n);
		if (cbasin.gamma[upSt - 1] < 0) {
			allCalibrated = 0;
			cerr << "station " << cbasin.name[upSt - 1]<< " is not calibrated! \n";
		}
	}

	if (!allCalibrated) {
		logFile << "Not all upstream stations are calibrated!\n";
		resultFile << "NOT_YET_CALIBRATED\n";
		cerr << "Not all upstream stations are calibrated!\n";
		logFile.close();
		resultFile.close();
		exit(-1);
	} else {
		// write file header
		resultFile << "# 1. column: value of gamma\n";
		resultFile << "# 2. column: efficiency criterion\n";
		resultFile << "# 3. column: sum of diff. betw. sim. and meas. values\n";
		resultFile << "# 4. column: lower limit for gamma\n";
		resultFile << "# 5. column: upper limit for gamma\n";
		resultFile << "# 6. column: average of measured values\n";
		resultFile << "# 7. column: number of years with measured data\n";
		resultFile << "# 8. column: resulting correction factor\n";
		resultFile << "# 9. column: average inflow from upstream\n";
		resultFile << "#10. column: average satisfied water use of the basin\n";
		resultFile << "#11. column: (obsolete)\n";
		resultFile << "#12. column: simulated runoff at station\n";
		resultFile << "#13. column: runoff generated within the basin\n";
		resultFile << "#14. column: new correction factor (basinwide average)\n";
	}

	logFile.close();
	resultFile.close();

	// check start value of gamma
	if (gammaFromFile < 0)
		*gamma = 1;
	else {
		if ((gammaFromFile > gammaUpperLimit)
			|| (gammaFromFile < gammaLowerLimit))
			*gamma = 1;
		else
			*gamma = gammaFromFile;
	}

	return calibStationNumber;
}

void calibGammaClass::readObservedData()
{
	// read RIVER.DAT 
	// file with two columns: year, measured runoff (m3/s) 
	// years do not necessarily need to be continuously
	char filename[250];

	sprintf(filename, "RIVER.DAT");
	ifstream observedDataFile(filename);

	if (!observedDataFile) {
		cerr << "Error while opening file " << filename << endl;
		exit(-1);
	} else {
		int i;
		float f;

		while (!observedDataFile.eof()) {
			observedDataFile >> i;
			observedDataFile >> f;
			measuredRunoff[i - options.evalStartYear]
				= f * 365 * 24 * 60 * 60 / 1000000000.0;
			// conversion from m3/s to km3/year
		}
	}
}

void calibGammaClass::init()
{
	for (int i = 0; i <= (options.end_year - options.evalStartYear); i++) {
		years[i] = options.evalStartYear + i;
		measuredRunoff[i] = -99;
		simulatedRunoff[i] = -99;
		simulatedInflow[i] = -99;
		simulatedWaterUse[i] = -99;
	}
	readObservedData();
}

void calibGammaClass::setRunoff(int year, float value)
{
	simulatedRunoff[year - options.evalStartYear] = value;
}

void calibGammaClass::setUpstInflow(int year, float value)
{
	simulatedInflow[year - options.evalStartYear] = value;
}

void calibGammaClass::setWaterUse(int year, float value)
{
	simulatedWaterUse[year - options.evalStartYear] = value;
}

float calibGammaClass::findNewGamma(float gamma_old)
{
	// compare simulated and measured values and return a value
	// for gamma, which will be used for the next calibration step

	ofstream logFile(logFileName, ios::app);
	ofstream resultFile(resultFileName, ios::app);

	int i;
	float measuredRunoffAvg;
        float measuredRunoffAvgAdapt;
	static float sumOfDifferences_old = 0;
	float sum1, sum2;
        static short callCounter = 0;

	float gamma;
	short evalStartYear, end_year;
	static float gamma_low = -99;
	static float gamma_high = -99;

        cout << "Executing calibration routine for gamma..." << endl;
	callCounter++;
        CallCounterNo = callCounter;

	evalStartYear = options.evalStartYear;
	end_year = options.end_year;

	float sumOfDifferences = 0;
	float measuredRunoffSum = 0;
	float simInflowSum = 0;
	float simWaterUseSum = 0;
	float simRunoffSum = 0;
	int numberOfYears = 0;

	logFile << "# Year\tObs. Runoff\tSim. Runoff\tUpst. Inflow" << "\tWater Use[km3/year]\n";
	for (i = 0; i <= (end_year - evalStartYear); i++) {
		logFile << years[i] << '\t'
			<< measuredRunoff[i] << '\t'
			<< simulatedRunoff[i] << '\t'
			<< simulatedInflow[i] << '\t' << simulatedWaterUse[i] << endl;

		if (measuredRunoff[i] > -1) {
			numberOfYears++;
			sumOfDifferences += simulatedRunoff[i] - measuredRunoff[i];
			measuredRunoffSum += measuredRunoff[i];
			simRunoffSum += simulatedRunoff[i];
			simInflowSum += simulatedInflow[i];
			simWaterUseSum += simulatedWaterUse[i];
		}
	}

        measuredRunoffAvg = measuredRunoffSum / numberOfYears;

	if (1 == callCounter) {
            // if already here everything is fine, quit calibration.
                    // termination condition:
                    // is difference between simulated and observed discharged
                    // less than one percent of observed discharge?
            logFile << "1% criterion: sumOfDifferences: " << sumOfDifferences << "\n";
            logFile << "1% criterion: noYears * Qobs: " << numberOfYears * measuredRunoffAvg << "\n";
            logFile << "1% criterion: calc: " << (fabs(sumOfDifferences / (numberOfYears * measuredRunoffAvg))) << "\n";
                    if ((fabs(sumOfDifferences / (numberOfYears * measuredRunoffAvg))) < 0.01) {
                            logFile << "Termination condition is fulfilled with +-1% condition.\n";
                            gamma = -99;
                            gammaCond = 1; // gamma value found, quit calibration.
                            calibStatus = 1;
                    }
                    else { // if 1% criterion is not fulfilled
		// when this routine is called for the first time, 
		// an upper and lower limit for gamma is defined.
		if (sumOfDifferences > 0) {
                        logFile << "First call: Value of gamma is " << gamma_old << " and to small.\n";
			if (gamma_old >= gammaUpperLimit) {
				gamma = -99;
                                logFile << "gamma is set to: " << gamma << "\n";
			} else {
				gamma_low = gamma_old;
				gamma = gammaUpperLimit;
				//gamma = gamma_old * 2.0;
                                logFile << "gamma is set to: " << gamma << "\n";
                                logFile << "gamma_low is set to: " << gamma_low << "\n";
                                logFile << "gamma_high is still: " << gamma_high << "\n";
			}
		} else {
                        logFile << "First call: Value of gamma is " << gamma_old << " and to great." << "\n";
			if (gamma_old <= gammaLowerLimit) {
				gamma = -99;
                                logFile << "gamma is set to: " << gamma << "\n";
			} else {
				gamma_high = gamma_old;
				//gamma = gamma_old/2.0;
				gamma = gammaLowerLimit;
                                logFile << "gamma is set to: " << gamma << "\n";
                                logFile << "gamma_high is set to: " << gamma_high << "\n";
                                logFile << "gamma_low is still: " << gamma_low << "\n";
                        }
			}
		}


        } else if (2 == callCounter) {
            // for second time, gamma is in every case either max or min.
            // this does not mean that gamma cannot be in between.
            // check first, if 1% condition is fulfilled.

            // if already here everything is fine, quit calibration.
                    // termination condition:
                    // is difference between simulated and observed discharged
                    // less than one percent of observed discharge?
                logFile << "1% criterion: sumOfDifferences: " << sumOfDifferences << "\n";
                logFile << "1% criterion: noYears * Qobs: " << numberOfYears * measuredRunoffAvg << "\n";
                logFile << "1% criterion: calc: " << (fabs(sumOfDifferences / (numberOfYears * measuredRunoffAvg))) << "\n";
                    if ((fabs(sumOfDifferences / (numberOfYears * measuredRunoffAvg))) < 0.01) {
                            logFile << "Termination condition is fulfilled with +-1% condition.\n";
                            gamma = -99;
                            gammaCond = 1; // gamma value found, quit calibration.
                            calibStatus = 1;
                    }
                    else { // modify gamma
                // when this is called the interval of the previous step
                // is devided into two halfs and we continue with a smaller
                // interval.
                if (sumOfDifferences > 0) {
                        logFile << "Value of gamma is to small.\n";
                        logFile << "Second call: Value of gamma is " << gamma_old << " and to small." << "\n";
                        //cout << "Value of gamma is to small." << endl;
                        if (gamma_old == gammaUpperLimit) {
                                gamma = gamma_old;
                                //gamma_low = gammaUpperLimit;
                                gamma_high = gammaUpperLimit;
                                //gamma_high = -99;
                                logFile << "but not yet terminating calibration.\n";
                                logFile << "gamma remains at " << gamma << "\n";
                                logFile << "but not yet terminating calibration." << "\n";

                        } else if (gamma_old < gammaUpperLimit){
                                if (gamma_high < 0) {
                                        // no upper limit has been found up to now
                                        gamma_low = gamma_old;
                                        gamma = gamma_old * 2.0;
                                        logFile << "gamma is set to: " << gamma << "\n";
                                        logFile << "gamma_high is still: " << gamma_high << "\n";
                                        logFile << "gamma_low is set to: " << gamma_low << "\n";
                                } else {
                                        gamma_low = gamma_old;
                                        if (gamma != gammaUpperLimit)
                                            gamma = (gamma_low + gamma_high) / 2.0;
                                        logFile << "gamma is set to: " << gamma << "\n";
                                        logFile << "gamma_high is still: " << gamma_high << "\n";
                                        logFile << "gamma_low is set to: " << gamma_low << "\n";
                                }
                        }
                } else {
                        logFile << "Value of gamma is to great.\n";
                        logFile << "Second call: Value of gamma is " << gamma_old << " and to great." << "\n";
                        if (gamma_old == gammaLowerLimit) {
                                gamma = gamma_old;
                                gamma_low = gammaLowerLimit;
                                logFile << "but not yet terminating calibration.\n";
                                logFile << "gamma remains at " << gamma << "\n";
                                logFile << "but not yet terminating calibration." << "\n";
                        } else if (gamma_old > gammaLowerLimit){
                        if (gamma_low < 0) {
                                gamma = gamma_old / 2.0;
                                gamma_high = gamma_old;
                                logFile << "gamma is set to: " << gamma << "\n";
                                logFile << "gamma_high is set to: " << gamma_high << "\n";
                                logFile << "gamma_low is still: " << gamma_low << "\n";
	} else {
                                gamma_high = gamma_old;
                                if (gamma != gammaLowerLimit)
                                    gamma = (gamma_low + gamma_high) / 2.0;
                                logFile << "gamma is set to: " << gamma << "\n";
                                logFile << "gamma_high is set to: " << gamma_high << "\n";
                                logFile << "gamma_low is still: " << gamma_low << "\n";
                        }
                }
                        }
                        }
        }    else { // callCounter > 2
            // if gamma is still max or min, terminate calibration and modify Qobs

            // if already here everything is fine, quit calibration.
                    // termination condition:
                    // is difference between simulated and observed discharged
                    // less than one percent of observed discharge?
            logFile << "1% criterion: sumOfDifferences: " << sumOfDifferences << "\n";
            logFile << "1% criterion: noYears * Qobs: " << numberOfYears * measuredRunoffAvg << "\n";
            logFile << "1% criterion: calc: " << (fabs(sumOfDifferences / (numberOfYears * measuredRunoffAvg))) << "\n";
                    if ((fabs(sumOfDifferences / (numberOfYears * measuredRunoffAvg))) < 0.01) {
                            logFile << "Termination condition is fulfilled with +-1% condition.\n";
                            gamma = -99;
                            gammaCond = 1; // gamma value found, quit calibration.
                            calibStatus = 1;
                    }
                    else {
		// when this is called the interval of the previous step
		// is devided into two halfs and we continue with a smaller
		// interval.
		if (sumOfDifferences > 0) {
			logFile << "Value of gamma is to small.\n";
                                   logFile << callCounter << "th call: Value of gamma is " << gamma_old << " and to small." << "\n";
			if (gamma_old == gammaUpperLimit) {
				gamma = -99;
                                        //gamma_high = -99;
                                        logFile << "still gamma to small, therefore terminating calibration.\n";
                                        logFile << "gamma gets " << gamma << "\n";
                                        logFile << "and search for gamma is terminating." << "\n";
                                        logFile << "now try if Qsim is within 10% of Qobs." << "\n";
                                        gammaCond = 2; // modify Qsim by 10%
                                } else if (gamma_old < gammaUpperLimit){
				if (gamma_high < 0) {
					// no upper limit has been found up to now
					gamma_low = gamma_old;
					gamma = gamma_old * 2.0;
                                                logFile << "gamma is set to: " << gamma << "\n";
                                                logFile << "gamma_high is still: " << gamma_high << "\n";
                                                logFile << "gamma_low is set to: " << gamma_low << "\n";
				} else {
					gamma_low = gamma_old;
					gamma = (gamma_low + gamma_high) / 2.0;
                                                logFile << "gamma is set to: " << gamma << "\n";
                                                logFile << "gamma_high is still: " << gamma_high << "\n";
                                                logFile << "gamma_low is set to: " << gamma_low << "\n";
				}
			}
		} else {
			logFile << "Value of gamma is to great.\n";
                                logFile << callCounter << "th call: Value of gamma is " << gamma_old << " and to great." << "\n";
                                if (gamma_old == gammaLowerLimit) {
                                    gamma = -99;
                                    //gamma_high = -99;
                                    logFile << "still gamma to great, therefore terminating calibration.\n";
                                    logFile << "gamma gets " << gamma << "\n";
                                    logFile << "and search for gamma is terminating." << "\n";
                                    logFile << "now try if Qsim is within 10% of Qobs." << "\n";
                                    gammaCond = 2; // modify Qsim by 10%
                            } else if (gamma_old > gammaLowerLimit){

			if (gamma_low < 0) {
				gamma = gamma_old / 2.0;
				gamma_high = gamma_old;
                                        logFile << "gamma is set to: " << gamma << "\n";
                                        logFile << "gamma_high is set to: " << gamma_high << "\n";
                                        logFile << "gamma_low is still: " << gamma_low << "\n";
			} else {
				gamma_high = gamma_old;
				gamma = (gamma_low + gamma_high) / 2.0;
                                        logFile << "gamma is set to: " << gamma << "\n";
                                        logFile << "gamma_high is set to: " << gamma_high << "\n";
                                        logFile << "gamma_low is still: " << gamma_low << "\n";
                                }
			}
		}
	}

                    }

        //measuredRunoffAvg = measuredRunoffSum / numberOfYears;

	float NashSutclCoeff;

	if (end_year - evalStartYear > 0) {
		// Calculation of 'efficiency criterion': Nash-Sutcliff-Coefficient
		// This can only be used if more than one year is simulated
		sum1 = 0;
		sum2 = 0;
		for (i = 0; i <= (end_year - evalStartYear); i += 1) {
			if (measuredRunoff[i] > -1) {
				sum1 += (measuredRunoff[i] - measuredRunoffAvg)
					* (measuredRunoff[i] - measuredRunoffAvg);
				sum2 += (simulatedRunoff[i] - measuredRunoff[i])
					* (simulatedRunoff[i] - measuredRunoff[i]);
			}
		}
		NashSutclCoeff = (sum1 - sum2) / sum1;
	} else
		NashSutclCoeff = -99;




        // if gamma reaches upper or lower limit (still after callCounter 1) and 1% condition is not fulfilled,
        // vary Qobs by 10% and check again
        // in case that measuredRunoffAvg is not adapted, store value (for resultFile)
        measuredRunoffAvgAdapt = measuredRunoffAvg;
        logFile << "simRunoffSum / numberOfYears" << simRunoffSum / numberOfYears << ".\n";
        logFile << "callCounter: " << callCounter << "\n";
        if ((gamma_old == gammaUpperLimit) && (gammaCond == 2)) {
                        measuredRunoffAvgAdapt = measuredRunoffAvg * 1.1;
                        logFile << "measuredRunoffAvg is multiplied by 1.1 (+10%).\n";
                        if ((simRunoffSum / numberOfYears) < measuredRunoffAvgAdapt) {
                            logFile << "10% criterion is fulfilled, CFA is still " << cellCorrFactor << ".\n";
                            calibStatus = 2;
                        }
                        else {
                            cellCorrFactorInd = 99.;
                            logFile << "10% criterion cannot be fulfilled, CFA is calculated.\n";
                        }
                }
                if ((gamma_old == gammaLowerLimit) && (gammaCond == 2)){
                        measuredRunoffAvgAdapt = measuredRunoffAvg * 0.9;
                        logFile << "measuredRunoffAvg is multiplied by 0.9 (-10%).\n";
                        if ((simRunoffSum / numberOfYears) > measuredRunoffAvgAdapt) {
                            logFile << "10% criterion is fulfilled, CFA is still " << cellCorrFactor << ".\n";
                            calibStatus = 2;
                        }
                        else {
                            cellCorrFactorInd = 99.;
                            logFile << "10% criterion cannot be fulfilled, CFA is calculated.\n";
		}
	}

	float runoffGeneratedInBasin;
// calculate CFA only if is set to 99 (10% criterion is not enough) using the adapted Qobs
logFile << "old cellCorrFactor " << cellCorrFactor << "\n";
logFile << "cellCorrFactorInd " << cellCorrFactorInd << "\n";
logFile << "measuredRunoffAvgAdapt " << measuredRunoffAvgAdapt << "\n";

        if (cellCorrFactorInd == 99.) {
            logFile << "CFA is calculatd here: " << "\n";
            cellCorrFactor = (measuredRunoffAvgAdapt + (simWaterUseSum - simInflowSum) / numberOfYears)
		/ ((simRunoffSum + simWaterUseSum - simInflowSum) / numberOfYears);
	runoffGeneratedInBasin = (simRunoffSum / numberOfYears)
		+ (simWaterUseSum / numberOfYears)
		- (simInflowSum / numberOfYears);
            logFile << "CFA is calculated to: " << cellCorrFactor << "\n";
            logFile << "runoffGeneratedInBasin " << runoffGeneratedInBasin << "\n";
	/* if (((simInflowSum/numberOfYears) > measuredRunoffAvg) ||  
	 * (runoffGeneratedInBasin <= 0)) cellCorrFactor = 1;
	 * else cellCorrFactor = cellCorrFactor_1; 
	 */
        }

	resultFile << gamma_old << '\t'
		<< NashSutclCoeff << '\t'
		<< sumOfDifferences << '\t'
		<< gamma_low << '\t'
		<< gamma_high << '\t'
                << measuredRunoffAvgAdapt << '\t'
		<< numberOfYears << '\t'
                << (measuredRunoffAvgAdapt
			/ (simRunoffSum / numberOfYears)) << '\t'
		<< simInflowSum / numberOfYears << '\t'
		<< simWaterUseSum / numberOfYears << '\t'
		<< 0 << '\t'
		<< simRunoffSum / numberOfYears << '\t'
		<< runoffGeneratedInBasin << '\t' << cellCorrFactor << '\t' << endl;

        if ((gamma > 0) && (gamma < (0.99 * gammaLowerLimit)) && (callCounter > 2)) {
		logFile << "Terminating calibration.\n";
		logFile << "Gamma has reached lower limit.\n";
		gamma = -99;
	}
        if ((gamma > 0) && //(gamma < 0.001) && // HMS 2015-08-05 for rare cases where gamma search would take too much time
		(fabs((sumOfDifferences_old / sumOfDifferences) - 1) < 0.0001)) {
		logFile << "Difference between obs. and sim. discharge has not changed\n";
		logFile << "since the last calibration step.\n";
                logFile << "Terminating therefore.\n";
		gamma = -99;
	}
	sumOfDifferences_old = sumOfDifferences;


	if ((gamma < 0) && (fabs(cellCorrFactor - 1) > 0.01)) {
                logFile << "Creating new grid with correction factors for cells" << "\n";
		createCorrectionGrid(evalStartYear, end_year,
                                                         simRunoffSum / numberOfYears, measuredRunoffAvgAdapt);
                calibStatus = 3;
        }

	logFile.close();
	resultFile.close();
	gammaOfPreviousRun = gamma_old;

	return (gamma);
}

float calibGammaClass::getCellCorrFactor()
{
	return cellCorrFactor;
}

short calibGammaClass::getCallCounter()
{
        return CallCounterNo;
}

int calibGammaClass::getCalibStatus()
{
        return calibStatus;
}

void calibGammaClass::writeCalibStatus(int calibStatus)
{
        ofstream statusFile(statusFileName, ios::app);
        statusFile << calibStatus << endl;
        statusFile.close();

}
void calibGammaClass::writeCorrFactors(float gamma, int cellCorrFactInd)
// here, gamma and CFA are needed to decide if simulatedRunoff needs to be multiplied by a factor or not
{
	float measuredRunoffSum = 0;
	float simRunoffSum = 0;
	int numberOfYears = 0;
        float cfs = 1;

	for (short i = 0; i <= (options.end_year - options.evalStartYear); i++) {
		if (measuredRunoff[i] > -1) {
			numberOfYears++;
			measuredRunoffSum += measuredRunoff[i];
			simRunoffSum += simulatedRunoff[i];
		}
	}
	        ofstream logFile(logFileName, ios::app);
        // adapt only, if gamma at limits and 10% criterion is not enough
        logFile << "gamma for CFS calculation is: " << gamma << "\n";
		logFile << "cellCorrFactInd for CFS calculation is: " << cellCorrFactInd << "\n";
        if ((gamma == gammaUpperLimit) && (cellCorrFactInd == 99)) {
            measuredRunoffSum = measuredRunoffSum * 1.1;
			logFile << "multiply measuredRunoffSum by 1.1" << "\n";
        }
        if ((gamma == gammaLowerLimit) && (cellCorrFactInd == 99)) {
            measuredRunoffSum = measuredRunoffSum * 0.9;
			logFile << "multiply measuredRunoffSum by 0.9" << "\n";
        }
        //calculate CFS if 10% criterium not fulfilled
        if (cellCorrFactInd == 99) {
            cfs = measuredRunoffSum / simRunoffSum;
			logFile << "measuredRunoffSum: " << measuredRunoffSum << "\n";
			logFile << "simRunoffSum: " << simRunoffSum << "\n";
		}
        else
            cfs = 1.0;
        logFile << "cfs1 " << cfs;
        //and set to 1.0 if +-1% deviation (same criterion like gamma in first try)
        if (((cfs > 1.0) && (cfs < 1.01)) || ((cfs < 1.0) && (cfs > 0.99))) {
            logFile << "cfs is not far away from 1.0: " << cfs << "\n";
            cfs = 1.0;
            logFile << "and is set to: " << cfs << "\n";
            }
        if ((cfs > 1.0) || (cfs < 1.0))
            calibStatus = 4;
        logFile << "cfs2 " << cfs;
        logFile.close();

	char corrFactFileName[250];

	sprintf(corrFactFileName, "STAT_CORR_FACTOR.OUT");
	ofstream corrFactFile(corrFactFileName);

	if (corrFactFile) {
		corrFactFile << "# Measured runoff\tSimulated runoff\tGamma\t"
			<< "Cell correction factor\tStation correction factor" << endl;
		corrFactFile << measuredRunoffSum / numberOfYears << '\t'
			<< simRunoffSum / numberOfYears << '\t'
			<< gammaOfPreviousRun << '\t'
                        << cellCorrFactor << '\t' << cfs << endl;
		corrFactFile.close();
	}
}

void calibGammaClass::createCorrectionGrid(short startYear, short endYear,
										   float simulatedDischarge, float measuredDischarge)
{
	extern gridioClass gridIO;
	extern dailyWaterBalanceClass dailyWaterBalance;
	extern signed short G_sbasin[ng];


	// create a grid with average values of
	// total cell runoff (includes vertical water balance of open 
	// water surfaces)
	char filename[250];
    //float G_meanCellRunoff[ng];
    //float G_annualCellRunoff[ng];
    float G_meanPotCellRunoff[ng]; // HMS including CFA (and the previous concept of potCellRunoff)
    float G_annualPotCellRunoff[ng];

	int n;

	short numberOfYears = endYear - startYear + 1;
	
	// initialize local variable G_meanCellRunoff
    //for (n = 0; n < ng; n++)
    //	G_meanCellRunoff[n] = 0.;

    //for (short year = startYear; year <= endYear; year++) {
    //	sprintf(filename, "%s/G_CELL_RUNOFF_%d.UNF0", options.output_dir, year);
    //	gridIO.readUnfFile(filename, ng, G_annualCellRunoff);
    //	for (n = 0; n < ng; n++)
    //		G_meanCellRunoff[n] += (G_annualCellRunoff[n]/numberOfYears);
    //}
    //sprintf(filename, "G_CELL_RUNOFF_MEAN.UNF0");
    //gridIO.writeUnfFile(filename, ng, G_meanCellRunoff);

    for (n = 0; n <= ng - 1; n++)
        G_meanPotCellRunoff[n] = 0.0;
    for (short year = startYear; year <= endYear; year++) {
        sprintf(filename, "%s/G_POT_CELL_RUNOFF_%d.UNF0", options.output_dir, year);
        gridIO.readUnfFile(filename, ng, G_annualPotCellRunoff);
        for (n = 0; n <= ng - 1; n++)
            G_meanPotCellRunoff[n] += G_annualPotCellRunoff[n];
    }

    for (n = 0; n <= ng - 1; n++)
        G_meanPotCellRunoff[n] /= numberOfYears;
    sprintf(filename, "G_POT_CELL_RUNOFF_MEAN.UNF0");
    gridIO.writeUnfFile(filename, ng, G_meanPotCellRunoff);


	float sumOfAbsValuesForBasin = 0;

	for (n = 0; n <= ng - 1; n++) {
		if (calibStationNumber == G_sbasin[n])
            //sumOfAbsValuesForBasin += fabs(G_meanCellRunoff[n]);
            sumOfAbsValuesForBasin += fabs(G_meanPotCellRunoff[n]);
	}
	cout << "Sum of absolute values of total cell runoff: " << sumOfAbsValuesForBasin << endl;
	
	// limitation of correction factor to 0.5 - 1.5 (2.1f): sumOfAbsValuesForBasin has to be at least twice
	// as big as simulatedDischarge - measuredDischarge

	// limit correction factor to changes of +/-50% (maxmium range of
	// correction factor (0.5-1.5).
	// -> changed for 2.1f. In former versions limits were 0 and 2
	// otherwise (correction factors below 0)
	// we would result in negative local runoff on land area
	// and the sign of runoff would change, which means: an increase in 
	// precipitation would lead to a decrease in water in the river system
	for (int n = 0; n < ng; ++n) {
		if (calibStationNumber == G_sbasin[n]) {
			dailyWaterBalance.G_cellCorrFact[n]
                //= 1 - ((sign(G_meanCellRunoff[n])
                = 1 - ((sign(G_meanPotCellRunoff[n])
						* (simulatedDischarge - measuredDischarge))
					   / sumOfAbsValuesForBasin);
			if (dailyWaterBalance.G_cellCorrFact[n] > 1.5) {
				dailyWaterBalance.G_cellCorrFact[n] = 1.5;
				cout << "Limiting correction factor to range 0.5 - 1.5. Limiting to 1.5" << endl;
			}

			if (dailyWaterBalance.G_cellCorrFact[n] < 0.5) {
				dailyWaterBalance.G_cellCorrFact[n] = 0.5;
				cout << "Limiting correction factor to range 0.5 - 1.5. Limiting to 0.5" << endl;
			}
				
		}
	} // end of for(n)

	sprintf(filename, "G_CORR_FACTOR.UNF0");
	gridIO.writeUnfFile(filename, ng, dailyWaterBalance.G_cellCorrFact);
}
