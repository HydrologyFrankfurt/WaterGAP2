
/***********************************************************************
*
* see former changes at file calibration.h.versioninfos.txt
*
***********************************************************************/
#if !defined (_calibration_h_)
#define _calibration_h_

class calibGammaClass {
  public:
	calibGammaClass(void);
	int cellCorrFactorInd = 1;
	int gammaCond = 0;
	int calibStatus = 0; // 1 only gamma, 2 gamma and 10 perc Qobsmod, 3 gamma, 10 perc Qobsmod and CFA, 4 gamma, 10 perc Qobsmod, CFA and CFS
	void init();
	void prepareFiles();
	short checkStations(float *gamma);
	void readObservedData();
	void setRunoff(int year, float value);
	void setUpstInflow(int year, float value);
	void setWaterUse(int year, float value);
	float findNewGamma(float gamma_old);
	void writeCorrFactors(float gamma, int cellCorrFactorInd);
	void writeCalibStatus(int calibStatus);
	float getCellCorrFactor();
    short getCallCounter();
    int getCalibStatus();
	void createCorrectionGrid(short startYear, short endYear,
							  float simulatedDischarge, float measuredDischarge);

    short CallCounterNo;

  private:
	char logFileName[250];
	char resultFileName[250];
	char statusFileName[250];

	short calibStationNumber;
	short years[100];
	float measuredRunoff[100];
	float simulatedRunoff[100];
	float simulatedInflow[100];
	float simulatedWaterUse[100];

	float gammaUpperLimit, gammaLowerLimit;
	float gammaOfPreviousRun;
    float cellCorrFactor = 1.;

};

template < class T > int sgn(T value)
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

#endif
