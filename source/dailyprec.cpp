#include <iostream>
#include <cstdlib>
#include <cstdio>
#include "dailyprec.h"
#include "option.h"

using namespace std;

float dailyPrecNrdClass::ran1(int *idum)
{
	// (C) Copr. 1986-92 Numerical Recipes Software 7
	const int IA = 16807;
	const int IM = 2147483647;
	const double AM = (1.0 / IM);
	const int IQ = 127773;
	const int IR = 2836;
	const short NTAB = 32;
	const double NDIV = (1 + (IM - 1.0) / NTAB);
	const double EPS = 1.2e-7;
	const double RNMX = (1.0 - EPS);

	int j;
	int k;
	static int iy = 0;
	static int iv[NTAB];
	float temp;

	if (*idum <= 0 || !iy) {
		if (-(*idum) < 1)
			*idum = 1;
		else
			*idum = -(*idum);
		for (j = NTAB + 7; j >= 0; j--) {
			k = (*idum) / IQ;
			*idum = IA * (*idum - k * IQ) - IR * k;
			if (*idum < 0)
				*idum += IM;
			if (j < NTAB)
				iv[j] = *idum;
		}
		iy = iv[0];
	}
	k = (*idum) / IQ;
	*idum = IA * (*idum - k * IQ) - IR * k;
	if (*idum < 0)
		*idum += IM;
	j = (int) (iy / NDIV);
	iy = iv[j];
	iv[j] = *idum;
	if ((temp = AM * iy) > RNMX)
		return RNMX;
	else
		return temp;
}

void dailyPrecNrdClass::init(int idum, const char *output_dir)
{
	extern optionClass options;

	char rain_days[31][31];

	float unirand[31], frac, probwd, probww;
	//float precip, precip_per_day; Never used
	short wet[31];	// dry: 0, wet: 1
	short nwet, mdays=-1, i, j, k, testwet;

	FILE *file_ptr;
	char filename[250];

	for (short x = 1; x <= 3; x++) {
		if (1 == x)
			mdays = 28;
		if (2 == x)
			mdays = 30;
		if (3 == x)
			mdays = 31;

		/* if number of wet days is equal to 1     */
		/* select the mid of the month as rain day */
		for (i = 0; i <= mdays - 1; i++)
			rain_days[i][0] = 0;
		rain_days[mdays / 2][0] = 1;
		/* if number of wet days is equal to mdays */
		/* all the days are raindays               */
		for (i = 0; i <= mdays - 1; i++)
			rain_days[i][mdays - 1] = 1;

		/* all the other days */
		for (k = 2; k < mdays; k++) {
			testwet = 0;
			nwet = k;
			frac = (float) nwet / (float) mdays;
			/* probability that a wet day follows a dry day (GENG et al., 1986) */
			probwd = 0.75 * frac;
			/* probability that a wet day follows a wet day (GENG et al., 1986) */
			probww = 0.25 + probwd;

			/*printf(" nwet= %d, probwd= %f, probww= %f\n", nwet, probwd, probww); */
			for (i = 0; i < mdays; i++)
				wet[i] = 0;

			if (frac > 0.5)
				wet[0] = 1;

			/*run the same month until realization with nwet days produced */

			while (testwet != nwet) {
				testwet = wet[0];
				for (i = 1; i <= mdays - 1; i++) {
					unirand[i] = ran1(&idum);

					if (wet[i - 1] == 1) {
						if (unirand[i] < probww)
							wet[i] = 1;
						else
							wet[i] = 0;
					} else {
						if (unirand[i] < probwd)
							wet[i] = 1;
						else
							wet[i] = 0;
					}
					testwet += wet[i];
				}
				if (testwet == nwet) {
					/*printf(" run no.  %d\n", j); */
					for (i = 0; i <= mdays - 1; i++) {
						/*printf(" %d: %f %d \n", i+1, unirand[i], wet[i]); */
						rain_days[i][nwet - 1] = wet[i];
					}
				}
			}
		}

		// write results to disk
		if (options.outRainDays) { // new output options 2.1f
			sprintf(filename, "%s/RAIN_DAYS%d.DAT", output_dir, mdays);
			file_ptr = fopen(filename, "w");
			if (file_ptr != NULL) {
				// printf("Writing %s\n",filename);
				for (i = 0; i <= mdays - 1; i++) {
					for (j = 0; j <= mdays - 1; j++) {
						fprintf(file_ptr, "%d ", rain_days[i][j]);
					}
					fprintf(file_ptr, "\n");
				}
			} else {
				printf("Can not write file %s !\n", filename);
				// simulation does not depend on this
				// exit(1);
			}
			fclose(file_ptr);
		}

		// copy results into the 'private' arrays of the class.
		if (28 == mdays) {
			for (i = 0; i <= mdays - 1; i++) {
				for (j = 0; j <= mdays - 1; j++) {
					rain_days28[i][j] = rain_days[i][j];
				}
			}
		}
		if (30 == mdays) {
			for (i = 0; i <= mdays - 1; i++) {
				for (j = 0; j <= mdays - 1; j++) {
					rain_days30[i][j] = rain_days[i][j];
				}
			}
		}
		if (31 == mdays) {
			for (i = 0; i <= mdays - 1; i++) {
				for (j = 0; j <= mdays - 1; j++) {
					rain_days31[i][j] = rain_days[i][j];
				}
			}
		}
	}
}

float dailyPrecNrdClass::getDailyPrec(const short day_in_month,
									  const short month,
									  const char numberOfRainDays, const short monthlyPrecipitation)
{

	const short int n_days_per_month[12]
	= { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };
	// const short int first_day_per_month[12] 
	//  = {1,32,60,90,121,151,182,213,243,274,304,334}; 

	if (0 == numberOfRainDays)
		return 0;

	if (28 == n_days_per_month[month]) {
		if (0 == rain_days28[day_in_month - 1][numberOfRainDays - 1])
			return 0;
		else
			return monthlyPrecipitation / (float) numberOfRainDays;
	}
	if (30 == n_days_per_month[month]) {
		if (0 == rain_days30[day_in_month - 1][numberOfRainDays - 1])
			return 0;
		else
			return monthlyPrecipitation / (float) numberOfRainDays;
	}
	if (31 == n_days_per_month[month]) {
		if (0 == rain_days31[day_in_month - 1][numberOfRainDays - 1])
			return 0;
		else
			return monthlyPrecipitation / (float) numberOfRainDays;
	}
	cerr << n_days_per_month[month]
		<< "is an illegal selection for the number of days in a month.\n";
	exit(1);
	
	return -1.0; //Never reached, remove compiler warning
}

//void main(void) {
//  dailyPrecNrdClass dailyRainNrd;
//  
//  dailyPrecNrd.init(-421112,"OUT_TEST");
//  
//  float a;
//  for (short day=1;day<=31;day++){
//    a = dailyPrecNrd.getDailyPrec(day,0,5,8);
//    printf("%d %f\n",day,a);
//  }
//}
