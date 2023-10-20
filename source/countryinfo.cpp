
/***********************************************************************
*
* nothing changed since 2.1f
*
***********************************************************************/
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <cctype>
#include <vector>
#include <algorithm>
#include <iostream>
#include "countryinfo.h"


void countryInfoClass::init(const char *input_dir, const unsigned short maxLengthOfCountryName)
{
	// in the default case only the first 15 characters 
	// of the country names are used.

	char filename[250];

	sprintf(filename, "%s/COUNTRY_NAMES.DAT", input_dir);
	ifstream countryFile(filename);

	if (!countryFile) {
		cerr << "Error while opening file " << filename << "." << endl;
		exit(1);
	} else {
		char *name;
		char string[250];
		short length, j, i = 0;

		// read commentlines indicated by # 
		while (!countryFile.eof() && countryFile.peek() == '#') {
			countryFile.getline(string, sizeof(string));
		}

		// set data type to 'decimal'
		// otherwise some compiler assume that numbers with leading 'zeros'
		// are octal 
		cin.setf(ios::dec, ios::basefield);

		while (countryFile >> j >> ws) {
			// read iso number and following whitespaces

			countryIsoNumber.push_back(j);

			// read country name.
			// country names might include space characters 
			// therefore 'getline' is used.  
			countryFile.getline(string, sizeof(string));
			length = strlen(string);

			// allocate memory and copy string with the desired length into it
			if (length >= maxLengthOfCountryName) {
				name = new char[maxLengthOfCountryName + 1];

				length = maxLengthOfCountryName;
				strncpy(name, string, length);
			} else {
				name = new char[length + 1];

				strncpy(name, string, length);
				name[length] = '\0';
			}
			countryName.push_back(name);
			i++;
		}
		countryFile.close();

		numberOfCountries = i;

		// avoid unused memory
		countryIsoNumber.resize(numberOfCountries);
		countryName.resize(numberOfCountries);
	}
}

short countryInfoClass::getNumberOfCountries(void)
{
	return numberOfCountries;
}

short countryInfoClass::getIsoNumber(const short n)
{
	if ((n >= 0) && (n <= numberOfCountries - 1))
		return countryIsoNumber[n];
	else {
		cerr << "Country number out of range: "
			<< n << " (0.." << numberOfCountries - 1 << ")." << endl;
		exit(-1);
	}
	
	//should never be reached
	return -1;	
}

short countryInfoClass::getArrayPosition(const short IsoNumber)
{
	vector < short >::iterator it
		= find(countryIsoNumber.begin(), countryIsoNumber.end(), IsoNumber);
	if (it != countryIsoNumber.end())
		return (it - countryIsoNumber.begin());
	else {
		cerr << "No country found with ISO number " << IsoNumber << endl;
		exit(-1);
	}
	
	//should never be reached
	return -1;
}

char *countryInfoClass::getCountryName(const short n)
{
	if ((n >= 0) && (n <= numberOfCountries - 1))
		return countryName[n];
	else {
		cerr << "Country number out of range: "
			<< n << "(0.." << numberOfCountries - 1 << ")." << endl;
		exit(-1);
	}

	//should never be reached
	//return '\0';		
}
