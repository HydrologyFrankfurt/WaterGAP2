
/***********************************************************************
*
* see former changes at file gridio.cpp.versioninfos.txt
*
***********************************************************************/
#include "gridio.h"

void error(const char *string1, const char *string2, const char *string3)
{
	cerr << "ERROR: " << string1 << ' ' << string2 << ' ' << string3 << endl;
	exit(1);
}

gridioClass::gridioClass(void)
{
	// detect endian type of the system:
	// write two one-byte numbers (unsigned char) into memory
	// and use a two-byte pointer (unsigned short) to detect
	// which byte order is used
	unsigned char c[2] = { 0, 255 };
	unsigned short *p;

	p = (unsigned short *) &c[0];
	unsigned short i;

	i = *p;
	if (65280 == i)
		systemEndianType = 1;	// x86
	if (255 == i)
		systemEndianType = 2;	// sparc

	// no assumtion about endian type of files has been made: 
	byteSwapFlag = false;
}

void gridioClass::setFileEndianType(short fileEndianType)
{
	if (fileEndianType == systemEndianType)
		byteSwapFlag = false;
	else
		byteSwapFlag = true;
}

short gridioClass::getEndianType()
{
	return systemEndianType;
}

void gridioClass::memSwap(char *array,
						  const unsigned short dataTypeSize, const unsigned int arrayByteSize)
{
	// swap bytes within a given array 
	// (needed to reformat between Intel and SUN
	// binary-data format).
	// Example:
	// short int tmpArray[3]={1,2,3};
	// memSwap((char*)tmpArray,sizeof(short int),
	//            sizeof(tmpArray));

	unsigned int sourcePosition, destPosition;

	//----- don't swap bytes
	if (dataTypeSize == 1)
		return;

	char tmpByteValue;

	for (unsigned int i = 0; i < arrayByteSize; i += dataTypeSize) {
		for (unsigned short j = 0; j < dataTypeSize / 2; j++) {
			sourcePosition = i + j;
			destPosition = (unsigned int) i + dataTypeSize - j - 1;
			tmpByteValue = array[sourcePosition];
			array[sourcePosition] = array[destPosition];
			array[destPosition] = tmpByteValue;
		}
	}
}

/********** READ FILES **********/
void gridioClass::readUnfFile(char *file, int n_values, double * grid)
{
	float *temp_grid=new float[n_values];

	readUnfFile(file, n_values, temp_grid);

	//Convert floats to double
	for (long i=0; i<n_values; ++i) {
		grid[i]=static_cast<double>(temp_grid[i]);
	}

	delete[] temp_grid; temp_grid = 0;
}

void gridioClass::writeUnfFile(char *file, int n_values, double * grid)
{
	float *temp_grid=new float[n_values];

	//Convert floats to double
	for (long i=0; i<n_values; ++i) {
		temp_grid[i]=static_cast<float>(grid[i]);
	}

	writeUnfFile(file, n_values, temp_grid);


	delete[] temp_grid; temp_grid = 0;
}
