
/***********************************************************************
*
* see former changes at file gridio.h.versioninfos.txt
*
***********************************************************************/

#include <iostream>
#include <cstdlib>
#include <cstdio>
using namespace std;

//#define debug

#if !defined (_gridio_h_)
#define _gridio_h_

void error(const char *string1, const char *string2 = "", const char *string3 = "");

// template functions have to be directly defined in this .h-file
template < class T > void read_grid(char *file, int n_values, T * grid)
{
	int x;
	FILE *file_ptr;

	file_ptr = fopen(file, "rb");
	if (file_ptr != NULL) {
		x = fread(grid, sizeof(T), n_values, file_ptr);
#ifdef debug
		cout << x << " objects successfully read from file " << file << endl;
#endif
		if (x != n_values) 
			error("Unable to read the required number of objects from file: ", file);
		fclose(file_ptr);
	} else
		error("Unable to open", file, "for reading!");
}

template < class T > void write_grid(char *file, int n_values, T * grid)
{
	int x;
	FILE *file_ptr;

	file_ptr = fopen(file, "wb");
	if (file_ptr != NULL) {
		x = fwrite(grid, sizeof(T), n_values, file_ptr);
#ifdef debug
		cout << x << " objects successfully written to file " << file << endl;
#endif
		if (x != n_values)
			error("Unable to write the required number of objects to file: ", file);
		fclose(file_ptr);
	} else
		error("Unable to open", file, "for writing!");
}

class gridioClass {
  private:
	bool byteSwapFlag;
	short systemEndianType;

  public:
	 gridioClass();
	void setFileEndianType(short endianType);
	short getEndianType();
	bool getSwapFlag();

	void memSwap(char *array, const unsigned short dataTypeSize, const unsigned int arrayByteSize);

	 template < class T > void readUnfFile(char *file, int n_values, T * grid);
	
	void readUnfFile(char *file, int n_values, double * grid);

	 template < class T > void writeUnfFile(char *file, int n_values, T * grid);
	
	void writeUnfFile(char *file, int n_values, double * grid);

};


template < class T > void gridioClass::readUnfFile(char *file, int n_values, T * grid)
{
	int x;
	FILE *file_ptr;

	file_ptr = fopen(file, "rb");
	if (file_ptr != NULL) {
		x = fread(grid, sizeof(T), n_values, file_ptr);
#ifdef debug
		cout << x << " objects successfully read from file " << file << endl;
#endif
		if (x != n_values) {
			cout << "124 x: " << x << "n_values: " << n_values << endl;
			error("Unable to read the required number of objects from file: ", file);
		}
		fclose(file_ptr);
	} else
		error("Unable to open", file, "for reading!");

	if (byteSwapFlag)
		memSwap((char *) grid, sizeof(T), n_values * sizeof(T));
}

template < class T > void gridioClass::writeUnfFile(char *file, int n_values, T * grid)
{
	int x;
	FILE *file_ptr;

	file_ptr = fopen(file, "wb");
	if (file_ptr != NULL) {
		if (byteSwapFlag)
			memSwap((char *) grid, sizeof(T), n_values * sizeof(T));
		x = fwrite(grid, sizeof(T), n_values, file_ptr);
		if (byteSwapFlag)
			memSwap((char *) grid, sizeof(T), n_values * sizeof(T));
		// memSwap has to be done before and after writing the file
		// so that the array contains the same contents again

#ifdef debug
		cout << x << " objects successfully written to file " << file << endl;
#endif
		if (x != n_values)
			error("Unable to write the required number of objects to file: ", file);
		fclose(file_ptr);
	} else
		error("Unable to open", file, "for writing!");
}

inline bool gridioClass::getSwapFlag()
{
	return byteSwapFlag;
}

#endif
