
/***********************************************************************
*
* nothing changed since 2.1f
*
***********************************************************************/
#if !defined (_countryinfo_h_)
#define _countryinfo_h_

#include <vector>

using namespace std;

class countryInfoClass {
	public:
		void init(const char *input_dir, const unsigned short maxLengthOfCountryName = 15);
		short getIsoNumber(const short n);
		short getArrayPosition(const short IsoNumber);
		short getNumberOfCountries(void);
		char *getCountryName(const short n);

	private:
		vector < short >countryIsoNumber;
		vector < char *>countryName;
		short numberOfCountries;
};
#endif
