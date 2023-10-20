#if !defined (_clcl_h_)
#define _clcl_h_


#include <vector>
#include "def.h"
#include "grid.h"

using namespace std;

class clclClass {
	public:
		clclClass();
		void calcKoep();
		Grid<float> clcl_alpha; 	// alpha factor for PT-calculations
		Grid<short> cls; 			// koeppen climate classification
		void init();
		void cleanup();
		void alphaKoep();
};


#endif
