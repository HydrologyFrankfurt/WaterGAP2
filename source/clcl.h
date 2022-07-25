#if !defined (_clcl_h_)
#define _clcl_h_


#include <vector>
#include "def.h"
//#include "gridio.h"

using namespace std;

class clclClass {
	public:
	clclClass();
	void calcKoep();  
	float (*clcl_alpha); 	// alpha factor for PT-calculations
	short (*cls); 			// koeppen climate classification
	void init();
	void cleanup();
	void alphaKoep();
	
	private:

};


#endif
