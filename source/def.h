#if !defined (_def_h_)
#define _def_h_

			/* must fit into int, for watergap.cpp and routing.cpp */
// for Standard-landmask (e.g. H08) (landmaskOpt = 0) comment in & dont forget do do changes in routing.h:
//  #define ng 66896		/* number of grid cells */
//  #define ng_land 66663  /* number of land cells */ HMS 2015-05-04 not anymore separated

// for Watch-CRU-landmask (e.g. WFD) (landmaskOpt = 1) comment in & dont forget to do changes in routing.h:
  #define ng 67420		/* number of grid cells in WATCH */
//  #define ng_land 67327  /* number of land cells in WATCH */

//#define ng     66636             // land and freshwater (lake) cells
//#define ng_land 66559             // land cells only

#define ng_climate 70412//number of grid cells in climate land mask // Reminder: Possible extension for future climate data sets other than WFD_bc_WFDEI)

//#define ng21 59831
//#define nlct 20		/* number of land cover types: GLCT */
#define nlct 18		/* number of land cover types: GLCC (later: CLC & GLCT) */
//#define nreg 18		/* number of image regions */

#define reservoir_dsc 5  // number of reservoir downstream cells for new reservoir algorithm; only for irrigation reservoirs
#define wateruse_dsc  5  // number of downstream cells for allocation of consumptive uses

#define ng_glolakcells 0 //GFZ:GLOLAK 
//#define ng_glolakcells 1170 /*no of lakecells in wg_big_lakes_area_frac.txt*/
			    /*for no allocation of glolakestorage define to zero!*/


#endif
