/**********************************************************************
* Modification (MH20160601N001):
* The arguments of the program was not working properly. During this
* modification the bug is fixed temporarily (in a rather orthodox manner).
*
* This modification is required for the parallel implementation of model
* calibration using Borg-MOEA 1.8.
*
* Modified by: H.M. Mehedi Hasan, Dated on: June 01, 2016
***********************************************************************/
/**********************************************************************
* Modifications by: Felix T. Portmann (FP), Institute of Physical Geography, Goethe University Frankfurt am Main
*
* FP20160915N001:
* Arguments of program from update MH20160601N001:
* Now in error message also arguments of single cells output are included
* (watergap.cpp, options.cpp)
*
* FP20160915N002: scoutSurfaceRunoff -> options.scoutCellSRunoff
* Correction of initialization of G_daily365CellSurfaceRunoff: options.scoutCellSRunoff instead of erroneous scoutSurfaceRunoff
* (routing.cpp)
*
* FP20160926N001: changed for daily output .31 and .365 arrays: replaced "+=" by consistent "="
* for several variables e.g. river availability, transported volume, river velocity
* (routing.cpp)
*
* FP20161005N001: double G_landAreaFracPrevTimestep__at__n = routing.G_landAreaFracPrevTimestep[n];
* (daily.cpp)
*
* FP20161010N001 Replacement: geo.G_contfreq[n] instead (geo.G_landfreq[n] + geo.G_fwaterfreq[n])
* Class geo.writeContFreq for writing resulting grid to file GCONTFREQ.UNF0
* (geo.cpp/h, daily.cpp, routing.cpp, watergap.cpp)
*
* FP20161010N002: remove internal selection "if (calibrun == 0)"; now selection before calling "if (geo.calibrun == 0)"
* (geo.cpp)
*
* FP20161011N001: Introduce monthly rainfall output (as monthly snowfall output is also written)
* (daily.cpp)
*
* FP20161012N001 Enabling treatment of maximum slope class in G_SLOPE_CLASS.UNF1 (currently in file max. 69)
* FP20161012N002 Corrected error treatment (for values outside class limits) for G_TEXTURE.UNF1
* (gw_frac.cpp)
*
* FP20161018N001 Evaluation of arguments since first argument
* FP20161018N002 Reservoir operation start years
* FP20161018N003 Enable reading input data from climate land mask
* FP20161018N004 Net use / net abstraction for reservoir algorithm: new filename
* (climate.cpp, def.h, geo.cpp, option.h/cpp, routing.h/cpp, watergap.cpp)
***********************************************************************/

/***********************************************************************
 *Since 2.1f (around the year 2006) the development of WaterGAP (version with
 *spatial resolution of 0.5°) at CESR and IPG resulted in different model
 *versions and features. At the WaterGAP coordination meeting in 2011 it was
 *decided to create a harmonized model version of WaterGAP 2 integrating all
 *current developments. This model version (WaterGAP 2.2 refvers) is the basis
 *for further code development. The old headers (version control infos) were
 *saved in a textfile and replaced by comments on changes since 2.1f. Names of
 *persons who contributed to the code changes were added in order to get the
 *right person in touch for questions. Where possible, references were added.
 *Within the most of the .cpp headers, comments were added (where the changes
 *occur). Here, a short summary is done. More detailled documentation is in
 *preparation.
 *
 * WaterGAP 2.2 refvers has a lot of new features, including: handling of daily
 * climate input (as days per month and as days per year); variable flow
 * velocity including bankfull flow; water use can be taken out of surface
 * water or ground water sources (water_use.cpp is the basis for the
 * preprocessing tool GWSWUSE); a new land cover classification scheme (IGBP)
 * and data basis (combination of GLCT & Corine) as well as revised land use
 * attributes are used; arid/humid definition can be calculated at Köppen
 * climates; several new radiation (e.g. land use class based emissivity value
 * for calculating outgoing long wave radiation) and PET options (e.g. crop
 * coefficients) were included; distribution of permafrost can be calculated;
 * the growing period is based on precipitation and not on PET anymore;
 * precipitation correction is done after Adam & Lettenmaier; some user friedly
 * things were added (e.g. MODEL_SETTINGS_OUT.DAT), the original sensitivity
 * analysis was deleted, several input files were updated and many output files
 * were added.
 *
 *WaterGAP 2.2 refvers was calibrated within CESR at 1xxx GRDC stations using
 *WATCH Forcing Data as climate input; IPG used 1323 GRDC stations and CRU 3.1
 *climate input (with GPCC v6 precip). The calibration period was set to
 *1971-2000 (if discharge data were available at this time span).
 *
 * Stephanie Eisner & Hannes Müller Schmied, Kassel and Frankfurt in June 2012
 *
 * see former changes at file watergap.cpp.versioninfos.txt
 *
 ***********************************************************************/
// WaterGAP
// Water - Global Assessment and Prognosis
//
// original developed by Frank Kaspar, Kassel, 1997-2000
// Center for Environmental Systems Reserach
// University of Kassel
// kaspar@usf.uni-kassel.de

// since that the following people were involved in code development and
//provided model input files (alphabetical order): Adam, Linda; Aus der Beek,
//Tim; Döll, Petra; Eisner, Stephanie; Fiedler (now Wilkinson), Kristina; Flörke,
//Martina; Hoffmann-Dobrev, Heike; Hunger, Martin; Müller Schmied, Hannes = HMS;
//Portmann, Felix = FP; Schneider, Christof; Schulze (now Verzano), Kerstin; Siebert,
//Stefan; Teichner (now Kynast), Ellen; Voß, Frank; Weiß, Martina; Zhang, Jing (tbc.)

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>

#ifdef _OPENMP
#include <omp.h>
#endif
#include <sys/time.h>
using namespace std;

// OE: New include files
#include "common.h"
#include "initializeWGHM.h"
#include "integrateWGHM.h"

int main(int argc, char* argv[]) {
    try{
        // AEMS: Evaluate command line:
        std::string fileName;
        switch(argc)
        {
            case 2:
                fileName = argv[1];
                break;
            default:
                std::cerr<<"USAGE: watergap <configFile>"<<std::endl;
                std::cerr<<"configFile: path to config file"<<std::endl;
                return -1;
        }

        // Variable declaration
        std::string progName = "OL"; // standard configuration for WaterGap without CDA
        WghmStateFile *wghmState;
        WghmStateFile *wghmMean;
        calibParamClass *calParam;
        AdditionalOutputInputFile *additionalOutIn;
        SnowInElevationFile *snow_in_elevation;
        ConfigFile *configFile;

        std::string path_wghmMean;
        long year;long month;
        long total_steps; long step;
        const char * fullname = fileName.c_str();
        const char * prog = progName.c_str();
        const char * path = path_wghmMean.c_str();

        initialize_wghm(fullname, wghmState,calParam, additionalOutIn, snow_in_elevation, &year, &month, prog, path, wghmMean);

        integrate_wghm(fullname,configFile,wghmState,calParam, additionalOutIn,snow_in_elevation,&step,&total_steps,&year,&month,prog);

    }
    catch(std::exception &e)
    {
        throw(Exception(std::string("In watergap::main():\n")+e.what()));
    }

}
