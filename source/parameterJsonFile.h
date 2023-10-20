#include "common.h"
#include "timestring.h"
#include "calib_param.h"
#include "def.h"//OE included
#include "matrix.h"//OE included
using namespace json11;

/*
This class is based on "wg_param2jsonfile.cpp" (Frankfurt University). It saves a WaterGAP calibration parameter set as JSON file. 
By Kerstin Schulze (02/2019) - schulze@geod.uni-bonn.de
OE,KS 12.03.2020: allow ALL calPar being spatially variable
*/

class parameterJsonFile
{ // Calibration parameters
    public:
    
    parameterJsonFile();
    parameterJsonFile(std::string filenameJson);
    void parameterJsonFile_cda(Matrix &, int *groupmatrixindex,int nr_cda_unit, int ids);// OE included
    void setCFAandCFSfromJson(std::string filename);
    void setCFAandCFStoConstantValue(double value);
    void setGammaHBVfromJson(std::string filename);
    void save(std::string filenameInput, std::string filenameInputArcID);
    void save_cda_txt(std::string,double value, Matrix &); // OE included
    std::vector<string> getNamesCalPar();// OE included


    // Getter
    std::vector<double> get_gammaHBV(){return _gammaHBV;};
    std::vector<double> get_CFA(){return _CFA;};
    std::vector<double> get_CFS(){return _CFS;};
    std::vector<double> get_rootDepthMult(){return _rootDepthMult;};
    std::vector<double> get_riverRoughnessCoeffMult(){return _riverRoughnessCoeffMult;};
    std::vector<double> get_lakeDepth(){return _lakeDepth;};
    std::vector<double> get_wetlandDepth(){return _wetlandDepth;};
    std::vector<double> get_swOutflowCoeff(){return _swOutflowCoeff;};
    std::vector<double> get_evapoRedFactExpMult(){return _evapoRedFactExpMult;};
    std::vector<double> get_netRadiationMult(){return _netRadiationMult;};
    std::vector<double> get_PTcoeffHumid(){return _PTcoeffHumid;};
    std::vector<double> get_PTcoeffArid(){return _PTcoeffArid;};
    std::vector<double> get_maxDailyPET(){return _maxDailyPET;};
    std::vector<double> get_mcwh(){return _mcwh;};
    std::vector<double> get_laiMult(){return _rootDepthMult;};
    std::vector<double> get_snowFreezeTemp(){return _snowFreezeTemp;};
    std::vector<double> get_snowMeltTemp(){return _snowMeltTemp;};
    std::vector<double> get_degreeDayFactorMult(){return _degreeDayFactorMult;};
    std::vector<double> get_temperatureGradient(){return _temperatureGradient;};
    std::vector<double> get_gwFactorMult(){return _gwFactorMult;};
    std::vector<double> get_rgMaxMult(){return _rgMaxMult;};
    std::vector<double> get_pCritAridGW(){return _pCritAridGW;};
    std::vector<double> get_gwOutflowCoeff(){return _gwOutflowCoeff;};
    std::vector<double> get_netAbstractionSWMult(){return _netAbstractionSWMult;};
    std::vector<double> get_netAbstractionGWMult(){return _netAbstractionGWMult;};
    std::vector<double> get_precipMult(){return _precipMult;};

    // Setter
    void set_gammaHBV(std::vector<double> x){_gammaHBV=x;};
    void set_CFA(std::vector<double> x){_CFA=x;};
    void set_CFS(std::vector<double> x){_CFS=x;};
    void set_rootDepthMult(std::vector<double> x){_rootDepthMult=x;};
    void set_riverRoughnessCoeffMult(std::vector<double> x){_riverRoughnessCoeffMult=x;};
    void set_lakeDepth(std::vector<double> x){_lakeDepth=x;};
    void set_wetlandDepth(std::vector<double> x){_wetlandDepth=x;};
    void set_swOutflowCoeff(std::vector<double> x){_swOutflowCoeff=x;};
    void set_evapoRedFactExpMult(std::vector<double> x){_evapoRedFactExpMult=x;};
    void set_netRadiationMult(std::vector<double> x){_netRadiationMult=x;};
    void set_PTcoeffHumid(std::vector<double> x){_PTcoeffHumid=x;};
    void set_PTcoeffArid(std::vector<double> x){_PTcoeffArid=x;};
    void set_maxDailyPET(std::vector<double> x){_maxDailyPET=x;};
    void set_mcwh(std::vector<double> x){_mcwh=x;};
    void set_laiMult(std::vector<double> x){_rootDepthMult=x;};
    void set_snowFreezeTemp(std::vector<double> x){_snowFreezeTemp=x;};
    void set_snowMeltTemp(std::vector<double> x){_snowMeltTemp=x;};
    void set_degreeDayFactorMult(std::vector<double> x){_degreeDayFactorMult=x;};
    void set_temperatureGradient(std::vector<double> x){_temperatureGradient=x;};
    void set_gwFactorMult(std::vector<double> x){_gwFactorMult=x;};
    void set_rgMaxMult(std::vector<double> x){_rgMaxMult=x;};
    void set_pCritAridGW(std::vector<double> x){_pCritAridGW=x;};
    void set_gwOutflowCoeff(std::vector<double> x){_gwOutflowCoeff=x;};
    void set_netAbstractionSWMult(std::vector<double> x){_netAbstractionSWMult=x;};
    void set_netAbstractionGWMult(std::vector<double> x){_netAbstractionGWMult=x;};
    void set_precipMult(std::vector<double> x){_precipMult=x;};

    private:
    std::vector<double> _gammaHBV;
    std::vector<double> _CFA;
    std::vector<double> _CFS;
    std::vector<double> _rootDepthMult;
    std::vector<double> _riverRoughnessCoeffMult;
    std::vector<double> _lakeDepth;
    std::vector<double> _wetlandDepth;
    std::vector<double> _swOutflowCoeff;
    std::vector<double> _evapoRedFactExpMult;
    std::vector<double> _netRadiationMult;
    std::vector<double> _PTcoeffHumid;
    std::vector<double> _PTcoeffArid;
    std::vector<double> _maxDailyPET;
    std::vector<double> _mcwh;
    std::vector<double> _laiMult;
    std::vector<double> _snowFreezeTemp;
    std::vector<double> _snowMeltTemp;
    std::vector<double> _degreeDayFactorMult;
    std::vector<double> _temperatureGradient;
    std::vector<double> _gwFactorMult;
    std::vector<double> _rgMaxMult;
    std::vector<double> _pCritAridGW;
    std::vector<double> _gwOutflowCoeff;
    std::vector<double> _netAbstractionSWMult;
    std::vector<double> _netAbstractionGWMult;
    std::vector<double> _precipMult;
};
