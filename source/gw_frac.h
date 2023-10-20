#if !defined (_gw_frac_h_)
#define _gw_frac_h_

#include "grid.h"
#include <string>
#include "calib_param.h"

class groundwaterFactorClass {
  public:
    groundwaterFactorClass();
    void createGrids(const std::string input_dir, const std::string output_dir, calibParamClass &calParam);

    Grid<signed char> G_texture;

    void setRgmax(int cell, short value);
    short getRgmax(int cell);

    void setgwFactor(int cell, float value);
    float getgwFactor(int cell);
    void setpermaFactor(int cell, float value);
    float getpermaFactor(int cell);

  private:
    Grid<short> G_Rgmax;
    Grid<float> G_textureFactor;
    Grid<float> G_aquiferFactor;
    Grid<float> G_gwFactor;
    Grid<float> G_gwFactorCorr; // introduced in WG2.2b to reduce fg in the Mississippi Embayment Regional Aquifer from 0.8-0.9 to 0.1
    Grid<float> G_permaFactor;

};

#endif
