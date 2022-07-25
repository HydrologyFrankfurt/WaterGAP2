class groundwaterFactorClass {
  public:
	groundwaterFactorClass();
	void createGrids(const char *input_dir, const char *output_dir);

        //signed char G_texture[ng_land]; //++++++++++ new for GW_RECHARGE (semi)arid
        signed char G_texture[ng]; //++++++++++ new for GW_RECHARGE (semi)arid

	void setRgmax(int cell, short value);
	short getRgmax(int cell);

	void setgwFactor(int cell, float value);
	float getgwFactor(int cell);
	void setpermaFactor(int cell, float value);
	float getpermaFactor(int cell);

  private:
        //short G_Rgmax[ng_land];
        short G_Rgmax[ng];
        //float G_gwFactor[ng_land];
        float G_gwFactor[ng];
        //float G_gwFactorCorr[ng_land]; //CR 2014-07-31: introduced in WG2.2b to reduce fg in the Mississippi Embayment Regional Aquifer from 0.8-0.9 to 0.1
        float G_gwFactorCorr[ng]; //CR 2014-07-31: introduced in WG2.2b to reduce fg in the Mississippi Embayment Regional Aquifer from 0.8-0.9 to 0.1
        //float G_permaFactor[ng_land];
        float G_permaFactor[ng];

};
