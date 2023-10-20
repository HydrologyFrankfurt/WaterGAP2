#include <iostream>
#include "option.h"
#include "glacierYear.h"

using namespace std;

extern optionClass options;

glacierYearClass::glacierYearClass(void)
{
}

void glacierYearClass::init()
{
}

void glacierYearClass::cleanup() 
{
}

void glacierYearClass::read_initial_glacier_area(short start_year)
{
    string glacierDir = string(options.glacier_dir);

    G_glacier_area.read(glacierDir + "/G_GLACIER_AREA_km2_" + to_string(start_year) + ".UNF0");
}

void glacierYearClass::read_glacier_data_daily_per_year(short actual_year)
{
    string glacierDir = string(options.glacier_dir);

    G_precipitation_on_glacier_d365.read(glacierDir + "/G_GLACIER_PRECIPITATION_km3_" + to_string(actual_year) + ".365.UNF0");
    G_glacier_mass_change_d365.read(glacierDir + "/G_GLACIER_MASS_CHANGE_km3_" + to_string(actual_year) + ".365.UNF0");
}

void glacierYearClass::read_glacier_data_yearly_per_year(short actual_year) 
{
    string glacierDir = string(options.glacier_dir);

    G_glacier_area.read(glacierDir + "/G_GLACIER_AREA_km2_" + to_string(actual_year) + ".UNF0");
}


glacierYearClass::~glacierYearClass(void)
{
}