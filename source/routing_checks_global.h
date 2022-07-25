#ifndef WATERGAP_CHECKS_GLOBAL_H
#define WATERGAP_CHECKS_GLOBAL_H
// Check Felix Portmann 2015 for invalid data
// START Definitions
extern short chk_year;
extern short chk_month;
extern short chk_day_in_month;
extern long chk_n_1;
extern short chk_year_1;
extern short chk_month_1;
extern short chk_day_in_month_1;
extern long chk_n_2;
extern short chk_year_2;
extern short chk_month_2;
extern short chk_day_in_month_2;
extern long chk_n_3;
extern short chk_year_3;
extern short chk_month_3;
extern short chk_day_in_month_3;

extern void check_value_short_LT0(short value);
extern void check_value_short_NaN(short value);
extern void check_value_double_LT0(double value);
extern void check_value_double_NaN(double value);

extern void check_array_short_LT0_YYMMDD(short* array, long n_cells, short year, short month, short day);
extern void check_array_short_NaN_YYMMDD(short* array, long n_cells, short year, short month, short day);
extern void check_array_double_LT0_YYMMDD(double* array, long n_cells, short year, short month, short day);
extern void check_array_double_NaN_YYMMDD(double* array, long n_cells, short year, short month, short day);

extern void check_array_double_LT0(double* array, long n_cells);
extern void check_array_double_NaN(double* array, long n_cells);

extern void check_array2Dng12_double_LT0(double array[ng][12]);
extern void check_array2Dng12_double_NaN(double array[ng][12]);

extern void check_array2Dng365_double_LT0(double array[ng][365]);
extern void check_array2Dng365_double_NaN(double array[ng][365]);

extern void list_array2Dng365_double_LT0__strgmsg(double array[ng][365], string const &s);
extern void list_array2Dng12_double_LT0__strgmsg(double array[ng][12], string const &s);

extern void list_array2Dng365_double_NaN__strgmsg(double array[ng][365], string const &s);
extern void list_array2Dng12_double_NaN__strgmsg(double array[ng][12], string const &s);

//
// Check & throw an exception
//
void routingClass::check_LT0__365TWS() {
    cout << "check_LT0__365TWS() started" << endl;
    if ((( 5 == options.grid_store) || (6 == options.grid_store)) && ((options.outTotalWaterInStoragesDaily_km3) || (options.outTotalWaterInStoragesDaily_mm)||(((options.scoutTWSkm3)||(options.scoutTWSmm))&&(2==options.day_store)))){
        check_array2Dng365_double_LT0(G_daily365TotalWaterInStorages_km3);
        if (options.outTotalWaterInStoragesDaily_mm || options.scoutTWSmm) {
            check_array2Dng365_double_LT0(G_daily365TotalWaterInStorages_mm);
        }
    }
    //endif options
}
// end check_LT0__365TWS

void routingClass::check_NaN__365TWS() {
    cout << "check_NaN__365TWS() started" << endl;
    if ((( 5 == options.grid_store) || (6 == options.grid_store)) && ((options.outTotalWaterInStoragesDaily_km3) || (options.outTotalWaterInStoragesDaily_mm)||(((options.scoutTWSkm3)||(options.scoutTWSmm))&&(2==options.day_store)))){
        check_array2Dng365_double_NaN(G_daily365TotalWaterInStorages_km3);
        if (options.outTotalWaterInStoragesDaily_mm || options.scoutTWSmm) {
            check_array2Dng365_double_NaN(G_daily365TotalWaterInStorages_mm);
        }
    }
    //endif options
}
// end check_NaN__365TWS


void routingClass::check_LT0__365AET() {
    cout << "check_LT0__365AET() started" << endl;
    check_array2Dng365_double_LT0(G_daily365CellAET);
}
// end check_LT0__365AET

void routingClass::check_NaN__365AET() {
    cout << "check_NaN__365AET() started" << endl;
    check_array2Dng365_double_NaN(G_daily365CellAET);
}
// end check_NaN__365AET


void routingClass::check_LT0__annualInit() {

    cout << "check_LT0__annualInit() started" << endl;

        check_array_double_LT0(G_UnsatisfiedUsePrevYear, ng);
        check_array_double_LT0(G_totalUnsatisfiedUse, ng);
            // may be negative = exclude
            // G_satisfiedUse[n];
            // G_actualUse[n];
            // G_potCellRunoff[n];

    if ((2 == options.grid_store) || (4 == options.grid_store) || (6 == options.grid_store)) {

        check_array2Dng12_double_LT0(G_monthlyConsistentPrecip);
        check_array2Dng12_double_LT0(G_monthlyRiverAvail);
            // may be negative = exclude
            // G_monthlyCellSurfaceRunoff[n][m];
        check_array2Dng12_double_LT0(G_monthlyOpenWaterEvap);
        check_array2Dng12_double_LT0(G_monthlyVelocity);
            // may be negative = exclude
            // G_monthlySurfStor[n][m];
            // G_monthlySurfStor_mm[n][m];
        check_array2Dng12_double_LT0(G_monthlyMinRiverAvail);
        check_array2Dng12_double_LT0(G_monthlyMaxRiverAvail);
        check_array2Dng12_double_LT0(G_monthlyGwrSwb);
        check_array2Dng12_double_LT0(G_monthlyLandAreaFrac);
        check_array2Dng12_double_LT0(G_monthlyFswb);
        check_array2Dng12_double_LT0(G_monthlyGwRunoff);
            // may be negative = exclude
            // G_monthlyPotCellRunoff[n][m];
            // G_monthlyPotCellRunoffDeficit[n][m];
            // G_monthlyCellRunoff[n][m];
            // G_monthlyCellAET[n][m];

        if (options.outRiverPET == 1) {
            check_array2Dng12_double_LT0(G_monthlyRiverAreaFrac);
            check_array2Dng12_double_LT0(G_monthlyRiverPET);
        }

        if (options.resOpt == 1) {
            check_array2Dng12_double_LT0(G_monthlyResStorage);
            check_array2Dng12_double_LT0(G_monthlyResStorageRatio);
            check_array2Dng12_double_LT0(G_monthlyResInflow);
            check_array2Dng12_double_LT0(G_monthlyResOutflow);
        }

        if (options.outRiverInUpstream) {
            check_array2Dng12_double_LT0(G_monthlyRiverInUpstream);
        }

        check_array2Dng12_double_LT0(G_monthlyRiverStorage);
            // may be negative = exclude
            //G_monthlyGwStorage[n][m];
            //G_monthlyLocLakeStorage[n][m];
            //G_monthlyGloLakeStorage[n][m];
        check_array2Dng12_double_LT0(G_monthlyLocWetlStorage);
        check_array2Dng12_double_LT0(G_monthlyGloWetlStorage);

            // may be negative = exclude
            // G_monthlySatisfiedUse[n][m];
            // G_monthlyActualUse[n][m];

    }
    // endif options

}
// end check_LT0__annualInit


void routingClass::check_NaN__annualInit() {

    cout << "check_NaN__annualInit() started" << endl;

        check_array_double_NaN(G_UnsatisfiedUsePrevYear, ng);
        check_array_double_NaN(G_totalUnsatisfiedUse, ng);
            // may be negative 3x
            check_array_double_NaN(G_satisfiedUse, ng);
            check_array_double_NaN(G_actualUse, ng);
            check_array_double_NaN(G_potCellRunoff, ng);


    if ((2 == options.grid_store) || (4 == options.grid_store) || (6 == options.grid_store)) {

        check_array2Dng12_double_NaN(G_monthlyConsistentPrecip);
        check_array2Dng12_double_NaN(G_monthlyRiverAvail);
            // may be negative 1x
            check_array2Dng12_double_NaN(G_monthlyCellSurfaceRunoff);
        check_array2Dng12_double_NaN(G_monthlyOpenWaterEvap);
        check_array2Dng12_double_NaN(G_monthlyVelocity);
            // may be negative 2x
            check_array2Dng12_double_NaN(G_monthlySurfStor);
            check_array2Dng12_double_NaN(G_monthlySurfStor_mm);
        check_array2Dng12_double_NaN(G_monthlyMinRiverAvail);
        check_array2Dng12_double_NaN(G_monthlyMaxRiverAvail);
        check_array2Dng12_double_NaN(G_monthlyGwrSwb);
        check_array2Dng12_double_NaN(G_monthlyLandAreaFrac);
        check_array2Dng12_double_NaN(G_monthlyFswb);
        check_array2Dng12_double_NaN(G_monthlyGwRunoff);
            // may be negative 4x
            check_array2Dng12_double_NaN(G_monthlyPotCellRunoff);
            check_array2Dng12_double_NaN(G_monthlyPotCellRunoffDeficit);
            check_array2Dng12_double_NaN(G_monthlyCellRunoff);
            check_array2Dng12_double_NaN(G_monthlyCellAET);

        if (options.outRiverPET == 1) {
            check_array2Dng12_double_NaN(G_monthlyRiverAreaFrac);
            check_array2Dng12_double_NaN(G_monthlyRiverPET);
        }

        if (options.resOpt == 1) {
            check_array2Dng12_double_NaN(G_monthlyResStorage);
            check_array2Dng12_double_NaN(G_monthlyResStorageRatio);
            check_array2Dng12_double_NaN(G_monthlyResInflow);
            check_array2Dng12_double_NaN(G_monthlyResOutflow);
        }

        if (options.outRiverInUpstream) {
            check_array2Dng12_double_NaN(G_monthlyRiverInUpstream);
        }

        check_array2Dng12_double_NaN(G_monthlyRiverStorage);
            // may be negative 3x
            check_array2Dng12_double_NaN(G_monthlyGwStorage);
            check_array2Dng12_double_NaN(G_monthlyLocLakeStorage);
            check_array2Dng12_double_NaN(G_monthlyGloLakeStorage);

        check_array2Dng12_double_NaN(G_monthlyLocWetlStorage);
        check_array2Dng12_double_NaN(G_monthlyGloWetlStorage);
        check_array2Dng12_double_NaN(G_monthlySatisfiedUse);
        check_array2Dng12_double_NaN(G_monthlyActualUse);

    }
    // endif options

}
// end check_NaN__annualInit


//
// Check & list interesting values
//
void routingClass::list_LT0__365TWS() {
    cout << "list_LT0__365TWS() started" << endl;
    if ((( 5 == options.grid_store) || (6 == options.grid_store)) && ((options.outTotalWaterInStoragesDaily_km3) || (options.outTotalWaterInStoragesDaily_mm)||(((options.scoutTWSkm3)||(options.scoutTWSmm))&&(2==options.day_store)))){
        string s_temp = string("G_daily365TotalWaterInStorages_km3");
        list_array2Dng365_double_LT0__strgmsg(G_daily365TotalWaterInStorages_km3, s_temp);
        if (options.outTotalWaterInStoragesDaily_mm || options.scoutTWSmm) {
            s_temp = string("G_daily365TotalWaterInStorages_mm");
            list_array2Dng365_double_LT0__strgmsg(G_daily365TotalWaterInStorages_mm, s_temp);
        }
    }
    //endif options
}
// end list_LT0__365TWS

void routingClass::list_NaN__365TWS() {
    cout << "list_NaN__365TWS() started" << endl;
    if ((( 5 == options.grid_store) || (6 == options.grid_store)) && ((options.outTotalWaterInStoragesDaily_km3) || (options.outTotalWaterInStoragesDaily_mm)||(((options.scoutTWSkm3)||(options.scoutTWSmm))&&(2==options.day_store)))){
        string s_temp = string("G_daily365TotalWaterInStorages_km3");
        list_array2Dng365_double_NaN__strgmsg(G_daily365TotalWaterInStorages_km3, s_temp);
        if (options.outTotalWaterInStoragesDaily_mm || options.scoutTWSmm) {
            s_temp = string("G_daily365TotalWaterInStorages_mm");
            list_array2Dng365_double_NaN__strgmsg(G_daily365TotalWaterInStorages_mm, s_temp);
        }
    }
    //endif options
}
// end list_NaN__365TWS

void routingClass::list_LT0__365AET() {
    cout << "list_LT0__365AET() started" << endl;
    string s_temp = string("G_daily365CellAET");
    list_array2Dng365_double_LT0__strgmsg(G_daily365CellAET, s_temp);
}
// end list_LT0__365AET

void routingClass::list_NaN__365AET() {
    cout << "list_NaN__365AET() started" << endl;
    string s_temp = string("G_daily365CellAET");
    list_array2Dng365_double_NaN__strgmsg(G_daily365CellAET, s_temp);
}
// end list_NaN__365AET


// List negative values of monthly output grids (selected variables from annualInit)
void routingClass::list_LT0__annualInit() {

    cout << "list_LT0__annualInit() started (selected variables)" << endl;

    string s_temp;

    if ((2 == options.grid_store) || (4 == options.grid_store) || (6 == options.grid_store)) {

        // discharge
        // (only positive values allowed)
        s_temp = string("G_monthlyRiverAvail");
        list_array2Dng12_double_LT0__strgmsg(G_monthlyRiverAvail, s_temp);

        // monthly cell AET
        // (a small number of cells with negative values expected)
        s_temp = string("G_monthlyCellAET");
        list_array2Dng12_double_LT0__strgmsg(G_monthlyCellAET, s_temp);

            // may be negative = exclude
            // s_temp = string("G_monthlyPotCellRunoff");
            // list_array2Dng12_double_LT0__strgmsg(G_monthlyPotCellRunoff, s_temp);
            // s_temp = string("G_monthlyPotCellRunoffDeficit");
            // list_array2Dng12_double_LT0__strgmsg(G_monthlyPotCellRunoffDeficit, s_temp);
            // s_temp = string("G_monthlyCellRunoff");
            // list_array2Dng12_double_LT0__strgmsg(G_monthlyCellRunoff, s_temp);

    }
    // endif options
}
// end list_LT0__annualInit

// List NaN values of monthly output grids (selected variables from annualInit)
void routingClass::list_NaN__annualInit() {

    cout << "list_NaN__annualInit() started (selected variables)" << endl;

    string s_temp;

    if ((2 == options.grid_store) || (4 == options.grid_store) || (6 == options.grid_store)) {

        // discharge
        s_temp = string("G_monthlyRiverAvail");
        list_array2Dng12_double_NaN__strgmsg(G_monthlyRiverAvail, s_temp);

        // monthly cell AET
        s_temp = string("G_monthlyCellAET");
        list_array2Dng12_double_NaN__strgmsg(G_monthlyCellAET, s_temp);

        s_temp = string("G_monthlyPotCellRunoff");
        list_array2Dng12_double_NaN__strgmsg(G_monthlyPotCellRunoff, s_temp);

        s_temp = string("G_monthlyPotCellRunoffDeficit");
        list_array2Dng12_double_NaN__strgmsg(G_monthlyPotCellRunoffDeficit, s_temp);

        s_temp = string("G_monthlyCellRunoff");
        list_array2Dng12_double_NaN__strgmsg(G_monthlyCellRunoff, s_temp);

    }
    // endif options
}
// end list_NaN__annualInit

// Check Felix Portmann 2015 for invalid data
// END Definitions


#endif // WATERGAP_CHECKS_GLOBAL_H
