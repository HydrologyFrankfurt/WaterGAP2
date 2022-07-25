#ifndef WATERGAP_CHECKS_GLOBAL_H
#define WATERGAP_CHECKS_GLOBAL_H
// Check Felix Portmann 2015: date (year/month/day_in_month) to check for invalid data
// global checks
short chk_year;
short chk_month;
short chk_day_in_month;

//in classes use:
//extern short chk_year;
//extern short chk_month;
//extern short chk_day_in_month;

// local checks
long chk_n_1;
short chk_year_1;
short chk_month_1;
short chk_day_in_month_1;
long chk_n_2;
short chk_year_2;
short chk_month_2;
short chk_day_in_month_2;
long chk_n_3;
short chk_year_3;
short chk_month_3;
short chk_day_in_month_3;

// Functions included by Felix Portmann 2015-06-03 to check for invalid data
/**
 * Exception when values less than zero are encountered
 */
struct ExceptionValue : public std::exception {
       std::string s;
       ExceptionValue(std::string ss) : s(ss) {}
       ~ExceptionValue() throw () {}
       const char* what() const throw() { return s.c_str();}
};

struct ExceptionValueLT0__val : public std::exception {
       std::string s;
       std::string blank_str = " ";
       std::string hyphen_str = "-";
       const double value = -1; // initialization only in standard c++11; compile in makefile with "CC = g++ -std=c++11"
       ExceptionValueLT0__val(std::string ss, double value) : s(ss), value(value) {}
       ~ExceptionValueLT0__val() throw () {}
       const char* what() const throw() { return (s + blank_str + string(" Value encountered: ") + to_string(value)).c_str();}
};

struct ExceptionValueNaN : public std::exception {
       std::string s;
       ExceptionValueNaN(std::string ss) : s(ss) {}
       ~ExceptionValueNaN() throw () {}
       const char* what() const throw() { return s.c_str();}
};

struct ExceptionValueLT0_array__n : public std::exception {
       std::string s;
       std::string blank_str = " ";
       std::string hyphen_str = "-";
       const long n = -1; // initialization only in standard c++11; compile in makefile with "CC = g++ -std=c++11"
       const short year = -1; // initialization only in standard c++11; compile in makefile with "CC = g++ -std=c++11"
       const short month = -1; // initialization only in standard c++11; compile in makefile with "CC = g++ -std=c++11"
       const short day = -1; // initialization only in standard c++11; compile in makefile with "CC = g++ -std=c++11"
//       ExceptionValueLT0(std::string ss) : s(ss) {}
       ExceptionValueLT0_array__n(std::string ss, long n) : s(ss), n(n) {}
       ~ExceptionValueLT0_array__n() throw () {}
       const char* what() const throw() { return (s + blank_str + string("array index: ") + to_string(n)).c_str();}
};

struct ExceptionValueLT0_array__n_date : public std::exception {
       std::string s;
       std::string blank_str = " ";
       std::string hyphen_str = "-";
       const long n = -1; // initialization only in standard c++11; compile in makefile with "CC = g++ -std=c++11"
       const short year = -1; // initialization only in standard c++11; compile in makefile with "CC = g++ -std=c++11"
       const short month = -1; // initialization only in standard c++11; compile in makefile with "CC = g++ -std=c++11"
       const short day = -1; // initialization only in standard c++11; compile in makefile with "CC = g++ -std=c++11"
//       ExceptionValueLT0(std::string ss) : s(ss) {}
       ExceptionValueLT0_array__n_date(std::string ss, long n, short year, short month, short day) : s(ss), n(n), year(year), month(month), day(day) {}
       ~ExceptionValueLT0_array__n_date() throw () {}
       const char* what() const throw() { return (s + blank_str + string("gcrc n: ") + to_string(n) + string(" YYYY-MM-DD: ") + to_string(year) + hyphen_str + to_string(month) + hyphen_str + to_string(day)).c_str();}
};

struct ExceptionValueLT0_array2Dng12__n_month : public std::exception {
       std::string s;
       std::string blank_str = " ";
       std::string hyphen_str = "-";
       const long n = -1; // initialization only in standard c++11; compile in makefile with "CC = g++ -std=c++11"
       const short month = -1; // initialization only in standard c++11; compile in makefile with "CC = g++ -std=c++11"
       ExceptionValueLT0_array2Dng12__n_month(std::string ss, long n, short month) : s(ss), n(n), month(month) {}
       ~ExceptionValueLT0_array2Dng12__n_month() throw () {}
       const char* what() const throw() { return (s + blank_str + string("array index: ") + to_string(n) + string(" month index: ") + to_string(month)).c_str();}
};

struct ExceptionValueLT0_array2Dng365__n_day : public std::exception {
       std::string s;
       std::string blank_str = " ";
       std::string hyphen_str = "-";
       const long n = -1; // initialization only in standard c++11; compile in makefile with "CC = g++ -std=c++11"
       const short day = -1; // initialization only in standard c++11; compile in makefile with "CC = g++ -std=c++11"
       ExceptionValueLT0_array2Dng365__n_day(std::string ss, long n, short day) : s(ss), n(n), day(day) {}
       ~ExceptionValueLT0_array2Dng365__n_day() throw () {}
       const char* what() const throw() { return (s + blank_str + string("array index: ") + to_string(n) + string(" index 365day: ") + to_string(day)).c_str();}
};



struct ExceptionValueNaN_array__n : public std::exception {
       std::string s;
       std::string blank_str = " ";
       std::string hyphen_str = "-";
       const long n = -1; // initialization only in standard c++11; compile in makefile with "CC = g++ -std=c++11"
       ExceptionValueNaN_array__n(std::string ss, long n) : s(ss), n(n) {}
       ~ExceptionValueNaN_array__n() throw () {}
       const char* what() const throw() { return (s + blank_str + string("array index: ") + to_string(n)).c_str();}
};

struct ExceptionValueNaN_array__n_date : public std::exception {
       std::string s;
       std::string blank_str = " ";
       std::string hyphen_str = "-";
       const long n = -1; // initialization only in standard c++11; compile in makefile with "CC = g++ -std=c++11"
       const short year = -1; // initialization only in standard c++11; compile in makefile with "CC = g++ -std=c++11"
       const short month = -1; // initialization only in standard c++11; compile in makefile with "CC = g++ -std=c++11"
       const short day = -1; // initialization only in standard c++11; compile in makefile with "CC = g++ -std=c++11"
       ExceptionValueNaN_array__n_date(std::string ss, long n, short year, short month, short day) : s(ss), n(n), year(year), month(month), day(day) {}
       ~ExceptionValueNaN_array__n_date() throw () {}
       const char* what() const throw() { return (s + blank_str + string("gcrc n: ") + to_string(n) + string(" YYYY-MM-DD: ") + to_string(year) + hyphen_str + to_string(month) + hyphen_str + to_string(day)).c_str();}
};

struct ExceptionValueNaN_array2Dng12__n_month : public std::exception {
       std::string s;
       std::string blank_str = " ";
       std::string hyphen_str = "-";
       const long n = -1; // initialization only in standard c++11; compile in makefile with "CC = g++ -std=c++11"
       const short month = -1; // initialization only in standard c++11; compile in makefile with "CC = g++ -std=c++11"
       ExceptionValueNaN_array2Dng12__n_month(std::string ss, long n, short month) : s(ss), n(n), month(month) {}
       ~ExceptionValueNaN_array2Dng12__n_month() throw () {}
       const char* what() const throw() { return (s + blank_str + string("array index: ") + to_string(n) + string(" month index: ") + to_string(month)).c_str();}
};

struct ExceptionValueNaN_array2Dng365__n_day : public std::exception {
       std::string s;
       std::string blank_str = " ";
       std::string hyphen_str = "-";
       const long n = -1; // initialization only in standard c++11; compile in makefile with "CC = g++ -std=c++11"
       const short day = -1; // initialization only in standard c++11; compile in makefile with "CC = g++ -std=c++11"
       ExceptionValueNaN_array2Dng365__n_day(std::string ss, long n, short day) : s(ss), n(n), day(day) {}
       ~ExceptionValueNaN_array2Dng365__n_day() throw () {}
       const char* what() const throw() { return (s + blank_str + string("array index: ") + to_string(n) + string(" index 365day: ") + to_string(day)).c_str();}
};
// Check functions



// (1) Scalar values
// (1.1) Less Than (LT) 0

void check_value_double_LT0(double value) {
    if (value < 0) {
            throw ExceptionValueLT0__val("check_value_double_LT0: value less than zero 0!", value);
    }
}

void check_value_float_LT0(float value) {
    if (value < 0) {
            throw ExceptionValueLT0__val("check_value_float_LT0: value less than zero 0!", value);
    }
}

void check_value_short_LT0(short value) {
    if (value < 0) {
            throw ExceptionValueLT0__val("check_value_short_LT0: value less than zero 0!", value);
    }
}


// (1.1) Not a Number (NaN)
void check_value_double_NaN(double value) {
    if (value != value) {
            throw ExceptionValue("check_value_double_NaN: value NaN encountered!");
    }
}

void check_value_float_NaN(float value) {
    if (value != value) {
            throw ExceptionValue("check_value_float_NaN: value NaN encountered!");
    }
}

void check_value_short_NaN(short value) {
    if (value != value) {
            throw ExceptionValue("check_value_short_NaN: value NaN encountered!");
    }
}



// (2) Checking arrays (1D)
// (2.1) Less Than (LT) 0
void check_array_double_LT0(double* array, long n_cells) {
    for (long i=0; i<n_cells; i++) {
        if (array[i] < 0) {
                throw ExceptionValueLT0_array__n("check_array_double_LT0: value less than zero 0 first at cell index (e.g. gcrc): ", i);
        }
    }
}

void check_array_float_LT0(float* array, long n_cells) {
    for (long i=0; i<n_cells; i++) {
        if (array[i] < 0) {
                throw ExceptionValueLT0_array__n("check_array_float_LT0: value less than zero 0 first at cell index (e.g. gcrc): ", i);
        }
    }
}

void check_array_short_LT0(short* array, long n_cells) {
    for (long i=0; i<n_cells; i++) {
        if (array[i] < 0) {
                throw ExceptionValueLT0_array__n("check_array_short_LT0: value less than zero 0 first at cell index (e.g. gcrc): ", i);
        }
    }
}


void check_array_double_LT0_YYMMDD(double* array, long n_cells, short year, short month, short day) {
    for (long i=0; i<n_cells; i++) {
        if (array[i] < 0) {
                throw ExceptionValueLT0_array__n_date("check_array_double_LT0: Array value LT 0 first at cell index (e.g. gcrc): ", i, year, month, day);
        }
    }
}

void check_array_float_LT0_YYMMDD(float* array, long n_cells, short year, short month, short day) {
    for (long i=0; i<n_cells; i++) {
        if (array[i] < 0) {
                throw ExceptionValueLT0_array__n_date("check_array_float_LT0: Array value LT 0 first at cell index (e.g. gcrc): ", i, year, month, day);
        }
    }
}

void check_array_short_LT0_YYMMDD(short* array, long n_cells, short year, short month, short day) {
    for (long i=0; i<n_cells; i++) {
        if (array[i] < 0) {
                throw ExceptionValueLT0_array__n_date("check_array_short_LT0: Array value LT 0 first at cell index (e.g. gcrc): ", i, year, month, day);
        }
    }
}


// (2.2) Not a Number (NaN)
void check_array_double_NaN(double* array, long n_cells) {
    for (long i=0; i<n_cells; i++) {
        if (array[i] != array[i]) {
                throw ExceptionValueNaN_array__n("check_array_double_NaN: value NaN first at cell index (e.g. gcrc): ", i);
        }
    }
}

void check_array_float_NaN(float* array, long n_cells) {
    for (long i=0; i<n_cells; i++) {
        if (array[i] != array[i]) {
                throw ExceptionValueNaN_array__n("check_array_float_NaN: value NaN first at cell index (e.g. gcrc): ", i);
        }
    }
}

void check_array_short_NaN(short* array, long n_cells) {
    for (long i=0; i<n_cells; i++) {
        if (array[i] != array[i]) {
                throw ExceptionValueNaN_array__n("check_array_short_NaN: value NaN first at cell index (e.g. gcrc): ", i);
        }
    }
}


void check_array_double_NaN_YYMMDD(double* array, long n_cells, short year, short month, short date) {
    for (long i=0; i<n_cells; i++) {
        if (array[i] != array[i]) {
                throw ExceptionValueNaN_array__n_date("check_array_double_NaN: Array value NaN first at cell index (e.g. gcrc): ", i, year, month, date);
        }
    }
}

void check_array_float_NaN_YYMMDD(float* array, long n_cells, short year, short month, short date) {
    for (long i=0; i<n_cells; i++) {
        if (array[i] != array[i]) {
                throw ExceptionValueNaN_array__n_date("check_array_float_NaN: Array value NaN first at cell index (e.g. gcrc): ", i, year, month, date);
        }
    }
}

void check_array_short_NaN_YYMMDD(short* array, long n_cells, short year, short month, short date) {
    for (long i=0; i<n_cells; i++) {
        if (array[i] != array[i]) {
                throw ExceptionValueNaN_array__n_date("check_array_short_NaN: Array value NaN first at cell index (e.g. gcrc): ", i, year, month, date);
        }
    }
}


// (3) Checking arrays 2D = ng*12 or ng*365
// 12 months
// (3.1.1) Less Than (LT) 0
void check_array2Dng12_double_LT0(double array[ng][12]) {
    for (long i=0; i<ng; i++) {
        for (long m=0; m<11; m++) {
            if (array[i][m] < 0) {
                throw ExceptionValueLT0_array2Dng12__n_month("check_array2Dng12_double_LT0: value less than zero 0 first at cell index (e.g. gcrc) & month (1-12): ", i, m+1);
            }
        }
    }
}
// (3.1.2) Not a Number (NaN)
void check_array2Dng12_double_NaN(double array[ng][12]) {
    for (long i=0; i<ng; i++) {
        for (long m=0; m<11; m++) {
            if (array[i][m] != array[i][m]) {
                throw ExceptionValueNaN_array2Dng12__n_month("check_array2Dng12_double_NaN: 2D array value NaN first at cell index (e.g. gcrc) & month (1-12): ", i, m+1);
            }
        }
    }
}
//365 days
// (3.2.1) Less Than (LT) 0
void check_array2Dng365_double_LT0(double array[ng][365]) {
    for (long i=0; i<ng; i++) {
        for (long d=0; d<365; d++) {
            if (array[i][d] < 0) {
                throw ExceptionValueLT0_array2Dng365__n_day("check_array2Dng365_double_LT0: value less than zero 0 first at cell index (e.g. gcrc) & day (1-365): ", i, d+1);
            }
        }
    }
}
// (3.2.2) Not a Number (NaN)
void check_array2Dng365_double_NaN(double array[ng][365]) {
    for (long i=0; i<ng; i++) {
        for (long d=0; d<365; d++) {
            if (array[i][d] != array[i][d]) {
                throw ExceptionValueNaN_array2Dng365__n_day("check_array2Dng365_double_NaN: 2D array value NaN first at cell index (e.g. gcrc) & day (1-365): ", i, d+1);
            }
        }
    }
}


// (4) Listing values LT0 for arrays 2D = ng*12 or ng*365
// 12 months
void list_array2Dng12_double_LT0__strgmsg(double array[ng][12], string const &s) {
//    cout << "list_array2Dng12_double_LT0__strgmsg started" << endl;
    for (long i=0; i<ng; i++) {
        for (long m=0; m<11; m++) {
            // check
//            if (500 == i && 9 == m) {
//                cout << "i: " << i << " m: " << m << " array[i][m]: " << array[i][m] << endl;
//            }
            if (array[i][m] < 0) {
                cout << s << " - value less than zero 0: cell index " << i << " month m+1 (1-12): " << m+1 << " array[i][m]: " << array[i][m] << endl;
            }
        }
    }
}
//365 days
void list_array2Dng365_double_LT0__strgmsg(double array[ng][365], string const &s) {
//    cout << "list_array2Dng365_double_LT0__strgmsg started" << endl;
    for (long i=0; i<ng; i++) {
        for (long d=0; d<365; d++) {
            // check
//            if (500 == i && 302 == d) {
//                cout << "i: " << i << " d: " << d << " array[i][d]: " << array[i][d] << endl;
//            }
            if (array[i][d] < 0) {
                cout << s << " - value less than zero 0: cell index " << i << " day d+1 (1-365): " << d+1 << " array[i][d]: " << array[i][d] << endl;
            }
        }
    }
}

// (5) Listing NaN LT0 for arrays 2D = ng*12 or ng*365
// 12 months
void list_array2Dng12_double_NaN__strgmsg(double array[ng][12], string const &s) {
//    cout << "list_array2Dng12_double_NaN__strgmsg started" << endl;
    for (long i=0; i<ng; i++) {
        for (long m=0; m<11; m++) {
            // check
//            if (500 == i && 9 == m) {
//                cout << "i: " << i << " m: " << m << " array[i][m]: " << array[i][m] << endl;
//            }
            if (array[i][m] != array[i][m]) {
                cout << s << " - value NaN: cell index " << i << " month m+1 (1-12): " << m+1 << " array[i][m]: " << array[i][m] << endl;
            }
        }
    }
}
//365 days
void list_array2Dng365_double_NaN__strgmsg(double array[ng][365], string const &s) {
    //    cout << "list_array2Dng365_double_NaN__strgmsg started" << endl;
    for (long i=0; i<ng; i++) {
        for (long d=0; d<365; d++) {
            // check
//            if (500 == i && 302 == d) {
//                cout << "i: " << i << " d: " << d << " array[i][d]: " << array[i][d] << endl;
//            }
            if (array[i][d] != array[i][d]) {
                cout << s << " - value NaN: cell index " << i << " day d+1 (1-365): " << d+1 << " array[i][d]: " << array[i][d] << endl;
            }
        }
    }
}

#endif // WATERGAP_CHECKS_GLOBAL_H
