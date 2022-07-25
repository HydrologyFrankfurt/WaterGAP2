class dailyPrecNrdClass {
  public:
	void init(int idum, const char *output_dir);
	float getDailyPrec(const short day_in_month,
					   const short month,
					   const char numberOfRainDays, const short monthlyPrecipitation);

  private:
	// random number generator
	float ran1(int *idum);

	// tables which contain a series of 0 and 1 for 
	// each possible combination of 'number of raindays per month'
	// and 'number of days per month'
	// 0: if day is not a rainday
	// 1: if day is a rainday
	// first counter is for the number of the day within the month
	// second counter is for the different values of 'raindays per month'
	char rain_days28[28][28];
	char rain_days30[30][30];
	char rain_days31[31][31];
};
