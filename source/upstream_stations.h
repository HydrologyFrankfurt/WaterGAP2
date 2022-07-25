class upstreamStationClass {
  private:
	void deleteAllListElements();
	void findStations(char *routing_dir, short allStationOption);

	struct listElement {
		short station;
		listElement *next;
	};

	listElement **upstreamStationList;
	short *numberOfUpstreamStations;
	int numberOfBasins;

  public:
	 upstreamStationClass();
	~upstreamStationClass();
	void findAllStations(char *routing_dir);
	void findDirectStations(char *routing_dir);
	short getNumberOfUpstreamStations(int n);
	short getUpstreamStation(int n, short i);
	void writeListToFile(char *filename, char *filename_no);

        void replaceGridValues(int st, float newValue, double *Grid, float *upValueList);
};
