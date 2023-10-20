#if !defined (_upstream_stations_h_)
#define _upstream_stations_h_

#include <string>
#include <vector>

class upstreamStationClass {
  	private:
		void deleteAllListElements();
		void findStations(const std::string routing_dir, short allStationOption);

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
		void findAllStations(const std::string routing_dir);
		void findDirectStations(const std::string routing_dir);
		short getNumberOfUpstreamStations(int n);
		short getUpstreamStation(int n, short i);
		void writeListToFile(const std::string filename, const std::string filename_no);
        void replaceGridValues(int st, float newValue, Grid<> & Grid, std::vector<float> upValueList);
};
#endif
