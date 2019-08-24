/*
 * P2PTRTiedEgoCentricNetworkData.h
 *
 *  Created on: Jul 23, 2010
 *      Author: duyvu
 */

#ifndef P2PTRTIEDEGOCENTRICNETWORKDATA_H_
#define P2PTRTIEDEGOCENTRICNETWORKDATA_H_

#include "TiedEgoCentricNetworkData.h"

namespace ndip {

class P2PTRTiedEgoCentricNetworkData: public ndip::TiedEgoCentricNetworkData {
protected:
	unsigned int windowSizeOfRenewalStatistic;
	void updateNodalNetworkStatisticsCache(const EdgeEvents& edgeEvents,
			double nextCitingTime, bool isCached,
			ListOfEdgeEvents& listOfActiveEvents);
public:
	// num of nodal network statistics
	// index 2 for brokee
	// index 3 for broker
	// index 4 for blocker
	// index 5 for first out-degree
	// index 6 for first out-degree
	// index 7 for renewal
	static const int NUMBER_OF_NODAL_NETWORK_STATISTICS = 8;
	static const int WINDOW_SIZE_OF_RENEWAL_STATISTIC = 180;

	P2PTRTiedEgoCentricNetworkData();
	P2PTRTiedEgoCentricNetworkData(
			const ExposureTimeVectorType& _nodeJoiningTimes,
			const VectorOfEdgeEvents& _vectorOfEdgeEvents,
			double _observationTimeStart, double _observationTimeEnd);
	virtual ~P2PTRTiedEgoCentricNetworkData();

	// roll a new edge and update nodal network statistics using cache structures
	void updateNodalNetworkStatistics(const EdgeEvents& edgeEvents);
	int getFirstOrderPopularityIndex();

	void setWindowSizeOfRenewalStatistic(unsigned int size) {
		windowSizeOfRenewalStatistic = size;
		cout << "Renewal statistic window size is "
				<< windowSizeOfRenewalStatistic << endl;
	}
	unsigned int getWindowSizeOfRenewalStatistic() {
		return windowSizeOfRenewalStatistic;
	}
};

}

#endif /* P2PTRTIEDEGOCENTRICNETWORKDATA_H_ */
