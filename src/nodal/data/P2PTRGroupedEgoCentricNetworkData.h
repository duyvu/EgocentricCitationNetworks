/*
 * P2PTRGroupedEgoCentricNetworkData.h
 *
 *  Created on: Jan 19, 2011
 *      Author: duyvu
 */

#ifndef P2PTRGROUPEDEGOCENTRICNETWORKDATA_H_
#define P2PTRGROUPEDEGOCENTRICNETWORKDATA_H_

#include "GroupedEgoCentricNetworkData.h"

namespace ndip {

class P2PTRGroupedEgoCentricNetworkData: public ndip::GroupedEgoCentricNetworkData {
protected:
	unsigned int windowSizeOfRenewalStatistic;
	void updateNodalNetworkStatisticsCache(
			const DiscreteTimeVectorOfEdgeEvents& edgeEvents,
			double nextCitingTime, bool isCached,
			ListOfDiscreteTimeVectorOfEdgeEvents& listOfActiveEvents);
public:
	// num of nodal network statistics
	// index 1 for first in-degree
	// index 1 for second in-degree
	// index 2 for seller
	// index 3 for broker
	// index 4 for buyer
	// index 5 for first out-degree
	// index 6 for second out-degree
	// index 7 for recency
	static const int NUMBER_OF_NODAL_NETWORK_STATISTICS = 8;
	static const int WINDOW_SIZE_OF_RENEWAL_STATISTIC = 180;

	P2PTRGroupedEgoCentricNetworkData();
	P2PTRGroupedEgoCentricNetworkData(
			const ExposureTimeVectorType& _nodeJoiningTimes,
			const VectorOfDiscreteTimeVectorOfEdgeEvents& _vectorOfEdgeEvents,
			double _observationTimeStart, double _observationTimeEnd);
	virtual ~P2PTRGroupedEgoCentricNetworkData();

	// roll a new edge and update nodal network statistics using cache structures
	void updateNodalNetworkStatistics(
			const DiscreteTimeVectorOfEdgeEvents& edgeEvents);
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

#endif /* P2PTRGROUPEDEGOCENTRICNETWORKDATA_H_ */
