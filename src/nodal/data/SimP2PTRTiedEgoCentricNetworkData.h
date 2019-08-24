/*
 * SimP2PTRTiedEgoCentricNetworkData.h
 *
 *  Created on: Aug 14, 2010
 *      Author: duyvu
 */

#ifndef SIMP2PTRTIEDEGOCENTRICNETWORKDATA_H_
#define SIMP2PTRTIEDEGOCENTRICNETWORKDATA_H_

#include "TiedEgoCentricNetworkData.h"

namespace ndip {

class SimP2PTRTiedEgoCentricNetworkData: public ndip::TiedEgoCentricNetworkData {
protected:
	ListOfEdgeEvents listOfActiveEvents;
	unsigned int windowSizeOfRenewalStatistic;
	void updateNetworkStatisticsWithoutCache(const EdgeEvents& edgeEvents);
public:
	static const int NUMBER_OF_NODAL_NETWORK_STATISTICS = 8;
	static const int WINDOW_SIZE_OF_RENEWAL_STATISTIC = 180;

	SimP2PTRTiedEgoCentricNetworkData();
	SimP2PTRTiedEgoCentricNetworkData(
			const ExposureTimeVectorType& _nodeJoiningTimes,
			const VectorOfEdgeEvents& _vectorOfEdgeEvents,
			double _observationTimeStart, double _observationTimeEnd);
	virtual ~SimP2PTRTiedEgoCentricNetworkData();

	void updateNodalNetworkStatistics(const EdgeEvents& edgeEvents);
	int getFirstOrderPopularityIndex() {
		return 0;
	}

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

#endif /* SIMP2PTRTIEDEGOCENTRICNETWORKDATA_H_ */
