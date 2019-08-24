/*
 * TextualP2PTRTiedEgoCentricNetworkData.h
 *
 *  Created on: Aug 17, 2010
 *      Author: duyvu
 */

#ifndef TEXTUALP2PTRTIEDEGOCENTRICNETWORKDATA_H_
#define TEXTUALP2PTRTIEDEGOCENTRICNETWORKDATA_H_

#include "TiedEgoCentricNetworkData.h"

namespace ndip {

class TextualP2PTRTiedEgoCentricNetworkData: public ndip::TiedEgoCentricNetworkData {
protected:
	unsigned int windowSizeOfRenewalStatistic;
	void updateNodalNetworkStatisticsCache(const EdgeEvents& edgeEvents,
			double nextCitingTime, bool isCached,
			ListOfEdgeEvents& listOfActiveEvents);
public:
	static const int NUMBER_OF_NODAL_NETWORK_STATISTICS = 9;
	static const int WINDOW_SIZE_OF_RENEWAL_STATISTIC = 180;

	TextualP2PTRTiedEgoCentricNetworkData();
	TextualP2PTRTiedEgoCentricNetworkData(
			const ExposureTimeVectorType& _nodeJoiningTimes,
			const VectorOfEdgeEvents& _vectorOfEdgeEvents,
			double _observationTimeStart, double _observationTimeEnd,
			const std::vector<string>& listOfNodalDataFiles);
	virtual ~TextualP2PTRTiedEgoCentricNetworkData();

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

#endif /* TEXTUALP2PTRTIEDEGOCENTRICNETWORKDATA_H_ */
