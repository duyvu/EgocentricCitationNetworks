/*
 * P2PTTiedEgoCentricNetworkData.h
 *
 *  Created on: Jul 21, 2010
 *      Author: duyvu
 */

#ifndef P2PTTIEDEGOCENTRICNETWORKDATA_H_
#define P2PTTIEDEGOCENTRICNETWORKDATA_H_

#include "TiedEgoCentricNetworkData.h"

namespace ndip {

class P2PTTiedEgoCentricNetworkData: public ndip::TiedEgoCentricNetworkData {
protected:
	void updateNodalNetworkStatisticsCache(const EdgeEvents& edgeEvents,
			bool isCached);
public:
	// num of nodal network statistics
	// index 1 for first in-degree
	// index 1 for second in-degree
	// index 2 for seller
	// index 3 for broker
	// index 4 for buyer
	// index 5 for first out-degree
	// index 6 for second out-degree
	static const int NUMBER_OF_NODAL_NETWORK_STATISTICS = 7;

	P2PTTiedEgoCentricNetworkData();
	P2PTTiedEgoCentricNetworkData(
			const ExposureTimeVectorType& _nodeJoiningTimes,
			const VectorOfEdgeEvents& _vectorOfEdgeEvents,
			double _observationTimeStart, double _observationTimeEnd);
	virtual ~P2PTTiedEgoCentricNetworkData();

	// roll a new edge and update nodal network statistics using cache structures
	void updateNodalNetworkStatistics(const EdgeEvents& edgeEvents);
	int getFirstOrderPopularityIndex();
};

}

#endif /* P2PTTIEDEGOCENTRICNETWORKDATA_H_ */
