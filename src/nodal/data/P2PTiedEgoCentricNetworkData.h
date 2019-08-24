/*
 * P2PTiedEgoCentricNetworkData.h
 *
 *  Created on: Jul 20, 2010
 *      Author: duyvu
 */

#ifndef P2PTIEDEGOCENTRICNETWORKDATA_H_
#define P2PTIEDEGOCENTRICNETWORKDATA_H_

#include "TiedEgoCentricNetworkData.h"

namespace ndip {

class P2PTiedEgoCentricNetworkData: public ndip::TiedEgoCentricNetworkData {
protected:
	void updateNodalNetworkStatisticsCache(const EdgeEvents& edgeEvents,
			bool isCached);
public:
	// num of nodal network statistics
	static const int NUMBER_OF_NODAL_NETWORK_STATISTICS = 2;

	P2PTiedEgoCentricNetworkData();
	P2PTiedEgoCentricNetworkData(
			const ExposureTimeVectorType& _nodeJoiningTimes,
			const VectorOfEdgeEvents& _vectorOfEdgeEvents,
			double _observationTimeStart, double _observationTimeEnd);
	virtual ~P2PTiedEgoCentricNetworkData();

	// roll a new edge and update nodal network statistics using cache structures
	void updateNodalNetworkStatistics(const EdgeEvents& edgeEvents);
	int getFirstOrderPopularityIndex();
};

}

#endif /* P2PTIEDEGOCENTRICNETWORKDATA_H_ */
