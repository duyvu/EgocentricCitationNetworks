/*
 * P2PTGroupedEgoCentricNetworkData.h
 *
 *  Created on: Jan 18, 2011
 *      Author: duyvu
 */

#ifndef P2PTGROUPEDEGOCENTRICNETWORKDATA_H_
#define P2PTGROUPEDEGOCENTRICNETWORKDATA_H_

#include "GroupedEgoCentricNetworkData.h"

namespace ndip {

class P2PTGroupedEgoCentricNetworkData: public ndip::GroupedEgoCentricNetworkData {
protected:
	void updateNodalNetworkStatisticsCache(
			const DiscreteTimeVectorOfEdgeEvents& edgeEvents, bool isCached);
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

	P2PTGroupedEgoCentricNetworkData();
	P2PTGroupedEgoCentricNetworkData(
			const ExposureTimeVectorType& _nodeJoiningTimes,
			const VectorOfDiscreteTimeVectorOfEdgeEvents& _vectorOfEdgeEvents,
			double _observationTimeStart, double _observationTimeEnd);
	virtual ~P2PTGroupedEgoCentricNetworkData();

	// roll a new edge and update nodal network statistics using cache structures
	void updateNodalNetworkStatistics(
			const DiscreteTimeVectorOfEdgeEvents& edgeEvents);
	int getFirstOrderPopularityIndex();
};

}

#endif /* P2PTGROUPEDEGOCENTRICNETWORKDATA_H_ */
