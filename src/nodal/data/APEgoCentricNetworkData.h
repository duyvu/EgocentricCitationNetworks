/*
 * APEgoCentricNetworkData.h
 *
 *  Created on: May 26, 2010
 *      Author: duyvu
 */

#ifndef APEGOCENTRICNETWORKDATA_H_
#define APEGOCENTRICNETWORKDATA_H_

#include "EgoCentricNetworkData.h"

namespace ndip {

class APEgoCentricNetworkData: public ndip::EgoCentricNetworkData {
protected:
	void updateNodalNetworkStatisticsCache(const EdgeTime& edgeTime,
			bool isCached);
public:
	// num of nodal network statistics
	static const int NUMBER_OF_NODAL_NETWORK_STATISTICS = 4;

	APEgoCentricNetworkData();
	APEgoCentricNetworkData(const ExposureTimeVectorType& _exposureTimes,
			const EdgeTimeVectorType& _edgeTimes, double _observationTime);

	virtual ~APEgoCentricNetworkData();

	// roll a new edge and update nodal network statistics using cache structures
	void updateNodalNetworkStatistics(const EdgeTime& edgeTime);
	int getFirstOrderPopularityIndex();
};

}

#endif /* APEGOCENTRICNETWORKDATA_H_ */
