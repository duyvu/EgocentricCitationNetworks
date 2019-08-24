/*
 * TAPEgoCentricNetworkData.h
 *
 *  Created on: Jun 1, 2010
 *      Author: duyvu
 */

#ifndef TAPEGOCENTRICNETWORKDATA_H_
#define TAPEGOCENTRICNETWORKDATA_H_

#include "EgoCentricNetworkData.h"

namespace ndip {

class TAPEgoCentricNetworkData: public ndip::EgoCentricNetworkData {
protected:
	void updateNodalNetworkStatisticsCache(const EdgeTime& edgeTime,
			bool isCached);
public:
	// num of nodal network statistics
	static const int NUMBER_OF_NODAL_NETWORK_STATISTICS = 6;

	TAPEgoCentricNetworkData();
	TAPEgoCentricNetworkData(const ExposureTimeVectorType& _exposureTimes,
			const EdgeTimeVectorType& _edgeTimes, double _observationTime);

	virtual ~TAPEgoCentricNetworkData();

	// roll a new edge and update nodal network statistics using cache structures
	void updateNodalNetworkStatistics(const EdgeTime& edgeTime);
	int getFirstOrderPopularityIndex();
};

}

#endif /* TAPEGOCENTRICNETWORKDATA_H_ */
