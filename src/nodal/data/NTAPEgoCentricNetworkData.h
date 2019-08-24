/*
 * NTAPEgoCentricNetworkData.h
 *
 *  Created on: Jun 20, 2010
 *      Author: duyvu
 */

#ifndef NTAPEGOCENTRICNETWORKDATA_H_
#define NTAPEGOCENTRICNETWORKDATA_H_

#include "TAPEgoCentricNetworkData.h"

namespace ndip {

class NTAPEgoCentricNetworkData: public ndip::EgoCentricNetworkData {
protected:
	void updateNodalNetworkStatisticsCache(const EdgeTime& edgeTime,
			bool isCached);
public:
	// num of nodal network statistics
	static const int NUMBER_OF_NODAL_NETWORK_STATISTICS = 10;

	NTAPEgoCentricNetworkData();
	NTAPEgoCentricNetworkData(const ExposureTimeVectorType& _exposureTimes,
			const EdgeTimeVectorType& _edgeTimes, double _observationTime,
			const std::vector<string>& listOfNodalDataFiles);

	virtual ~NTAPEgoCentricNetworkData();

	// roll a new edge and update nodal network statistics using cache structures
	void updateNodalNetworkStatistics(const EdgeTime& edgeTime);
	int getFirstOrderPopularityIndex();
};

}

#endif /* NTAPEGOCENTRICNETWORKDATA_H_ */
