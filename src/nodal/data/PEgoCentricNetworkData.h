/*
 * PEgoCentricNetworkData.h
 *
 *  Created on: Jun 15, 2010
 *      Author: duyvu
 */

#ifndef PEGOCENTRICNETWORKDATA_H_
#define PEGOCENTRICNETWORKDATA_H_

#include "EgoCentricNetworkData.h"

namespace ndip {

class PEgoCentricNetworkData: public ndip::EgoCentricNetworkData {
public:
	// num of nodal network statistics
	static const int NUMBER_OF_NODAL_NETWORK_STATISTICS = 1;

	PEgoCentricNetworkData();
	PEgoCentricNetworkData(const ExposureTimeVectorType& _exposureTimes,
			const EdgeTimeVectorType& _edgeTimes, double _observationTime);

	virtual ~PEgoCentricNetworkData();

	// roll a new edge and update nodal network statistics using cache structures
	void updateNodalNetworkStatistics(const EdgeTime& edgeTime);
	int getFirstOrderPopularityIndex();
};

}

#endif /* PEGOCENTRICNETWORKDATA_H_ */
