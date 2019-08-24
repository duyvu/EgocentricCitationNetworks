/*
 * SimTAPEgoCentricNetworkData.h
 *
 *  Created on: Aug 10, 2010
 *      Author: duyvu
 */

#ifndef SIMTAPEGOCENTRICNETWORKDATA_H_
#define SIMTAPEGOCENTRICNETWORKDATA_H_

#include "EgoCentricNetworkData.h"

namespace ndip {

class SimTAPEgoCentricNetworkData: public ndip::EgoCentricNetworkData {
protected:
	void updateNetworkStatisticsWithoutCache(const EdgeTime& edgeTime);
public:
	// num of nodal network statistics
	static const int NUMBER_OF_NODAL_NETWORK_STATISTICS = 6;

	SimTAPEgoCentricNetworkData();
	SimTAPEgoCentricNetworkData(const ExposureTimeVectorType& _exposureTimes,
			const EdgeTimeVectorType& _edgeTimes, double _observationTime);

	virtual ~SimTAPEgoCentricNetworkData();

	// roll a new edge and update nodal network statistics using cache structures
	void updateNodalNetworkStatistics(const EdgeTime& edgeTime);
	int getFirstOrderPopularityIndex();
};

}

#endif /* SIMTAPEGOCENTRICNETWORKDATA_H_ */
