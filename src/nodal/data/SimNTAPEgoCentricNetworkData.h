/*
 * SimNTAPEgoCentricNetworkData.h
 *
 *  Created on: Jul 6, 2010
 *      Author: duyvu
 */

#ifndef SIMNTAPEGOCENTRICNETWORKDATA_H_
#define SIMNTAPEGOCENTRICNETWORKDATA_H_

#include "EgoCentricNetworkData.h"

namespace ndip {

class SimNTAPEgoCentricNetworkData: public ndip::EgoCentricNetworkData {
protected:
	void updateNetworkStatisticsWithoutCache(const EdgeTime& edgeTime);
public:
	// num of nodal network statistics
	static const int NUMBER_OF_NODAL_NETWORK_STATISTICS = 10;

	SimNTAPEgoCentricNetworkData();
	SimNTAPEgoCentricNetworkData(const ExposureTimeVectorType& _exposureTimes,
			const EdgeTimeVectorType& _edgeTimes, double _observationTime,
			const std::vector<string>& listOfNodalDataFiles);

	virtual ~SimNTAPEgoCentricNetworkData();

	// roll a new edge and update nodal network statistics using cache structures
	void updateNodalNetworkStatistics(const EdgeTime& edgeTime);
	int getFirstOrderPopularityIndex();
};

}

#endif /* SIMNTAPEGOCENTRICNETWORKDATA_H_ */
