/*
 * PGTiedEgoCentricNetworkData.h
 *
 *  Created on: Jan 22, 2011
 *      Author: duyvu
 */

#ifndef PGTIEDEGOCENTRICNETWORKDATA_H_
#define PGTIEDEGOCENTRICNETWORKDATA_H_

#include "TiedEgoCentricNetworkData.h"

namespace ndip {

class PGTiedEgoCentricNetworkData: public ndip::TiedEgoCentricNetworkData {
public:
	// num of nodal network statistics
	static const int NUMBER_OF_NODAL_NETWORK_STATISTICS = 2;

	PGTiedEgoCentricNetworkData();
	PGTiedEgoCentricNetworkData(
			const ExposureTimeVectorType& _nodeJoiningTimes,
			const VectorOfEdgeEvents& _vectorOfEdgeEvents,
			double _observationTimeStart, double _observationTimeEnd);
	virtual ~PGTiedEgoCentricNetworkData();

	// roll a new edge and update nodal network statistics using cache structures
	void updateNodalNetworkStatistics(const EdgeEvents& edgeEvents);
	int getFirstOrderPopularityIndex() {
		return 0;
	}

};

}

#endif /* PGTIEDEGOCENTRICNETWORKDATA_H_ */
