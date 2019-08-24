/*
 * PTiedEgoCentricNetworkData.h
 *
 *  Created on: Jul 18, 2010
 *      Author: duyvu
 */

#ifndef PTIEDEGOCENTRICNETWORKDATA_H_
#define PTIEDEGOCENTRICNETWORKDATA_H_

#include "TiedEgoCentricNetworkData.h"

namespace ndip {

class PTiedEgoCentricNetworkData: public ndip::TiedEgoCentricNetworkData {
public:
	// num of nodal network statistics
	static const int NUMBER_OF_NODAL_NETWORK_STATISTICS = 1;

	PTiedEgoCentricNetworkData();
	PTiedEgoCentricNetworkData(const ExposureTimeVectorType& _nodeJoiningTimes,
			const VectorOfEdgeEvents& _vectorOfEdgeEvents,
			double _observationTimeStart, double _observationTimeEnd);
	virtual ~PTiedEgoCentricNetworkData();

	// roll a new edge and update nodal network statistics using cache structures
	void updateNodalNetworkStatistics(const EdgeEvents& edgeEvents);
	int getFirstOrderPopularityIndex();
};

}

#endif /* PTIEDEGOCENTRICNETWORKDATA_H_ */
