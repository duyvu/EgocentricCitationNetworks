/*
 * PGroupedEgoCentricNetworkData.h
 *
 *  Created on: Jan 17, 2011
 *      Author: duyvu
 */

#ifndef PGROUPEDEGOCENTRICNETWORKDATA_H_
#define PGROUPEDEGOCENTRICNETWORKDATA_H_

#include "GroupedEgoCentricNetworkData.h"

namespace ndip {

class PGroupedEgoCentricNetworkData: public ndip::GroupedEgoCentricNetworkData {
public:
	// num of nodal network statistics
	static const int NUMBER_OF_NODAL_NETWORK_STATISTICS = 1;

	PGroupedEgoCentricNetworkData();
	PGroupedEgoCentricNetworkData(
			const ExposureTimeVectorType& _nodeJoiningTimes,
			const VectorOfDiscreteTimeVectorOfEdgeEvents& _vectorOfEdgeEvents,
			double _observationTimeStart, double _observationTimeEnd);
	virtual ~PGroupedEgoCentricNetworkData();

	// roll a new edge and update nodal network statistics using cache structures
	void updateNodalNetworkStatistics(
			const DiscreteTimeVectorOfEdgeEvents& edgeEvents);
	int getFirstOrderPopularityIndex();
};

}

#endif /* PGROUPEDEGOCENTRICNETWORKDATA_H_ */
