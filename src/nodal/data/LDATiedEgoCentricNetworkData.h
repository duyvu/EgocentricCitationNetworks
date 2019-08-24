/*
 * LDATiedEgoCentricNetworkData.h
 *
 *  Created on: Jan 21, 2011
 *      Author: duyvu
 */

#ifndef LDATIEDEGOCENTRICNETWORKDATA_H_
#define LDATIEDEGOCENTRICNETWORKDATA_H_

#include "TiedEgoCentricNetworkData.h"

namespace ndip {

class LDATiedEgoCentricNetworkData: public ndip::TiedEgoCentricNetworkData {
protected:
	DoubleMatrix LDA;
public:
	static const unsigned int NUMBER_OF_LDA_TOPICS = 50;
	static const unsigned int NUMBER_OF_NODAL_NETWORK_STATISTICS =
			NUMBER_OF_LDA_TOPICS;

	LDATiedEgoCentricNetworkData();
	LDATiedEgoCentricNetworkData(
			const ExposureTimeVectorType& _nodeJoiningTimes,
			const VectorOfEdgeEvents& _vectorOfEdgeEvents,
			double _observationTimeStart, double _observationTimeEnd,
			const std::vector<string>& listOfNodalDataFiles);
	virtual ~LDATiedEgoCentricNetworkData();

	void updateNodalNetworkStatistics(const EdgeEvents& edgeEvents) {
		cout << "This is not supported in LDA data" << endl;
	}
	void updateNodalNetworkStatistics(unsigned int indexOfEdgeEvents);

	int getFirstOrderPopularityIndex();

};

}

#endif /* LDATIEDEGOCENTRICNETWORKDATA_H_ */
