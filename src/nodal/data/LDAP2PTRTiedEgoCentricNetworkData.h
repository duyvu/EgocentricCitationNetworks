/*
 * LDAP2PTRTiedEgoCentricNetworkData.h
 *
 *  Created on: Jan 13, 2011
 *      Author: duyvu
 */

#ifndef LDAP2PTRTIEDEGOCENTRICNETWORKDATA_H_
#define LDAP2PTRTIEDEGOCENTRICNETWORKDATA_H_

#include "TiedEgoCentricNetworkData.h"

namespace ndip {

class LDAP2PTRTiedEgoCentricNetworkData: public ndip::TiedEgoCentricNetworkData {
protected:
	unsigned int windowSizeOfRenewalStatistic;
	void updateNodalNetworkStatisticsCache(const EdgeEvents& edgeEvents,
			double nextCitingTime, bool isCached,
			ListOfEdgeEvents& listOfActiveEvents);
	DoubleMatrix LDA;
public:
	static const int NUMBER_OF_LDA_TOPICS = 50;
	static const int NUMBER_OF_NODAL_NETWORK_STATISTICS = 8
			+ NUMBER_OF_LDA_TOPICS;
	static const int WINDOW_SIZE_OF_RENEWAL_STATISTIC = 180;

	LDAP2PTRTiedEgoCentricNetworkData();
	LDAP2PTRTiedEgoCentricNetworkData(
			const ExposureTimeVectorType& _nodeJoiningTimes,
			const VectorOfEdgeEvents& _vectorOfEdgeEvents,
			double _observationTimeStart, double _observationTimeEnd,
			const std::vector<string>& listOfNodalDataFiles);
	virtual ~LDAP2PTRTiedEgoCentricNetworkData();

	void updateNodalNetworkStatistics(const EdgeEvents& edgeEvents) {
		cout << "This is not supported in LDA data" << endl;
	}
	void updateNodalNetworkStatistics(unsigned int indexOfEdgeEvents);

	int getFirstOrderPopularityIndex();

	void setWindowSizeOfRenewalStatistic(unsigned int size) {
		windowSizeOfRenewalStatistic = size;
		cout << "Renewal statistic window size is "
				<< windowSizeOfRenewalStatistic << endl;
	}
	unsigned int getWindowSizeOfRenewalStatistic() {
		return windowSizeOfRenewalStatistic;
	}

};

}

#endif /* LDAP2PTRTIEDEGOCENTRICNETWORKDATA_H_ */
