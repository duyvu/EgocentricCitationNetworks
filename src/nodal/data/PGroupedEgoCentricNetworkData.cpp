/*
 * PGroupedEgoCentricNetworkData.cpp
 *
 *  Created on: Jan 17, 2011
 *      Author: duyvu
 */

#include "PGroupedEgoCentricNetworkData.h"

namespace ndip {

PGroupedEgoCentricNetworkData::PGroupedEgoCentricNetworkData() {
	// TODO Auto-generated constructor stub

}

PGroupedEgoCentricNetworkData::PGroupedEgoCentricNetworkData(
		const ExposureTimeVectorType& _nodeJoiningTimes,
		const VectorOfDiscreteTimeVectorOfEdgeEvents& _vectorOfEdgeEvents,
		double _observationTimeStart, double _observationTimeEnd) :
	GroupedEgoCentricNetworkData(_nodeJoiningTimes, _vectorOfEdgeEvents,
			_observationTimeStart, _observationTimeEnd) {

	cout << "START the constructor of PGroupedEgoCentricNetworkData" << endl;

	// Initialize the current nodal network statistics
	numOfNodalNetworkStatistics
			= PGroupedEgoCentricNetworkData::NUMBER_OF_NODAL_NETWORK_STATISTICS;
	currentNodalNetworkStatistics = UnsignedIntMatrix(numOfVertices,
			numOfNodalNetworkStatistics);
	currentNodalNetworkStatistics.clear();

	// Initialize the recent node contacted times
	currentRecentNodeCitedTimes = ExposureTimeVectorType(numOfVertices);
	for (unsigned int i = 0; i < numOfVertices; i++)
		currentRecentNodeCitedTimes[i] = nodeJoiningTimes[i];

	// Roll the current nodal network statistics
	// and the graph up to the observation time
	for (unsigned int i = 0; i < beginningEdgeIndex; i++) {
		//cout << i << endl;
		DiscreteTimeVectorOfEdgeEvents edgeEvents =
				vectorOfDiscreteTimeVectorOfEdgeEvents[i];
		this->updateNodalNetworkStatistics(edgeEvents);
	}

	// Copy these current nodal network statistics and the graph up to the observation time
	// to save the time for iterative calls later
	beginningNodalNetworkStatistics = UnsignedIntMatrix(
			currentNodalNetworkStatistics);
	beginningRecentNodeCitedTimes = ExposureTimeVectorType(
			currentRecentNodeCitedTimes);

	// initialize the cache structure
	int cacheSize = endingEdgeIndex - beginningEdgeIndex + 1;
	cout << "The size of cache is " << cacheSize << endl;
	vectorOfNodalUpdateMaps = VectorOfNodalUpdateMaps(cacheSize);

	// this corresponds to the first edge join the network at the observation time
	// vectorOfComingVertices[0]: the list of nodes join the network from the observation time (the first edge time) up to
	// the time of next edge.
	// vectorOfNodalUpdateMaps[0]: the list of nodes updated after entering the first edge to the network.
	// roll up to the end and cache data
	currentCacheIndex = 0;
	for (unsigned int i = beginningEdgeIndex; i <= endingEdgeIndex; i++, currentCacheIndex++) {
		//cout << i << endl;
		DiscreteTimeVectorOfEdgeEvents edgeEvents =
				vectorOfDiscreteTimeVectorOfEdgeEvents[i];
		vectorOfNodalUpdateMaps(currentCacheIndex) = NodalUpdateMap();
		for (unsigned int k = 0; k < edgeEvents.second.size(); k++) {
			for (unsigned int i = 0; i < edgeEvents.second[k].second.size(); i++) {
				unsigned int citedNode = edgeEvents.second[k].second[i];
				NodalUpdateMap::iterator iter = vectorOfNodalUpdateMaps(
						currentCacheIndex).find(citedNode);
				if (iter != vectorOfNodalUpdateMaps(currentCacheIndex).end()) {
					NodalUpdates myNodalUpdates = ((NodalUpdates) iter->second);
					myNodalUpdates.push_back(StatValue(0,
							currentNodalNetworkStatistics(citedNode, 0)));
					vectorOfNodalUpdateMaps(currentCacheIndex)[citedNode]
							= myNodalUpdates;
				} else {
					NodalUpdates myNodalUpdates = NodalUpdates();
					myNodalUpdates.push_back(StatValue(0,
							currentNodalNetworkStatistics(citedNode, 0)));
					vectorOfNodalUpdateMaps(currentCacheIndex).insert(
							make_pair(citedNode, myNodalUpdates));
				}
			}
		}
	}

	// reset everything for the next use
	dirty = true;
	finish();

	cout << "FINISH the constructor of PGroupedEgoCentricNetworkData" << endl;
}

PGroupedEgoCentricNetworkData::~PGroupedEgoCentricNetworkData() {
	// TODO Auto-generated destructor stub
}

// roll a new edge and update nodal network statistics using cache structures
void PGroupedEgoCentricNetworkData::updateNodalNetworkStatistics(
		const DiscreteTimeVectorOfEdgeEvents& edgeEvents) {
	double citingTime = edgeEvents.first;
	for (unsigned int k = 0; k < edgeEvents.second.size(); k++) {
		for (unsigned int i = 0; i < edgeEvents.second[k].second.size(); i++) {
			unsigned int citedNode = edgeEvents.second[k].second[i];
			currentRecentNodeCitedTimes[citedNode] = citingTime;
			currentNodalNetworkStatistics(citedNode, 0) += 1;
		}
	}
	currentCacheIndex++;
}

int PGroupedEgoCentricNetworkData::getFirstOrderPopularityIndex() {
	return 0;
}

}
