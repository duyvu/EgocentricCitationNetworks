/*
 * PTiedEgoCentricNetworkData.cpp
 *
 *  Created on: Jul 18, 2010
 *      Author: duyvu
 */

#include "PTiedEgoCentricNetworkData.h"

namespace ndip {

PTiedEgoCentricNetworkData::PTiedEgoCentricNetworkData() {
	// TODO Auto-generated constructor stub

}

PTiedEgoCentricNetworkData::PTiedEgoCentricNetworkData(
		const ExposureTimeVectorType& _nodeJoiningTimes,
		const VectorOfEdgeEvents& _vectorOfEdgeEvents,
		double _observationTimeStart, double _observationTimeEnd) :
	TiedEgoCentricNetworkData(_nodeJoiningTimes, _vectorOfEdgeEvents,
			_observationTimeStart, _observationTimeEnd) {

	cout << "START the constructor of PTiedEgoCentricNetworkData" << endl;

	// Initialize the current nodal network statistics
	numOfNodalNetworkStatistics
			= PTiedEgoCentricNetworkData::NUMBER_OF_NODAL_NETWORK_STATISTICS;
	currentNodalNetworkStatistics = SparseUnsignedIntMatrix(numOfVertices,
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
		EdgeEvents edgeEvents = vectorOfEdgeEvents[i];
		this->updateNodalNetworkStatistics(edgeEvents);
	}

	// Copy these current nodal network statistics and the graph up to the observation time
	// to save the time for iterative calls later
	beginningNodalNetworkStatistics = SparseUnsignedIntMatrix(
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
		EdgeEvents edgeEvents = vectorOfEdgeEvents[i];
		vectorOfNodalUpdateMaps(currentCacheIndex) = NodalUpdateMap();
		for (unsigned int i = 0; i < edgeEvents.second.size(); i++) {
			unsigned int citedNode = edgeEvents.second[i];
			NodalUpdates myNodalUpdates = NodalUpdates();
			myNodalUpdates.push_back(StatValue(0,
					currentNodalNetworkStatistics(citedNode, 0)));
			vectorOfNodalUpdateMaps(currentCacheIndex).insert(make_pair(
					citedNode, myNodalUpdates));
		}
	}

	// reset everything for the next use
	dirty = true;
	finish();

	cout << "FINISH the constructor of PTiedEgoCentricNetworkData" << endl;
}

PTiedEgoCentricNetworkData::~PTiedEgoCentricNetworkData() {
	// TODO Auto-generated destructor stub
}

// roll a new edge and update nodal network statistics using cache structures
void PTiedEgoCentricNetworkData::updateNodalNetworkStatistics(
		const EdgeEvents& edgeEvents) {

	unsigned int citingNode = edgeEvents.first;
	double citingTime = nodeJoiningTimes[citingNode];

	for (unsigned int i = 0; i < edgeEvents.second.size(); i++) {
		unsigned int citedNode = edgeEvents.second[i];
		currentRecentNodeCitedTimes[citedNode] = citingTime;
		currentNodalNetworkStatistics(citedNode, 0) += 1;
	}

	currentCacheIndex++;
}

int PTiedEgoCentricNetworkData::getFirstOrderPopularityIndex() {
	return 0;
}

}
