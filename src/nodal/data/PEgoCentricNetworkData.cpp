/*
 * PEgoCentricNetworkData.cpp
 *
 *  Created on: Jun 15, 2010
 *      Author: duyvu
 */

#include "PEgoCentricNetworkData.h"

namespace ndip {

PEgoCentricNetworkData::PEgoCentricNetworkData() {
	// TODO Auto-generated constructor stub

}

PEgoCentricNetworkData::PEgoCentricNetworkData(
		const ExposureTimeVectorType& _exposureTimes,
		const EdgeTimeVectorType& _edgeTimes, double _observationTime) :
	EgoCentricNetworkData(_exposureTimes, _edgeTimes, _observationTime) {

	// Initialize the current nodal network statistics
	numOfNodalNetworkStatistics
			= PEgoCentricNetworkData::NUMBER_OF_NODAL_NETWORK_STATISTICS;
	currentNodalNetworkStatistics = SparseUnsignedIntMatrix(numOfVertices,
			numOfNodalNetworkStatistics);
	currentNodalNetworkStatistics.clear();

	// Initialize the recent node contacting times
	currentRecentNodeContactingTimes = ExposureTimeVectorType(numOfVertices);
	for (unsigned int i = 0; i < numOfVertices; i++)
		currentRecentNodeContactingTimes[i] = nodeJoiningTimes[i];

	// Initialize the recent node contacted times
	currentRecentNodeContactedTimes = ExposureTimeVectorType(numOfVertices);
	for (unsigned int i = 0; i < numOfVertices; i++)
		currentRecentNodeContactedTimes[i] = nodeJoiningTimes[i];

	// Roll the current nodal network statistics
	// and the graph up to the observation time
	for (unsigned int i = 0; i < beginningEdgeIndex; i++) {
		//cout << i << endl;
		EdgeTime edgeTime = edgeJoiningTimes[i];
		this->updateNodalNetworkStatistics(edgeTime);
	}

	// Copy these current nodal network statistics and the graph up to the observation time
	// to save the time for iterative calls later
	beginningNodalNetworkStatistics = SparseUnsignedIntMatrix(
			currentNodalNetworkStatistics);
	beginningRecentNodeContactingTimes = ExposureTimeVectorType(
			currentRecentNodeContactingTimes);
	beginningRecentNodeContactedTimes = ExposureTimeVectorType(
			currentRecentNodeContactedTimes);

	// initialize the cache structure
	int cacheSize = numOfEdges - beginningEdgeIndex;
	cout << "cacheSize " << cacheSize << endl;
	vectorOfNodalUpdateMaps = VectorOfNodalUpdateMaps(cacheSize);

	// this corresponds to the first edge join the network at the observation time
	// vectorOfComingVertices[0]: the list of nodes join the network from the observation time (the first edge time) up to
	// the time of next edge.
	// vectorOfNodalUpdateMaps[0]: the list of nodes updated after entering the first edge to the network.
	// roll up to the end and cache data
	for (unsigned int i = beginningEdgeIndex, currentCacheIndex = 0; i
			< numOfEdges; i++, currentCacheIndex++) {
		unsigned int nodeB = edgeJoiningTimes[i].first.second;
		NodalUpdates myNodalUpdates = NodalUpdates();
		myNodalUpdates.push_back(StatValue(0, currentNodalNetworkStatistics(
				nodeB, 0)));
		vectorOfNodalUpdateMaps(currentCacheIndex) = NodalUpdateMap();
		vectorOfNodalUpdateMaps(currentCacheIndex).insert(make_pair(nodeB,
				myNodalUpdates));
	}

	// reset everything for the next use
	dirty = true;
	finish();
}

PEgoCentricNetworkData::~PEgoCentricNetworkData() {
	// TODO Auto-generated destructor stub
}

// roll a new edge and update nodal network statistics using cache structures
void PEgoCentricNetworkData::updateNodalNetworkStatistics(
		const EdgeTime& edgeTime) {
	unsigned int nodeA = edgeTime.first.first;
	unsigned int nodeB = edgeTime.first.second;
	currentRecentNodeContactingTimes[nodeA] = edgeTime.second;
	currentRecentNodeContactedTimes[nodeB] = edgeTime.second;
	currentNodalNetworkStatistics(nodeB, 0) += 1;
	currentCacheIndex++;
}

int PEgoCentricNetworkData::getFirstOrderPopularityIndex() {
	return 0;
}

}
