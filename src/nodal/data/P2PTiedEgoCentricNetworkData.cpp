/*
 * P2PTiedEgoCentricNetworkData.cpp
 *
 *  Created on: Jul 20, 2010
 *      Author: duyvu
 */

#include "P2PTiedEgoCentricNetworkData.h"

namespace ndip {

P2PTiedEgoCentricNetworkData::P2PTiedEgoCentricNetworkData() {
	// TODO Auto-generated constructor stub

}

P2PTiedEgoCentricNetworkData::P2PTiedEgoCentricNetworkData(
		const ExposureTimeVectorType& _nodeJoiningTimes,
		const VectorOfEdgeEvents& _vectorOfEdgeEvents,
		double _observationTimeStart, double _observationTimeEnd) :
	TiedEgoCentricNetworkData(_nodeJoiningTimes, _vectorOfEdgeEvents,
			_observationTimeStart, _observationTimeEnd) {

	cout << "START the constructor of P2PTiedEgoCentricNetworkData" << endl;

	// Initialize the current nodal network statistics
	numOfNodalNetworkStatistics
			= P2PTiedEgoCentricNetworkData::NUMBER_OF_NODAL_NETWORK_STATISTICS;
	currentNodalNetworkStatistics = SparseUnsignedIntMatrix(numOfVertices,
			numOfNodalNetworkStatistics);
	currentNodalNetworkStatistics.clear();

	// Initialize the current graph
	currentGraph = Graph(numOfVertices);

	// Initialize the recent node contacted times
	currentRecentNodeCitedTimes = ExposureTimeVectorType(numOfVertices);
	for (unsigned int i = 0; i < numOfVertices; i++)
		currentRecentNodeCitedTimes[i] = nodeJoiningTimes[i];

	// Roll the current nodal network statistics
	// and the graph up to the observation time
	for (unsigned int e = 0; e < beginningEdgeIndex; e++) {
		EdgeEvents edgeEvents = vectorOfEdgeEvents[e];
		this->updateNodalNetworkStatisticsCache(edgeEvents, false);
	}

	// Copy these current nodal network statistics and the graph up to the observation time
	// to save the time for iterative calls later
	beginningNodalNetworkStatistics = SparseUnsignedIntMatrix(
			currentNodalNetworkStatistics);
	beginningGraph = Graph(currentGraph);
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
	currentCacheIndex = 0;
	// roll up to the end and cache data
	for (unsigned int e = beginningEdgeIndex; e <= endingEdgeIndex; e++) {
		EdgeEvents edgeEvents = vectorOfEdgeEvents[e];
		this->updateNodalNetworkStatisticsCache(edgeEvents, true);
	}

	// reset everything for the next use
	dirty = true;
	finish();

	cout << "FINISH the constructor of P2PTiedEgoCentricNetworkData" << endl;
}

P2PTiedEgoCentricNetworkData::~P2PTiedEgoCentricNetworkData() {
	// TODO Auto-generated destructor stub
}

void P2PTiedEgoCentricNetworkData::updateNodalNetworkStatisticsCache(
		const EdgeEvents& edgeEvents, bool isCached) {

	std::set<unsigned int> setOf2PathUpdatedNodes;
	if (isCached) {
		vectorOfNodalUpdateMaps(currentCacheIndex) = NodalUpdateMap();
		setOf2PathUpdatedNodes = std::set<unsigned int>();
	}

	unsigned int citingNode = edgeEvents.first;
	double citingTime = nodeJoiningTimes[citingNode];
	for (unsigned int i = 0; i < edgeEvents.second.size(); i++) {
		unsigned int citedNode = edgeEvents.second[i];

		// update current cited time
		currentRecentNodeCitedTimes[citedNode] = citingTime;

		// update the first in-degree statistic
		currentNodalNetworkStatistics(citedNode, 0) += 1;
		if (isCached) {
			NodalUpdates myNodalUpdates = NodalUpdates();
			myNodalUpdates.push_back(StatValue(0,
					currentNodalNetworkStatistics(citedNode, 0)));
			vectorOfNodalUpdateMaps(currentCacheIndex).insert(make_pair(
					citedNode, myNodalUpdates));
		}

		// citingNode A
		// citedNode B
		// A -> B, B -> C then add A into the second in-degree statistic of C
		GraphTraits::out_edge_iterator out_B, out_B_end;
		for (tie(out_B, out_B_end) = out_edges(citedNode, currentGraph); out_B
				!= out_B_end; ++out_B) {
			unsigned int nodeC = target(*out_B, currentGraph);
			// add citingNode A into the second order popularity of C
			currentNodalNetworkStatistics(nodeC, 1) += 1;
			// cache the update of nodeC
			if (isCached)
				setOf2PathUpdatedNodes.insert(nodeC);
		}

		// add the edge to the current network
		add_edge(citingNode, citedNode, currentGraph);
	}

	if (isCached) {
		for (std::set<unsigned int>::iterator setIterator =
				setOf2PathUpdatedNodes.begin(); setIterator
				!= setOf2PathUpdatedNodes.end(); setIterator++) {
			unsigned int nodeC = *setIterator;
			NodalUpdateMap::iterator iter = vectorOfNodalUpdateMaps(
					currentCacheIndex).find(nodeC);
			if (iter != vectorOfNodalUpdateMaps(currentCacheIndex).end()) {
				NodalUpdates myNodalUpdates = ((NodalUpdates) iter->second);
				myNodalUpdates.push_back(StatValue(1,
						currentNodalNetworkStatistics(nodeC, 1)));
				vectorOfNodalUpdateMaps(currentCacheIndex)[nodeC]
						= myNodalUpdates;
			} else {
				NodalUpdates myNodalUpdates = NodalUpdates();
				myNodalUpdates.push_back(StatValue(1,
						currentNodalNetworkStatistics(nodeC, 1)));
				vectorOfNodalUpdateMaps(currentCacheIndex).insert(make_pair(
						nodeC, myNodalUpdates));
			}
		}
	}

	// for the next vectorOfNodalUpdateMaps
	currentCacheIndex++;
}

void P2PTiedEgoCentricNetworkData::updateNodalNetworkStatistics(
		const EdgeEvents& edgeEvents) {

	unsigned int citingNode = edgeEvents.first;
	double citingTime = nodeJoiningTimes[citingNode];
	for (unsigned int i = 0; i < edgeEvents.second.size(); i++) {
		unsigned int citedNode = edgeEvents.second[i];
		currentRecentNodeCitedTimes[citedNode] = citingTime;
	}

	//unsigned int updatedNodeCount = 0;
	//unsigned int updatedStatisticsCount = 0;
	for (NodalUpdateMap::iterator it = vectorOfNodalUpdateMaps(
			currentCacheIndex).begin(); it != vectorOfNodalUpdateMaps(
			currentCacheIndex).end(); ++it) {
		unsigned int updatedNode = it->first;
		//updatedNodeCount++;
		//cout << updatedNodeCount << " updating the node " << myNode << ": ";
		NodalUpdates myNodeUpdates = ((NodalUpdates) it->second);
		for (NodalUpdates::iterator it1 = myNodeUpdates.begin(); it1
				!= myNodeUpdates.end(); it1++) {
			StatValue value = *it1;
			currentNodalNetworkStatistics(updatedNode, value.first)
					= value.second;
			//updatedStatisticsCount++;
			//cout << "(" << value.first << ", " << value.second << ") ; ";
		}
		//cout << endl;
	}

	//cout << "The number of updated nodes is " << updatedNodeCount << endl;
	//cout << "The number of updated statistics is " << updatedStatisticsCount << endl;

	// add the edge to the current network; however, if it is NOT necessary
	// because we already cached all statistics
	// add_edge(nodeA, nodeB, currentGraph);

	// for the next VectorOfUpdateElements
	currentCacheIndex++;
}

int P2PTiedEgoCentricNetworkData::getFirstOrderPopularityIndex() {
	return 0;
}

}
