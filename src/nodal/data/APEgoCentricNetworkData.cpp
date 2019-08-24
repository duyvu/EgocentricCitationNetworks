/*
 * APEgoCentricNetworkData.cpp
 *
 *  Created on: May 26, 2010
 *      Author: duyvu
 */

#include "APEgoCentricNetworkData.h"

namespace ndip {

APEgoCentricNetworkData::APEgoCentricNetworkData() {
	// TODO Auto-generated constructor stub

}

APEgoCentricNetworkData::APEgoCentricNetworkData(
		const ExposureTimeVectorType& _exposureTimes,
		const EdgeTimeVectorType& _edgeTimes, double _observationTime) :
	EgoCentricNetworkData(_exposureTimes, _edgeTimes, _observationTime) {

	// Initialize the current nodal network statistics
	numOfNodalNetworkStatistics
			= APEgoCentricNetworkData::NUMBER_OF_NODAL_NETWORK_STATISTICS;
	currentNodalNetworkStatistics = SparseUnsignedIntMatrix(numOfVertices,
			numOfNodalNetworkStatistics);

	// Initialize the current graph
	currentGraph = Graph(numOfVertices);

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
		this->updateNodalNetworkStatisticsCache(edgeTime, false);
	}

	// Copy these current nodal network statistics and the graph up to the observation time
	// to save the time for iterative calls later
	beginningNodalNetworkStatistics = SparseUnsignedIntMatrix(
			currentNodalNetworkStatistics);
	beginningGraph = Graph(currentGraph);
	beginningRecentNodeContactingTimes = ExposureTimeVectorType(
				currentRecentNodeContactingTimes);
	beginningRecentNodeContactedTimes = ExposureTimeVectorType(
			currentRecentNodeContactedTimes);

	// run through the whole edge times to construct the cache
#ifdef DEBUG_AP_EGO_CENTRIC_NETWORK_DATA
	cout << "Caching..." << endl;
#endif

	// initialize the cache structure
	int cacheSize = numOfEdges - beginningEdgeIndex;
	cout << "cacheSize " << cacheSize << endl;
	vectorOfNodalUpdateMaps = VectorOfNodalUpdateMaps(cacheSize);

	// this corresponds to the first edge join the network at the observation time
	// vectorOfComingVertices[0]: the list of nodes join the network from the observation time (the first edge time) up to
	// the time of next edge.
	// vectorOfNodalUpdateMaps[0]: the list of nodes updated after entering the first edge to the network.
	currentCacheIndex = 0;

	// roll up to the end and cache data
	for (unsigned int e = beginningEdgeIndex; e < numOfEdges; e++) {
		//cout << "Caching the edge index " << e << endl;
		updateNodalNetworkStatisticsCache(edgeJoiningTimes[e], true);
	}

	// reset everything for the next use
	dirty = true;
	finish();
}

APEgoCentricNetworkData::~APEgoCentricNetworkData() {
	// TODO Auto-generated destructor stub
}

// prepare cache structures
void APEgoCentricNetworkData::updateNodalNetworkStatisticsCache(
		const EdgeTime& edgeTime, bool isCached) {

	unsigned int nodeA = edgeTime.first.first;
	unsigned int nodeB = edgeTime.first.second;
	currentRecentNodeContactingTimes[nodeA] = edgeTime.second;
	currentRecentNodeContactedTimes[nodeB] = edgeTime.second;

	if (isCached) {
		vectorOfNodalUpdateMaps(currentCacheIndex) = NodalUpdateMap();
	}

	/* update first order activity */
	currentNodalNetworkStatistics(nodeA, 0) += 1;
	/* update second order activity */
	// A -> B, B -> C then add C into A second order activity list
	currentNodalNetworkStatistics(nodeA, 1) += currentNodalNetworkStatistics(
			nodeB, 0);

	/* update first order popularity */
	currentNodalNetworkStatistics(nodeB, 2) += 1;
	/* update second order popularity */
	// C -> A, A -> B then add C into B second order popularity list
	currentNodalNetworkStatistics(nodeB, 3) += currentNodalNetworkStatistics(
			nodeA, 2);

	// C -> A, A -> B then add B into the second order activity of C
	GraphTraits::in_edge_iterator in_A, in_A_end;
	for (tie(in_A, in_A_end) = in_edges(nodeA, currentGraph); in_A != in_A_end; ++in_A) {
		unsigned int nodeC = source(*in_A, currentGraph);

		// The new A -> B IS a reciprocity edge
		if (nodeC == nodeB) {
			// over-count the second order popularity of B by 1
			currentNodalNetworkStatistics(nodeB, 3) -= 1;
			// do not count self second order activity
		} else {

			// The new A -> B is NOT a reciprocity edge

			// add B to C second order activity list
			currentNodalNetworkStatistics(nodeC, 1) += 1;

			// cache the update of nodeC
			if (isCached) {
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
					vectorOfNodalUpdateMaps(currentCacheIndex).insert(
							make_pair(nodeC, myNodalUpdates));
				}
			}
		}
	}

	// A -> B, B -> C then add A into the second order popularity of C
	GraphTraits::out_edge_iterator out_B, out_B_end;
	for (tie(out_B, out_B_end) = out_edges(nodeB, currentGraph); out_B
			!= out_B_end; ++out_B) {
		unsigned int nodeC = target(*out_B, currentGraph);

		// The new A -> B IS a reciprocity edge
		if (nodeC == nodeA) {
			// over-count the second order activity of A by 1
			currentNodalNetworkStatistics(nodeA, 1) -= 1;
			// do not count self second order popularity
		} else {

			// The new A -> B is NOT a reciprocity edge

			// add A into the second order popularity of C
			currentNodalNetworkStatistics(nodeC, 3) += 1;

			// cache the update of nodeC
			if (isCached) {
				NodalUpdateMap::iterator iter = vectorOfNodalUpdateMaps(
						currentCacheIndex).find(nodeC);
				if (iter != vectorOfNodalUpdateMaps(currentCacheIndex).end()) {
					NodalUpdates myNodalUpdates = ((NodalUpdates) iter->second);
					myNodalUpdates.push_back(StatValue(3,
							currentNodalNetworkStatistics(nodeC, 3)));
					vectorOfNodalUpdateMaps(currentCacheIndex)[nodeC]
							= myNodalUpdates;
				} else {
					NodalUpdates myNodalUpdates = NodalUpdates();
					myNodalUpdates.push_back(StatValue(3,
							currentNodalNetworkStatistics(nodeC, 3)));
					vectorOfNodalUpdateMaps(currentCacheIndex).insert(
							make_pair(nodeC, myNodalUpdates));
				}
			}
		}
	}

	// Finally cache nodeA and nodeB. Has to wait to the end to avoid the over-counting of
	// the second order activity of A and the second order popularity of B
	if (isCached) {
		// Cache nodeA
		NodalUpdateMap::iterator iter = vectorOfNodalUpdateMaps(
				currentCacheIndex).find(nodeA);
		if (iter != vectorOfNodalUpdateMaps(currentCacheIndex).end()) {
			NodalUpdates myNodalUpdates = ((NodalUpdates) iter->second);
			myNodalUpdates.push_back(StatValue(0,
					currentNodalNetworkStatistics(nodeA, 0)));
			myNodalUpdates.push_back(StatValue(1,
					currentNodalNetworkStatistics(nodeA, 1)));
			vectorOfNodalUpdateMaps(currentCacheIndex)[nodeA] = myNodalUpdates;
		} else {
			NodalUpdates myNodalUpdates = NodalUpdates();
			myNodalUpdates.push_back(StatValue(0,
					currentNodalNetworkStatistics(nodeA, 0)));
			myNodalUpdates.push_back(StatValue(1,
					currentNodalNetworkStatistics(nodeA, 1)));
			vectorOfNodalUpdateMaps(currentCacheIndex).insert(make_pair(nodeA,
					myNodalUpdates));
		}

		// Cache nodeB
		iter = vectorOfNodalUpdateMaps(currentCacheIndex).find(nodeB);
		if (iter != vectorOfNodalUpdateMaps(currentCacheIndex).end()) {
			NodalUpdates myNodalUpdates = ((NodalUpdates) iter->second);
			myNodalUpdates.push_back(StatValue(2,
					currentNodalNetworkStatistics(nodeB, 2)));
			myNodalUpdates.push_back(StatValue(3,
					currentNodalNetworkStatistics(nodeB, 3)));
			vectorOfNodalUpdateMaps(currentCacheIndex)[nodeB] = myNodalUpdates;
		} else {
			NodalUpdates myNodalUpdates = NodalUpdates();
			myNodalUpdates.push_back(StatValue(2,
					currentNodalNetworkStatistics(nodeB, 2)));
			myNodalUpdates.push_back(StatValue(3,
					currentNodalNetworkStatistics(nodeB, 3)));
			vectorOfNodalUpdateMaps(currentCacheIndex).insert(make_pair(nodeB,
					myNodalUpdates));
		}
	}

	// add the edge to the current network
	add_edge(nodeA, nodeB, currentGraph);

	// for the next vectorOfNodalUpdateMaps
	currentCacheIndex++;
}

// roll a new edge and update nodal network statistics using cache structures
void APEgoCentricNetworkData::updateNodalNetworkStatistics(
		const EdgeTime& edgeTime) {

	unsigned int nodeA = edgeTime.first.first;
	unsigned int nodeB = edgeTime.first.second;
	currentRecentNodeContactingTimes[nodeA] = edgeTime.second;
	currentRecentNodeContactedTimes[nodeB] = edgeTime.second;

	//unsigned int updatedNodeCount = 0;
	//unsigned int updatedStatisticsCount = 0;

	for (NodalUpdateMap::iterator it = vectorOfNodalUpdateMaps(
			currentCacheIndex).begin(); it != vectorOfNodalUpdateMaps(
			currentCacheIndex).end(); ++it) {

		unsigned int myNode = it->first;

		//updatedNodeCount++;
		//cout << updatedNodeCount << " updating the node " << myNode << ": ";

		NodalUpdates myNodeUpdates = ((NodalUpdates) it->second);
		for (NodalUpdates::iterator it1 = myNodeUpdates.begin(); it1
				!= myNodeUpdates.end(); it1++) {
			StatValue value = *it1;
			currentNodalNetworkStatistics(myNode, value.first) = value.second;

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

int APEgoCentricNetworkData::getFirstOrderPopularityIndex() {
	return 2;
}

}
