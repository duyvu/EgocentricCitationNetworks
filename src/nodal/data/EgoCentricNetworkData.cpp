/*
 * EgoCentricNetworkData.cpp
 *
 *  Created on: May 26, 2010
 *      Author: duyvu
 */

#include "EgoCentricNetworkData.h"

namespace ndip {

EgoCentricNetworkData::EgoCentricNetworkData() {
	// TODO Auto-generated constructor stub
}

EgoCentricNetworkData::EgoCentricNetworkData(
		const ExposureTimeVectorType& _nodeTimes,
		const EdgeTimeVectorType& _edgeTimes, double _observationTime) {

	// set network data structures
	nodeJoiningTimes = _nodeTimes;
	numOfVertices = nodeJoiningTimes.size();
	edgeJoiningTimes = _edgeTimes;
	numOfEdges = edgeJoiningTimes.size();
	observationTime = _observationTime;

	// set the index of the first edge starting from observation time
	beginningEdgeIndex = 0;
	for (; beginningEdgeIndex < edgeJoiningTimes.size(); beginningEdgeIndex++) {
		if (edgeJoiningTimes[beginningEdgeIndex].second >= observationTime)
			break;
	}

	// set the lists of vertices coming in [t_e, t_{e+1}) for events
	// after the observation time.
	// HERE WE ASSUME DATA OF BOTH NODES AND EDGES ARE INTO TIME ORDER
	vectorOfComingVertices = VectorOfListsOfVertices(numOfEdges
			- beginningEdgeIndex); // redundant by 2 elements
	double currentEdgeTime = observationTime;
	unsigned int currentVertex = 0;
	unsigned int i = beginningEdgeIndex;
	for (; i < numOfEdges && currentVertex < numOfVertices; i++) {

		cout << "i = " << i << endl;
		cout << "currentVertex = " << currentVertex << "; ";
		cout << "edge [" << edgeJoiningTimes[i].first.first << ", "
				<< edgeJoiningTimes[i].first.second << "]: ";

		vectorOfComingVertices(i - beginningEdgeIndex) = ListOfVertices();
		double nextEdgeTime = edgeJoiningTimes[i + 1].second;

		// Bypass nodes that join before the previous edge time, i.e. including this time
		while (currentVertex < numOfVertices && nodeJoiningTimes[currentVertex]
				< currentEdgeTime) {
			currentVertex++;
		}

		// Record nodes that join from the previous time to the time of the next edge
		while (currentVertex < numOfVertices && nodeJoiningTimes[currentVertex]
				< nextEdgeTime) {
			vectorOfComingVertices(i - beginningEdgeIndex).push_back(
					currentVertex);
			cout << currentVertex << " ; ";
			currentVertex++;
		}

		cout << endl;
		currentEdgeTime = nextEdgeTime;
	}

	// empty for the rest
	for (; i < numOfEdges; i++) {
		cout << "i = " << i << endl;
		vectorOfComingVertices(i - beginningEdgeIndex) = ListOfVertices();
	}
}

EgoCentricNetworkData::~EgoCentricNetworkData() {
	// TODO Auto-generated destructor stub
}

void EgoCentricNetworkData::start() {

#ifdef DEBUG_EGO_CENTRIC_NETWORK_DATA
	cout << "Start Rolling Network Data" << endl;
#endif

	if (dirty) {
		// restart the nodal network statistics
		currentNodalNetworkStatistics = SparseUnsignedIntMatrix(
				beginningNodalNetworkStatistics);
		// restart the graph
		currentGraph = Graph(beginningGraph);
		// restart the recent node contact times
		currentRecentNodeContactingTimes = ExposureTimeVectorType(
						beginningRecentNodeContactingTimes);
		currentRecentNodeContactedTimes = ExposureTimeVectorType(
				beginningRecentNodeContactedTimes);
		// restart cache structure which is used to update nodal network statistics rather
		// recompute them
		currentCacheIndex = 0;
		// the data container is clean now and can start the rolling update
		dirty = false;
	}

}

void EgoCentricNetworkData::finish() {

#ifdef DEBUG_EGO_CENTRIC_NETWORK_DATA
	cout << "Finish Rolling Network Data" << endl;
#endif

	dirty = true;

}

void EgoCentricNetworkData::summaryCache() {

	unsigned int cacheSize = vectorOfNodalUpdateMaps.size();
	cout << "The number of cached maps is " << cacheSize << endl;

	int numUpdatedNodes = 0;
	int numUpdatedStatistics = 0;
	for (unsigned int i = 0; i < cacheSize; i++) {
		numUpdatedNodes += vectorOfNodalUpdateMaps(i).size();

		for (NodalUpdateMap::iterator it = vectorOfNodalUpdateMaps(i).begin(); it
				!= vectorOfNodalUpdateMaps(i).end(); ++it) {

			NodalUpdates myNodeUpdates = ((NodalUpdates) it->second);
			for (NodalUpdates::iterator it1 = myNodeUpdates.begin(); it1
					!= myNodeUpdates.end(); it1++) {
				numUpdatedStatistics++;
			}
		}
	}

	cout << "The total number of updated nodes is " << numUpdatedNodes << endl;
	cout << "The average number of updated nodes per new edge is " << 1.0
			* numUpdatedNodes / cacheSize << endl;

	cout << "The total number of updated statistics is "
			<< numUpdatedStatistics << endl;
	cout << "The average number of updated statistics per new edge is " << 1.0
			* numUpdatedStatistics / cacheSize << endl;
}

}
