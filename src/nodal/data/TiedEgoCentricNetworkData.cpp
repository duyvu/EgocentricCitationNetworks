/*
 * TEgoCentricNetworkData.cpp
 *
 *  Created on: Jul 18, 2010
 *      Author: duyvu
 */

#include "TiedEgoCentricNetworkData.h"

using namespace std;

namespace ndip {

TiedEgoCentricNetworkData::TiedEgoCentricNetworkData() {
	// TODO Auto-generated constructor stub

}

TiedEgoCentricNetworkData::TiedEgoCentricNetworkData(
		const ExposureTimeVectorType& _nodeJoiningTimes,
		const VectorOfEdgeEvents& _vectorOfEdgeEvents,
		double _observationTimeStart, double _observationTimeEnd) {

	cout << "START the constructor of TiedEgoCentricNetworkData" << endl;

	// set network data structures
	nodeJoiningTimes = _nodeJoiningTimes;
	numOfVertices = nodeJoiningTimes.size();
	vectorOfEdgeEvents = _vectorOfEdgeEvents;
	numOfDistinctEdgeEvents = vectorOfEdgeEvents.size();
	observationTimeStart = _observationTimeStart;
	observationTimeEnd = _observationTimeEnd;

	// set the index of the first edge event starting from observation time
	beginningEdgeIndex = 0;
	for (; beginningEdgeIndex < vectorOfEdgeEvents.size(); beginningEdgeIndex++) {
		if (nodeJoiningTimes[vectorOfEdgeEvents[beginningEdgeIndex].first]
				>= observationTimeStart)
			break;
	}
	//cout << "beginningEdgeIndex = " << beginningEdgeIndex << endl;

	// set the index of the last edge event at the end of observation process
	for (endingEdgeIndex = beginningEdgeIndex; endingEdgeIndex
			< vectorOfEdgeEvents.size(); endingEdgeIndex++) {
		if (nodeJoiningTimes[vectorOfEdgeEvents[endingEdgeIndex].first]
				> observationTimeEnd)
			break;
	}
	endingEdgeIndex--;
	//cout << "endingEdgeIndex = " << endingEdgeIndex << endl;

	// set the lists of vertices coming in [t_e, t_{e+1}) for events
	// after the observation time.
	vectorOfComingVertices = VectorOfListsOfVertices(endingEdgeIndex
			- beginningEdgeIndex + 1);
	double currentEventTime = observationTimeStart;
	unsigned int currentVertex = 0;
	unsigned int i = beginningEdgeIndex;
	for (; i <= endingEdgeIndex && currentVertex < numOfVertices; i++) {

		vectorOfComingVertices(i - beginningEdgeIndex) = ListOfVertices();
		double nextEventTime =
				nodeJoiningTimes[vectorOfEdgeEvents[i + 1].first];

		// Bypass nodes that join before the previous edge time, i.e. including this time
		while (currentVertex < numOfVertices && nodeJoiningTimes[currentVertex]
				< currentEventTime) {
			currentVertex++;
		}

		// Record nodes that join from the previous time to the time of the next edge
		while (currentVertex < numOfVertices && nodeJoiningTimes[currentVertex]
				< nextEventTime) {
			vectorOfComingVertices(i - beginningEdgeIndex).push_back(
					currentVertex);
			//cout << currentVertex << " ; ";
			currentVertex++;
		}

		//cout << endl;
		currentEventTime = nextEventTime;
	}

	// empty for the rest
	for (; i <= endingEdgeIndex; i++)
		vectorOfComingVertices(i - beginningEdgeIndex) = ListOfVertices();

	cout << "FINISH the constructor of TiedEgoCentricNetworkData" << endl;
}

TiedEgoCentricNetworkData::~TiedEgoCentricNetworkData() {
	// TODO Auto-generated destructor stub
}

// start rolling the network data edge by edge
void TiedEgoCentricNetworkData::start() {
	if (dirty) {
		// restart the nodal network statistics
		//currentNodalNetworkStatistics = SparseUnsignedIntMatrix(
		//		beginningNodalNetworkStatistics);
		currentNodalNetworkStatistics = beginningNodalNetworkStatistics;
		// restart the graph
		//currentGraph = Graph(beginningGraph);
		currentGraph = beginningGraph;
		// restart the recent node contact times
		//currentRecentNodeCitedTimes = ExposureTimeVectorType(
		//		beginningRecentNodeCitedTimes);
		currentRecentNodeCitedTimes = beginningRecentNodeCitedTimes;
		// restart cache structure which is used to update nodal network statistics rather
		// recompute them
		currentCacheIndex = 0;
		// the data container is clean now and can start the rolling update
		dirty = false;
	}
}

// finish rolling the network data edge by edge
void TiedEgoCentricNetworkData::finish() {
	dirty = true;
}

void TiedEgoCentricNetworkData::summaryCache() {
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
