/*
 * SimP2PTRTiedEgoCentricNetworkData.cpp
 *
 *  Created on: Aug 14, 2010
 *      Author: duyvu
 */

#include "SimP2PTRTiedEgoCentricNetworkData.h"

namespace ndip {

SimP2PTRTiedEgoCentricNetworkData::SimP2PTRTiedEgoCentricNetworkData() {
	// TODO Auto-generated constructor stub

}

SimP2PTRTiedEgoCentricNetworkData::SimP2PTRTiedEgoCentricNetworkData(
		const ExposureTimeVectorType& _nodeJoiningTimes,
		const VectorOfEdgeEvents& _vectorOfEdgeEvents,
		double _observationTimeStart, double _observationTimeEnd) :
	TiedEgoCentricNetworkData(_nodeJoiningTimes, _vectorOfEdgeEvents,
			_observationTimeStart, _observationTimeEnd) {

	cout << "START the constructor of SimP2PTRTiedEgoCentricNetworkData"
			<< endl;

	// Revise this hard code later
	windowSizeOfRenewalStatistic = WINDOW_SIZE_OF_RENEWAL_STATISTIC;

	// Initialize the current nodal network statistics
	numOfNodalNetworkStatistics
			= SimP2PTRTiedEgoCentricNetworkData::NUMBER_OF_NODAL_NETWORK_STATISTICS;
	currentNodalNetworkStatistics = SparseUnsignedIntMatrix(numOfVertices,
			numOfNodalNetworkStatistics);
	currentNodalNetworkStatistics.clear();

	// Initialize the current graph
	currentGraph = Graph(numOfVertices);

	// Initialize the recent node contacted times
	currentRecentNodeCitedTimes = ExposureTimeVectorType(numOfVertices);
	for (unsigned int i = 0; i < numOfVertices; i++)
		currentRecentNodeCitedTimes[i] = nodeJoiningTimes[i];

	listOfActiveEvents = ListOfEdgeEvents();
	for (unsigned int e = 0; e < beginningEdgeIndex; e++) {
		EdgeEvents edgeEvents = vectorOfEdgeEvents[e];
		this->updateNetworkStatisticsWithoutCache(edgeEvents);
	}

	// Copy these current nodal network statistics and the graph up to the observation time
	// to save the time for iterative calls later
	beginningNodalNetworkStatistics = SparseUnsignedIntMatrix(
			currentNodalNetworkStatistics);
	beginningGraph = Graph(currentGraph);
	beginningRecentNodeCitedTimes = ExposureTimeVectorType(
			currentRecentNodeCitedTimes);

	// reset everything for the next use
	dirty = true;
	finish();
}

SimP2PTRTiedEgoCentricNetworkData::~SimP2PTRTiedEgoCentricNetworkData() {
	// TODO Auto-generated destructor stub
}

void SimP2PTRTiedEgoCentricNetworkData::updateNetworkStatisticsWithoutCache(
		const EdgeEvents& edgeEvents) {

	unsigned int citingNode = edgeEvents.first;
	double citingTime = nodeJoiningTimes[citingNode];
	unsigned int numOfSecondOutPathsFromCitingNode = 0;
	for (unsigned int i = 0; i < edgeEvents.second.size(); i++) {
		unsigned int citedNode = edgeEvents.second[i];

		// update current cited time
		currentRecentNodeCitedTimes[citedNode] = citingTime;

		// update the first in-degree statistic
		currentNodalNetworkStatistics(citedNode, 0) += 1;

		// citingNode A
		// citedNode B
		// A -> B, B -> C then add A into the second in-degree statistic of C
		GraphTraits::out_edge_iterator out_B, out_B_end;
		for (tie(out_B, out_B_end) = out_edges(citedNode, currentGraph); out_B
				!= out_B_end; ++out_B) {
			unsigned int nodeC = target(*out_B, currentGraph);
			numOfSecondOutPathsFromCitingNode++;
			// add citingNode A into the second order popularity of C
			currentNodalNetworkStatistics(nodeC, 1) += 1;
			// cache the update of nodeC
		}

		// triangles
		for (unsigned int j = 0; j < i; j++) {
			unsigned int citedNodeB = edgeEvents.second[j];
			if (edgeExist(citedNode, citedNodeB)) {
				currentNodalNetworkStatistics(citedNode, 3) += 1;
				currentNodalNetworkStatistics(citedNodeB, 2) += 1;
				currentNodalNetworkStatistics(citingNode, 4) += 1;
			} else if (edgeExist(citedNodeB, citedNode)) {
				currentNodalNetworkStatistics(citedNode, 2) += 1;
				currentNodalNetworkStatistics(citedNodeB, 3) += 1;
				currentNodalNetworkStatistics(citingNode, 4) += 1;
			}
		}

		// update the renewal statistic
		currentNodalNetworkStatistics(citedNode, 7) += 1;

		// add the edge to the current network
		add_edge(citingNode, citedNode, currentGraph);
	}

	// update the first and second out-degree statistics
	currentNodalNetworkStatistics(citingNode, 5) += edgeEvents.second.size();
	currentNodalNetworkStatistics(citingNode, 6)
			+= numOfSecondOutPathsFromCitingNode;

	while (listOfActiveEvents.size() > 0) {
		EdgeEvents myEvents = listOfActiveEvents.front();
		unsigned int myCitingNode = myEvents.first;
		double myCitingTime = nodeJoiningTimes[myCitingNode];
		if (myCitingTime < (citingTime - windowSizeOfRenewalStatistic)) {
			for (unsigned int j = 0; j < myEvents.second.size(); j++) {
				unsigned int myCitedNode = myEvents.second[j];
				currentNodalNetworkStatistics(myCitedNode, 7) -= 1;
			}
			listOfActiveEvents.pop_front();
		} else
			break;
	}
	// push the new edge event to the list of active events
	listOfActiveEvents.push_back(edgeEvents);
}

void SimP2PTRTiedEgoCentricNetworkData::updateNodalNetworkStatistics(
		const EdgeEvents& edgeEvents) {
	updateNetworkStatisticsWithoutCache(edgeEvents);
}

}
