/*
 * P2PTRTiedEgoCentricNetworkData.cpp
 *
 *  Created on: Jul 23, 2010
 *      Author: duyvu
 */

#include "P2PTRTiedEgoCentricNetworkData.h"

namespace ndip {

P2PTRTiedEgoCentricNetworkData::P2PTRTiedEgoCentricNetworkData() {
	// TODO Auto-generated constructor stub

}

P2PTRTiedEgoCentricNetworkData::P2PTRTiedEgoCentricNetworkData(
		const ExposureTimeVectorType& _nodeJoiningTimes,
		const VectorOfEdgeEvents& _vectorOfEdgeEvents,
		double _observationTimeStart, double _observationTimeEnd) :
	TiedEgoCentricNetworkData(_nodeJoiningTimes, _vectorOfEdgeEvents,
			_observationTimeStart, _observationTimeEnd) {

	cout << "START the constructor of P2PTRTiedEgoCentricNetworkData" << endl;

	// Revise this hard code later
	windowSizeOfRenewalStatistic = WINDOW_SIZE_OF_RENEWAL_STATISTIC;

	// Initialize the current nodal network statistics
	numOfNodalNetworkStatistics
			= P2PTRTiedEgoCentricNetworkData::NUMBER_OF_NODAL_NETWORK_STATISTICS;
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
	ListOfEdgeEvents listOfActiveEvents = ListOfEdgeEvents();
	for (unsigned int e = 0; e < beginningEdgeIndex; e++) {
		EdgeEvents edgeEvents = vectorOfEdgeEvents[e];
		double nextCitingTime =
				nodeJoiningTimes[vectorOfEdgeEvents[e + 1].first];
		this->updateNodalNetworkStatisticsCache(edgeEvents, nextCitingTime,
				false, listOfActiveEvents);
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
		double nextCitingTime;
		if (e < endingEdgeIndex)
			nextCitingTime = nodeJoiningTimes[vectorOfEdgeEvents[e + 1].first];
		else
			// not really important since we do not use the covariate matrix after the last event.
			nextCitingTime = nodeJoiningTimes[vectorOfEdgeEvents[e].first] + .1;
		this->updateNodalNetworkStatisticsCache(edgeEvents, nextCitingTime,
				true, listOfActiveEvents);
	}

	// reset everything for the next use
	dirty = true;
	finish();

	cout << "FINISH the constructor of P2PTRTiedEgoCentricNetworkData" << endl;
}

P2PTRTiedEgoCentricNetworkData::~P2PTRTiedEgoCentricNetworkData() {
	// TODO Auto-generated destructor stub
}

void P2PTRTiedEgoCentricNetworkData::updateNodalNetworkStatisticsCache(
		const EdgeEvents& edgeEvents, double nextCitingTime, bool isCached,
		ListOfEdgeEvents& listOfActiveEvents) {

	std::set<unsigned int> setOf2PathUpdatedNodes;
	std::set<unsigned int> setOfBrokeeUpdatedNodes;
	std::set<unsigned int> setOfBrokerUpdatedNodes;
	bool willUpdateCitingNodeBlockStatistic = false;
	std::set<unsigned int> setOfRenewalStatisticUpdatedNodes;
	if (isCached) {
		vectorOfNodalUpdateMaps(currentCacheIndex) = NodalUpdateMap();
		setOf2PathUpdatedNodes = std::set<unsigned int>();
		setOfBrokeeUpdatedNodes = std::set<unsigned int>();
		setOfBrokerUpdatedNodes = std::set<unsigned int>();
		setOfRenewalStatisticUpdatedNodes = std::set<unsigned int>();
	}

	unsigned int citingNode = edgeEvents.first;
	double citingTime = nodeJoiningTimes[citingNode];
	unsigned int numOfSecondOutPathsFromCitingNode = 0;
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
			numOfSecondOutPathsFromCitingNode++;
			// add citingNode A into the second order popularity of C
			currentNodalNetworkStatistics(nodeC, 1) += 1;
			// cache the update of nodeC
			if (isCached)
				setOf2PathUpdatedNodes.insert(nodeC);
		}

		// triangles
		for (unsigned int j = 0; j < i; j++) {
			unsigned int citedNodeB = edgeEvents.second[j];
			if (edgeExist(citedNode, citedNodeB)) {
				currentNodalNetworkStatistics(citedNode, 3) += 1;
				currentNodalNetworkStatistics(citedNodeB, 2) += 1;
				currentNodalNetworkStatistics(citingNode, 4) += 1;
				if (isCached) {
					setOfBrokerUpdatedNodes.insert(citedNode);
					setOfBrokeeUpdatedNodes.insert(citedNodeB);
					willUpdateCitingNodeBlockStatistic = true;
				}
			} else if (edgeExist(citedNodeB, citedNode)) {
				currentNodalNetworkStatistics(citedNode, 2) += 1;
				currentNodalNetworkStatistics(citedNodeB, 3) += 1;
				currentNodalNetworkStatistics(citingNode, 4) += 1;
				if (isCached) {
					setOfBrokeeUpdatedNodes.insert(citedNode);
					setOfBrokerUpdatedNodes.insert(citedNodeB);
					willUpdateCitingNodeBlockStatistic = true;
				}
			}
		}

		// update the renewal statistic
		currentNodalNetworkStatistics(citedNode, 7) += 1;
		if (isCached)
			setOfRenewalStatisticUpdatedNodes.insert(citedNode);

		// add the edge to the current network
		add_edge(citingNode, citedNode, currentGraph);
	}
	// push the new edge event to the list of active events
	listOfActiveEvents.push_back(edgeEvents);

	// update the first and second out-degree statistics
	currentNodalNetworkStatistics(citingNode, 5) += edgeEvents.second.size();
	currentNodalNetworkStatistics(citingNode, 6)
			+= numOfSecondOutPathsFromCitingNode;

	// update the renewal statistic
	// drop out-of-date citations by the time of the next event
	while (listOfActiveEvents.size() > 0) {
		EdgeEvents myEvents = listOfActiveEvents.front();
		unsigned int myCitingNode = myEvents.first;
		double myCitingTime = nodeJoiningTimes[myCitingNode];
		// please remove me if I am obsolete
		if (myCitingTime < (nextCitingTime - windowSizeOfRenewalStatistic)) {
			for (unsigned int j = 0; j < myEvents.second.size(); j++) {
				unsigned int myCitedNode = myEvents.second[j];
				currentNodalNetworkStatistics(myCitedNode, 7) -= 1;
				if (isCached) {
					setOfRenewalStatisticUpdatedNodes.insert(myCitedNode);
				}
			}
			listOfActiveEvents.pop_front();
		} else
			break;
	}

	if (isCached) {
		// update 2-path statistics
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
		// update brokee nodes
		for (std::set<unsigned int>::iterator setIterator =
				setOfBrokeeUpdatedNodes.begin(); setIterator
				!= setOfBrokeeUpdatedNodes.end(); setIterator++) {
			unsigned int nodeC = *setIterator;
			NodalUpdateMap::iterator iter = vectorOfNodalUpdateMaps(
					currentCacheIndex).find(nodeC);
			if (iter != vectorOfNodalUpdateMaps(currentCacheIndex).end()) {
				NodalUpdates myNodalUpdates = ((NodalUpdates) iter->second);
				myNodalUpdates.push_back(StatValue(2,
						currentNodalNetworkStatistics(nodeC, 2)));
				vectorOfNodalUpdateMaps(currentCacheIndex)[nodeC]
						= myNodalUpdates;
			} else {
				NodalUpdates myNodalUpdates = NodalUpdates();
				myNodalUpdates.push_back(StatValue(2,
						currentNodalNetworkStatistics(nodeC, 2)));
				vectorOfNodalUpdateMaps(currentCacheIndex).insert(make_pair(
						nodeC, myNodalUpdates));
			}
		}
		// update broker nodes
		for (std::set<unsigned int>::iterator setIterator =
				setOfBrokerUpdatedNodes.begin(); setIterator
				!= setOfBrokerUpdatedNodes.end(); setIterator++) {
			unsigned int nodeC = *setIterator;
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
				vectorOfNodalUpdateMaps(currentCacheIndex).insert(make_pair(
						nodeC, myNodalUpdates));
			}
		}
		// update block node
		if (willUpdateCitingNodeBlockStatistic) {
			NodalUpdates myNodalUpdates = NodalUpdates();
			myNodalUpdates.push_back(StatValue(4,
					currentNodalNetworkStatistics(citingNode, 4)));
			vectorOfNodalUpdateMaps(currentCacheIndex).insert(make_pair(
					citingNode, myNodalUpdates));
		}
		// update out-degree
		NodalUpdateMap::iterator iter = vectorOfNodalUpdateMaps(
				currentCacheIndex).find(citingNode);
		if (iter != vectorOfNodalUpdateMaps(currentCacheIndex).end()) {
			NodalUpdates myNodalUpdates = ((NodalUpdates) iter->second);
			myNodalUpdates.push_back(StatValue(5,
					currentNodalNetworkStatistics(citingNode, 5)));
			myNodalUpdates.push_back(StatValue(6,
					currentNodalNetworkStatistics(citingNode, 6)));
			vectorOfNodalUpdateMaps(currentCacheIndex)[citingNode]
					= myNodalUpdates;
		} else {
			NodalUpdates myNodalUpdates = NodalUpdates();
			myNodalUpdates.push_back(StatValue(5,
					currentNodalNetworkStatistics(citingNode, 5)));
			myNodalUpdates.push_back(StatValue(6,
					currentNodalNetworkStatistics(citingNode, 6)));
			vectorOfNodalUpdateMaps(currentCacheIndex).insert(make_pair(
					citingNode, myNodalUpdates));
		}
		// update renewal nodes
		for (std::set<unsigned int>::iterator setIterator =
				setOfRenewalStatisticUpdatedNodes.begin(); setIterator
				!= setOfRenewalStatisticUpdatedNodes.end(); setIterator++) {
			unsigned int nodeC = *setIterator;
			NodalUpdateMap::iterator iter = vectorOfNodalUpdateMaps(
					currentCacheIndex).find(nodeC);
			if (iter != vectorOfNodalUpdateMaps(currentCacheIndex).end()) {
				NodalUpdates myNodalUpdates = ((NodalUpdates) iter->second);
				myNodalUpdates.push_back(StatValue(7,
						currentNodalNetworkStatistics(nodeC, 7)));
				vectorOfNodalUpdateMaps(currentCacheIndex)[nodeC]
						= myNodalUpdates;
			} else {
				NodalUpdates myNodalUpdates = NodalUpdates();
				myNodalUpdates.push_back(StatValue(7,
						currentNodalNetworkStatistics(nodeC, 7)));
				vectorOfNodalUpdateMaps(currentCacheIndex).insert(make_pair(
						nodeC, myNodalUpdates));
			}
		}
	}

	// for the next vectorOfNodalUpdateMaps
	currentCacheIndex++;
}

void P2PTRTiedEgoCentricNetworkData::updateNodalNetworkStatistics(
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

int P2PTRTiedEgoCentricNetworkData::getFirstOrderPopularityIndex() {
	return 0;
}

}
