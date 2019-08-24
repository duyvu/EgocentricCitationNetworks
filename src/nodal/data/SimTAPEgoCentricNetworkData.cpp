/*
 * SimTAPEgoCentricNetworkData.cpp
 *
 *  Created on: Aug 10, 2010
 *      Author: duyvu
 */

#include "SimTAPEgoCentricNetworkData.h"

namespace ndip {

SimTAPEgoCentricNetworkData::SimTAPEgoCentricNetworkData() {
	// TODO Auto-generated constructor stub

}

SimTAPEgoCentricNetworkData::SimTAPEgoCentricNetworkData(
		const ExposureTimeVectorType& _exposureTimes,
		const EdgeTimeVectorType& _edgeTimes, double _observationTime) :

	EgoCentricNetworkData(_exposureTimes, _edgeTimes, _observationTime) {
	// Initialize the current nodal network statistics
	numOfNodalNetworkStatistics
			= SimTAPEgoCentricNetworkData::NUMBER_OF_NODAL_NETWORK_STATISTICS;
	currentNodalNetworkStatistics = SparseUnsignedIntMatrix(numOfVertices,
			numOfNodalNetworkStatistics);
	currentNodalNetworkStatistics.clear();

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
		this->updateNetworkStatisticsWithoutCache(edgeTime);
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

	// reset everything for the next use
	dirty = true;
	finish();
}

SimTAPEgoCentricNetworkData::~SimTAPEgoCentricNetworkData() {
	// TODO Auto-generated destructor stub
}

void SimTAPEgoCentricNetworkData::updateNetworkStatisticsWithoutCache(
		const EdgeTime& edgeTime) {
	unsigned int nodeA = edgeTime.first.first;
	unsigned int nodeB = edgeTime.first.second;
	currentRecentNodeContactingTimes[nodeA] = edgeTime.second;
	currentRecentNodeContactedTimes[nodeB] = edgeTime.second;

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
		}
	}

	// update the successful transitive activity count
	// this statistic measures how many times a node receives contacts through transitive activity.
	// Is it more likely that an actor will be contacted more if she is successful in this kind of
	// configuration?
	// Different from the transitive triangle count, this statistic takes into account the order of
	// edge formation as well as the pivotal node (the node B below).
	// A -> B to B -> C and C -> A => B gets one more successful transitive activity.
	for (tie(in_A, in_A_end) = in_edges(nodeA, currentGraph); in_A != in_A_end; ++in_A) {
		unsigned int nodeC = source(*in_A, currentGraph);
		if (nodeC == nodeA || nodeC == nodeB)
			continue;
		for (tie(out_B, out_B_end) = out_edges(nodeB, currentGraph); out_B
				!= out_B_end; ++out_B) {
			if (nodeC == target(*out_B, currentGraph))
				currentNodalNetworkStatistics(nodeB, 4) += 1;
		}
	}

	// update the successful transitive popularity count
	// this statistic measures how many times a node receives contacts through transitive popularity.
	// Is it more likely that she will be contacted more if she is successful in this kind of
	// configuration?
	// Different from the transitive triangle count, this statistic takes into account the order of
	// edge formation and the pivotal node (the node B below).
	//A -> B to A -> C and C -> B => B gets one more successful transitive popularity.
	GraphTraits::in_edge_iterator in_B, in_B_end;
	GraphTraits::out_edge_iterator out_A, out_A_end;
	for (tie(out_A, out_A_end) = out_edges(nodeA, currentGraph); out_A
			!= out_A_end; ++out_A) {
		unsigned int nodeC = target(*out_A, currentGraph);
		if (nodeC == nodeA || nodeC == nodeB)
			continue;
		for (tie(in_B, in_B_end) = in_edges(nodeB, currentGraph); in_B
				!= in_B_end; ++in_B) {
			if (nodeC == source(*in_B, currentGraph))
				currentNodalNetworkStatistics(nodeB, 5) += 1;
		}
	}

	// add the edge to the current network
	add_edge(nodeA, nodeB, currentGraph);
}

void SimTAPEgoCentricNetworkData::updateNodalNetworkStatistics(
		const EdgeTime& edgeTime) {

	unsigned int nodeA = edgeTime.first.first;
	unsigned int nodeB = edgeTime.first.second;
	currentRecentNodeContactingTimes[nodeA] = edgeTime.second;
	currentRecentNodeContactedTimes[nodeB] = edgeTime.second;

	// update nodal network statistics
	// and add the edge to the current network. This step is NECCESSARY because
	// network statistics are not cached in this implementation
	updateNetworkStatisticsWithoutCache(edgeTime);
}

int SimTAPEgoCentricNetworkData::getFirstOrderPopularityIndex() {
	return 2;
}

}
