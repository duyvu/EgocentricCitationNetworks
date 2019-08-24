/*
 * TEgoCentricNetworkData.h
 *
 *  Created on: Jul 18, 2010
 *      Author: duyvu
 */

#ifndef TEGOCENTRICNETWORKDATA_H_
#define TEGOCENTRICNETWORKDATA_H_

#include "DataTypes.h"

namespace ndip {

class TiedEgoCentricNetworkData {
protected:
	// the number of nodes
	unsigned int numOfVertices;
	// a time-ordered list of nodes and corresponding join times. It is used to compute the age time.
	ExposureTimeVectorType nodeJoiningTimes;

	// the number of distinct edge events
	unsigned int numOfDistinctEdgeEvents;
	// a time-ordered list of edges and corresponding join times
	VectorOfEdgeEvents vectorOfEdgeEvents;

	// assume the observation time is the time of the first edge event
	// recorded on the continuous-time scale. Otherwise, the first element of vectorOfComingVertices
	// is the list of vertices coming from the observation time to the first edge after this time.
	// The assumption guarantees the synchronization of the new node lists and the update lists.
	double observationTimeStart;
	// assume the observation time is the time of the last edge event
	double observationTimeEnd;

	// The graph before the observation time
	Graph beginningGraph;

	// The nodal network statistics before the observation time
	//SparseUnsignedIntMatrix beginningNodalNetworkStatistics;
	DoubleMatrix beginningNodalNetworkStatistics;

	// a list of recent contacted times to nodes after the observation time
	// it will be used to extract the gap times
	ExposureTimeVectorType beginningRecentNodeCitedTimes;

	// Other nodal covariates can be added here as a full matrix.

	// This model is dirty, need to be restarted before use.
	bool dirty;

	// The current cache index. The 0 index corresponds to updates from [observationTimes, the next edge).
	// We will set observationTimes equal to the time of the last edge when we start collecting on continuous-time scale
	unsigned int currentCacheIndex;
public:
	// The number of network statistics
	unsigned int numOfNodalNetworkStatistics;

	// The first edge event starting from the observation time
	unsigned int beginningEdgeIndex;
	// The last edge event at the end of the observation process
	unsigned int endingEdgeIndex;

	// The network up to now
	Graph currentGraph;

	// The nodal network statistics up to now
	//SparseUnsignedIntMatrix currentNodalNetworkStatistics;
	DoubleMatrix currentNodalNetworkStatistics;

	// a list of recent contacted times to nodes. It is used to compute the gap time.
	ExposureTimeVectorType currentRecentNodeCitedTimes;

	// These below cache structures are computed starting from the observation time.

	// Between times of two edge events e and e+1, some new vertices join the network
	// Each element of this vector is a list of vertices coming in [t_e, t_{e+1}).
	VectorOfListsOfVertices vectorOfComingVertices;

	// Data structure for caching the network statistics update
	// Each element of this vector is a map of nodes updated after adding the edge t_e
	// Each node is mapped to a list of updated statistics <nodalNetworkStatIndex, newValue>.
	// IMPORTANT NODE for retrieve the update set:
	// The nodes updated after a new kth edge is added include nodes (keys) of the map which is
	// kth element of the vector.
	VectorOfNodalUpdateMaps vectorOfNodalUpdateMaps;

	TiedEgoCentricNetworkData();
	TiedEgoCentricNetworkData(const ExposureTimeVectorType& _nodeJoiningTimes,
			const VectorOfEdgeEvents& _vectorOfEdgeEvents,
			double _observationTimeStart, double _observationTimeEnd);
	virtual ~TiedEgoCentricNetworkData();

	// start rolling the network data edge by edge
	void start();

	// roll a new edge and update nodal network statistics using cache structures
	virtual void updateNodalNetworkStatistics(const EdgeEvents& edgeEvents) = 0;
	virtual void updateNodalNetworkStatistics(unsigned int indexOfEdgeEvents) {
		updateNodalNetworkStatistics(vectorOfEdgeEvents[indexOfEdgeEvents]);
	}

	virtual int getFirstOrderPopularityIndex() = 0;

	// finish rolling the network data edge by edge
	void finish();

	void summaryCache();

	void getCurrentGraph(Graph& graph) {
		graph = currentGraph;
	}

	unsigned int getCurrentCacheIndex() {
		return currentCacheIndex;
	}

	bool edgeExist(unsigned int nodeA, unsigned int nodeB) {
		GraphTraits::out_edge_iterator out_A, out_A_end;
		for (tie(out_A, out_A_end) = out_edges(nodeA, currentGraph); out_A
				!= out_A_end; ++out_A) {
			unsigned int nodeC = target(*out_A, currentGraph);
			if (nodeC == nodeB)
				return true;
		}
		return false;
	}
};
}

#endif /* TEGOCENTRICNETWORKDATA_H_ */
