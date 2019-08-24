/*
 * EgoCentricNetworkData.h
 *
 *  Created on: May 26, 2010
 *      Author: duyvu
 */

#ifndef EGOCENTRICNETWORKDATA_H_
#define EGOCENTRICNETWORKDATA_H_

#include "DataTypes.h"

namespace ndip {

class EgoCentricNetworkData {
protected:
	// the number of nodes
	unsigned int numOfVertices;
	// a time-ordered list of nodes and corresponding join times. It is used to compute the age time.
	ExposureTimeVectorType nodeJoiningTimes;

	// the number of edges
	unsigned int numOfEdges;
	// a time-ordered list of edges and corresponding join times
	EdgeTimeVectorType edgeJoiningTimes;

	// assume the observation time is the time of the first edge
	// recorded on the continuous-time scale. Otherwise, the first element of vectorOfComingVertices
	// is the list of vertices coming from the observation time to the first edge after this time.
	// The assumption guarantees the synchronization of the new node lists and the update lists.
	double observationTime;

	// The graph before the observation time
	Graph beginningGraph;

	// The nodal network statistics before the observation time
	SparseUnsignedIntMatrix beginningNodalNetworkStatistics;

	// a list of recent contacted times to nodes after the observation time
	// it will be used to extract the gap times
	ExposureTimeVectorType beginningRecentNodeContactedTimes;

	// a list of recent contacting times from nodes after the observation time
	// it will be used to extract the gap times
	ExposureTimeVectorType beginningRecentNodeContactingTimes;

	// Other nodal covariates can be added here as a full matrix.

	// This model is dirty, need to be restarted before use.
	bool dirty;

	// The current cache index. The 0 index corresponds to updates from [observationTimes, the next edge).
	// We will set observationTimes equal to the time of the last edge when we start collecting on continuous-time scale
	unsigned int currentCacheIndex;
public:

	// The number of network statistics
	unsigned int numOfNodalNetworkStatistics;

	// The first edge starting from the observation time
	unsigned int beginningEdgeIndex;

	// The network up to now
	Graph currentGraph;

	// The nodal network statistics up to now
	SparseUnsignedIntMatrix currentNodalNetworkStatistics;

	// a list of recent contacted times to nodes. It is used to compute the gap time.
	ExposureTimeVectorType currentRecentNodeContactedTimes;

	// a list of recent contacting times from nodes. It is used to compute the gap time.
	ExposureTimeVectorType currentRecentNodeContactingTimes;

	// These below cache structures are computed starting from the observation time.

	// Between times of two edges e and e+1, some new vertices join the network
	// Each element of this vector is a list of vertices coming in [t_e, t_{e+1}).
	VectorOfListsOfVertices vectorOfComingVertices;

	// Data structure for caching the network statistics update
	// Each element of this vector is a map of nodes updated after adding the edge t_e
	// Each node is mapped to a list of updated statistics <nodalNetworkStatIndex, newValue>.
	// IMPORTANT NODE for retrieve the update set:
	// The nodes updated after a new kth edge is added include nodes (keys) of the map which is
	// kth element of the vector.
	VectorOfNodalUpdateMaps vectorOfNodalUpdateMaps;

	EgoCentricNetworkData();
	EgoCentricNetworkData(const ExposureTimeVectorType& _exposureTimes,
			const EdgeTimeVectorType& _edgeTimes, double _observationTime);
	virtual ~EgoCentricNetworkData();

	// start rolling the network data edge by edge
	void start();

	// roll a new edge and update nodal network statistics using cache structures
	virtual void updateNodalNetworkStatistics(const EdgeTime& edgeTime) = 0;
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

#endif /* EGOCENTRICNETWORKDATA_H_ */
