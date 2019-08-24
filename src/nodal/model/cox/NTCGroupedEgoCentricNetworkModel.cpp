/*
 * NTCGroupedEgoCentricNetworkModel.cpp
 *
 *  Created on: Jan 18, 2011
 *      Author: duyvu
 */

#include "NTCGroupedEgoCentricNetworkModel.h"

#include "util/ranker.h"

#include <iostream>
#include <fstream>

namespace ndip {

NTCGroupedEgoCentricNetworkModel::NTCGroupedEgoCentricNetworkModel() {
	// TODO Auto-generated constructor stub

}

NTCGroupedEgoCentricNetworkModel::~NTCGroupedEgoCentricNetworkModel() {
	// TODO Auto-generated destructor stub
}

double NTCGroupedEgoCentricNetworkModel::computeLogLikelihood(
		const DoubleVector& beta, bool isUsingLBF) {
	if (isUsingLBF) {
		return computeLogLikelihoodWithLBF(beta);
	} else {
		cout << "This feature has not been supported!" << endl;
	}
}

double NTCGroupedEgoCentricNetworkModel::computeLogLikelihoodWithLBF(
		const DoubleVector& beta) {

	//cout << "computing LogLikelihood with LBF..." << endl;

	double logLLH = 0;

	// start rolling edges
	data->start();

	// the term that will be forwarded between two iterations
	double kappa = 0;

	unsigned int e = data->beginningEdgeIndex;

	// START - The first event

	// add those cited nodes
	DiscreteTimeVectorOfEdgeEvents edgeEvents =
			vectorOfDiscreteTimeVectorOfEdgeEvents[e];
	double currentCitingTime = edgeEvents.first;
	unsigned int numOfTiedEvents = 0;
	for (unsigned int k = 0; k < edgeEvents.second.size(); k++) {
		numOfTiedEvents += edgeEvents.second[k].second.size();
		for (unsigned int i = 0; i < edgeEvents.second[k].second.size(); i++) {
			unsigned int citedNode = edgeEvents.second[k].second[i];
			logLLH += computeSumNodalStatistics(citedNode, beta);
		}
	}

	// compute kappa(1)
	for (unsigned int i = 0; i < numOfVertices; i++) {
		// only consider those nodes are already in the network
		if ((nodeJoiningTimes[i] < currentCitingTime)) {
			// sum over statistics of the node i
			double nodalStatSum = computeSumNodalStatistics(i, beta);
			kappa += exp(nodalStatSum);
		}
	}

	//cout << "kappa(" << e << "): " << kappa << endl;
	logLLH -= numOfTiedEvents * log(kappa);

	// compute the lookahead term
	double lookaheadTerm = computeKappaOfNodalUpdateMap(
			data->getCurrentCacheIndex(), beta);

	// Update statistics and add the edge to the current graph
	data->updateNodalNetworkStatistics(edgeEvents);
	e++;

	// END - kappa(1)

	// START - kappa(k) for k = 2,...,M
	for (; e <= data->endingEdgeIndex; e++) {
		// add those cited nodes
		DiscreteTimeVectorOfEdgeEvents edgeEvents =
				vectorOfDiscreteTimeVectorOfEdgeEvents[e];
		double currentCitingTime = edgeEvents.first;
		unsigned int numOfTiedEvents = 0;
		for (unsigned int k = 0; k < edgeEvents.second.size(); k++) {
			numOfTiedEvents += edgeEvents.second[k].second.size();
			for (unsigned int i = 0; i < edgeEvents.second[k].second.size(); i++) {
				unsigned int citedNode = edgeEvents.second[k].second[i];
				logLLH += computeSumNodalStatistics(citedNode, beta);
			}
		}

		// Here we update kappa by adding the difference
		double forwardTerm = 0;
		forwardTerm += computeKappaOfNodalUpdateMap(
				data->getCurrentCacheIndex() - 1, beta);
		forwardTerm += computeKappaOfComingVertices(
				data->getCurrentCacheIndex() - 1, beta);
		kappa = kappa - lookaheadTerm + forwardTerm;

		//cout << "kappa(" << e << "): " << kappa << endl;
		logLLH -= numOfTiedEvents * log(kappa);

		lookaheadTerm = computeKappaOfNodalUpdateMap(
				data->getCurrentCacheIndex(), beta);

		// Update statistics and add the edge to the current graph
		data->updateNodalNetworkStatistics(edgeEvents);
	}
	// END - kappa(k)

	// finish rolling edges
	data->finish();
	return logLLH;
}

double NTCGroupedEgoCentricNetworkModel::computeSumNodalStatistics(
		unsigned int nodeID, const DoubleVector& beta) {
	double nodalStatSum = 0;
	for (unsigned int statIdx = 0; statIdx < data->numOfNodalNetworkStatistics; statIdx++)
		nodalStatSum += beta(statIdx) * data->currentNodalNetworkStatistics(
				nodeID, statIdx);
	return nodalStatSum;
}

double NTCGroupedEgoCentricNetworkModel::computeKappaOfNodalUpdateMap(
		int cacheIndex, const DoubleVector& beta) {
	double result = 0;
	for (NodalUpdateMap::iterator it =
			data->vectorOfNodalUpdateMaps(cacheIndex).begin(); it
			!= data->vectorOfNodalUpdateMaps(cacheIndex).end(); ++it) {
		unsigned int updatedNode = it->first;
		double nodalStatSum = computeSumNodalStatistics(updatedNode, beta);
		result += exp(nodalStatSum);
	}
	return result;
}

double NTCGroupedEgoCentricNetworkModel::computeKappaOfComingVertices(
		int cacheIndex, const DoubleVector& beta) {
	double result = 0;
	for (ListOfVertices::iterator it =
			data->vectorOfComingVertices(cacheIndex).begin(); it
			!= data->vectorOfComingVertices(cacheIndex).end(); ++it) {
		unsigned int updatedNode = *it;
		double nodalStatSum = computeSumNodalStatistics(updatedNode, beta);
		result += exp(nodalStatSum);
	}
	return result;
}

void NTCGroupedEgoCentricNetworkModel::computeLogLikelihoodOverEvents(
		DoubleVector& cummulativeLogLLH, const DoubleVector& beta,
		bool isUsingLBF) {
}

void NTCGroupedEgoCentricNetworkModel::computeLogLikelihoodOverEdges(
		ofstream& byEdgesOutputFile, const DoubleVector& beta, bool isUsingLBF) {

	if (!isUsingLBF)
		cout << "This feature has not implemented yet" << endl;

	// start rolling edges
	data->start();

	// the term that will be forwarded between two iterations
	double kappa = 0;

	unsigned int e = data->beginningEdgeIndex;

	// START - The first event

	// add those cited nodes
	DiscreteTimeVectorOfEdgeEvents edgeEvents =
			vectorOfDiscreteTimeVectorOfEdgeEvents[e];
	double currentCitingTime = edgeEvents.first;

	// compute kappa(1)
	for (unsigned int i = 0; i < numOfVertices; i++) {
		// only consider those nodes are already in the network
		if ((nodeJoiningTimes[i] < currentCitingTime)) {
			// sum over statistics of the node i
			double nodalStatSum = computeSumNodalStatistics(i, beta);
			kappa += exp(nodalStatSum);
		}
	}

	for (unsigned int k = 0; k < edgeEvents.second.size(); k++) {
		unsigned int citingNode = edgeEvents.second[k].first;
		for (unsigned int i = 0; i < edgeEvents.second[k].second.size(); i++) {
			unsigned int citedNode = edgeEvents.second[k].second[i];
			double logLLH = computeSumNodalStatistics(citedNode, beta) - log(
					kappa);
			byEdgesOutputFile << citingNode << "\t" << citedNode << "\t"
					<< currentCitingTime << "\t" << logLLH << endl;
		}
	}

	// compute the lookahead term
	double lookaheadTerm = computeKappaOfNodalUpdateMap(
			data->getCurrentCacheIndex(), beta);

	// Update statistics and add the edge to the current graph
	data->updateNodalNetworkStatistics(edgeEvents);
	e++;

	// END - kappa(1)

	// START - kappa(k) for k = 2,...,M
	for (; e <= data->endingEdgeIndex; e++) {
		// add those cited nodes
		DiscreteTimeVectorOfEdgeEvents edgeEvents =
				vectorOfDiscreteTimeVectorOfEdgeEvents[e];
		double currentCitingTime = edgeEvents.first;

		// Here we update kappa by adding the difference
		double forwardTerm = 0;
		forwardTerm += computeKappaOfNodalUpdateMap(
				data->getCurrentCacheIndex() - 1, beta);
		forwardTerm += computeKappaOfComingVertices(
				data->getCurrentCacheIndex() - 1, beta);
		kappa = kappa - lookaheadTerm + forwardTerm;

		for (unsigned int k = 0; k < edgeEvents.second.size(); k++) {
			unsigned int citingNode = edgeEvents.second[k].first;
			for (unsigned int i = 0; i < edgeEvents.second[k].second.size(); i++) {
				unsigned int citedNode = edgeEvents.second[k].second[i];
				double logLLH = computeSumNodalStatistics(citedNode, beta)
						- log(kappa);
				byEdgesOutputFile << citingNode << "\t" << citedNode << "\t"
						<< currentCitingTime << "\t" << logLLH << endl;
			}
		}

		lookaheadTerm = computeKappaOfNodalUpdateMap(
				data->getCurrentCacheIndex(), beta);

		// Update statistics and add the edge to the current graph
		data->updateNodalNetworkStatistics(edgeEvents);
	}
	// END - kappa(k)

	// finish rolling edges
	data->finish();
}

void NTCGroupedEgoCentricNetworkModel::computeScoreVectorAndInformationMatrix(
		DoubleVector& u, DoubleMatrix& I, const DoubleVector& beta) {
	u.clear();
	I.clear();

	// start rolling edges
	data->start();

	// the terms that will be forwarded between two iterations
	double kappa = 0;
	DoubleVector deltaKappa(data->numOfNodalNetworkStatistics);
	deltaKappa.clear();
	DoubleMatrix squaredDeltaKappa(data->numOfNodalNetworkStatistics,
			data->numOfNodalNetworkStatistics);
	squaredDeltaKappa.clear();

	unsigned int e = data->beginningEdgeIndex;

	// START - The first event

	// add those cited nodes
	DiscreteTimeVectorOfEdgeEvents edgeEvents =
			vectorOfDiscreteTimeVectorOfEdgeEvents[e];
	double currentCitingTime = edgeEvents.first;
	unsigned int numOfTiedEvents = 0;
	NodeSet citingNodes;
	for (unsigned int k = 0; k < edgeEvents.second.size(); k++) {
		citingNodes.insert(edgeEvents.second[k].first);
		numOfTiedEvents += edgeEvents.second[k].second.size();
		for (unsigned int i = 0; i < edgeEvents.second[k].second.size(); i++) {
			unsigned int citedNode = edgeEvents.second[k].second[i];
			// update u by statistics of the cited node
			for (unsigned int h = 0; h < data->numOfNodalNetworkStatistics; h++)
				u(h) += data->currentNodalNetworkStatistics(citedNode, h);
		}
	}

	// collect kappa, deltaKappa, and squaredDeltaKappa
	for (unsigned int i = 0; i < numOfVertices; i++) {
		// only consider those nodes are already in the network
		if ((nodeJoiningTimes[i] < currentCitingTime))
			computeNodalDerivatives(i, beta, kappa, deltaKappa,
					squaredDeltaKappa);
	}

	// update u and I
	updateUandI(u, I, kappa, deltaKappa, squaredDeltaKappa, numOfTiedEvents);

	// compute the lookahead terms
	double LA_Kappa;
	DoubleVector LA_DeltaKappa(data->numOfNodalNetworkStatistics);
	DoubleMatrix LA_SquaredDeltaKappa(data->numOfNodalNetworkStatistics,
			data->numOfNodalNetworkStatistics);
	computeKappaDerivativesOfNodalUpdateMap(data->getCurrentCacheIndex(), beta,
			LA_Kappa, LA_DeltaKappa, LA_SquaredDeltaKappa, true, citingNodes);

	// Update statistics and add the edge to the current graph
	data->updateNodalNetworkStatistics(edgeEvents);
	e++;

	// END - kappa(1)

	// START - kappa(k) for k = 2,...,M
	for (; e <= data->endingEdgeIndex; e++) {
		DiscreteTimeVectorOfEdgeEvents edgeEvents =
				vectorOfDiscreteTimeVectorOfEdgeEvents[e];
		double currentCitingTime = edgeEvents.first;
		unsigned int numOfTiedEvents = 0;
		NodeSet citingNodes;
		for (unsigned int k = 0; k < edgeEvents.second.size(); k++) {
			citingNodes.insert(edgeEvents.second[k].first);
			numOfTiedEvents += edgeEvents.second[k].second.size();
			for (unsigned int i = 0; i < edgeEvents.second[k].second.size(); i++) {
				unsigned int citedNode = edgeEvents.second[k].second[i];
				// update u by statistics of the cited node
				for (unsigned int h = 0; h < data->numOfNodalNetworkStatistics; h++)
					u(h) += data->currentNodalNetworkStatistics(citedNode, h);
			}
		}

		// compute the forward terms
		double FW_Kappa;
		DoubleVector FW_DeltaKappa(data->numOfNodalNetworkStatistics);
		DoubleMatrix FW_SquaredDeltaKappa(data->numOfNodalNetworkStatistics,
				data->numOfNodalNetworkStatistics);
		computeKappaDerivativesOfNodalUpdateMap(data->getCurrentCacheIndex()
				- 1, beta, FW_Kappa, FW_DeltaKappa, FW_SquaredDeltaKappa, true,
				citingNodes);
		computeKappaDerivativesOfComingVertices(data->getCurrentCacheIndex()
				- 1, beta, FW_Kappa, FW_DeltaKappa, FW_SquaredDeltaKappa, false);

		// here we update kappa, deltaKappa, squaredDeltaKappa by adding the difference
		kappa = kappa - LA_Kappa + FW_Kappa;
		for (unsigned int h = 0; h < data->numOfNodalNetworkStatistics; h++)
			deltaKappa(h) = deltaKappa(h) - LA_DeltaKappa(h) + FW_DeltaKappa(h);
		for (unsigned int h = 0; h < data->numOfNodalNetworkStatistics; h++)
			for (unsigned int g = 0; g <= h; g++)
				squaredDeltaKappa(h, g) = squaredDeltaKappa(h, g)
						- LA_SquaredDeltaKappa(h, g) + FW_SquaredDeltaKappa(h,
						g);

		// update u and I
		updateUandI(u, I, kappa, deltaKappa, squaredDeltaKappa, numOfTiedEvents);

		// look ahead
		computeKappaDerivativesOfNodalUpdateMap(data->getCurrentCacheIndex(),
				beta, LA_Kappa, LA_DeltaKappa, LA_SquaredDeltaKappa, true,
				citingNodes);

		// Update statistics and add the edge to the current graph
		data->updateNodalNetworkStatistics(edgeEvents);
	}
	// END - kappa(k)

	// finish rolling edges
	data->finish();

	// reflex the information matrix since it is symmetric.
	for (unsigned int h = 0; h < data->numOfNodalNetworkStatistics; h++)
		for (unsigned int g = 0; g < h; g++) {
			I(g, h) = I(h, g);
		}
}

void NTCGroupedEgoCentricNetworkModel::computeNodalDerivatives(
		unsigned int nodeID, const DoubleVector& beta, double& kappa,
		DoubleVector& deltaKappa, DoubleMatrix& squaredDeltaKappa) {
	// sum over statistics of the node i
	double nodalStatSum = computeSumNodalStatistics(nodeID, beta);

	// exp(sumOfStatistics)
	double temp = exp(nodalStatSum);

	// kappa
	kappa += temp;

	// deltaKappa
	for (unsigned int h = 0; h < data->numOfNodalNetworkStatistics; h++)
		deltaKappa(h) += data->currentNodalNetworkStatistics(nodeID, h) * temp;

	// squaredDelataKappa
	for (unsigned int h = 0; h < data->numOfNodalNetworkStatistics; h++)
		for (unsigned int g = 0; g <= h; g++)
			squaredDeltaKappa(h, g) += data->currentNodalNetworkStatistics(
					nodeID, h) * data->currentNodalNetworkStatistics(nodeID, g)
					* temp;
}

void NTCGroupedEgoCentricNetworkModel::computeKappaDerivativesOfNodalUpdateMap(
		int cacheIndex, const DoubleVector& beta, double& kappa,
		DoubleVector& deltaKappa, DoubleMatrix& squaredDeltaKappa,
		bool resetKappas, const NodeSet& citingNodes) {
	if (resetKappas) {
		kappa = 0;
		deltaKappa.clear();
		squaredDeltaKappa.clear();
	}

	for (NodalUpdateMap::iterator it =
			data->vectorOfNodalUpdateMaps(cacheIndex).begin(); it
			!= data->vectorOfNodalUpdateMaps(cacheIndex).end(); ++it) {
		// get the updated node
		unsigned int updatedNode = it->first;
		if (citingNodes.find(updatedNode) != citingNodes.end()) {
			//cout << "Here we go!!!!" << endl;
			continue;
		}
		computeNodalDerivatives(updatedNode, beta, kappa, deltaKappa,
				squaredDeltaKappa);
	}
}

void NTCGroupedEgoCentricNetworkModel::computeKappaDerivativesOfComingVertices(
		int cacheIndex, const DoubleVector& beta, double& kappa,
		DoubleVector& deltaKappa, DoubleMatrix& squaredDeltaKappa,
		bool resetKappas) {
	if (resetKappas) {
		kappa = 0;
		deltaKappa.clear();
		squaredDeltaKappa.clear();
	}

	for (ListOfVertices::iterator it =
			data->vectorOfComingVertices(cacheIndex).begin(); it
			!= data->vectorOfComingVertices(cacheIndex).end(); ++it) {
		// get the new node
		unsigned int newNode = *it;
		computeNodalDerivatives(newNode, beta, kappa, deltaKappa,
				squaredDeltaKappa);
	}
}

void NTCGroupedEgoCentricNetworkModel::updateUandI(DoubleVector& u,
		DoubleMatrix& I, const double& kappa, const DoubleVector& deltaKappa,
		const DoubleMatrix& squaredDeltaKappa, unsigned int numOfTiedEvents) {
	// update u
	for (unsigned int h = 0; h < data->numOfNodalNetworkStatistics; h++) {
		u(h) -= numOfTiedEvents * deltaKappa(h) / kappa;
	}
	// update I
	for (unsigned int h = 0; h < data->numOfNodalNetworkStatistics; h++)
		for (unsigned int g = 0; g <= h; g++)
			I(h, g) += numOfTiedEvents * (squaredDeltaKappa(h, g) / kappa
					- (deltaKappa(h) / kappa) * (deltaKappa(g) / kappa));
}

void NTCGroupedEgoCentricNetworkModel::getSumRanks(DoubleVector& sumRanks,
		const DoubleVector& beta, const std::string& tieMethod) {

}

void NTCGroupedEgoCentricNetworkModel::getRanks(VectorOfDoubleVector& ranks,
		const DoubleVector& beta, const std::string& tieMethod) {

}

void NTCGroupedEgoCentricNetworkModel::getRanks(ofstream& byEdgesOutputFile,
		ofstream& byEventsOutputFile, ofstream& accumlativeOutputFile,
		const DoubleVector& beta, const std::string& tieMethod) {
	// start rolling edges
	data->start();

	double accumRanks = 0;
	for (unsigned int e = data->beginningEdgeIndex, index = 0; e
			<= data->endingEdgeIndex; e++, index++) {
		if (index % 100 == 0)
			cout << "Working up to " << index << endl;

		DiscreteTimeVectorOfEdgeEvents edgeEvents =
				vectorOfDiscreteTimeVectorOfEdgeEvents[e];
		double currentCitingTime = edgeEvents.first;

		std::vector<double> intensities = std::vector<double>(numOfVertices);
		for (unsigned int i = 0; i < numOfVertices; i++) {
			if (nodeJoiningTimes[i] < currentCitingTime) {
				intensities[i] = exp(computeSumNodalStatistics(i, beta));
			} else
				intensities[i] = 0;
		}
		std::vector<double> currentNodeRanks = std::vector<double>(
				numOfVertices);
		ndip::rankhigh(intensities, currentNodeRanks, tieMethod);

		unsigned long eventAccumRanks = 0;
		unsigned int eventNumOfCitedNodes = 0;
		for (unsigned int k = 0; k < edgeEvents.second.size(); k++) {
			unsigned int citingNode = edgeEvents.second[k].first;
			for (unsigned int i = 0; i < edgeEvents.second[k].second.size(); i++) {
				unsigned int citedNode = edgeEvents.second[k].second[i];
				byEdgesOutputFile << citingNode << "\t" << citedNode << "\t"
						<< currentCitingTime << "\t"
						<< currentNodeRanks[citedNode] << endl;
				eventAccumRanks += currentNodeRanks[citedNode];
				eventNumOfCitedNodes++;
			}
		}
		byEventsOutputFile << currentCitingTime << "\t" << eventAccumRanks
				* 1.0 / eventNumOfCitedNodes << endl;

		accumRanks += eventAccumRanks;
		accumlativeOutputFile << currentCitingTime << "\t" << accumRanks
				<< endl;

		// Update statistics and add the edge to the current graph
		data->updateNodalNetworkStatistics(edgeEvents);
	}

	// finish rolling edges
	data->finish();
}
}
