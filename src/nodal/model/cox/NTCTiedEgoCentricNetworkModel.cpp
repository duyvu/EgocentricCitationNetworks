/*
 * NTCTiedEgoCentricNetworkModel.cpp
 *
 *  Created on: Jul 19, 2010
 *      Author: duyvu
 */

#include "NTCTiedEgoCentricNetworkModel.h"

#include "util/ranker.h"

#include <iostream>
#include <fstream>

namespace ndip {

NTCTiedEgoCentricNetworkModel::NTCTiedEgoCentricNetworkModel() {
	// TODO Auto-generated constructor stub

}

NTCTiedEgoCentricNetworkModel::~NTCTiedEgoCentricNetworkModel() {
	// TODO Auto-generated destructor stub
}

double NTCTiedEgoCentricNetworkModel::computeLogLikelihood(
		const DoubleVector& beta, bool isUsingLBF) {
	if (isUsingLBF) {
		return computeLogLikelihoodWithLBF(beta);
	} else {
		return computeLogLikelihoodWithoutLBF(beta);
	}
}

void NTCTiedEgoCentricNetworkModel::computeLogLikelihoodOverEvents(
		DoubleVector& cummulativeLogLLH, const DoubleVector& beta,
		bool isUsingLBF) {
	if (isUsingLBF) {
		computeLogLikelihoodOverEventsWithLBF(cummulativeLogLLH, beta);
	} else
		cout << "This feature has not implemented yet" << endl;
}

double NTCTiedEgoCentricNetworkModel::computeLogLikelihoodWithLBF(
		const DoubleVector& beta) {

	double logLLH = 0;

	// start rolling edges
	data->start();

	// the term that will be forwarded between two iterations
	double kappa = 0;

	unsigned int e = data->beginningEdgeIndex;

	// START - The first event

	// add those cited nodes
	unsigned int citingNode = vectorOfEdgeEvents[e].first;
	double currentCitingTime = nodeJoiningTimes[citingNode];
	unsigned int numOfTiedEvents = vectorOfEdgeEvents[e].second.size();
	for (unsigned int i = 0; i < numOfTiedEvents; i++) {
		unsigned int citedNode = vectorOfEdgeEvents[e].second[i];
		logLLH += computeSumNodalStatistics(citedNode, beta);
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
	data->updateNodalNetworkStatistics(vectorOfEdgeEvents[e]);
	e++;

	// END - kappa(1)

	// START - kappa(k) for k = 2,...,M
	for (; e <= data->endingEdgeIndex; e++) {
		// add those cited nodes
		unsigned int citingNode = vectorOfEdgeEvents[e].first;
		double currentCitingTime = nodeJoiningTimes[citingNode];
		unsigned int numOfTiedEvents = vectorOfEdgeEvents[e].second.size();
		for (unsigned int i = 0; i < numOfTiedEvents; i++) {
			unsigned int citedNode = vectorOfEdgeEvents[e].second[i];
			logLLH += computeSumNodalStatistics(citedNode, beta);
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
		data->updateNodalNetworkStatistics(vectorOfEdgeEvents[e]);
	}
	// END - kappa(k)

	// finish rolling edges
	data->finish();
	return logLLH;
}

void NTCTiedEgoCentricNetworkModel::computeLogLikelihoodOverEventsWithLBF(
		DoubleVector& cummulativeLogLLH, const DoubleVector& beta) {

	cummulativeLogLLH.resize(data->endingEdgeIndex - data->beginningEdgeIndex
			+ 1, false);
	cummulativeLogLLH.clear();

	// start rolling edges
	data->start();

	// the term that will be forwarded between two iterations
	double kappa = 0;

	unsigned int e = data->beginningEdgeIndex;
	unsigned int index = 0;
	// START - The first event

	// add those cited nodes
	unsigned int citingNode = vectorOfEdgeEvents[e].first;
	double currentCitingTime = nodeJoiningTimes[citingNode];
	unsigned int numOfTiedEvents = vectorOfEdgeEvents[e].second.size();
	for (unsigned int i = 0; i < numOfTiedEvents; i++) {
		unsigned int citedNode = vectorOfEdgeEvents[e].second[i];
		cummulativeLogLLH[index] += computeSumNodalStatistics(citedNode, beta);
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
	cummulativeLogLLH[index] -= numOfTiedEvents * log(kappa);

	// compute the lookahead term
	double lookaheadTerm = computeKappaOfNodalUpdateMap(
			data->getCurrentCacheIndex(), beta);

	// Update statistics and add the edge to the current graph
	data->updateNodalNetworkStatistics(vectorOfEdgeEvents[e]);
	e++;
	index++;

	// END - kappa(1)

	// START - kappa(k) for k = 2,...,M
	for (; e <= data->endingEdgeIndex; e++, index++) {
		// add those cited nodes
		unsigned int citingNode = vectorOfEdgeEvents[e].first;
		double currentCitingTime = nodeJoiningTimes[citingNode];
		unsigned int numOfTiedEvents = vectorOfEdgeEvents[e].second.size();
		for (unsigned int i = 0; i < numOfTiedEvents; i++) {
			unsigned int citedNode = vectorOfEdgeEvents[e].second[i];
			cummulativeLogLLH[index] += computeSumNodalStatistics(citedNode,
					beta);
		}

		// Here we update kappa by adding the difference
		double forwardTerm = 0;
		forwardTerm += computeKappaOfNodalUpdateMap(
				data->getCurrentCacheIndex() - 1, beta);
		forwardTerm += computeKappaOfComingVertices(
				data->getCurrentCacheIndex() - 1, beta);
		kappa = kappa - lookaheadTerm + forwardTerm;

		//cout << "kappa(" << e << "): " << kappa << endl;
		cummulativeLogLLH[index] -= numOfTiedEvents * log(kappa);

		// accummulate
		cummulativeLogLLH[index] += cummulativeLogLLH[index - 1];

		lookaheadTerm = computeKappaOfNodalUpdateMap(
				data->getCurrentCacheIndex(), beta);

		// Update statistics and add the edge to the current graph
		data->updateNodalNetworkStatistics(vectorOfEdgeEvents[e]);
	}
	// END - kappa(k)

	// finish rolling edges
	data->finish();
}

void NTCTiedEgoCentricNetworkModel::computeLogLikelihoodOverEdges(
		ofstream& byEdgesOutputFile, const DoubleVector& beta, bool isUsingLBF) {
	if (!isUsingLBF)
		cout << "This feature has not implemented yet" << endl;

	// start rolling edges
	data->start();

	// the term that will be forwarded between two iterations
	double kappa = 0;

	unsigned int e = data->beginningEdgeIndex;
	unsigned int index = 0;
	// START - The first event


	unsigned int citingNode = vectorOfEdgeEvents[e].first;
	double currentCitingTime = nodeJoiningTimes[citingNode];

	// compute kappa(1)
	for (unsigned int i = 0; i < numOfVertices; i++) {
		// only consider those nodes are already in the network
		if ((nodeJoiningTimes[i] < currentCitingTime)) {
			// sum over statistics of the node i
			double nodalStatSum = computeSumNodalStatistics(i, beta);
			kappa += exp(nodalStatSum);
		}
	}

	// add those cited nodes
	unsigned int numOfTiedEvents = vectorOfEdgeEvents[e].second.size();
	for (unsigned int i = 0; i < numOfTiedEvents; i++) {
		unsigned int citedNode = vectorOfEdgeEvents[e].second[i];
		byEdgesOutputFile << citingNode << "\t" << citedNode << "\t"
				<< currentCitingTime << "\t" << computeSumNodalStatistics(
				citedNode, beta) - log(kappa) << endl;
	}

	// compute the lookahead term
	double lookaheadTerm = computeKappaOfNodalUpdateMap(
			data->getCurrentCacheIndex(), beta);

	// Update statistics and add the edge to the current graph
	data->updateNodalNetworkStatistics(vectorOfEdgeEvents[e]);
	e++;
	index++;

	// END - kappa(1)

	// START - kappa(k) for k = 2,...,M
	for (; e <= data->endingEdgeIndex; e++, index++) {

		unsigned int citingNode = vectorOfEdgeEvents[e].first;
		double currentCitingTime = nodeJoiningTimes[citingNode];

		// Here we update kappa by adding the difference
		double forwardTerm = 0;
		forwardTerm += computeKappaOfNodalUpdateMap(
				data->getCurrentCacheIndex() - 1, beta);
		forwardTerm += computeKappaOfComingVertices(
				data->getCurrentCacheIndex() - 1, beta);
		kappa = kappa - lookaheadTerm + forwardTerm;

		lookaheadTerm = computeKappaOfNodalUpdateMap(
				data->getCurrentCacheIndex(), beta);

		// add those cited nodes
		unsigned int numOfTiedEvents = vectorOfEdgeEvents[e].second.size();
		for (unsigned int i = 0; i < numOfTiedEvents; i++) {
			unsigned int citedNode = vectorOfEdgeEvents[e].second[i];
			byEdgesOutputFile << citingNode << "\t" << citedNode << "\t"
					<< currentCitingTime << "\t" << computeSumNodalStatistics(
					citedNode, beta) - log(kappa) << endl;
		}

		// Update statistics and add the edge to the current graph
		data->updateNodalNetworkStatistics(vectorOfEdgeEvents[e]);
	}
	// END - kappa(k)

	// finish rolling edges
	data->finish();
}

double NTCTiedEgoCentricNetworkModel::computeSumNodalStatistics(
		unsigned int nodeID, const DoubleVector& beta) {
	double nodalStatSum = 0;
	for (unsigned int statIdx = 0; statIdx < data->numOfNodalNetworkStatistics; statIdx++)
		nodalStatSum += beta(statIdx) * data->currentNodalNetworkStatistics(
				nodeID, statIdx);
	return nodalStatSum;
}

double NTCTiedEgoCentricNetworkModel::computeKappaOfNodalUpdateMap(
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

double NTCTiedEgoCentricNetworkModel::computeKappaOfComingVertices(
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

double NTCTiedEgoCentricNetworkModel::computeLogLikelihoodWithoutLBF(
		const DoubleVector& beta) {

	double logLLH = 0;

	// start rolling edges
	data->start();

	// cache the exposure indicator for better performance
	// checking a bool indicator is faster than comparing
	// two double values of the current time and the joining time of a node
	bool exposureIndicators[numOfVertices];
	prepareExposureIndicators(exposureIndicators);

	// iteratively calculate the likelihood edge by edge
	// HERE we ASSUME that the observation time coincides with the time of the first edge
	// under consideration
	for (unsigned int e = data->beginningEdgeIndex; e <= data->endingEdgeIndex; e++) {
		// add those cited nodes
		unsigned int citingNode = vectorOfEdgeEvents[e].first;
		double currentCitingTime = nodeJoiningTimes[citingNode];
		unsigned int numOfTiedEvents = vectorOfEdgeEvents[e].second.size();
		for (unsigned int i = 0; i < numOfTiedEvents; i++) {
			unsigned int citedNode = vectorOfEdgeEvents[e].second[i];
			logLLH += computeSumNodalStatistics(citedNode, beta);
		}

		double kappa = 0;
		for (unsigned int i = 0; i < numOfVertices; i++) {
			// only consider those nodes are already in the network
			if (exposureIndicators[i] || (nodeJoiningTimes[i]
					< currentCitingTime)) {
				// in case the exposure of i has not been set up, we set it
				// this happens when i joins the network in [previousEdgeTime, currentEdgeTime)
				if (!exposureIndicators[i])
					exposureIndicators[i] = true;
				double nodalStatSum = computeSumNodalStatistics(i, beta);
				kappa += exp(nodalStatSum);
			}
		}

		//cout << "kappa(" << e << "): " << kappa << endl;
		logLLH -= numOfTiedEvents * log(kappa);

		// Update statistics and add the edge to the current graph
		data->updateNodalNetworkStatistics(vectorOfEdgeEvents[e]);
	}

	// finish rolling edges
	data->finish();

	return logLLH;
}

void NTCTiedEgoCentricNetworkModel::computeScoreVectorAndInformationMatrix(
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
	unsigned int citingNode = vectorOfEdgeEvents[e].first;
	double currentCitingTime = nodeJoiningTimes[citingNode];
	unsigned int numOfTiedEvents = vectorOfEdgeEvents[e].second.size();
	for (unsigned int i = 0; i < numOfTiedEvents; i++) {
		unsigned int citedNode = vectorOfEdgeEvents[e].second[i];
		// update u by statistics of the cited node
		for (unsigned int h = 0; h < data->numOfNodalNetworkStatistics; h++)
			u(h) += data->currentNodalNetworkStatistics(citedNode, h);
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
			LA_Kappa, LA_DeltaKappa, LA_SquaredDeltaKappa, true, citingNode);

	// Update statistics and add the edge to the current graph
	data->updateNodalNetworkStatistics(vectorOfEdgeEvents[e]);
	e++;

	// END - kappa(1)

	// START - kappa(k) for k = 2,...,M
	for (; e <= data->endingEdgeIndex; e++) {

		unsigned int citingNode = vectorOfEdgeEvents[e].first;
		double currentCitingTime = nodeJoiningTimes[citingNode];
		unsigned int numOfTiedEvents = vectorOfEdgeEvents[e].second.size();
		for (unsigned int i = 0; i < numOfTiedEvents; i++) {
			unsigned int citedNode = vectorOfEdgeEvents[e].second[i];
			// update u by statistics of the cited node
			for (unsigned int h = 0; h < data->numOfNodalNetworkStatistics; h++)
				u(h) += data->currentNodalNetworkStatistics(citedNode, h);
		}

		// compute the forward terms
		double FW_Kappa;
		DoubleVector FW_DeltaKappa(data->numOfNodalNetworkStatistics);
		DoubleMatrix FW_SquaredDeltaKappa(data->numOfNodalNetworkStatistics,
				data->numOfNodalNetworkStatistics);
		computeKappaDerivativesOfNodalUpdateMap(data->getCurrentCacheIndex()
				- 1, beta, FW_Kappa, FW_DeltaKappa, FW_SquaredDeltaKappa, true,
				citingNode);
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
				citingNode);

		// Update statistics and add the edge to the current graph
		data->updateNodalNetworkStatistics(vectorOfEdgeEvents[e]);
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

void NTCTiedEgoCentricNetworkModel::computeNodalDerivatives(
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

void NTCTiedEgoCentricNetworkModel::computeKappaDerivativesOfNodalUpdateMap(
		int cacheIndex, const DoubleVector& beta, double& kappa,
		DoubleVector& deltaKappa, DoubleMatrix& squaredDeltaKappa,
		bool resetKappas, unsigned int citingNode) {

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
		if (updatedNode == citingNode) {
			//cout << "Here we go!!!!" << endl;
			continue;
		}
		computeNodalDerivatives(updatedNode, beta, kappa, deltaKappa,
				squaredDeltaKappa);
	}
}

void NTCTiedEgoCentricNetworkModel::computeKappaDerivativesOfComingVertices(
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

void NTCTiedEgoCentricNetworkModel::updateUandI(DoubleVector& u,
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

double NTCTiedEgoCentricNetworkModel::computeLogLikelihoodAndScoreVector(
		DoubleVector& u, const DoubleVector& beta) {

	double logLLH = 0;
	u.clear();

	// start rolling edges
	data->start();

	// the term that will be forwarded between two iterations
	double kappa = 0;
	DoubleVector deltaKappa(data->numOfNodalNetworkStatistics);
	deltaKappa.clear();

	unsigned int e = data->beginningEdgeIndex;

	// START - The first event

	// add those cited nodes
	unsigned int citingNode = vectorOfEdgeEvents[e].first;
	double currentCitingTime = nodeJoiningTimes[citingNode];
	unsigned int numOfTiedEvents = vectorOfEdgeEvents[e].second.size();
	for (unsigned int i = 0; i < numOfTiedEvents; i++) {
		unsigned int citedNode = vectorOfEdgeEvents[e].second[i];
		logLLH += computeSumNodalStatistics(citedNode, beta);
		// update u by statistics of the cited node
		for (unsigned int h = 0; h < data->numOfNodalNetworkStatistics; h++)
			u(h) += data->currentNodalNetworkStatistics(citedNode, h);
	}

	// compute kappa(1)
	for (unsigned int i = 0; i < numOfVertices; i++) {
		// only consider those nodes are already in the network
		if ((nodeJoiningTimes[i] < currentCitingTime))
			computeNodalGradients(i, beta, kappa, deltaKappa);
	}

	//cout << "kappa(" << e << "): " << kappa << endl;
	// update logLLH
	logLLH -= numOfTiedEvents * log(kappa);
	// update u
	updateU(u, kappa, deltaKappa, numOfTiedEvents);

	// compute the lookahead terms
	double LA_Kappa;
	DoubleVector LA_DeltaKappa(data->numOfNodalNetworkStatistics);
	computeKappaGradientsOfNodalUpdateMap(data->getCurrentCacheIndex(), beta,
			LA_Kappa, LA_DeltaKappa, true, citingNode);

	// Update statistics and add the edge to the current graph
	data->updateNodalNetworkStatistics(vectorOfEdgeEvents[e]);
	e++;

	// END - kappa(1)

	// START - kappa(k) for k = 2,...,M
	for (; e <= data->endingEdgeIndex; e++) {
		// add those cited nodes
		unsigned int citingNode = vectorOfEdgeEvents[e].first;
		double currentCitingTime = nodeJoiningTimes[citingNode];
		unsigned int numOfTiedEvents = vectorOfEdgeEvents[e].second.size();
		for (unsigned int i = 0; i < numOfTiedEvents; i++) {
			unsigned int citedNode = vectorOfEdgeEvents[e].second[i];
			logLLH += computeSumNodalStatistics(citedNode, beta);
			// update u by statistics of the cited node
			for (unsigned int h = 0; h < data->numOfNodalNetworkStatistics; h++)
				u(h) += data->currentNodalNetworkStatistics(citedNode, h);
		}

		// compute the forward terms
		double FW_Kappa;
		DoubleVector FW_DeltaKappa(data->numOfNodalNetworkStatistics);
		computeKappaGradientsOfNodalUpdateMap(data->getCurrentCacheIndex() - 1,
				beta, FW_Kappa, FW_DeltaKappa, true, citingNode);
		computeKappaGradientsOfComingVertices(data->getCurrentCacheIndex() - 1,
				beta, FW_Kappa, FW_DeltaKappa, false);

		// here we update kappa, deltaKappa by adding the difference
		kappa = kappa - LA_Kappa + FW_Kappa;
		for (unsigned int h = 0; h < data->numOfNodalNetworkStatistics; h++)
			deltaKappa(h) = deltaKappa(h) - LA_DeltaKappa(h) + FW_DeltaKappa(h);

		//cout << "kappa(" << e << "): " << kappa << endl;
		// update logLLH
		logLLH -= numOfTiedEvents * log(kappa);
		// update u
		updateU(u, kappa, deltaKappa, numOfTiedEvents);

		// look ahead
		computeKappaGradientsOfNodalUpdateMap(data->getCurrentCacheIndex(),
				beta, LA_Kappa, LA_DeltaKappa, true, citingNode);

		// Update statistics and add the edge to the current graph
		data->updateNodalNetworkStatistics(vectorOfEdgeEvents[e]);
	}
	// END - kappa(k)

	// finish rolling edges
	data->finish();
	return logLLH;
}

void NTCTiedEgoCentricNetworkModel::computeNodalGradients(unsigned int nodeID,
		const DoubleVector& beta, double& kappa, DoubleVector& deltaKappa) {
	// sum over statistics of the node i
	double nodalStatSum = computeSumNodalStatistics(nodeID, beta);

	// exp(sumOfStatistics)
	double temp = exp(nodalStatSum);

	// kappa
	kappa += temp;

	// deltaKappa
	for (unsigned int h = 0; h < data->numOfNodalNetworkStatistics; h++)
		deltaKappa(h) += data->currentNodalNetworkStatistics(nodeID, h) * temp;
}

void NTCTiedEgoCentricNetworkModel::computeKappaGradientsOfNodalUpdateMap(
		int cacheIndex, const DoubleVector& beta, double& kappa,
		DoubleVector& deltaKappa, bool resetKappas, unsigned int citingNode) {
	if (resetKappas) {
		kappa = 0;
		deltaKappa.clear();
	}

	for (NodalUpdateMap::iterator it =
			data->vectorOfNodalUpdateMaps(cacheIndex).begin(); it
			!= data->vectorOfNodalUpdateMaps(cacheIndex).end(); ++it) {
		// get the updated node
		unsigned int updatedNode = it->first;
		if (updatedNode == citingNode) {
			//cout << "Here we go!!!!" << endl;
			continue;
		}
		computeNodalGradients(updatedNode, beta, kappa, deltaKappa);
	}

}

void NTCTiedEgoCentricNetworkModel::computeKappaGradientsOfComingVertices(
		int cacheIndex, const DoubleVector& beta, double& kappa,
		DoubleVector& deltaKappa, bool resetKappas) {
	if (resetKappas) {
		kappa = 0;
		deltaKappa.clear();
	}

	for (ListOfVertices::iterator it =
			data->vectorOfComingVertices(cacheIndex).begin(); it
			!= data->vectorOfComingVertices(cacheIndex).end(); ++it) {
		// get the new node
		unsigned int newNode = *it;
		computeNodalGradients(newNode, beta, kappa, deltaKappa);
	}

}

void NTCTiedEgoCentricNetworkModel::updateU(DoubleVector& u,
		const double& kappa, const DoubleVector& deltaKappa,
		unsigned int numOfTiedEvents) {
	// update u
	for (unsigned int h = 0; h < data->numOfNodalNetworkStatistics; h++) {
		u(h) -= numOfTiedEvents * deltaKappa(h) / kappa;
	}
}

void NTCTiedEgoCentricNetworkModel::computeOnlineGradientAscentMLE(
		DoubleVector& beta, double learningRate) {

	DoubleVector stocGrad = DoubleVector(beta.size());
	stocGrad.clear();

	// start rolling edges
	data->start();

	// the term that will be forwarded between two iterations
	double kappa = 0;
	DoubleVector deltaKappa(data->numOfNodalNetworkStatistics);
	deltaKappa.clear();

	unsigned int e = data->beginningEdgeIndex;

	// START - The first event

	// add those cited nodes
	unsigned int citingNode = vectorOfEdgeEvents[e].first;
	double currentCitingTime = nodeJoiningTimes[citingNode];
	unsigned int numOfTiedEvents = vectorOfEdgeEvents[e].second.size();
	for (unsigned int i = 0; i < numOfTiedEvents; i++) {
		unsigned int citedNode = vectorOfEdgeEvents[e].second[i];
		// update u by statistics of the cited node
		for (unsigned int h = 0; h < data->numOfNodalNetworkStatistics; h++)
			stocGrad(h) += data->currentNodalNetworkStatistics(citedNode, h);
	}

	// compute kappa(1)
	for (unsigned int i = 0; i < numOfVertices; i++) {
		// only consider those nodes are already in the network
		if ((nodeJoiningTimes[i] < currentCitingTime))
			computeNodalGradients(i, beta, kappa, deltaKappa);
	}

	// update u
	updateU(stocGrad, kappa, deltaKappa, numOfTiedEvents);

	// update beta
	//	cout << "stocGrad = ";
	//	for (unsigned int k = 0; k < stocGrad.size(); k++)
	//		cout << stocGrad[k] << " ";
	//	cout << endl;
	beta += learningRate * stocGrad;
	//	cout << "beta = ";
	//	for (unsigned int k = 0; k < beta.size(); k++)
	//		cout << beta[k] << " ";
	//	cout << endl;

	// compute the lookahead terms
	double LA_Kappa;
	DoubleVector LA_DeltaKappa(data->numOfNodalNetworkStatistics);
	computeKappaGradientsOfNodalUpdateMap(data->getCurrentCacheIndex(), beta,
			LA_Kappa, LA_DeltaKappa, true, citingNode);

	// Update statistics and add the edge to the current graph
	data->updateNodalNetworkStatistics(vectorOfEdgeEvents[e]);
	e++;

	// END - kappa(1)

	// START - kappa(k) for k = 2,...,M
	for (; e <= data->endingEdgeIndex; e++) {

		stocGrad.clear();

		// add those cited nodes
		unsigned int citingNode = vectorOfEdgeEvents[e].first;
		double currentCitingTime = nodeJoiningTimes[citingNode];
		unsigned int numOfTiedEvents = vectorOfEdgeEvents[e].second.size();
		for (unsigned int i = 0; i < numOfTiedEvents; i++) {
			unsigned int citedNode = vectorOfEdgeEvents[e].second[i];
			// update u by statistics of the cited node
			for (unsigned int h = 0; h < data->numOfNodalNetworkStatistics; h++)
				stocGrad(h)
						+= data->currentNodalNetworkStatistics(citedNode, h);
		}

		// compute the forward terms
		double FW_Kappa;
		DoubleVector FW_DeltaKappa(data->numOfNodalNetworkStatistics);
		computeKappaGradientsOfNodalUpdateMap(data->getCurrentCacheIndex() - 1,
				beta, FW_Kappa, FW_DeltaKappa, true, citingNode);
		computeKappaGradientsOfComingVertices(data->getCurrentCacheIndex() - 1,
				beta, FW_Kappa, FW_DeltaKappa, false);

		// here we update kappa, deltaKappa by adding the difference
		kappa = kappa - LA_Kappa + FW_Kappa;
		for (unsigned int h = 0; h < data->numOfNodalNetworkStatistics; h++)
			deltaKappa(h) = deltaKappa(h) - LA_DeltaKappa(h) + FW_DeltaKappa(h);

		// update u
		updateU(stocGrad, kappa, deltaKappa, numOfTiedEvents);

		// update beta
		//		cout << "stocGrad = ";
		//		for (unsigned int k = 0; k < stocGrad.size(); k++)
		//			cout << stocGrad[k] << " ";
		//		cout << endl;
		beta += learningRate * stocGrad;
		//		cout << "beta = ";
		//		for (unsigned int k = 0; k < beta.size(); k++)
		//			cout << beta[k] << " ";
		//		cout << endl;

		// look ahead
		computeKappaGradientsOfNodalUpdateMap(data->getCurrentCacheIndex(),
				beta, LA_Kappa, LA_DeltaKappa, true, citingNode);

		// Update statistics and add the edge to the current graph
		data->updateNodalNetworkStatistics(vectorOfEdgeEvents[e]);
	}
	// END - kappa(k)

	// finish rolling edges
	data->finish();
}

void NTCTiedEgoCentricNetworkModel::computeCummulativeBaselineHazard(
		VectorOfCummulativePoints& cumBH, const DoubleVector& beta) {

	cumBH.resize(data->endingEdgeIndex - data->beginningEdgeIndex + 1, false);

	// start rolling edges
	data->start();

	// the term that will be forwarded between two iterations
	double kappa = 0;

	unsigned int e = data->beginningEdgeIndex;

	// START - The first event
	unsigned int citingNode = vectorOfEdgeEvents[e].first;
	double currentCitingTime = nodeJoiningTimes[citingNode];
	unsigned int numOfTiedEvents = vectorOfEdgeEvents[e].second.size();
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
	cumBH(0) = make_pair(currentCitingTime, numOfTiedEvents / kappa);

	// compute the lookahead term
	double lookaheadTerm = computeKappaOfNodalUpdateMap(
			data->getCurrentCacheIndex(), beta);

	// Update statistics and add the edge to the current graph
	data->updateNodalNetworkStatistics(vectorOfEdgeEvents[e]);
	e++;

	// END - kappa(1)

	// START - kappa(k) for k = 2,...,M
	for (int index = 1; e <= data->endingEdgeIndex; e++, index++) {
		// The next event
		unsigned int citingNode = vectorOfEdgeEvents[e].first;
		double currentCitingTime = nodeJoiningTimes[citingNode];
		unsigned int numOfTiedEvents = vectorOfEdgeEvents[e].second.size();

		// Here we update kappa by adding the difference
		double forwardTerm = 0;
		forwardTerm += computeKappaOfNodalUpdateMap(
				data->getCurrentCacheIndex() - 1, beta);
		forwardTerm += computeKappaOfComingVertices(
				data->getCurrentCacheIndex() - 1, beta);
		kappa = kappa - lookaheadTerm + forwardTerm;

		//cout << "kappa(" << e << "): " << kappa << endl;
		cumBH(index) = make_pair(currentCitingTime, cumBH(index - 1).second
				+ numOfTiedEvents / kappa);

		lookaheadTerm = computeKappaOfNodalUpdateMap(
				data->getCurrentCacheIndex(), beta);

		// Update statistics and add the edge to the current graph
		data->updateNodalNetworkStatistics(vectorOfEdgeEvents[e]);
	}
	// END - kappa(k)

	// finish rolling edges
	data->finish();
}

void NTCTiedEgoCentricNetworkModel::computeSchoenfeldResiduals(
		DoubleMatrix& schoenfeldResiduals, const DoubleVector& beta,
		bool isPlusBeta) {

	unsigned int numOfCovariates = data->numOfNodalNetworkStatistics;
	schoenfeldResiduals.resize(data->endingEdgeIndex - data->beginningEdgeIndex
			+ 1, numOfCovariates, false);
	schoenfeldResiduals.clear();
	DoubleVector residuals = DoubleVector(numOfCovariates);
	residuals.clear();

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
	residuals.clear();
	unsigned int citingNode = vectorOfEdgeEvents[e].first;
	double currentCitingTime = nodeJoiningTimes[citingNode];
	unsigned int numOfTiedEvents = vectorOfEdgeEvents[e].second.size();
	for (unsigned int i = 0; i < numOfTiedEvents; i++) {
		unsigned int citedNode = vectorOfEdgeEvents[e].second[i];
		// update u by statistics of the cited node
		for (unsigned int h = 0; h < data->numOfNodalNetworkStatistics; h++)
			residuals(h) += data->currentNodalNetworkStatistics(citedNode, h);
	}

	// collect kappa, deltaKappa, and squaredDeltaKappa
	for (unsigned int i = 0; i < numOfVertices; i++) {
		// only consider those nodes are already in the network
		if ((nodeJoiningTimes[i] < currentCitingTime))
			computeNodalDerivatives(i, beta, kappa, deltaKappa,
					squaredDeltaKappa);
	}

	// compute the residuals
	computeSchoenfeldResiduals(residuals, kappa, deltaKappa, squaredDeltaKappa,
			numOfTiedEvents);
	// add to the matrix result
	for (unsigned int k = 0; k < numOfCovariates; k++)
		schoenfeldResiduals(0, k) = residuals(k);

	// compute the lookahead terms
	double LA_Kappa;
	DoubleVector LA_DeltaKappa(data->numOfNodalNetworkStatistics);
	DoubleMatrix LA_SquaredDeltaKappa(data->numOfNodalNetworkStatistics,
			data->numOfNodalNetworkStatistics);
	computeKappaDerivativesOfNodalUpdateMap(data->getCurrentCacheIndex(), beta,
			LA_Kappa, LA_DeltaKappa, LA_SquaredDeltaKappa, true, citingNode);

	// Update statistics and add the edge to the current graph
	data->updateNodalNetworkStatistics(vectorOfEdgeEvents[e]);
	e++;

	// END - kappa(1)

	// START - kappa(k) for k = 2,...,M
	for (int index = 1; e <= data->endingEdgeIndex; e++, index++) {
		residuals.clear();
		unsigned int citingNode = vectorOfEdgeEvents[e].first;
		double currentCitingTime = nodeJoiningTimes[citingNode];
		unsigned int numOfTiedEvents = vectorOfEdgeEvents[e].second.size();
		for (unsigned int i = 0; i < numOfTiedEvents; i++) {
			unsigned int citedNode = vectorOfEdgeEvents[e].second[i];
			// update u by statistics of the cited node
			for (unsigned int h = 0; h < data->numOfNodalNetworkStatistics; h++)
				residuals(h) += data->currentNodalNetworkStatistics(citedNode,
						h);
		}

		// compute the forward terms
		double FW_Kappa;
		DoubleVector FW_DeltaKappa(data->numOfNodalNetworkStatistics);
		DoubleMatrix FW_SquaredDeltaKappa(data->numOfNodalNetworkStatistics,
				data->numOfNodalNetworkStatistics);
		computeKappaDerivativesOfNodalUpdateMap(data->getCurrentCacheIndex()
				- 1, beta, FW_Kappa, FW_DeltaKappa, FW_SquaredDeltaKappa, true,
				citingNode);
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

		// compute the residuals
		computeSchoenfeldResiduals(residuals, kappa, deltaKappa,
				squaredDeltaKappa, numOfTiedEvents);
		for (unsigned int k = 0; k < numOfCovariates; k++)
			schoenfeldResiduals(index, k) = residuals(k);

		computeKappaDerivativesOfNodalUpdateMap(data->getCurrentCacheIndex(),
				beta, LA_Kappa, LA_DeltaKappa, LA_SquaredDeltaKappa, true,
				citingNode);

		// Update statistics and add the edge to the current graph
		data->updateNodalNetworkStatistics(vectorOfEdgeEvents[e]);
	}
	// END - kappa(k)

	// finish rolling edges
	data->finish();

	if (isPlusBeta)
		for (unsigned int i = 0; i < schoenfeldResiduals.size1(); i++)
			for (unsigned int k = 0; k < schoenfeldResiduals.size2(); k++)
				schoenfeldResiduals(i, k) += beta(k);
}

void NTCTiedEgoCentricNetworkModel::computeSchoenfeldResiduals(
		DoubleVector& residuals, double kappa, const DoubleVector& deltaKappa,
		const DoubleMatrix& squaredDeltaKappa, unsigned int numOfTiedEvents) {
	// update residuals
	unsigned int numOfCovariates = residuals.size();
	for (unsigned int h = 0; h < numOfCovariates; h++) {
		residuals(h) -= numOfTiedEvents * deltaKappa(h) / kappa;
	}
	//	cout << "Residual vector: " << endl;
	//	for (int h = 0; h < residuals.size(); h++)
	//		cout << residuals(h) << "  ";
	//	cout << endl;

	// compute V
	DoubleMatrix V = DoubleMatrix(numOfCovariates, numOfCovariates);
	for (unsigned int h = 0; h < numOfCovariates; h++)
		for (unsigned int g = 0; g <= h; g++) {
			V(h, g) = numOfTiedEvents * (squaredDeltaKappa(h, g) / kappa
					- (deltaKappa(h) / kappa) * (deltaKappa(g) / kappa));
			V(g, h) = V(h, g);
		}
	//	cout << "V vector: " << endl;
	//	for (int k = 0; k < V.size1(); k++) {
	//		for (int h = 0; h < V.size2(); h++)
	//			cout << V(k, h) << "  ";
	//		cout << endl;
	//	}

	// V^{-1} * u, standardize residuals
	boost::numeric::ublas::permutation_matrix<> pm(numOfCovariates);
	lu_factorize(V, pm);
	lu_substitute(V, pm, residuals);
}

void NTCTiedEgoCentricNetworkModel::getSumRanks(DoubleVector& sumRanks,
		const DoubleVector& beta, const std::string& tieMethod) {
	sumRanks.resize(data->endingEdgeIndex - data->beginningEdgeIndex + 1, false);
	sumRanks.clear();

	// start rolling edges
	data->start();

	for (unsigned int e = data->beginningEdgeIndex, index = 0; e
			<= data->endingEdgeIndex; e++, index++) {
		if (index % 100 == 0)
			cout << "Working up to " << index << endl;

		unsigned int citingNode = vectorOfEdgeEvents[e].first;
		double currentCitingTime = nodeJoiningTimes[citingNode];

		std::vector<double> intensities = std::vector<double>(numOfVertices);
		for (unsigned int i = 0; i < numOfVertices; i++) {
			if (nodeJoiningTimes[i] < currentCitingTime) {
				double nodeStatSum = 0;
				for (unsigned int h = 0; h < data->numOfNodalNetworkStatistics; h++)
					nodeStatSum += beta(h)
							* data->currentNodalNetworkStatistics(i, h);
				intensities[i] = exp(nodeStatSum);
			} else
				intensities[i] = 0;
		}

		std::vector<double> currentNodeRanks = std::vector<double>(
				numOfVertices);
		ndip::rankhigh(intensities, currentNodeRanks, tieMethod);
		unsigned int numOfCitedNodes = vectorOfEdgeEvents[e].second.size();
		for (unsigned int i = 0; i < numOfCitedNodes; i++) {
			unsigned int citedNode = vectorOfEdgeEvents[e].second[i];
			sumRanks(index) += currentNodeRanks[citedNode];
		}

		// cummulative
		if (index > 0)
			sumRanks(index) += sumRanks(index - 1);

		// Update statistics and add the edge to the current graph
		data->updateNodalNetworkStatistics(vectorOfEdgeEvents[e]);
	}

	// finish rolling edges
	data->finish();
}

void NTCTiedEgoCentricNetworkModel::getRanks(VectorOfDoubleVector& ranks,
		const DoubleVector& beta, const std::string& tieMethod) {
	ranks.resize(data->endingEdgeIndex - data->beginningEdgeIndex + 1, false);

	// start rolling edges
	data->start();

	for (unsigned int e = data->beginningEdgeIndex, index = 0; e
			<= data->endingEdgeIndex; e++, index++) {
		if (index % 100 == 0)
			cout << "Working up to " << index << endl;

		unsigned int citingNode = vectorOfEdgeEvents[e].first;
		double currentCitingTime = nodeJoiningTimes[citingNode];

		std::vector<double> intensities = std::vector<double>(numOfVertices);
		for (unsigned int i = 0; i < numOfVertices; i++) {
			if (nodeJoiningTimes[i] < currentCitingTime) {
				double nodeStatSum = 0;
				for (unsigned int h = 0; h < data->numOfNodalNetworkStatistics; h++)
					nodeStatSum += beta(h)
							* data->currentNodalNetworkStatistics(i, h);
				intensities[i] = exp(nodeStatSum);
			} else
				intensities[i] = 0;
		}

		std::vector<double> currentNodeRanks = std::vector<double>(
				numOfVertices);
		ndip::rankhigh(intensities, currentNodeRanks, tieMethod);
		unsigned int numOfCitedNodes = vectorOfEdgeEvents[e].second.size();
		ranks(index) = DoubleVector(numOfCitedNodes);
		for (unsigned int i = 0; i < numOfCitedNodes; i++) {
			unsigned int citedNode = vectorOfEdgeEvents[e].second[i];
			ranks(index)(i) = currentNodeRanks[citedNode];
		}

		// Update statistics and add the edge to the current graph
		data->updateNodalNetworkStatistics(vectorOfEdgeEvents[e]);
	}

	// finish rolling edges
	data->finish();
}

void NTCTiedEgoCentricNetworkModel::getRanks(ofstream& byEdgesOutputFile,
		ofstream& byEventsOutputFile, ofstream& accumlativeOutputFile,
		const DoubleVector& beta, const std::string& tieMethod) {
	// start rolling edges
	data->start();

	double accumRanks = 0;
	for (unsigned int e = data->beginningEdgeIndex, index = 0; e
			<= data->endingEdgeIndex; e++, index++) {
		if (index % 100 == 0)
			cout << "Working up to " << index << endl;

		unsigned int citingNode = vectorOfEdgeEvents[e].first;
		double currentCitingTime = nodeJoiningTimes[citingNode];

		std::vector<double> intensities = std::vector<double>(numOfVertices);
		for (unsigned int i = 0; i < numOfVertices; i++) {
			if (nodeJoiningTimes[i] < currentCitingTime) {
				double nodeStatSum = 0;
				for (unsigned int h = 0; h < data->numOfNodalNetworkStatistics; h++)
					nodeStatSum += beta(h)
							* data->currentNodalNetworkStatistics(i, h);
				intensities[i] = exp(nodeStatSum);
			} else
				intensities[i] = 0;
		}

		std::vector<double> currentNodeRanks = std::vector<double>(
				numOfVertices);
		ndip::rankhigh(intensities, currentNodeRanks, tieMethod);

		unsigned long eventAccumRanks = 0;
		unsigned int numOfCitedNodes = vectorOfEdgeEvents[e].second.size();
		for (unsigned int i = 0; i < numOfCitedNodes; i++) {
			unsigned int citedNode = vectorOfEdgeEvents[e].second[i];
			byEdgesOutputFile << citingNode << "\t" << citedNode << "\t"
					<< currentCitingTime << "\t" << currentNodeRanks[citedNode]
					<< endl;
			eventAccumRanks += currentNodeRanks[citedNode];
		}
		byEventsOutputFile << currentCitingTime << "\t" << eventAccumRanks
				* 1.0 / numOfCitedNodes << endl;

		accumRanks += eventAccumRanks;
		accumlativeOutputFile << currentCitingTime << "\t" << accumRanks
				<< endl;

		// Update statistics and add the edge to the current graph
		data->updateNodalNetworkStatistics(vectorOfEdgeEvents[e]);
	}

	// finish rolling edges
	data->finish();
}

void NTCTiedEgoCentricNetworkModel::simulateNetwork(Graph& graph,
		DoubleMatrix& covariates, const DoubleVector& beta) {
}

}
