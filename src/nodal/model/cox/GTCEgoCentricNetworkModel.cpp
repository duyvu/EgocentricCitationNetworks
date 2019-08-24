/*
 * GTCEgoCentricNetworkModel.cpp
 *
 *  Created on: Jun 1, 2010
 *      Author: duyvu
 */

#include "GTCEgoCentricNetworkModel.h"

#include "util/ranker.h"

#include <cfloat>

namespace ndip {

GTCEgoCentricNetworkModel::GTCEgoCentricNetworkModel() {
	// TODO Auto-generated constructor stub

}

GTCEgoCentricNetworkModel::~GTCEgoCentricNetworkModel() {
	// TODO Auto-generated destructor stub
}

double GTCEgoCentricNetworkModel::computeLogLikelihood(
		const DoubleVector& beta, bool isUsingLBF) {
	if (isUsingLBF) {
		return computeLogLikelihoodWithLBF(beta);
	} else {
		return computeLogLikelihoodWithoutLBF(beta);
	}
}

double GTCEgoCentricNetworkModel::computeKappaOfNodalUpdateMap(int cacheIndex,
		const DoubleVector& beta, double currentEdgeTime) {
	double result = 0;
	for (NodalUpdateMap::iterator it =
			data->vectorOfNodalUpdateMaps(cacheIndex).begin(); it
			!= data->vectorOfNodalUpdateMaps(cacheIndex).end(); ++it) {
		unsigned int updatedNode = it->first;
		// sum over statistics of the contacted nodeB
		double nodeStatSum = 0;
		for (unsigned int h = 0; h < data->numOfNodalNetworkStatistics; h++)
			nodeStatSum += beta(h) * data->currentNodalNetworkStatistics(
					updatedNode, h);
		nodeStatSum += beta(gapTimeBetaIndex) * (currentEdgeTime
				- data->currentRecentNodeContactedTimes[updatedNode]);
		result += exp(nodeStatSum);
	}
	return result;
}

double GTCEgoCentricNetworkModel::computeKappaOfComingVertices(int cacheIndex,
		const DoubleVector& beta, double currentEdgeTime) {
	double result = 0;
	for (ListOfVertices::iterator it =
			data->vectorOfComingVertices(cacheIndex).begin(); it
			!= data->vectorOfComingVertices(cacheIndex).end(); ++it) {
		unsigned int updatedNode = *it;
		// sum over statistics of the contacted nodeB
		double nodeStatSum = 0;
		for (unsigned int h = 0; h < data->numOfNodalNetworkStatistics; h++)
			nodeStatSum += beta(h) * data->currentNodalNetworkStatistics(
					updatedNode, h);
		nodeStatSum += beta(gapTimeBetaIndex) * (currentEdgeTime
				- data->currentRecentNodeContactedTimes[updatedNode]);
		result += exp(nodeStatSum);
	}
	return result;
}

double GTCEgoCentricNetworkModel::computeLogLikelihoodWithLBF(
		const DoubleVector& beta) {
	double logLLH = 0;

	// start rolling edges
	data->start();

	// the term that will be forwarded between two iterations
	double kappa = 0;

	// START - kappa(1)
	unsigned int e = data->beginningEdgeIndex;
	unsigned int nodeB = edgeJoiningTimes[e].first.second;
	double currentEdgeTime = edgeJoiningTimes[e].second;

	// sum over statistics of the contacted nodeB
	for (unsigned int statIdx = 0; statIdx < data->numOfNodalNetworkStatistics; statIdx++)
		logLLH += beta(statIdx) * data->currentNodalNetworkStatistics(nodeB,
				statIdx);
	logLLH += beta(gapTimeBetaIndex) * (currentEdgeTime
			- data->currentRecentNodeContactedTimes[nodeB]);

	for (unsigned int i = 0; i < numOfVertices; i++) {
		// only consider those nodes are already in the network
		if ((nodeJoiningTimes[i] < currentEdgeTime)) {
			// sum over statistics of the node i
			double nodeStatSum = 0;
			for (unsigned int statIdx = 0; statIdx
					< data->numOfNodalNetworkStatistics; statIdx++)
				nodeStatSum += beta(statIdx)
						* data->currentNodalNetworkStatistics(i, statIdx);
			nodeStatSum += beta(gapTimeBetaIndex) * (currentEdgeTime
					- data->currentRecentNodeContactedTimes[i]);
			kappa += exp(nodeStatSum);
		}
	}

	//cout << "kappa(" << e << "): " << kappa << endl;

	logLLH -= log(kappa);

	// compute the lookahead term
	double lookaheadTerm = computeKappaOfNodalUpdateMap(
			data->getCurrentCacheIndex(), beta, currentEdgeTime);

	// Update statistics and add the edge to the current graph
	data->updateNodalNetworkStatistics(edgeJoiningTimes[e]);
	// update time
	double previousEdgeTime = currentEdgeTime;
	e++;

	// END - kappa(1)

	// START - kappa(k) for k = 2,...,M
	for (; e < numOfEdges; e++) {
		int nodeB = edgeJoiningTimes[e].first.second;
		double currentEdgeTime = edgeJoiningTimes[e].second;

		// sum over statistics of the contacted nodeB
		for (unsigned int statIdx = 0; statIdx
				< data->numOfNodalNetworkStatistics; statIdx++)
			logLLH += beta(statIdx) * data->currentNodalNetworkStatistics(
					nodeB, statIdx);
		logLLH += beta(gapTimeBetaIndex) * (currentEdgeTime
				- data->currentRecentNodeContactedTimes[nodeB]);

		// Here we update kappa by adding the difference
		double forwardTerm = 0;
		forwardTerm += computeKappaOfNodalUpdateMap(
				data->getCurrentCacheIndex() - 1, beta, currentEdgeTime);
		forwardTerm += computeKappaOfComingVertices(
				data->getCurrentCacheIndex() - 1, beta, currentEdgeTime);
		double scaler = exp(beta(gapTimeBetaIndex) * (currentEdgeTime
				- previousEdgeTime));
		kappa = (kappa - lookaheadTerm) * scaler + forwardTerm;

		//cout << "kappa(" << e << "): " << kappa << endl;

		logLLH -= log(kappa);

		lookaheadTerm = computeKappaOfNodalUpdateMap(
				data->getCurrentCacheIndex(), beta, currentEdgeTime);

		// update statistics and add the edge to the current graph
		data->updateNodalNetworkStatistics(edgeJoiningTimes[e]);
		// update time
		previousEdgeTime = currentEdgeTime;
	}
	// END - kappa(k)

	// finish rolling edges
	data->finish();
	return logLLH;
}

double GTCEgoCentricNetworkModel::computeLogLikelihoodWithoutLBF(
		const DoubleVector& beta) {
	double logLLH = 0;

	return logLLH;
}

double GTCEgoCentricNetworkModel::computeLogLikelihoodAndScoreVector(
		DoubleVector& u, const DoubleVector& beta) {
	double logLLH = 0;

	return logLLH;
}

void GTCEgoCentricNetworkModel::computeKappaDerivativesOfNodalUpdateMap(
		int cacheIndex, const DoubleVector& beta, double& kappa,
		DoubleVector& deltaKappa, DoubleMatrix& squaredDeltaKappa,
		double currentEdgeTime, bool resetKappas) {

	if (resetKappas) {
		kappa = 0;
		deltaKappa.clear();
		squaredDeltaKappa.clear();
	}

	for (NodalUpdateMap::iterator it =
			data->vectorOfNodalUpdateMaps(cacheIndex).begin(); it
			!= data->vectorOfNodalUpdateMaps(cacheIndex).end(); ++it) {

		// get the node
		unsigned int updatedNode = it->first;

		// sum over statistics of the contacted nodeB
		double nodeStatSum = 0;
		for (unsigned int h = 0; h < data->numOfNodalNetworkStatistics; h++)
			nodeStatSum += beta(h) * data->currentNodalNetworkStatistics(
					updatedNode, h);
		nodeStatSum += beta(gapTimeBetaIndex) * (currentEdgeTime
				- data->currentRecentNodeContactedTimes[updatedNode]);

		double temp = exp(nodeStatSum);

		// kappa
		kappa += temp;

		// score kappa
		for (unsigned int h = 0; h < data->numOfNodalNetworkStatistics; h++)
			deltaKappa(h)
					+= data->currentNodalNetworkStatistics(updatedNode, h)
							* temp;
		deltaKappa(gapTimeBetaIndex) += (currentEdgeTime
				- data->currentRecentNodeContactedTimes[updatedNode]) * temp;

		// information kappa
		for (unsigned int h = 0; h < data->numOfNodalNetworkStatistics; h++) {
			for (unsigned int g = 0; g <= h; g++)
				squaredDeltaKappa(h, g) += data->currentNodalNetworkStatistics(
						updatedNode, h) * data->currentNodalNetworkStatistics(
						updatedNode, g) * temp;
			squaredDeltaKappa(gapTimeBetaIndex, h)
					+= data->currentNodalNetworkStatistics(updatedNode, h)
							* (currentEdgeTime
									- data->currentRecentNodeContactedTimes[updatedNode])
							* temp;
		}
		squaredDeltaKappa(gapTimeBetaIndex, gapTimeBetaIndex)
				+= (currentEdgeTime
						- data->currentRecentNodeContactedTimes[updatedNode])
						* (currentEdgeTime
								- data->currentRecentNodeContactedTimes[updatedNode])
						* temp;
	}
}

void GTCEgoCentricNetworkModel::computeKappaDerivativesOfComingVertices(
		int cacheIndex, const DoubleVector& beta, double& kappa,
		DoubleVector& deltaKappa, DoubleMatrix& squaredDeltaKappa,
		double currentEdgeTime, bool resetKappas) {

	if (resetKappas) {
		kappa = 0;
		deltaKappa.clear();
		squaredDeltaKappa.clear();
	}

	for (ListOfVertices::iterator it =
			data->vectorOfComingVertices(cacheIndex).begin(); it
			!= data->vectorOfComingVertices(cacheIndex).end(); ++it) {

		// get the node
		unsigned int updatedNode = *it;

		// sum over statistics of the contacted nodeB
		double nodeStatSum = 0;
		for (unsigned int h = 0; h < data->numOfNodalNetworkStatistics; h++)
			nodeStatSum += beta(h) * data->currentNodalNetworkStatistics(
					updatedNode, h);
		nodeStatSum += beta(gapTimeBetaIndex) * (currentEdgeTime
				- data->currentRecentNodeContactedTimes[updatedNode]);

		double temp = exp(nodeStatSum);

		// kappa
		kappa += temp;

		// score kappa
		for (unsigned int h = 0; h < data->numOfNodalNetworkStatistics; h++)
			deltaKappa(h)
					+= data->currentNodalNetworkStatistics(updatedNode, h)
							* temp;
		deltaKappa(gapTimeBetaIndex) += (currentEdgeTime
				- data->currentRecentNodeContactedTimes[updatedNode]) * temp;

		// information kappa
		for (unsigned int h = 0; h < data->numOfNodalNetworkStatistics; h++) {
			for (unsigned int g = 0; g <= h; g++)
				squaredDeltaKappa(h, g) += data->currentNodalNetworkStatistics(
						updatedNode, h) * data->currentNodalNetworkStatistics(
						updatedNode, g) * temp;
			squaredDeltaKappa(gapTimeBetaIndex, h)
					+= data->currentNodalNetworkStatistics(updatedNode, h)
							* (currentEdgeTime
									- data->currentRecentNodeContactedTimes[updatedNode])
							* temp;
		}
		squaredDeltaKappa(gapTimeBetaIndex, gapTimeBetaIndex)
				+= (currentEdgeTime
						- data->currentRecentNodeContactedTimes[updatedNode])
						* (currentEdgeTime
								- data->currentRecentNodeContactedTimes[updatedNode])
						* temp;
	}

}

void GTCEgoCentricNetworkModel::updateUandI(DoubleVector& u, DoubleMatrix& I,
		double kappa, const DoubleVector& deltaKappa,
		const DoubleMatrix& squaredDeltaKappa) {
	unsigned int numOfCovariates = data->numOfNodalNetworkStatistics + 1;
	// update u
	for (unsigned int h = 0; h < numOfCovariates; h++) {
		u(h) -= deltaKappa(h) / kappa;
	}
	// update I
	for (unsigned int h = 0; h < numOfCovariates; h++)
		for (unsigned int g = 0; g <= h; g++) {
			I(h, g) += (squaredDeltaKappa(h, g) / kappa - (deltaKappa(h)
					/ kappa) * (deltaKappa(g) / kappa));
		}
}

void GTCEgoCentricNetworkModel::computeScoreVectorAndInformationMatrix(
		DoubleVector& u, DoubleMatrix& I, const DoubleVector& beta) {

	unsigned int numOfCovariates = data->numOfNodalNetworkStatistics + 1;

	u.clear();
	I.clear();

	// start rolling edges
	data->start();

	// the terms that will be forwarded between two iterations
	double kappa = 0;
	DoubleVector deltaKappa(numOfCovariates);
	deltaKappa.clear();
	DoubleMatrix squaredDeltaKappa(numOfCovariates, numOfCovariates);
	squaredDeltaKappa.clear();

	unsigned int e = data->beginningEdgeIndex;

	// START - kappa(1)

	unsigned int nodeB = edgeJoiningTimes[e].first.second;
	double currentEdgeTime = edgeJoiningTimes[e].second;

	// update u by statistics of the contacted nodeB
	for (unsigned int h = 0; h < data->numOfNodalNetworkStatistics; h++)
		u(h) += data->currentNodalNetworkStatistics(nodeB, h);
	u(gapTimeBetaIndex) += (currentEdgeTime
			- data->currentRecentNodeContactedTimes[nodeB]);

	// collect kappa, deltaKappa, and squaredDeltaKappa
	for (unsigned int i = 0; i < numOfVertices; i++) {
		// only consider those nodes are already in the network
		if ((nodeJoiningTimes[i] < currentEdgeTime)) {
			// sum over statistics of the node i
			double nodeStatSum = 0;
			for (unsigned int h = 0; h < data->numOfNodalNetworkStatistics; h++)
				nodeStatSum += beta(h) * data->currentNodalNetworkStatistics(i,
						h);
			nodeStatSum += beta(gapTimeBetaIndex) * (currentEdgeTime
					- data->currentRecentNodeContactedTimes[i]);

			double temp = exp(nodeStatSum);

			// kappa
			kappa += temp;

			// deltaKappa
			for (unsigned int h = 0; h < data->numOfNodalNetworkStatistics; h++)
				deltaKappa(h) += data->currentNodalNetworkStatistics(i, h)
						* temp;
			deltaKappa(gapTimeBetaIndex) += (currentEdgeTime
					- data->currentRecentNodeContactedTimes[i]) * temp;

			// squaredDelataKappa
			for (unsigned int h = 0; h < data->numOfNodalNetworkStatistics; h++) {
				for (unsigned int g = 0; g <= h; g++)
					squaredDeltaKappa(h, g)
							+= data->currentNodalNetworkStatistics(i, h)
									* data->currentNodalNetworkStatistics(i, g)
									* temp;
				squaredDeltaKappa(gapTimeBetaIndex, h)
						+= data->currentNodalNetworkStatistics(i, h)
								* (currentEdgeTime
										- data->currentRecentNodeContactedTimes[i])
								* temp;
			}
			squaredDeltaKappa(gapTimeBetaIndex, gapTimeBetaIndex)
					+= (currentEdgeTime
							- data->currentRecentNodeContactedTimes[i])
							* (currentEdgeTime
									- data->currentRecentNodeContactedTimes[i])
							* temp;
		}
	}

	// update the score vector and the Hessian matrix
	updateUandI(u, I, kappa, deltaKappa, squaredDeltaKappa);

	// compute the lookahead terms
	double LA_Kappa;
	DoubleVector LA_DeltaKappa(numOfCovariates);
	DoubleMatrix LA_SquaredDeltaKappa(numOfCovariates, numOfCovariates);
	computeKappaDerivativesOfNodalUpdateMap(data->getCurrentCacheIndex(), beta,
			LA_Kappa, LA_DeltaKappa, LA_SquaredDeltaKappa, currentEdgeTime,
			true);

	// Update statistics and add the edge to the current graph
	data->updateNodalNetworkStatistics(edgeJoiningTimes[e]);
	e++;

	// END - kappa(1)

	// START - kappa(k) for k = 2,...,M
	for (; e < numOfEdges; e++) {
		int nodeB = edgeJoiningTimes[e].first.second;
		double currentEdgeTime = edgeJoiningTimes[e].second;

		// update u by statistics of the contacted nodeB
		for (unsigned int h = 0; h < data->numOfNodalNetworkStatistics; h++)
			u(h) += data->currentNodalNetworkStatistics(nodeB, h);
		u(gapTimeBetaIndex) += (currentEdgeTime
				- data->currentRecentNodeContactedTimes[nodeB]);

		// compute the forward terms
		double FW_Kappa;
		DoubleVector FW_DeltaKappa(numOfCovariates);
		DoubleMatrix FW_SquaredDeltaKappa(numOfCovariates, numOfCovariates);
		computeKappaDerivativesOfNodalUpdateMap(data->getCurrentCacheIndex()
				- 1, beta, FW_Kappa, FW_DeltaKappa, FW_SquaredDeltaKappa,
				currentEdgeTime, true);
		computeKappaDerivativesOfComingVertices(data->getCurrentCacheIndex()
				- 1, beta, FW_Kappa, FW_DeltaKappa, FW_SquaredDeltaKappa,
				currentEdgeTime, false);

		// here we update kappa, deltaKappa, squaredDeltaKappa by adding the difference
		double previousEdgeTime = edgeJoiningTimes[e - 1].second;
		double scaler = exp(beta(gapTimeBetaIndex) * (currentEdgeTime
				- previousEdgeTime));
		kappa = (kappa - LA_Kappa) * scaler + FW_Kappa;
		for (unsigned int h = 0; h < data->numOfNodalNetworkStatistics; h++)
			deltaKappa(h) = (deltaKappa(h) - LA_DeltaKappa(h)) * scaler
					+ FW_DeltaKappa(h);
		deltaKappa(gapTimeBetaIndex) = (deltaKappa(gapTimeBetaIndex)
				- LA_DeltaKappa(gapTimeBetaIndex)) * scaler + (currentEdgeTime
				- previousEdgeTime) * (kappa - LA_Kappa) * scaler
				+ FW_DeltaKappa(gapTimeBetaIndex);
		for (unsigned int h = 0; h < data->numOfNodalNetworkStatistics; h++) {
			for (unsigned int g = 0; g <= h; g++)
				squaredDeltaKappa(h, g) = (squaredDeltaKappa(h, g)
						- LA_SquaredDeltaKappa(h, g)) * scaler
						+ FW_SquaredDeltaKappa(h, g);
			squaredDeltaKappa(gapTimeBetaIndex, h) = (squaredDeltaKappa(
					gapTimeBetaIndex, h) - LA_SquaredDeltaKappa(
					gapTimeBetaIndex, h)) * scaler + (currentEdgeTime
					- previousEdgeTime) * (deltaKappa(h) - LA_DeltaKappa(h))
					* scaler + FW_SquaredDeltaKappa(gapTimeBetaIndex, h);
		}
		squaredDeltaKappa(gapTimeBetaIndex, gapTimeBetaIndex)
				= (squaredDeltaKappa(gapTimeBetaIndex, gapTimeBetaIndex)
						- LA_SquaredDeltaKappa(gapTimeBetaIndex,
								gapTimeBetaIndex)) * scaler + (currentEdgeTime
						- previousEdgeTime) * (deltaKappa(gapTimeBetaIndex)
						- LA_DeltaKappa(gapTimeBetaIndex)) * scaler
						+ (currentEdgeTime - previousEdgeTime)
								* (currentEdgeTime - previousEdgeTime) * (kappa
								- LA_Kappa) * scaler + FW_SquaredDeltaKappa(
						gapTimeBetaIndex, gapTimeBetaIndex);

		//cout << "kappa(" << e << "): " << kappa << endl;

		// update the score vector U and the Hessian matrix I
		updateUandI(u, I, kappa, deltaKappa, squaredDeltaKappa);

		computeKappaDerivativesOfNodalUpdateMap(data->getCurrentCacheIndex(),
				beta, LA_Kappa, LA_DeltaKappa, LA_SquaredDeltaKappa,
				currentEdgeTime, true);

		// Update statistics and add the edge to the current graph
		data->updateNodalNetworkStatistics(edgeJoiningTimes[e]);
	}
	// END - kappa(k)

	// finish rolling edges
	data->finish();

	// reflex the information matrix since it is symmetric.
	for (unsigned int h = 0; h < numOfCovariates; h++)
		for (unsigned int g = 0; g < h; g++) {
			I(g, h) = I(h, g);
		}
}

void GTCEgoCentricNetworkModel::getNodalCovariateMatrix(
		DoubleMatrix& covariates, double time) {
	covariates.resize(numOfVertices, data->numOfNodalNetworkStatistics + 1,
			false);
	// start rolling edges
	data->start();

	double currentEdgeTime = 0;
	for (unsigned int e = data->beginningEdgeIndex; e < numOfEdges; e++) {
		currentEdgeTime = edgeJoiningTimes[e].second;
		if (currentEdgeTime < time)
			data->updateNodalNetworkStatistics(edgeJoiningTimes[e]);
	}

	for (unsigned int i = 0; i < numOfVertices; i++) {
		for (unsigned int k = 0; k < data->numOfNodalNetworkStatistics; k++)
			covariates(i, k) = data->currentNodalNetworkStatistics(i, k);
		covariates(i, data->numOfNodalNetworkStatistics) = currentEdgeTime
				- data->currentRecentNodeContactedTimes[i];
	}

	// finish rolling edges
	data->finish();
}

void GTCEgoCentricNetworkModel::computeCummulativeBaselineHazard(
		VectorOfCummulativePoints& cumBH, const DoubleVector& beta) {

	cumBH.resize(numOfEdges - data->beginningEdgeIndex, false);

	// start rolling edges
	data->start();

	// the term that will be forwarded between two iterations
	double kappa = 0;

	// START - kappa(1)
	unsigned int e = data->beginningEdgeIndex;
	double currentEdgeTime = edgeJoiningTimes[e].second;
	for (unsigned int i = 0; i < numOfVertices; i++) {
		// only consider those nodes are already in the network
		if ((nodeJoiningTimes[i] < currentEdgeTime)) {
			// sum over statistics of the node i
			double nodeStatSum = 0;
			for (unsigned int statIdx = 0; statIdx
					< data->numOfNodalNetworkStatistics; statIdx++)
				nodeStatSum += beta(statIdx)
						* data->currentNodalNetworkStatistics(i, statIdx);
			nodeStatSum += beta(gapTimeBetaIndex) * (currentEdgeTime
					- data->currentRecentNodeContactedTimes[i]);
			kappa += exp(nodeStatSum);
		}
	}

	//cout << "kappa(" << e << "): " << kappa << endl;
	cumBH(0) = make_pair(currentEdgeTime, 1.0 / kappa);

	// compute the lookahead term
	double lookaheadTerm = computeKappaOfNodalUpdateMap(
			data->getCurrentCacheIndex(), beta, currentEdgeTime);

	// Update statistics and add the edge to the current graph
	data->updateNodalNetworkStatistics(edgeJoiningTimes[e]);
	// update time
	double previousEdgeTime = currentEdgeTime;
	e++;

	// END - kappa(1)

	// START - kappa(k) for k = 2,...,M
	for (int index = 1; e < numOfEdges; e++, index++) {
		double currentEdgeTime = edgeJoiningTimes[e].second;
		// Here we update kappa by adding the difference
		double forwardTerm = 0;
		forwardTerm += computeKappaOfNodalUpdateMap(
				data->getCurrentCacheIndex() - 1, beta, currentEdgeTime);
		forwardTerm += computeKappaOfComingVertices(
				data->getCurrentCacheIndex() - 1, beta, currentEdgeTime);
		double scaler = exp(beta(gapTimeBetaIndex) * (currentEdgeTime
				- previousEdgeTime));
		kappa = (kappa - lookaheadTerm) * scaler + forwardTerm;

		//cout << "kappa(" << e << "): " << kappa << endl;
		cumBH(index) = make_pair(currentEdgeTime, cumBH(index - 1).second + 1.0
				/ kappa);

		lookaheadTerm = computeKappaOfNodalUpdateMap(
				data->getCurrentCacheIndex(), beta, currentEdgeTime);

		// update statistics and add the edge to the current graph
		data->updateNodalNetworkStatistics(edgeJoiningTimes[e]);
		// update time
		previousEdgeTime = currentEdgeTime;
	}
	// END - kappa(k)

	// finish rolling edges
	data->finish();
}

void GTCEgoCentricNetworkModel::computeCoxSnellResiduals(
		VectorOfCummulativePoints& csResiduals,
		const DoubleVector& beta) {

	VectorOfCummulativePoints cumBH =
			VectorOfCummulativePoints(numOfEdges
					- data->beginningEdgeIndex);
	csResiduals.resize(numOfEdges - data->beginningEdgeIndex, false);

	// start rolling edges
	data->start();

	// the term that will be forwarded between two iterations
	double kappa = 0;

	// START - kappa(1)
	unsigned int e = data->beginningEdgeIndex;
	unsigned int nodeB = edgeJoiningTimes[e].first.second;
	double currentEdgeTime = edgeJoiningTimes[e].second;

	// sum over statistics of the contacted nodeB
	double nodeBR = 0;
	for (unsigned int statIdx = 0; statIdx < data->numOfNodalNetworkStatistics; statIdx++)
		nodeBR += beta(statIdx) * data->currentNodalNetworkStatistics(nodeB,
				statIdx);
	nodeBR += beta(gapTimeBetaIndex) * (currentEdgeTime
			- data->currentRecentNodeContactedTimes[nodeB]);
	nodeBR = exp(nodeBR);

	for (unsigned int i = 0; i < numOfVertices; i++) {
		// only consider those nodes are already in the network
		if ((nodeJoiningTimes[i] < currentEdgeTime)) {
			// sum over statistics of the node i
			double nodeStatSum = 0;
			for (unsigned int statIdx = 0; statIdx
					< data->numOfNodalNetworkStatistics; statIdx++)
				nodeStatSum += beta(statIdx)
						* data->currentNodalNetworkStatistics(i, statIdx);
			nodeStatSum += beta(gapTimeBetaIndex) * (currentEdgeTime
					- data->currentRecentNodeContactedTimes[i]);
			kappa += exp(nodeStatSum);
		}
	}

	//cout << "kappa(" << e << "): " << kappa << endl;
	cumBH(0) = make_pair(currentEdgeTime, 1.0 / kappa);
	//csResiduals(0) = make_pair(cumBH(0).second * nodeBR, currentEdgeTime);
	unsigned int popIndex = data->getFirstOrderPopularityIndex();
	csResiduals(0) = make_pair(cumBH(0).second * nodeBR,
			data->currentNodalNetworkStatistics(nodeB, popIndex));

	// compute the lookahead term
	double lookaheadTerm = computeKappaOfNodalUpdateMap(
			data->getCurrentCacheIndex(), beta, currentEdgeTime);

	// Update statistics and add the edge to the current graph
	data->updateNodalNetworkStatistics(edgeJoiningTimes[e]);
	// update time
	double previousEdgeTime = currentEdgeTime;
	e++;

	// END - kappa(1)

	// START - kappa(k) for k = 2,...,M
	for (int index = 1; e < numOfEdges; e++, index++) {
		int nodeB = edgeJoiningTimes[e].first.second;
		double currentEdgeTime = edgeJoiningTimes[e].second;

		// sum over statistics of the contacted nodeB
		double nodeBR = 0;
		for (unsigned int statIdx = 0; statIdx
				< data->numOfNodalNetworkStatistics; statIdx++)
			nodeBR += beta(statIdx) * data->currentNodalNetworkStatistics(
					nodeB, statIdx);
		nodeBR += beta(gapTimeBetaIndex) * (currentEdgeTime
				- data->currentRecentNodeContactedTimes[nodeB]);
		nodeBR = exp(nodeBR);

		// Here we update kappa by adding the difference
		double forwardTerm = 0;
		forwardTerm += computeKappaOfNodalUpdateMap(
				data->getCurrentCacheIndex() - 1, beta, currentEdgeTime);
		forwardTerm += computeKappaOfComingVertices(
				data->getCurrentCacheIndex() - 1, beta, currentEdgeTime);
		double scaler = exp(beta(gapTimeBetaIndex) * (currentEdgeTime
				- previousEdgeTime));
		kappa = (kappa - lookaheadTerm) * scaler + forwardTerm;

		//cout << "kappa(" << e << "): " << kappa << endl;
		cumBH(index) = make_pair(currentEdgeTime, cumBH(index - 1).second + 1.0
				/ kappa);
		csResiduals(index) = make_pair(nodeBR / kappa,
				data->currentNodalNetworkStatistics(nodeB, popIndex));

		lookaheadTerm = computeKappaOfNodalUpdateMap(
				data->getCurrentCacheIndex(), beta, currentEdgeTime);

		// update statistics and add the edge to the current graph
		data->updateNodalNetworkStatistics(edgeJoiningTimes[e]);
		// update time
		previousEdgeTime = currentEdgeTime;
	}
	// END - kappa(k)

	// finish rolling edges
	data->finish();
}

void GTCEgoCentricNetworkModel::computeMartingaleResiduals(
		VectorOfCummulativePoints& martingaleResiduals,
		VectorOfCummulativePoints& timeMartingaleResiduals,
		const DoubleVector& beta) {

	martingaleResiduals.resize(numOfVertices, false);
	timeMartingaleResiduals.resize(numOfEdges - data->beginningEdgeIndex, false);
	unsigned int popIndex = data->getFirstOrderPopularityIndex();

	// start rolling edges
	data->start();

	DoubleVector beginningNumberOfEvents = DoubleVector(numOfVertices);
	for (unsigned int i = 0; i < numOfVertices; i++)
		beginningNumberOfEvents(i) = data->currentNodalNetworkStatistics(i,
				popIndex);

	DoubleVector terms = DoubleVector(numOfVertices);
	DoubleVector accums = DoubleVector(numOfVertices);
	terms.clear();
	accums.clear();
	for (unsigned int e = data->beginningEdgeIndex, index = 0; e < numOfEdges; e++, index++) {
		if (index % 1000 == 0)
			cout << "Working up to " << index << endl;

		double currentEdgeTime = edgeJoiningTimes[e].second;

		double kappa = 0;
		for (unsigned int i = 0; i < numOfVertices; i++) {
			// only consider those nodes are already in the network
			if ((nodeJoiningTimes[i] < currentEdgeTime)) {
				// sum over statistics of the node i
				double nodeStatSum = 0;
				for (unsigned int statIdx = 0; statIdx
						< data->numOfNodalNetworkStatistics; statIdx++)
					nodeStatSum += beta(statIdx)
							* data->currentNodalNetworkStatistics(i, statIdx);
				nodeStatSum += beta(gapTimeBetaIndex) * (currentEdgeTime
						- data->currentRecentNodeContactedTimes[i]);
				terms(i) = exp(nodeStatSum);
				kappa += terms(i);
			}
		}

		list<double> values = list<double> ();
		if (e == (numOfEdges - 1)) {
			for (unsigned int i = 0; i < numOfVertices; i++) {
				if ((nodeJoiningTimes[i] < currentEdgeTime)) {
					accums(i) += terms(i) / kappa;
					double numOfEvents = data->currentNodalNetworkStatistics(i,
							popIndex) - beginningNumberOfEvents(i);
					double residual = numOfEvents - accums(i);
					values.push_back(residual);
					martingaleResiduals(i) = make_pair(numOfEvents, residual);
				}
			}
		} else {
			for (unsigned int i = 0; i < numOfVertices; i++) {
				if ((nodeJoiningTimes[i] < currentEdgeTime)) {
					accums(i) += terms(i) / kappa;
					double numOfEvents = data->currentNodalNetworkStatistics(i,
							popIndex) - beginningNumberOfEvents(i);
					double residual = numOfEvents - accums(i);
					values.push_back(residual);
				}
			}
		}

		double mean;
		double std;
		computeMeanStandardDeviation(mean, std, values);
		timeMartingaleResiduals(index) = make_pair(mean, std);
		//	cout << "e = " << e << endl;
		//	cout << "mean = " << mean << endl;
		//	cout << "std = " << std << endl;

		// update statistics and add the edge to the current graph
		data->updateNodalNetworkStatistics(edgeJoiningTimes[e]);
	}
	// END - kappa(k)

	// finish rolling edges
	data->finish();
}

void GTCEgoCentricNetworkModel::computeSchoenfeldResiduals(
		DoubleMatrix& schoenfeldResiduals, const DoubleVector& beta,
		bool isPlusBeta) {
	unsigned int numOfCovariates = data->numOfNodalNetworkStatistics + 1;
	schoenfeldResiduals.resize(numOfEdges - data->beginningEdgeIndex,
			numOfCovariates, false);
	DoubleVector residuals = DoubleVector(numOfCovariates);

	// start rolling edges
	data->start();

	// the terms that will be forwarded between two iterations
	double kappa = 0;
	DoubleVector deltaKappa(numOfCovariates);
	deltaKappa.clear();
	DoubleMatrix squaredDeltaKappa(numOfCovariates, numOfCovariates);
	squaredDeltaKappa.clear();

	unsigned int e = data->beginningEdgeIndex;

	// START - kappa(1)

	unsigned int nodeB = edgeJoiningTimes[e].first.second;
	double currentEdgeTime = edgeJoiningTimes[e].second;

	// update u by statistics of the contacted nodeB
	residuals.clear();
	for (unsigned int h = 0; h < data->numOfNodalNetworkStatistics; h++)
		residuals(h) = data->currentNodalNetworkStatistics(nodeB, h);
	residuals(gapTimeBetaIndex) = (currentEdgeTime
			- data->currentRecentNodeContactedTimes[nodeB]);

	// collect kappa, deltaKappa, and squaredDeltaKappa
	for (unsigned int i = 0; i < numOfVertices; i++) {
		// only consider those nodes are already in the network
		if ((nodeJoiningTimes[i] < currentEdgeTime)) {
			// sum over statistics of the node i
			double nodeStatSum = 0;
			for (unsigned int h = 0; h < data->numOfNodalNetworkStatistics; h++)
				nodeStatSum += beta(h) * data->currentNodalNetworkStatistics(i,
						h);
			nodeStatSum += beta(gapTimeBetaIndex) * (currentEdgeTime
					- data->currentRecentNodeContactedTimes[i]);

			double temp = exp(nodeStatSum);

			// kappa
			kappa += temp;

			// deltaKappa
			for (unsigned int h = 0; h < data->numOfNodalNetworkStatistics; h++)
				deltaKappa(h) += data->currentNodalNetworkStatistics(i, h)
						* temp;
			deltaKappa(gapTimeBetaIndex) += (currentEdgeTime
					- data->currentRecentNodeContactedTimes[i]) * temp;

			// squaredDelataKappa
			for (unsigned int h = 0; h < data->numOfNodalNetworkStatistics; h++) {
				for (unsigned int g = 0; g <= h; g++)
					squaredDeltaKappa(h, g)
							+= data->currentNodalNetworkStatistics(i, h)
									* data->currentNodalNetworkStatistics(i, g)
									* temp;
				squaredDeltaKappa(gapTimeBetaIndex, h)
						+= data->currentNodalNetworkStatistics(i, h)
								* (currentEdgeTime
										- data->currentRecentNodeContactedTimes[i])
								* temp;
			}
			squaredDeltaKappa(gapTimeBetaIndex, gapTimeBetaIndex)
					+= (currentEdgeTime
							- data->currentRecentNodeContactedTimes[i])
							* (currentEdgeTime
									- data->currentRecentNodeContactedTimes[i])
							* temp;
		}
	}

	// compute the residuals
	computeSchoenfeldResiduals(residuals, kappa, deltaKappa, squaredDeltaKappa);
	// add to the matrix result
	unsigned int index = 0;
	for (unsigned int k = 0; k < numOfCovariates; k++)
		schoenfeldResiduals(index, k) = residuals(k);
	index++;

	// compute the lookahead terms
	double LA_Kappa;
	DoubleVector LA_DeltaKappa(numOfCovariates);
	DoubleMatrix LA_SquaredDeltaKappa(numOfCovariates, numOfCovariates);
	computeKappaDerivativesOfNodalUpdateMap(data->getCurrentCacheIndex(), beta,
			LA_Kappa, LA_DeltaKappa, LA_SquaredDeltaKappa, currentEdgeTime,
			true);

	// Update statistics and add the edge to the current graph
	data->updateNodalNetworkStatistics(edgeJoiningTimes[e]);
	e++;

	// END - kappa(1)

	// START - kappa(k) for k = 2,...,M
	for (; e < numOfEdges; e++) {
		int nodeB = edgeJoiningTimes[e].first.second;
		double currentEdgeTime = edgeJoiningTimes[e].second;

		// update u by statistics of the contacted nodeB
		residuals.clear();
		for (unsigned int h = 0; h < data->numOfNodalNetworkStatistics; h++)
			residuals(h) = data->currentNodalNetworkStatistics(nodeB, h);
		residuals(gapTimeBetaIndex) = (currentEdgeTime
				- data->currentRecentNodeContactedTimes[nodeB]);

		// compute the forward terms
		double FW_Kappa;
		DoubleVector FW_DeltaKappa(numOfCovariates);
		DoubleMatrix FW_SquaredDeltaKappa(numOfCovariates, numOfCovariates);
		computeKappaDerivativesOfNodalUpdateMap(data->getCurrentCacheIndex()
				- 1, beta, FW_Kappa, FW_DeltaKappa, FW_SquaredDeltaKappa,
				currentEdgeTime, true);
		computeKappaDerivativesOfComingVertices(data->getCurrentCacheIndex()
				- 1, beta, FW_Kappa, FW_DeltaKappa, FW_SquaredDeltaKappa,
				currentEdgeTime, false);

		// here we update kappa, deltaKappa, squaredDeltaKappa by adding the difference
		double previousEdgeTime = edgeJoiningTimes[e - 1].second;
		double scaler = exp(beta(gapTimeBetaIndex) * (currentEdgeTime
				- previousEdgeTime));
		kappa = (kappa - LA_Kappa) * scaler + FW_Kappa;
		for (unsigned int h = 0; h < data->numOfNodalNetworkStatistics; h++)
			deltaKappa(h) = (deltaKappa(h) - LA_DeltaKappa(h)) * scaler
					+ FW_DeltaKappa(h);
		deltaKappa(gapTimeBetaIndex) = (deltaKappa(gapTimeBetaIndex)
				- LA_DeltaKappa(gapTimeBetaIndex)) * scaler + (currentEdgeTime
				- previousEdgeTime) * (kappa - LA_Kappa) * scaler
				+ FW_DeltaKappa(gapTimeBetaIndex);
		for (unsigned int h = 0; h < data->numOfNodalNetworkStatistics; h++) {
			for (unsigned int g = 0; g <= h; g++)
				squaredDeltaKappa(h, g) = (squaredDeltaKappa(h, g)
						- LA_SquaredDeltaKappa(h, g)) * scaler
						+ FW_SquaredDeltaKappa(h, g);
			squaredDeltaKappa(gapTimeBetaIndex, h) = (squaredDeltaKappa(
					gapTimeBetaIndex, h) - LA_SquaredDeltaKappa(
					gapTimeBetaIndex, h)) * scaler + (currentEdgeTime
					- previousEdgeTime) * (deltaKappa(h) - LA_DeltaKappa(h))
					* scaler + FW_SquaredDeltaKappa(gapTimeBetaIndex, h);
		}
		squaredDeltaKappa(gapTimeBetaIndex, gapTimeBetaIndex)
				= (squaredDeltaKappa(gapTimeBetaIndex, gapTimeBetaIndex)
						- LA_SquaredDeltaKappa(gapTimeBetaIndex,
								gapTimeBetaIndex)) * scaler + (currentEdgeTime
						- previousEdgeTime) * (deltaKappa(gapTimeBetaIndex)
						- LA_DeltaKappa(gapTimeBetaIndex)) * scaler
						+ (currentEdgeTime - previousEdgeTime)
								* (currentEdgeTime - previousEdgeTime) * (kappa
								- LA_Kappa) * scaler + FW_SquaredDeltaKappa(
						gapTimeBetaIndex, gapTimeBetaIndex);

		//cout << "kappa(" << e << "): " << kappa << endl;

		// compute the residuals
		computeSchoenfeldResiduals(residuals, kappa, deltaKappa,
				squaredDeltaKappa);
		for (unsigned int k = 0; k < numOfCovariates; k++)
			schoenfeldResiduals(index, k) = residuals(k);
		index++;

		computeKappaDerivativesOfNodalUpdateMap(data->getCurrentCacheIndex(),
				beta, LA_Kappa, LA_DeltaKappa, LA_SquaredDeltaKappa,
				currentEdgeTime, true);

		// Update statistics and add the edge to the current graph
		data->updateNodalNetworkStatistics(edgeJoiningTimes[e]);
	}
	// END - kappa(k)

	// finish rolling edges
	data->finish();

	if (isPlusBeta)
		for (unsigned int i = 0; i < schoenfeldResiduals.size1(); i++)
			for (unsigned int k = 0; k < schoenfeldResiduals.size2(); k++)
				schoenfeldResiduals(i, k) += beta(k);
}

void GTCEgoCentricNetworkModel::computeSchoenfeldResiduals(
		DoubleVector& residuals, double kappa, const DoubleVector& deltaKappa,
		const DoubleMatrix& squaredDeltaKappa) {
	// update residuals
	unsigned int numOfCovariates = residuals.size();
	for (unsigned int h = 0; h < numOfCovariates; h++) {
		residuals(h) -= deltaKappa(h) / kappa;
	}
	// compute V
	DoubleMatrix V = DoubleMatrix(numOfCovariates, numOfCovariates);
	for (unsigned int h = 0; h < numOfCovariates; h++)
		for (unsigned int g = 0; g <= h; g++) {
			V(h, g) = (squaredDeltaKappa(h, g) / kappa
					- (deltaKappa(h) / kappa) * (deltaKappa(g) / kappa));
			V(g, h) = V(h, g);
		}
	// V^{-1} * u, standardized residuals
	boost::numeric::ublas::permutation_matrix<> pm(numOfCovariates);
	lu_factorize(V, pm);
	lu_substitute(V, pm, residuals);
}

void GTCEgoCentricNetworkModel::getTrainingRanks(DoubleVector& ranks,
		const DoubleVector& beta, const std::string& tieMethod) {
	ranks.resize(numOfEdges - data->beginningEdgeIndex, false);

	// start rolling edges
	data->start();

	for (unsigned int e = data->beginningEdgeIndex, index = 0; e < numOfEdges; e++, index++) {
		if (index % 1000 == 0)
			cout << "Working up to " << index << endl;

		int nodeB = edgeJoiningTimes[e].first.second;
		double currentEdgeTime = edgeJoiningTimes[e].second;

		std::vector<double> intensities = std::vector<double>(numOfVertices);
		for (unsigned int i = 0; i < numOfVertices; i++) {
			if (nodeJoiningTimes[i] < currentEdgeTime) {
				double nodeStatSum = 0;
				for (unsigned int h = 0; h < data->numOfNodalNetworkStatistics; h++)
					nodeStatSum += beta(h)
							* data->currentNodalNetworkStatistics(i, h);
				nodeStatSum += beta(gapTimeBetaIndex) * (currentEdgeTime
						- data->currentRecentNodeContactedTimes[i]);
				intensities[i] = exp(nodeStatSum);
			} else {
				intensities[i] = 0;
			}
		}

		std::vector<double> currentNodeRanks = std::vector<double>(
				numOfVertices);
		ndip::rankhigh(intensities, currentNodeRanks, tieMethod);
		ranks(index) = currentNodeRanks[nodeB];
		//cout << "node " << nodeB << " has rank " << currentNodeRanks[nodeB]
		//		<< endl;

		// update statistics and add the edge to the current graph
		data->updateNodalNetworkStatistics(edgeJoiningTimes[e]);
	}

	// finish rolling edges
	data->finish();
}

void GTCEgoCentricNetworkModel::simulateNetwork(Graph& graph,
		DoubleMatrix& covariates, const DoubleVector& beta) {
}

void GTCEgoCentricNetworkModel::simulateNetworkWithBothEnds(Graph& graph,
		DoubleMatrix& covariates, const DoubleVector& betaOut,
		const DoubleVector& betaIn) {
}

}
