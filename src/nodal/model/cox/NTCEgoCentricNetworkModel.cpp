/*
 * NTCEgoCentricNetworkModel.cpp
 *
 *  Created on: May 27, 2010
 *      Author: duyvu
 */

#include "NTCEgoCentricNetworkModel.h"

#include "util/ranker.h"

namespace ndip {

NTCEgoCentricNetworkModel::NTCEgoCentricNetworkModel() {
	// TODO Auto-generated constructor stub

}

NTCEgoCentricNetworkModel::~NTCEgoCentricNetworkModel() {
	// TODO Auto-generated destructor stub
}

double NTCEgoCentricNetworkModel::computeLogLikelihood(
		const DoubleVector& beta, bool isUsingLBF) {
	// update this assert when time covariates such as age, gap or other nodal
	// covariates are added.
	assert(beta.size() == data->numOfNodalNetworkStatistics);

	if (isUsingLBF) {
		return computeLogLikelihoodWithLBF(beta);
	} else {
		return computeLogLikelihoodWithoutLBF(beta);
	}
}

double NTCEgoCentricNetworkModel::computeKappaOfNodalUpdateMap(int cacheIndex,
		const DoubleVector& beta) {
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
		result += exp(nodeStatSum);
	}
	return result;
}

double NTCEgoCentricNetworkModel::computeKappaOfComingVertices(int cacheIndex,
		const DoubleVector& beta) {
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
		result += exp(nodeStatSum);
	}
	return result;
}

void NTCEgoCentricNetworkModel::computeKappaDerivativesOfNodalUpdateMap(
		int cacheIndex, const DoubleVector& beta, double& kappa,
		DoubleVector& deltaKappa, DoubleMatrix& squaredDeltaKappa,
		bool resetKappas) {

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

		double temp = exp(nodeStatSum);

		// kappa
		kappa += temp;

		// score kappa
		for (unsigned int h = 0; h < data->numOfNodalNetworkStatistics; h++)
			deltaKappa(h)
					+= data->currentNodalNetworkStatistics(updatedNode, h)
							* temp;
		// information kappa
		for (unsigned int h = 0; h < data->numOfNodalNetworkStatistics; h++)
			for (unsigned int g = 0; g <= h; g++)
				squaredDeltaKappa(h, g) += data->currentNodalNetworkStatistics(
						updatedNode, h) * data->currentNodalNetworkStatistics(
						updatedNode, g) * temp;
	}

}

void NTCEgoCentricNetworkModel::computeKappaDerivativesOfComingVertices(
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

		// get the node
		unsigned int updatedNode = *it;

		// sum over statistics of the contacted nodeB
		double nodeStatSum = 0;
		for (unsigned int h = 0; h < data->numOfNodalNetworkStatistics; h++)
			nodeStatSum += beta(h) * data->currentNodalNetworkStatistics(
					updatedNode, h);

		double temp = exp(nodeStatSum);

		// kappa
		kappa += temp;

		// score kappa
		for (unsigned int h = 0; h < data->numOfNodalNetworkStatistics; h++)
			deltaKappa(h)
					+= data->currentNodalNetworkStatistics(updatedNode, h)
							* temp;
		// information kappa
		for (unsigned int h = 0; h < data->numOfNodalNetworkStatistics; h++)
			for (unsigned int g = 0; g <= h; g++)
				squaredDeltaKappa(h, g) += data->currentNodalNetworkStatistics(
						updatedNode, h) * data->currentNodalNetworkStatistics(
						updatedNode, g) * temp;
	}

}

double NTCEgoCentricNetworkModel::computeLogLikelihoodWithLBF(
		const DoubleVector& beta) {
	double logLLH = 0;

	// start rolling edges
	data->start();

	// the term that will be forwarded between two iterations
	double kappa = 0;

	unsigned int e = data->beginningEdgeIndex;

	// START - kappa(1)

	unsigned int nodeB = edgeJoiningTimes[e].first.second;
	double currentEdgeTime = edgeJoiningTimes[e].second;

	// sum over statistics of the contacted nodeB
	for (unsigned int statIdx = 0; statIdx < data->numOfNodalNetworkStatistics; statIdx++)
		logLLH += beta(statIdx) * data->currentNodalNetworkStatistics(nodeB,
				statIdx);

	for (unsigned int i = 0; i < numOfVertices; i++) {
		// only consider those nodes are already in the network
		if ((nodeJoiningTimes[i] < currentEdgeTime)) {
			// sum over statistics of the node i
			double nodeStatSum = 0;
			for (unsigned int statIdx = 0; statIdx
					< data->numOfNodalNetworkStatistics; statIdx++)
				nodeStatSum += beta(statIdx)
						* data->currentNodalNetworkStatistics(i, statIdx);
			kappa += exp(nodeStatSum);
		}
	}

	//cout << "kappa(" << e << "): " << kappa << endl;

	logLLH -= log(kappa);

	// compute the lookahead term
	double lookaheadTerm = computeKappaOfNodalUpdateMap(
			data->getCurrentCacheIndex(), beta);

	// Update statistics and add the edge to the current graph
	data->updateNodalNetworkStatistics(edgeJoiningTimes[e]);
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

		// Here we update kappa by adding the difference
		double forwardTerm = 0;
		forwardTerm += computeKappaOfNodalUpdateMap(
				data->getCurrentCacheIndex() - 1, beta);
		forwardTerm += computeKappaOfComingVertices(
				data->getCurrentCacheIndex() - 1, beta);
		kappa = kappa - lookaheadTerm + forwardTerm;

		//cout << "kappa(" << e << "): " << kappa << endl;

		logLLH -= log(kappa);

		lookaheadTerm = computeKappaOfNodalUpdateMap(
				data->getCurrentCacheIndex(), beta);

		// Update statistics and add the edge to the current graph
		data->updateNodalNetworkStatistics(edgeJoiningTimes[e]);
	}
	// END - kappa(k)

	// finish rolling edges
	data->finish();
	return logLLH;
}

double NTCEgoCentricNetworkModel::computeLogLikelihoodWithoutLBF(
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
	for (unsigned int e = data->beginningEdgeIndex; e < numOfEdges; e++) {
		int nodeB = edgeJoiningTimes[e].first.second;
		double currentEdgeTime = edgeJoiningTimes[e].second;

		// sum over statistics of the contacted nodeB
		for (unsigned int statIdx = 0; statIdx
				< data->numOfNodalNetworkStatistics; statIdx++)
			logLLH += beta(statIdx) * data->currentNodalNetworkStatistics(
					nodeB, statIdx);

		double kappa = 0;
		for (unsigned int i = 0; i < numOfVertices; i++) {
			// only consider those nodes are already in the network
			if (exposureIndicators[i]
					|| (nodeJoiningTimes[i] < currentEdgeTime)) {
				// in case the exposure of i has not been set up, we set it
				// this happens when i joins the network in [previousEdgeTime, currentEdgeTime)
				if (!exposureIndicators[i])
					exposureIndicators[i] = true;
				// sum over statistics of the node i
				double nodeStatSum = 0;
				for (unsigned int statIdx = 0; statIdx
						< data->numOfNodalNetworkStatistics; statIdx++)
					nodeStatSum += beta(statIdx)
							* data->currentNodalNetworkStatistics(i, statIdx);
				kappa += exp(nodeStatSum);
			}
		}

		//cout << "kappa(" << e << "): " << kappa << endl;

		logLLH -= log(kappa);

		// Update statistics and add the edge to the current graph
		data->updateNodalNetworkStatistics(edgeJoiningTimes[e]);
	}

	// finish rolling edges
	data->finish();

	return logLLH;
}

void NTCEgoCentricNetworkModel::computeScoreVectorAndInformationMatrix(
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

	// START - kappa(1)

	unsigned int nodeB = edgeJoiningTimes[e].first.second;
	double currentEdgeTime = edgeJoiningTimes[e].second;

	// update u by statistics of the contacted nodeB
	for (unsigned int h = 0; h < data->numOfNodalNetworkStatistics; h++)
		u(h) += data->currentNodalNetworkStatistics(nodeB, h);

	// collect kappa, deltaKappa, and squaredDeltaKappa
	for (unsigned int i = 0; i < numOfVertices; i++) {
		// only consider those nodes are already in the network
		if ((nodeJoiningTimes[i] < currentEdgeTime)) {
			// sum over statistics of the node i
			double nodeStatSum = 0;
			for (unsigned int h = 0; h < data->numOfNodalNetworkStatistics; h++)
				nodeStatSum += beta(h) * data->currentNodalNetworkStatistics(i,
						h);

			double temp = exp(nodeStatSum);

			// kappa
			kappa += temp;

			// deltaKappa
			for (unsigned int h = 0; h < data->numOfNodalNetworkStatistics; h++)
				deltaKappa(h) += data->currentNodalNetworkStatistics(i, h)
						* temp;

			// squaredDelataKappa
			for (unsigned int h = 0; h < data->numOfNodalNetworkStatistics; h++)
				for (unsigned int g = 0; g <= h; g++)
					squaredDeltaKappa(h, g)
							+= data->currentNodalNetworkStatistics(i, h)
									* data->currentNodalNetworkStatistics(i, g)
									* temp;
		}
	}

	// update u
	for (unsigned int h = 0; h < data->numOfNodalNetworkStatistics; h++) {
		u(h) -= deltaKappa(h) / kappa;
	}

	// update I
	for (unsigned int h = 0; h < data->numOfNodalNetworkStatistics; h++)
		for (unsigned int g = 0; g <= h; g++) {
			I(h, g) += (squaredDeltaKappa(h, g) / kappa - (deltaKappa(h)
					/ kappa) * (deltaKappa(g) / kappa));

			//cout << "I(h,g) " << I(h, g) << endl;
			//cout << "kappa " << kappa << endl;
			//cout << "deltaKappa(h) " << deltaKappa(h) << endl;
			//cout << "deltaKappa(g) " << deltaKappa(g) << endl;
		}

	// compute the lookahead terms
	double LA_Kappa;
	DoubleVector LA_DeltaKappa(data->numOfNodalNetworkStatistics);
	DoubleMatrix LA_SquaredDeltaKappa(data->numOfNodalNetworkStatistics,
			data->numOfNodalNetworkStatistics);
	computeKappaDerivativesOfNodalUpdateMap(data->getCurrentCacheIndex(), beta,
			LA_Kappa, LA_DeltaKappa, LA_SquaredDeltaKappa, true);

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

		// compute the forward terms
		double FW_Kappa;
		DoubleVector FW_DeltaKappa(data->numOfNodalNetworkStatistics);
		DoubleMatrix FW_SquaredDeltaKappa(data->numOfNodalNetworkStatistics,
				data->numOfNodalNetworkStatistics);
		computeKappaDerivativesOfNodalUpdateMap(data->getCurrentCacheIndex()
				- 1, beta, FW_Kappa, FW_DeltaKappa, FW_SquaredDeltaKappa, true);
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

		//cout << "kappa(" << e << "): " << kappa << endl;

		// update u
		for (unsigned int h = 0; h < data->numOfNodalNetworkStatistics; h++) {
			u(h) -= deltaKappa(h) / kappa;
		}

		// update I
		for (unsigned int h = 0; h < data->numOfNodalNetworkStatistics; h++)
			for (unsigned int g = 0; g <= h; g++) {
				I(h, g) += (squaredDeltaKappa(h, g) / kappa - (deltaKappa(h)
						/ kappa) * (deltaKappa(g) / kappa));

				//cout << "I(h,g) " << I(h, g) << endl;
				//cout << "kappa " << kappa << endl;
				//cout << "deltaKappa(h) " << deltaKappa(h) << endl;
				//cout << "deltaKappa(g) " << deltaKappa(g) << endl;
			}

		computeKappaDerivativesOfNodalUpdateMap(data->getCurrentCacheIndex(),
				beta, LA_Kappa, LA_DeltaKappa, LA_SquaredDeltaKappa, true);

		// Update statistics and add the edge to the current graph
		data->updateNodalNetworkStatistics(edgeJoiningTimes[e]);
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

double NTCEgoCentricNetworkModel::computeLogLikelihoodAndScoreVector(
		DoubleVector& u, const DoubleVector& beta) {
	double logLLH = 0;

	return logLLH;
}

void NTCEgoCentricNetworkModel::getNodalCovariateMatrix(
		DoubleMatrix& covariates, double time) {

}

void NTCEgoCentricNetworkModel::computeCummulativeBaselineHazard(
		VectorOfCummulativePoints& cumBH, const DoubleVector& beta) {

	cumBH.resize(numOfEdges - data->beginningEdgeIndex, false);

	// start rolling edges
	data->start();

	// the term that will be forwarded between two iterations
	double kappa = 0;

	unsigned int e = data->beginningEdgeIndex;

	// START - kappa(1)
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
			kappa += exp(nodeStatSum);
		}
	}

	//cout << "kappa(" << e << "): " << kappa << endl;
	cumBH(0) = make_pair(currentEdgeTime, 1.0 / kappa);

	// compute the lookahead term
	double lookaheadTerm = computeKappaOfNodalUpdateMap(
			data->getCurrentCacheIndex(), beta);

	// Update statistics and add the edge to the current graph
	data->updateNodalNetworkStatistics(edgeJoiningTimes[e]);
	e++;

	// END - kappa(1)

	// START - kappa(k) for k = 2,...,M
	for (int index = 1; e < numOfEdges; e++, index++) {
		double currentEdgeTime = edgeJoiningTimes[e].second;
		// Here we update kappa by adding the difference
		double forwardTerm = 0;
		forwardTerm += computeKappaOfNodalUpdateMap(
				data->getCurrentCacheIndex() - 1, beta);
		forwardTerm += computeKappaOfComingVertices(
				data->getCurrentCacheIndex() - 1, beta);
		kappa = kappa - lookaheadTerm + forwardTerm;

		//cout << "kappa(" << e << "): " << kappa << endl;
		cumBH(index) = make_pair(currentEdgeTime, cumBH(index - 1).second + 1.0
				/ kappa);

		lookaheadTerm = computeKappaOfNodalUpdateMap(
				data->getCurrentCacheIndex(), beta);

		// Update statistics and add the edge to the current graph
		data->updateNodalNetworkStatistics(edgeJoiningTimes[e]);
	}
	// END - kappa(k)

	// finish rolling edges
	data->finish();
}

void NTCEgoCentricNetworkModel::computeCoxSnellResiduals(
		VectorOfCummulativePoints& csResiduals,
		const DoubleVector& beta) {
}

void NTCEgoCentricNetworkModel::computeMartingaleResiduals(
		VectorOfCummulativePoints& martingaleResiduals,
		VectorOfCummulativePoints& timeMartingaleResiduals,
		const DoubleVector& beta) {

}

void NTCEgoCentricNetworkModel::computeSchoenfeldResiduals(
		DoubleMatrix& schoenfeldResiduals, const DoubleVector& beta,
		bool isPlusBeta) {

}

void NTCEgoCentricNetworkModel::getTrainingRanks(DoubleVector& ranks,
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

void NTCEgoCentricNetworkModel::simulateNetwork(Graph& graph,
		DoubleMatrix& covariates, const DoubleVector& beta) {

}

void NTCEgoCentricNetworkModel::simulateNetworkWithBothEnds(Graph& graph,
		DoubleMatrix& covariates, const DoubleVector& betaOut,
		const DoubleVector& betaIn) {

}

}
