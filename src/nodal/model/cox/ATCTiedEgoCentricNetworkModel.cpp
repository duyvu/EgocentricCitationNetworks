// underflow issue

/*
 * ATCTiedEgoCentricNetworkModel.cpp
 *
 *  Created on: Oct 28, 2010
 *      Author: duyvu
 */

#include "ATCTiedEgoCentricNetworkModel.h"

namespace ndip {

ATCTiedEgoCentricNetworkModel::ATCTiedEgoCentricNetworkModel() {
	// TODO Auto-generated constructor stub

}

ATCTiedEgoCentricNetworkModel::~ATCTiedEgoCentricNetworkModel() {
	// TODO Auto-generated destructor stub
}

double ATCTiedEgoCentricNetworkModel::computeLogLikelihood(
		const DoubleVector& beta, bool isUsingLBF) {
	if (isUsingLBF) {
		return computeLogLikelihoodWithLBF(beta);
	} else {
		return computeLogLikelihoodWithoutLBF(beta);
	}
}

void ATCTiedEgoCentricNetworkModel::computeLogLikelihoodOverEvents(
		DoubleVector& cummulativeLogLLH, const DoubleVector& beta,
		bool isUsingLBF) {

}

double ATCTiedEgoCentricNetworkModel::computeLogLikelihoodWithLBF(
		const DoubleVector& beta) {

	double logLLH = 0;

	// start rolling edges
	data->start();

	// the term that will be forwarded between two iterations
	double kappa = 0;

	unsigned int e = data->beginningEdgeIndex;

	// START - kappa(1)

	unsigned int citingNode = vectorOfEdgeEvents[e].first;
	double currentCitingTime = nodeJoiningTimes[citingNode];
	unsigned int numOfTiedEvents = vectorOfEdgeEvents[e].second.size();
	for (unsigned int i = 0; i < numOfTiedEvents; i++) {
		unsigned int citedNode = vectorOfEdgeEvents[e].second[i];
		logLLH += computeSumNodalStatistics(citedNode, beta, currentCitingTime);
	}

	for (unsigned int i = 0; i < numOfVertices; i++) {
		// only consider those nodes are already in the network
		if ((nodeJoiningTimes[i] < currentCitingTime)) {
			// sum over statistics of the node i
			double nodalStatSum = computeSumNodalStatistics(i, beta,
					currentCitingTime);
			kappa += exp(nodalStatSum);
		}
	}

	//cout << "kappa(" << e << "): " << kappa << endl;
	logLLH -= numOfTiedEvents * log(kappa);

	// compute the lookahead term
	double lookaheadTerm = computeKappaOfNodalUpdateMap(
			data->getCurrentCacheIndex(), beta, currentCitingTime);

	// update statistics and add the edge to the current graph
	data->updateNodalNetworkStatistics(vectorOfEdgeEvents[e]);
	// update time
	double previousCitingTime = currentCitingTime;
	e++;

	// END - kappa(1)

	// START - kappa(k) for k = 2,...,M
	for (; e <= data->endingEdgeIndex; e++) {
		// add those cited nodes
		unsigned int citingNode = vectorOfEdgeEvents[e].first;
		currentCitingTime = nodeJoiningTimes[citingNode];
		unsigned int numOfTiedEvents = vectorOfEdgeEvents[e].second.size();
		for (unsigned int i = 0; i < numOfTiedEvents; i++) {
			unsigned int citedNode = vectorOfEdgeEvents[e].second[i];
			logLLH += computeSumNodalStatistics(citedNode, beta,
					currentCitingTime);
		}

		// Here we update kappa by adding the difference
		double forwardTerm = 0;
		forwardTerm += computeKappaOfNodalUpdateMap(
				data->getCurrentCacheIndex() - 1, beta, currentCitingTime);
		forwardTerm += computeKappaOfComingVertices(
				data->getCurrentCacheIndex() - 1, beta, currentCitingTime);
		double scaler = exp(beta(ageTimeBetaIndex) * (currentCitingTime
				- previousCitingTime));
		kappa = (kappa - lookaheadTerm) * scaler + forwardTerm;

		//cout << "kappa(" << e << "): " << kappa << endl;
		logLLH -= numOfTiedEvents * log(kappa);

		lookaheadTerm = computeKappaOfNodalUpdateMap(
				data->getCurrentCacheIndex(), beta, currentCitingTime);

		// update statistics and add the edge to the current graph
		data->updateNodalNetworkStatistics(vectorOfEdgeEvents[e]);
		// update time
		previousCitingTime = currentCitingTime;

		cout << logLLH << endl;
	}
	// END - kappa(k)

	// finish rolling edges
	data->finish();

	return logLLH;
}

void ATCTiedEgoCentricNetworkModel::computeLogLikelihoodOverEdges(
		ofstream& byEdgesOutputFile, const DoubleVector& beta, bool isUsingLBF) {

}

double ATCTiedEgoCentricNetworkModel::computeSumNodalStatistics(
		unsigned int nodeID, const DoubleVector& beta, double currentCitingTime) {
	// sum over statistics of the cited node
	double nodalStatSum = 0;
	for (unsigned int h = 0; h < data->numOfNodalNetworkStatistics; h++)
		nodalStatSum += beta(h)
				* data->currentNodalNetworkStatistics(nodeID, h);
	nodalStatSum += beta(ageTimeBetaIndex) * (currentCitingTime
			- nodeJoiningTimes[nodeID]);
}

double ATCTiedEgoCentricNetworkModel::computeKappaOfNodalUpdateMap(
		int cacheIndex, const DoubleVector& beta, double currentCitingTime) {
	double result = 0;
	for (NodalUpdateMap::iterator it =
			data->vectorOfNodalUpdateMaps(cacheIndex).begin(); it
			!= data->vectorOfNodalUpdateMaps(cacheIndex).end(); ++it) {
		unsigned int updatedNode = it->first;
		// sum over statistics of the cited node
		double nodalStatSum = computeSumNodalStatistics(updatedNode, beta,
				currentCitingTime);
		result += exp(nodalStatSum);
	}
	return result;
}

double ATCTiedEgoCentricNetworkModel::computeKappaOfComingVertices(
		int cacheIndex, const DoubleVector& beta, double currentCitingTime) {
	double result = 0;
	for (ListOfVertices::iterator it =
			data->vectorOfComingVertices(cacheIndex).begin(); it
			!= data->vectorOfComingVertices(cacheIndex).end(); ++it) {
		unsigned int newNode = *it;
		// sum over statistics of the new node
		double nodalStatSum = computeSumNodalStatistics(newNode, beta,
				currentCitingTime);
		result += exp(nodalStatSum);
	}
	return result;
}

double ATCTiedEgoCentricNetworkModel::computeLogLikelihoodWithoutLBF(
		const DoubleVector& beta) {

}

void ATCTiedEgoCentricNetworkModel::computeScoreVectorAndInformationMatrix(
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

	// add those cited nodes
	unsigned int citingNode = vectorOfEdgeEvents[e].first;
	double currentCitingTime = nodeJoiningTimes[citingNode];
	//cout << "currentCitingTime = " << currentCitingTime << endl;
	unsigned int numOfTiedEvents = vectorOfEdgeEvents[e].second.size();
	for (unsigned int i = 0; i < numOfTiedEvents; i++) {
		unsigned int citedNode = vectorOfEdgeEvents[e].second[i];
		// update u by statistics of the cited node
		for (unsigned int h = 0; h < data->numOfNodalNetworkStatistics; h++)
			u(h) += data->currentNodalNetworkStatistics(citedNode, h);
		u(ageTimeBetaIndex)
				+= (currentCitingTime - nodeJoiningTimes[citedNode]);
		//cout << "citedTime = " << nodeJoiningTimes[citedNode]
		//		<< endl;
	}

	// collect kappa, deltaKappa, and squaredDeltaKappa
	for (unsigned int i = 0; i < numOfVertices; i++) {
		// only consider those nodes are already in the network
		if ((nodeJoiningTimes[i] < currentCitingTime))
			computeNodalDerivatives(i, beta, kappa, deltaKappa,
					squaredDeltaKappa, currentCitingTime);
	}

	// update the score vector and the Hessian matrix
	updateUandI(u, I, kappa, deltaKappa, squaredDeltaKappa, numOfTiedEvents);

	// compute the lookahead terms
	double LA_Kappa;
	DoubleVector LA_DeltaKappa(numOfCovariates);
	DoubleMatrix LA_SquaredDeltaKappa(numOfCovariates, numOfCovariates);
	computeKappaDerivativesOfNodalUpdateMap(data->getCurrentCacheIndex(), beta,
			LA_Kappa, LA_DeltaKappa, LA_SquaredDeltaKappa, currentCitingTime,
			true);

	// Update statistics and add the edge to the current graph
	data->updateNodalNetworkStatistics(vectorOfEdgeEvents[e]);
	double previousCitingTime = currentCitingTime;
	e++;

	// END - kappa(1)

	// START - kappa(k) for k = 2,...,M
	for (; e <= data->endingEdgeIndex; e++) {

		unsigned int citingNode = vectorOfEdgeEvents[e].first;
		currentCitingTime = nodeJoiningTimes[citingNode];
		//cout << "currentCitingTime = " << currentCitingTime << endl;
		unsigned int numOfTiedEvents = vectorOfEdgeEvents[e].second.size();
		for (unsigned int i = 0; i < numOfTiedEvents; i++) {
			unsigned int citedNode = vectorOfEdgeEvents[e].second[i];
			// update u by statistics of the cited node
			for (unsigned int h = 0; h < data->numOfNodalNetworkStatistics; h++)
				u(h) += data->currentNodalNetworkStatistics(citedNode, h);
			u(ageTimeBetaIndex) += (currentCitingTime
					- nodeJoiningTimes[citedNode]);
			//cout << "citedTime = "
			//		<< nodeJoiningTimes[citedNode] << endl;
		}

		// compute the forward terms
		double FW_Kappa;
		DoubleVector FW_DeltaKappa(numOfCovariates);
		DoubleMatrix FW_SquaredDeltaKappa(numOfCovariates, numOfCovariates);
		computeKappaDerivativesOfNodalUpdateMap(data->getCurrentCacheIndex()
				- 1, beta, FW_Kappa, FW_DeltaKappa, FW_SquaredDeltaKappa,
				currentCitingTime, true);
		computeKappaDerivativesOfComingVertices(data->getCurrentCacheIndex()
				- 1, beta, FW_Kappa, FW_DeltaKappa, FW_SquaredDeltaKappa,
				currentCitingTime, false);

		// here we update kappa, deltaKappa, squaredDeltaKappa by adding the difference
		double scaler = exp(beta(ageTimeBetaIndex) * (currentCitingTime
				- previousCitingTime));
		// kappa
		kappa = (kappa - LA_Kappa) * scaler + FW_Kappa;
		// deltaKappa
		for (unsigned int h = 0; h < data->numOfNodalNetworkStatistics; h++)
			deltaKappa(h) = (deltaKappa(h) - LA_DeltaKappa(h)) * scaler
					+ FW_DeltaKappa(h);
		deltaKappa(ageTimeBetaIndex) = (deltaKappa(ageTimeBetaIndex)
				- LA_DeltaKappa(ageTimeBetaIndex)) * scaler
				+ (currentCitingTime - previousCitingTime) * (kappa - LA_Kappa)
						* scaler + FW_DeltaKappa(ageTimeBetaIndex);
		// squaredDeltaKappa
		for (unsigned int h = 0; h < data->numOfNodalNetworkStatistics; h++) {
			for (unsigned int g = 0; g <= h; g++)
				squaredDeltaKappa(h, g) = (squaredDeltaKappa(h, g)
						- LA_SquaredDeltaKappa(h, g)) * scaler
						+ FW_SquaredDeltaKappa(h, g);
			squaredDeltaKappa(ageTimeBetaIndex, h) = (squaredDeltaKappa(
					ageTimeBetaIndex, h) - LA_SquaredDeltaKappa(
					ageTimeBetaIndex, h)) * scaler + (currentCitingTime
					- previousCitingTime) * (deltaKappa(h) - LA_DeltaKappa(h))
					* scaler + FW_SquaredDeltaKappa(ageTimeBetaIndex, h);
		}
		squaredDeltaKappa(ageTimeBetaIndex, ageTimeBetaIndex)
				= (squaredDeltaKappa(ageTimeBetaIndex, ageTimeBetaIndex)
						- LA_SquaredDeltaKappa(ageTimeBetaIndex,
								ageTimeBetaIndex)) * scaler
						+ (currentCitingTime - previousCitingTime)
								* (deltaKappa(ageTimeBetaIndex)
										- LA_DeltaKappa(ageTimeBetaIndex))
								* scaler + (currentCitingTime
						- previousCitingTime) * (currentCitingTime
						- previousCitingTime) * (kappa - LA_Kappa) * scaler
						+ FW_SquaredDeltaKappa(ageTimeBetaIndex,
								ageTimeBetaIndex);

		// update the score vector U and the Hessian matrix I
		updateUandI(u, I, kappa, deltaKappa, squaredDeltaKappa, numOfTiedEvents);

		computeKappaDerivativesOfNodalUpdateMap(data->getCurrentCacheIndex(),
				beta, LA_Kappa, LA_DeltaKappa, LA_SquaredDeltaKappa,
				currentCitingTime, true);

		// Update statistics and add the edge to the current graph
		data->updateNodalNetworkStatistics(vectorOfEdgeEvents[e]);
		previousCitingTime = currentCitingTime;
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

void ATCTiedEgoCentricNetworkModel::computeNodalDerivatives(
		unsigned int nodeID, const DoubleVector& beta, double& kappa,
		DoubleVector& deltaKappa, DoubleMatrix& squaredDeltaKappa,
		double currentCitingTime) {

	// sum over statistics of the contacted nodeB
	double nodalStatSum = computeSumNodalStatistics(nodeID, beta,
			currentCitingTime);

	double temp = exp(nodalStatSum);

	// kappa
	kappa += temp;

	// score kappa
	for (unsigned int h = 0; h < data->numOfNodalNetworkStatistics; h++)
		deltaKappa(h) += data->currentNodalNetworkStatistics(nodeID, h) * temp;
	deltaKappa(ageTimeBetaIndex) += (currentCitingTime
			- nodeJoiningTimes[nodeID]) * temp;

	// information kappa
	for (unsigned int h = 0; h < data->numOfNodalNetworkStatistics; h++) {
		for (unsigned int g = 0; g <= h; g++)
			squaredDeltaKappa(h, g) += data->currentNodalNetworkStatistics(
					nodeID, h) * data->currentNodalNetworkStatistics(nodeID, g)
					* temp;
		squaredDeltaKappa(ageTimeBetaIndex, h)
				+= data->currentNodalNetworkStatistics(nodeID, h)
						* (currentCitingTime - nodeJoiningTimes[nodeID]) * temp;
	}
	squaredDeltaKappa(ageTimeBetaIndex, ageTimeBetaIndex) += (currentCitingTime
			- nodeJoiningTimes[nodeID]) * (currentCitingTime
			- nodeJoiningTimes[nodeID]) * temp;
}

void ATCTiedEgoCentricNetworkModel::computeKappaDerivativesOfNodalUpdateMap(
		int cacheIndex, const DoubleVector& beta, double& kappa,
		DoubleVector& deltaKappa, DoubleMatrix& squaredDeltaKappa,
		double currentCitingTime, bool resetKappas) {

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
		computeNodalDerivatives(updatedNode, beta, kappa, deltaKappa,
				squaredDeltaKappa, currentCitingTime);
	}
}

void ATCTiedEgoCentricNetworkModel::computeKappaDerivativesOfComingVertices(
		int cacheIndex, const DoubleVector& beta, double& kappa,
		DoubleVector& deltaKappa, DoubleMatrix& squaredDeltaKappa,
		double currentCitingTime, bool resetKappas) {

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
				squaredDeltaKappa, currentCitingTime);
	}
}

void ATCTiedEgoCentricNetworkModel::updateUandI(DoubleVector& u,
		DoubleMatrix& I, const double& kappa, const DoubleVector& deltaKappa,
		const DoubleMatrix& squaredDeltaKappa, unsigned int numOfTiedEvents) {
	unsigned int numOfCovariates = data->numOfNodalNetworkStatistics + 1;
	// update u
	for (unsigned int h = 0; h < numOfCovariates; h++) {
		u(h) -= numOfTiedEvents * deltaKappa(h) / kappa;
	}
	// update I
	for (unsigned int h = 0; h < numOfCovariates; h++)
		for (unsigned int g = 0; g <= h; g++) {
			I(h, g) += numOfTiedEvents * (squaredDeltaKappa(h, g) / kappa
					- (deltaKappa(h) / kappa) * (deltaKappa(g) / kappa));
		}
}

double ATCTiedEgoCentricNetworkModel::computeLogLikelihoodAndScoreVector(
		DoubleVector& u, const DoubleVector& beta) {

}

void ATCTiedEgoCentricNetworkModel::computeOnlineGradientAscentMLE(
		DoubleVector& beta, double learningRate) {

}

void ATCTiedEgoCentricNetworkModel::computeCummulativeBaselineHazard(
		VectorOfCummulativePoints& cumBH, const DoubleVector& beta) {

}

void ATCTiedEgoCentricNetworkModel::computeSchoenfeldResiduals(
		DoubleMatrix& schoenfeldResiduals, const DoubleVector& beta,
		bool isPlusBeta) {

}

void ATCTiedEgoCentricNetworkModel::getSumRanks(DoubleVector& sumRanks,
		const DoubleVector& beta, const std::string& tieMethod) {

}

void ATCTiedEgoCentricNetworkModel::getRanks(VectorOfDoubleVector& ranks,
		const DoubleVector& beta, const std::string& tieMethod) {

}

void ATCTiedEgoCentricNetworkModel::simulateNetwork(Graph& graph,
		DoubleMatrix& covariates, const DoubleVector& beta) {

}

}
