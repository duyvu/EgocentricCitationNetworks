/*
 * DenseNTCTiedEgoCentricNetworkModel.cpp
 *
 *  Created on: Jan 13, 2011
 *      Author: duyvu
 */

#include "DenseNTCTiedEgoCentricNetworkModel.h"

#include "util/ranker.h"

#include <iostream>
#include <fstream>

namespace ndip {

DenseNTCTiedEgoCentricNetworkModel::DenseNTCTiedEgoCentricNetworkModel() {
	// TODO Auto-generated constructor stub

}

DenseNTCTiedEgoCentricNetworkModel::~DenseNTCTiedEgoCentricNetworkModel() {
	// TODO Auto-generated destructor stub
}

double DenseNTCTiedEgoCentricNetworkModel::computeLogLikelihood(
		const DoubleVector& beta, bool isUsingLBF) {
	if (isUsingLBF) {
		cout
				<< "This LBF feature is not supported in this class! Full iteration will be used!!!"
				<< endl;
	}
	return computeLogLikelihoodWithoutLBF(beta);
}

double DenseNTCTiedEgoCentricNetworkModel::computeLogLikelihoodWithoutLBF(
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

		if (e % 100 == 0)
			cout << "Working up to " << e << endl;

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
		data->updateNodalNetworkStatistics(e);
	}

	// finish rolling edges
	data->finish();

	return logLLH;
}

void DenseNTCTiedEgoCentricNetworkModel::computeLogLikelihoodOverEdges(
		ofstream& byEdgesOutputFile, const DoubleVector& beta, bool isUsingLBF) {

	if (!isUsingLBF)
		cout << "This feature has not implemented yet" << endl;

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

		if (e % 100 == 0)
			cout << "Working up to " << e << endl;

		// add those cited nodes
		unsigned int citingNode = vectorOfEdgeEvents[e].first;
		double currentCitingTime = nodeJoiningTimes[citingNode];

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

		unsigned int numOfTiedEvents = vectorOfEdgeEvents[e].second.size();
		for (unsigned int i = 0; i < numOfTiedEvents; i++) {
			unsigned int citedNode = vectorOfEdgeEvents[e].second[i];
			byEdgesOutputFile << citingNode << "\t" << citedNode << "\t"
					<< currentCitingTime << "\t" << computeSumNodalStatistics(
					citedNode, beta) - log(kappa) << endl;
		}

		// Update statistics and add the edge to the current graph
		data->updateNodalNetworkStatistics(e);
	}

	// finish rolling edges
	data->finish();
}

double DenseNTCTiedEgoCentricNetworkModel::computeSumNodalStatistics(
		unsigned int nodeID, const DoubleVector& beta) {
	double nodalStatSum = 0;
	for (unsigned int statIdx = 0; statIdx < data->numOfNodalNetworkStatistics; statIdx++)
		nodalStatSum += beta(statIdx) * data->currentNodalNetworkStatistics(
				nodeID, statIdx);
	return nodalStatSum;
}

void DenseNTCTiedEgoCentricNetworkModel::computeLogLikelihoodOverEvents(
		DoubleVector& cummulativeLogLLH, const DoubleVector& beta,
		bool isUsingLBF) {
	if (isUsingLBF) {
		cout << "This feature is not supported in this class!" << endl;
	}

	cummulativeLogLLH.resize(data->endingEdgeIndex - data->beginningEdgeIndex
			+ 1, false);
	cummulativeLogLLH.clear();

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
	for (unsigned int e = data->beginningEdgeIndex, index = 0; e
			<= data->endingEdgeIndex; e++, index++) {

		if (e % 100 == 0)
			cout << "Working up to " << e << endl;

		// add those cited nodes
		unsigned int citingNode = vectorOfEdgeEvents[e].first;
		double currentCitingTime = nodeJoiningTimes[citingNode];
		unsigned int numOfTiedEvents = vectorOfEdgeEvents[e].second.size();
		for (unsigned int i = 0; i < numOfTiedEvents; i++) {
			unsigned int citedNode = vectorOfEdgeEvents[e].second[i];
			cummulativeLogLLH[index] += computeSumNodalStatistics(citedNode,
					beta);
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
		cummulativeLogLLH[index] -= numOfTiedEvents * log(kappa);

		// Update statistics and add the edge to the current graph
		data->updateNodalNetworkStatistics(e);
	}

	// finish rolling edges
	data->finish();
}

void DenseNTCTiedEgoCentricNetworkModel::computeScoreVectorAndInformationMatrix(
		DoubleVector& u, DoubleMatrix& I, const DoubleVector& beta) {

	u.clear();
	I.clear();

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
	double kappa;
	DoubleVector deltaKappa(data->numOfNodalNetworkStatistics);
	DoubleMatrix squaredDeltaKappa(data->numOfNodalNetworkStatistics,
			data->numOfNodalNetworkStatistics);
	for (unsigned int e = data->beginningEdgeIndex; e <= data->endingEdgeIndex; e++) {

		if (e % 10 == 0)
			cout << "Working up to " << e << endl;

		unsigned int citingNode = vectorOfEdgeEvents[e].first;
		double currentCitingTime = nodeJoiningTimes[citingNode];
		unsigned int numOfTiedEvents = vectorOfEdgeEvents[e].second.size();
		for (unsigned int i = 0; i < numOfTiedEvents; i++) {
			unsigned int citedNode = vectorOfEdgeEvents[e].second[i];
			// update u by statistics of the cited node
			for (unsigned int h = 0; h < data->numOfNodalNetworkStatistics; h++)
				u(h) += data->currentNodalNetworkStatistics(citedNode, h);
		}

		// here we update kappa, deltaKappa, squaredDeltaKappa
		kappa = 0;
		deltaKappa.clear();
		squaredDeltaKappa.clear();
		for (unsigned int i = 0; i < numOfVertices; i++) {
			// only consider those nodes are already in the network
			if (exposureIndicators[i] || (nodeJoiningTimes[i]
					< currentCitingTime)) {
				// in case the exposure of i has not been set up, we set it
				// this happens when i joins the network in [previousEdgeTime, currentEdgeTime)
				if (!exposureIndicators[i])
					exposureIndicators[i] = true;
				computeNodalDerivatives(i, beta, kappa, deltaKappa,
						squaredDeltaKappa);
			}
		}

		// update u and I
		updateUandI(u, I, kappa, deltaKappa, squaredDeltaKappa, numOfTiedEvents);

		// Update statistics and add the edge to the current graph
		data->updateNodalNetworkStatistics(e);
	}

	// finish rolling edges
	data->finish();

	// reflex the information matrix since it is symmetric.
	for (unsigned int h = 0; h < data->numOfNodalNetworkStatistics; h++)
		for (unsigned int g = 0; g < h; g++) {
			I(g, h) = I(h, g);
		}
}

void DenseNTCTiedEgoCentricNetworkModel::computeNodalDerivatives(
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

void DenseNTCTiedEgoCentricNetworkModel::updateUandI(DoubleVector& u,
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

double DenseNTCTiedEgoCentricNetworkModel::computeLogLikelihoodAndScoreVector(
		DoubleVector& u, const DoubleVector& beta) {
	return 0;
}

void DenseNTCTiedEgoCentricNetworkModel::computeOnlineGradientAscentMLE(
		DoubleVector& beta, double learningRate) {

}

void DenseNTCTiedEgoCentricNetworkModel::computeCummulativeBaselineHazard(
		VectorOfCummulativePoints& cumBH, const DoubleVector& beta) {

	cumBH.resize(data->endingEdgeIndex - data->beginningEdgeIndex + 1, false);

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
	for (unsigned int e = data->beginningEdgeIndex, index = 0; e
			<= data->endingEdgeIndex; e++, index++) {

		if (e % 100 == 0)
			cout << "Working up to " << e << endl;

		// add those cited nodes
		unsigned int citingNode = vectorOfEdgeEvents[e].first;
		double currentCitingTime = nodeJoiningTimes[citingNode];
		unsigned int numOfTiedEvents = vectorOfEdgeEvents[e].second.size();

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

		if (index > 0)
			cumBH(index) = make_pair(currentCitingTime, cumBH(index - 1).second
					+ numOfTiedEvents / kappa);
		else
			cumBH(index)
					= make_pair(currentCitingTime, numOfTiedEvents / kappa);

		// Update statistics and add the edge to the current graph
		data->updateNodalNetworkStatistics(e);
	}

	// finish rolling edges
	data->finish();
}

void DenseNTCTiedEgoCentricNetworkModel::computeSchoenfeldResiduals(
		DoubleMatrix& schoenfeldResiduals, const DoubleVector& beta,
		bool isPlusBeta) {

}

void DenseNTCTiedEgoCentricNetworkModel::getSumRanks(DoubleVector& sumRanks,
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
		data->updateNodalNetworkStatistics(e);
	}

	// finish rolling edges
	data->finish();
}

void DenseNTCTiedEgoCentricNetworkModel::getRanks(VectorOfDoubleVector& ranks,
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
		data->updateNodalNetworkStatistics(e);
	}

	// finish rolling edges
	data->finish();
}

void DenseNTCTiedEgoCentricNetworkModel::getRanks(ofstream& byEdgesOutputFile,
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
		data->updateNodalNetworkStatistics(e);
	}

	// finish rolling edges
	data->finish();
}

void DenseNTCTiedEgoCentricNetworkModel::simulateNetwork(Graph& graph,
		DoubleMatrix& covariates, const DoubleVector& beta) {

}

}
