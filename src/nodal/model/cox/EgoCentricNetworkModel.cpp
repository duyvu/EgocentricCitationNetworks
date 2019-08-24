/*
 * EgoCentricNetworkModel.cpp
 *
 *  Created on: May 27, 2010
 *      Author: duyvu
 */

#include "EgoCentricNetworkModel.h"
#include "nodal/data/PEgoCentricNetworkData.h"
#include "nodal/data/APEgoCentricNetworkData.h"
#include "nodal/data/TAPEgoCentricNetworkData.h"
#include "nodal/data/NTAPEgoCentricNetworkData.h"
#include "nodal/data/SimTAPEgoCentricNetworkData.h"
#include "nodal/data/SimNTAPEgoCentricNetworkData.h"

#include "util/StringTokenizer.h"
#include <iostream>
#include <fstream>
#include <math.h>

namespace ndip {

EgoCentricNetworkModel::EgoCentricNetworkModel() {
	// TODO Auto-generated constructor stub
	data = NULL;
}

EgoCentricNetworkModel::~EgoCentricNetworkModel() {
	// TODO Auto-generated destructor stub
	if (data != NULL)
		delete data;
}

void EgoCentricNetworkModel::setNetworkData(char* exposureTimeFile,
		char* edgeTimeFile) {
	// Enter the node joining time
	cout << "The node joining time file is " << exposureTimeFile << endl;
	string line;
	ifstream expoFile(exposureTimeFile);
	if (expoFile.is_open()) {
		// The first line is the number of vertices
		if (!expoFile.eof()) {
			getline(expoFile, line);
			StringTokenizer strtok = StringTokenizer(line, "\t");
			nodeJoiningTimes.reserve(strtok.nextIntToken());
		}
		while (!expoFile.eof()) {
			getline(expoFile, line);
			StringTokenizer strtok = StringTokenizer(line, "\t");
			if (strtok.countTokens() == 2) {
				strtok.nextIntToken();
				double time = strtok.nextDoubleToken(); // scale up the time
				nodeJoiningTimes.push_back(time);
			}
		}
		numOfVertices = nodeJoiningTimes.size();
		expoFile.close();
	} else
		cout << "Unable to open file " << exposureTimeFile << endl;

	// Enter the edge joining time
	cout << "The edge joining time file is " << edgeTimeFile << endl;
	ifstream edgeFile(edgeTimeFile);
	if (edgeFile.is_open()) {
		// The first line is the number of edges
		if (!edgeFile.eof()) {
			getline(edgeFile, line);
			StringTokenizer strtok = StringTokenizer(line, "\t");
			edgeJoiningTimes.reserve(strtok.nextIntToken());
		}
		while (!edgeFile.eof()) {
			getline(edgeFile, line);
			StringTokenizer strtok = StringTokenizer(line, "\t");
			if (strtok.countTokens() == 3) {
				int nodeA = strtok.nextIntToken();
				int nodeB = strtok.nextIntToken();
				double time = strtok.nextDoubleToken(); // scale up the time
				edgeJoiningTimes.push_back(EdgeTime(Edge(nodeA, nodeB), time));
			}
		}
		numOfEdges = edgeJoiningTimes.size();
		edgeFile.close();
	} else
		cout << "Unable to open file " << edgeTimeFile << endl;
}

void EgoCentricNetworkModel::setModel(int _dataType, double _observationTime,
		const std::vector<string>& listOfNodalDataFiles) {
	observationTime = _observationTime;
	if (data != NULL)
		delete data;
	switch (_dataType) {
	case PREFERENTIAL_ATTACHMENT_EGO_CENTRIC_MODEL:
		cout << "here we go" << endl;
		data = new PEgoCentricNetworkData(nodeJoiningTimes, edgeJoiningTimes,
				observationTime);
		initializeMore();
		break;
	case ACTIVITY_POPULARITY_EGO_CENTRIC_MODEL:
		data = new APEgoCentricNetworkData(nodeJoiningTimes, edgeJoiningTimes,
				observationTime);
		initializeMore();
		break;
	case TRANSITIVITY_ACTIVITY_POPULARITY_EGO_CENTRIC_MODEL:
		data = new TAPEgoCentricNetworkData(nodeJoiningTimes, edgeJoiningTimes,
				observationTime);
		initializeMore();
		break;
	case NODAL_TRANSITIVITY_ACTIVITY_POPULARITY_EGO_CENTRIC_MODEL:
		data = new NTAPEgoCentricNetworkData(nodeJoiningTimes,
				edgeJoiningTimes, observationTime, listOfNodalDataFiles);
		initializeMore();
		break;
	case SIMULATING_TRANSITIVITY_ACTIVITY_POPULARITY_EGO_CENTRIC_MODEL:
		data = new SimTAPEgoCentricNetworkData(nodeJoiningTimes,
				edgeJoiningTimes, observationTime);
		initializeMore();
		break;
	case SIMULATING_NODAL_TRANSITIVITY_ACTIVITY_POPULARITY_EGO_CENTRIC_MODEL:
		data = new SimNTAPEgoCentricNetworkData(nodeJoiningTimes,
				edgeJoiningTimes, observationTime, listOfNodalDataFiles);
		initializeMore();
		break;
	default:
		data = new APEgoCentricNetworkData(nodeJoiningTimes, edgeJoiningTimes,
				observationTime);
		initializeMore();
	}
}

void EgoCentricNetworkModel::prepareExposureIndicators(bool* exposureIndicators) {
	for (int i = 0; i < numOfVertices; i++) {
		// a node is exposed if it joined before the observation time
		// not include this time
		exposureIndicators[i] = (nodeJoiningTimes[i] < observationTime);
	}
}

void EgoCentricNetworkModel::getMLE(DoubleVector& beta, DoubleMatrix& cov,
		double firstPercentageStepLength, int maxNRSteps, double scoreTolerance) {

	DoubleVector u(beta.size());
	DoubleMatrix I(beta.size(), beta.size());
	double currentLogLLH = computeLogLikelihood(beta, true);
	for (int k = 0; k < maxNRSteps; k++) {

		cout << "Step " << k << endl;

		// obtain the score vector and the information matrix
		computeScoreVectorAndInformationMatrix(u, I, beta);

		cout << "u: " << u << endl;
		//		cout << "I: " << endl;
		//		for (int k = 0; k < cov.size1(); k++) {
		//			for (int h = 0; h < cov.size2(); h++)
		//				cout << I(k, h) << "  ";
		//			cout << endl;
		//		}

		// obtain the covariance estimate
		InvertMatrix(I, cov);

		cout << "Beta Cov: " << endl;
		for (unsigned int k = 0; k < cov.size1(); k++) {
			for (unsigned int h = 0; h < cov.size2(); h++)
				cout << cov(k, h) << "  ";
			cout << endl;
		}

		// check the score vector to stop if all of elements
		// are smaller than the tolerance
		double maxError = fabs(u(0));
		for (unsigned int k = 1; k < u.size(); k++)
			if (maxError <= fabs(u(k)))
				maxError = fabs(u(k));

		boost::numeric::ublas::permutation_matrix<> pm(beta.size());
		lu_factorize(I, pm);
		lu_substitute(I, pm, u);

		// search for a good step length
		double stepLength = 1.0;
		DoubleVector newBeta = DoubleVector(beta.size());
		newBeta = beta + stepLength * u;
		double newLogLLH = computeLogLikelihood(newBeta, true);
		unsigned searchCount = 1;
		if (isnan(newLogLLH) || newLogLLH < currentLogLLH) {
			while (isnan(newLogLLH) || newLogLLH < currentLogLLH) {
				stepLength /= 2.0;
				newBeta = beta + stepLength * u;
				newLogLLH = computeLogLikelihood(newBeta, true);
				searchCount++;
			}
		}

		cout << "Number of function evaluations for step length search "
				<< searchCount << endl;

		// update beta
		beta = newBeta;
		cout << "New Beta Estimates: " << endl;
		for (unsigned int h = 0; h < beta.size(); h++)
			cout << beta(h) << "  ";
		cout << endl;

		// update currentlogLLH
		currentLogLLH = newLogLLH;
		cout << "LogLLH: " << currentLogLLH << endl;

		if (fabs(maxError) < scoreTolerance)
			break;
	}

}

void EgoCentricNetworkModel::computeMeanStandardDeviation(double& mean,
		double& stDeviation, list<double>& values) {
	mean = 0;
	for (list<double>::iterator it = values.begin(); it != values.end(); it++)
		mean += *it;
	mean /= values.size();
	stDeviation = 0;
	for (list<double>::iterator it = values.begin(); it != values.end(); it++)
		stDeviation += (*it - mean) * (*it - mean);
	stDeviation /= (values.size() - 1);
	stDeviation = sqrt(stDeviation);
}

}
