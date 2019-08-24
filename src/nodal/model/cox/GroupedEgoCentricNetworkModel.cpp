/*
 * GroupedEgoCentricNetworkModel.cpp
 *
 *  Created on: Jan 17, 2011
 *      Author: duyvu
 */

#include "GroupedEgoCentricNetworkModel.h"

#include "nodal/data/GroupedEgoCentricNetworkData.h"
#include "nodal/data/PGroupedEgoCentricNetworkData.h"
#include "nodal/data/P2PTGroupedEgoCentricNetworkData.h"
#include "nodal/data/P2PTRGroupedEgoCentricNetworkData.h"

#include "util/ranker.h"

#include "util/StringTokenizer.h"
#include <iostream>
#include <fstream>
#include <math.h>
#include <float.h>

namespace ndip {

GroupedEgoCentricNetworkModel::GroupedEgoCentricNetworkModel() {
	// TODO Auto-generated constructor stub

}

GroupedEgoCentricNetworkModel::~GroupedEgoCentricNetworkModel() {
	// TODO Auto-generated destructor stub
	if (data != NULL)
		delete data;
}

// read data of nodes and edges from files
void GroupedEgoCentricNetworkModel::setNetworkData(
		char* nodeJoiningTimesFileName, char* edgeEventsFileName) {

	// Enter the node joining times
	cout << "The node joining times file is " << nodeJoiningTimesFileName
			<< endl;
	string line;
	ifstream nodeJoiningTimesFile(nodeJoiningTimesFileName);
	if (nodeJoiningTimesFile.is_open()) {
		// The first line is the number of vertices
		if (!nodeJoiningTimesFile.eof()) {
			getline(nodeJoiningTimesFile, line);
			StringTokenizer strtok = StringTokenizer(line, "\t");
			numOfVertices = strtok.nextIntToken();
			nodeJoiningTimes = ExposureTimeVectorType(numOfVertices);
		}
		int currentIndex = 0;
		while (!nodeJoiningTimesFile.eof()) {
			getline(nodeJoiningTimesFile, line);
			StringTokenizer strtok = StringTokenizer(line, "\t");
			strtok.nextIntToken(); // the IDs must be from 0 to N - 1; therefore we do not to read them in
			double time = strtok.nextDoubleToken(); // by days
			nodeJoiningTimes[currentIndex++] = time;
		}
		nodeJoiningTimesFile.close();
	} else
		cout << "Unable to open file " << nodeJoiningTimesFileName << endl;

	// Enter the edge joining time
	cout << "The edge events file is " << edgeEventsFileName << endl;
	ifstream edgeEventsFile(edgeEventsFileName);
	if (edgeEventsFile.is_open()) {
		// The first line is the number of discrete times
		if (!edgeEventsFile.eof()) {
			getline(edgeEventsFile, line);
			StringTokenizer strtok = StringTokenizer(line, "\t");
			numOfDistinctEdgeEvents = strtok.nextIntToken();
			vectorOfDiscreteTimeVectorOfEdgeEvents
					= VectorOfDiscreteTimeVectorOfEdgeEvents(
							numOfDistinctEdgeEvents);
		}
		cout << "Number of event times is " << numOfDistinctEdgeEvents << endl;
		unsigned int totalCitingNodes = 0;
		unsigned int totalCitedNodes = 0;
		for (unsigned int k = 0; k < numOfDistinctEdgeEvents; k++) {
			getline(edgeEventsFile, line);
			StringTokenizer strtok = StringTokenizer(line, "\t");
			double eventTime = strtok.nextDoubleToken();
			int numOfCitingNodes = strtok.nextIntToken();
			totalCitingNodes += numOfCitingNodes;
			VectorOfEdgeEvents vectorOfEdgeEvents = VectorOfEdgeEvents(
					numOfCitingNodes);
			for (unsigned int i = 0; i < numOfCitingNodes; i++) {
				getline(edgeEventsFile, line);
				StringTokenizer strtok = StringTokenizer(line, "\t");
				int citingNode = strtok.nextIntToken();
				NodeVector nodeVector = NodeVector();
				while (strtok.hasMoreTokens())
					nodeVector.push_back(strtok.nextIntToken());
				vectorOfEdgeEvents[i] = make_pair(citingNode, nodeVector);
				totalCitedNodes += nodeVector.size();
			}
			vectorOfDiscreteTimeVectorOfEdgeEvents[k] = make_pair(eventTime,
					vectorOfEdgeEvents);
		}
		cout << "Number of citing papers is " << totalCitingNodes << endl;
		cout << "Number of cited papers is " << totalCitedNodes << endl;
		edgeEventsFile.close();
	} else
		cout << "Unable to open file " << edgeEventsFileName << endl;
}

// choose which ego centric model to use
// for a set of network statistics of interest,
// we will build an corresponding ego centric network data
// the likelhood and ML estimation code will be independent with
// the selected ego centric network data
void GroupedEgoCentricNetworkModel::setModel(int _modelType,
		double _observationTimeStart, double _observationTimeEnd,
		const std::vector<string>& listOfNodalDataFiles) {

	observationTimeStart = _observationTimeStart;
	observationTimeEnd = _observationTimeEnd;

	if (data != NULL)
		delete data;

	switch (_modelType) {
	case PREFERENTIAL_ATTACHMENT_EGO_CENTRIC_MODEL:
		cout << "PREFERENTIAL ATTACHMENT EGO CENTRIC MODEL is in use!" << endl;
		data = new PGroupedEgoCentricNetworkData(nodeJoiningTimes,
				vectorOfDiscreteTimeVectorOfEdgeEvents, observationTimeStart,
				observationTimeEnd);
		initializeMore();
		break;
	case P2PT_EGO_CENTRIC_MODEL:
		cout << "P2PT EGO CENTRIC MODEL is in use!" << endl;
		data = new P2PTGroupedEgoCentricNetworkData(nodeJoiningTimes,
				vectorOfDiscreteTimeVectorOfEdgeEvents, observationTimeStart,
				observationTimeEnd);
		initializeMore();
		break;
	case P2PTR_EGO_CENTRIC_MODEL:
		cout << "P2PTR EGO CENTRIC MODEL is in use!" << endl;
		data = new P2PTRGroupedEgoCentricNetworkData(nodeJoiningTimes,
				vectorOfDiscreteTimeVectorOfEdgeEvents, observationTimeStart,
				observationTimeEnd);
		initializeMore();
		break;
	default:
		cout
				<< "BY DEFAULT!!! PREFERENTIAL ATTACHMENT EGO CENTRIC MODEL is in use!"
				<< endl;
		data = new PGroupedEgoCentricNetworkData(nodeJoiningTimes,
				vectorOfDiscreteTimeVectorOfEdgeEvents, observationTimeStart,
				observationTimeEnd);
		initializeMore();
	}
}

void GroupedEgoCentricNetworkModel::getMLE(DoubleVector& beta,
		DoubleMatrix& cov, double logLLHTolerance, int maxNRSteps,
		double scoreTolerance) {

	DoubleVector u(beta.size());
	DoubleMatrix I(beta.size(), beta.size());
	double currentLogLLH = computeLogLikelihood(beta, true);
	cout << "Initial logLLH = " << currentLogLLH << endl;
	for (int k = 0; k < maxNRSteps; k++) {

		cout << "Step " << k << endl;

		// obtain the score vector and the information matrix
		computeScoreVectorAndInformationMatrix(u, I, beta);

		cout << "u: " << u << endl;

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
		double relLogLLHImprovement = fabs((currentLogLLH - newLogLLH)
				/ currentLogLLH);
		currentLogLLH = newLogLLH;
		cout << "LogLLH: " << currentLogLLH << endl;

		if (fabs(maxError) < scoreTolerance)
			break;
		if (relLogLLHImprovement < logLLHTolerance)
			break;
	}
}

void GroupedEgoCentricNetworkModel::prepareExposureIndicators(
		bool* exposureIndicators) {
	for (unsigned int i = 0; i < numOfVertices; i++) {
		// a node is exposed if it joined before the observation time
		// not include this time
		exposureIndicators[i] = (nodeJoiningTimes[i] < observationTimeStart);
	}
}

void GroupedEgoCentricNetworkModel::getRecencyBaselineRanks(
		ofstream& byEdgesOutputFile, ofstream& byEventsOutputFile,
		ofstream& accumlativeOutputFile, const std::string& tieMethod) {
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
				intensities[i] = currentCitingTime
						- data->currentRecentNodeCitedTimes[i];
			} else
				intensities[i] = DBL_MAX;
		}
		std::vector<double> currentNodeRanks = std::vector<double>(
				numOfVertices);
		ndip::rank(intensities, currentNodeRanks, tieMethod);

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
