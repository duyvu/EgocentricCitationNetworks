/*
 * TiedEgoCentricNetworkModel.cpp
 *
 *  Created on: Jul 18, 2010
 *      Author: duyvu
 */

#include "TiedEgoCentricNetworkModel.h"
#include "nodal/data/PTiedEgoCentricNetworkData.h"
#include "nodal/data/PGTiedEgoCentricNetworkData.h"
#include "nodal/data/P2PTiedEgoCentricNetworkData.h"
#include "nodal/data/P2PTTiedEgoCentricNetworkData.h"
#include "nodal/data/P2PTRTiedEgoCentricNetworkData.h"
#include "nodal/data/TextualP2PTRTiedEgoCentricNetworkData.h"
#include "nodal/data/LDAP2PTRTiedEgoCentricNetworkData.h"
#include "nodal/data/LDATiedEgoCentricNetworkData.h"
#include "nodal/data/SimP2PTRTiedEgoCentricNetworkData.h"

#include "util/ranker.h"
#include "util/StringTokenizer.h"

#include <iostream>
#include <fstream>
#include <float.h>

namespace ndip {

TiedEgoCentricNetworkModel::TiedEgoCentricNetworkModel() {
	// TODO Auto-generated constructor stub
	data = NULL;
}

TiedEgoCentricNetworkModel::~TiedEgoCentricNetworkModel() {
	// TODO Auto-generated destructor stub
	if (data != NULL) {
		delete data;
		data = NULL;
	}
}

// read data of nodes and edges from files
void TiedEgoCentricNetworkModel::setNetworkData(char* nodeJoiningTimesFileName,
		char* edgeEventsFileName) {

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
			cout << "numOfVertices = " << numOfVertices << endl;
			nodeJoiningTimes = ExposureTimeVectorType(numOfVertices);
		}
		int currentIndex = 0;
		while (!nodeJoiningTimesFile.eof()) {
			getline(nodeJoiningTimesFile, line);
			StringTokenizer strtok = StringTokenizer(line, "\t");
			strtok.nextIntToken(); // the IDs must be from 0 to N - 1; therefore we do not to read them in
			double time = strtok.nextDoubleToken(); // by days
			nodeJoiningTimes.push_back(time);
			currentIndex++;
		}
		nodeJoiningTimesFile.close();
	} else
		cout << "Unable to open file " << nodeJoiningTimesFileName << endl;

	// Enter the edge joining time
	cout << "The edge events file is " << edgeEventsFileName << endl;
	ifstream edgeEventsFile(edgeEventsFileName);
	if (edgeEventsFile.is_open()) {
		// The first line is the number of edge events
		if (!edgeEventsFile.eof()) {
			getline(edgeEventsFile, line);
			StringTokenizer strtok = StringTokenizer(line, "\t");
			numOfDistinctEdgeEvents = strtok.nextIntToken();
			cout << "numOfDistinctEdgeEvents = " << numOfDistinctEdgeEvents
					<< endl;
		}
		int currentIndex = 0;
		while (!edgeEventsFile.eof()) {
			getline(edgeEventsFile, line);
			StringTokenizer strtok = StringTokenizer(line, "\t");
			int citingNode = strtok.nextIntToken();
			NodeVector nodeVector = NodeVector();
			while (strtok.hasMoreTokens())
				nodeVector.push_back(strtok.nextIntToken());
			vectorOfEdgeEvents.push_back(make_pair(citingNode, nodeVector));
			currentIndex++;
		}
		edgeEventsFile.close();
	} else
		cout << "Unable to open file " << edgeEventsFileName << endl;

	//	cout << "vectorOfEdgeEvents size = " << vectorOfEdgeEvents.size() << endl;
	//	cout << "vectorOfEdgeEvents[0].first = " << vectorOfEdgeEvents[0].first
	//			<< endl;
	//	cout << "vectorOfEdgeEvents[25003].first = "
	//			<< vectorOfEdgeEvents[25003].first << endl;

	cout << "Done with setNetworkData(...)" << endl;
}

// choose which ego centric model to use
// for a set of network statistics of interest,
// we will build an corresponding ego centric network data
// the likelhood and ML estimation code will be independent with
// the selected ego centric network data
void TiedEgoCentricNetworkModel::setModel(int _modelType,
		double _observationTimeStart, double _observationTimeEnd,
		const std::vector<string>& listOfNodalDataFiles) {

	observationTimeStart = _observationTimeStart;
	observationTimeEnd = _observationTimeEnd;

	if (data != NULL)
		delete data;

	switch (_modelType) {
	case PREFERENTIAL_ATTACHMENT_EGO_CENTRIC_MODEL:
		cout << "PREFERENTIAL ATTACHMENT EGO CENTRIC MODEL is in use!" << endl;
		data = new PTiedEgoCentricNetworkData(nodeJoiningTimes,
				vectorOfEdgeEvents, observationTimeStart, observationTimeEnd);
		initializeMore();
		break;
	case P2P_EGO_CENTRIC_MODEL:
		cout << "P2P EGO CENTRIC MODEL is in use!" << endl;
		data = new P2PTiedEgoCentricNetworkData(nodeJoiningTimes,
				vectorOfEdgeEvents, observationTimeStart, observationTimeEnd);
		initializeMore();
		break;
	case P2PT_EGO_CENTRIC_MODEL:
		cout << "P2PT EGO CENTRIC MODEL is in use!" << endl;
		data = new P2PTTiedEgoCentricNetworkData(nodeJoiningTimes,
				vectorOfEdgeEvents, observationTimeStart, observationTimeEnd);
		initializeMore();
		break;
	case P2PTR_EGO_CENTRIC_MODEL:
		cout << "P2PTR EGO CENTRIC MODEL is in use!" << endl;
		data = new P2PTRTiedEgoCentricNetworkData(nodeJoiningTimes,
				vectorOfEdgeEvents, observationTimeStart, observationTimeEnd);
		initializeMore();
		break;
	case TextualP2PTR_EGO_CENTRIC_MODEL:
		cout << "Textual P2PTR EGO CENTRIC MODEL is in use!" << endl;
		data = new TextualP2PTRTiedEgoCentricNetworkData(nodeJoiningTimes,
				vectorOfEdgeEvents, observationTimeStart, observationTimeEnd,
				listOfNodalDataFiles);
		initializeMore();
		break;
	case LDA_P2PTR_EGO_CENTRIC_MODEL:
		cout << "LDA P2PTR EGO CENTRIC MODEL is in use!" << endl;
		data = new LDAP2PTRTiedEgoCentricNetworkData(nodeJoiningTimes,
				vectorOfEdgeEvents, observationTimeStart, observationTimeEnd,
				listOfNodalDataFiles);
		initializeMore();
		break;
	case LDA_EGO_CENTRIC_MODEL:
		cout << "LDA EGO CENTRIC MODEL is in use!" << endl;
		data = new LDATiedEgoCentricNetworkData(nodeJoiningTimes,
				vectorOfEdgeEvents, observationTimeStart, observationTimeEnd,
				listOfNodalDataFiles);
		initializeMore();
		break;
	case SIMULATING_P2PTR_EGO_CENTRIC_MODEL:
		cout << "Simulating P2PTR EGO CENTRIC MODEL is in use!" << endl;
		data = new SimP2PTRTiedEgoCentricNetworkData(nodeJoiningTimes,
				vectorOfEdgeEvents, observationTimeStart, observationTimeEnd);
		initializeMore();
		break;
	case PA_GAP_EGO_CENTRIC_MODEL:
		cout << "PA GAP EGO CENTRIC MODEL is in use!" << endl;
		data = new PGTiedEgoCentricNetworkData(nodeJoiningTimes,
				vectorOfEdgeEvents, observationTimeStart, observationTimeEnd);
		initializeMore();
		break;
	default:
		data = new PTiedEgoCentricNetworkData(nodeJoiningTimes,
				vectorOfEdgeEvents, observationTimeStart, observationTimeEnd);
		initializeMore();
		break;
	}

}

void TiedEgoCentricNetworkModel::getMLE(DoubleVector& beta, DoubleMatrix& cov,
		double logLLHTolerance, int maxNRSteps, double scoreTolerance) {

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
		double relLogLLHImprovement = fabs(
				(currentLogLLH - newLogLLH) / currentLogLLH);
		currentLogLLH = newLogLLH;
		cout << "LogLLH: " << currentLogLLH << endl;

		if (fabs(maxError) < scoreTolerance)
			break;
		if (relLogLLHImprovement < logLLHTolerance)
			break;
	}
}

void TiedEgoCentricNetworkModel::getRecencyBaselineRanks(
		ofstream& byEdgesOutputFile, ofstream& byEventsOutputFile,
		ofstream& accumlativeOutputFile, const std::string& tieMethod) {
	// start rolling edges
	data->start();

	double accumRanks = 0;
	for (unsigned int e = data->beginningEdgeIndex, index = 0;
			e <= data->endingEdgeIndex; e++, index++) {
		if (index % 100 == 0)
			cout << "Working up to " << index << endl;

		unsigned int citingNode = vectorOfEdgeEvents[e].first;
		double currentCitingTime = nodeJoiningTimes[citingNode];

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
		unsigned int numOfCitedNodes = vectorOfEdgeEvents[e].second.size();
		for (unsigned int i = 0; i < numOfCitedNodes; i++) {
			unsigned int citedNode = vectorOfEdgeEvents[e].second[i];
			byEdgesOutputFile << citingNode << "\t" << citedNode << "\t"
					<< currentCitingTime << "\t" << currentNodeRanks[citedNode]
					<< endl;
			eventAccumRanks += currentNodeRanks[citedNode];
		}
		byEventsOutputFile << currentCitingTime << "\t"
				<< eventAccumRanks * 1.0 / numOfCitedNodes << endl;

		accumRanks += eventAccumRanks;
		accumlativeOutputFile << currentCitingTime << "\t" << accumRanks
				<< endl;

		// Update statistics and add the edge to the current graph
		data->updateNodalNetworkStatistics(vectorOfEdgeEvents[e]);
	}

	// finish rolling edges
	data->finish();
}

void TiedEgoCentricNetworkModel::prepareExposureIndicators(
		bool* exposureIndicators) {
	for (unsigned int i = 0; i < numOfVertices; i++) {
		// a node is exposed if it joined before the observation time
		// not include this time
		exposureIndicators[i] = (nodeJoiningTimes[i] < observationTimeStart);
	}
}

void TiedEgoCentricNetworkModel::recordNetworkStatistics(unsigned int nodeID,
		ofstream& outputFile) {

	// start rolling edges
	data->start();

	for (unsigned int e = data->beginningEdgeIndex, index = 0;
			e <= data->endingEdgeIndex; e++, index++) {
		unsigned int citingNode = vectorOfEdgeEvents[e].first;
		double currentCitingTime = nodeJoiningTimes[citingNode];

		if (currentCitingTime > nodeJoiningTimes[nodeID]) {
			outputFile << currentCitingTime << "\t";
			for (unsigned int k = 0; k < data->numOfNodalNetworkStatistics; k++)
				outputFile << data->currentNodalNetworkStatistics(nodeID, k)
						<< "\t";
			outputFile << endl;
		}

		data->updateNodalNetworkStatistics(e);
	}

	// finish rolling edges
	data->finish();
}

}
