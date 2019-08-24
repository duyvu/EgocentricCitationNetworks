/*
 * LDATiedEgoCentricNetworkData.cpp
 *
 *  Created on: Jan 21, 2011
 *      Author: duyvu
 */

#include "LDATiedEgoCentricNetworkData.h"

#include "util/StringTokenizer.h"

#include <iostream>
#include <fstream>

namespace ndip {

LDATiedEgoCentricNetworkData::LDATiedEgoCentricNetworkData() {
	// TODO Auto-generated constructor stub

}

LDATiedEgoCentricNetworkData::LDATiedEgoCentricNetworkData(
		const ExposureTimeVectorType& _nodeJoiningTimes,
		const VectorOfEdgeEvents& _vectorOfEdgeEvents,
		double _observationTimeStart, double _observationTimeEnd,
		const std::vector<string>& listOfNodalDataFiles) :
	TiedEgoCentricNetworkData(_nodeJoiningTimes, _vectorOfEdgeEvents,
			_observationTimeStart, _observationTimeEnd) {

	cout << "START the constructor of LDATiedEgoCentricNetworkData" << endl;

	// Initialize the current nodal network statistics
	numOfNodalNetworkStatistics
			= LDATiedEgoCentricNetworkData::NUMBER_OF_NODAL_NETWORK_STATISTICS;
	currentNodalNetworkStatistics = DoubleMatrix(numOfVertices,
			numOfNodalNetworkStatistics);
	currentNodalNetworkStatistics.clear();

	// Initialize the current graph
	currentGraph = Graph(numOfVertices);

	// Initialize the recent node contacted times
	currentRecentNodeCitedTimes = ExposureTimeVectorType(numOfVertices);
	for (unsigned int i = 0; i < numOfVertices; i++)
		currentRecentNodeCitedTimes[i] = nodeJoiningTimes[i];

	// Read LDA matrix
	const char* nodalDataFile = listOfNodalDataFiles[0].c_str();
	cout << "LDA file is " << nodalDataFile << endl;
	string line;
	ifstream nodalFile(nodalDataFile);
	LDA = DoubleMatrix(numOfVertices, NUMBER_OF_LDA_TOPICS);
	LDA.clear();
	if (nodalFile.is_open()) {
		getline(nodalFile, line); // ignore number of topics
		while (!nodalFile.eof()) {
			getline(nodalFile, line);
			StringTokenizer strtok = StringTokenizer(line, "\t");
			if (strtok.countTokens() == 1 + NUMBER_OF_LDA_TOPICS) {
				int nodeID = strtok.nextIntToken();
				//cout << nodeID << ": ";
				for (unsigned int k = 0; k < NUMBER_OF_LDA_TOPICS; k++) {
					double ldaProp = strtok.nextDoubleToken();
					LDA(nodeID, k) = ldaProp;
					//cout << ldaProp << "\t";
				}
				//cout << endl;
			} else {
				int nodeID = strtok.nextIntToken();
				//cout << nodeID << ": ";
				cout << "Some LDA covariates " << nodeID
						<< " are missing!!! I will assume them zero" << endl;
				unsigned int k = 0;
				while (strtok.hasMoreTokens()) {
					double ldaProp = strtok.nextDoubleToken();
					LDA(nodeID, k++) = ldaProp;
					//cout << ldaProp << "\t";
				}
				//cout << endl;
			}
		}
		nodalFile.close();
	} else
		cout << "Unable to open file " << nodalDataFile << endl;

	// Roll the current nodal network statistics
	// and the graph up to the observation time
	for (unsigned int indexOfEdgeEvents = 0; indexOfEdgeEvents
			< beginningEdgeIndex; indexOfEdgeEvents++) {
		unsigned int citingNode = vectorOfEdgeEvents[indexOfEdgeEvents].first;
		double citingTime = nodeJoiningTimes[citingNode];
		for (unsigned int i = 0; i
				< vectorOfEdgeEvents[indexOfEdgeEvents].second.size(); i++) {
			unsigned int citedNode =
					vectorOfEdgeEvents[indexOfEdgeEvents].second[i];
			currentRecentNodeCitedTimes[citedNode] = citingTime;
		}
	}
	// update the LDA covariates
	unsigned int nextCitingNode = vectorOfEdgeEvents[beginningEdgeIndex].first;
	double nextCitingTime = nodeJoiningTimes[nextCitingNode];
	for (unsigned int i = 0; i < numOfVertices; i++) {
		if (nodeJoiningTimes[i] < nextCitingTime) {
			for (unsigned int ldaIndex = 0; ldaIndex < NUMBER_OF_LDA_TOPICS; ldaIndex++)
				currentNodalNetworkStatistics(i, ldaIndex) = LDA(i, ldaIndex)
						* LDA(nextCitingNode, ldaIndex);
		} else
			break;
	}

	// Copy these current nodal network statistics and the graph up to the observation time
	// to save the time for iterative calls later
	beginningNodalNetworkStatistics = DoubleMatrix(
			currentNodalNetworkStatistics);
	beginningGraph = Graph(currentGraph);
	beginningRecentNodeCitedTimes = ExposureTimeVectorType(
			currentRecentNodeCitedTimes);

	cout << "No cache is in use" << endl;

	// reset everything for the next use
	dirty = true;
	finish();

	cout << "FINISH the constructor of LDATiedEgoCentricNetworkData" << endl;
}

LDATiedEgoCentricNetworkData::~LDATiedEgoCentricNetworkData() {
	// TODO Auto-generated destructor stub
}

void LDATiedEgoCentricNetworkData::updateNodalNetworkStatistics(
		unsigned int indexOfEdgeEvents) {

	unsigned int citingNode = vectorOfEdgeEvents[indexOfEdgeEvents].first;
	double citingTime = nodeJoiningTimes[citingNode];
	for (unsigned int i = 0; i
			< vectorOfEdgeEvents[indexOfEdgeEvents].second.size(); i++) {
		unsigned int citedNode =
				vectorOfEdgeEvents[indexOfEdgeEvents].second[i];
		currentRecentNodeCitedTimes[citedNode] = citingTime;
	}

	// update the LDA covariates
	unsigned int nextCitingNode =
			vectorOfEdgeEvents[indexOfEdgeEvents + 1].first;
	double nextCitingTime = nodeJoiningTimes[nextCitingNode];
	for (unsigned int i = 0; i < numOfVertices; i++) {
		if (nodeJoiningTimes[i] < nextCitingTime) {
			for (unsigned int ldaIndex = 0; ldaIndex < NUMBER_OF_LDA_TOPICS; ldaIndex++)
				currentNodalNetworkStatistics(i, ldaIndex) = LDA(i, ldaIndex)
						* LDA(nextCitingNode, ldaIndex);
		} else
			break;
	}

	// for the next VectorOfUpdateElements
	currentCacheIndex++;
}

int LDATiedEgoCentricNetworkData::getFirstOrderPopularityIndex() {
	return -1;
}

}
