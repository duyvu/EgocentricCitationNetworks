/*
 * icmlArXivHepTh.h
 *
 *  Created on: Jan 20, 2011
 *      Author: duyvu
 */

#ifndef ICMLARXIVHEPTH_H_
#define ICMLARXIVHEPTH_H_

#include "util/StringTokenizer.h"

#include "nodal/data/TiedEgoCentricNetworkData.h"
#include "nodal/data/PTiedEgoCentricNetworkData.h"
#include "nodal/data/PGTiedEgoCentricNetworkData.h"
#include "nodal/data/P2PTTiedEgoCentricNetworkData.h"
#include "nodal/data/P2PTRTiedEgoCentricNetworkData.h"
#include "nodal/data/LDAP2PTRTiedEgoCentricNetworkData.h"
#include "nodal/data/LDATiedEgoCentricNetworkData.h"

#include "nodal/model/cox/TiedEgoCentricNetworkModel.h"
#include "nodal/model/cox/NTCTiedEgoCentricNetworkModel.h"
#include "nodal/model/cox/DenseNTCTiedEgoCentricNetworkModel.h"

#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <ctime>

using namespace std;
using namespace ndip;

void trainTiedEgoCentricNetworkModel(unsigned int modelType,
		double trainObservationTimeStart, double trainObservationTimeEnd,
		const char * signature) {

	cout << "TRAINING TiedEgoCentricNetworkModel" << endl;
	cout.precision(16);

	TiedEgoCentricNetworkModel* traingModel;
	if (modelType == 3 || modelType == 4 || modelType == 5)
		traingModel = new DenseNTCTiedEgoCentricNetworkModel();
	else
		traingModel = new NTCTiedEgoCentricNetworkModel();

	char* nodeJoiningTimesFileName = "./data/CitationHepTh/Cit-HepTh-nodes.txt";
	char* edgeEventsFileName = "./data/CitationHepTh/Cit-HepTh-edges.txt";
	traingModel->setNetworkData(nodeJoiningTimesFileName, edgeEventsFileName);
	cout << "Done with model->setNetworkData()" << endl;

	std::vector<string> trainNodalFiles;
	int numOfNodalNetworkStatistics = 0;
	cout << "Train Observation Time Start is " << trainObservationTimeStart
			<< endl;
	cout << "Train Observation Time End is " << trainObservationTimeEnd << endl;
	ofstream betaOutputFile;
	betaOutputFile.precision(16);
	std::stringstream postfix;
	postfix << "./output/";
	switch (modelType) {
	case 0:
		cout << "PA model in use" << endl;
		traingModel->setModel(PREFERENTIAL_ATTACHMENT_EGO_CENTRIC_MODEL,
				trainObservationTimeStart, trainObservationTimeEnd,
				trainNodalFiles);
		numOfNodalNetworkStatistics
				= PTiedEgoCentricNetworkData::NUMBER_OF_NODAL_NETWORK_STATISTICS;
		postfix << "PA-Beta-" << signature << ".txt";
		betaOutputFile.open(postfix.str().c_str());
		betaOutputFile << "PA " << trainObservationTimeStart << "\t"
				<< trainObservationTimeEnd << endl;
		break;
	case 1:
		cout << "P2PT model in use" << endl;
		traingModel->setModel(P2PT_EGO_CENTRIC_MODEL,
				trainObservationTimeStart, trainObservationTimeEnd,
				trainNodalFiles);
		numOfNodalNetworkStatistics
				= P2PTTiedEgoCentricNetworkData::NUMBER_OF_NODAL_NETWORK_STATISTICS;
		postfix << "P2PT-Beta-" << signature << ".txt";
		betaOutputFile.open(postfix.str().c_str());
		betaOutputFile << "P2PT " << trainObservationTimeStart << "\t"
				<< trainObservationTimeEnd << endl;
		break;
	case 2:
		cout << "P2PTR model in use" << endl;
		traingModel->setModel(P2PTR_EGO_CENTRIC_MODEL,
				trainObservationTimeStart, trainObservationTimeEnd,
				trainNodalFiles);
		numOfNodalNetworkStatistics
				= P2PTRTiedEgoCentricNetworkData::NUMBER_OF_NODAL_NETWORK_STATISTICS;
		postfix << "P2PTR-Beta-" << signature << ".txt";
		betaOutputFile.open(postfix.str().c_str());
		betaOutputFile << "P2PTR " << trainObservationTimeStart << "\t"
				<< trainObservationTimeEnd << endl;
		break;
	case 3:
		cout << "LDA_P2PTR model in use" << endl;
		trainNodalFiles = std::vector<string>(1);
		trainNodalFiles[0] = "./data/CitationHepTh/Cit-HepTh-LDA-50-mapped.txt";
		traingModel->setModel(LDA_P2PTR_EGO_CENTRIC_MODEL,
				trainObservationTimeStart, trainObservationTimeEnd,
				trainNodalFiles);
		numOfNodalNetworkStatistics
				= LDAP2PTRTiedEgoCentricNetworkData::NUMBER_OF_NODAL_NETWORK_STATISTICS;
		postfix << "LDA-P2PTR-Beta-" << signature << ".txt";
		betaOutputFile.open(postfix.str().c_str());
		betaOutputFile << "LDA_P2PTR " << trainObservationTimeStart << "\t"
				<< trainObservationTimeEnd << endl;
		break;
	case 4:
		cout << "LDA model in use" << endl;
		trainNodalFiles = std::vector<string>(1);
		trainNodalFiles[0] = "./data/CitationHepTh/Cit-HepTh-LDA-50-mapped.txt";
		traingModel->setModel(LDA_EGO_CENTRIC_MODEL, trainObservationTimeStart,
				trainObservationTimeEnd, trainNodalFiles);
		numOfNodalNetworkStatistics
				= LDATiedEgoCentricNetworkData::NUMBER_OF_NODAL_NETWORK_STATISTICS;
		postfix << "LDA-Beta-" << signature << ".txt";
		betaOutputFile.open(postfix.str().c_str());
		betaOutputFile << "LDA " << trainObservationTimeStart << "\t"
				<< trainObservationTimeEnd << endl;
		break;
	case 5:
		cout << "PA and GAP Time model in use" << endl;
		traingModel->setModel(PA_GAP_EGO_CENTRIC_MODEL,
				trainObservationTimeStart, trainObservationTimeEnd,
				trainNodalFiles);
		numOfNodalNetworkStatistics
				= PGTiedEgoCentricNetworkData::NUMBER_OF_NODAL_NETWORK_STATISTICS;
		postfix << "PA-GAP-Beta-" << signature << ".txt";
		betaOutputFile.open(postfix.str().c_str());
		betaOutputFile << "PA-GAP " << trainObservationTimeStart << "\t"
				<< trainObservationTimeEnd << endl;
		break;
	default:
		cout << "PA model in use" << endl;
		traingModel->setModel(PREFERENTIAL_ATTACHMENT_EGO_CENTRIC_MODEL,
				trainObservationTimeStart, trainObservationTimeEnd,
				trainNodalFiles);
		numOfNodalNetworkStatistics
				= PTiedEgoCentricNetworkData::NUMBER_OF_NODAL_NETWORK_STATISTICS;
		postfix << "PA-Beta-" << signature << ".txt";
		betaOutputFile.open(postfix.str().c_str());
		betaOutputFile << "PA " << trainObservationTimeStart << "\t"
				<< trainObservationTimeEnd << endl;
		break;
	}
	cout << "Done with model->setModel()" << endl;

	clock_t start;

	DoubleVector beta(numOfNodalNetworkStatistics);
	beta.clear();

	double tolLogLLH = 1e-100;

	DoubleMatrix cov(beta.size(), beta.size());
	cov.clear();

	cout << "Running MLE" << endl;

	start = std::clock();
	traingModel->getMLE(beta, cov, tolLogLLH, 50, 1e-3);

	cout << "Done with MLE" << endl;

	cout << "----------------FINAL RESULT----------------" << endl;

	cout << "Time of Newton-Raphson with LBF: " << (clock() - start) / 1e6
			<< " seconds" << endl;

	cout << "Computing training log-likelihood" << endl;
	start = std::clock();
	double trainingLogLLH = traingModel->computeLogLikelihood(beta, true);
	cout << "Time of training log-likelihood computation with LBF: "
			<< (clock() - start) / 1e6 << " seconds" << endl;
	cout << "training logLLH = " << trainingLogLLH << endl;
	betaOutputFile << trainingLogLLH << endl;

	cout << "Beta Estimates: " << endl;
	for (int h = 0; h < beta.size(); h++) {
		cout << beta(h) << " ";
		betaOutputFile << beta(h) << "\t";
	}
	cout << endl;
	betaOutputFile << endl;

	cout << "Beta Covariance: " << endl;
	for (int k = 0; k < cov.size1(); k++) {
		for (int h = 0; h < cov.size2(); h++) {
			cout << cov(k, h) << "  ";
			betaOutputFile << cov(k, h) << "\t";
		}
		cout << endl;
		betaOutputFile << endl;
	}

	betaOutputFile.close();
	if (traingModel != NULL)
		delete traingModel;
}

void testTiedEgoCentricNetworkModel(unsigned int modelType,
		double testObservationTimeStart, double testObservationTimeEnd,
		const char * signature) {

	cout << "TESTING TiedEgoCentricNetworkModel" << endl;
	cout.precision(16);

	TiedEgoCentricNetworkModel* testModel;
	if (modelType == 3 || modelType == 4 || modelType == 5)
		testModel = new DenseNTCTiedEgoCentricNetworkModel();
	else
		testModel = new NTCTiedEgoCentricNetworkModel();

	char* nodeJoiningTimesFileName = "./data/CitationHepTh/Cit-HepTh-nodes.txt";
	char* edgeEventsFileName = "./data/CitationHepTh/Cit-HepTh-edges.txt";
	testModel->setNetworkData(nodeJoiningTimesFileName, edgeEventsFileName);
	cout << "Done with model->setNetworkData()" << endl;

	std::vector<string> testNodalFiles;

	int numOfNodalNetworkStatistics = 0;
	std::stringstream betaPostfix;
	betaPostfix << "./data/CitationHepTh/icml/";
	DoubleVector beta;
	ifstream betaInputFile;
	string line;
	StringTokenizer tokenizer = StringTokenizer("", "");
	unsigned int betaIndex = 0;

	cout << "Test Observation Time Start is " << testObservationTimeStart
			<< endl;
	cout << "Test Observation Time End is " << testObservationTimeEnd << endl;

	ofstream testRankByEdgesOutputFile;
	testRankByEdgesOutputFile.precision(16);
	std::stringstream testRankByEdgesPostfix;
	testRankByEdgesPostfix << "./data/CitationHepTh/icml/";

	ofstream testRankByEventsOutputFile;
	testRankByEventsOutputFile.precision(16);
	std::stringstream testRankByEventsPostfix;
	testRankByEventsPostfix << "./data/CitationHepTh/icml/";

	ofstream testRankAccumOutputFile;
	testRankAccumOutputFile.precision(16);
	std::stringstream testRankAccumPostfix;
	testRankAccumPostfix << "./data/CitationHepTh/icml/";

	ofstream testLikelihoodByEdgesOutputFile;
	testLikelihoodByEdgesOutputFile.precision(16);
	std::stringstream testLikelihoodByEdgesPostfix;
	testLikelihoodByEdgesPostfix << "./data/CitationHepTh/icml/";

	switch (modelType) {
	case 0:
		cout << "PA model in use" << endl;
		testModel->setModel(PREFERENTIAL_ATTACHMENT_EGO_CENTRIC_MODEL,
				testObservationTimeStart, testObservationTimeEnd,
				testNodalFiles);
		numOfNodalNetworkStatistics
				= PTiedEgoCentricNetworkData::NUMBER_OF_NODAL_NETWORK_STATISTICS;

		betaPostfix << "PA-Beta-" << signature << ".txt";

		testRankByEdgesPostfix << "PA-test-rank-by-edges-" << signature
				<< ".txt";
		testRankByEventsPostfix << "PA-test-rank-by-events-" << signature
				<< ".txt";
		testRankAccumPostfix << "PA-test-rank-accum-" << signature << ".txt";

		testLikelihoodByEdgesPostfix << "PA-test-likelihood-by-edges-"
				<< signature << ".txt";
		break;
	case 1:
		cout << "P2PT model in use" << endl;
		testModel->setModel(P2PT_EGO_CENTRIC_MODEL, testObservationTimeStart,
				testObservationTimeEnd, testNodalFiles);
		numOfNodalNetworkStatistics
				= P2PTTiedEgoCentricNetworkData::NUMBER_OF_NODAL_NETWORK_STATISTICS;

		betaPostfix << "P2PT-Beta-" << signature << ".txt";

		testRankByEdgesPostfix << "P2PT-test-rank-by-edges-" << signature
				<< ".txt";
		testRankByEventsPostfix << "P2PT-test-rank-by-events-" << signature
				<< ".txt";
		testRankAccumPostfix << "P2PT-test-rank-accum-" << signature << ".txt";
		testLikelihoodByEdgesPostfix << "P2PT-test-likelihood-by-edges-"
				<< signature << ".txt";
		break;
	case 2:
		cout << "P2PTR model in use" << endl;
		testModel->setModel(P2PTR_EGO_CENTRIC_MODEL, testObservationTimeStart,
				testObservationTimeEnd, testNodalFiles);
		numOfNodalNetworkStatistics
				= P2PTRTiedEgoCentricNetworkData::NUMBER_OF_NODAL_NETWORK_STATISTICS;

		betaPostfix << "P2PTR-Beta-" << signature << ".txt";

		testRankByEdgesPostfix << "P2PTR-test-rank-by-edges-" << signature
				<< ".txt";
		testRankByEventsPostfix << "P2PTR-test-rank-by-events-" << signature
				<< ".txt";
		testRankAccumPostfix << "P2PTR-test-rank-accum-" << signature << ".txt";
		testLikelihoodByEdgesPostfix << "P2PTR-test-likelihood-by-edges-"
				<< signature << ".txt";
		break;
	case 3:
		cout << "LDA_P2PTR model in use" << endl;
		testNodalFiles = std::vector<string>(1);
		testNodalFiles[0] = "./data/CitationHepTh/Cit-HepTh-LDA-50-mapped.txt";
		testModel->setModel(LDA_P2PTR_EGO_CENTRIC_MODEL,
				testObservationTimeStart, testObservationTimeEnd,
				testNodalFiles);
		numOfNodalNetworkStatistics
				= LDAP2PTRTiedEgoCentricNetworkData::NUMBER_OF_NODAL_NETWORK_STATISTICS;

		betaPostfix << "LDA-P2PTR-Beta-" << signature << ".txt";

		testRankByEdgesPostfix << "LDA-P2PTR-test-rank-by-edges-" << signature
				<< ".txt";
		testRankByEventsPostfix << "LDA-P2PTR-test-rank-by-events-"
				<< signature << ".txt";
		testRankAccumPostfix << "LDA-P2PTR-test-rank-accum-" << signature
				<< ".txt";
		testLikelihoodByEdgesPostfix << "LDA-P2PTR-test-likelihood-by-edges-"
				<< signature << ".txt";
		break;
	case 4:
		cout << "LDA model in use" << endl;
		testNodalFiles = std::vector<string>(1);
		testNodalFiles[0] = "./data/CitationHepTh/Cit-HepTh-LDA-50-mapped.txt";
		testModel->setModel(LDA_EGO_CENTRIC_MODEL, testObservationTimeStart,
				testObservationTimeEnd, testNodalFiles);
		numOfNodalNetworkStatistics
				= LDATiedEgoCentricNetworkData::NUMBER_OF_NODAL_NETWORK_STATISTICS;

		betaPostfix << "LDA-Beta-" << signature << ".txt";

		testRankByEdgesPostfix << "LDA-test-rank-by-edges-" << signature
				<< ".txt";
		testRankByEventsPostfix << "LDA-test-rank-by-events-" << signature
				<< ".txt";
		testRankAccumPostfix << "LDA-test-rank-accum-" << signature << ".txt";
		testLikelihoodByEdgesPostfix << "LDA-test-likelihood-by-edges-"
				<< signature << ".txt";
		break;
	case 5:
		cout << "PA GAP model in use" << endl;
		testModel->setModel(PA_GAP_EGO_CENTRIC_MODEL, testObservationTimeStart,
				testObservationTimeEnd, testNodalFiles);
		numOfNodalNetworkStatistics
				= PGTiedEgoCentricNetworkData::NUMBER_OF_NODAL_NETWORK_STATISTICS;

		betaPostfix << "PA-GAP-Beta-" << signature << ".txt";

		testRankByEdgesPostfix << "PA-GAP-test-rank-by-edges-" << signature
				<< ".txt";
		testRankByEventsPostfix << "PA-GAP-test-rank-by-events-" << signature
				<< ".txt";
		testRankAccumPostfix << "PA-GAP-test-rank-accum-" << signature
				<< ".txt";
		testLikelihoodByEdgesPostfix << "PA-GAP-test-likelihood-by-edges-"
				<< signature << ".txt";
		break;
	default:
		cout << "Please enter a correct model type!" << endl;
	}
	cout << "Done with model->setModel()" << endl;

	beta = DoubleVector(numOfNodalNetworkStatistics);
	cout << "Beta file is " << betaPostfix.str() << endl;
	betaInputFile.open(betaPostfix.str().c_str());
	getline(betaInputFile, line);
	getline(betaInputFile, line);
	getline(betaInputFile, line);
	tokenizer = StringTokenizer(line, "\t");
	cout << "Beta: ";
	while (tokenizer.hasMoreTokens()) {
		beta(betaIndex++) = tokenizer.nextDoubleToken();
		cout << beta(betaIndex - 1) << "\t";
	}
	cout << endl;
	betaInputFile.close();

	testRankByEdgesOutputFile.open(testRankByEdgesPostfix.str().c_str());
	testRankByEdgesOutputFile << testObservationTimeStart << "\t"
			<< testObservationTimeEnd << endl;

	testRankByEventsOutputFile.open(testRankByEventsPostfix.str().c_str());
	testRankByEventsOutputFile << testObservationTimeStart << "\t"
			<< testObservationTimeEnd << endl;

	testRankAccumOutputFile.open(testRankAccumPostfix.str().c_str());
	testRankAccumOutputFile << testObservationTimeStart << "\t"
			<< testObservationTimeEnd << endl;

	string tieMethod = "average";
	testModel->getRanks(testRankByEdgesOutputFile, testRankByEventsOutputFile,
			testRankAccumOutputFile, beta, tieMethod);

	testRankByEdgesOutputFile.close();
	testRankByEventsOutputFile.close();
	testRankAccumOutputFile.close();

	clock_t start;

	//	cout << "Computing test log-likelihood" << endl;
	//	start = std::clock();
	//	double testLogLLH = testModel->computeLogLikelihood(beta, true);
	//	cout << "Time of test log-likelihood computation with LBF: " << (clock()
	//			- start) / 1e6 << " seconds" << endl;
	//	cout << "Test logLLH = " << testLogLLH << endl;

	testLikelihoodByEdgesOutputFile.open(
			testLikelihoodByEdgesPostfix.str().c_str());
	testLikelihoodByEdgesOutputFile << testObservationTimeStart << "\t"
			<< testObservationTimeEnd << endl;

	testModel->computeLogLikelihoodOverEdges(testLikelihoodByEdgesOutputFile,
			beta, true);

	testLikelihoodByEdgesOutputFile.close();

	ofstream testRankByEdgesBaselineOutputFile;
	testRankByEdgesBaselineOutputFile.open(
			"./data/CitationHepTh/icml/baseline-test-rank-by-edges.txt");
	testRankByEdgesBaselineOutputFile.precision(16);
	ofstream testRankByEventsBaselineOutputFile;
	testRankByEventsBaselineOutputFile.open(
			"./data/CitationHepTh/icml/baseline-test-rank-by-events.txt");
	testRankByEventsBaselineOutputFile.precision(16);
	ofstream testRankAccumBaselineOutputFile;
	testRankAccumBaselineOutputFile.open(
			"./data/CitationHepTh/icml/baseline-test-rank-accum.txt");
	testRankAccumBaselineOutputFile.precision(16);

	testModel->getRecencyBaselineRanks(testRankByEdgesBaselineOutputFile,
			testRankByEventsBaselineOutputFile,
			testRankAccumBaselineOutputFile, tieMethod);

	testRankByEdgesBaselineOutputFile.close();
	testRankByEventsBaselineOutputFile.close();
	testRankAccumBaselineOutputFile.close();

	if (testModel != NULL)
		delete testModel;
}

void getPiecewiseCummulativeBaselines(unsigned int modelType,
		double observationTimeStart, double observationTimeEnd,
		const char * signature) {
	cout << "Baseline TiedEgoCentricNetworkModel" << endl;
	cout.precision(16);

	TiedEgoCentricNetworkModel* testModel;
	if (modelType == 3 || modelType == 4 || modelType == 5)
		testModel = new DenseNTCTiedEgoCentricNetworkModel();
	else
		testModel = new NTCTiedEgoCentricNetworkModel();

	char* nodeJoiningTimesFileName = "./data/CitationHepTh/Cit-HepTh-nodes.txt";
	char* edgeEventsFileName = "./data/CitationHepTh/Cit-HepTh-edges.txt";
	testModel->setNetworkData(nodeJoiningTimesFileName, edgeEventsFileName);
	cout << "Done with model->setNetworkData()" << endl;

	std::vector<string> testNodalFiles;

	int numOfNodalNetworkStatistics = 0;
	std::stringstream betaPostfix;
	betaPostfix << "./data/CitationHepTh/icml/";
	DoubleVector beta;
	ifstream betaInputFile;
	string line;
	StringTokenizer tokenizer = StringTokenizer("", "");
	unsigned int betaIndex = 0;

	cout << "Test Observation Time Start is " << observationTimeStart << endl;
	cout << "Test Observation Time End is " << observationTimeEnd << endl;

	cout << "LDA_P2PTR model in use" << endl;
	testNodalFiles = std::vector<string>(1);
	testNodalFiles[0] = "./data/CitationHepTh/Cit-HepTh-LDA-50-mapped.txt";
	testModel->setModel(LDA_P2PTR_EGO_CENTRIC_MODEL, observationTimeStart,
			observationTimeEnd, testNodalFiles);
	numOfNodalNetworkStatistics
			= LDAP2PTRTiedEgoCentricNetworkData::NUMBER_OF_NODAL_NETWORK_STATISTICS;

	betaPostfix << "LDA-P2PTR-Beta-" << signature << ".txt";

	cout << "Done with model->setModel()" << endl;

	beta = DoubleVector(numOfNodalNetworkStatistics);
	cout << "Beta file is " << betaPostfix.str() << endl;
	betaInputFile.open(betaPostfix.str().c_str());
	getline(betaInputFile, line);
	getline(betaInputFile, line);
	getline(betaInputFile, line);
	tokenizer = StringTokenizer(line, "\t");
	cout << "Beta: ";
	while (tokenizer.hasMoreTokens()) {
		beta(betaIndex++) = tokenizer.nextDoubleToken();
		cout << beta(betaIndex - 1) << "\t";
	}
	cout << endl;
	betaInputFile.close();

	VectorOfCummulativePoints cumBH = VectorOfCummulativePoints();
	testModel->computeCummulativeBaselineHazard(cumBH, beta);
	ofstream outputFile;
	outputFile.open("./data/CitationHepTh/icml/cumBH.txt");
	outputFile.precision(20);
	for (VectorOfCummulativePoints::iterator it = cumBH.begin(); it
			!= cumBH.end(); it++) {
		CummulativePoint point = *it;
		outputFile << point.first << "\t" << point.second << endl;
	}
	outputFile.close();

}

void mineSomePatterns(unsigned int nodeID, double observationTimeStart,
		double observationTimeEnd) {
	cout << "MINING TiedEgoCentricNetworkModel" << endl;
	cout.precision(16);

	TiedEgoCentricNetworkModel* testModel;
	//testModel = new NTCTiedEgoCentricNetworkModel();
	testModel = new DenseNTCTiedEgoCentricNetworkModel();

	char* nodeJoiningTimesFileName = "./data/CitationHepTh/Cit-HepTh-nodes.txt";
	char* edgeEventsFileName = "./data/CitationHepTh/Cit-HepTh-edges.txt";
	testModel->setNetworkData(nodeJoiningTimesFileName, edgeEventsFileName);

	//std::vector<string> nodalFiles;
	//testModel->setModel(P2PTR_EGO_CENTRIC_MODEL, observationTimeStart,
	//		observationTimeEnd, nodalFiles);

	std::vector<string> nodalFiles = std::vector<string>(1);
	nodalFiles[0] = "./data/CitationHepTh/Cit-HepTh-LDA-50-mapped.txt";
	testModel->setModel(LDA_P2PTR_EGO_CENTRIC_MODEL, observationTimeStart,
			observationTimeEnd, nodalFiles);

	std::stringstream outputPostfix;
	outputPostfix << "./data/CitationHepTh/icml/NodeStat_";
	outputPostfix << nodeID << ".txt";
	ofstream outputFile;
	outputFile.open(outputPostfix.str().c_str());
	outputFile.precision(16);
	testModel->recordNetworkStatistics(nodeID, outputFile);
	outputFile.close();

	if (testModel != NULL)
		delete testModel;
}

#endif /* ICMLARXIVHEPTH_H_ */

