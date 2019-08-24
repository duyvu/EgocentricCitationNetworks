/*
 * TiedEgoCentricNetworkModel.h
 *
 *  Created on: Jul 18, 2010
 *      Author: duyvu
 */

#ifndef TIEDEGOCENTRICNETWORKMODEL_H_
#define TIEDEGOCENTRICNETWORKMODEL_H_

#include "nodal/data/TiedEgoCentricNetworkData.h"

namespace ndip {

class TiedEgoCentricNetworkModel {
protected:
	// the number of nodes
	unsigned int numOfVertices;
	// a time-ordered list of nodes and corresponding join times. It is used to compute the age time.
	ExposureTimeVectorType nodeJoiningTimes;

	// the number of distinct edge events
	unsigned int numOfDistinctEdgeEvents;
	// a time-ordered list of edges and corresponding join times
	VectorOfEdgeEvents vectorOfEdgeEvents;

	// assume the observation time is the time of the first edge
	// recorded on the continuous-time scale. Otherwise, the first element of vectorOfComingVertices
	// is the list of vertices coming from the observation time to the first edge after this time.
	// The assumption guarantees the synchronization of the new node lists and the update lists.
	double observationTimeStart;

	double observationTimeEnd;

	TiedEgoCentricNetworkData* data;

	void prepareExposureIndicators(bool* exposureIndicators);

	// for late initialization, i.e. it is used to initialize those
	// member variables which are only present in subclasses
	virtual void initializeMore() = 0;
public:
	TiedEgoCentricNetworkModel();
	virtual ~TiedEgoCentricNetworkModel();

	// read data of nodes and edges from files
	void setNetworkData(char* nodeJoiningTimesFileName,
			char* edgeEventsFileName);

	// choose which ego centric model to use
	// for a set of network statistics of interest,
	// we will build an corresponding ego centric network data
	// the likelhood and ML estimation code will be independent with
	// the selected ego centric network data
	void setModel(int _modelType, double _observationTimeStart,
			double _observationTimeEnd,
			const std::vector<string>& listOfNodalDataFiles);

	// timing covariates such as age or gap will be handled in
	// different subclasses of EgoCentricNetworkModel
	// Each subclass will compute this likelihood in different ways
	// depending on age or gap time is used
	// this method will be used for free-gradient optimization methods
	// such as Nelder-Mead
	virtual double computeLogLikelihood(const DoubleVector& beta,
			bool isUsingLBF) = 0;

	virtual void computeLogLikelihoodOverEvents(
			DoubleVector& cummulativeLogLLH, const DoubleVector& beta,
			bool isUsingLBF) = 0;

	virtual void computeLogLikelihoodOverEdges(ofstream& byEdgesOutputFile,
			const DoubleVector& beta, bool isUsingLBF) = 0;

	// this method will be used for Newton-Raphson optimization methods
	virtual void computeScoreVectorAndInformationMatrix(DoubleVector& u,
			DoubleMatrix& I, const DoubleVector& beta) = 0;

	// this method will be used for gradient optimization methods
	virtual double computeLogLikelihoodAndScoreVector(DoubleVector& u,
			const DoubleVector& beta) = 0;

	virtual void computeOnlineGradientAscentMLE(DoubleVector& beta,
			double learningRate) = 0;

	// get ML estimates of parameters
	void getMLE(DoubleVector& beta, DoubleMatrix& cov, double logLLHTolerance,
			int maxNRSteps, double scoreTolerance);

	// a helper function for Nelder-Mead optimization method
	// in case we only use function values to get ML estimates
	static double computeMinusLogLLH(int n, const double *x, double *grad,
			void *my_func_data) {

		TiedEgoCentricNetworkModel* model =
				(TiedEgoCentricNetworkModel*) my_func_data;

		DoubleVector beta(n);
		for (int k = 0; k < n; k++)
			beta(k) = x[k];

		if (grad) {
			DoubleVector u(n);
			double logLLH = model->computeLogLikelihoodAndScoreVector(u, beta);
			logLLH *= -1.0;
			for (int k = 0; k < n; k++)
				grad[k] = -u[k];

			cout << "x = ";
			for (int k = 0; k < n; k++)
				cout << x[k] << ", ";
			cout << endl;
			cout << "log-LLH = " << logLLH << endl;
			cout << "grad = ";
			for (int k = 0; k < n; k++)
				cout << grad[k] << ", ";
			cout << endl;

			return logLLH;
		} else {
			double logLLH = model->computeLogLikelihood(beta, true);
			logLLH *= -1.0;
			cout << "log-LLH = " << logLLH << endl;
			return logLLH;
		}
	}

	TiedEgoCentricNetworkData* getNetworkDataObject() {
		return data;
	}

	virtual void computeCummulativeBaselineHazard(
			VectorOfCummulativePoints& cumBH, const DoubleVector& beta) = 0;

	virtual void computeSchoenfeldResiduals(DoubleMatrix& schoenfeldResiduals,
			const DoubleVector& beta, bool isPlusBeta) = 0;

	virtual void getSumRanks(DoubleVector& sumRanks, const DoubleVector& beta,
			const std::string& tieMethod) = 0;
	virtual void getRanks(VectorOfDoubleVector& ranks,
			const DoubleVector& beta, const std::string& tieMethod) = 0;
	virtual void getRanks(ofstream& byEdgesOutputFile,
			ofstream& byEventsOutputFile, ofstream& accumlativeOutputFile,
			const DoubleVector& beta, const std::string& tieMethod) = 0;
	virtual void getRecencyBaselineRanks(ofstream& byEdgesOutputFile,
			ofstream& byEventsOutputFile, ofstream& accumlativeOutputFile,
			const std::string& tieMethod);

	virtual void simulateNetwork(Graph& graph, DoubleMatrix& covariates,
			const DoubleVector& beta) = 0;

	unsigned int getNumOfVertices() {
		return numOfVertices;
	}

	void recordNetworkStatistics(unsigned int nodeID, ofstream& outputFile);
};

}

#endif /* TIEDEGOCENTRICNETWORKMODEL_H_ */
