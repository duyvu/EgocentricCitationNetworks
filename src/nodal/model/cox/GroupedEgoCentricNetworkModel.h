/*
 * GroupedEgoCentricNetworkModel.h
 *
 *  Created on: Jan 17, 2011
 *      Author: duyvu
 */

#ifndef GROUPEDEGOCENTRICNETWORKMODEL_H_
#define GROUPEDEGOCENTRICNETWORKMODEL_H_

#include "nodal/data/GroupedEgoCentricNetworkData.h"

namespace ndip {

class GroupedEgoCentricNetworkModel {
protected:
	// the number of nodes
	unsigned int numOfVertices;
	// a time-ordered list of nodes and corresponding join times. It is used to compute the age time.
	ExposureTimeVectorType nodeJoiningTimes;

	// the number of distinct edge events
	unsigned int numOfDistinctEdgeEvents;
	// a time-ordered list of edges and corresponding join times
	VectorOfDiscreteTimeVectorOfEdgeEvents
			vectorOfDiscreteTimeVectorOfEdgeEvents;

	// assume the observation time is the time of the first edge event
	// recorded on the continuous-time scale. Otherwise, the first element of vectorOfComingVertices
	// is the list of vertices coming from the observation time to the first edge after this time.
	// The assumption guarantees the synchronization of the new node lists and the update lists.
	double observationTimeStart;
	// assume the observation time is the time of the last edge event
	double observationTimeEnd;

	GroupedEgoCentricNetworkData* data;

	void prepareExposureIndicators(bool* exposureIndicators);

	// for late initialization, i.e. it is used to initialize those
	// member variables which are only present in subclasses
	virtual void initializeMore() = 0;
public:
	GroupedEgoCentricNetworkModel();
	virtual ~GroupedEgoCentricNetworkModel();

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

	virtual void computeLogLikelihoodOverEdges(ofstream& outputFile,
				const DoubleVector& beta, bool isUsingLBF) = 0;

	// this method will be used for Newton-Raphson optimization methods
	virtual void computeScoreVectorAndInformationMatrix(DoubleVector& u,
			DoubleMatrix& I, const DoubleVector& beta) = 0;

	// get ML estimates of parameters
	void getMLE(DoubleVector& beta, DoubleMatrix& cov, double logLLHTolerance,
			int maxNRSteps, double scoreTolerance);

	GroupedEgoCentricNetworkData* getNetworkDataObject() {
		return data;
	}

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
	unsigned int getNumOfVertices() {
		return numOfVertices;
	}
};

}

#endif /* GROUPEDEGOCENTRICNETWORKMODEL_H_ */
