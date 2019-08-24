/*
 * EgoCentricNetworkModel.h
 *
 *  Created on: May 27, 2010
 *      Author: duyvu
 */

#ifndef EGOCENTRICNETWORKMODEL_H_
#define EGOCENTRICNETWORKMODEL_H_

#include "DataTypes.h"
#include "nodal/data/EgoCentricNetworkData.h"

namespace ndip {

class EgoCentricNetworkModel {
protected:
	// the number of nodes
	unsigned int numOfVertices;
	// a time-ordered list of nodes and corresponding join times. It is used to compute the age time.
	ExposureTimeVectorType nodeJoiningTimes;

	// the number of edges
	unsigned int numOfEdges;
	// a time-ordered list of edges and corresponding join times
	EdgeTimeVectorType edgeJoiningTimes;

	// assume the observation time is the time of the first edge
	// recorded on the continuous-time scale. Otherwise, the first element of vectorOfComingVertices
	// is the list of vertices coming from the observation time to the first edge after this time.
	// The assumption guarantees the synchronization of the new node lists and the update lists.
	double observationTime;

	EgoCentricNetworkData* data;

	void prepareExposureIndicators(bool* exposureIndicators);
	void computeMeanStandardDeviation(double& mean, double& stDeviation, list<
			double>& values);

	virtual void initializeMore() = 0;
public:
	EgoCentricNetworkModel();
	virtual ~EgoCentricNetworkModel();

	// read data of nodes and edges from files
	void setNetworkData(char* exposureTimeFile, char* edgeTimeFile);

	// choose which ego centric model to use
	// for a set of network statistics of interest,
	// we will build an corresponding ego centric network data
	// the likelhood and ML estimation code will be independent with
	// the selected ego centric network data
	void setModel(int _modelType, double _observationTime, const std::vector<
			string>& listOfNodalDataFiles);

	// timing covariates such as age or gap will be handled in
	// different subclasses of EgoCentricNetworkModel
	// Each subclass will compute this likelihood in different ways
	// depending on age or gap time is used
	// this method will be used for free-gradient optimization methods
	// such as Nelder-Mead
	virtual double computeLogLikelihood(const DoubleVector& beta,
			bool isUsingLBF) = 0;
	// this method will be used for gradient-based optimization methods
	virtual double computeLogLikelihoodAndScoreVector(DoubleVector& u,
			const DoubleVector& beta) = 0;
	// this method will be used for Newton-Raphson optimization methods
	virtual void computeScoreVectorAndInformationMatrix(DoubleVector& u,
			DoubleMatrix& I, const DoubleVector& beta) = 0;
	// get ML estimates of parameters
	void getMLE(DoubleVector& beta, DoubleMatrix& cov,
			double firstPercentageStepLength, int maxNRSteps,
			double scoreTolerance);
	// a helper function for Nelder-Mead optimization method
	// in case we only use function values to get ML estimates
	static double computeMinusLogLLH(int n, const double *x, double *grad,
			void *my_func_data) {
		EgoCentricNetworkModel* model = (EgoCentricNetworkModel*) my_func_data;
		DoubleVector beta(n);
		for (int k = 0; k < n; k++)
			beta(k) = x[k];
		if (grad) {
			DoubleVector u(n);
			DoubleMatrix I(n, n);
			model->computeScoreVectorAndInformationMatrix(u, I, beta);
			for (int k = 0; k < n; k++)
				grad[k] = -u[k];
		}
		double result = model->computeLogLikelihood(beta, true);

		cout << "x = ";
		for (int k = 0; k < n; k++)
			cout << x[k] << ", ";
		cout << endl;

		return -result; // nlopt works on minimizing problems
	}

	virtual void
			getNodalCovariateMatrix(DoubleMatrix& covariates, double time) = 0;

	virtual void
			computeCummulativeBaselineHazard(
					VectorOfCummulativePoints& cumBH,
					const DoubleVector& beta) = 0;

	virtual void computeCoxSnellResiduals(
			VectorOfCummulativePoints& csResiduals,
			const DoubleVector& beta) = 0;
	virtual void computeMartingaleResiduals(
			VectorOfCummulativePoints& martingaleResiduals,
			VectorOfCummulativePoints& timeMartingaleResiduals,
			const DoubleVector& beta) = 0;

	virtual void computeSchoenfeldResiduals(DoubleMatrix& schoenfeldResiduals,
			const DoubleVector& beta, bool isPlusBeta) = 0;

	virtual void getTrainingRanks(DoubleVector& ranks,
			const DoubleVector& beta, const std::string& tieMethod) = 0;

	virtual void simulateNetwork(Graph& graph, DoubleMatrix& covariates,
			const DoubleVector& beta) = 0;
	virtual void simulateNetworkWithBothEnds(Graph& graph,
			DoubleMatrix& covariates, const DoubleVector& betaOut,
			const DoubleVector& betaIn) = 0;

	unsigned int getNumOfVertices() {
		return numOfVertices;
	}
};

}

#endif /* EGOCENTRICNETWORKMODEL_H_ */
